############################################################################
# ReDGAP - Structural Variant (SV) filter code for Whole Genome and Panel Gene Analysis with files annotated by TCAG SickKids pipeline 1.7.3.
############################################################################

messages <- c(
  "You are currently using the Structural Variant (SV) filter code for Whole Genome and Panel Gene Analysis annotated by TCAG Sickkids pipeline 1.7.3. ",
  "Select your patient file, and have a panel of genes as a text file if you want to use one."
)
for (msg in messages) { cat(msg, "\n"); Sys.sleep(1) }

install_and_load <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

cran_packages <- c("data.table", 
                   "conflicted", 
                   "tidyverse", 
                   "writexl",
                   "dplyr", 
                   "readr", 
                   "readxl", 
                   "openxlsx")
lapply(cran_packages, install_and_load)

###
# HELPER FUNCTIONS
###

# Running log of step -> row count, becomes the SV_Analysis_Pipeline sheet
Counter <- list(Step = character(0), Count = numeric(0))
log_step <- function(step, n) {
  Counter$Step  <<- c(Counter$Step, step)
  Counter$Count <<- c(Counter$Count, n)
}

printFilter <- function(df, param, numLines) {
  cat("Number of SVs removed for", param, "is:", numLines - nrow(df), "\n")
  cat(nrow(df), "SVs remain.\n")
  nrow(df)
}

# Text -> number, quietly (anything non-numeric becomes NA)
num <- function(x) suppressWarnings(as.numeric(x))

# TRUE where a value passes a threshold rule. Missing values PASS, because an SV with no frequency/overlap information is treated as rare, not removed.
passes_rule <- function(x, op, cutoff) {
  val  <- num(x)
  keep <- if (op == "<=") val <= cutoff else val < cutoff
  is.na(val) | keep
}

# ClinVar helpers (used for both scoring and for building the clinical list)
clinvar_path <- function(df) {
  !is.na(df$ClinVar_ClinSign) &
    grepl("Pathogenic|Likely_pathogenic", df$ClinVar_ClinSign, ignore.case = TRUE) &
    !grepl("Conflicting_classifications_of_pathogenicity", df$ClinVar_ClinSign, ignore.case = TRUE)
}
clinvar_vus <- function(df) {
  !is.na(df$ClinVar_ClinSign) & grepl("Uncertain", df$ClinVar_ClinSign, ignore.case = TRUE)
}

# TRUE where the SV is closer than `cutoff` bp to the nearest exon boundary
near_exon <- function(df, col, cutoff) {
  d <- abs(num(df[[col]]))
  !is.na(d) & d < cutoff
}

############################################################################
# IMPORT DATA
############################################################################
print('Input patient file')
fileName <- file.choose()
if (!tolower(tools::file_ext(fileName)) %in% c('tsv', 'csv', 'txt', 'dsv')) {
  stop("Unsupported file format", call. = FALSE)
}
outName  <- sub("\\..*$", "", basename(fileName))
repTimes <- length(strsplit(readLines(fileName, n = 1), ',')[[1]])
dfall <- fread(fileName, colClasses = rep('character', repTimes),
                na.strings = "NA", data.table = FALSE)
names(dfall) <- sapply(strsplit(names(dfall), ':'), function(i) i[[length(i)]])

numLines <- nrow(dfall)
cat("Initial row count:", numLines, "\n")
log_step("Initial number of SVs", numLines)

# Convert empty value to NA
dfall[] <- lapply(dfall, function(x) {x[x==""] <- NA; x})

# Variant Score Column
dfall$VScore <- 0   # Create Variant score column that will be added to later

############################################################################
# OPTIONS (both questions are asked now, so the rest of the script can run without waiting)
############################################################################
if (interactive()) {
  ans2 <- tolower(trimws(readline("Do you want to remove SVs overlapping polymorphic regions by >= 70% (dirtyRegion)? [yes/no]: ")))
  ans1 <- tolower(trimws(readline("Do you want to use a panel? [yes/no]: ")))
}

############################################################################
# ADD GENE PANEL (Optional)
############################################################################
panel         <- NULL
panel_pattern <- NULL

if (exists("ans1") && ans1 %in% c("y", "yes")) {
  if (!"gene_symbol" %in% colnames(dfall)) {
    cat("No gene_symbol column found, so the panel will not be used.\n")
  } else {
    print('Input gene panel file')
    panelFile <- file.choose()
    print(panelFile)
    if (!tolower(tools::file_ext(panelFile)) %in% c('tsv', 'csv', 'txt', 'dsv')) {
      stop("Unsupported file format", call. = FALSE)
    }
    panel <- read.table(panelFile, header = TRUE, colClasses = "character")
    cat("Number of panel genes:", nrow(panel), "\n")

    # Gene names are in the first column; match them whether separated by | or ;
    panel_genes   <- panel[[1]]
    panel_pattern <- paste(paste0("(^|[;|])", panel_genes, "($|[;|])"), collapse = "|")
  }
}

############################################################################
# VARIANT SCORING
# All rows are scored BEFORE any filtering or subsetting, so every sheet shows the same final VScore and duplicate rows can be removed reliably later
############################################################################

###
# SCORE BASED ON VARIANTS
# EXAMPLE: list(need = "<column that must exist>", points = <how much to add>, test = function(df) <TRUE/FALSE per row>)
###
simple_score_rules <- list(
  list(need = "ClinGen_SV_pctOvlp_max", points = 10, test = function(df)
    num(df$ClinGen_SV_pctOvlp_max) >= 30),
  list(need = "decipher_SV_pctOvlp_max", points = 10, test = function(df)
    num(df$decipher_SV_pctOvlp_max) >= 30),
  list(need = "ClinVar_ClinSign", points = 10, test = function(df) clinvar_path(df)),
  list(need = "ClinVar_ClinSign", points = 2,  test = function(df) clinvar_vus(df)),
  list(need = "cds_symbol", points = 20, test = function(df)
    !is.na(df$cds_symbol)),
  list(need = "numberOfGeneSymbols", points = 10, test = function(df)
    num(df$numberOfGeneSymbols) >= 3),
  list(need = "numberOfGeneSymbols", points = 5, test = function(df)
    num(df$numberOfGeneSymbols) == 2),
  list(need = "Fantom5DPI", points = 5, test = function(df)
    !is.na(df$Fantom5DPI) & df$Fantom5DPI == "1"),
  list(need = "Fantom5DPI", points = 10, test = function(df)
    !is.na(df$Fantom5DPI) & df$Fantom5DPI == "2"),
  # Near an exon boundary: +2 within 40 bp, and a further +3 within 10 bp (so +5 in total)
  list(need = "nearestLeftExonDistance",  points = 2, test = function(df) near_exon(df, "nearestLeftExonDistance", 40)),
  list(need = "nearestLeftExonDistance",  points = 3, test = function(df) near_exon(df, "nearestLeftExonDistance", 11)),
  list(need = "nearestRightExonDistance", points = 2, test = function(df) near_exon(df, "nearestRightExonDistance", 40)),
  list(need = "nearestRightExonDistance", points = 3, test = function(df) near_exon(df, "nearestRightExonDistance", 11))
)

for (rule in simple_score_rules) {
  if (rule$need %in% colnames(dfall)) {
    idx <- rule$test(dfall)
    idx[is.na(idx)] <- FALSE
    dfall$VScore[idx] <- dfall$VScore[idx] + rule$points
  }
}

###
# ADD SCORE BASED ON GENE PANEL
###
if (!is.null(panel_pattern)) {
  idpanel <- grepl(panel_pattern, dfall$gene_symbol)
  dfall$VScore[idpanel] <- dfall$VScore[idpanel] + 50
}

dfall <- dfall[order(-dfall$VScore), ]

############################################################################
# CLINICAL DATABASES (ClinGen, DECIPHER, ClinVar)
# Built BEFORE the frequency filters, so clinically relevant SVs are never lost to them
############################################################################
clinical_defs <- list(
  list(need = "ClinGen_SV_pctOvlp_max",
       label = "Number of SVs overlapping 30%+ ClinGen SVs",
       test = function(df) num(df$ClinGen_SV_pctOvlp_max) >= 30),
  list(need = "decipher_SV_pctOvlp_max",
       label = "Number of SVs overlapping 30%+ decipher SVs",
       test = function(df) num(df$decipher_SV_pctOvlp_max) >= 30),
  list(need = "ClinVar_ClinSign",
       label = "Number of SVs overlapping ClinVar P, LP, VUS",
       test = function(df) clinvar_path(df) | clinvar_vus(df))
)

clinical <- NULL
for (d in clinical_defs) {
  if (d$need %in% colnames(dfall)) {
    idx <- d$test(dfall)
    idx[is.na(idx)] <- FALSE
    part <- dfall[idx, ]
    cat(d$label, ":", nrow(part), "\n")
    log_step(d$label, nrow(part))
    clinical <- rbind(clinical, part)
  }
}
if (!is.null(clinical)) {
  clinical <- clinical[!duplicated(clinical), ]
  clinical <- clinical[order(-clinical$VScore), ]
  cat("Unique SVs overlapping ClinGen, decipher, or ClinVar:", nrow(clinical), "\n")
  log_step("Number of unique SVs overlapping 30%+ ClinGen, decipher, or ClinVar SVs", nrow(clinical))
}

clinvar_path_SV <- if ("ClinVar_ClinSign" %in% colnames(dfall)) 
  dfall[clinvar_path(dfall), ] else NULL       

############################################################################
# FILTERING FOR ALLELE FREQUENCY
# Frequency columns from the TCAG pipeline are in percent, except gnomAD (fraction)
############################################################################
cat("Keeping SVs with frequency lower than 0.5% in large databases (1% in NCBI nstd186).\n")

# Written as a tribble so it is easy to change the rules by row (more human readable)
af_filters <- tibble::tribble(
  ~col,                            ~op,   ~cutoff, ~label,
  "gnomAD_AF_MAX_max90pctRecOvlp", "<=",  0.005,   "Removed SVs with MAF > 0.5% in gnomADv4.1",
  "otgMantaPctFreq_90pctRecOvlp",  "<=",  0.5,     "Removed SVs with MAF > 0.5% in 1000g+ called with Manta.",
  "otgDellyPctFreq_90pctRecOvlp",  "<=",  0.5,     "Removed SVs with MAF > 0.5% in 1000g+ called with Delly.",
  "otgDgnPctFreq_90pctRecOvlp",    "<=",  0.5,     "Removed SVs with AF > 0.5% in 1000g+ called with Dragen.",
  "all_DgnPctFreq_90pctRecOvlp",   "<=",  0.5,     "Removed SVs with AF > 0.5% in an internal database of ~10600 unrelated samples called with Dragen.",
  "nstd186_pctOvlp",               "<=",  1,       "Removed SVs with MAF > 1% in NCBI curated common SV-nstd186.",
  "DGV_commonPct",                 "<",   80,      "Removed SVs with 80%+ overlap with SVs with a frequency > 1% in DGV."
)


for (i in seq_len(nrow(af_filters))) {
  r <- af_filters[i, ]
  if (r$col %in% colnames(dfall)) {
    dfall <- dfall[passes_rule(dfall[[r$col]], r$op, r$cutoff), ]
    numLines <- printFilter(dfall, r$col, numLines)
    log_step(r$label, nrow(dfall))
  }
}

############################################################################
# QUALITY FILTER (Optional)
############################################################################
if (exists("ans2") && ans2 %in% c("y", "yes") && "dirtyRegion_pctOvlp" %in% colnames(dfall)) {
  dfall <- dfall[passes_rule(dfall$dirtyRegion_pctOvlp, "<", 70), ]
  numLines <- printFilter(dfall, "dirtyRegion_pctOvlp", numLines)
  log_step("Removed SVs that overlap polymorphic regions > 70% of their size.", nrow(dfall))
}

############################################################################
# CREATE GENE PANEL SHEET
############################################################################
panelSV <- NULL
if (!is.null(panel_pattern)) {
  panelSV <- dfall[grepl(panel_pattern, dfall$gene_symbol), ]
  cat(nrow(panelSV), "SVs from genes in the panel were saved.\n")
  log_step("SVs from panel were saved", nrow(panelSV))
}

############################################################################
# MAKING SUB DATAFRAMES
############################################################################
dfall <- dfall[order(-dfall$VScore), ]
rare  <- dfall

###
# Keep genic variants (map to at least one gene)
###
genic <- dfall
if ("gene_symbol" %in% colnames(dfall)) {
  genic <- dfall[!is.na(dfall$gene_symbol), ]
  numLines <- printFilter(genic, "non-genic SVs", nrow(dfall))
  log_step("Rare SVs mapped in genes", nrow(genic))
}

###
# Coding region (CDS) list
###
CDS <- NULL
if ("cds_symbol" %in% colnames(genic)) {
  CDS <- genic[!is.na(genic$cds_symbol), ]
  cat(nrow(CDS), "rare SVs overlapping CDS were saved.\n")
  log_step("Rare SVs overlapping CDS", nrow(CDS))
}

###
# Splice list (within 20 bp of a canonical exon boundary)
###
spli <- NULL
exon_cols <- intersect(c("nearestLeftExonDistance", "nearestRightExonDistance"), colnames(genic))
if (length(exon_cols) > 0) {
  idx  <- Reduce(`|`, lapply(exon_cols, function(cl) near_exon(genic, cl, 20)))
  spli <- genic[idx, ]
  cat(nrow(spli), "rare SVs within 20 bp of canonical splicing sites were saved.\n")
  log_step("Rare SVs <20 bp of canonical splicing sites", nrow(spli))
}

###
# Promoter list (FANTOM5)
###
F5 <- NULL
if ("Fantom5DPI" %in% colnames(genic)) {
  F5 <- genic[!is.na(genic$Fantom5DPI) & genic$Fantom5DPI != "0", ]
  cat(nrow(F5), "rare SVs in promoters (FANTOM5 != 0) were saved.\n")
  log_step("Rare SVs in promoters (FANTOM5 != 0)", nrow(F5))
}

############################################################################
# ASSEMBLE THE HIGH IMPACT SV LIST
############################################################################
strong_parts <- list(CDS, spli, F5, clinvar_path_SV)
strong_parts <- strong_parts[!sapply(strong_parts, is.null)]

if (length(strong_parts) > 0) {
  strong_SV <- do.call(rbind, strong_parts)
  strong_SV <- strong_SV[!duplicated(strong_SV), ]
  strong_SV <- strong_SV[order(-strong_SV$VScore), ]
} else {
  strong_SV <- dfall[0, ]
}
log_step("Strong SVs: SVs in CDS, splicing or regulatory regions, or P/LP in ClinVar", nrow(strong_SV))

############################################################################
# WRITE OUT THE DATA
############################################################################
FilterRows <- data.frame(Step = Counter$Step, Count = Counter$Count, stringsAsFactors = FALSE)

sheets <- list(
  SV_Analysis_Pipeline = FilterRows,
  High_impact_SV       = strong_SV,
  rare_SV              = rare,
  Genic_SV             = genic,
  clin_db_SV           = clinical    # NULL (and skipped) if no ClinGen/decipher/ClinVar columns exist
)
if (!is.null(panelSV)) {
  sheets$Panel_SV   <- panelSV
  sheets$Panel_list <- panel
}
sheets <- sheets[!sapply(sheets, is.null)]

wb <- openxlsx::createWorkbook()
for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  writeData(wb, sheet = nm, x = sheets[[nm]])
}

final_file <- paste0(outName, "_High_impact_SV.xlsx")
if (file.exists(final_file)) {
  counter <- 1
  dup_file <- paste0(outName, "_", counter, "_High_impact_SV.xlsx")
  while (file.exists(dup_file)) {
    counter <- counter + 1
    dup_file <- paste0(outName, "_", counter, "_High_impact_SV.xlsx")
  }
  saveWorkbook(wb, file = dup_file)
  cat(paste0("File saved as '", dup_file, "'.\n"))
} else {
  saveWorkbook(wb, file = final_file)
  cat(paste0(final_file, " saved successfully!\n"))
}
