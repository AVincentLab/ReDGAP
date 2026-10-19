############################################################################
# ReDGAP - Single Nucleotide Variant (SNV) filter code for Whole Genome and Panel Gene Analysis with files annotated by TCAG SickKids pipeline 28.5.
############################################################################

messages <- c(
  "You are currently using the Single Nucleotide Variant (SNV) filter code for Whole Genome and Panel Gene Analysis annotated by TCAG Sickkids pipeline 28.5. "
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

cran_packages <- c("data.table", "conflicted", "tidyverse", "writexl",
                    "dplyr", "readr", "readxl", "openxlsx")
lapply(cran_packages, install_and_load)

###
# PHENOTYPE TERMS (used by both MPO and HPO scoring rules)
###
PHENO_TERMS <- c("vision", "eye", "retina", "macula")

###
# HELPER FUNCTIONS
###

# Running log of step -> row count, becomes the SNV_Filter_History sheet
Counter <- list(Step = character(0), Count = numeric(0))
log_step <- function(step, n) {
  Counter$Step  <<- c(Counter$Step, step)
  Counter$Count <<- c(Counter$Count, n)
}

printFilter <- function(df, param, numLines) {
  cat("Number of variants removed for", param, "is:", numLines - nrow(df), "\n")
  cat(nrow(df), "variants remain.\n")
  nrow(df)
}

############################################################################
# OPTIONS FOR QUALITY FILTERS (questions are asked now, so the rest of the script can run without waiting)
############################################################################
if (interactive()) {
  ans2 <- tolower(trimws(readline("Do you want to keep only variants with FILTER == PASS? [yes/no]: ")))
  ans3 <- tolower(trimws(readline("Do you want to keep only variants with read depth between 5-120? [yes/no]: ")))
}

############################################################################
# IMPORT PATIENT DATA
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
log_step("Initial Number of Variants", numLines)

# Convert empty value to NA (dfall[dfall==""] <- NA)
dfall[] <- lapply(dfall, function(x) {x[x==""] <- NA; x})

# Create splice AI delta column
splice_cols <- c("spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL")
if (all(splice_cols %in% colnames(dfall))) {
  dfall$sAI_delta <- rowSums(
    sapply(splice_cols, function(cl) suppressWarnings(as.numeric(dfall[[cl]]))),
    na.rm = TRUE
  )
}
# Create promoter AI column
if ("promoterAI_score" %in% colnames(dfall)) {
  dfall$promoterAInum <- suppressWarnings(as.numeric(dfall$promoterAI_score))
  dfall$promoterAInum[is.na(dfall$promoterAInum)] <- 0
}
# Create new DP column
if (all(c("AD_ALT", "AD_REF") %in% colnames(dfall))) {
  dfall$newDP <- suppressWarnings(as.numeric(dfall$AD_ALT) + as.numeric(dfall$AD_REF))
}

# Variant Score Column
dfall$VScore <- 0   # Create Variant score column that will be added to later

############################################################################
# ADD GENE PANEL (Optional)
############################################################################
panel_pattern <- NULL

if (interactive()) {
  ans1 <- tolower(trimws(readline("Do you want to use a panel? [yes/no]: ")))
}
if (exists("ans1") && ans1 %in% c("y", "yes") && "gene_symbol" %in% colnames(dfall)) {
  print('Input gene panel file')
  panelFile <- file.choose()
  print(panelFile)
  if (!tolower(tools::file_ext(panelFile)) %in% c('tsv', 'csv', 'txt', 'dsv')) {
    stop("Unsupported file format", call. = FALSE)
  }
  panel <- read.table(panelFile, header = TRUE, colClasses = "character")
  cat("Number of panel genes:", nrow(panel), "\n")

  panel_genes   <- panel[[1]]
  panel_pattern <- paste(paste0("(^|[;|])", panel_genes, "($|[;|])"), collapse = "|")
}

############################################################################
# FILTERING FOR ALLELE FREQUENCY
############################################################################
cat("Keeping Allele freq lower than 0.5% in large databases and 1% in internal databases.\n")

af_filters <- list(
  list(col = "A1000g_freq_max",                 max = 0.005),
  list(col = "gnomadV4_exome_fafmax_faf95_max", max = 0.005),
  list(col = "gnomadV4_exome_AF_grpmax",        max = 0.005),
  list(col = "gnomadV4_genome_fafmax_faf95_max",max = 0.005),
  list(col = "gnomadV4_genome_AF_grpmax",       max = 0.005),
  list(col = "hsint_AF_all",                    max = 0.01),
  list(col = "hsint_AF_EUR",                    max = 0.01),
  list(col = "gnomadV4_exome_AF_remaining",     max = 0.005),
  list(col = "gnomadV4_genome_AF_remaining",    max = 0.005),
  list(col = "gnomadV4_genome_nhomalt",         max = 7),
  list(col = "gnomadV4_exome_nhomalt",          max = 7)
)

for (f in af_filters) {
  if (f$col %in% colnames(dfall)) {
    val  <- suppressWarnings(as.numeric(dfall[[f$col]]))
    keep <- is.na(val) | val <= f$max
    dfall <- dfall[keep, ]
    numLines <- printFilter(dfall, f$col, numLines)
    log_step(paste0("Removed SNVs failing ", f$col, " threshold (<= ", f$max, ")"), nrow(dfall))
  }
}

############################################################################
# QUALITY FILTERS (Optional)
############################################################################
if (exists("ans2") && ans2 %in% c("y", "yes") && "FILTER" %in% colnames(dfall)) {
  dfall <- dfall[dfall$FILTER %in% "PASS", ]
  log_step("Variants with FILTER != PASS were removed", nrow(dfall))
}

if (exists("ans3") && ans3 %in% c("y", "yes") && "newDP" %in% colnames(dfall)) {
  dp <- suppressWarnings(as.numeric(dfall$newDP))
  dfall <- dfall[!is.na(dp) & dp >= 5 & dp <= 120, ]
  log_step("Quality Filter: 5 <= Depth <= 120", nrow(dfall))
}

############################################################################
# VARIANT SCORING
############################################################################

###
# PREDICTOR VOTE (in-silico predictors)
# +1 is added to VScore for passing any threshold AND later a subset will be made by passing two or more thresholds
###

# Written as a tribble so it is easy to change the predictor rules by row (more human readable)
predictor_rules <- tibble::tribble(
  ~col,                   ~direction, ~cutoff, ~description,
  "phylopMam_avg",        ">=",        1.5,    "Mammalian conservation (phyloP)",
  "phylopVert100_avg",    ">=",        4.0,    "100-way vertebrate conservation (phyloP)",
  "PROVEAN_score",        "<=",       -2.5,    "PROVEAN (more negative = more damaging)",
  "CADD_phred",           ">=",       15.0,    "CADD Phred-scaled score",
  "alpha_missense_score", ">=",        0.8,    "AlphaMissense pathogenicity probability",
  "primateai3d_score",    ">=",        0.8,    "PrimateAI-3D score",
  "MPC_score",            ">=",        2.0,    "Missense badness/constraint (MPC)",
  "REVEL_score",          ">=",        0.5,    "REVEL ensemble score"
)

# Sanity check the table itself, so a typo (e.g. ">" instead of ">=") fails
# loudly right away instead of silently doing the wrong thing later.
stopifnot(
  "direction must be '>=' or '<='" = all(predictor_rules$direction %in% c(">=", "<=")),
  "cutoff must be numeric"         = is.numeric(predictor_rules$cutoff)
)

compute_predictor_votes <- function(df, rules) {
  present <- rules[rules$col %in% colnames(df), , drop = FALSE]
  if (nrow(present) == 0) return(rep(0, nrow(df)))
  votes <- vapply(seq_len(nrow(present)), function(i) {
    x <- suppressWarnings(as.numeric(df[[present$col[i]]]))
    if (present$direction[i] == ">=") x >= present$cutoff[i] else x <= present$cutoff[i]
  }, logical(nrow(df)))
  if (is.null(dim(votes))) votes <- matrix(votes, ncol = 1)
  rowSums(votes, na.rm = TRUE)
}

dfall$predictor_votes <- compute_predictor_votes(dfall, predictor_rules)
dfall$VScore[dfall$predictor_votes >= 1] <- dfall$VScore[dfall$predictor_votes >= 1] + 1

###
# SCORE BASED ON VARIANTS
# EXAMPLE: list(need = "<column that must exist>", points = <how much to add>, test = function(df) <TRUE/FALSE per row>)
###
simple_score_rules <- list(
  list(need = "Clinvar_SIG", points = 30, test = function(df)
    !is.na(df$Clinvar_SIG) &
      grepl("Pathogenic|Likely_pathogenic", df$Clinvar_SIG, ignore.case = TRUE) &
      !grepl("Conflicting_classifications_of_pathogenicity", df$Clinvar_SIG, ignore.case = TRUE)),
  list(need = "Clinvar_SIG", points = 2, test = function(df)
    !is.na(df$Clinvar_SIG) & grepl("Uncertain", df$Clinvar_SIG, ignore.case = TRUE)),
  list(need = "effect_priority", points = 15, test = function(df)
    !is.na(df$effect_priority) &
      df$effect_priority %in% c('frameshift deletion', 'frameshift insertion', 'stopgain')),
  list(need = "effect_priority", points = 5, test = function(df)
    !is.na(df$effect_priority) &
      df$effect_priority %in% c('nonframeshift insertion', 'nonframeshift deletion', 'startloss', 'stoploss')),
  list(need = "sAI_delta", points = 1, test = function(df)
    !is.na(df$sAI_delta) & df$sAI_delta >= 0.2),
  list(need = "sAI_delta", points = 0.5, test = function(df)
    !is.na(df$sAI_delta) & df$sAI_delta >= 0.1 & df$sAI_delta < 0.2),
  list(need = "MPO", points = 0.5, test = function(df)
    !is.na(df$MPO) & grepl(paste(PHENO_TERMS, collapse = "|"), df$MPO, ignore.case = TRUE)),
  list(need = "HPO", points = 0.5, test = function(df)
    !is.na(df$HPO) & grepl(paste(PHENO_TERMS, collapse = "|"), df$HPO, ignore.case = TRUE)),
  list(need = "Fantom5DPI", points = 0.5, test = function(df)
    !is.na(df$Fantom5DPI) & df$Fantom5DPI == "2"),
  list(need = "promoterAInum", points = 0.5, test = function(df)
    !is.na(df$promoterAInum) & abs(df$promoterAInum) > 0.1)
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
if (!is.null(panel_pattern) && "gene_symbol" %in% colnames(dfall)) {
  idpanel <- grepl(panel_pattern, dfall$gene_symbol)
  dfall$VScore[idpanel] <- dfall$VScore[idpanel] + 50
}

dfall <- dfall[order(-dfall$VScore), ]

############################################################################
# CLINVAR (create a dataframe with clinically relevant variants)
############################################################################
mClinVar <- NULL
if ("Clinvar_SIG" %in% colnames(dfall)) {
  idcv1_early <- !is.na(dfall$Clinvar_SIG) &
    grepl("Pathogenic|Likely_pathogenic", dfall$Clinvar_SIG, ignore.case = TRUE) &
    !grepl("Conflicting_classifications_of_pathogenicity", dfall$Clinvar_SIG, ignore.case = TRUE)
  idcv2_early <- !is.na(dfall$Clinvar_SIG) &
    grepl("Uncertain", dfall$Clinvar_SIG, ignore.case = TRUE)

  mClinVar <- dfall[idcv1_early | idcv2_early, ]
  mClinVar <- mClinVar[!duplicated(mClinVar), ]

  log_step("Number of SNVs overlapping ClinVar P, LP, VUS", nrow(mClinVar))
}

############################################################################
# CREATE GENE PANEL SHEET
############################################################################
panelSNV <- NULL
if (!is.null(panel_pattern)) {
  panelSNV <- dfall[grepl(panel_pattern, dfall$gene_symbol), ]
  cat(nrow(panelSNV), "SNVs from genes in the panel were saved.\n")
  log_step("SNVs from panel were saved", nrow(panelSNV))
}

############################################################################
# MAKING SUB DATAFRAMES
############################################################################

# Column changes to be applied to full set so row matching later removes duplicates
dfall$dbscSNV_ADA_SCORE <- suppressWarnings(as.numeric(dfall$dbscSNV_ADA_SCORE))
dfall$dbscSNV_ADA_SCORE[is.na(dfall$dbscSNV_ADA_SCORE)] <- 0
dfall$dbscSNV_RF_SCORE  <- suppressWarnings(as.numeric(dfall$dbscSNV_RF_SCORE))
dfall$dbscSNV_RF_SCORE[is.na(dfall$dbscSNV_RF_SCORE)] <- 0

###
# ClinVar List (Pathogenic/Likely pathogenic)
###
ClinVar_strong <- NULL
if ("Clinvar_SIG" %in% colnames(dfall)) {
  idcv1 <- !is.na(dfall$Clinvar_SIG) &
    grepl("Pathogenic|Likely_pathogenic", dfall$Clinvar_SIG, ignore.case = TRUE) &
    !grepl("Conflicting_classifications_of_pathogenicity", dfall$Clinvar_SIG, ignore.case = TRUE)
  ClinVar_strong <- dfall[idcv1, ]
}

###
# Promoter and Enhancer Lists (Regulatory/non-coding -> promoters and enhancers are often upstream)
###
promoterAI <- NULL
if ("promoterAInum" %in% colnames(dfall)) {
  promoterAI <- dfall[!is.na(dfall$promoterAInum) & abs(dfall$promoterAInum) > 0.1, ]
  cat(nrow(promoterAI), "SNVs with absolute promoterAI values > 0.1 were saved.\n")
  log_step("SNVs with absolute promoterAI values > 0.1 were saved", nrow(promoterAI))
}

FANTOM5 <- NULL
if ("Fantom5DPI" %in% colnames(dfall)) {
  FANTOM5 <- dfall[!is.na(dfall$Fantom5DPI) & dfall$Fantom5DPI == "2", ]
  cat(nrow(FANTOM5), "SNVs overlapping robust FANTOM5 motifs were saved.\n")
  log_step("SNVs overlapping robust FANTOM5 motifs were saved", nrow(FANTOM5))
}

###
# Remove the non-genic regions from dfall
###
RemoveForGenic <- c('downstream', 'intergenic', 'ncRNA_exonic', 'ncRNA_intronic',
                     'ncRNA_splicing', 'ncRNA_UTR3', 'ncRNA_UTR5', 'upstream',
                     'upstream;downstream', 'UTR3', 'UTR5', 'UTR5;UTR3',
                     'ncRNA_exonic;ncRNA_splicing')
genic_dfall <- dfall[!dfall$typeseq_priority %in% RemoveForGenic, ]
log_step("Genic Variants", nrow(genic_dfall))

###
# Splice List
###
cat("Splice AI filtering - cutoff is 0.1.\n")
spliceAI_data <- genic_dfall[!is.na(genic_dfall$sAI_delta) & genic_dfall$sAI_delta >= 0.1, ]

dbscSNV <- genic_dfall[genic_dfall$dbscSNV_ADA_SCORE >= 0.6 | genic_dfall$dbscSNV_RF_SCORE >= 0.6, ]

Splice <- rbind(spliceAI_data, dbscSNV)
Splice <- Splice[!duplicated(Splice), ]
cat(nrow(Splice), "variants passed SpliceAI (>=0.1) or dbscSNV ADA/RF (>=0.6) splice filtering.\n")
log_step("Kept SpliceAI sum-of-delta >0.1 and dbscSNV (ADA or RF) >=0.6", nrow(Splice))

###
# Effect-priority List (loss-of-function and other qualifying nonframeshift indels)
###
LOFkeep <- c('frameshift deletion', 'frameshift insertion', 'stoploss', 'stopgain', 'startloss', 'startgain')
LOF <- genic_dfall[genic_dfall$effect_priority %in% LOFkeep, ]

keep2 <- c('nonframeshift insertion', 'nonframeshift deletion')
LOF2 <- genic_dfall[genic_dfall$effect_priority %in% keep2, ]
LOF2$dbsnp_region[is.na(LOF2$dbsnp_region)] <- 0
LOF2 <- LOF2[str_count(as.character(LOF2$dbsnp_region), ',') <= 4, ]
LOF2$dbsnp_region[LOF2$dbsnp_region == 0] <- NA

effect_data <- rbind(LOF, LOF2)
effect_data <- effect_data[!duplicated(effect_data), ]
log_step("Saved LOF (stopgain, frameshift), inframe, startloss, stoploss, startgain", nrow(effect_data))

###
# Predictor vote subset (>=2 of the 8 predictors clear cutoff)
###
pathogenicfilter_data <- genic_dfall[genic_dfall$predictor_votes >= 2, ]
log_step("Saved genic variants with two predictive scores over set threshold", nrow(pathogenicfilter_data))

############################################################################
# ASSEMBLE THE STRONG SNV LIST
############################################################################
strong_parts <- list(ClinVar_strong, effect_data, pathogenicfilter_data, FANTOM5, Splice, promoterAI)
strong_parts <- strong_parts[!sapply(strong_parts, is.null)]

strong_SNV <- do.call(rbind, strong_parts)

strong_SNV <- strong_SNV[!duplicated(strong_SNV), ]
#strong_SNV <- strong_SNV[order(strong_SNV$gene_symbol), ]
strong_SNV <- strong_SNV[order(-strong_SNV$VScore), ]
log_step("Predicted as damaging, deleterious, conserved or splicing", nrow(strong_SNV))

############################################################################
# WRITE OUT THE DATA
############################################################################
FilterRows <- data.frame(Step = Counter$Step, Count = Counter$Count, stringsAsFactors = FALSE)

use_panel <- !is.null(panelSNV)

wb <- openxlsx::createWorkbook()
sheet_name_snv        <- "SNV_Filter_History"
sheet_name_all        <- "strong_SNV"
sheet_name_ClinVar    <- "ClinVar"
sheet_name_panel      <- "Panel"
sheet_name_panelgenes <- "Panel_Gene_lt"

addWorksheet(wb, sheet_name_snv)
addWorksheet(wb, sheet_name_all)
addWorksheet(wb, sheet_name_ClinVar)
if (use_panel) {
  addWorksheet(wb, sheet_name_panel)
  addWorksheet(wb, sheet_name_panelgenes)
}

final_file <- paste0(outName, "_Pipeline_SNVs.xlsx")
writeData(wb, sheet = sheet_name_snv, x = FilterRows)
writeData(wb, sheet = sheet_name_all, x = strong_SNV)
writeData(wb, sheet = sheet_name_ClinVar, x = mClinVar)
if (use_panel) {
  writeData(wb, sheet = sheet_name_panel, x = panelSNV)
  writeData(wb, sheet = sheet_name_panelgenes, x = panel)
}

if (file.exists(final_file)) {
  counter <- 1
  dup_file <- paste0(outName, "_", counter, "_Pipeline_SNVs.xlsx")
  while (file.exists(dup_file)) {
    counter <- counter + 1
    dup_file <- paste0(outName, "_", counter, "_Pipeline_SNVs.xlsx")
  }
  saveWorkbook(wb, file = dup_file)
  cat(paste0("File saved as '", dup_file, "'.\n"))
} else {
  saveWorkbook(wb, file = final_file)
  cat(paste0(final_file, " saved successfully!\n"))
}

