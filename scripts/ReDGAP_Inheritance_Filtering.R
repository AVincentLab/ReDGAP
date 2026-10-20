############################################################################
# ReDGAP - Inheritance filter code (recessive / dominant) for SNV, SV and CNV files from the ReDGAP pipeline scripts.
############################################################################
# This script uses gene names to match variants, so every file needs a column called "gene_symbol". If your column has a different name, edit your data before running this.
# Use | to separate gene names for variants that overlap multiple genes.
# It can also be used to find genes shared between affected family members. The shared genes are saved in their entirety, so check the resulting rows for the specific variant matches yourself.

if (!interactive()) {
  stop("This script asks questions as it runs, so run it in an interactive R session (e.g. RStudio).", call. = FALSE)
}

messages <- c(
  "You are currently using the Inheritance Filtering script.",
  "This script is designed to utilize gene names for matching variants.",
  "Note that you will need a column called gene_symbol for this code to work; if you do not have this column name, edit your data before running this."
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
# HELPER FUNCTIONS
###

# Text -> number, quietly (anything non-numeric becomes NA)
num <- function(x) suppressWarnings(as.numeric(x))

# Ask a multiple choice question, repeat until a valid number is typed. Returns the number chosen.
ask_choice <- function(question, options) {
  cat(question, "\n")
  for (i in seq_along(options)) cat("[", i, "] ", options[i], "\n", sep = "")
  repeat {
    ans <- num(trimws(readline("Enter a number: ")))
    if (!is.na(ans) && ans %in% seq_along(options)) return(ans)
    cat("Please pick from the given numbers.\n")
  }
}

printFilter <- function(original_data, modified_data) {
  cat("This action removed", nrow(original_data) - nrow(modified_data), "rows.\n")
  cat("There are currently", nrow(modified_data), "rows remaining.\n")
}

# Make a name safe to use in a file name (Sample ID) or an Excel sheet name
clean_name <- function(x) gsub("[^A-Za-z0-9._-]", "_", gsub(" ", "", x))

# Excel sheet names: max 31 characters, none of [ ] : * ? / \, and no repeats (Excel ignores upper/lower case).
# A repeat gets _2, _3, ... instead of crashing or being silently cut to the same name.
safe_sheet_name <- function(name, existing = character(0)) {
  name <- substr(gsub("[\\[\\]:*?/\\\\]", "_", name, perl = TRUE), 1, 31)
  base <- name
  k <- 1
  while (tolower(name) %in% tolower(existing)) {
    k <- k + 1
    suffix <- paste0("_", k)
    name <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
  }
  name
}

add_sheet <- function(wb, name, data) {
  name <- safe_sheet_name(name, as.character(names(wb)))
  addWorksheet(wb, name)
  writeData(wb, sheet = name, x = data)
  invisible(name)
}

# Save with the usual numbering if the file already exists: <stem>.xlsx, <stem>_1.xlsx, <stem>_2.xlsx ...
save_unique <- function(wb, stem) {
  final_file <- paste0(stem, ".xlsx")
  k <- 0
  while (file.exists(final_file)) {
    k <- k + 1
    final_file <- paste0(stem, "_", k, ".xlsx")
  }
  saveWorkbook(wb, file = final_file)
  cat(paste0("File saved as '", final_file, "'.\n"))
  final_file
}

# One row per (variant, gene) pair. Genes can be separated by | ; or , and blanks/NA are ignored.
# row_id points back to the row of the original data.
gene_map <- function(df) {
  empty <- data.frame(row_id = integer(0), gene = character(0), stringsAsFactors = FALSE)
  if (!"gene_symbol" %in% colnames(df) || nrow(df) == 0) return(empty)
  parts <- strsplit(as.character(df$gene_symbol), "[|;,]")
  map <- data.frame(
    row_id = rep(seq_along(parts), lengths(parts)),
    gene   = trimws(as.character(unlist(parts))),
    stringsAsFactors = FALSE
  )
  map <- map[!is.na(map$gene) & map$gene != "", ]
  map[!duplicated(map), ]
}

all_genes <- function(df) unique(gene_map(df)$gene)

# Excel cells hold at most 32767 characters, so a very long gene list is cut off instead of corrupting the file
genes_to_text <- function(genes) {
  txt <- paste(unique(genes), collapse = ", ")
  if (nchar(txt) > 32000) txt <- paste0(substr(txt, 1, 32000), " ... (cut off, too long for one Excel cell)")
  txt
}

# Running table of what happened to each sheet, becomes the Analysis_Statistics sheet
stats_data <- data.frame()
add_to_stats <- function(sheet_name, sheet_type, total_variants, biallelic_count,
                         genes_involved, genes_retained = NULL, notes = "") {
  new_row <- data.frame(
    Sheet_Name                = sheet_name,
    Sheet_Type                = sheet_type,
    Total_Variants            = total_variants,
    Biallelic_Recessive_Count = biallelic_count,
    Genes_Involved            = genes_to_text(genes_involved),
    Gene_Count_Total          = length(unique(genes_involved)),
    Gene_Count_Retained       = if (is.null(genes_retained)) NA else length(unique(genes_retained)),
    Notes                     = notes,
    stringsAsFactors = FALSE
  )
  stats_data <<- rbind(stats_data, new_row)
}

# Read one file (csv/tsv/txt/xlsx/xls). Returns a named list of data frames, one per sheet chosen.
# Everything is read as text, like the pipeline scripts do, so sheets from different files can be compared.
load_file <- function(taken_names) {
  cat("Please select a file to begin.\n")
  chosen <- file.choose()
  ext    <- tolower(tools::file_ext(chosen))
  out    <- list()

  add_loaded <- function(df, label, fallback) {
    label <- clean_name(label)
    if (label == "") label <- clean_name(fallback)
    label <- safe_sheet_name(label, c(taken_names, names(out)))
    if (!"gene_symbol" %in% colnames(df)) {
      cat("Warning:", label, "has no gene_symbol column, so it can not be used for gene matching.\n")
    }
    if ("VScore" %in% colnames(df)) df$VScore <- num(df$VScore)
    out[[label]] <<- df
    cat("Loaded", label, "with", nrow(df), "rows.\n")
  }

  if (ext %in% c("csv", "tsv", "txt", "dsv")) {
    df <- fread(chosen, colClasses = "character", na.strings = c("NA", ""), data.table = FALSE)
    label <- readline("Please write the file name -i.e the variant file type- (e.g. SNV, SV, CNV): ")
    add_loaded(df, label, sub("\\..*$", "", basename(chosen)))

  } else if (ext %in% c("xlsx", "xls")) {
    sheetlist <- excel_sheets(chosen)
    cat("Choose the sheets you want to load in; each number corresponds with the adjacent sheet:\n")
    for (i in seq_along(sheetlist)) cat(i, " - ", sheetlist[i], "\n")

    keep <- integer(0)
    repeat {
      input <- trimws(readline("Enter a number (or 'q' to finish): "))
      if (tolower(input) == "q") break
      picked <- num(strsplit(input, "[ ,]+")[[1]])
      if (length(picked) > 0 && all(!is.na(picked) & picked >= 1 & picked <= length(sheetlist))) {
        keep <- unique(c(keep, picked))
      } else {
        cat("Invalid input. Please enter a valid number between 1 and", length(sheetlist), ".\n")
      }
    }
    if (length(keep) == 0) {
      cat("No sheets were chosen from this file.\n")
    } else {
      cat("These are the sheets we're keeping:", sheetlist[keep], "\n")
      for (k in keep) {
        df <- as.data.frame(read_xlsx(chosen, sheet = sheetlist[k], col_types = "text"))
        label <- readline(paste0("Please write the file name -i.e the variant file type- for sheet '",
                                 sheetlist[k], "' (e.g. SNV, SV, CNV; blank = sheet name): "))
        add_loaded(df, label, sheetlist[k])
      }
    }

  } else {
    cat("Please pick a file in tabular format!\n")
  }
  out
}

############################################################################
# IMPORT DATA
############################################################################
sample_id <- clean_name(readline("Please write the Sample ID: "))

sheets <- list()   # one data frame per loaded file/sheet, named by what you typed for it
repeat {
  sheets <- c(sheets, load_file(names(sheets)))
  if (length(sheets) == 0) {
    cat("Nothing has been loaded yet, please pick a file.\n")
    next
  }
  if (ask_choice("Do you want to add more files and sheets?", c("Yes", "No")) == 2) {
    cat("Okay, no files will be added.\n")
    break
  }
  cat("OK, adding some more...\n")
}

############################################################################
# CHOOSE INHERITANCE MODE
############################################################################
inheritance_mode <- ask_choice("Choose inheritance mode:", c("Recessive", "Dominant"))

if (inheritance_mode == 1) {

  ##########################################################################
  # RECESSIVE
  ##########################################################################
  wb <- createWorkbook()
  all_genes_saved <- data.frame(Comparison = character(), gene_symbol = character(), stringsAsFactors = FALSE)
  analyzed_genes  <- character(0)
  retained_genes  <- character(0)

  ###
  # STEP 1: RECESSIVE CANDIDATES WITHIN EACH SHEET
  # A variant is kept if (in order of priority):
  #   - its gene has more than one variant in the sheet (possible compound heterozygote), or
  #   - it is homozygous (Zygosity == hom-alt, or GT 1/1 or 1|1 when there is no Zygosity column)
  ###
  cat("Step 1: Finding recessive candidates within each sheet...\n")

  for (nm in names(sheets)) {
    df <- sheets[[nm]]
    cat("Processing", nm, "for recessive candidates...\n")

    reason <- rep(NA_character_, nrow(df))
    genes  <- character(0)
    map    <- gene_map(df)

    # 1. Genes with multiple variants
    if (nrow(map) > 0) {
      gene_counts <- table(map$gene)
      dup_genes   <- names(gene_counts)[gene_counts > 1]
      if (length(dup_genes) > 0) {
        cat("Found", length(dup_genes), "genes with multiple variants in", nm, "\n")
        reason[unique(map$row_id[map$gene %in% dup_genes])] <- "duplicate_gene"
        genes <- c(genes, dup_genes)
      }
    }

    # 2. Homozygous variants (hom-alt in Zygosity), or 3. 1/1 or 1|1 in the GT column [SV]
    if ("Zygosity" %in% colnames(df)) {
      hom <- which(!is.na(df$Zygosity) & df$Zygosity == "hom-alt")
      hom_label <- "hom_alt_zygosity"
      cat("Found", length(hom), "homozygous variants (hom-alt) in", nm, "\n")
    } else {
      GTname <- grep("(^|[_.])GT([_.]|$)", colnames(df), value = TRUE)
      hom <- integer(0)
      if (length(GTname) > 0) {
        hom <- which(!is.na(df[[GTname[1]]]) & df[[GTname[1]]] %in% c("1/1", "1|1"))
        cat("Found", length(hom), "homozygous variants (1/1 or 1|1) in", nm, "\n")
      }
      hom_label <- "homozygous_GT"
    }
    hom_new <- hom[is.na(reason[hom])]   # rows already kept as duplicate genes keep that reason
    reason[hom_new] <- hom_label
    genes <- c(genes, map$gene[map$row_id %in% hom])

    # Gene-level scores (only if the sheet has a VScore column).
    # A variant in several genes gets the highest score of its genes.
    candidates <- df[!is.na(reason), , drop = FALSE]
    candidates$recessive_reason <- reason[!is.na(reason)]
    if ("VScore" %in% colnames(candidates) && nrow(candidates) > 0) {
      cat("Calculating gene-level scores for recessive genes in", nm, "...\n")
      candidates$VScore[is.na(candidates$VScore)] <- 0
      cmap <- gene_map(candidates)
      cmap$VScore <- candidates$VScore[cmap$row_id]
      avg_all  <- tapply(cmap$VScore, cmap$gene, mean)
      avg_top2 <- tapply(cmap$VScore, cmap$gene, function(v) mean(sort(v, decreasing = TRUE)[seq_len(min(2, length(v)))]))
      by_row   <- split(cmap$gene, factor(cmap$row_id, levels = seq_len(nrow(candidates))))
      best <- function(avg) vapply(by_row, function(g) if (length(g) == 0) NA_real_ else max(avg[g]), numeric(1))
      candidates$GeneScore_AvgAll  <- unname(best(avg_all))
      candidates$GeneScore_AvgTop2 <- unname(best(avg_top2))
    }

    # Statistics for the original sheet
    sheet_genes    <- all_genes(df)
    analyzed_genes <- c(analyzed_genes, sheet_genes)
    add_to_stats(sheet_name = nm, sheet_type = "Original", total_variants = nrow(df), biallelic_count = 0,
                 genes_involved = sheet_genes, notes = "Original input data")

    recessive_sheet_name <- paste0(nm, "_recessive")
    if (nrow(candidates) > 0) {
      add_sheet(wb, recessive_sheet_name, candidates)
      kept_genes     <- all_genes(candidates)
      retained_genes <- c(retained_genes, kept_genes)
      add_to_stats(sheet_name = recessive_sheet_name, sheet_type = "Recessive", total_variants = nrow(candidates),
                   biallelic_count = nrow(candidates), genes_involved = genes, genes_retained = kept_genes,
                   notes = "Biallelic/recessive candidates")
      cat("Total recessive candidates found in", nm, ":", nrow(candidates), "\n")
    } else {
      cat("No recessive candidates found in", nm, "\n")
      add_to_stats(sheet_name = recessive_sheet_name, sheet_type = "Recessive", total_variants = 0,
                   biallelic_count = 0, genes_involved = character(0), notes = "No recessive candidates found")
    }
  }

  ###
  # STEP 2: GENES SHARED BETWEEN SHEETS (every pair of sheets is compared)
  # The full rows of every variant in a shared gene are saved, from both sheets.
  ###
  cat("\nStep 2: Finding matches between sheets...\n")

  if (length(sheets) < 2) {
    cat("No other sheet to compare!\n")
  } else {
    for (i in seq_len(length(sheets) - 1)) {
      for (j in (i + 1):length(sheets)) {
        name_i <- names(sheets)[i]
        name_j <- names(sheets)[j]
        comparison_name <- paste0(name_i, "_vs_", name_j)
        cat("Comparing sheet", name_i, "and", name_j, ".\n")

        map_i   <- gene_map(sheets[[i]])
        map_j   <- gene_map(sheets[[j]])
        matches <- base::intersect(map_i$gene, map_j$gene)

        if (length(matches) == 0) {
          cat("There are no matches between", name_i, "and", name_j, ".\n")
          add_to_stats(sheet_name = comparison_name, sheet_type = "Comparison", total_variants = 0,
                       biallelic_count = NA, genes_involved = character(0), notes = "No overlapping genes found")
          next
        }
        cat("There are", length(matches), "gene matches between", name_i, "and", name_j, ".\n")

        sheet_i_matches <- sheets[[i]][unique(map_i$row_id[map_i$gene %in% matches]), , drop = FALSE]
        sheet_j_matches <- sheets[[j]][unique(map_j$row_id[map_j$gene %in% matches]), , drop = FALSE]
        sheet_i_matches$SourceSheet <- name_i
        sheet_j_matches$SourceSheet <- name_j

        gene_symbol_Data <- data.frame(Comparison = comparison_name, gene_symbol = matches, stringsAsFactors = FALSE)
        add_sheet(wb, paste0("Genes_", comparison_name), data.frame(gene_symbol = matches))
        add_sheet(wb, paste0(name_i, "_", comparison_name), sheet_i_matches)
        add_sheet(wb, paste0(name_j, "_", comparison_name), sheet_j_matches)

        all_genes_saved <- rbind(all_genes_saved, gene_symbol_Data)
        add_to_stats(sheet_name = comparison_name, sheet_type = "Comparison",
                     total_variants = nrow(sheet_i_matches) + nrow(sheet_j_matches), biallelic_count = NA,
                     genes_involved = matches, genes_retained = matches,
                     notes = paste0("Overlap between ", name_i, " and ", name_j))
      }
    }
  }

  ###
  # STEP 3: STATISTICS
  ###
  cat("\nStep 3: Creating statistics summary...\n")
  add_sheet(wb, "Analysis_Statistics", stats_data)

  summary_stats <- data.frame(
    Metric = c(
      "Total Original Sheets",
      "Total Recessive Sheets Created",
      "Total Comparison Sheets Created",
      "Total Unique Genes Analyzed",
      "Total Unique Genes Retained",
      "Total Biallelic/Recessive Variants Found"
    ),
    Value = c(
      sum(stats_data$Sheet_Type == "Original"),
      sum(stats_data$Sheet_Type == "Recessive" & stats_data$Total_Variants > 0),
      sum(stats_data$Sheet_Type == "Comparison" & stats_data$Total_Variants > 0),
      length(unique(analyzed_genes)),
      length(unique(c(retained_genes, all_genes_saved$gene_symbol))),
      sum(stats_data$Biallelic_Recessive_Count, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  add_sheet(wb, "Summary_Statistics", summary_stats)

  if (nrow(all_genes_saved) > 0) {
    add_sheet(wb, "All_Unique_Genes", unique(all_genes_saved))
  }

  ###
  # WRITE OUT THE DATA
  ###
  save_unique(wb, paste0(sample_id, "_Inheritance_Filter_Recessive"))

  cat("\nAnalysis complete!\n")
  cat("All results saved in a single Excel file with multiple sheets:\n")
  cat("- Analysis_Statistics: Detailed statistics for each sheet\n")
  cat("- Summary_Statistics: Overall analysis summary\n")
  cat("- Individual recessive candidate sheets\n")
  cat("- Gene overlap comparison sheets\n")

} else {

  ##########################################################################
  # DOMINANT
  # The SNV file (if there is one) is filtered for rare variants, then every sheet is saved.
  ##########################################################################
  cat("Here are the loaded sheets:\n")
  for (i in seq_along(sheets)) cat(i, " - ", names(sheets)[i], "\n")
  snv_choice <- num(trimws(readline("Select the SNV file by number, or type '0' if none: ")))

  # Written as a tribble so it is easy to change the rules by row (more human readable)
  dominant_af_filters <- tibble::tribble(
    ~col,               ~max,    ~label,
    "A1000g_freq_max",  0.0003,  "1000 Genomes",
    "freq_max",         0.0003,  "gnomAD"
  )

  if (!is.na(snv_choice) && snv_choice >= 1 && snv_choice <= length(sheets)) {
    snv_data <- sheets[[snv_choice]]
    present  <- dominant_af_filters[dominant_af_filters$col %in% colnames(snv_data), ]

    if (nrow(present) == 0) {
      cat("The selected file does not contain known population frequency columns.\n")
    } else {
      keep <- rep(TRUE, nrow(snv_data))
      for (i in seq_len(nrow(present))) {
        val  <- num(snv_data[[present$col[i]]])
        keep <- keep & (is.na(val) | val <= present$max[i])   # missing frequency = rare, so it is kept
      }
      filtered_snv <- snv_data[keep, , drop = FALSE]
      sheets[[snv_choice]] <- filtered_snv
      cat("SNV file filtered for allele frequency (<= 0.03% in:", paste(present$label, collapse = ", "), ")\n")
      printFilter(snv_data, filtered_snv)
    }
  } else {
    cat("No SNV file selected or invalid input.\n")
  }

  ###
  # WRITE OUT THE DATA
  ###
  wb <- createWorkbook()
  for (nm in names(sheets)) add_sheet(wb, nm, sheets[[nm]])

  stats_data <- data.frame(
    Sheet_Name        = names(sheets),
    Variant_Count     = vapply(sheets, nrow, numeric(1)),
    Unique_Gene_Count = vapply(sheets, function(d) length(all_genes(d)), numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  summary_stats <- data.frame(
    Metric = c("Total Sheets Processed", "Total Variants Across All Sheets", "Total Unique Genes Across All Sheets"),
    Value  = c(nrow(stats_data),
               sum(stats_data$Variant_Count, na.rm = TRUE),
               length(unique(unlist(lapply(sheets, all_genes))))),
    stringsAsFactors = FALSE
  )
  add_sheet(wb, "Summary_Statistics", stats_data)
  writeData(wb, sheet = "Summary_Statistics", x = summary_stats, startRow = nrow(stats_data) + 3)

  save_unique(wb, paste0(sample_id, "_Inheritance_Filter_Dominant"))
}
