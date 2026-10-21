# How to use ReDGAP for filtering SNV, SVs, and CNVs

This walkthrough further explains the scripts for analyzing **single nucleotide variants (SNVs)**, **structural variants (SVs)**, and **copy number variants (CNVs)**, and later applying **inheritance filtering**. These tools are designed to facilitate variant discovery in inherited retinal dystrophies (IRDs) and other rare diseases. Included in this folder are 3 sample input files: 'Sample_SNV_Input', 'Sample_SV_Input', and 'Sample_CNV_Input'. 

Note: The ReDGAP scripts were written to filter through variants files annotated by **The Centre for Applied Genomics (TCAG)** at The Hospital for Sick Children (SickKids). These scripts can be used to filter through  other annotated files — for example, variants called in-house, through another sequencing core, or by independent analysis — either by renaming required columns in the annotated variant files, or by modifying the scripts.

---

## **SNV Pipeline** 

The SNV pipeline filters, scores, and ranks small variants (SNVs and indels). The script is written to filter and find strong SNVs (i.e. predicted as damaging, deleterious, conserved or splicing), SNVs found on ClinVar, and SNVs that match an optional panel list of genes given. 

### Input Requirements  
Files must be in a tabular format and supported extensions are: `tsv`, `csv`, `txt`, `dsv`. 

**Required Columns**
| Column Name | General Purpose |
| :--- | :--- |
| **typeseq_priority** | Used to distinguish genic vs. non-genic variants and remove intergenic, UTR, upstream/downstream, and ncRNA variants before downstream filtering. |
| **dbscSNV_ADA_SCORE** | Splice prediction score from dbscSNV. Used for splice variant prioritization (ADA score ≥ 0.6). |
| **dbscSNV_RF_SCORE** | Splice prediction score from dbscSNV. Used for splice variant prioritization (Random Forest score ≥ 0.6). |
| **effect_priority** | Functional consequence annotation. Used to identify loss-of-function, frameshift, stopgain, stoploss, startloss, startgain, and nonframeshift variants. |
| **dbsnp_region** | Used to filter nonframeshift insertion/deletion variants based on the number of dbSNP entries in the region. |

*Note: The script as written should be modified if the above columns are missing from the annotated variant file*

**Optional Columns**

| Column Name | General Purpose |
| :--- | :--- |
| **gene_symbol** | Gene name(s) affected by the variant. Used for gene panel matching and prioritization. |
| **FILTER** | Variant quality status (e.g., PASS). Optional quality filtering. |
| **AD_ALT** | Alternate allele read depth. Combined with AD_REF to calculate total depth (newDP). |
| **AD_REF** | Reference allele read depth. Combined with AD_ALT to calculate total depth (newDP). |
| **spliceAI_DS_AG** | SpliceAI acceptor gain delta score. Contributes to splice impact score (sAI_delta). |
| **spliceAI_DS_AL** | SpliceAI acceptor loss delta score. Contributes to splice impact score (sAI_delta). |
| **spliceAI_DS_DG** | SpliceAI donor gain delta score. Contributes to splice impact score (sAI_delta). |
| **spliceAI_DS_DL** | SpliceAI donor loss delta score. Contributes to splice impact score (sAI_delta). |
| **promoterAI_score** | Predicted promoter regulatory impact score. Used for scoring and promoter variant output. |
| **A1000g_freq_max** | Maximum allele frequency in the 1000 Genomes dataset. Used for rare variant filtering. |
| **gnomadV4_exome_fafmax_faf95_max** | Maximum credible filtering allele frequency (gnomAD exomes). Used for rare variant filtering. |
| **gnomadV4_exome_AF_grpmax** | Highest observed population allele frequency in gnomAD exomes. Used for rare variant filtering. |
| **gnomadV4_genome_fafmax_faf95_max** | Maximum credible filtering allele frequency (gnomAD genomes). Used for rare variant filtering. |
| **gnomadV4_genome_AF_grpmax** | Highest observed population allele frequency in gnomAD genomes. Used for rare variant filtering. |
| **gnomadV4_exome_AF_remaining** | Remaining exome population allele frequency. Used for rare variant filtering. |
| **gnomadV4_genome_AF_remaining** | Remaining genome population allele frequency. Used for rare variant filtering. |
| **gnomadV4_genome_nhomalt** | Number of homozygous alternate individuals in gnomAD genomes. Used for rarity filtering. |
| **gnomadV4_exome_nhomalt** | Number of homozygous alternate individuals in gnomAD exomes. Used for rarity filtering. |
| **hsint_AF_all** | Internal cohort allele frequency. Used for rarity filtering. |
| **hsint_AF_EUR** | Internal European cohort allele frequency. Used for rarity filtering. |
| **phylopMam_avg** | Mammalian conservation score. Predictor vote. |
| **phylopVert100_avg** | Vertebrate conservation score. Predictor vote. |
| **PROVEAN_score** | In-silico deleteriousness prediction. Predictor vote. |
| **CADD_phred** | CADD pathogenicity score. Predictor vote. |
| **alpha_missense_score** | AlphaMissense pathogenicity score. Predictor vote. |
| **primateai3d_score** | PrimateAI-3D pathogenicity score. Predictor vote. |
| **MPC_score** | Missense badness/constraint score. Predictor vote. |
| **REVEL_score** | Ensemble missense pathogenicity score. Predictor vote. |
| **Clinvar_SIG** | Clinical significance from ClinVar. Used for scoring, ClinVar reports, and strong variant identification. |
| **MPO** | Mammalian Phenotype Ontology annotation. Used for phenotype relevance scoring. |
| **HPO** | Human Phenotype Ontology annotation. Used for phenotype relevance scoring. |
| **Fantom5DPI** | Indicates overlap with robust FANTOM5 regulatory/promoter elements. Used for scoring and output lists. |


**Example SNV Input (excerpt):**  

```
gene_symbol    FILTER    AD_REF    AD_ALT    1000g_all    typeseq_priority    spliceAI_DS_AG    spliceAI_DS_AL    spliceAI_DS_DG    spliceAI_DS_DL    Clinvar_SIG    Clinvar_ReviewStatus    effect_priority    sift_score    polyphen_score    ma_score    spx_dpsi    phylopMam_avg    phylopVert100_avg    PROVEAN_score    CADD_phred    REVEL_score    MPO    HPO
ABCA4          PASS      35        12        0.0001       missense_variant    0.02              0.00              0.00              0.01              Pathogenic     reviewed_by_expert      missense           0.01          0.95             1.2        -0.05        1.75            2.41                -2.3             25.6          0.86           MP:0005390    HP:0000548
```

### Running the Script
1. Open `ReDGAP_SNV_GRCh38_TCAG_pipeline28_5.R` in an **interactive** R session (e.g. RStudio) and run it. Missing packages are installed automatically.
2. Select the patient variant file when the file browser opens.
3. Answer the prompts (type `yes` or `no`):
   - Keep only variants with `FILTER == PASS`?
   - Keep only variants with read depth (`newDP`) between 5 and 120?
   - Use a panel gene list? If yes, you will be asked to select the panel `txt` file.
4. The script prints how many variants each step removed or kept, then saves the workbook in your working directory.

### Output  
#### Output Format
For all variants in the output there are additional columns to the ones listed above: 

```
gene_symbol    FILTER    AD_REF    AD_ALT    ...    sAI_delta    newDP    VariantScore
ABCA4          PASS      35        12        ...    0.42         47       6.5
```

- `sAI_delta` — combined SpliceAI prediction score
- `newDP` — normalized read depth across reference and alternate alleles  
- `VariantScore` — cumulative variant score used for ranking

The output file is named `<input file name>_Pipeline_SNVs.xlsx` (a number is added, e.g. `_1_`, if the file already exists).

#### Output Workbooks
The SNV script will output an excel workbook with three to five sheets. There are two *optional* sheets added when the user decides to provide a panel gene list.

| Sheet Name | Description |
| :--- | :--- |
| **SNV_Filter_History** | A summary sheet of the steps in the filtering and the number of variants removed or kept|
| **strong_SNV** | The variants that are predicted as damaging, deleterious, conserved or splicing. A group list of variants that are included at least one of the following categories: Pathogenic or likely pathogenic in ClinVar, have some effect (i.e. "stopgain", "frameshift", etc.), pass at least two in-silico predictors, have a FANTOM5 score = 2, have a spiceAI score >= 0.1 OR dbscSNV (ADA or RF) >= 0.6, OR have a promoterAI score > 0.1 |
| **ClinVar** | The variants that have a ClinVar designation as pathogenic, likely pathogenic, or uncertain significance |
| **Panel** | *OPTIONAL* The variants that match the the gene panel provided by user |
| **Panel_Gene** | *OPTIONAL* List of the panel genes used |

### Panel Gene List
Users can decide to add a list of panel genes which are genes of interest that you want to keep variants for. This list should be a `txt` file with a column header of `gene_symbol` followed by a gene per line. 

### Allele Frequency Filtering
Variants are filtered to subset for rarity. A variant with a missing value in a column is treated as rare and is kept.

| Column | Database / Resource | Filter Threshold |
| :--- | :--- | :--- |
| **A1000g_freq_max** | 1000 Genomes Project maximum allele frequency | Keep variants with MAF ≤ 0.5% (0.005) |
| **gnomadV4_exome_fafmax_faf95_max** | gnomAD v4 Exomes FAF95 maximum filtering allele frequency | Keep variants with FAF ≤ 0.5% (0.005) |
| **gnomadV4_exome_AF_grpmax** | gnomAD v4 Exomes maximum population-specific allele frequency (grpmax) | Keep variants with AF ≤ 0.5% (0.005) |
| **gnomadV4_genome_fafmax_faf95_max** | gnomAD v4 Genomes FAF95 maximum filtering allele frequency | Keep variants with FAF ≤ 0.5% (0.005) |
| **gnomadV4_genome_AF_grpmax** | gnomAD v4 Genomes maximum population-specific allele frequency (grpmax) | Keep variants with AF ≤ 0.5% (0.005) |
| **hsint_AF_all** | SickKids Internal Database (all populations) allele frequency | Keep variants with AF ≤ 1.0% (0.01) |
| **hsint_AF_EUR** | SickKids Internal Database (European subset) allele frequency | Keep variants with AF ≤ 1.0% (0.01) |
| **gnomadV4_exome_AF_remaining** | gnomAD v4 Exomes allele frequency in remaining populations | Keep variants with AF ≤ 0.5% (0.005) |
| **gnomadV4_genome_AF_remaining** | gnomAD v4 Genomes allele frequency in remaining populations | Keep variants with AF ≤ 0.5% (0.005) |
| **gnomadV4_genome_nhomalt** | gnomAD v4 Genomes homozygous alternate carrier count | Keep variants with ≤ 7 homozygous individuals |
| **gnomadV4_exome_nhomalt** | gnomAD v4 Exomes homozygous alternate carrier count | Keep variants with ≤ 7 homozygous individuals |


### Scoring Variants
Variants are scored (**VScore**) to help prioritise which ones to explore first. 

<img width="1350" height="727" alt="SNV_VScore_Chart" src="https://github.com/user-attachments/assets/c5f5dab7-a7bc-4e95-b6d7-7b8da066b5f3" />


---

## **SV & CNV Pipelines**
The scripts for filtering structural variants (SV) and copy number variants (CNV) are very similar.

### Input Requirements
Files must be in a tabular format and supported extensions are: `tsv`, `csv`, `txt`, `dsv`. 

**Columns**
The scripts checks for column existence before running them, but using all columns leads to best filtering and scoring.

| Column Name | General Purpose |
| :--- | :--- |
| **gene_symbol** | Maps SVs to genes. Required for panel-gene filtering, creation of Genic SV lists, and panel-specific outputs. |
| **cds_symbol** | Identifies SVs overlapping coding sequences (CDS) and contributes to variant scoring and the High Impact SV list. |
| **ClinVar_ClinSign** | Provides ClinVar clinical significance annotations used for pathogenic/VUS scoring and clinical database outputs. |
| **ClinGen_SV_pctOvlp_max** | Percentage overlap with ClinGen pathogenic/curated SVs. Used for scoring and clinical database prioritization. |
| **decipher_SV_pctOvlp_max** | Percentage overlap with DECIPHER SVs. Used for scoring and clinical database prioritization. |
| **numberOfGeneSymbols** | Number of genes affected by the SV. Used to increase variant priority scores. |
| **nearestLeftExonDistance** | Distance to the nearest left exon boundary. Used for splice-site proximity scoring and splice-region identification. |
| **nearestRightExonDistance** | Distance to the nearest right exon boundary. Used for splice-site proximity scoring and splice-region identification. |
| **Fantom5DPI** | Indicates overlap with FANTOM5 promoter regions. Used for promoter scoring and regulatory SV prioritization. |
| **gnomAD_AF_MAX_max90pctRecOvlp** | Maximum population allele frequency of a matching structural variant in gnomAD (calculated using a strict 90% reciprocal overlap threshold). Used to remove common SVs. |
| **otgMantaPctFreq_90pctRecOvlp** | Frequency of matching SVs in the 1000 Genomes+ Manta callset. Used for rarity filtering. |
| **otgDellyPctFreq_90pctRecOvlp** | Frequency of matching SVs in the 1000 Genomes+ Delly callset. Used for rarity filtering. |
| **otgDgnPctFreq_90pctRecOvlp** | Frequency of matching SVs in the 1000 Genomes+ DRAGEN callset. Used for rarity filtering. |
| **all_DgnPctFreq_90pctRecOvlp** | Frequency of matching SVs in the internal DRAGEN structural variant database. Used for rarity filtering. |
| **nstd186_pctOvlp** | Overlap with NCBI common SV set (nstd186). Used to remove common variants. |
| **DGV_commonPct** | Percentage overlap with common structural variants in the Database of Genomic Variants (DGV). Used to remove common variants. |

### Running the Scripts
1. Open `ReDGAP_SV_GRCh38_TCAG_pipeline1_7_3.R` or `ReDGAP_CNV_GRCh38_TCAG_pipeline1_7_3.R` in an **interactive** R session (e.g. RStudio) and run it.
2. Select the patient variant file when the file browser opens.
3. Answer the prompts (type `yes` or `no`):
   - Use a panel gene list? If yes, you will be asked to select the panel `txt` file.
   - Remove variants overlapping polymorphic regions (`dirtyRegion_pctOvlp`) by ≥ 70%?
5. The output is saved as `<input file name>_High_impact_SV.xlsx` or `<input file name>_High_impact_CNVs.xlsx`.

### Output Workbook
The SV and CNV scripts will output an excel workbook with five to seven sheets. There are two *optional* sheets added when the user decides to provide a panel gene list.

| Sheet Name | Description |
| :--- | :--- |
| **X_Analysis_Pipeline** | A summary sheet of the steps in the filtering and the number of variants removed or kept |
| **Strong_X** | The variants that are predicted to have a high impact. The list of variants that are included in at least one of the following categories: Pathogenic or likely pathogenic in ClinVar, are present in the coding region of genes (`cds_symbol`), are near a canonical splice site (<20bp), or have a non-zero FANTOM5 |
| **rare_X** | The variants that have passed through the frequency filtering |
| **Genic_X** | The variants that include genes in `gene_symbol` |
| **clin_db_X** | The variants that have a ClinVar designation as pathogenic, likely pathogenic, or uncertain significance OR overlap with variants from either ClinGen or DECIPHER databases by 30% |
| **Panel_X** | *OPTIONAL* The variants that match the the gene panel provided by user |
| **Panel_list** | *OPTIONAL* List of the panel genes used |

*NOTE: The "X" in the Sheet Name of the workbook will be either "SV" or "CNV" depending on the script run.*

### Panel Gene List
Users can decide to add a list of panel genes which are genes of interest that you want to keep variants for. This list should be a `txt` file with a column header of `gene_symbol` followed by a gene per line. 

### Allele Frequency Filtering
Variants are filtered to subset for rarity. A variant with a missing value in a column is treated as rare and is kept. The SV script uses the 90% reciprocal overlap columns listed below; the CNV script uses the equivalent 50% reciprocal overlap columns (e.g. `gnomAD_AF_MAX_max50pctRecOvlp`, `otgErdsPctFreq_50pctRecOvlp`, `otgCnvnPctFreq_50pctRecOvlp`, `otgDgnPctFreq_50pctRecOvlp`, `all_DgnPctFreq_50pctRecOvlp`) with the same thresholds.

| Filter | Database / Resource | Filter Threshold |
| :--- | :--- | :--- |
| **gnomAD_AF_MAX_max90pctRecOvlp** | Maximum population allele frequency of a matching structural variant in gnomAD v4.1 (calculated using a strict 90% reciprocal overlap threshold) | Keep SVs with MAF ≤ 0.5% (0.005) |
| **otgMantaPctFreq_90pctRecOvlp** | atching structural variant in the 1000 Genomes Project callset generated with Manta (90% reciprocal overlap threshold) | Keep SVs with MAF ≤ 0.5% (0.005) |
| **otgDellyPctFreq_90pctRecOvlp** | Matching structural variant in the 1000 Genomes Project callset generated with Delly (90% reciprocal overlap threshold) | Keep SVs with MAF ≤ 0.5% (0.005) |
| **otgDgnPctFreq_90pctRecOvlp** | Allele frequency of a matching structural variant in the 1000 Genomes Project callset generated with DRAGEN (90% reciprocal overlap threshold) | Keep SVs with AF ≤ 0.5% (0.005) |
| **all_DgnPctFreq_90pctRecOvlp** | Allele frequency of a matching structural variant in an internal database of ~10,600 unrelated samples called with DRAGEN (90% reciprocal overlap threshold) | Keep SVs with AF ≤ 0.5% (0.005) |
| **nstd186_pctOvlp** | Minor allele frequency of a matching structural variant in the NCBI curated common SV dataset nstd186 | Keep SVs with MAF ≤ 1.0% (0.01) |
| **DGV_commonPct** | Reciprocal overlap with common structural variants (>1% frequency) from the Database of Genomic Variants (DGV) | Keep SVs with <80% overlap with common DGV variants |

### Scoring Variants
Variants are scored (**VScore**) to help prioritise which ones to explore first. 
<img width="1350" height="576" alt="SV_VScore_Chart" src="https://github.com/user-attachments/assets/344df0d0-cd86-4afe-bdee-341a3c8cb110" />


---
## **Inheritance Filtering**
The inheritance filtering script (`ReDGAP_Inheritance_Filtering.R`) is run **after** the SNV, SV and/or CNV scripts. It takes one, or any combination, of the pipeline outputs and looks for genes that fit a **recessive** or **dominant** model. It matches variants **by gene name**, so it can be used for a single patient or to find genes shared between affected family members.

### Input Requirements
- Files can be `tsv`, `csv`, `txt`, `dsv`, `xlsx` or `xls`. For Excel files you choose which sheets to load (e.g. `strong_SNV` and `Strong_SV` from the pipeline workbooks).
- Every sheet needs a column called **`gene_symbol`**. Variants that overlap several genes should separate gene names with `|` (`;` and `,` are also understood). Sheets without this column load with a warning and can not be used for gene matching.
- For **recessive** mode, homozygous variants are recognised from a `Zygosity` column (value `hom-alt`) or, if there is no `Zygosity` column, from a genotype column whose name contains `GT` (value `1/1` or `1|1`).
- A `VScore` column (created by the pipeline scripts) is optional, but if present it is used to score genes.
- Everything is read as text, so sheets from different files can be compared.

### Running the Script
The script asks questions as it runs, so it must be run in an **interactive** R session (e.g. RStudio).

1. Write the **Sample ID** (used in the output file name).
2. Select a file. For an Excel file, enter the number of each sheet you want to load (separated by spaces or commas) and type `q` to finish.
3. Give each loaded sheet a short label for its variant type (e.g. `SNV`, `SV`, `CNV`). These labels are used in the output sheet names.
4. Choose whether to add more files and sheets. Repeat until finished.
5. Choose the inheritance mode: `[1] Recessive` or `[2] Dominant`.

### Recessive Mode
**Step 1: Recessive candidates within each loaded sheet.** A variant is kept if:
- its gene has **more than one variant** in that sheet (possible compound heterozygote), or
- it is **homozygous** (see Input Requirements).

Kept variants get a `recessive_reason` column (`duplicate_gene`, `hom_alt_zygosity` or `homozygous_GT`). If the sheet has a `VScore` column, two gene-level scores are added: `GeneScore_AvgAll` (mean VScore of the gene's candidate variants) and `GeneScore_AvgTop2` (mean of the gene's two highest VScores). A variant in several genes takes the highest score of its genes.

**Step 2: Genes shared between sheets.** Every pair of loaded sheets is compared (e.g. `SNV` vs `SV`, `SNV` vs `CNV`, `SV` vs `CNV`). For each pair with genes in common, the genes are listed and **all rows** of the original sheets that include a shared gene are saved. Check these rows yourself for the specific variant matches.

**Step 3: Statistics.** Summary sheets are added to the workbook.

The output is saved as `<SampleID>_Inheritance_Filter_Recessive.xlsx` (a number is added if the file already exists).

#### Recessive Output Workbook
Example below uses sheets labelled `SNV`, `SV` and `CNV`. Comparison sheets are abbreviated as `A_vs_B`, where `A` and `B` stand for any two loaded sheet labels; one set of three sheets is made for each pair that has shared genes. Excel sheet names are limited to 31 characters (long names are shortened and repeated names get `_2`, `_3`...).

| Sheet Name | Created When | Description |
| :--- | :--- | :--- |
| **SNV_recessive** | Each loaded sheet (`X_recessive`) has at least one candidate | Variants from the `SNV` sheet in genes with >1 variant or that are homozygous, with `recessive_reason` (and `GeneScore_AvgAll` / `GeneScore_AvgTop2` if VScore exists). |
| **SV_recessive** | As above | Same as above for the `SV` sheet. |
| **CNV_recessive** | As above | Same as above for the `CNV` sheet. |
| **Genes_A_vs_B** | Pair `A`, `B` share at least one gene | The list of genes (`gene_symbol`) shared between sheet `A` and sheet `B`. |
| **A_A_vs_B** | As above | All rows of sheet `A` with shared genes, with a `SourceSheet` column. |
| **B_A_vs_B** | As above | All rows of sheet `B` with shared genes, with a `SourceSheet` column. |
| **Analysis_Statistics** | Always | One row per original sheet, recessive sheet and comparison: `Sheet_Name`, `Sheet_Type` (Original / Recessive / Comparison), `Total_Variants`, `Biallelic_Recessive_Count`, `Genes_Involved`, `Gene_Count_Total`, `Gene_Count_Retained`, `Notes`. |
| **Summary_Statistics** | Always | Totals: original sheets, recessive sheets created, comparison sheets created, unique genes analyzed, unique genes retained, biallelic/recessive variants found. |
| **All_Unique_Genes** | At least one comparison found genes | Every shared gene and the comparison it came from (`Comparison`, `gene_symbol`). |

Example with three loaded sheets (`SNV`, `SV`, `CNV`) where all pairs share genes:

| Sheet Name |
| :--- |
| SNV_recessive, SV_recessive, CNV_recessive |
| Genes_SNV_vs_SV, SNV_SNV_vs_SV, SV_SNV_vs_SV |
| Genes_SNV_vs_CNV, SNV_SNV_vs_CNV, CNV_SNV_vs_CNV |
| Genes_SV_vs_CNV, SV_SV_vs_CNV, CNV_SV_vs_CNV |
| Analysis_Statistics, Summary_Statistics, All_Unique_Genes |

### Dominant Mode
1. The loaded sheets are listed. Enter the number of the **SNV sheet**, or `0` if there is none.
2. The selected SNV sheet is filtered for very rare variants: **≤ 0.03% (0.0003)** in `A1000g_freq_max` and `freq_max` (gnomAD). Variants with a missing frequency are kept. If neither column exists in the sheet, no frequency filter is applied.
3. All loaded sheets (SV and CNV sheets are not filtered) are saved to one workbook.

The output is saved as `<SampleID>_Inheritance_Filter_Dominant.xlsx`.

#### Dominant Output Workbook

| Sheet Name | Description |
| :--- | :--- |
| **Each loaded sheet** (e.g. `SNV`, `SV`, `CNV`) | The sheet as loaded. The chosen SNV sheet has been filtered for allele frequency ≤ 0.03%; other sheets are unchanged. |
| **Summary_Statistics** | A table with `Sheet_Name`, `Variant_Count` and `Unique_Gene_Count` for each sheet, followed (a few rows below) by totals: total sheets processed, total variants across all sheets, and total unique genes across all sheets. |

---
## **Walkthrough**
This example follows one patient with SNV, SV and CNV files.

### 1. Prepare your files
- Annotated variant files from TCAG for SNVs (pipeline 28.5), and SVs and CNVs (pipeline 1.7.3), saved as `tsv`, `csv`, `txt` or `dsv`.
- *Optional:* a panel gene list as a `txt` file with the header `gene_symbol` and one gene per line.
- R (with RStudio recommended). Packages are installed automatically the first time a script is run.

### 2. Run the SNV script
1. Open and run `ReDGAP_SNV_GRCh38_TCAG_pipeline28_5.R`.
2. Select the SNV file, then answer the PASS filter, read depth and panel questions.
3. Output: `<name>_Pipeline_SNVs.xlsx` with `SNV_Filter_History`, `strong_SNV`, `ClinVar` (and `Panel`, `Panel_Gene_lt`).

### 3. Run the SV script
1. Open and run `ReDGAP_SV_GRCh38_TCAG_pipeline1_7_3.R`.
2. Select the SV file, then answer the panel and polymorphic region questions.
3. Output: `<name>_High_impact_SV.xlsx` with `SV_Analysis_Pipeline`, `Strong_SV`, `rare_SV`, `Genic_SV`, `clin_db_SV` (and `Panel_SV`, `Panel_list`).

### 4. Run the CNV script
1. Open and run `ReDGAP_CNV_GRCh38_TCAG_pipeline1_7_3.R`.
2. Select the CNV file, then answer the same questions as for the SV script.
3. Output: `<name>_High_impact_CNVs.xlsx` with `CNV_Analysis_Pipeline`, `Strong_CNV`, `rare_CNV`, `Genic_CNV`, `clin_db_CNV` (and `Panel_CNV`, `Panel_list`).

You do not need to run all three. Run the ones you have data for and use any combination of the outputs in the next step.

### 5. Run the inheritance filter
1. Open and run `ReDGAP_Inheritance_Filtering.R`.
2. Enter the Sample ID.
3. Load the sheets you want to compare. Choose the output workbooks one at a time, select the sheets (e.g. `strong_SNV`, `Strong_SV`, `Strong_CNV`), and label each one (`SNV`, `SV`, `CNV`). A sheet from another family member can be loaded the same way to look for shared genes.
4. Choose **Recessive** or **Dominant**.
5. Open `<SampleID>_Inheritance_Filter_Recessive.xlsx` or `<SampleID>_Inheritance_Filter_Dominant.xlsx`.

### 6. Interpret the results
- **Recessive:** look at the `X_recessive` sheets for genes with multiple variants or homozygous variants (see `recessive_reason`), ranked by `GeneScore_AvgTop2` / `GeneScore_AvgAll`. Genes appearing in the `Genes_A_vs_B` sheets are hit by more than one variant type (e.g. an SNV on one allele and a CNV on the other). The comparison sheets contain every row in a shared gene, so confirm the specific variants and their phase/segregation yourself.
- **Dominant:** look at the rare SNVs together with the SV and CNV sheets, and review each for segregation with the phenotype.

---
