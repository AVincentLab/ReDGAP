# How to use ReDGAP for filtering SNV, SVs, and CNVs

This walkthrough further explains the scripts for analyzing **single nucleotide variants (SNVs)**, **structural variants (SVs)**, and **copy number variants (CNVs)**, and later applying **inheritance filtering**. These tools are designed to facilitate variant discovery in inherited retinal dystrophies (IRDs) and other rare diseases. Included in this folder are 3 sample input files: 'Sample_SNV_Input', 'Sample_SV_Input', and 'Sample_CNV_Input'. 

Note: The ReDGAP scripts were written to filter through variants files annotated by **The Centre for Applied Genomics (TCAG)** at The Hospital for Sick Children (SickKids). These scripts can be used to filter through  other annotated files — for example, variants called in-house, through another sequencing core, or by independent analysis — either by renaming required columns in the annotated variant files, or by modifying the scripts.

---

## 1. SNV Pipeline 

The SNV pipeline filters, scores, and ranks small variants (SNVs and indels). The script is written to filter and find strong SNVs (i.e. predicted as damaging, deleterious, conserved or splicing), SNVs found on ClinVar, and SNVs that match an optional panel list of genes given. 

### Input Requirements  
Files must be in a tabular format and supported extensions are: `tsv`, `csv`, `txt`, `dsv`. 

**Required Columns**
| Column Name | General Purpose |
| :--- | :--- |
| **seq_priority** | Used to distinguish genic vs. non-genic variants and remove intergenic, UTR, upstream/downstream, and ncRNA variants before downstream filtering. |
| **dbscSNV_ADA_SCORE** | Splice prediction score from dbscSNV. Used for splice variant prioritization (ADA score ≥ 0.6). |
| **dbscSNV_RF_SCORE** | Splice prediction score from dbscSNV. Used for splice variant prioritization (Random Forest score ≥ 0.6). |
| **effect_priority** | Functional consequence annotation. Used to identify loss-of-function, frameshift, stopgain, stoploss, startloss, startgain, and nonframeshift variants. |
| **dbsnp_region** | Used to filter nonframeshift insertion/deletion variants based on the number of dbSNP entries in the region. |
Note: The script as written should be modified if the above columns are missing from the annotated variant file

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
Variants are filtered to subset for rarity.

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




## 2. SV Pipeline

## 3. CNV Pipeline

## 4. Inheritance Filtering

---
