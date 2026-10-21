# ReDGAP
**Retinal Degeneration Genome Analysis Pipeline**

ReDGAP is a phenotype-guided, semi-automated genome analysis pipeline written in R for the prioritization of single nucleotide variants (SNVs), structural variants (SVs), and copy number variants (CNVs). The pipeline integrates variant filtering, scoring, inheritance modeling, and phenotype-driven prioritization to accelerate rare disease variant discovery, with a particular focus on inherited retinal dystrophies (IRDs).

ReDGAP supports the joint analysis of multiple variant classes and produces concise, review-ready outputs that facilitate downstream interpretation and candidate gene discovery. The framework was developed and validated using whole-genome sequencing data from IRD cohorts and achieved a 100% re-identification rate in previously solved cases while resolving 80% of a novel unsolved cohort.

---

## Publication

If you use or adapt ReDGAP in your research, please cite:

> Ahmed L, Tavares E, Li JM, Ahmed K, Mehta M, Eileen C, Ah-Sen G, Osman R, Green-Sanderson K, Dvaladze A, et al. *A novel phenotype-guided genome analysis pipeline for variant discovery*. **npj Genomic Medicine** (2026). DOI: 10.1038/s41525-026-00557-0.

**Paper:** https://doi.org/10.1038/s41525-026-00557-0 

---

## Key Features

- SNV, SV, and CNV prioritization workflows
- Interactive R-based analysis with minimal setup
- Phenotype-guided panel gene analysis (PGA)
- Whole genome analysis (WGA) mode
- Flexible and customizable filtering criteria
- Evidence-weighted variant scoring system
- Recessive and dominant inheritance filtering
- Joint analysis across variant classes
- Excel workbook outputs for rapid review
- Adaptable to non-IRD rare disease studies 

---

## Pipeline Overview

ReDGAP consists of four major workflows:

### 1. SNV Analysis
Filters, scores, and ranks SNVs and small indels using:

- Population allele frequencies
- Variant consequence annotations
- Splice prediction tools
- Conservation metrics
- Multiple pathogenicity predictors
- ClinVar annotations
- Phenotype information (HPO/MPO)
- Optional gene panel filtering 

### 2. Structural Variant (SV) Analysis

Identifies and prioritizes:

- Deletions
- Duplications
- Insertions
- Inversions
- Other structural rearrangements

SVs are filtered against multiple population databases and prioritized using gene overlap, ClinVar, ClinGen, DECIPHER, splice proximity, and regulatory annotations. 

### 3. Copy Number Variant (CNV) Analysis

Performs filtering and prioritization of CNVs using a workflow analogous to the SV pipeline while incorporating frequency and functional impact information.

### 4. Inheritance Filtering

Combines outputs from SNV, SV, and CNV analysis to identify candidate genes under:

- **Autosomal recessive inheritance**
  - Homozygous variants
  - Potential compound heterozygous events
  - Cross-variant-class gene matches

- **Autosomal dominant inheritance**
  - Rare heterozygous candidate variants
  - Stringent frequency filtering

This workflow enables identification of genes impacted by multiple variant classes, such as an SNV on one allele and an SV or CNV on the other.

---

## Example Workflow

Detailed instructions, required columns, filtering thresholds, scoring rules, output descriptions, and a complete end-to-end example are provided in [HowToUse.md]([Sample Files and Walkthrough/HowToUse.md](https://github.com/AVincentLab/ReDGAP/blob/main/Sample%20Files%20and%20Walkthrough/HowToUse.md)). A brief overview of the [workflow](https://github.com/AVincentLab/ReDGAP/blob/main/Sample%20Files%20and%20Walkthrough/HowToUse.md#walkthrough) is also included.

---
## Repository Structure

ReDGAP/
├── scripts/
│   ├── ReDGAP_SNV_GRCh38_TCAG_pipeline28_5.R
│   ├── ReDGAP_SV_GRCh38_TCAG_pipeline1_7_3.R
│   ├── ReDGAP_CNV_GRCh38_TCAG_pipeline1_7_3.R
│   └── ReDGAP_Inheritance_Filtering.R  
├── Sample Files and Walkthrough/
│   ├── HowToUse.md
│   ├── Sample_SNV_Input.tsv
│   ├── Sample_SV_Input.*
│   └── Sample_CNV_Input.* 
├── LICENSE
└── README.md

---
## License

Please see the repository license for usage terms.

---
## Contact

For questions or contributions, please feel free to get in touch with the corresponding author.
