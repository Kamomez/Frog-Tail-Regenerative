# Frog Tail Regeneration - ROC Identification

Single-cell RNA-seq analysis for identifying the Regenerative Organizing Cell (ROC) in the Xenopus tail.

Reference Papers:

https://www.science.org/doi/full/10.1126/science.aav9996

https://www.nature.com/articles/s41587-025-02694-w

https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aav9996&file=aav9996_aztekin_sm.pdf

---

## Start
Dataset cleaned_processed_frogtail.h5ad available to download from the journal website

recommend to pip install requirement.txt first

#### run local file step by step (Recommended)
local code/01_data_preprocessing.py

local code/02_data_denoising.py

local code/03_batch_correction_clustering.py

local code/04_biomarker_analysis.py

#### run local ipynb
local code/frog_tail_analysis_local.ipynb

#### run colab ipynb
colab code/frog_tail_analysis_colab.ipynb


---
## Overview

**Dataset**: 13,199 cells 31,535 genes, 4 time points (days 0, 1, 2, 4)  
**Key Finding**: Cluster 2 identified as ROC with **99.8% early enrichment** (days 0-1)  

---

## Results

### Top 10 Marker Genes (Logistic Regression)
| Rank | Gene | Importance Score | Notes |
|------|------|-----------------|-------|
| 1 | Xelaev18032448m.g | 0.574 | Wound healing response |
| 2 | akr1c2.L | 0.530 | Oxidative stress regulation |
| 3 | Xelaev18014353m.g | 0.449 | Cell signaling |
| 4 | nos2.L | 0.448 | Immune/inflammatory response |
| 5 | loc100492852.L | 0.365 | Metabolic regulation |
| 6 | gstp1.L | 0.363 | Detoxification enzyme |
| 7 | pc.L | 0.353 | Gluconeogenesis |
| 8 | Xelaev18022490m.g | 0.347 | Cell proliferation |
| 9 | sst.L | 0.326 | Growth hormone inhibitor |
| 10 | loc100135371.L | 0.308 | Developmental regulation |

### Top 10 Marker Genes (Wilcoxon Test)
| Rank | Gene | P-value | Log2 Fold Change |
|------|------|---------|------------------|
| 1 | itln1.L | 2.49e-261 | 2.69 |
| 2 | loc100127564.L | 1.76e-255 | 2.27 |
| 3 | Xelaev18019895m.g | 2.39e-238 | 2.01 |
| 4 | Xelaev18034081m.g | 1.37e-234 | 2.11 |
| 5 | Xelaev18026753m.g | 1.77e-234 | 2.15 |
| 6 | akr1c2.L | 4.46e-219 | 2.61 |
| 7 | loc100492852.L | 1.23e-218 | 2.50 |
| 8 | fbp1.L | 2.08e-218 | 2.18 |
| 9 | pc.L | 6.78e-216 | 2.40 |
| 10 | pck1.L | 2.98e-211 | 1.93 |

### High-Confidence Markers (Validated by Both Methods)
8 genes identified by both Logistic Regression AND Wilcoxon test:
- `akr1c2.L` - Oxidative stress response (top in both methods)
- `pc.L` - Gluconeogenesis enzyme
- `gstp1.L` - Glutathione S-transferase
- `loc100492852.L` - Metabolic regulation
- `krt.S` - Keratin structural protein
- `grhl3.S` - Grainyhead-like transcription factor
- `Xelaev18022528m.g` - Novel gene
- `Xelaev18020205m.g` - Novel gene

---

## Clustering Quality

| Method | Denoising | Batch Correction | Silhouette | ARI |
|--------|-----------|------------------|------------|-----|
| **Best Quality** | MAGIC | Scanorama | **0.485** | 0.486 |
| Baseline | - | Harmony | 0.451 | 0.512 |
| Standard | MAGIC | - | 0.461 | 0.546 |
| Baseline | - | BBKNN | 0.451 | 0.512 |
| Standard | ALRA | - | 0.451 | 0.480 |
| **Fastest** | - | Linear Regression | 0.419 | 0.523 |
