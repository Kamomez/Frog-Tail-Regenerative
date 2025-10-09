# Frog Tail Regeneration - ROC Identification

Single-cell RNA-seq analysis for identifying the Regenerative Organizing Cell (ROC) in the Xenopus tail.

Reference Papers:
https://www.science.org/doi/full/10.1126/science.aav9996
https://www.nature.com/articles/s41587-025-02694-w
https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aav9996&file=aav9996_aztekin_sm.pdf

---

## Quick Start
Dataset cleaned_processed_frogtail.h5ad available to download from the journal website
```bash
# run local file step by step (Recommended)
python "local code/01_data_preprocessing.py"
python "local code/02_data_denoising.py"
python "local code/03_batch_correction_clustering.py"
python "local code/04_biomarker_analysis.py"
```

---
## Overview

**Dataset**: 13,199 cells  31,535 genes, 4 time points (days 0, 1, 2, 4)  
**Key Finding**: Cluster 2 identified as ROC with **99.8% early enrichment** (days 0-1)  

### Top Markers
| Gene | Score | Function |
|------|-------|----------|
| Xelaev18032448m.g | 0.574 | Wound healing |
| akr1c2.L | 0.530 | Oxidative stress |
| nos2.L | 0.448 | Immune response |

### Performance
Best silhouette score and lowest runtime
| Method | Silhouette | Runtime |
|--------|-----------|---------|
| MAGIC+Scanorama | **0.485** | 165s |
| Baseline+Harmony | 0.451 | **27s** |

