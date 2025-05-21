
# METABRIC Breast Cancer Analysis

A machine learning project analyzing the METABRIC breast cancer dataset, focusing on:
- Gene mutations
- Clinical features
- Survival analysis
- Subtype classification

---

## Dataset

This project uses the [METABRIC breast cancer dataset](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric), which includes:

* 31 clinical attributes
* mRNA z-score expression levels for 331 genes
* Mutation status for 175 genes
* Data from 1,904 breast cancer patients

---


## 📂 **Project Structure**
```
METABRIC_Project/
│
├── data/                           # Dataset files
│   └── METABRIC_RNA_Mutation.csv   # Main dataset
│
├── scripts/                        # Python scripts for analysis
│   └── analysis.py                
│
├── results/                        # Model outputs and visualizations
│   └── plots/                      # Plots generated during analysis
│
├── images/                         # Plots for documentation
│   ├── age_distribution.png
│   ├── tumor_stage_distribution.png
│   ├── mutation_heatmap.png
│   └── survival_distribution.png
│
└── README.md                       # Project documentation
```
## 📊 Visualizations

- **Age at Diagnosis Distribution**
  ![Age at Diagnosis](images/age_distribution.png)

- **Type of Breast Surgery**
  ![Type of Breast Surgery](images/type_of_breast_surgery.png)

- **Cancer Type Distribution**
  ![Cancer Type Distribution](images/cancer_type_distribution.png)

