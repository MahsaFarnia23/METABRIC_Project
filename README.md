
# 🔬 Survival Analysis of Breast Cancer Patients Based on Gene Mutation Status (METABRIC Dataset)

This project conducts a **survival analysis** using the **METABRIC RNA + Mutation dataset**, focusing on how mutations in key genes affect patient outcomes. It includes **descriptive statistics**, **data visualization**, **mutation profiling**, and **Kaplan-Meier survival analysis** with **log-rank testing**.


## Dataset
* **Name**: `METABRIC_RNA_Mutation.csv`
* **Source**: [METABRIC breast cancer dataset](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)
* **Content**:
- 31 clinical metadata
- mRNA z-score expression levels for 331 genes
- Gene mutation status for 175 genes
- Survival outcome data for over 1,900 breast cancer patient


## 📂 Dataset

* **Name**: `METABRIC_RNA_Mutation.csv`  
* **Source**: [METABRIC breast cancer dataset](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric)  
* **Content**:  
  - 🧾 31 clinical metadata features  
  - 📈 mRNA z-score expression levels for 331 genes  
  - 🧬 Gene mutation status for 175 genes  
  - 🕒 Survival outcome data for over 1,900 breast cancer patients  

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



---




## ⚙️ Workflow Summary

### 1️⃣ Data Loading and Inspection

* Load and inspect dataset structure
* Detect and report missing values
* Explore unique values for each categorical column

### 2️⃣ Descriptive Statistics

* Numerical summary using `.describe()`
* Categorical summary using `.describe(include='object')`

---

## 📊 Exploratory Visualizations

To understand the clinical landscape, we visualized important patient features:

### 🎂 Age at Diagnosis

A histogram showing the distribution of age among patients.

```python
sns.histplot(df['age_at_diagnosis'], bins=30, kde=True)
```

### 🏥 Type of Breast Surgery

A countplot showing frequencies of surgery types.

```python
sns.countplot(x='type_of_breast_surgery', data=df)
```

### 🧬 Cancer Type Distribution

A countplot for the distribution of different cancer types.

```python
sns.countplot(x='cancer_type', data=df)
```

These visualizations provide a clinical overview and help contextualize the mutation analysis.

---

## 🧬 Mutation-Based Survival Analysis

### 🔎 Mutation Processing

* Mutation columns ending in `_mut` are extracted.
* For each gene:

  * A binary column is created (`mutated` vs `wildtype`).
  * Patients are grouped accordingly.

### 📉 Survival Modeling

* **Kaplan-Meier Estimator** is used to model survival curves for each gene.
* **Log-rank tests** are performed to compare survival distributions.
* Genes with significant differences are identified.

### ✅ Top Genes Output

Top 10 genes with the lowest **p-values** are reported as potentially significant in influencing survival.

### 📈 Kaplan-Meier Plots

The survival curves for the top 3 genes are plotted to visually demonstrate the impact of mutations.

```python
kmf.fit(group[duration_col], group[event_col], label=label_name)
kmf.plot_survival_function()
```

---

## 🔍 Notable Mutations: Example Interpretation

> Based on output, interpretations might include:

* **TP53**: Frequently mutated; associated with aggressive tumors and **poorer prognosis**.
* **PIK3CA**: Involved in PI3K pathway; sometimes linked to **better prognosis** in certain subtypes.
* **GATA3**: Regulates differentiation; mutations may affect **luminal subtype** development.

---

## 📦 Dependencies

```bash
pip install pandas numpy matplotlib seaborn lifelines
```

---

## ▶️ How to Run

1. Save the dataset as `METABRIC_RNA_Mutation.csv` in your project directory.
2. Run the Python script:

```bash
python your_script_name.py
```

---

## 🧠 Future Directions

* Apply **multivariate Cox regression** models
* Visualize gene-gene interactions or mutational burden
* Correlate mutations with PAM50 or other molecular subtypes

---

## 📜 License

This project is for academic and non-commercial research purposes.

---

Would you like me to export this as a `.md` file ready for upload to your GitHub repo?


