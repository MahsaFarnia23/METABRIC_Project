import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
from lifelines import CoxPHFitter
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report, RocCurveDisplay
from sklearn.model_selection import train_test_split
import shap



warnings.filterwarnings("ignore")

# Set styles for plots
sns.set(style='whitegrid')

# Load the data
df = pd.read_csv('data/METABRIC_RNA_Mutation.csv')

# Quick check
print("Shape of the data:", df.shape)

# Check missing values
missing = df.isnull().sum()
missing = missing[missing > 0].sort_values(ascending=False)
if not missing.empty:
    print("Missing values found in the following columns:")
    print(missing)
else:
    print("No missing values!")

# Summary stats
print(df.describe())
print(df.dtypes.value_counts())

# Unique values in object columns
object_cols = df.select_dtypes(include='object').columns
for col in object_cols:
    print(f"{col}: {df[col].nunique()} unique values")

# Visualization: Age at Diagnosis
plt.figure(figsize=(6, 4))
sns.histplot(df['age_at_diagnosis'].dropna(), bins=30, kde=True, color='skyblue')
plt.title('Age at Diagnosis Distribution')
plt.xlabel('Age at Diagnosis')
plt.ylabel('Frequency')
plt.tight_layout()
plt.show()

# Visualization: Type of Breast Surgery
plt.figure(figsize=(6, 4))
sns.countplot(data=df, x='type_of_breast_surgery', palette='Set2')
plt.title('Type of Breast Surgery Distribution')
plt.xlabel('Surgery Type')
plt.ylabel('Count')
plt.tight_layout()
plt.show()

# Visualization: Cancer Type
plt.figure(figsize=(6, 4))
sns.countplot(data=df, x='cancer_type', palette='Set1')
plt.title('Cancer Type Distribution')
plt.xlabel('Cancer Type')
plt.ylabel('Count')
plt.tight_layout()
plt.show()


# Identify survival columns
duration_col = 'overall_survival_months' if 'overall_survival_months' in df.columns else 'os_months'
event_col = 'overall_survival' if 'overall_survival' in df.columns else 'death_from_cancer'

# Convert event column to binary (1=event/death, 0=censored)
df[event_col] = df[event_col].apply(lambda x: 1 if str(x).lower() in ['yes', '1', 'true', 'dead'] else 0)

# Find all mutation columns that end in '_mut'
mutation_cols = [col for col in df.columns if col.endswith('_mut')]

# Filter survival-ready dataframe
surv_ready_df = df[[duration_col, event_col] + mutation_cols].dropna()

# Create binary flags for each mutation and run log-rank tests
results = []
for gene in mutation_cols:
    # Binary column: 0 if '0', 1 otherwise
    surv_ready_df[f'{gene}_bin'] = surv_ready_df[gene].apply(lambda x: 0 if x == '0' else 1)
    
    mutated = surv_ready_df[surv_ready_df[f'{gene}_bin'] == 1]
    wildtype = surv_ready_df[surv_ready_df[f'{gene}_bin'] == 0]

    if len(mutated) > 10 and len(wildtype) > 10:  # Avoid small groups
        res = logrank_test(
            wildtype[duration_col], mutated[duration_col],
            event_observed_A=wildtype[event_col], event_observed_B=mutated[event_col]
        )
        results.append({
            'gene': gene,
            'p_value': res.p_value,
            'test_statistic': res.test_statistic
        })

# Create result dataframe
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('p_value')
print("\nTop genes by statistical significance:")
print(results_df.head(10))

# Plot top 3 genes
top_genes = results_df.head(3)['gene']
kmf = KaplanMeierFitter()

for gene in top_genes:
    plt.figure(figsize=(8, 6))
    surv_ready_df[f'{gene}_bin'] = surv_ready_df[gene].apply(lambda x: 0 if x == '0' else 1)
    for label, group in surv_ready_df.groupby(f'{gene}_bin'):
        label_name = f"{gene} Mutated" if label == 1 else f"{gene} Wildtype"
        kmf.fit(group[duration_col], group[event_col], label=label_name)
        kmf.plot_survival_function()
    plt.title(f'Survival Analysis by {gene} Mutation Status')
    plt.xlabel('Time (months)')
    plt.ylabel('Survival Probability')
    plt.grid(True)
    plt.show()


#Add Cox Proportional Hazards Model 
# Load dataset
df = pd.read_csv("data/METABRIC_RNA_Mutation.csv")

# Identify survival columns
duration_col = 'overall_survival_months' if 'overall_survival_months' in df.columns else 'os_months'
event_col = 'overall_survival' if 'overall_survival' in df.columns else 'death_from_cancer'

# Convert event to binary
df[event_col] = df[event_col].apply(lambda x: 1 if str(x).lower() in ['yes', '1', 'true', 'dead'] else 0)

# Select mutation columns of interest
selected_genes = ['tp53_mut', 'map3k1_mut', 'tg_mut', 'arid1b_mut']

# Create binary indicators
for gene in selected_genes:
    df[f"{gene}_bin"] = df[gene].apply(lambda x: 0 if x == '0' else 1)

# Prepare input for Cox model
cox_df = df[[duration_col, event_col] + [f"{gene}_bin" for gene in selected_genes]].dropna()

# Fit model
cph = CoxPHFitter()
cph.fit(cox_df, duration_col=duration_col, event_col=event_col)

# Print results
cph.print_summary()

#Gene Expression Feature Engineering

# Step 1: Select numeric columns (assumed to be gene expression + mutation)
df_numeric = df.select_dtypes(include=[np.number])

# Step 2: Remove known clinical numeric columns (keep gene-related features only)
clinical_cols = ['age_at_diagnosis', 'chemotherapy', 'cohort']
df_expression = df_numeric.drop(columns=[col for col in clinical_cols if col in df_numeric.columns])

# Step 3: Compute average gene expression per patient
df['GeneExpression_Average'] = df_expression.mean(axis=1)

# Step 4: Perform PCA to extract top 5 principal components
pca = PCA(n_components=5)
pca_components = pca.fit_transform(df_expression.fillna(0))

# Step 5: Add PCA features to the original DataFrame
for i in range(5):
    df[f'PCA_{i+1}'] = pca_components[:, i]

# Step 6 (optional): Save to a new CSV file
df[['patient_id', 'GeneExpression_Average', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']].to_csv("data/METABRIC_features.csv", index=False)

# Preview new features
print(df[['patient_id', 'GeneExpression_Average', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']].head())


# -------------------------------
# 1. Load and Merge Data
# -------------------------------
df = pd.read_csv(r'data/METABRIC_RNA_Mutation.csv')

features_df = pd.read_csv(r'data/METABRIC_features.csv')

# Ensure patient_id is same type
df["patient_id"] = df["patient_id"].astype(str)
features_df["patient_id"] = features_df["patient_id"].astype(str)

# Merge features
df = df.merge(features_df, on="patient_id", how="inner")

# -------------------------------
# 2. Create Survival Label
# -------------------------------
# Binary label: 1 = survived ≥ 60 months, 0 = short survival
df['Survival_Group'] = df['overall_survival_months'].apply(lambda x: 1 if x >= 60 else 0)

# -------------------------------
# 3. Prepare Features & Target
# -------------------------------
feature_cols = ['GeneExpression_Average', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']
X = df[feature_cols]
y = df['Survival_Group']

# Drop missing values
valid_idx = y.notna() & X.notna().all(axis=1)
X = X[valid_idx]
y = y[valid_idx]

# -------------------------------
# 4. Train/Test Split
# -------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# -------------------------------
# 5. Train Random Forest
# -------------------------------
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# -------------------------------
# 6. Evaluate Model
# -------------------------------
y_pred = rf.predict(X_test)
y_prob = rf.predict_proba(X_test)[:, 1]

# Metrics
roc_auc = roc_auc_score(y_test, y_prob)
conf_matrix = confusion_matrix(y_test, y_pred)
report = classification_report(y_test, y_pred)

# Print evaluation
print(f"ROC AUC Score: {roc_auc:.3f}")
print("Confusion Matrix:")
print(conf_matrix)
print("Classification Report:")
print(report)

# ROC Curve Plot
RocCurveDisplay.from_estimator(rf, X_test, y_test)
plt.title(f"Random Forest ROC Curve (AUC = {roc_auc:.2f})")
plt.grid(True)
plt.show()



# -------------------------------
# 7. SHAP Interpretation (Sample-limited for memory efficiency)
# -------------------------------

# Take a subset of the test set
X_test_sample = X_test.sample(n=100, random_state=42)

# Create the SHAP explainer and compute values
explainer = shap.Explainer(rf, X_train)
shap_values = explainer(X_test_sample)

# Save SHAP summary (bar plot)
shap.summary_plot(shap_values, X_test_sample, plot_type="bar", show=False)
plt.savefig("shap_summary_bar.png", bbox_inches="tight")
plt.clf()

# Save SHAP beeswarm plot
shap.summary_plot(shap_values, X_test_sample, show=False)
plt.savefig("shap_beeswarm.png", bbox_inches="tight")
plt.clf()
