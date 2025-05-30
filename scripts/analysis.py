import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
from lifelines import CoxPHFitter
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

# Set styles for plots
sns.set(style='whitegrid')

# Load the data
df = pd.read_csv(r'C:\Users\farm2103\project\Breast Cancer Gene Expression Profiles (METABRIC)\METABRIC_RNA_Mutation.csv')

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
df = pd.read_csv("METABRIC_RNA_Mutation.csv")

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
df[['patient_id', 'GeneExpression_Average', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']].to_csv("METABRIC_features.csv", index=False)

# Preview new features
print(df[['patient_id', 'GeneExpression_Average', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'PCA_5']].head())
