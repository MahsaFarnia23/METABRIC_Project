import pandas as pd
from sklearn.preprocessing import LabelEncoder

# 1️⃣ Load the Dataset
data = pd.read_csv('METABRIC_RNA_Mutation.csv')

# 2️⃣ Display basic information
print("Number of Rows and Columns:", data.shape)
print("\nColumn Names:\n", data.columns)
print("\nFirst Few Rows:\n", data.head())
print("\nNumber of Missing Values:\n", data.isnull().sum().sum())

# 3️⃣ Check Data Types
print("\nData Types:\n", data.dtypes)

# 4️⃣ Statistical Summary for Numerical and Categorical Columns
print("\nStatistical Summary for Numerical Columns:\n", data.describe())

# Get only the categorical columns
categorical_columns = data.select_dtypes(include=['object']).columns
print("\nStatistical Summary for Categorical Columns:\n", data[categorical_columns].describe())

# 5️⃣ Identify columns with missing values
missing_values = data.isnull().sum()
print("\nColumns with Missing Values:\n", missing_values[missing_values > 0])

# 6️⃣ Fix mixed types: Convert problematic columns to string
# These columns were giving DtypeWarnings
mixed_type_columns = [678, 688, 690, 692]
for col in mixed_type_columns:
    data.iloc[:, col] = data.iloc[:, col].astype(str)

# 7️⃣ Fill Missing Values
# For numerical columns, fill with the mean
numerical_cols = ['age_at_diagnosis', 'tumor_size', 'mutation_count']
for col in numerical_cols:
    data[col] = data[col].fillna(data[col].mean())

# For categorical columns, fill with the mode (most frequent value)
categorical_cols = [
    'type_of_breast_surgery', 'cancer_type_detailed', 'cellularity',
    'er_status_measured_by_ihc', 'neoplasm_histologic_grade',
    'tumor_other_histologic_subtype', 'primary_tumor_laterality',
    'oncotree_code', '3-gene_classifier_subtype', 'tumor_stage'
]

for col in categorical_cols:
    data[col] = data[col].fillna(data[col].mode()[0])

# Special cases
data['death_from_cancer'] = data['death_from_cancer'].fillna(data['death_from_cancer'].mode()[0])
data['tumor_stage'] = data['tumor_stage'].fillna(data['tumor_stage'].mode()[0])

# ✅ Verify if there are still missing values
print("\nMissing values after cleaning:\n", data.isnull().sum().sum())

# 8️⃣ Label Encoding for Binary Columns
# For columns that are binary (2 categories)
binary_cols = ['type_of_breast_surgery', 'cancer_type']
le = LabelEncoder()

for col in binary_cols:
    data[col] = le.fit_transform(data[col])

# 9️⃣ One-Hot Encoding for Multi-Category Columns
# For columns with more than 2 categories, use One-Hot Encoding
multi_category_cols = [
    'cancer_type_detailed', 'cellularity', 'er_status_measured_by_ihc',
    'neoplasm_histologic_grade', 'tumor_other_histologic_subtype',
    'primary_tumor_laterality', 'oncotree_code', '3-gene_classifier_subtype', 'tumor_stage'
]

# Apply One-Hot Encoding
data = pd.get_dummies(data, columns=multi_category_cols, drop_first=True)

# 🔄 Display the first few rows and shape of the DataFrame
print("\nFirst Few Rows After Encoding:\n", data.head())
print("\nData Shape After Encoding:", data.shape)
