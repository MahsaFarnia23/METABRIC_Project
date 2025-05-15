import pandas as pd

# Load the dataset
data = pd.read_csv('METABRIC_RNA_Mutation.csv')

# Display basic information
print("Number of Rows and Columns:", data.shape)
print("\nColumn Names:\n", data.columns)
print("\nFirst Few Rows:\n", data.head())
print("\nNumber of Missing Values:\n", data.isnull().sum().sum())
# Check data types
print("\nData Types:\n", data.dtypes)

# Get basic statistics for numerical columns
print("\nStatistical Summary for Numerical Columns:\n", data.describe())

# Get basic statistics for categorical columns
categorical_columns = data.select_dtypes(include=['object']).columns
print("\nStatistical Summary for Categorical Columns:\n", data[categorical_columns].describe())
# Check columns with missing values
missing_values = data.isnull().sum()
print("Columns with missing values:\n", missing_values[missing_values > 0])
# Convert problematic columns to string to prevent issues
mixed_type_columns = [678, 688, 690, 692]
for col in mixed_type_columns:
    data.iloc[:, col] = data.iloc[:, col].astype(str)

numerical_cols = ['age_at_diagnosis', 'tumor_size', 'mutation_count']
for col in numerical_cols:
    data[col].fillna(data[col].mean(), inplace=True)

categorical_cols = [
    'type_of_breast_surgery', 'cancer_type_detailed', 'cellularity',
    'er_status_measured_by_ihc', 'neoplasm_histologic_grade',
    'tumor_other_histologic_subtype', 'primary_tumor_laterality',
    'oncotree_code', '3-gene_classifier_subtype', 'tumor_stage'
]

for col in categorical_cols:
    data[col].fillna(data[col].mode()[0], inplace=True)

data['death_from_cancer'].fillna(data['death_from_cancer'].mode()[0], inplace=True)
data['tumor_stage'].fillna(data['tumor_stage'].mode()[0], inplace=True)

# Check if there are still missing values
print("Missing values after cleaning:\n", data.isnull().sum().sum())

