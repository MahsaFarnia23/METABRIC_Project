import pandas as pd

# Load the dataset
data = pd.read_csv('METABRIC_RNA_Mutation.csv')

# Display basic information
print("Number of Rows and Columns:", data.shape)
print("\nColumn Names:\n", data.columns)
print("\nFirst Few Rows:\n", data.head())
print("\nNumber of Missing Values:\n", data.isnull().sum().sum())
