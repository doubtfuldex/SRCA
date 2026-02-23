import pandas as pd

# Load the data
df = pd.read_csv('movingCSP.csv')

# Initialize an empty list to store the results
results = []

# Define your new parameters
window_size = 5  # change this to your preferred window size
min_residue_types = 3  # change this to your preferred minimum number of different types of residues

# Create an empty dictionary to store the counts
residue_type_counts = {'Hydrophobic': 0, 'Polar charged': 0, 'Polar noncharged': 0, 'Glycine': 0}

# Loop over each set of consecutive residues with size = window_size
for i in range(len(df) - (window_size - 1)):
    subset = df.iloc[i:i+window_size]
    
    # Check if the subset contains at least min_residue_types different types of residues
    if subset['Type of residue'].nunique() >= min_residue_types:
        # Find the residue with the highest CSP
        max_csp_residue = subset.loc[subset['CSPs'].idxmax()]
        
        # Append the result to the list
        results.append((max_csp_residue['Residue'], max_csp_residue['Type of residue'], max_csp_residue['CSPs']))

        # Count the unique residue types in the chunk and update the dictionary
        unique_types = subset['Type of residue'].unique()
        for residue_type in unique_types:
            residue_type_counts[residue_type] += 1

# Convert the results to a DataFrame
results_df = pd.DataFrame(results, columns=['Residue', 'Type of residue', 'CSP'])

# Print the results
print(results_df)
print(residue_type_counts)

# Count the number of times each type of residue appears in the results
counts = results_df['Type of residue'].value_counts()

# Print the counts
print(counts)

