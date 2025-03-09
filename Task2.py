import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Step 1: Read the dataset into a DataFrame
Data_2 = pd.read_excel(r'C:\Users\fraga\OneDrive\Desktop\Statistics in genetics\Project\Data_2.xlsx')
#print(Data_2.columns)


# Summing male and female expression values across all flies for each gene
male_columns = [col for col in Data_2.columns if 'M1' in col or 'M2' in col]  # Male columns
female_columns = [col for col in Data_2.columns if 'F1' in col or 'F2' in col]  # Female columns

# Ensure the data is numeric for male and female expression
Data_2[male_columns] = Data_2[male_columns].apply(pd.to_numeric, errors='coerce')
Data_2[female_columns] = Data_2[female_columns].apply(pd.to_numeric, errors='coerce')

# Sum the expression values for males and females separately
Data_2['Male_expression'] = Data_2[male_columns].sum(axis=1)
Data_2['Female_expression'] = Data_2[female_columns].sum(axis=1)

# Calculate the log-fold change (LFC)
Data_2['LFC'] = np.log2(Data_2['Male_expression'] / Data_2['Female_expression'])

# Perform t-test for each gene to compare male vs female expression and extract significant genes
significant_genes = []

for index, row in Data_2.iterrows():
    # Ensure the row data for both male and female columns are numeric and not NaN
    male_data = row[male_columns].dropna().apply(pd.to_numeric, errors='coerce')
    female_data = row[female_columns].dropna().apply(pd.to_numeric, errors='coerce')

    # Perform t-test if data is valid (numeric)
    if not male_data.empty and not female_data.empty:
        t_stat_val, p_val_val = stats.ttest_ind(male_data, female_data)
        
        # Filter out genes with LFC under 0.8 or -0.8 (but not low expression)
        if p_val_val < 0.05 and (abs(row['LFC']) >= 1):  # Only include genes with LFC ≥ 0.8 or ≤ -0.9
            Tissue_Specificity = row['Tissue Specificity']  # Get the tissue specificity value
            significant_genes.append((row['name'], row['LFC'], p_val_val, Tissue_Specificity))

    else:
        print(f"name: {row['name']} - Invalid data for t-test")

# Convert significant genes to a DataFrame for further processing
significant_df = pd.DataFrame(significant_genes, columns=['name', 'LFC', 'p_value', 'Tissue Specificity'])

# Display the results of significant genes
#print(significant_df[['name', 'LFC', 'p_value', 'Tissue Specificity']])

num_significant_genes = len(significant_genes)

# Print the number of genes
#print(f'Number of significant genes: {num_significant_genes}')


# Plotting the distribution of LFC values for significant genes
plt.figure(figsize=(10, 6))
plt.hist(significant_df['LFC'], bins=50, color='skyblue', edgecolor='black')
plt.title('Distribution of Log-Fold Change (LFC) for Significant Genes')
plt.xlabel('Log-Fold Change (LFC)')
plt.ylabel('Number of Genes')
plt.grid(True)
plt.show()


# Strip invisible characters or whitespace from all column names
significant_df.columns = significant_df.columns.str.strip()

# Now proceed with the analysis
males_df = significant_df[significant_df['LFC'] > 0]  # Positive LFC genes
females_df = significant_df[significant_df['LFC'] < 0]  # Negative LFC genes

# Pearson and Spearman for males (positive LFC)
spearman_males, p_spearman_males = stats.spearmanr(males_df['LFC'], males_df['Tissue Specificity'])
pearson_males, p_pearson_males = stats.pearsonr(males_df['LFC'], males_df['Tissue Specificity'])

# Pearson and Spearman for females (negative LFC)
spearman_females, p_spearman_females = stats.spearmanr(females_df['LFC'], females_df['Tissue Specificity'])
pearson_females, p_pearson_females = stats.pearsonr(females_df['LFC'], females_df['Tissue Specificity'])


# Print the results
print(f"Spearman Correlation (Males - Positive LFC): {spearman_males}, p-value: {p_spearman_males}")
print(f"Pearson Correlation (Males - Positive LFC): {pearson_males}, p-value: {p_pearson_males}")

print(f"Spearman Correlation (Females - Negative LFC): {spearman_females}, p-value: {p_spearman_females}")
print(f"Pearson Correlation (Females - Negative LFC): {pearson_females}, p-value: {p_pearson_females}")

# Print the mean value of Tissue Specificity for males and females
Specificity_males = males_df['Tissue Specificity'].mean()
Specificity_females = females_df['Tissue Specificity'].mean()

print(f"Mean Tissue Specificity for Males: {Specificity_males}")
print(f"Mean Tissue Specificity for Females: {Specificity_females}")

# Optionally, plotting the correlation data
plt.figure(figsize=(12, 6))

# Plotting for Males (positive LFC)
plt.subplot(1, 2, 1)
plt.scatter(males_df['LFC'], males_df['Tissue Specificity'], color='blue', label="Males (Positive LFC)")
plt.title('LFC vs Tissue Specificity (Males - Positive LFC)')
plt.xlabel('LFC')
plt.ylabel('Tissue Specificity')
plt.legend()

# Plotting for Females (negative LFC)
plt.subplot(1, 2, 2)
plt.scatter(females_df['LFC'], females_df['Tissue Specificity'], color='red', label="Females (Negative LFC)")
plt.title('LFC vs Tissue Specificity (Females - Negative LFC)')
plt.xlabel('LFC')
plt.ylabel('Tissue Specificity')
plt.legend()

plt.tight_layout()
plt.show()
