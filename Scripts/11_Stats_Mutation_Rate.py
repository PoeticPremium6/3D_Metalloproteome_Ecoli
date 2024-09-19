import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Specify the necessary columns for analysis
metal_columns = ['UniProt_ID', 'Gene_Length', 'Metal_type']  # Add any other relevant columns you need
non_metal_columns = ['UniProt_ID', 'Gene_Length']  # Adjust for non-metal data

# Load the datasets in chunks and process them
chunk_size = 100000  # Adjust chunk size based on memory capacity
metal_data_chunks = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_metals_final.csv',
                                usecols=metal_columns, chunksize=chunk_size, low_memory=True)
non_metal_data_chunks = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\refined_integrated_dataset_nonmetal.csv',
                                    usecols=non_metal_columns, chunksize=chunk_size, low_memory=True)

# Concatenate chunks into DataFrames
metal_data = pd.concat(metal_data_chunks, ignore_index=True)
non_metal_data = pd.concat(non_metal_data_chunks, ignore_index=True)

# Load mutation data, ensuring only necessary columns are loaded
mutation_data = pd.read_csv('C:\\Users\\jonat\\OneDrive - University of Glasgow\\Metalloproteome\\Submission\\categorized_dataset.csv',
                            usecols=['Entry', 'Grantham_Category'], low_memory=True)

# Rename 'Entry' to 'UniProt_ID' for consistency
mutation_data.rename(columns={'Entry': 'UniProt_ID'}, inplace=True)

# Convert columns to memory-efficient data types
metal_data['UniProt_ID'] = metal_data['UniProt_ID'].astype('category')
non_metal_data['UniProt_ID'] = non_metal_data['UniProt_ID'].astype('category')
mutation_data['UniProt_ID'] = mutation_data['UniProt_ID'].astype('category')
mutation_data['Grantham_Category'] = mutation_data['Grantham_Category'].astype('category')

# Count the number of mutations for each UniProt ID
mutation_counts = mutation_data['UniProt_ID'].value_counts().reset_index()
mutation_counts.columns = ['UniProt_ID', 'Number_of_Mutations']

# Merge mutation counts with metal and non-metal data
metal_mutations = pd.merge(metal_data, mutation_counts, on='UniProt_ID', how='left')
non_metal_mutations = pd.merge(non_metal_data, mutation_counts, on='UniProt_ID', how='left')

# Merge Grantham_Category from the mutation data
metal_mutations = pd.merge(metal_mutations, mutation_data[['UniProt_ID', 'Grantham_Category']], on='UniProt_ID', how='left')
non_metal_mutations = pd.merge(non_metal_mutations, mutation_data[['UniProt_ID', 'Grantham_Category']], on='UniProt_ID', how='left')

# Calculate mutation rates (mutations per gene length)
metal_mutations['Mutation_Rate'] = metal_mutations['Number_of_Mutations'] / metal_mutations['Gene_Length']
non_metal_mutations['Mutation_Rate'] = non_metal_mutations['Number_of_Mutations'] / non_metal_mutations['Gene_Length']

# Group by Grantham_Category and Metal_Type, then calculate average mutation rate
mutation_rate_metal = metal_mutations.groupby(['Grantham_Category', 'Metal_type'])['Mutation_Rate'].mean().unstack(fill_value=0)

# Add the non-metal data as a separate column
non_metal_avg_rate = non_metal_mutations.groupby('Grantham_Category')['Mutation_Rate'].mean()
mutation_rate_metal['Non-Metal'] = non_metal_avg_rate

# Apply log scale transformation
mutation_rate_metal_log = np.log1p(mutation_rate_metal)

# Flip the y-axis by sorting the index in descending order
mutation_rate_metal_log = mutation_rate_metal_log.sort_index(ascending=False)

# Save the result to a file
mutation_rate_metal_log.to_csv('mutation_rate_metal_log.csv')

# Plot (optional)
plt.figure(figsize=(12, 8))
sns.heatmap(mutation_rate_metal_log, cmap='viridis_r', linewidths=0.5)
plt.title('Log-scaled Mutation Rates by Grantham Category and Metal Type')
plt.ylabel('Grantham Category')
plt.xlabel('Metal Type')
plt.show()

import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Assuming 'mutation_rate_metal_log' is the log-scaled DataFrame you have calculated earlier
# We will melt the DataFrame into long form to make it suitable for statistical testing

# Melt the mutation_rate_metal_log dataframe
mutation_rate_long = mutation_rate_metal_log.reset_index().melt(id_vars='Grantham_Category', var_name='Metal', value_name='Mutation_Rate')

# Drop any missing values (NaNs) from the dataset
mutation_rate_long = mutation_rate_long.dropna()

# Step 1: ANOVA to check for overall differences across metals
anova_result = stats.f_oneway(*[mutation_rate_long[mutation_rate_long['Metal'] == metal]['Mutation_Rate'] for metal in mutation_rate_long['Metal'].unique()])

print(f"ANOVA result: F-statistic = {anova_result.statistic}, p-value = {anova_result.pvalue}")

# Step 2: If ANOVA shows significant differences, proceed to post-hoc test
if anova_result.pvalue < 0.05:
    print("ANOVA is significant. Proceeding with post-hoc tests...")

    # Step 3: Tukey's HSD test for pairwise comparison
    tukey_result = pairwise_tukeyhsd(endog=mutation_rate_long['Mutation_Rate'], groups=mutation_rate_long['Metal'], alpha=0.05)
    
    # Print the summary of Tukey HSD test
    print(tukey_result.summary())
    
    # Optional: Visualize the results of Tukey HSD test
    tukey_result.plot_simultaneous()
else:
    print("ANOVA is not significant. No further testing required.")
