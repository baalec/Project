import pandas as pd
import numpy as np

file_path = r'C:\Users\fraga\OneDrive\Desktop\Statistics in genetics\Project\Data_2.xlsx'
Data_2 = pd.read_excel(file_path)
tissue_specificity = Data_2['Tissue Specificity']
threshold_ts = tissue_specificity.quantile(0.75)
top_genes = Data_2[tissue_specificity >= threshold_ts]
expression_columns = [col for col in Data_2.columns if col.endswith('mean')]
trachea_genes = top_genes[top_genes[expression_columns].idxmax(axis=1).str.contains('trachea', case=False)]
eye_genes = top_genes[top_genes[expression_columns].idxmax(axis=1).str.contains('eye', case=False)]
print(f'Number of genes most expressed in trachea: {len(trachea_genes)}')
print(f'Number of genes most expressed in eye: {len(eye_genes)}')
trachea_mean = trachea_genes[expression_columns].mean(axis=1)
eye_mean = eye_genes[expression_columns].mean(axis=1)

def create_binary_coexpression_matrix(mean_series, threshold_ratio=0.8):
    n = len(mean_series)
    matrix = np.zeros((n, n), dtype=int)
    values = mean_series.values
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i, j] = 0
            else:
                ratio = min(values[i], values[j]) / max(values[i], values[j])
                matrix[i, j] = 1 if ratio >= threshold_ratio else 0
    return pd.DataFrame(matrix, index=mean_series.index, columns=mean_series.index)

trachea_coexp = create_binary_coexpression_matrix(trachea_mean, threshold_ratio=0.8)
eye_coexp = create_binary_coexpression_matrix(eye_mean, threshold_ratio=0.8)
print("Trachea Co-expression Binary Matrix shape:", trachea_coexp.shape)
print("Eye Co-expression Binary Matrix shape:", eye_coexp.shape)
trachea_coexp.to_excel('trachea_coexpression_binary.xlsx')
eye_coexp.to_excel('eye_coexpression_binary.xlsx')

def normalize(matrix):
    col_sum = matrix.sum(axis=0)
    col_sum[col_sum == 0] = 1
    return matrix / col_sum

def expansion(matrix, power):
    return np.linalg.matrix_power(matrix, power)

def inflation(matrix, inflation_factor):
    inflated = np.power(matrix, inflation_factor)
    return normalize(inflated)

def mcl(adj_matrix, expansion_power=2, inflation_power=2, max_iterations=100, convergence_check=1e-6):
    M = adj_matrix.astype(float).copy()
    np.fill_diagonal(M, 1)
    M = normalize(M)
    for i in range(max_iterations):
        last = M.copy()
        M = expansion(M, expansion_power)
        M = inflation(M, inflation_power)
        if np.abs(M - last).sum() < convergence_check:
            break
    return M

def get_clusters(mcl_matrix, gene_names, tissue_type, threshold=0.0001):
    clusters = {}
    for i, row in enumerate(mcl_matrix):
        sig = tuple(np.where(row > threshold)[0])
        sig = tuple(int(x) for x in sig)
        if len(sig) > 1:
            clusters.setdefault(sig, []).append(i)
    print(f"{tissue_type} Clusters:")
    for idx, (signature, rows) in enumerate(clusters.items(), start=1):
        signature_gene_names = tuple(gene_names[i] for i in signature)  # mapping signature indices to gene names
        print(f"{tissue_type} Cluster {idx}: {signature_gene_names} size: {len(signature_gene_names)}")  # printing gene names
    return clusters

def export_clusters_for_david(clusters, tissue, gene_names):
    output_filename = f"{tissue}_clusters.david"
    
    with open(output_filename, 'w') as f:
        for cluster in clusters:
            # Get the gene names for the current cluster
            cluster_gene_names = [gene_names[i] for i in cluster]
            
            # Write the indices (optional) and gene names to the file
            f.write(f"({', '.join(map(str, cluster))}):\n")
            for gene_name in cluster_gene_names:
                f.write(f"{gene_name}\n")
            f.write("\n")  # Empty line after each cluster
            
trachea_adj = trachea_coexp.values.astype(float)
eye_adj = eye_coexp.values.astype(float)
trachea_mcl = mcl(trachea_adj, expansion_power=2, inflation_power=2, max_iterations=100)
eye_mcl = mcl(eye_adj, expansion_power=2, inflation_power=2, max_iterations=100)
# Create gene_names list using the order from the co-expression matrix index to ensure proper mapping
trachea_gene_names = [trachea_genes.loc[i, 'name'] for i in trachea_coexp.index]  # corrected gene_names order
eye_gene_names = [eye_genes.loc[i, 'name'] for i in eye_coexp.index]  # corrected gene_names order
trachea_clusters = get_clusters(trachea_mcl, trachea_gene_names, tissue_type="Trachea", threshold=0.0001)
eye_clusters = get_clusters(eye_mcl, eye_gene_names, tissue_type="Eye", threshold=0.0001)
export_clusters_for_david(trachea_clusters, "Trachea", trachea_gene_names)  # pass gene names
export_clusters_for_david(eye_clusters, "Eye", eye_gene_names)  # pass gene names