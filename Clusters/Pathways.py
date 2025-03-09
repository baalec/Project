import os

# Define the directory where your files are stored
directory = r"C:\Users\fraga\OneDrive\Desktop\Statistics in genetics\Project\Clusters"  # Your actual directory path

# Function to read files and extract pathway names
def extract_significant_pathways(file_path):
    significant_pathways = []
    
    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that are headers or don't contain relevant data
            if line.startswith("Sublist") or not line.strip():
                continue
            
            # Split the line into columns
            columns = line.split("\t")
            
            # Extract the p-value and count % from the line (assuming they are in the right columns)
            try:
                p_value = float(columns[7])  # Column with p-value (adjust index if needed)
                count_percent = float(columns[6])  # Column with count % (adjust index if needed)
                
                # Check if the conditions are met
                if p_value < 0.01 and count_percent > 15:
                    pathway_name = columns[2]  # Assuming the pathway name is in column 2 (adjust if needed)
                    significant_pathways.append(pathway_name)
            except (ValueError, IndexError):
                # Handle any lines that don't have the expected format
                continue

    return significant_pathways

# Iterate through all files in the directory
for filename in os.listdir(directory):
    # Check if the file matches the expected pattern (eye_cluster or trachea_cluster)
    if filename.startswith("eye_cluster") or filename.startswith("trachea_cluster"):
        file_path = os.path.join(directory, filename)
        
        # Extract significant pathways from the file
        pathways = extract_significant_pathways(file_path)
        
        # Print the results for each file
        if pathways:
            print(f"Significant pathways in {filename}:")
            for pathway in pathways:
                print(f"  - {pathway}")
        else:
            print(f"No significant pathways found in {filename}.")
