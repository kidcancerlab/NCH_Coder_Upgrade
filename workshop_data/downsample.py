import h5py
import pandas as pd
import os
import gzip
import numpy as np
import json

# Taking in input/scRNA/mouse_brain_ref_data.hdf5 and
# input/scRNA/mouse_brain_metadata.csv
# and downsampling to 1000 cells per cell type
# cell type is in subclass_label column
file_path = 'input/scRNA/mouse_brain_ref_data.hdf5'
with h5py.File(file_path, 'r') as hdf:
    data_group = hdf['data']

    # Get gene names
    genes = data_group['gene'][:]

    # Get sample names
    samples = data_group['samples'][:]
    samples = [s.decode('utf-8') for s in samples]

    # Access counts within the 'data' group
    dataset = data_group['counts']
    data = dataset[:]

meta_data = pd.read_csv('input/scRNA/mouse_brain_metadata.csv', index_col=0)

# Randomly sample up to 1000 cells per cell type, and keep these indices
sample_barcodes = []
max_n_cells = 500
for cell_type in meta_data['subclass_label'].unique():
    cell_type_sub = \
        meta_data[meta_data['subclass_label'] == cell_type]
    sampled_sub = \
        cell_type_sub.sample(
            n=min(max_n_cells, len(cell_type_sub)),
            random_state=0
        )
    sample_barcodes.extend(
        sampled_sub.index.tolist()
    )

new_meta_data = meta_data.loc[sample_barcodes]

counts_indices = [samples.index(s) for s in sample_barcodes]

new_samples = [samples[i] for i in counts_indices]
new_counts = data[:, counts_indices]

output_dir = "cellranger_output"
# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# In counts matrix features should be rows, samples are columns

def write_matrix_10x(matrix, output_file, metadata=None):
    """
    Writes a NumPy matrix to the 10X Cell Ranger count output format.
    Args:
        matrix (np.ndarray): The matrix to write.
        output_file (str): The path to the output file.
        metadata (dict, optional): Optional metadata to include in the header.
    """
    rows, cols = matrix.shape
    nnz = np.count_nonzero(matrix)  # Count non-zero entries
    with open(output_file, 'w') as f:
        # Write the header
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        # Write metadata JSON
        if metadata:
            metadata_str = json.dumps(metadata)
            f.write(f"%metadata_json: {metadata_str}\n")
        else:
            f.write("%metadata_json: {}\n")
        # Write the dimensions and number of non-zero entries
        f.write(f"{rows} {cols} {nnz}\n")
        # Write the non-zero entries in coordinate format (1-based indexing)
        for i in range(rows):
            for j in range(cols):
                if matrix[i, j] != 0:
                    f.write(f"{i + 1} {j + 1} {matrix[i, j]}\n")

# Write the matrix.mtx file in the format that Cell Ranger expects
write_matrix_10x(new_counts, os.path.join(output_dir, "matrix.mtx"))
# Compress matrix.mtx
os.system(f"pigz {os.path.join(output_dir, 'matrix.mtx')}")

# Write the features.tsv file
with gzip.open(os.path.join(output_dir, "features.tsv.gz"), "wt") as f:
    for gene in genes:
        f.write(
            gene.decode('utf-8') + '\t' \
            + gene.decode('utf-8') + '\t' \
            + "Gene Expression\n"
        )

# Write the barcodes.tsv file
with gzip.open(os.path.join(output_dir, "barcodes.tsv.gz"), "wt") as f:
    for barcode in new_samples:
        f.write(barcode + "\n")

# Write out the downsampled metadata
with gzip.open(os.path.join(output_dir, "metadata.tsv.gz"), "wt") as f:
    new_meta_data.to_csv(f, sep='\t', index=True)
