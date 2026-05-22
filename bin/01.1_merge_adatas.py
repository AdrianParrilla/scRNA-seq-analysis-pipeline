#!/opt/env/bin/python

import os
import scanpy as sc
import pandas as pd
import anndata as ad
from glob import glob
import argparse

def parse_metadata(adata, metadata_path, sample_key="sample"):
    """
    Merges metadata from a CSV into an AnnData object.
    """
    print(f"Parsing metadata from {metadata_path}...")

    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Metadata file not found at: {os.path.abspath(metadata_path)}")
    
    # Load metadata
    samples_metadata = pd.read_csv(metadata_path)
    
    # Ensure the join key exists in both objects
    if sample_key not in adata.obs.columns:
        raise KeyError(f"'{sample_key}' not found in adata.obs")
    if sample_key not in samples_metadata.columns:
        raise KeyError(f"'{sample_key}' not found in metadata file")

    # Setting index on metadata allows for a cleaner join
    samples_metadata = samples_metadata.set_index(sample_key)
    
    adata.obs = adata.obs.join(samples_metadata, on=sample_key, how="left")

    # Redefine cell types
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].astype('category')

    # remove samples that are not in the metadata (NA)
    adata = adata[~adata.obs[sample_key].isna()]

    print(f'removed {adata.obs[sample_key].isna().sum()} barcodes with NA values in {sample_key}')

    # reorder sample_key to match metadata order
    desired_order = samples_metadata.index.drop_duplicates().tolist()
    adata.obs[sample_key] = pd.Categorical(adata.obs[sample_key], categories=desired_order, ordered=True)
    sorted_barcodes = adata.obs.sort_values(sample_key).index
    adata = adata[sorted_barcodes]
    
    print("Metadata successfully integrated.")
    return adata



def main(adatas_dir, metadata_path, filename, sample_key="sample"):
    
    search_path = os.path.join(adatas_dir, "*_QC.h5ad")
    files = glob(search_path)

    if not files:
        raise FileNotFoundError(f"No .h5ad files found in {adatas_dir}")

    print(f"Found {len(files)} files. Concatenating...")

    adatas_dict = { 
        os.path.basename(f).split('_QC.h5ad')[0]: sc.read_h5ad(f)  # get sample name
        for f in files 
    }

    # Concatenate all
    adata_concat = ad.concat(adatas_dict, label=sample_key, join="outer")
    adata_concat.obs_names_make_unique()

    # remove genes express in less than 3 cells
    sc.pp.filter_genes(adata_concat, min_cells=3)
    adata_concat = adata_concat[:, ~adata_concat.var.index.str.contains('DEPRECATED_')]

    # setting cell names as column
    adata_concat.obs['cell_name'] = adata_concat.obs.index.values.copy()

    #Add metadata
    adata_concat = parse_metadata(adata_concat, metadata_path, sample_key=sample_key)

    # Set counts layer
    adata_concat.layers['counts'] = adata_concat.X.copy()

    print("Saving adata after QC...")
    adata_concat.write(f'01_{filename}_concat_QC.h5ad')
    print("QC completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC of scRNA seq data')
    parser.add_argument('--adatas_dir', required=True, help='Path to the pkl file containing a dictionary sample:Anndata')
    parser.add_argument('--metadata_path', required=True, help='Path to the sample metadata')
    parser.add_argument('--sample_key', required=True, help='Metadata field containing sample names')
    parser.add_argument('-f', '--filename', type=str, required=False, default='adata', help='Filename for the output files')

    
    args = parser.parse_args()
    
    main(adatas_dir=args.adatas_dir,
            metadata_path=args.metadata_path,
            filename=args.filename,
            sample_key=args.sample_key)
