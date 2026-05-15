#!/opt/env/bin/python
    
import scanpy as sc
import pandas as pd
import anndata as ad
from glob import glob

def parse_metadata(adata, metadata_path, sample_key="sample"):
    """
    Merges metadata from a CSV into an AnnData object.
    """
    print(f"Parsing metadata from {metadata_path}...")
    
    try:
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

        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                adata.obs[col] = adata.obs[col].astype('category')
        
        print("Metadata successfully integrated.")
        return adata

    except FileNotFoundError:
        print(f"Error: The file at {metadata_path} was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


def main(adatas_dir, metadata_path, filename, sample_key="sample"):
    
    files = glob("samples_QC/*.h5ad")
    adatas = [sc.read_h5ad(f) for f in files]

    # Concatenate all
    adata_concat = ad.concat(adatas, label="sample", join="outer")
    adata_concat.obs_names_make_unique()

    # remove genes express in less than 3 cells
    sc.pp.filter_genes(adata_concat, min_cells=3)
    adata_concat = adata_concat[:, ~adata_concat.var.index.str.contains('DEPRECATED_')]

    # setting cell names as column
    adatas_concat.obs['cell_name'] = adatas_concat.obs.index.values.copy()

    #Add metadata
    adata = parse_metadata(adatas_concat, metadata_path, sample_key=sample_key)

    # Set counts layer
    adata.layers['counts'] = adata.X

    print("Saving adata after QC...")
    adatas_concat.write(f'01_{filename}_concat_QC.h5ad')
    print("QC completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC of scRNA seq data')
    parser.add_argument('--adatas_dir', required=True, help='Path to the pkl file containing a dictionary sample:Anndata')
    parser.add_argument('--metadata_path', required=True, help='Path to the sample metadata')
    parser.add_argument('--sample_key', required=True, help='Metadata field containing sample names')
    parser.add_argument('-f', '--filename', type=str, required=False, default='adata', help='Filename for the output files')

    
    args = parser.parse_args()
    
    main(args.data_pkl,
        args.metadata_path,
        args.sample_key,
        args.filename)
