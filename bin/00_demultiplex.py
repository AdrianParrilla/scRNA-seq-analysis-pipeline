#!/opt/env/bin/python

import os
from glob import glob
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import argparse



def collect_barcodes(alignment_dir) -> dict[str, list[str]]:
    """
    Get the barcodes from a pooled library aligned with cellranger
    input:
    -aligment_dir: path to cellranger output directory

    output:
    - a dictionary with sample names and their barcodes
    """

    samples_pooled = glob(f'{alignment_dir}/outs/per_sample_outs/*')

    sample_barcodes = {}

    for sample in samples_pooled:
        sample_name = os.path.basename(sample)
        pool_basename = os.path.basename(alignment_dir)
        sample_name = sample_name + '-' + pool_basename.split('_')[1] # add pool suffix to sample name

        try:
            barcodes_path = f'{sample}/count/sample_filtered_barcodes.csv'

            barcodes_df = pd.read_csv(barcodes_path, header=None)
            barcodes = barcodes_df[1].values.tolist()

        except:
            print(f'{barcodes_path} does not exists')
            continue

        sample_barcodes[sample_name] = barcodes

    return sample_barcodes


def demultiplex(cb_matrix, sample_barcodes):

    adata_cb = sc.read_10x_h5(cb_matrix)
    adata_cb.var_names_make_unique()

    os.makedirs("samples_demultiplexed", exist_ok=True)

    for sample, barcodes in sample_barcodes.items():

        missing_barcodes = set(barcodes) - set(adata_cb.obs.index.tolist()) # check that all barcodes are included

        if missing_barcodes:
            print(f"Number of missing barcodes for sample {sample}:", len(missing_barcodes))

        valid_barcodes = [b for b in barcodes if b in adata_cb.obs_names]
        adata_sample = adata_cb[valid_barcodes].copy() # subset the unfiltered matrix with the barcodes from cellranger

        short_name = sample.rsplit('-', 1)[0]
        adata_sample.obs_names = [f"{short_name}-{b}" for b in adata_sample.obs_names]

        # Save individually to be collected later
        adata_sample.write_h5ad(f"samples_demultiplexed/{sample}.h5ad")
  

def main(alignment_dir, cb_matrix_raw_path):

    sample_barcodes = collect_barcodes(alignment_dir)

    demultiplex(cb_matrix_raw_path, sample_barcodes)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Demultiplex cellbender results')
    parser.add_argument('--alignment_dir', required=True, help='Path to Cellranger output')
    parser.add_argument('--cb_matrix', required=True, help='Path to Cellbdender output raw h5 matrix')
    
    args = parser.parse_args()
    
    main(args.alignment_dir,
        args.cb_matrix
        )





