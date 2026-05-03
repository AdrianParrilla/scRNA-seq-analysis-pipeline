#!/opt/conda/bin/python

import os 
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scib
import argparse
import torch
import jax


def integration_metrics(adata, batch_key, label_key, filename):

    print('\n>>> Calculating integration metrics', flush=True)

    valid_embeds = [r for r in adata.obsm.keys() if r not in ['X_pca', 'X_umap']] # original is called unintegrated

    print(f"Detected embeddings: {valid_embeds}")

    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=valid_embeds,
        bio_conservation_metrics=BioConservation(isolated_labels=True, 
            nmi_ari_cluster_labels_kmeans=True, 
            silhouette_label=True, 
            clisi_knn=True),
        batch_correction_metrics=BatchCorrection(bras=True,
            ilisi_knn=True, 
            kbet_per_label=False, 
            graph_connectivity=True, 
            pcr_comparison=True),
        pre_integrated_embedding_obsm_key='Unintegrated', 
        n_jobs=4
        )

    bm.benchmark()

    metrics = bm.get_results(min_max_scale=False)

    metrics.to_csv(f"{filename}_benchmark_results.csv")

    bm.plot_results_table(
        min_max_scale=False, 
        show=False,
        save_dir=".")

    print('\n>>> Metrics calculation done!', flush=True)


def main(adata_dir, label_key='cell_type', batch_key='batch', filename = 'adata'):

    # Detect PyTorch device
    device_pt = "gpu" if torch.cuda.is_available() else "cpu"
    print(f">>> PyTorch is using device: {device_pt}", flush=True)

    # Detect JAX device and backend
    print(f">>> JAX default backend: {jax.default_backend()}", flush=True)
    print(f">>> JAX available devices: {jax.devices()}", flush=True)

    print('\n>>> Loading adata...', flush=True)
    adata = sc.read_h5ad(adata_dir)
    print('\n>>> adata succesfully loaded', flush=True)

    integration_metrics(adata, batch_key, label_key, filename)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Integrate scRNA seq data from different batches using scVI')
    parser.add_argument('-i', '--adata_dir', type=str, required=True, help='Input directory of the h5ad file containing the anndata object.')
    parser.add_argument('-l', '--label_key', type=str, required=False, default='cell_type', help='Column from adata.obs where the cell types are specified. Default: cell_type')
    parser.add_argument('-b', '--batch_key', type=str, required=False, default='batch', help='Column from adata.obs where the sample batch is specified. Default: batch')
    parser.add_argument('-f', '--filename', type=str, required=False, default='adata', help='Filename for the output files')

    args = parser.parse_args()

    print(
        f"adata_dir: {args.adata_dir}\n"
        f"label_key: {args.label_key}\n"
        f"batch_key: {args.batch_key}\n"
        f"filename: {args.filename}"
    )

    main(
     adata_dir= args.adata_dir,
     label_key= args.label_key,
     batch_key= args.batch_key, 
     filename=args.filename
     )
