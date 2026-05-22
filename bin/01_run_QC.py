#!/opt/env/bin/python

import os
import argparse
from glob import glob
import warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import binom
from itertools import product
import pickle
import logging 
import anndata2ri
from scipy.stats import median_abs_deviation as mad
import rpy2.rinterface_lib.callbacks as rcb
from rpy2.robjects.conversion import localconverter
import rpy2.robjects as ro

warnings.simplefilter("ignore", FutureWarning)

rcb.logger.setLevel(logging.ERROR)


# Load R libraries
ro.r('''
suppressMessages({
    library(Seurat)
    library(scater)
    library(scDblFinder)
    library(SingleCellExperiment)
    library(BiocParallel)
})
''')




def mad_outlier(adata, metric: str, nmads: int, upper_only: bool = False):
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    return (M > np.median(M) + nmads * mad(M))


def flag_outliers(adata):

    bool_vector = (
        mad_outlier(adata, "log1p_total_counts", 5)
        | mad_outlier(adata, "log1p_n_genes_by_counts", 5)
        | mad_outlier(adata, "pct_counts_in_top_20_genes", 5)
        | mad_outlier(adata, "pct_counts_mt", 3, upper_only=True)
        | (adata.obs["pct_counts_mt"] > 5) # setting maximum tolerated mitochondrial percentage
        | (adata.obs["n_genes_by_counts"] < 100) # flagging cells with less than 100 genes
    )

    adata.obs['outlier'] = bool_vector.astype(bool)
    

    return adata


def detect_doublets(adata):
    '''
    Perform doublet detection with scDblFinder

    returns: The input AnnData object with two new .obs columns:
        - 'scDblFinder_score'
        - 'scDblFinder_class'
    '''

    data_mat = adata.X.T

    with localconverter(anndata2ri.converter):
        ro.globalenv["data_mat"] = data_mat
        data_mat = ro.r('as(data_mat, "dgCMatrix")')

    

    ro.r(f'''
    set.seed({123})
    sce <- scDblFinder(
        SingleCellExperiment(list(counts=data_mat))
    )
    ''')

    doublet_score = np.array(ro.r('sce$scDblFinder.score'))
    doublet_class = np.array(ro.r('sce$scDblFinder.class'))

    # Add to AnnData.obs
    adata.obs['scDblFinder_score'] = doublet_score
    adata.obs['scDblFinder_class'] = doublet_class

    adata.obs['scDblFinder_class'] = pd.Categorical(adata.obs['scDblFinder_class']).rename_categories({1: 'singlet', 2: 'doublet'})

    return adata


def qc(adata):
    """
    Perform QC and detect doublets on a given Anndata
    """

    adata.var_names_make_unique()

    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )

    remove = ['total_counts_mt', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb']   

    adata.obs.drop(columns=remove, inplace=True)

    adata = flag_outliers(adata)

    adata = detect_doublets(adata)

    return adata



def main(sample_id, h5ad_file):

    print(f'Loading adata for sample {sample_id}',flush=True)
    adata = sc.read_h5ad(h5ad_file)
    print('adata succesfully loaded', flush=True)

    print('>>> Running QC...', flush=True)
    adata = qc(adata)


    print("Saving adata after QC...")
    adata.write(f'{sample_id}_QC.h5ad')
    print("QC completed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC of scRNA seq data')
    parser.add_argument('-f', '--sample_id', type=str, required=False, default='adata', help='Sample_id for the output files')
    parser.add_argument('--h5ad_file', required=True, help='Path to the h5ad file')

    
    args = parser.parse_args()
    
    main(args.sample_id,
    args.h5ad_file
    )