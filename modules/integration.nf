process integration { 
    label "integration"
    
    publishDir path: { "${out_dir}/integration/" }, mode: 'copy', overwrite: true
    
    input:
        path adata
        val filename
        val n_genes
        val label_key
        val batch_key
        val out_dir

    output:
        path "*.h5ad",        emit: adata_integrated
        path "umap_plots/*_UMAP_*.png",  emit: umap_plots

    script:
        """
        03_integration.py \\
            --adata_dir ${adata} \\
            --filename ${filename} \\
            --n_genes ${n_genes} \\
            --label_key ${label_key} \\
            --batch_key ${batch_key}
        """
        }