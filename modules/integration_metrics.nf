process integration_metrics { 
    label "integration_metrics"
    
    publishDir path: { "${out_dir}/integration/" }, mode: 'copy', overwrite: true
    
    input:
        path adata
        val filename
        val label_key
        val batch_key
        val out_dir

    output:
        path "*_benchmark_results.csv"             , emit: metrics_csv
        path "*.svg"                               , emit: scib_plot

    script:
        """
        04_integration_metrics.py \\
            --adata_dir ${adata} \\
            --filename ${filename} \\
            --label_key ${label_key} \\
            --batch_key ${batch_key}
        """
        }