process plot_QC { 
    label "plot_QC" 

    publishDir path: { "${out_dir}/QC_results/" }, mode: 'copy', overwrite: true

    input:
        path adata
        val out_dir
        val sample_key
        val color_by
        val filename

    output:
        path "*_QC_violin_plots.png",      emit: qc_violins
        path "*_QC_doublet_count.png",     emit: qc_doublets
        path '*_QC_summary.csv',           emit: qc_summary
        path '*_QC_cell_counts.csv',       emit: qc_cell_counts
        path '*_QC_thresholds.csv',        emit: qc_thresholds
    

    script:
        """
        02_plot_QC.py --adata ${adata} --sample_key ${sample_key} --color_by ${color_by} --filename ${filename}
        """
    }