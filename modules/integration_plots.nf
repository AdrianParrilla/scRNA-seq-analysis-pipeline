process integration_plots { 
    label "integration_plots"
    
    publishDir path: { "${out_dir}/integration/" }, mode: 'copy', overwrite: true

    input:
        path metrics_csv
        val out_dir

    output:
        path "integration_dot_plot.png"    
        path "radar_plot.png"    

    script:
        """
        04.1_integration_plots.py \\
            --metrics_csv ${metrics_csv}
        """
        }