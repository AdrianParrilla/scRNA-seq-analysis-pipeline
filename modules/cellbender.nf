process cellbender { 
    label "cellbender"

    publishDir path: { "${out_dir}/cellbender/${dataset}" }, mode: 'copy', overwrite: true


    input:
        tuple val(dataset), path(alignment_dir), val(out_dir), val(lr), val(exp), val(total)
         

    output:
        tuple val(dataset), path ("*_cellbender_filtered.h5"), emit: cb_matrix
        tuple val(dataset), path ("*_cellbender.h5")         , emit: cb_matrix_raw
        tuple val(dataset), path ("*.csv")                   , emit: metrics
        tuple val(dataset), path ("*_cellbender.pdf")        , emit: report
        tuple val(dataset), path ("*_cellbender_report.html"), emit: html_report, optional: true
        tuple val(dataset), path ("*_cellbender.log")        , emit: log

    script:
        def input_mat = "${alignment_dir}/outs/multi/count/raw_feature_bc_matrix.h5"

        def lr_flag    = lr ? "--learning-rate ${lr}" : ""
        def exp_flag   = exp ? "--expected-cells ${exp}" : ""
        def total_flag = total ? "--total-droplets-included ${total}" : ""

        """
        cellbender remove-background \\
            --cuda \\
            --input ${input_mat} \\
            --output ${dataset}_cellbender \\
            ${lr_flag} \\
            ${exp_flag} \\
            ${total_flag}
        """
        } 