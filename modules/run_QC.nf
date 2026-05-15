process run_QC {
    label "run_QC"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(sample_h5ad), val(output_dir)

    output:
    path "${sample_id}_QC.h5ad", emit: qc_h5ad

    script:
        """
        01_run_QC.py --sample_id ${sample_id} --h5ad_file ${sample_h5ad} 
        """
}