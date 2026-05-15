process merge_adatas {
    label "merge_adatas"

    publishDir path: { "${out_dir}/checkpoints/" }, mode: 'copy', overwrite: true

    input:
    path "input_h5ads/*"
    path metadata
    val out_dir
    val sample_key
    val filename

    output:
    path "01_${filename}_concat_QC.h5ad", emit: merged_adata

    script:
    """
    01.1_merge_adatas.py \\
        --adatas_dir input_h5ads \\
        --metadata_path ${metadata} \\
        --sample_key ${sample_key} \\
        --filename ${filename}
    """
}