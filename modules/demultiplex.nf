process demultiplex { 
    label "demultiplex"
    tag "${dataset}"

    publishDir path: { "${out_dir}/checkpoints/" }, mode: 'copy', overwrite: true


    input:
        tuple val(dataset), path(alignment_dir), val(out_dir), path(cb_matrix_raw)

    output:
        path "samples_demultiplexed/*.h5ad", emit: sample_adatas


    script:
        """
        00_demultiplex.py --alignment_dir ${alignment_dir} --cb_matrix ${cb_matrix_raw}

        """
        } 