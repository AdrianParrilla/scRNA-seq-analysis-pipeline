process QC { 
    label "QC"

    input:
        path data_pkl
        path metadata 
        val sample_key
        val filename

    output:
        path "*.h5ad", emit: adata_QC

    script:
        """
        01_preprocessing.py --data_pkl ${data_pkl} --metadata_path ${metadata} --sample_key ${sample_key} --filename ${filename}
        """
        }