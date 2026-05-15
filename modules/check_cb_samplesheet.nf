process check_cb_samplesheet { 
    label "samplesheet"

    input:
        path samplesheet

    output:
    path samplesheet, emit: validated_csv

    script:
        """
        check_cb_samplesheet.py --samplesheet ${samplesheet}
        """
        }