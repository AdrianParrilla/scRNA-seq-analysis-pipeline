#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include modules
include { parse_cb_samplesheet } from './modules/parse_cb_samplesheet'
include { cellbender } from './modules/cellbender'
include { QC } from './modules/QC'
include { plot_QC } from './modules/plot_QC'
include { integration } from './modules/integration'
include { integration_metrics } from './modules/integration_metrics'

/*
workflow {
    main:
    input_mat = channel.fromPath(params.cellbender.input_mat)
    data_pkl = channel.fromPath(params.QC.data_pkl)
    metadata = channel.fromPath(params.metadata)
    adata = channel.fromPath(params.adata_QC)
    adata_integrated = channel.fromPath(params.adata_integrated)

    //cellbender_results = cellbender(input_mat)

    QC(data_pkl, metadata, params.sample_key)
    plot_QC(QC.out.adata_QC, params.sample_key, params.plot_QC.color_by)
    integration(QC.out.adata_QC, params.filename, params.integration.n_genes, params.integration.label_key, params.integration.batch_key)
    integration_metrics(integration.out.adata_integrated, params.filename, params.integration.label_key, params.integration.batch_key)

    }
    */


workflow QC_WORKFLOW {
    main:
        data_pkl = channel.fromPath(params.QC.data_pkl)
        metadata = channel.fromPath(params.metadata)
        adata = channel.fromPath(params.adata_QC)

        QC(data_pkl, metadata, params.sample_key, params.filename)
        plot_QC(QC.out.adata_QC, params.sample_key, params.plot_QC.color_by, params.filename)
        
}


workflow INTEGRATION_WORKFLOW {
    main:
        adata_annotated = channel.fromPath(params.adata_annotated)

        integration(adata_annotated, params.filename, params.integration.n_genes, params.integration.label_key, params.integration.batch_key)
        integration_metrics(integration.out.adata_integrated, params.filename, params.integration.label_key, params.integration.batch_key)
}


workflow {
    if (params.run_workflow == 'QC') {
        log.info "Executing QC Workflow..."
        QC_WORKFLOW()
        
    } else if (params.run_workflow == 'integration') {
        log.info "Executing Integration Workflow..."
        INTEGRATION_WORKFLOW()
        
    } else {
        // Failsafe in case of a typo in the config
        exit 1, "ERROR: Invalid workflow specified. Set params.run_workflow to 'QC' or 'integration' in your config."
    }
}
