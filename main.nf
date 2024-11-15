#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { Basecall } from './modules.nf'
include { Basecall_custom } from './modules.nf'
include { Demux } from './modules.nf'
include { Trim } from './modules.nf'

params.generate_db = false
params.DEMUX = false
params.STANDARD = false
params.CUSTOM_DEMUX = false
params.dorado_runoptions = ""
    workflow{
    if (params.STANDARD){
        Basecall(
            params.INPUT,
            params.MODEL,
            params.BASE
	)
    }
    if (params.DEMUX) {
        Demux(
            Basecall.out,
            file(params.SAMPLESHEET),
            params.KIT,
            params.BASE
        )
        Trim(
            Demux.out.flatten()
        )
    }
    if (params.CUSTOM_DEMUX) {
        Basecall_custom(
            params.INPUT,
            params.MODEL,
            params.BASE,
            params.BARCODE_ARRANGEMENT,
            params.BARCODE_SEQUENCES,
            params.dorado_runoptions
	)   
        
        }
    }
