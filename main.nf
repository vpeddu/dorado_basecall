#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { Basecall } from './modules.nf'
include { Demux } from './modules.nf'

params.generate_db = false
params.DEMUX = false
    workflow{
        Basecall( 
            params.INPUT,
            params.MODEL,
            params.BASE
	)
    if (params.DEMUX) {
        Demux(
            Basecall.out,
            file(params.SAMPLESHEET),
            params.KIT,
            params.BASE
        )
    }     
    }
