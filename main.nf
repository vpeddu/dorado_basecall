#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { Basecall } from './modules.nf'
include { Basecall_custom } from './modules.nf'
include { Demux_custom } from './modules.nf'
include { Demux } from './modules.nf'
include { Trim } from './modules.nf'
include { Chopper_live } from './modules.nf'
include { Basecall_RNA } from './modules.nf'
include { Chopper } from './modules.nf'
include { Modkit } from './modules.nf'
include { Align_genome } from './modules.nf'
include { Filter_RCS } from './modules.nf'
include { Align_genome_live } from './modules.nf'
include { Filter_RCS_live } from './modules.nf'
include { Accumulate_live } from './modules.nf'

params.generate_db = false
params.DEMUX = false
params.STANDARD = false
params.CUSTOM_DEMUX = false
params.dRNA = false
params.LIVE = false
params.dorado_runoptions = ""
params.chopper_quality = 10

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
        Demux_custom(
            Basecall_custom.out,
            file(params.SAMPLESHEET),
            params.BASE,
            params.BARCODE_ARRANGEMENT,
            params.BARCODE_SEQUENCES,
        ) 
        
        }
    if (params.LIVE){
    input_read_Ch = Channel.watchPath("${params.INPUT}*.fastq.gz")
    Filter_RCS_live(
        input_read_Ch,
        file(params.RCS)
        )
    Chopper_live(
        Filter_RCS_live.out,
        params.chopper_quality
        )
    Align_genome_live(
        Chopper_live.out,
        file(params.GENOME)
        )
    Accumulate_live.scan(
        Align_genome_live.out
        )
    }

    if (params.dRNA){
    Mods_ch = channel.from('', 'pseU', 'm6A_DRACH','inosine_m6A', 'm5C')
    Basecall_RNA(
        params.INPUT,
        params.ACCURACY,
        Mods_ch,
        params.BASE
        )
    Filter_RCS(
        Basecall_RNA.out,
        file(params.RCS)
        )
    Chopper(
        Filter_RCS.out,
        params.chopper_quality
        )
    Align_genome(
        Chopper.out,
        file(params.GENOME)
        )
    Modkit(
        Align_genome.out[2]
        )
    }
    }
