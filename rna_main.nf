#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.recursion=true

include { Basecall } from './modules.nf'
include { Basecall_RNA } from './modules.nf'
include { Align_genome } from './modules.nf'
include { Chopper } from './modules.nf'
include { Modkit } from './modules.nf'
// include { Basecall_m6a } from './modules.nf'
// include { Basecall_pseU } from './modules.nf'
include { Filter_RCS } from './modules.nf'
include { Demux } from './modules.nf'
include { Trim } from './modules.nf'
include { Filter_RCS_live } from './modules.nf'
include { Chopper_live } from './modules.nf'
include { Align_genome_live } from './modules.nf'
include { Accumulate_live } from './modules.nf'

params.generate_db = false
params.DEMUX = false
params.dRNA = false
params.LIVE = false
params.chopper_quality = 10

    workflow{
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
    Mods_ch = channel.from('', 'pseU', 'm6A')
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
    if (params.DEMUX) {
    Basecall(
            params.INPUT,
            params.MODEL,
            params.BASE
	)
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
    }
