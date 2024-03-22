process Basecall { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/${base}", mode: 'symlink', overwrite: true
container  "genomicpariscentre/dorado:0.5.3"
cpus 24
memory '128 GB'
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input: 
    path pod5
    path model
    val base
//    val modelname
output: 
    file "basecalled.fastq.gz"

script:
"""
#!/bin/bash

ls -lah 


#export CUDA_VISIBLE_DEVICES=0
#export NVIDIA_VISIBLE_DEVICES=0,1,2,3,4,5

ls ${pod5}
ls ${model}

/opt/dorado/bin/dorado basecaller ${model} ${pod5} ${params.dorado_runoptions}\
    -x "cuda:auto" \
    -r \
    --emit-fastq >basecalled.fastq

gzip basecalled.fastq

echo "finished basecalling"

"""

}
process Demux { 
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/demux/${base}", mode: 'symlink', overwrite: true
cpus 32
memory '128 GB'
container  "genomicpariscentre/dorado:0.5.3"
beforeScript 'chmod o+rw .'

input: 
    file  basecalled_fastq
    file samplesheet
    val kit
    val base
//    val modelname
output: 
    file  "${base}.barcoded_output"

script:
"""
#!/bin/bash

ls -lah 

echo "dorado version"
/opt/dorado/bin/dorado --version 2>&1

/opt/dorado/bin/dorado demux ${basecalled_fastq} \
    --emit-fastq \
    --output-dir ${base}.barcoded_output \
    --sample-sheet ${samplesheet} \
    --kit-name ${kit} \
    -t ${task.cpus} 
echo "finished Demux"

#dorado demux basecalled.fastq.gz \
#    --sample-sheet samplesheet.csv \
#    --emit-fastq \
#    --output-dir barcoded_output \
#    --kit-name SQK-NBD114-24 \
#    -t 32

"""
}
