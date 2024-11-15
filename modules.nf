process Basecall {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/", mode: 'symlink', overwrite: true
container  "genomicpariscentre/dorado:0.8.3"
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

/opt/dorado/bin/dorado basecaller ${model} ${pod5} ${dorado_runoptions}\
    -x "cuda:auto" \
    -r \
    --emit-fastq >basecalled.fastq

gzip basecalled.fastq

echo "finished basecalling"

"""

}


process Basecall_custom {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/", mode: 'symlink', overwrite: true
container  "genomicpariscentre/dorado:0.8.3"
cpus 24
memory '128 GB'
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input:
    path pod5
    path model
    val base
    path barcode_arrangement
    path barcode_sequences
    val dorado_runoptions
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

/opt/dorado/bin/dorado basecaller ${model} ${pod5} ${dorado_runoptions}\
    -x "cuda:auto" \
    -r \
    --barcode-arrangement ${barcode_arrangement} \
    --barcode-sequences ${barcode_sequences} \
    --emit-fastq >basecalled.fastq

gzip basecalled.fastq

echo "finished basecalling"

"""

}

process Demux {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/demux/", mode: 'symlink', overwrite: true
cpus 32
memory '128 GB'
container  "genomicpariscentre/dorado:0.8.3"
beforeScript 'chmod o+rw .'

input:
    file  basecalled_fastq
    file samplesheet
    val kit
    val base
//    val modelname
output:
    file  "${base}.barcoded_output/*.fastq.gz"

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

find ${base}.barcoded_output/ -name *.fastq | xargs -I {} gzip {}

#dorado demux basecalled.fastq.gz \
#    --sample-sheet samplesheet.csv \
#    --emit-fastq \
#    --output-dir barcoded_output \
#    --kit-name SQK-NBD114-24 \
#    -t 32

"""
}

process Trim {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/trim/", mode: 'symlink', overwrite: true
cpus 32
memory '32 GB'
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    file demuxed_fastq
//    val modelname
output:
    file  "*.trim.fastq.gz"

script:
"""
#!/bin/bash

ls -lah

echo "dorado version"
/opt/dorado/bin/dorado --version 2>&1

sname=\$(basename -s .fastq.gz ${demuxed_fastq})

            porechop -i ${demuxed_fastq} \
                --format auto \
                -t ${task.cpus} \
                -o \$sname.trim.fastq.gz
"""
}
