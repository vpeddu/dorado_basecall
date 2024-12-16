process Basecall {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/", mode: 'copy', overwrite: true
container  "genomicpariscentre/dorado:0.7.2"
cpus 24
memory '128 GB'
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input:
    path pod5
    val model
    val base
//    val modelname
output:
    file "*.basecalled.bam"

script:
"""
#!/bin/bash

ls -lah

ls ${pod5}
ls ${model}

/opt/dorado/bin/dorado basecaller ${model} ${pod5} \
    -x "cuda:auto" \
    -r  > ${base}.basecalled.bam


echo "finished basecalling"

"""
}

process Basecall_RNA {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/basecall/", mode: 'copy', overwrite: true
container  "genomicpariscentre/dorado:0.7.2"
cpus 24
memory '128 GB'
beforeScript 'chmod o+rw .'
label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
input:
    path pod5
    val accuracy 
    val mod
    val base
//    val modelname
output:
    tuple val(base), file("*.basecalled.bam"), val(mod), val(accuracy)

script:
"""
#!/bin/bash

ls -lah

if [ -z "${mod}" ]
then
    /opt/dorado/bin/dorado basecaller ${accuracy} ${pod5} \
    -x "cuda:auto" \
    --estimate-poly-a \
    -r  > ${base}.${accuracy}.basecalled.bam

else
    /opt/dorado/bin/dorado basecaller ${accuracy},${mod} ${pod5} \
    -x "cuda:auto" \
    --estimate-poly-a \
    -r  > ${base}.${mod}.${accuracy}.basecalled.bam

fi

echo "finished basecalling"
"""
}

process Filter_RCS {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/RCS_removed/", mode: 'copy', overwrite: true
cpus 32
memory '64 GB'
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(bam), val(mod), val(accuracy)
    file RCS
//    val modelname
output:
        tuple val(base), file("*.dcsremoved.fastq.gz"), val(mod), val(accuracy)

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

if [ -z "${mod}" ]
then
samtools fastq -TML,BC,MM,pt ${bam} | gzip  > ${base}.${accuracy}.fastq.gz

minimap2 -ax map-ont \
    -y \
	-t ${task.cpus} \
	${RCS} \
    ${base}.${accuracy}.fastq.gz | \
	samtools view -Sb -f 4 - | \
	samtools fastq -TML,BC,MM,pt - | \
	gzip > ${base}.${accuracy}.dcsremoved.fastq.gz

else
samtools fastq -TML,BC,MM,pt ${bam} | gzip  > ${base}.${mod}.${accuracy}.fastq.gz

minimap2 -ax map-ont \
    -y \
	-t ${task.cpus} \
	${RCS} \
    ${base}.${mod}.${accuracy}.fastq.gz | \
	samtools view -Sb -f 4 - | \
	samtools fastq -TML,BC,MM,pt - | \
	gzip > ${base}.${mod}.${accuracy}.dcsremoved.fastq.gz

fi

"""
}

process Chopper {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Chopper/", mode: 'copy', overwrite: true
cpus 32
memory '64 GB'
container  "quay.io/biocontainers/chopper:0.8.0--hdcf5f25_0"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(fastq), val(mod), val(accuracy)
    val chopqual
output:
        tuple val(base), file("*.filtered.fastq.gz"), val(mod), val(accuracy)

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

if [ -z "${mod}" ]
then

    zcat ${fastq} | \
        /usr/local/bin/chopper --quality ${chopqual} \
                    --threads ${task.cpus} | \
                    gzip > ${base}.${accuracy}.filtered.fastq.gz

else
    zcat ${fastq} | \
        /usr/local/bin/chopper --quality ${chopqual} \
                --threads ${task.cpus} | \
                gzip > ${base}.${mod}.${accuracy}.filtered.fastq.gz

fi

"""
}

process Align_genome {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Alignment/", mode: 'copy', overwrite: true
cpus 32
memory '64 GB'
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(rcs_trimmed_fastq), val(mod), val(accuracy)
    file GENOME
//    val modelname
output:
    file  "*.aligned.bam"
    file  "*.aligned.bam.bai"
    tuple val(base), file("*.aligned.bam"), file("*.aligned.bam.bai"), val(mod), val(accuracy)

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

if [ -z "${mod}" ]
then

minimap2 -ax map-ont \
    -y \
	-t ${task.cpus} \
	${GENOME} \
    ${rcs_trimmed_fastq} | \
    samtools view  -Sb - |
    samtools sort -@ 4 - > ${base}.${accuracy}.aligned.bam

samtools index ${base}.${accuracy}.aligned.bam

else

minimap2 -ax map-ont \
    -y \
	-t ${task.cpus} \
	${GENOME} \
    ${rcs_trimmed_fastq} | \
    samtools view  -Sb - |
    samtools sort -@ 4 - > ${base}.${mod}.${accuracy}.aligned.bam

samtools index ${base}.${mod}.${accuracy}.aligned.bam
fi

"""
}

/////live 

process Filter_RCS_live {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/RCS_removed/", mode: 'copy', overwrite: true
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    file fastq
    file(RCS)
//    val modelname
output:
        tuple env(base), file("*.dcsremoved.fastq.gz")
script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

base=\$(basename ${fastq} .fastq.gz)


minimap2 -ax map-ont \
	-t ${task.cpus} \
	${RCS} \
    ${fastq} | \
	samtools view -Sb -f 4 - | \
	samtools fastq - | \
	gzip > \$base.dcsremoved.fastq.gz

"""
}

process Chopper_live {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Chopper/", mode: 'copy', overwrite: true
container  "quay.io/biocontainers/chopper:0.8.0--hdcf5f25_0"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(rcs_trimmed_fastq)
    val chopqual
output:
        tuple val(base), file("*.filtered.fastq.gz")

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam


    zcat ${rcs_trimmed_fastq} | \
        /usr/local/bin/chopper --quality ${chopqual} \
                    --threads ${task.cpus} | \
                    gzip > ${base}.filtered.fastq.gz

"""
}

process Align_genome_live {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Alignment/", mode: 'copy', overwrite: true
cpus 32
memory '64 GB'
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(rcs_trimmed_fastq)
    file GENOME
//    val modelname
output:
    //file  "*.aligned.bam"
    //file  "*.aligned.bam.bai"
    
    //tuple val(base), file("*.aligned.bam"), file("*.aligned.bam.bai")
    file "*.aligned.bam"

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

minimap2 -ax map-ont \
    -y \
	-t ${task.cpus} \
	${GENOME} \
    ${rcs_trimmed_fastq} | \
    samtools view  -Sb - |
    samtools sort -@ 4 - > ${base}.aligned.bam

samtools index ${base}.aligned.bam

"""
}

process Accumulate_live {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Alignment/", mode: 'copy', overwrite: true
cpus 4
memory '64 GB'
container  "alexdhill/complete-seq:latest"
beforeScript 'chmod o+rw .'

input:
    //tuple val(base), file(bam), file(bai)
    file bam
//    val modelname
output:
    file  "*.merged.bam"

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam
echo "iteration ${task.index}"
prev=\$(expr ${task.index} - 1)

echo "prev \$prev"
if [ \$prev -ne 0 ]; then
    samtools merge \
        -@ ${task.cpus} \
        -o ${task.index}.merged.bam \$prev.merged.bam *.aligned.bam 
else
    mv ${bam} ${task.index}.merged.bam
fi



"""
}
///end live


process Modkit {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/Modkit/", mode: 'copy', overwrite: true
cpus 32
memory '64 GB'
container  "ontresearch/modkit:sha849f94ddad41742d5b65b6667d48a532d7323c29"
beforeScript 'chmod o+rw .'

input:
    tuple val(base), file(bam), file(bai), val(mod), val(accuracy)
//    val modelname
output:
    file  "*.readcalls.tsv"

script:
"""
#!/bin/bash
###TODO modbam -> tagged fastq -> minimap2 -T -> filtered_modbam

if [ -z "${mod}" ]
then
    /opt/custflow/epi2meuser/conda/bin/modkit extract ${bam} ${base}.${accuracy}.readcalls.tsv
else
    /opt/custflow/epi2meuser/conda/bin/modkit extract ${bam} null --read-calls ${base}.${mod}.${accuracy}.readcalls.tsv 
fi

"""
}

process Demux {
//conda "${baseDir}/env/env.yml"
publishDir "${params.output}/demux/", mode: 'symlink', overwrite: true
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
    -t ${task.cpus} \
    --no-trim
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
