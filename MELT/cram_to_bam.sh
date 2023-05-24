#!/bin/bash
#SBATCH --job-name=cram_to_bam
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --output=%x-%j.out

CRAM=$1
BAM=`echo $CRAM | sed 's/.cram/.bam/'`

export REF_CACHE=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/REF_CACHE/GRCh37d5/%2s/%2s/%s:/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/REF_CACHE/hg19/%2s/%2s/%s
export REF_PATH=/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/REF_CACHE/GRCh37d5/%2s/%2s/%s:/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/REF_CACHE/hg19/%2s/%2s/%s

module load samtools/1.10

samtools view -bo $BAM $CRAM
samtools index $BAM


## does this show up