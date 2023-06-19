#!/bin/bash
#SBATCH --job-name=makeVCF_merged
#SBATCH --time=100:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out

PROJECT=166
SMAPLE=120948
BAM="/hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/166_120948.bam"
MELT="/hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/MELT.jar"
GENOME_REF="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa"
ELEMENT_REF="/hpf/largeprojects/ccmbio/swu/MELT_me_merged/transposon_file_list.txt"

module purge

java -jar ${MELT} MakeVCF \
-genotypingdir /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/Genotype_out/ \
-p /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/GroupAnalysis_out/ \
-h ${GENOME_REF} \
-t ${ELEMENT_REF} \
-w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/MakeVCF_out/test ; \

for me in "ALU" "LINE1" "SVA"; do 
    mv "$me".final_comp.vcf data/MakeVCF_out/test/{params.project}_"$me".final_comp.vcf 
    mv data/MakeVCF_out/test/"$me".hum.list data/MakeVCF_out/test/${PROJECT}_"$me".hum.list
done 