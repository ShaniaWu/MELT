import pandas as pd
import os

project = config["run"]["project"]
samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)

## Target Rule ---
rule all:
    input: 
        # [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam", project=config["run"]["project"], sample=s) for s in samples.index]
        expand("data/{project}/MakeVCF_out/{project}_ALU.final_comp.vcf".format(project=project)),
        expand("data/{project}/MakeVCF_out/{project}_LINE1.final_comp.vcf".format(project=project)),
        expand("data/{project}/MakeVCF_out/{project}_SVA.final_comp.vcf".format(project=project))


rule Preprocess:
    input:                        
        bam = "data/{project}/{project}_{sample}.bam"
    output: 
        disc = "data/{project}/{project}_{sample}.bam.disc",
        bai = "data/{project}/{project}_{sample}.bam.disc.bai",
        fq = "data/{project}/{project}_{sample}.bam.fq"
    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"]
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} Preprocess \
        -bamfile {input.bam} \
        -h {params.genome_ref}
        '''
# java -jar /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/MELT.jar Preprocess \
# -bamfile /home/shania.wu/MELT_me_merged/data/{project}/166_120948.bam \
# -h /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa


rule IndivAnalysis:
    input: 
        bam = "data/{project}/{project}_{sample}.bam",
        disc = "data/{project}/{project}_{sample}.bam.disc",
        bai = "data/{project}/{project}_{sample}.bam.disc.bai",
        fq = "data/{project}/{project}_{sample}.bam.fq"
    output: 
        final_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam",
        final_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam.bai",
        # pulled_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.pulled.sorted.bam",
        # pulled_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.pulled.sorted.bam.bai",
        breaks_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam",
        breaks_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam.bai",

        final_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam",
        final_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam.bai",
        # pulled_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.pulled.sorted.bam",
        # pulled_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.pulled.sorted.bam.bai",
        breaks_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam",
        breaks_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam.bai",

        final_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam",
        final_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam.bai",
        # pulled_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.pulled.sorted.bam",
        # pulled_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.pulled.sorted.bam.bai",
        breaks_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam",
        breaks_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam.bai"

    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"],
        element_ref = config["run"]["element_ref"]
    log:
        "logs/{project}_{sample}.log"
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} IndivAnalysis \
        -b hs37d5/NC_007605 \
        -c 40 \
        -h {params.genome_ref} \
        -bamfile {input.bam}\
        -t {params.element_ref}\
        -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/IndivAnalysis_out/
        '''
# java -jar /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/MELT.jar IndivAnalysis \
# -b hs37d5/NC_007605 \
# -c 40 \
# -h /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa \
# -bamfile /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/166_120948.bam \
# -t /hpf/largeprojects/ccmbio/swu/MELT_me_merged/transposon_file_list.txt \
# -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/IndivAnalysis_out/test/


rule GroupAnalysis:
    input: 
        # final_bam_ALU = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_ALU = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index],

        # final_bam_LINE1 = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_LINE1 = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index],

        # final_bam_SVA = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_SVA = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index]
        "data/{project}/IndivAnalysis_out/"

    output: 
        pre_geno_ALU = "data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv"
        
    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"],
        project = project,
        element_ref = config["run"]["element_ref"]
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} GroupAnalysis \
        -h {params.genome_ref} \
        -n /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed \
        -t {params.element_ref}\
        -discoverydir /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/IndivAnalysis_out/ \
        -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/GroupAnalysis_out/ 
        '''

## Note: cannot use line below as MELT "Genotype" step looks for specific output file name without prefix
# mv data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv data/{project}/GroupAnalysis_out/{params.project}_ALU.pre_geno.tsv \

# java -jar /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/MELT.jar GroupAnalysis \
#     -h /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa \
#     -n /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/add_bed_files/1KGP_Hg19/hg19.genes.bed \
#     -t /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/me_refs/1KGP_Hg19/ALU_MELT.zip \
#     -discoverydir /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/IndivAnalysis_out/ \
#     -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/dataIndivAnalysis_out/ ; \
#     mv data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv data/{project}/GroupAnalysis_out/{wildcards.project}.ALU.pre_geno.tsv



rule Genotype:
    input: 
        pre_geno_ALU = "data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv",
        bam = "data/{project}/{project}_{sample}.bam"
    output: 
        genotype_tsv_ALU = "data/{project}/Genotype_out/{project}_{sample}.ALU.tsv",
        genotype_tsv_LINE1 = "data/{project}/Genotype_out/{project}_{sample}.LINE1.tsv",
        genotype_tsv_SVA = "data/{project}/Genotype_out/{project}_{sample}.SVA.tsv"
    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"],
        element_ref = config["run"]["element_ref"]
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} Genotype \
        -bamfile {input.bam} \
        -h {params.genome_ref} \
        -t {params.element_ref}\
	    -p /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/GroupAnalysis_out/ \
        -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/Genotype_out/
        '''
# java -jar /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT//MELTv2.2.2/MELT.jar Genotype \
# 	-bamfile "/home/shania.wu/MELT_me_merged/data/{project}/166_120948.bam" \
# 	-h /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa  \
# 	-t /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/me_refs/1KGP_Hg19/ALU_MELT.zip \
# 	-p /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/GroupAnalysis_out/ \
# 	-w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/Genotype_out/


rule MakeVCF:
    input: 
        pre_geno_ALU = "data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv",
        genotype_tsv_ALU =  [expand("data/{project}/Genotype_out/{project}_{sample}.ALU.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_LINE1 =  [expand("data/{project}/Genotype_out/{project}_{sample}.LINE1.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_SVA =  [expand("data/{project}/Genotype_out/{project}_{sample}.SVA.tsv".format(project=project, sample=s)) for s in samples.index]
    output: 
        VCF_ALU = "data/{project}/MakeVCF_out/{project}_ALU.final_comp.vcf",
        VCF_LINE1 = "data/{project}/MakeVCF_out/{project}_LINE1.final_comp.vcf",
        VCF_SVA = "data/{project}/MakeVCF_out/{project}_SVA.final_comp.vcf"

    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"],
        project = project,
        element_ref = config["run"]["element_ref"]
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        java -jar {params.melt} MakeVCF \
        -genotypingdir /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/Genotype_out/ \
        -p /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/GroupAnalysis_out/ \
        -h {params.genome_ref} \
        -t {params.element_ref}\
        -w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{wildcards.project}/MakeVCF_out/ ; \

        for me in "ALU" "LINE1" "SVA"; do 
            mv "$me".final_comp.vcf data/{wildcards.project}/MakeVCF_out/{wildcards.project}_"$me".final_comp.vcf
            mv data/{wildcards.project}/MakeVCF_out/"$me".hum.list data/{wildcards.project}/MakeVCF_out/{wildcards.project}_"$me".hum.list
        done 
        '''

# 
# java -jar /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT//MELTv2.2.2/MELT.jar MakeVCF \
# 	-genotypingdir /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/Genotype_out/ \
# 	-h /hpf/largeprojects/ccm_dccforge/dccdipg/Common/genomes/GRCh37d5/GRCh37d5.fa  \
#     -t /hpf/largeprojects/ccmbio/ccmmarvin_shared/genomes/166/MELT/MELTv2.2.2/me_1refs/1KGP_Hg19/ALU_MELT.zip \
# 	-w /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/MakeVCF_out/ \
# 	-p /hpf/largeprojects/ccmbio/swu/MELT_me_merged/data/{project}/GroupAnalysis_out/ 