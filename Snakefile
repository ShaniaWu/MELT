import pandas as pd
import os

project = config["run"]["project"]
samples = pd.read_table(config["run"]["samples"]).set_index("sample", drop=False)

## Target Rule ---
rule all:
    input: 
        # [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam", project=config["run"]["project"], sample=s) for s in samples.index]
        # expand("data/{project}/MakeVCF_out/{project}_ALU.final_comp.vcf".format(project=project)),
        # expand("data/{project}/MakeVCF_out/{project}_LINE1.final_comp.vcf".format(project=project)),
        # expand("data/{project}/MakeVCF_out/{project}_SVA.final_comp.vcf".format(project=project))
        # expand("data/{project}/MakeVCF_out/{project}_MELT_no_ac0.vcf.gz".format(project=project))
        # expand("data/{project}/MakeVCF_out/{project}_MELT_no_symb_allele.vcf".format(project=project))
        expand("data/{project}/snpeff_out/{project}_MELT_snpeff.vcf".format(project=project))


rule Preprocess:
    input:                        
        bam = "data/{project}/{project}_{sample}.bam"
    output: 
        disc = temp("data/{project}/{project}_{sample}.bam.disc"),
        bai = temp("data/{project}/{project}_{sample}.bam.disc.bai"),
        fq = temp("data/{project}/{project}_{sample}.bam.fq")
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


rule IndivAnalysis:
    input: 
        bam = [expand("data/{project}/{project}_{sample}.bam".format(project=project, sample=s)) for s in samples.index],
        disc = [expand("data/{project}/{project}_{sample}.bam.disc".format(project=project, sample=s)) for s in samples.index],
        bai = [expand("data/{project}/{project}_{sample}.bam.disc.bai".format(project=project, sample=s)) for s in samples.index],
        fq = [expand("data/{project}/{project}_{sample}.bam.fq".format(project=project, sample=s)) for s in samples.index]
    output: 
        indiv_dir = temp(directory("data/{project}/IndivAnalysis_out/"))
        # final_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam",
        # final_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam.bai",
        # # pulled_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.pulled.sorted.bam",
        # # pulled_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.pulled.sorted.bam.bai",
        # breaks_bam_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam",
        # breaks_bai_ALU = "data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam.bai",

        # final_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam",
        # final_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam.bai",
        # # pulled_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.pulled.sorted.bam",
        # # pulled_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.pulled.sorted.bam.bai",
        # breaks_bam_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam",
        # breaks_bai_LINE1 = "data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam.bai",

        # final_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam",
        # final_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam.bai",
        # # pulled_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.pulled.sorted.bam",
        # # pulled_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.pulled.sorted.bam.bai",
        # breaks_bam_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam",
        # breaks_bai_SVA = "data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam.bai"

    params:
        melt = config["run"]["MELT"],
        genome_ref = config["run"]["genome_ref"],
        element_ref = config["run"]["element_ref"]
    log:
        "logs/{project}.log"
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
        -w /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/IndivAnalysis_out/
        '''


rule GroupAnalysis:
    input: 
        # indiv_dir = [expand("data/{project}/IndivAnalysis_out/".format(project=project))]
        indiv_dir = "data/{project}/IndivAnalysis_out/"
        # final_bam_ALU = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_ALU = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.ALU.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index],

        # final_bam_LINE1 = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_LINE1 = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.LINE1.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index],

        # final_bam_SVA = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.aligned.final.sorted.bam".format(project=project, sample=s)) for s in samples.index],
        # breaks_bam_SVA = [expand("data/{project}/IndivAnalysis_out/{project}_{sample}.SVA.hum_breaks.sorted.bam".format(project=project, sample=s)) for s in samples.index]

    output: 
        pre_geno_ALU = temp("data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv"),
        pre_geno_LINE1 = temp("data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv"),
        pre_geno_SVA = temp("data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv")
        
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
        -discoverydir {input.indiv_dir} \
        -w /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/GroupAnalysis_out/ 
        '''
## Note: cannot add project name as prefix to output name as MELT "Genotype" step looks for specific output file name without any prefixes


rule Genotype:
    input: 
        pre_geno_ALU = "data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv",
        bam = "data/{project}/{project}_{sample}.bam"
    output: 
        genotype_tsv_ALU = temp("data/{project}/Genotype_out/{project}_{sample}.ALU.tsv"),
        genotype_tsv_LINE1 = temp("data/{project}/Genotype_out/{project}_{sample}.LINE1.tsv"),
        genotype_tsv_SVA = temp("data/{project}/Genotype_out/{project}_{sample}.SVA.tsv")
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
	    -p /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/GroupAnalysis_out/ \
        -w /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/Genotype_out/
        '''


rule MakeVCF:
    input: 
        pre_geno_ALU = "data/{project}/GroupAnalysis_out/ALU.pre_geno.tsv",
        pre_geno_LINE1 = "data/{project}/GroupAnalysis_out/LINE1.pre_geno.tsv",
        pre_geno_SVA = "data/{project}/GroupAnalysis_out/SVA.pre_geno.tsv",
        genotype_tsv_ALU =  [expand("data/{project}/Genotype_out/{project}_{sample}.ALU.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_LINE1 =  [expand("data/{project}/Genotype_out/{project}_{sample}.LINE1.tsv".format(project=project, sample=s)) for s in samples.index],
        genotype_tsv_SVA =  [expand("data/{project}/Genotype_out/{project}_{sample}.SVA.tsv".format(project=project, sample=s)) for s in samples.index]
    output: 
        VCF_ALU = temp("data/{project}/MakeVCF_out/{project}_ALU.final_comp.vcf.gz"),
        VCF_LINE1 = temp("data/{project}/MakeVCF_out/{project}_LINE1.final_comp.vcf.gz"),
        VCF_SVA = temp("data/{project}/MakeVCF_out/{project}_SVA.final_comp.vcf.gz")

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
        -genotypingdir /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/Genotype_out/ \
        -p /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/GroupAnalysis_out/ \
        -h {params.genome_ref} \
        -t {params.element_ref}\
        -w /hpf/largeprojects/ccmbio/swu/MELT/data/{wildcards.project}/MakeVCF_out/ ; 

        for me in "ALU" "LINE1" "SVA"; do 
            mv "$me".final_comp.vcf data/{wildcards.project}/MakeVCF_out/{wildcards.project}_"$me".final_comp.vcf;
            mv data/{wildcards.project}/MakeVCF_out/"$me".hum.list data/{wildcards.project}/MakeVCF_out/{wildcards.project}_"$me".hum.list
        done;

        for file in data/{wildcards.project}/MakeVCF_out/*vcf; do
            bgzip $file; tabix $file.gz; 
        done
        '''


rule MergeVCFs:
    input: 
        VCF_ALU = "data/{project}/MakeVCF_out/{project}_ALU.final_comp.vcf.gz",
        VCF_LINE1 = "data/{project}/MakeVCF_out/{project}_LINE1.final_comp.vcf.gz",
        VCF_SVA = "data/{project}/MakeVCF_out/{project}_SVA.final_comp.vcf.gz"

    output:
        merged_VCF = temp("data/{project}/MakeVCF_out/{project}_MELT.vcf.gz"),
        index = temp("data/{project}/MakeVCF_out/{project}_MELT.vcf.gz.tbi")

    params:
        project = project
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        bcftools concat -a -Oz  -o data/{wildcards.project}/MakeVCF_out/{wildcards.project}_MELT.vcf.gz {input.VCF_ALU} {input.VCF_LINE1} {input.VCF_SVA};
        tabix data/{wildcards.project}/MakeVCF_out/{wildcards.project}_MELT.vcf.gz;
        '''


rule FilterVCF:
    input: 
        merged_VCF = "data/{project}/MakeVCF_out/{project}_MELT.vcf.gz",
        index = "data/{project}/MakeVCF_out/{project}_MELT.vcf.gz.tbi"

    output:
        filtered_VCF = temp("data/{project}/MakeVCF_out/{project}_MELT_no_ac0.vcf.gz")
    params:
        project = project
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        bcftools filter -i 'FILTER!="ac0"' {input.merged_VCF} -Oz -o data/{wildcards.project}/MakeVCF_out/{wildcards.project}_MELT_no_ac0.vcf.gz
        '''


rule ReplaceSymbAllele:
    input: 
        filtered_VCF = "data/{project}/MakeVCF_out/{project}_MELT_no_ac0.vcf.gz"
    output: 
        noSymbAllele = temp("data/{project}/MakeVCF_out/{project}_MELT_no_symb_allele.vcf.gz")
    params:
        project = project
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        PREFIX=`echo {input.filtered_VCF} | sed 's/no_ac0.vcf.gz//'`;

        zcat {input.filtered_VCF} | sed 's/<INS:ME:ALU>/N/g' | sed 's/<INS:ME:LINE1>/N/g' | sed 's/<INS:ME:SVA>/N/g' > "$PREFIX"no_symb_allele.vcf;
        bgzip "$PREFIX"no_symb_allele.vcf;
        '''


rule snpeff:
    input: 
        noSymbAllele = "data/{project}/MakeVCF_out/{project}_MELT_no_symb_allele.vcf.gz"
    output: 
        snpeff_VCF = "data/{project}/snpeff_out/{project}_MELT_snpeff.vcf" 
    params:
        project = project,
        snpeff = config["run"]["snpeff"],
        snpeff_data = config["run"]["snpeff_data"]
    conda: 
        "envs/melt.yaml"
    shell:
        '''
        PREFIX=$(echo {input.noSymbAllele} | sed 's/MakeVCF_out/snpeff_out/' | sed 's/no_symb_allele.vcf.gz//')

        {params.snpeff} \
            -v \
            -Xms750m -Xmx20g  \
            -i vcf -o vcf \
            -dataDir {params.snpeff_data} \
            GRCh37.75 \
            -s "data/{project}/snpeff_out/snpEff_summary.html" \
            {input.noSymbAllele} > "$PREFIX"snpeff.vcf
        '''

rule MELT_report:
    input: 
        snpeff_VCF = "data/{project}/snpeff_out/{project}_MELT_snpeff.vcf" 
    output: 
        report = "data/{project}/report_out/{project}_MELT_report.csv" 
    params:
        project = project,
        protein_coding_genes = config["run"]["protein_coding_genes"],
        hgmd_db = config["run"]["hgmd_db"],
        exon_bed = config["run"]["exon_bed"],
        exac = config["run"]["exac"],
        omim = config["run"]["omim"],
        gnomad = config["run"]["gnomad"],
        biomart = config["run"]["biomart"],
        hpo = config["run"]["hpo_dir"] + project + "_HPO.txt", ###########??
        sv_counts = config["run"]["sv_counts"]

    conda: 
        "envs/melt.yaml"

    script:
        "scripts/MELT_report.py"