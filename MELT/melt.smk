configfile: "config.yaml"

rule Preprocess:
    input: 
        bam = "data/{project}_{sample}_chr22.bam"
    output: 
        disc = "data/{project}_{sample}_chr22.bam.disc",
        bai = "data/{project}_{sample}_chr22.bam.disc.bai",
        fq = "data/{project}_{sample}_chr22.bam.disc.fq"
    params:
        melt = config[run][MELT]
        ref = config[run][ref]
    conda: 
        "../envs/melt.yaml"
    shell:
        java -jar {params.melt} Preprocess \
        -bamfile {input.bam} \
        -h {params.ref}


# rule IndivAnalysis:
#     input:
#     output:
#     params:
