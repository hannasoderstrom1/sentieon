# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2022, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"


rule sentieon_bwa_mem:
    input:
        reads=lambda wildcards: alignment_input(wildcards),
        #reference=config.get("sentieon", {}).get("reference", ""),
    output:
        bam = "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.bam",
        #"sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Align fastq files {input.reads} using Sentieon bwa mem against {params.reference}"
    shell:
        #"echo {input} > {output}"
        "{params.sentieon} bwa mem "
            "-M -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' "
            "-t {threads} {params.reference} {input.reads} "
        "| {params.sentieon} util sort -o {output.bam} -t {threads} --sam2bam -i -"

rule sentieon_dedup:
    input:
        lambda wildcards: [
            "sentieon/bwa_mem/{sample}_%s_%s_%s_{type}.bam" % (u.flowcell, u.lane, u.barcode)
            for u in get_units(units, wildcards)
        ],
    output:
        "sentieon/dedup/{sample}_{type}_DEDUP.bam",
        "sentieon/dedup/{sample}_{type}_DEDUP.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
    log:
        "sentieon/dedup/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/dedup/{sample}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Mark/remove duplicate reads in bam file {input} using Sentieon dedup algorithm"
    shell:
        #"echo {input} > {output}"
        "{params.sentieon} driver -t {threads} "
            "-i {input} "
            "--algo LocusCollector "
            "--fun score_info "
            "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt ;"
        "{params.sentieon} driver "
            "-t {threads} "
            "-i {input} "
            "--algo Dedup "
            "--rmdup "
            "--score_info sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP_score.txt "
            "--metrics sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.txt "
            "sentieon/dedup/{wildcards.sample}_{wildcards.type}_DEDUP.bam"


rule sentieon_realigner:
    input:
        bam="sentieon/dedup/{sample}_{type}_DEDUP.bam",
        #reference=config.get("sentieon", {}).get("reference", ""),
    output:
        "sentieon/realign/{sample}_{type}_REALIGNED.bam",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        mills=config.get("sentieon", {}).get("mills", ""),
    log:
        "sentieon/realign/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/realign/{sample}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Indel realignment of bam file {input.bam} using Sentieon realigner"
    shell:
        "{params.sentieon} driver -t {threads} -r {params.reference} -i {input.bam} --algo Realigner -k {params.mills} {output}"


rule sentieon_qualcal:
    input:
        bam="sentieon/realign/{sample}_{type}_REALIGNED.bam",
        #reference=config.get("sentieon", {}).get("reference", ""),
    output:
        "sentieon/qualcal/{sample}_{type}_RECAL_DATA.TABLE",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        mills=config.get("sentieon", {}).get("mills", ""),
        dbsnp=config.get("sentieon", {}).get("dbsnp", ""),
    log:
        "sentieon/qualcal/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "sentieon/qualcal/{sample}_{type}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("sentieon", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Calculate recalibration table of {input.bam} using Sentieon QualCal algorithm"
    shell:
        "{params.sentieon} driver -t {threads} -r {params.reference} -i {input.bam} --algo QualCal -k {params.mills} -k {params.dbsnp} {output}"


rule sentieon_tnscope:
    input:
        #tumorbam=expand("sentieon/realign/{sample}_T_REALIGNED.bam", sample=tumor_sample),
        tumorbam="sentieon/realign/{tumor_sample}_T_REALIGNED.bam",
        #normalbam="sentieon/realign/{normal_sample}_N_REALIGNED.bam",
        normalbam="sentieon/realign/%s_N_REALIGNED.bam" % (normal_sample),
        tumortable="sentieon/qualcal/{tumor_sample}_T_RECAL_DATA.TABLE",
        normaltable="sentieon/qualcal/%s_N_RECAL_DATA.TABLE" % (normal_sample),
    output:
        #tnscope = expand("sentieon/tnscope/{sample}_TNscope_tn.vcf", sample=tumor_sample),
        tnscope = "sentieon/tnscope/{tumor_sample}_TNscope_tn.vcf",
        tnscope_bam = "sentieon/tnscope/{tumor_sample}_REALIGNED_realignedTNscope.bam",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
        reference=config.get("sentieon", {}).get("reference", ""),
        sentieon=config.get("sentieon", {}).get("sentieon", ""),
        callsettings=config.get("sentieon", {}).get("tnscope_settings", ""),
    log:
        #expand("sentieon/tnscope/{sample}.output.log", sample=tumor_sample),
        "sentieon/tnscope/{tumor_sample}.output.log",
    benchmark:
        repeat(
            #expand("sentieon/tnscope/{sample}.output.benchmark.tsv", sample=tumor_sample),
            "sentieon/tnscope/{tumor_sample}.output.benchmark.tsv",
            config.get("sentieon", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("tnscope", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sentieon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sentieon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sentieon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sentieon", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sentieon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sentieon", {}).get("container", config["default_container"])
    conda:
        "../envs/sentieon.yaml"
    message:
        "{rule}: Call SNVs and structural variants in {input.tumorbam} using matched normal {input.normalbam} using Sentieon TNScope"
    shell:
        #"echo {input.tumorbam} > {output.tnscope}"
        "{params.sentieon} driver -t {threads} -r {params.reference} "
            "-i {input.tumorbam} -q {input.tumortable} -i {input.normalbam} -q {input.normaltable} "
            "--algo TNscope --tumor_sample {tumor_sample} --normal_sample {normal_sample} --bam_output {output.tnscope_bam} "
            "{params.callsettings} {output.tnscope}"

