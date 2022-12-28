# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Hanna Soderstrom"
__copyright__ = "Copyright 2022, Hanna Soderstrom"
__email__ = "hanna.soderstrom@gu.se"
__license__ = "GPL-3"


rule sentieon_bwa_mem:
    input:
        reads=lambda wildcards: alignment_input(wildcards),
    output:
        "sentieon/bwa_mem/{sample}_{flowcell}_{lane}_{barcode}_{type}.output.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
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
        "{rule}: Do stuff on sentieon/{rule}/{wildcards.sample}_{wildcards.type}.input"
    shell:
        "echo {input} > {output}"

rule sentieon_dedup:
    input:
        lambda wildcards: [
            "sentieon/bwa_mem/{sample}_%s_%s_%s_{type}.output.txt" % (u.flowcell, u.lane, u.barcode)
            for u in get_units(units, wildcards)
        ],
    output:
        "sentieon/dedup/{sample}_{type}.output.txt",
    params:
        extra=config.get("sentieon", {}).get("extra", ""),
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
        "{rule}: Do stuff on sentieon/{rule}/{wildcards.sample}_{wildcards.type}.input"
    shell:
        "echo {input} > {output}"
