import os
import pandas as pd

configfile: "config.yaml"

# Проверка обязательных параметров конфига
required_params = ["reference", "input_dir", "output_dir", "threads", "picard_jar", 
                  "reference_regions", "annotation", "windows"]
for param in required_params:
    if param not in config:
        raise ValueError(f"Required parameter '{param}' missing in config.yaml")

# Автоматическое определение образцов
SAMPLES = [f for f in os.listdir(config["input_dir"]) 
          if os.path.isdir(os.path.join(config["input_dir"], f)) and 
          os.path.exists(os.path.join(config["input_dir"], f, f"{f}_read1.fastq.gz")) and
          f != "ZHA261221"]  # Исключаем метагеномный образец

if not SAMPLES:
    raise ValueError("No samples found in input directory!")

rule all:
    input:
        expand(os.path.join(config["output_dir"], "metrics/raw_fastqc/{fastqc_sample}_read1_fastqc/fastqc_data.txt"), fastqc_sample=SAMPLES),
        expand(os.path.join(config["output_dir"], "reports/{report_sample}_final_report.txt"), report_sample=SAMPLES),
        expand(os.path.join(config["output_dir"], "metrics/{metric_sample}_quality_metrics.txt"), metric_sample=SAMPLES)

rule fastqc_raw:
    input:
        r1 = os.path.join(config['input_dir'], "{sample}/{sample}_read1.fastq.gz"),
        r2 = os.path.join(config['input_dir'], "{sample}/{sample}_read2.fastq.gz")
    output:
        data = os.path.join(config['output_dir'], "metrics/raw_fastqc/{sample}_read1_fastqc/fastqc_data.txt"),
        html = os.path.join(config['output_dir'], "metrics/raw_fastqc/{sample}_read1_fastqc.html")
    threads: config['threads']
    log:
        os.path.join(config['output_dir'], "logs/{sample}_fastqc.log")
    shell:
        """set -euo pipefail
        
        # Создаем директории
        mkdir -p {config[output_dir]}/metrics/raw_fastqc
        
        # Запускаем FastQC
        echo "Starting FastQC analysis for {wildcards.sample}" > {log}
        fastqc {input.r1} {input.r2} \
            -o {config[output_dir]}/metrics/raw_fastqc \
            --quiet \
            2>> {log} || exit 1
            
        # Ждем создания zip-файла
        echo "Waiting for FastQC zip file..." >> {log}
        zip_file="{config[output_dir]}/metrics/raw_fastqc/{wildcards.sample}_read1_fastqc.zip"
        while [ ! -f "$zip_file" ]; do
            sleep 1
            echo "Waiting for $zip_file..." >> {log}
        done
        
        # Распаковываем результаты
        echo "Extracting FastQC results..." >> {log}
        unzip -o "$zip_file" -d "{config[output_dir]}/metrics/raw_fastqc/" 2>> {log}
        
        # Проверяем создание файлов
        echo "Checking output files..." >> {log}
        if [ ! -f {output.data} ]; then
            echo "Error: FastQC data file not created" >> {log}
            exit 1
        fi
        
        if [ ! -f {output.html} ]; then
            echo "Error: FastQC HTML file not created" >> {log}
            exit 1
        fi
        
        echo "FastQC analysis completed successfully" >> {log}"""

rule bwa_alignment:
    input:
        r1 = rules.fastqc_raw.input.r1,
        r2 = rules.fastqc_raw.input.r2,
        fastqc = rules.fastqc_raw.output.data
    output:
        temp(os.path.join(config['output_dir'], "aligned/{sample}.sam"))
    threads: config['threads']
    log: os.path.join(config['output_dir'], "logs/{sample}_bwa.log")
    shell:
        """bwa mem -t {threads} \
        -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:lib1\\tPL:ILLUMINA' \
        {config[reference]} \
        {input.r1} {input.r2} > {output} 2> {log}"""

rule process_bam:
    input:
        rules.bwa_alignment.output
    output:
        bam = os.path.join(config['output_dir'], "aligned/{sample}.sorted.bam"),
        bai = os.path.join(config['output_dir'], "aligned/{sample}.sorted.bam.bai"),
        stats = os.path.join(config['output_dir'], "metrics/{sample}_bam_stats.txt")
    threads: config['threads']
    log: os.path.join(config['output_dir'], "logs/{sample}_samtools.log")
    shell:
        """samtools view -@ {threads} -Sb {input} 2>> {log} | \
        samtools sort -@ {threads} -o {output.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        samtools stats {output.bam} > {output.stats} 2>> {log}"""

rule mark_duplicates:
    input:
        rules.process_bam.output.bam
    output:
        bam = os.path.join(config['output_dir'], "aligned/{sample}.dedup.bam"),
        metrics = os.path.join(config['output_dir'], "metrics/{sample}.dups_metrics.txt")
    log:
        os.path.join(config['output_dir'], "logs/{sample}_mark_duplicates.log")
    threads: config['threads']
    shell:
        """set -e
        echo "Starting mark_duplicates for {wildcards.sample}" > {log}
        echo "Input file: {input}" >> {log}
        echo "Output files: {output.bam}, {output.metrics}" >> {log}
        
        # Проверяем входной файл
        if [ ! -f {input} ]; then
            echo "Error: Input file {input} does not exist" >> {log}
            exit 1
        fi
        
        # Создаем директории
        mkdir -p $(dirname {output.bam}) $(dirname {output.metrics})
        
        # Проверяем права доступа
        if [ ! -w $(dirname {output.bam}) ]; then
            echo "Error: No write permission in output directory" >> {log}
            exit 1
        fi
        
        # Запускаем Picard
        echo "Running Picard MarkDuplicates..." >> {log}
        java -jar {config[picard_jar]} MarkDuplicates \
        INPUT={input} \
        OUTPUT={output.bam} \
        METRICS_FILE={output.metrics} \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        2>> {log}
        
        # Проверяем результат
        if [ ! -f {output.bam} ]; then
            echo "Error: Output BAM file was not created" >> {log}
            exit 1
        fi
        
        echo "MarkDuplicates completed successfully" >> {log}"""

rule coverage_analysis:
    input:
        bam = rules.mark_duplicates.output.bam,
        reference_bed = config["reference_regions"]
    output:
        coverage = os.path.join(config['output_dir'], "coverage/{sample}_coverage.txt")
    log:
        os.path.join(config['output_dir'], "logs/{sample}_coverage.log")
    shell:
        """mkdir -p $(dirname {output.coverage})
        bedtools coverage -a {input.reference_bed} -b {input.bam} > {output.coverage} 2>> {log}"""

rule analyze_mitochondrial:
    input:
        bam = rules.mark_duplicates.output.bam,
        annotation = config["annotation"]
    output:
        mito_stats = os.path.join(config['output_dir'], "metrics/{sample}_mito_stats.txt")
    shell:
        """samtools view -h {input.bam} MT | \
        samtools stats | grep ^SN | cut -f 2- > {output.mito_stats}"""

rule analyze_exome:
    input:
        bam = rules.mark_duplicates.output.bam,
        annotation = config["annotation"]
    output:
        exome_stats = os.path.join(config['output_dir'], "metrics/{sample}_exome_stats.txt")
    shell:
        """bedtools coverage -a {input.annotation} -b {input.bam} | \
        grep exon > {output.exome_stats}"""

rule analyze_onco_panel:
    input:
        bam = rules.mark_duplicates.output.bam,
        annotation = config["annotation"]
    output:
        onco_stats = os.path.join(config['output_dir'], "metrics/{sample}_onco_stats.txt")
    shell:
        """bedtools coverage -a {input.annotation} -b {input.bam} | \
        grep -E 'BRCA1|BRCA2|TP53|EGFR|KRAS' > {output.onco_stats}"""

rule analyze_metagenome:
    input:
        bam = rules.mark_duplicates.output.bam,
        windows = config["windows"]
    output:
        meta_stats = os.path.join(config['output_dir'], "metrics/{sample}_meta_stats.txt")
    shell:
        """bedtools coverage -a {input.windows} -b {input.bam} > {output.meta_stats}"""

rule classify_sample:
    input:
        coverage = rules.coverage_analysis.output.coverage,
        annotation = config["annotation"]
    output:
        flag = os.path.join(config['output_dir'], "reports/{sample}_type.txt")
    params:
        output_file = os.path.join(config['output_dir'], "reports/{sample}_type.txt")
    script:
        "scripts/classify_sample.py"

rule analyze_sample:
    input:
        bam = rules.mark_duplicates.output.bam,
        metrics = rules.mark_duplicates.output.metrics,
        sample_type = rules.classify_sample.output.flag,
        mito_stats = rules.analyze_mitochondrial.output.mito_stats,
        exome_stats = rules.analyze_exome.output.exome_stats,
        onco_stats = rules.analyze_onco_panel.output.onco_stats,
        meta_stats = rules.analyze_metagenome.output.meta_stats
    output:
        report = os.path.join(config['output_dir'], "reports/{sample}_final_report.txt")
    log:
        os.path.join(config['output_dir'], "logs/{sample}_analyze.log")
    shell:
        "python scripts/generate_report.py {input.bam} {input.metrics} {output.report} {input.mito_stats} {input.exome_stats} {input.onco_stats} {input.meta_stats} 2>> {log}"

rule collect_quality_metrics:
    input:
        bam = rules.mark_duplicates.output.bam,
        metrics = rules.mark_duplicates.output.metrics,
        mito_stats = rules.analyze_mitochondrial.output.mito_stats,
        exome_stats = rules.analyze_exome.output.exome_stats,
        onco_stats = rules.analyze_onco_panel.output.onco_stats,
        meta_stats = rules.analyze_metagenome.output.meta_stats
    output:
        quality_metrics = os.path.join(config['output_dir'], "metrics/{sample}_quality_metrics.txt")
    shell:
        "python scripts/collect_metrics.py {input.bam} {input.metrics} {input.mito_stats} {input.exome_stats} {input.onco_stats} {input.meta_stats} {output.quality_metrics}"

rule collect_metrics:
    input:
        fastqc = rules.fastqc_raw.output.data,
        bam_stats = rules.process_bam.output.stats
    output:
        metrics = os.path.join(config['output_dir'], "metrics/{sample}_metrics.json")
    params:
        output_file = os.path.join(config['output_dir'], "metrics/{sample}_metrics.json")
    script:
        "scripts/collect_metrics.py"
