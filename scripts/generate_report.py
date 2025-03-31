#!/usr/bin/env python3
import sys
import os
import pysam
from pathlib import Path
import zipfile
import numpy as np
import yaml

def read_mito_stats(mito_stats_path):
    """Читает статистики митохондриального генома."""
    stats = {}
    try:
        with open(mito_stats_path) as f:
            for line in f:
                if 'reads mapped:' in line:
                    stats['mapped_reads'] = int(line.split('\t')[1])
                elif 'bases mapped:' in line:
                    stats['bases_mapped'] = int(line.split('\t')[1])
                elif 'average length:' in line:
                    stats['avg_length'] = float(line.split('\t')[1])
    except Exception as e:
        print(f"WARNING: Failed to read mitochondrial stats: {str(e)}", file=sys.stderr)
    return stats

def read_coverage_stats(stats_path):
    """Читает статистики покрытия."""
    try:
        with open(stats_path) as f:
            # Берем среднее покрытие из первой строки
            line = f.readline().strip().split('\t')
            return float(line[2]) if len(line) > 2 else 0.0
    except Exception as e:
        print(f"WARNING: Failed to read coverage stats: {str(e)}", file=sys.stderr)
        return 0.0

def main():
    if len(sys.argv) != 8:
        raise ValueError("Usage: generate_report.py <input.bam> <metrics.txt> <output.report> <mito_stats> <exome_stats> <onco_stats> <meta_stats>")
    
    bam_path = sys.argv[1]
    metrics_path = sys.argv[2]
    output_path = sys.argv[3]
    mito_stats_path = sys.argv[4]
    exome_stats_path = sys.argv[5]
    onco_stats_path = sys.argv[6]
    meta_stats_path = sys.argv[7]
    
    # Проверка файлов
    for path in [bam_path, metrics_path, mito_stats_path, exome_stats_path, onco_stats_path, meta_stats_path]:
        if not Path(path).exists():
            raise FileNotFoundError(f"File not found: {path}")

    # Получаем имя образца из пути к BAM
    sample_name = Path(bam_path).stem.split('.')[0]
    
    # Путь к FastQC
    fastqc_zip = os.path.join(
        os.path.dirname(bam_path),
        f"../../metrics/raw_fastqc/{sample_name}/{sample_name}_read1_fastqc.zip"
    )

    try:
        # 1. FastQC метрики
        fastqc_data = parse_fastqc_data(fastqc_zip) if Path(fastqc_zip).exists() else {}
        
        # 2. Статистики выравнивания
        total_reads, mapped_reads, non_unique = get_mapping_stats(bam_path)
        
        # 3. Размеры вставок
        mean_insert, median_insert, mode_insert = get_insert_size_stats(bam_path)
        
        # 4. Метрики дубликатов
        dup_percent = parse_picard_metrics(metrics_path)
        
        # 5. Дополнительные метрики
        mito_stats = read_mito_stats(mito_stats_path)
        exome_coverage = read_coverage_stats(exome_stats_path)
        onco_coverage = read_coverage_stats(onco_stats_path)
        meta_coverage = read_coverage_stats(meta_stats_path)
        
        # Генерация отчета
        with open(output_path, 'w') as report:
            report.write(f"Sample: {sample_name}\n")
            report.write("=" * 50 + "\n\n")
            
            # Базовые метрики
            report.write("Basic Metrics:\n")
            report.write("-" * 20 + "\n")
            report.write(f"Total reads: {total_reads:,}\n")
            report.write(f"Read length: {fastqc_data.get('sequence_length', 'N/A')} bp\n")
            report.write(f"Mapped reads: {mapped_reads:,} ({mapped_reads/total_reads*100:.2f}%)\n")
            report.write(f"Adapter content: {fastqc_data.get('adapter_percent', 0.0):.2f}%\n")
            report.write(f"Duplicate rate: {dup_percent:.2f}%\n")
            report.write(f"Non-unique alignments: {non_unique:,} ({non_unique/total_reads*100:.2f}%)\n")
            report.write(f"Overlapping reads (est.): {fastqc_data.get('overlapping_percent', 0.0):.2f}%\n")
            report.write(f"Insert size - Mean: {mean_insert}, Median: {median_insert}, Mode: {mode_insert}\n\n")
            
            # Митохондриальный геном
            report.write("Mitochondrial Genome:\n")
            report.write("-" * 20 + "\n")
            report.write(f"Mapped reads: {mito_stats.get('mapped_reads', 0):,}\n")
            report.write(f"Bases mapped: {mito_stats.get('bases_mapped', 0):,}\n")
            report.write(f"Average length: {mito_stats.get('avg_length', 0):.2f} bp\n\n")
            
            # Экзомная панель
            report.write("Exome Panel:\n")
            report.write("-" * 20 + "\n")
            report.write(f"Average coverage: {exome_coverage:.2f}x\n\n")
            
            # Онко-панель
            report.write("Oncology Panel:\n")
            report.write("-" * 20 + "\n")
            report.write(f"Average coverage: {onco_coverage:.2f}x\n\n")
            
            # Метагеном
            report.write("Metagenome:\n")
            report.write("-" * 20 + "\n")
            report.write(f"Average coverage: {meta_coverage:.2f}x\n\n")
            
            # Picard метрики
            report.write("Picard MarkDuplicates Metrics:\n")
            report.write("-" * 20 + "\n")
            with open(metrics_path) as f:
                report.write(f.read())
        
        # Добавляем тип исследования
        type_file = Path(bam_path).parent.parent / "reports" / f"{sample_name}_type.txt"
        if type_file.exists():
            with open(output_path, 'a') as report:
                research_type = type_file.read_text().strip()
                report.write(f"\nResearch Type: {research_type}\n")

    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 