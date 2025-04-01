#!/usr/bin/env python3
import sys
import os
import pysam
from pathlib import Path
import numpy as np
from collections import Counter

def get_bam_stats(bam_file):
    """Получает базовые статистики из BAM файла."""
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'duplicate_reads': 0,
        'read_lengths': [],
        'insert_sizes': []
    }
    
    bam = pysam.AlignmentFile(bam_file)
    stats['total_reads'] = bam.count()
    
    for read in bam:
        if read.query_length:
            stats['read_lengths'].append(read.query_length)
        if read.is_mapped:
            stats['mapped_reads'] += 1
        else:
            stats['unmapped_reads'] += 1
        if read.is_duplicate:
            stats['duplicate_reads'] += 1
        if read.is_paired and read.template_length:
            stats['insert_sizes'].append(abs(read.template_length))
    
    return stats

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
        print("Usage: collect_metrics.py <bam_file> <metrics_file> <mito_stats> <exome_stats> <onco_stats> <meta_stats> <output_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    metrics_file = sys.argv[2]
    mito_stats_path = sys.argv[3]
    exome_stats_path = sys.argv[4]
    onco_stats_path = sys.argv[5]
    meta_stats_path = sys.argv[6]
    output_file = sys.argv[7]
    
    # Проверка файлов
    for path in [bam_file, metrics_file, mito_stats_path, exome_stats_path, onco_stats_path, meta_stats_path]:
        if not Path(path).exists():
            raise FileNotFoundError(f"File not found: {path}")
    
    try:
        # 1. Базовые статистики BAM
        bam_stats = get_bam_stats(bam_file)
        
        # 2. Статистики митохондриального генома
        mito_stats = read_mito_stats(mito_stats_path)
        
        # 3. Статистики покрытия
        exome_coverage = read_coverage_stats(exome_stats_path)
        onco_coverage = read_coverage_stats(onco_stats_path)
        meta_coverage = read_coverage_stats(meta_stats_path)
        
        # 4. Расчет дополнительных метрик
        read_length_stats = {
            'mean': np.mean(bam_stats['read_lengths']) if bam_stats['read_lengths'] else 0,
            'median': np.median(bam_stats['read_lengths']) if bam_stats['read_lengths'] else 0,
            'mode': Counter(bam_stats['read_lengths']).most_common(1)[0][0] if bam_stats['read_lengths'] else 0
        }
        
        insert_size_stats = {
            'mean': np.mean(bam_stats['insert_sizes']) if bam_stats['insert_sizes'] else 0,
            'median': np.median(bam_stats['insert_sizes']) if bam_stats['insert_sizes'] else 0,
            'mode': Counter(bam_stats['insert_sizes']).most_common(1)[0][0] if bam_stats['insert_sizes'] else 0
        }
        
        # Записываем все метрики в файл
        with open(output_file, 'w') as f:
            # Базовые метрики
            f.write("Basic Metrics:\n")
            f.write(f"Total reads: {bam_stats['total_reads']:,}\n")
            f.write(f"Mapped reads: {bam_stats['mapped_reads']:,} ({bam_stats['mapped_reads']/bam_stats['total_reads']*100:.2f}%)\n")
            f.write(f"Unmapped reads: {bam_stats['unmapped_reads']:,} ({bam_stats['unmapped_reads']/bam_stats['total_reads']*100:.2f}%)\n")
            f.write(f"Duplicate reads: {bam_stats['duplicate_reads']:,} ({bam_stats['duplicate_reads']/bam_stats['total_reads']*100:.2f}%)\n\n")
            
            # Статистики длины ридов
            f.write("Read Length Statistics:\n")
            f.write(f"Mean: {read_length_stats['mean']:.2f} bp\n")
            f.write(f"Median: {read_length_stats['median']:.2f} bp\n")
            f.write(f"Mode: {read_length_stats['mode']:.2f} bp\n\n")
            
            # Статистики размера вставки
            f.write("Insert Size Statistics:\n")
            f.write(f"Mean: {insert_size_stats['mean']:.2f} bp\n")
            f.write(f"Median: {insert_size_stats['median']:.2f} bp\n")
            f.write(f"Mode: {insert_size_stats['mode']:.2f} bp\n\n")
            
            # Митохондриальный геном
            f.write("Mitochondrial Genome:\n")
            f.write(f"Mapped reads: {mito_stats.get('mapped_reads', 0):,}\n")
            f.write(f"Bases mapped: {mito_stats.get('bases_mapped', 0):,}\n")
            f.write(f"Average length: {mito_stats.get('avg_length', 0):.2f} bp\n\n")
            
            # Покрытие
            f.write("Coverage Statistics:\n")
            f.write(f"Exome coverage: {exome_coverage:.2f}x\n")
            f.write(f"Onco panel coverage: {onco_coverage:.2f}x\n")
            f.write(f"Metagenome coverage: {meta_coverage:.2f}x\n")
            
            # Picard метрики
            f.write("\nPicard MarkDuplicates Metrics:\n")
            with open(metrics_file) as mf:
                f.write(mf.read())
                
    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main() 