#!/usr/bin/env python3
import sys
import pandas as pd
import pysam
from collections import Counter

def get_bam_stats(bam_file):
    bam = pysam.AlignmentFile(bam_file)
    stats = {
        'total_reads': bam.count(),
        'read_lengths': [],
        'insert_sizes': [],
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'duplicate_reads': 0
    }
    
    for read in bam:
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

def calculate_metrics(stats):
    metrics = {
        'total_reads': stats['total_reads'],
        'average_read_length': sum(stats['read_lengths']) / len(stats['read_lengths']) if stats['read_lengths'] else 0,
        'mapped_percentage': (stats['mapped_reads'] / stats['total_reads'] * 100) if stats['total_reads'] > 0 else 0,
        'duplicate_percentage': (stats['duplicate_reads'] / stats['total_reads'] * 100) if stats['total_reads'] > 0 else 0
    }
    
    if stats['insert_sizes']:
        metrics.update({
            'average_insert_size': sum(stats['insert_sizes']) / len(stats['insert_sizes']),
            'median_insert_size': sorted(stats['insert_sizes'])[len(stats['insert_sizes'])//2],
            'mode_insert_size': Counter(stats['insert_sizes']).most_common(1)[0][0]
        })
    
    return metrics

def main():
    if len(sys.argv) != 8:
        print("Usage: collect_metrics.py <bam_file> <metrics_file> <mito_stats> <exome_stats> <onco_stats> <meta_stats> <output_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    metrics_file = sys.argv[2]
    mito_stats = sys.argv[3]
    exome_stats = sys.argv[4]
    onco_stats = sys.argv[5]
    meta_stats = sys.argv[6]
    output_file = sys.argv[7]
    
    # Собираем базовые метрики из BAM файла
    stats = get_bam_stats(bam_file)
    metrics = calculate_metrics(stats)
    
    # Читаем дополнительные метрики из файлов
    with open(mito_stats, 'r') as f:
        metrics['mitochondrial_coverage'] = float(f.readline().strip().split('\t')[2])
    
    with open(exome_stats, 'r') as f:
        metrics['exome_coverage'] = float(f.readline().strip().split('\t')[2])
    
    with open(onco_stats, 'r') as f:
        metrics['onco_panel_coverage'] = float(f.readline().strip().split('\t')[2])
    
    with open(meta_stats, 'r') as f:
        metrics['metagenome_coverage'] = float(f.readline().strip().split('\t')[2])
    
    # Записываем все метрики в выходной файл
    with open(output_file, 'w') as f:
        for key, value in metrics.items():
            f.write(f"{key}\t{value}\n")

if __name__ == "__main__":
    main() 