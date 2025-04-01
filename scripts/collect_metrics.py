#!/usr/bin/env python3
import sys
import os
import pysam
from pathlib import Path
import numpy as np
from collections import Counter
import json
import pandas as pd

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

def parse_fastqc_data(fastqc_data_file):
    """
    Парсит файл fastqc_data.txt и извлекает метрики.
    """
    metrics = {}
    current_section = None
    
    try:
        with open(fastqc_data_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>>'):
                    if line.startswith('>>END_MODULE'):
                        current_section = None
                    else:
                        current_section = line[2:].split('\t')[0]
                        continue
                
                if current_section == 'Basic Statistics':
                    if ':' in line:
                        key, value = line.split(':', 1)
                        key = key.strip()
                        value = value.strip()
                        
                        # Преобразуем числовые значения
                        if value.replace('.', '').replace(',', '').isdigit():
                            value = float(value.replace(',', ''))
                        elif value.endswith('%'):
                            value = float(value.rstrip('%')) / 100
                        elif value.endswith(' bp'):
                            value = int(value.split()[0])
                        
                        metrics[key] = value
                        
    except Exception as e:
        print(f"Warning: Failed to parse FastQC data: {str(e)}", file=sys.stderr)
    
    return metrics

def parse_bam_stats(bam_stats_file):
    """
    Парсит файл статистики BAM и извлекает метрики.
    """
    metrics = {}
    try:
        with open(bam_stats_file, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Преобразуем числовые значения
                    if value.replace('.', '').replace(',', '').isdigit():
                        value = float(value.replace(',', ''))
                    elif value.endswith('%'):
                        value = float(value.rstrip('%')) / 100
                    
                    metrics[key] = value
                    
    except Exception as e:
        print(f"Warning: Failed to parse BAM stats: {str(e)}", file=sys.stderr)
    
    return metrics

def collect_metrics(fastqc_data_file, bam_stats_file, output_file):
    """
    Собирает метрики из FastQC и BAM статистики.
    """
    try:
        # Парсим метрики из FastQC
        fastqc_metrics = parse_fastqc_data(fastqc_data_file)
        
        # Парсим метрики из BAM статистики
        bam_metrics = parse_bam_stats(bam_stats_file)
        
        # Объединяем метрики
        metrics = {
            'fastqc': fastqc_metrics,
            'bam': bam_metrics
        }
        
        # Сохраняем метрики в JSON
        with open(output_file, 'w') as f:
            json.dump(metrics, f, indent=2)
            
    except Exception as e:
        print(f"Error collecting metrics: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) != 4:
        print("Usage: collect_metrics.py <fastqc_data_file> <bam_stats_file> <output_file>")
        sys.exit(1)
    
    fastqc_data_file = sys.argv[1]
    bam_stats_file = sys.argv[2]
    output_file = sys.argv[3]
    
    collect_metrics(fastqc_data_file, bam_stats_file, output_file)

if __name__ == "__main__":
    main() 