#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import gffutils
from pathlib import Path

def get_oncogenes_from_gff3(gff3_file):
    """Извлекает список онкогенов из GFF3 файла."""
    oncogenes = set()
    try:
        # Создаем базу данных GFF3
        db = gffutils.create_db(gff3_file, ':memory:', force=True)
        
        # Список известных онкогенов
        known_oncogenes = {
            'BRCA1', 'BRCA2', 'TP53', 'EGFR', 'KRAS', 'NRAS', 'HRAS',
            'PIK3CA', 'PTEN', 'AKT1', 'BRAF', 'MYC', 'ERBB2', 'MET',
            'ALK', 'ROS1', 'RET', 'NTRK1', 'NTRK2', 'NTRK3'
        }
        
        # Ищем гены в GFF3
        for gene in db.features_of_type('gene'):
            gene_name = gene.attributes.get('Name', [''])[0]
            if gene_name in known_oncogenes:
                oncogenes.add(gene_name)
                
    except Exception as e:
        print(f"WARNING: Failed to parse GFF3 file: {str(e)}", file=sys.stderr)
    
    return oncogenes

def analyze_oncogene_coverage(coverage_file, oncogenes):
    """Анализирует покрытие онкогенов."""
    try:
        coverage_data = pd.read_csv(coverage_file, sep='\t', header=None)
        if len(coverage_data.columns) < 4:
            return 0.0, 0
        
        # Получаем среднее покрытие для онкогенов
        oncogene_coverage = coverage_data[coverage_data[0].isin(oncogenes)][3].mean()
        oncogene_count = len(coverage_data[coverage_data[0].isin(oncogenes)])
        
        return oncogene_coverage, oncogene_count
    except Exception as e:
        print(f"WARNING: Failed to analyze oncogene coverage: {str(e)}", file=sys.stderr)
        return 0.0, 0

def classify_sample(coverage_file, gff3_file):
    """
    Классифицирует тип образца на основе покрытия и наличия онкогенов.
    Возвращает один из типов:
    - mitochondrial
    - exome
    - onco_panel
    - metagenome
    """
    try:
        # Получаем список онкогенов из GFF3
        oncogenes = get_oncogenes_from_gff3(gff3_file)
        
        # Читаем файл покрытия
        coverage_data = pd.read_csv(coverage_file, sep='\t', header=None)
        if len(coverage_data.columns) < 4:
            return "unknown"
        
        # Получаем среднее покрытие
        mean_coverage = coverage_data[3].mean()
        
        # Анализируем покрытие онкогенов
        oncogene_coverage, oncogene_count = analyze_oncogene_coverage(coverage_file, oncogenes)
        
        # Классифицируем на основе покрытия и онкогенов
        if mean_coverage > 5000:  # Митохондриальный геном обычно имеет очень высокое покрытие (>5000x)
            return "mitochondrial"
        elif oncogene_coverage > 50 and oncogene_count >= 3:  # Онко-панели имеют хорошее покрытие онкогенов
            return "onco_panel"
        elif mean_coverage > 100:  # Экзомное секвенирование обычно имеет покрытие >100x
            return "exome"
        else:  # Метагеномное секвенирование обычно имеет низкое покрытие
            return "metagenome"
            
    except Exception as e:
        print(f"Error classifying sample: {str(e)}", file=sys.stderr)
        return "unknown"

def main():
    # Получаем параметры из Snakemake
    coverage_file = snakemake.input.coverage
    gff3_file = snakemake.input.annotation
    output_file = snakemake.output.flag
    
    # Классифицируем образец
    sample_type = classify_sample(coverage_file, gff3_file)
    
    # Записываем результат
    with open(output_file, 'w') as f:
        f.write(sample_type)

if __name__ == "__main__":
    main() 