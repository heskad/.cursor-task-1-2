#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import gffutils
from pathlib import Path

def get_oncogenes_from_gff3(gff3_file):
    """
    Извлекает список онкогенов из GFF3 файла.
    """
    # Список известных онкогенов
    known_oncogenes = {
        'BRCA1', 'BRCA2', 'TP53', 'EGFR', 'KRAS', 'NRAS', 'BRAF', 'PIK3CA',
        'PTEN', 'AKT1', 'MYC', 'ERBB2', 'MET', 'ALK', 'ROS1', 'RET',
        'PDGFRA', 'KIT', 'FLT3', 'JAK2', 'BCR', 'ABL1', 'NOTCH1', 'CTNNB1',
        'SMAD4', 'CDKN2A', 'RB1', 'VHL', 'NF1', 'NF2', 'TSC1', 'TSC2'
    }
    
    try:
        # Создаем базу данных GFF3 в памяти с отключенной проверкой дубликатов
        db = gffutils.create_db(
            gff3_file,
            ':memory:',
            force=True,
            keep_order=True,
            sort_attribute_values=True,
            merge_strategy='create_unique'
        )
        
        # Ищем гены, которые могут быть онкогенами
        oncogenes = set()
        for gene in db.features_of_type('gene'):
            gene_name = gene.attributes.get('gene_name', [''])[0]
            if gene_name in known_oncogenes:
                oncogenes.add(gene_name)
        
        return list(oncogenes)
        
    except Exception as e:
        print(f"Warning: Failed to parse GFF3 file: {str(e)}", file=sys.stderr)
        return list(known_oncogenes)  # Возвращаем базовый список онкогенов в случае ошибки

def analyze_oncogene_coverage(coverage_file, oncogenes):
    """
    Анализирует покрытие онкогенов в файле покрытия.
    """
    try:
        # Читаем файл покрытия
        coverage_data = pd.read_csv(coverage_file, sep='\t', header=None)
        if len(coverage_data.columns) < 4:
            return 0, 0
        
        # Получаем среднее покрытие и количество генов с хорошим покрытием
        mean_coverage = coverage_data[3].mean()
        genes_with_good_coverage = sum(1 for x in coverage_data[3] if x > 30)
        
        return mean_coverage, genes_with_good_coverage
        
    except Exception as e:
        print(f"Warning: Failed to analyze coverage: {str(e)}", file=sys.stderr)
        return 0, 0

def classify_sample(coverage_file, gff3_file):
    """
    Классифицирует тип образца на основе покрытия и наличия онкогенов.
    """
    try:
        # Получаем список онкогенов
        oncogenes = get_oncogenes_from_gff3(gff3_file)
        
        # Анализируем покрытие
        mean_coverage, genes_with_good_coverage = analyze_oncogene_coverage(coverage_file, oncogenes)
        
        # Классифицируем на основе покрытия и онкогенов
        if mean_coverage > 5000:  # Митохондриальный геном обычно имеет очень высокое покрытие
            return "mitochondrial"
        elif mean_coverage > 100:  # Экзомное секвенирование обычно имеет покрытие >100x
            return "exome"
        elif mean_coverage > 50 and genes_with_good_coverage > 5:  # Онко-панели обычно имеют покрытие 50-100x и содержат несколько онкогенов
            return "onco_panel"
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