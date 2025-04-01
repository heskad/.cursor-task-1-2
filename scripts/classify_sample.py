#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np

def classify_sample(coverage_file):
    """
    Классифицирует тип образца на основе покрытия.
    Возвращает один из типов:
    - mitochondrial
    - exome
    - onco_panel
    - metagenome
    """
    try:
        # Читаем файл покрытия
        coverage_data = pd.read_csv(coverage_file, sep='\t', header=None)
        if len(coverage_data.columns) < 4:
            return "unknown"
        
        # Получаем среднее покрытие
        mean_coverage = coverage_data[3].mean()
        
        # Классифицируем на основе покрытия
        if mean_coverage > 1000:  # Митохондриальный геном обычно имеет очень высокое покрытие
            return "mitochondrial"
        elif mean_coverage > 100:  # Экзомное секвенирование обычно имеет покрытие >100x
            return "exome"
        elif mean_coverage > 50:  # Онко-панели обычно имеют покрытие 50-100x
            return "onco_panel"
        else:  # Метагеномное секвенирование обычно имеет низкое покрытие
            return "metagenome"
            
    except Exception as e:
        print(f"Error classifying sample: {str(e)}", file=sys.stderr)
        return "unknown"

def main():
    if len(sys.argv) != 3:
        print("Usage: classify_sample.py <coverage_file> <output_file>")
        sys.exit(1)
    
    coverage_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Классифицируем образец
    sample_type = classify_sample(coverage_file)
    
    # Записываем результат
    with open(output_file, 'w') as f:
        f.write(sample_type)

if __name__ == "__main__":
    main() 