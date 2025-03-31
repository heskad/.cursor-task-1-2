#!/usr/bin/env python3
import sys
import os
import pysam
from pathlib import Path
import zipfile
import numpy as np
import yaml

def parse_fastqc_data(zip_path):
    """Извлекает метрики из fastqc_data.txt в ZIP-архиве FastQC."""
    metrics = {
        "total_sequences": 0,
        "sequence_length": 0,
        "adapter_percent": 0.0,
        "overlapping_percent": 0.0
    }
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            # Определяем имя папки внутри архива (формат: {sample}_read{1|2}_fastqc)
            base_name = Path(zip_path).stem
            data_file = f"{base_name}/fastqc_data.txt"
            
            with z.open(data_file) as f:
                data = f.read().decode('utf-8')
                
                # Парсим основные метрики
                for line in data.split('\n'):
                    if 'Total Sequences' in line:
                        metrics["total_sequences"] = int(line.split('\t')[-1])
                    elif 'Sequence length' in line:
                        metrics["sequence_length"] = int(line.split('\t')[-1])
                    elif '>>Adapter Content' in line:
                        # Берем максимальный процент адаптеров по всем позициям
                        adapter_section = data.split('>>Adapter Content')[1].split('>>END_MODULE')[0]
                        adapter_lines = [l.split('\t') for l in adapter_section.split('\n')[2:-1] if l]
                        max_percent = max(float(pos) for line in adapter_lines for pos in line[1:] if pos)
                        metrics["adapter_percent"] = max_percent
                    elif '>>Overrepresented sequences' in line:
                        # Оценка перекрывающихся ридов (примерная)
                        overrep_section = data.split('>>Overrepresented sequences')[1].split('>>END_MODULE')[0]
                        metrics["overlapping_percent"] = len(overrep_section.split('\n')) * 0.1  # Примерная оценка
                        
    except Exception as e:
        print(f"WARNING: FastQC parsing failed: {str(e)}", file=sys.stderr)
    
    return metrics

def get_insert_size_stats(bam_path, sample_size=10000):
    """Собирает статистики размера вставки (выборка для скорости)."""
    insert_sizes = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, read in enumerate(bam):
            if i >= sample_size * 2:  # Проверяем только часть данных
                break
            if read.is_proper_pair and read.is_read1 and read.template_length > 0:
                insert_sizes.append(abs(read.template_length))
    
    if not insert_sizes:
        return (0, 0, 0)
    
    return (
        round(np.mean(insert_sizes)),
        int(np.median(insert_sizes)),
        int(max(set(insert_sizes), key=insert_sizes.count))
    )

def get_mapping_stats(bam_path):
    """Собирает статистики выравнивания."""
    total = mapped = non_unique = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            total += 1
            if not read.is_unmapped:
                mapped += 1
            if read.has_tag('NH') and read.get_tag('NH') > 1:
                non_unique += 1
    return total, mapped, non_unique

def parse_picard_metrics(metrics_path):
    """Извлекает процент дубликатов из файла Picard."""
    with open(metrics_path) as f:
        for line in f:
            if line.startswith('PERCENT_DUPLICATION'):
                return float(line.strip().split('\t')[1]) * 100
    return 0.0

def main():
    if len(sys.argv) != 4:
        raise ValueError("Usage: generate_report.py <input.bam> <metrics.txt> <output.report>")
    
    bam_path = sys.argv[1]
    metrics_path = sys.argv[2]
    output_path = sys.argv[3]
    
    # Проверка файлов
    if not Path(bam_path).exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")
    if not Path(metrics_path).exists():
        raise FileNotFoundError(f"Metrics file not found: {metrics_path}")

    # Получаем имя образца из пути к BAM
    sample_name = Path(bam_path).stem.split('.')[0]
    
    # Путь к FastQC (предполагаем стандартную структуру папок)
    fastqc_zip = os.path.join(
        os.path.dirname(bam_path),
        f"../../metrics/raw_fastqc/{sample_name}/{sample_name}_read1_fastqc.zip"
    )

    # Собираем все метрики
    try:
        # 1. FastQC метрики
        fastqc_data = parse_fastqc_data(fastqc_zip) if Path(fastqc_zip).exists() else {}
        
        # 2. Статистики выравнивания
        total_reads, mapped_reads, non_unique = get_mapping_stats(bam_path)
        
        # 3. Размеры вставок
        mean_insert, median_insert, mode_insert = get_insert_size_stats(bam_path)
        
        # 4. Метрики дубликатов
        dup_percent = parse_picard_metrics(metrics_path)
        
        # Генерация отчета
        with open(output_path, 'w') as report:
            report.write(f"Sample: {sample_name}\n")
            report.write("=" * 50 + "\n")
            report.write(f"Total reads: {total_reads:,}\n")
            report.write(f"Read length: {fastqc_data.get('sequence_length', 'N/A')} bp\n")
            report.write(f"Mapped reads: {mapped_reads:,} ({mapped_reads/total_reads*100:.2f}%)\n")
            report.write(f"Adapter content: {fastqc_data.get('adapter_percent', 0.0):.2f}%\n")
            report.write(f"Duplicate rate: {dup_percent:.2f}%\n")
            report.write(f"Non-unique alignments: {non_unique:,} ({non_unique/total_reads*100:.2f}%)\n")
            report.write(f"Overlapping reads (est.): {fastqc_data.get('overlapping_percent', 0.0):.2f}%\n")
            report.write(f"Insert size - Mean: {mean_insert}, Median: {median_insert}, Mode: {mode_insert}\n\n")
            
            # Добавляем сырые данные из Picard
            report.write("Picard MarkDuplicates Metrics:\n")
            with open(metrics_path) as f:
                report.write(f.read())
        
        # Добавляем тип исследования (один блок)
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
