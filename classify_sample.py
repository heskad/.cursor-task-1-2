import sys
import pandas as pd

def classify(coverage_file):
    df = pd.read_csv(coverage_file, sep='\t', names=['chr', 'start', 'end', 'region', 'coverage'])
    
    # Инициализация значений покрытия
    mito_coverage = exome_coverage = oncopanel_coverage = amplicon_coverage = 0.0
    
    if not df.empty:
        # Митохондриальная ДНК
        if 'chrM' in df['chr'].values:
            mito_coverage = df[df['chr'] == 'chrM']['coverage'].mean()
        
        # Экзомные регионы
        if 'exome' in df['region'].values:
            exome_coverage = df[df['region'] == 'exome']['coverage'].mean()
        
        # Онко-панель (пример для BRCA1)
        if 'BRCA1' in df['region'].values:
            oncopanel_coverage = df[df['region'] == 'BRCA1']['coverage'].mean()
        
        # Ампликоны
        if 'amplicon' in df['region'].values:
            amplicon_coverage = df[df['region'] == 'amplicon']['coverage'].mean()
    
    # Логика классификации
    if mito_coverage > 1000:
        return "mitochondrial"
    elif exome_coverage > 50:
        return "exome_panel"
    elif amplicon_coverage > 500:
        return "amplicon_oncopanel"
    elif oncopanel_coverage > 100 and (df['coverage'].std() < 20):
        return "hybridization_oncopanel"
    else:
        return "metagenome"

if __name__ == "__main__":
    coverage_path = sys.argv[1]
    research_type = classify(coverage_path)
    with open(sys.argv[2], 'w') as f:
        f.write(research_type)
