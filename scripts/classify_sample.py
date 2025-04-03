import pandas as pd

def classify_sequencing(coverage_file, oncogenes):
    df = pd.read_csv(coverage_file, sep='\t', header=None, names=["chr", "start", "end", "reads", "coverage", "length", "normalized_coverage"])
    
    mean_cov = df["normalized_coverage"].mean()
    min_cov = df["normalized_coverage"].min()
    max_cov = df["normalized_coverage"].max()
    high_cov_chromosomes = df[df["normalized_coverage"] > 0.001]
    mitochondrial_cov = df[df["chr"] == "chrM"]["normalized_coverage"].values[0] if "chrM" in df["chr"].values else 0
    
    high_cov_genes = [gene for gene in oncogenes if gene in df["chr"].values and df[df["chr"] == gene]["coverage"].values[0] > 50]
    
    output = []
    output.append(f"Покрытие остальных хромосом в диапазоне {min_cov:.6f} - {max_cov:.6f}, что соответствует экзомному секвенированию или онко-панели.")
    
    if not high_cov_chromosomes.empty:
        for _, row in high_cov_chromosomes.iterrows():
            output.append(f"Высокое покрытие на {row['chr']} ({row['normalized_coverage']:.6f}) — возможно, это целенаправленная панель с онкогенами.")
    
    if mitochondrial_cov < 0.01:
        output.append("Покрытие chrM слишком низкое для митохондриального → это не митохондриальная ДНК.")
    
    output.append("\nВердикт")
    if high_cov_genes:
        output.append("Если в онкогенах есть хорошее покрытие (>50x) → это онко-панель.")
    if mean_cov > 0.0005 and mean_cov < 0.0016:
        output.append("Если покрытие экзонов выше 80% → экзомное секвенирование.")
    if mean_cov < 0.0005:
        output.append("Если много микробного контаминационного ридов → метагеном.")
    
    return "\n".join(output)