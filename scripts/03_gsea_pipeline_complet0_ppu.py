import pandas as pd
import gseapy as gp
import numpy as np
import os
import sys

# --- CONFIGURAÇÕES ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# 1. ARQUIVO GMT (Apontando para Pseudomonas)
ARQUIVO_GMT = os.path.join(BASE_DIR, "data", "external_data", "pseudomonas_kegg.gmt")

# 2. ARQUIVO DE ANOTAÇÃO (Enterobacter)
ARQUIVO_ANOTACAO = os.path.join(BASE_DIR, "data", "external_data", "GCF_026637775.1_ASM2663777v1_feature_table.txt.gz")

# --- LISTA DE TAREFAS (Note os sufixos _PPU nas pastas de saída) ---
TAREFAS = [
    {
        "nome": "Adaptação Inicial (D7 vs D1) - PPU",
        "input": "res_D7_vs_D1.tabular",
        "output_rnk": "D7_vs_D1.rnk", 
        "output_folder": "gsea_D7_vs_D1_PPU" # <--- Pasta nova!
    },
    {
        "nome": "Fase Tardia (D15 vs D1) - PPU",
        "input": "res_D15_vs_D1.tabular",
        "output_rnk": "D15_vs_D1.rnk",
        "output_folder": "gsea_D15_vs_D1_PPU" # <--- Pasta nova!
    },
    {
        "nome": "Transição (D15 vs D7) - PPU",
        "input": "res_D15_vs_D7.tabular",
        "output_rnk": "D15_vs_D7.rnk",
        "output_folder": "gsea_D15_vs_D7_PPU" # <--- Pasta nova!
    }
]

def processar_deseq2(arquivo_input, arquivo_output_rnk):
    """Lê o tabular, converte e gera RNK"""
    caminho_input = os.path.join(BASE_DIR, "data", "raw_data", arquivo_input)
    caminho_output = os.path.join(BASE_DIR, "data", "processed_data", arquivo_output_rnk)
    
    if not os.path.exists(caminho_input):
        print(f"   ⚠️ Input não achado: {arquivo_input}")
        return None

    # Carrega anotação
    df_annot = pd.read_csv(ARQUIVO_ANOTACAO, sep='\t', compression='gzip')
    df_annot = df_annot[['locus_tag', 'symbol']].dropna().drop_duplicates(subset='locus_tag')

    # Carrega DESeq2
    colunas_padrao = ['locus_tag', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    df_de = pd.read_csv(caminho_input, sep='\t', header=None, names=colunas_padrao)
    df_de['locus_tag'] = df_de['locus_tag'].astype(str).str.strip()

    # Calcula Ranking
    if df_de['stat'].isnull().any():
        df_de['pvalue'] = df_de['pvalue'].fillna(1)
        df_de['log2FoldChange'] = df_de['log2FoldChange'].fillna(0)
        df_de['stat'] = np.sign(df_de['log2FoldChange']) * -np.log10(df_de['pvalue'] + 1e-300)
    
    df_rank = df_de[['locus_tag', 'stat']].copy()
    
    # Cruza IDs
    df_merged = df_rank.merge(df_annot, on='locus_tag', how='inner')
    
    if len(df_merged) == 0: return None

    # Salva RNK
    os.makedirs(os.path.dirname(caminho_output), exist_ok=True)
    df_final = df_merged.groupby('symbol')['stat'].mean().reset_index()
    df_final = df_final.sort_values(by='stat', ascending=False)
    df_final.to_csv(caminho_output, sep='\t', index=False, header=False)
    
    return caminho_output

def rodar_gsea_task(arquivo_rnk, pasta_saida):
    caminho_saida = os.path.join(BASE_DIR, "results", pasta_saida)
    try:
        gp.prerank(rnk=arquivo_rnk,
                   gene_sets=ARQUIVO_GMT,
                   threads=4,
                   min_size=5,
                   max_size=1000,
                   permutation_num=1000,
                   outdir=caminho_saida,
                   seed=42,
                   verbose=False)
        return True
    except Exception as e:
        print(f"   ❌ Erro GSEA: {e}")
        return False

# --- EXECUÇÃO ---
print("="*50)
print("🚀 PIPELINE PSEUDOMONAS PUTIDA (PPU)")
print("="*50)

if not os.path.exists(ARQUIVO_GMT):
    print("❌ ERRO: Arquivo 'pseudomonas_kegg.gmt' não encontrado.")
    print("Rode o script 00_gerar_gmt_pseudomonas.py primeiro!")
    sys.exit()

for tarefa in TAREFAS:
    print(f"\n🔹 {tarefa['nome']} ...")
    rnk = processar_deseq2(tarefa['input'], tarefa['output_rnk'])
    if rnk:
        print(f"   ⏳ Rodando GSEA...")
        if rodar_gsea_task(rnk, tarefa['output_folder']):
            print(f"   🎉 Sucesso! Pasta: {tarefa['output_folder']}")
    else:
        print("   ⚠️ Erro na conversão.")

print("\n🏁 FIM.")