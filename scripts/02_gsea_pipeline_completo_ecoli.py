import pandas as pd
import gseapy as gp
import numpy as np
import os
import sys

# --- CONFIGURAÇÕES GERAIS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# Caminhos Fixos
ARQUIVO_ANOTACAO = os.path.join(BASE_DIR, "data", "external_data", "GCF_026637775.1_ASM2663777v1_feature_table.txt.gz")
ARQUIVO_GMT = os.path.join(BASE_DIR, "data", "external_data", "ecoli_kegg_oficial.gmt")

# --- LISTA DE TAREFAS (O que vamos rodar) ---
TAREFAS = [
    {
        "nome": "Adaptação Inicial (D7 vs D1)",
        "input": "res_D7_vs_D1.tabular",
        "output_rnk": "D7_vs_D1.rnk",
        "output_folder": "gsea_D7_vs_D1"
    },
    {
        "nome": "Fase Tardia (D15 vs D1)",
        "input": "res_D15_vs_D1.tabular",
        "output_rnk": "D15_vs_D1.rnk",
        "output_folder": "gsea_D15_vs_D1"
    },
    {
        "nome": "Transição (D15 vs D7)",
        "input": "res_D15_vs_D7.tabular",
        "output_rnk": "D15_vs_D7.rnk",
        "output_folder": "gsea_D15_vs_D7"
    }
]

def processar_deseq2(arquivo_input, arquivo_output_rnk):
    """Lê o tabular do Galaxy, converte IDs e gera RNK"""
    caminho_input = os.path.join(BASE_DIR, "data", "raw_data", arquivo_input)
    caminho_output = os.path.join(BASE_DIR, "data", "processed_data", arquivo_output_rnk)
    
    if not os.path.exists(caminho_input):
        print(f"   ⚠️ ARQUIVO NÃO ENCONTRADO: {arquivo_input}")
        return None

    # 1. Carregar Anotação
    df_annot = pd.read_csv(ARQUIVO_ANOTACAO, sep='\t', compression='gzip')
    df_annot = df_annot[['locus_tag', 'symbol']].dropna().drop_duplicates(subset='locus_tag')

    # 2. Carregar DESeq2 (Galaxy sem cabeçalho)
    colunas_padrao = ['locus_tag', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    df_de = pd.read_csv(caminho_input, sep='\t', header=None, names=colunas_padrao)
    df_de['locus_tag'] = df_de['locus_tag'].astype(str).str.strip()

    # 3. Calcular Ranking
    # Se stat for nulo, recalcula usando logFC e pvalue
    if df_de['stat'].isnull().any():
        df_de['pvalue'] = df_de['pvalue'].fillna(1)
        df_de['log2FoldChange'] = df_de['log2FoldChange'].fillna(0)
        df_de['stat'] = np.sign(df_de['log2FoldChange']) * -np.log10(df_de['pvalue'] + 1e-300)
    
    df_rank = df_de[['locus_tag', 'stat']].copy()

    # 4. Cruzar IDs
    df_merged = df_rank.merge(df_annot, on='locus_tag', how='inner')
    if len(df_merged) == 0:
        print("   ❌ Erro: Nenhum gene convertido.")
        return None

    # 5. Salvar RNK
    os.makedirs(os.path.dirname(caminho_output), exist_ok=True)
    df_final = df_merged.groupby('symbol')['stat'].mean().reset_index()
    df_final = df_final.sort_values(by='stat', ascending=False)
    df_final.to_csv(caminho_output, sep='\t', index=False, header=False)
    
    return caminho_output

def rodar_gsea_task(arquivo_rnk, pasta_saida):
    """Roda o GSEApy para um arquivo RNK"""
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
                   verbose=False) # Verbose False para limpar o terminal
        return True
    except Exception as e:
        print(f"   ❌ Erro no GSEA: {e}")
        return False

# --- EXECUÇÃO PRINCIPAL ---
print("="*50)
print("🚀 INICIANDO PIPELINE COMPLETO (3 COMPARAÇÕES)")
print("="*50)

if not os.path.exists(ARQUIVO_GMT):
    print("❌ ERRO CRÍTICO: Biblioteca KEGG não encontrada.")
    print("Rode o script 00_gerar_gmt_kegg_v3.py primeiro!")
    sys.exit()

for tarefa in TAREFAS:
    print(f"\n🔹 Processando: {tarefa['nome']} ...")
    
    # 1. Converter
    rnk_gerado = processar_deseq2(tarefa['input'], tarefa['output_rnk'])
    
    if rnk_gerado:
        print(f"   ✅ Conversão OK! ({tarefa['output_rnk']})")
        
        # 2. Rodar GSEA
        print(f"   ⏳ Rodando GSEA...")
        sucesso = rodar_gsea_task(rnk_gerado, tarefa['output_folder'])
        
        if sucesso:
            print(f"   🎉 GSEA Concluído! Resultados em: results/{tarefa['output_folder']}")
        else:
            print("   ⚠️ Falha no GSEA.")
    else:
        print("   ⚠️ Pulei o GSEA pois a conversão falhou.")

print("\n" + "="*50)
print("🏁 TUDO PRONTO! Verifique a pasta 'results'.")