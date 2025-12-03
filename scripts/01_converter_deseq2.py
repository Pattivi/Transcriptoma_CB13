import pandas as pd
import numpy as np
import sys
import os

# --- TRUQUE DE CAMINHOS ---
# Pega o diretório onde ESTE script está salvo (.../TRANSCRIPTOMA/scripts)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Pega o diretório raiz do projeto (.../TRANSCRIPTOMA)
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# --- CONFIGURAÇÕES DE ARQUIVOS (CAMINHOS RELATIVOS) ---
# Entrada: Onde você vai salvar o output do DESeq2
ARQUIVO_DESEQ2 = os.path.join(BASE_DIR, "data", "raw_data", "res_deseq2.csv") 

# Entrada: O arquivo do NCBI (já está na pasta certa)
ARQUIVO_ANOTACAO = os.path.join(BASE_DIR, "data", "external_data", "GCF_026637775.1_ASM2663777v1_feature_table.txt.gz")

# Saída: Vamos salvar o processado na pasta processed_data
ARQUIVO_SAIDA = os.path.join(BASE_DIR, "data", "processed_data", "lista_pronta_para_gsea.rnk")

# Configurações do arquivo (MANTENHA COMO ESTAVA)
SEPARADOR_DESEQ2 = ',' # Ou '\t' se for tabulação

print("--- INICIANDO CONVERSÃO (ESTRUTURA ORGANIZADA) ---")
print(f"Diretório base do projeto: {BASE_DIR}")

try:
    # Verifica se a pasta processed_data existe, se não, cria
    os.makedirs(os.path.dirname(ARQUIVO_SAIDA), exist_ok=True)

    # --- PASSO 1: Carregar Anotação ---
    print(f"> Lendo anotação: {os.path.basename(ARQUIVO_ANOTACAO)}")
    if not os.path.exists(ARQUIVO_ANOTACAO):
        raise FileNotFoundError(f"Arquivo não encontrado: {ARQUIVO_ANOTACAO}")

    df_annot = pd.read_csv(ARQUIVO_ANOTACAO, sep='\t', compression='gzip')
    df_annot = df_annot[['locus_tag', 'symbol']].dropna().drop_duplicates(subset='locus_tag')

    # --- PASSO 2: Carregar DESeq2 ---
    print(f"> Lendo DESeq2: {os.path.basename(ARQUIVO_DESEQ2)}")
    if not os.path.exists(ARQUIVO_DESEQ2):
        raise FileNotFoundError(f"Arquivo não encontrado: {ARQUIVO_DESEQ2}\n   -> Certifique-se de ter baixado o arquivo na pasta 'data/raw_data'")

    df_de = pd.read_csv(ARQUIVO_DESEQ2, sep=SEPARADOR_DESEQ2, index_col=0)
    df_de.index.name = 'locus_tag'
    df_de = df_de.reset_index()
    df_de['locus_tag'] = df_de['locus_tag'].astype(str).str.strip()

    # --- PASSO 3: Ranking ---
    if 'stat' in df_de.columns:
        df_rank = df_de[['locus_tag', 'stat']].copy()
    elif 'log2FoldChange' in df_de.columns and 'pvalue' in df_de.columns:
        print("  Calculando Rank manualmente...")
        df_de['pvalue'] = df_de['pvalue'].fillna(1)
        df_de['log2FoldChange'] = df_de['log2FoldChange'].fillna(0)
        df_de['stat'] = np.sign(df_de['log2FoldChange']) * -np.log10(df_de['pvalue'] + 1e-300)
        df_rank = df_de[['locus_tag', 'stat']].copy()
    else:
        raise ValueError("Colunas necessárias não encontradas no DESeq2.")

    # --- PASSO 4: Cruzamento ---
    df_merged = df_rank.merge(df_annot, on='locus_tag', how='inner')
    df_final = df_merged.groupby('symbol')['stat'].mean().reset_index()
    df_final = df_final.sort_values(by='stat', ascending=False)

    # --- PASSO 5: Salvar ---
    df_final.to_csv(ARQUIVO_SAIDA, sep='\t', index=False, header=False)

    print(f"\n✅ SUCESSO! Arquivo salvo em: data/processed_data/{os.path.basename(ARQUIVO_SAIDA)}")

except Exception as e:
    print(f"\n❌ ERRO: {e}")