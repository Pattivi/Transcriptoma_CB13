import gseapy as gp
import os
import matplotlib.pyplot as plt

# --- CAMINHOS DINÂMICOS ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

# Entrada: O arquivo .rnk que o script 01 gerou
ARQUIVO_RNK = os.path.join(BASE_DIR, "data", "processed_data", "lista_pronta_para_gsea.rnk")

# Saída: Pasta results na raiz
PASTA_SAIDA = os.path.join(BASE_DIR, "results", "gsea_analysis")

# Biblioteca (Enterobacter -> E. coli K12 model)
BIBLIOTECA = "KEGG_2019_Escherichia_coli_K-12_MG1655"

print(f"--- Rodando GSEA ---")
print(f"Input: {ARQUIVO_RNK}")
print(f"Output: {PASTA_SAIDA}")

if not os.path.exists(ARQUIVO_RNK):
    print(f"❌ ERRO: Arquivo .rnk não encontrado em processed_data. Rode o script 01 primeiro!")
    exit()

try:
    # Garante que a pasta results existe
    os.makedirs(PASTA_SAIDA, exist_ok=True)

    pre_res = gp.prerank(rnk=ARQUIVO_RNK,
                         gene_sets=BIBLIOTECA,
                         threads=4,
                         min_size=5,
                         max_size=1000,
                         permutation_num=1000,
                         outdir=PASTA_SAIDA,
                         seed=42,
                         verbose=True)

    print(f"\n✅ Análise concluída! Verifique a pasta 'results'.")

except Exception as e:
    print(f"Erro: {e}")