import requests
import os

# --- MUDANÇA AQUI ---
ORGANISMO = "ppu" # Código KEGG para Pseudomonas putida kt2440
# --------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
PASTA_DESTINO = os.path.join(BASE_DIR, "data", "external_data")
ARQUIVO_SAIDA = os.path.join(PASTA_DESTINO, "pseudomonas_kegg.gmt") # Nome novo

print(f"--- GERADOR GMT: PSEUDOMONAS PUTIDA ({ORGANISMO}) ---")

try:
    os.makedirs(PASTA_DESTINO, exist_ok=True)

    # 1. Baixar Genes
    print("> Baixando genes P. putida...")
    url_genes = f"http://rest.kegg.jp/list/{ORGANISMO}"
    resp_genes = requests.get(url_genes)
    gene_map = {}
    
    for linha in resp_genes.text.strip().split('\n'):
        if not linha: continue
        colunas = linha.split('\t')
        kegg_id = colunas[0] # ppu:PP_0001
        
        # Lógica de nome (igual ao script V3)
        info = ""
        for col in colunas:
            if ";" in col:
                info = col
                break
        if not info and len(colunas) > 1: info = colunas[-1]
        
        parts = info.split(';')
        symbol = parts[0].strip()
        
        # Limpeza
        if symbol in ["CDS", "tRNA", "rRNA", "Gene"] or symbol.startswith(f"{ORGANISMO}:"):
             symbol = kegg_id.replace(f"{ORGANISMO}:", "")
             
        gene_map[kegg_id] = symbol

    # 2. Baixar Vias
    print("> Baixando vias...")
    url_link = f"http://rest.kegg.jp/link/pathway/{ORGANISMO}"
    resp_link = requests.get(url_link)
    pathway_genes = {}
    
    for linha in resp_link.text.strip().split('\n'):
        cols = linha.split('\t')
        if len(cols) < 2: continue
        gene_id = cols[0]
        pathway_id = cols[1]
        gene_symbol = gene_map.get(gene_id, gene_id.replace(f"{ORGANISMO}:", ""))
        
        if pathway_id not in pathway_genes: pathway_genes[pathway_id] = []
        pathway_genes[pathway_id].append(gene_symbol)

    # 3. Nomes das Vias
    print("> Baixando nomes das vias...")
    url_names = f"http://rest.kegg.jp/list/pathway/{ORGANISMO}"
    resp_names = requests.get(url_names)
    pathway_names = {line.split('\t')[0]: line.split('\t')[1] for line in resp_names.text.strip().split('\n')}

    # 4. Salvar
    print(f"> Salvando {os.path.basename(ARQUIVO_SAIDA)}...")
    with open(ARQUIVO_SAIDA, 'w') as f:
        for path_id, genes in pathway_genes.items():
            if len(genes) < 3: continue
            nome_via = pathway_names.get(path_id, path_id)
            f.write(f"{nome_via}\t{path_id}\t" + "\t".join(genes) + "\n")

    print("\n✅ SUCESSO! Banco da Pseudomonas criado.")

except Exception as e:
    print(f"❌ Erro: {e}")