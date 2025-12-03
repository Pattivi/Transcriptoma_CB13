import requests
import os
import sys

# --- CONFIGURAÇÕES ---
ORGANISMO = "eco" # E. coli K-12
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
PASTA_DESTINO = os.path.join(BASE_DIR, "data", "external_data")
ARQUIVO_SAIDA = os.path.join(PASTA_DESTINO, "ecoli_kegg_oficial.gmt")

print(f"--- GERADOR DE GMT KEGG V3 (Correção de Colunas) ---")

try:
    os.makedirs(PASTA_DESTINO, exist_ok=True)

    # 1. Baixar dicionário de genes
    print("> Baixando lista de genes do KEGG...")
    url_genes = f"http://rest.kegg.jp/list/{ORGANISMO}"
    resp_genes = requests.get(url_genes)
    resp_genes.raise_for_status()

    gene_map = {}
    linhas = resp_genes.text.strip().split('\n')
    
    # Diagnóstico visual da primeira linha
    print(f"   Formato detectado: {linhas[0]}")

    for linha in linhas:
        if not linha: continue
        colunas = linha.split('\t')
        
        # Lógica inteligente para achar o nome
        # O KEGG pode mandar 2 ou 4 colunas. O nome sempre está na última que contém ';'
        # Exemplo: "eco:b0001  CDS  190..255  thrL; thr operon leader..."
        
        kegg_id = colunas[0] # eco:b0001
        
        # Procura a coluna que tem a descrição (contém ';')
        info = ""
        for col in colunas:
            if ";" in col:
                info = col
                break
        
        # Se não achou ';' em lugar nenhum, pega a última coluna disponível
        if not info and len(colunas) > 1:
            info = colunas[-1]

        # Extrai o símbolo (primeira palavra antes do ;)
        parts = info.split(';')
        symbol = parts[0].strip()
        
        # Limpeza final (caso ainda tenha sobrado lixo)
        if symbol in ["CDS", "tRNA", "rRNA", "Gene"] or symbol.startswith("eco:"):
             # Tenta usar o ID limpo (b0001) se não tiver nome real
             symbol = kegg_id.replace("eco:", "")

        gene_map[kegg_id] = symbol

    print(f"   Mapping criado. Exemplo: ('eco:b0001' -> '{gene_map.get('eco:b0001')}')")

    # 2. Baixar Vias (Igual ao anterior)
    print("> Baixando estrutura das vias...")
    url_link = f"http://rest.kegg.jp/link/pathway/{ORGANISMO}"
    resp_link = requests.get(url_link)
    
    pathway_genes = {}
    for linha in resp_link.text.strip().split('\n'):
        cols = linha.split('\t')
        if len(cols) < 2: continue
        
        gene_id = cols[0]
        pathway_id = cols[1]
        
        # Converte ID -> Símbolo (agora correto)
        gene_symbol = gene_map.get(gene_id, gene_id.replace("eco:", ""))
        
        # Ignora se ainda for "CDS" (improvável agora)
        if gene_symbol == "CDS": continue

        if pathway_id not in pathway_genes:
            pathway_genes[pathway_id] = []
        pathway_genes[pathway_id].append(gene_symbol)

    # 3. Baixar Nomes das Vias
    print("> Baixando nomes das vias...")
    url_names = f"http://rest.kegg.jp/list/pathway/{ORGANISMO}"
    resp_names = requests.get(url_names)
    pathway_names = {line.split('\t')[0]: line.split('\t')[1] for line in resp_names.text.strip().split('\n')}

    # 4. Salvar
    print(f"> Salvando {os.path.basename(ARQUIVO_SAIDA)}...")
    with open(ARQUIVO_SAIDA, 'w') as f:
        count = 0
        for path_id, genes in pathway_genes.items():
            if len(genes) < 3: continue
            nome_via = pathway_names.get(path_id, path_id)
            linha = f"{nome_via}\t{path_id}\t" + "\t".join(genes) + "\n"
            f.write(linha)
            count += 1

    print(f"\n✅ SUCESSO V3! Biblioteca gerada corretamente.")
    print(f"Total de vias: {count}")
    print("Pode rodar o GSEA agora!")

except Exception as e:
    print(f"❌ Erro: {e}")