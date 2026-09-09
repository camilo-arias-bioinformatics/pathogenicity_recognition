# %% Importar paquetes necesarios
import pandas as pd
from Bio import SeqIO

# %% Exploración inicial de la tabla de datos
data = pd.read_csv("data/raw/clinvar_result.txt", sep="\t")

for col in data.columns:
    print(f"\n{col}")
    print(data[col].unique())

# %% Generar tabla con registros patogénicos y por variación de un único nucleótido
data_sorted = data[data["Germline classification"] == "Pathogenic"]
data_sorted = data_sorted[data_sorted["Variant type"] == "single nucleotide variant"]
print(data_sorted.info())

# %% Eliminar registros con los que no se pueda trabajar
data_sorted[data_sorted["Gene(s)"].isna()]
data_sorted = data_sorted.dropna(subset=["Gene(s)"]).reset_index(drop=True)

# %% Exportar la tabla generada como un archivo csv
data_sorted.to_csv("data/processed/hcm_pathogenic_and_snp.csv", index=False)

# %% Conservar registros únicos para genes desde el nombre de la mutación
genes = []

for a in data_sorted["Name"]:
    end = a.find("(")
    genes.append(a[0:end])

genes = list(set(genes))
print(genes)
len(genes)

# %%Función para descargar secuencias de una lista de genes a un archivo fasta
def buscar_secuencias(accessions):
    from Bio import Entrez

    Entrez.email = "j.camilo.ariasospina@gmail.com"

    with open("results/gene_sequences.fasta", "w") as outfile:

        for acc in accessions:

            fetch = Entrez.efetch(
                db="nucleotide",
                id=acc,
                rettype="fasta",
                retmode="text"
            )

            outfile.write(fetch.read())
            outfile.write("\n")

            print(f"Downloaded {acc}")
buscar_secuencias(genes)

# %% Revisar el contenido del archivo fasta
fasta_file = "results/gene_sequences.fasta"

records = list(SeqIO.parse(fasta_file, "fasta"))

print(f"Number of sequences: {len(records)}\n")

for i, record in enumerate(records, start=1):
    print(f"{i:>4} | {record.id} | length = {len(record.seq)}")

# %% Importar las secuencias con One Hot Encoding

from one_hot_encoding import create_dataset_from_fasta

ids, dataset = create_dataset_from_fasta("results/gene_sequences.fasta")

# %% Crear una tabla de datos que contenga las etiquetas y secuencias

nucleotide_numbers = {"A":0, "C":1, "G":2, "T":3, "U":4}
label = [1] * (len(data_sorted))
gene = []
location = []
mutation_from = []
mutation_to = []

for a in data_sorted["Name"]:
    end = a.find("(")
    gene.append(a[0:end])

    start = a.find(":c.")
    end = a.find(">")
    location.append(a[(start+3):(end-1)])

    position = a.find(">")
    mutation_from.append(nucleotide_numbers[a[(position-1)]])
    mutation_to.append(nucleotide_numbers[a[(position+1)]])

mutations_dataset = pd.DataFrame({"gene": gene, "location": location, "mutation_from": mutation_from, "mutation_to": mutation_to})

print(mutations_dataset.head())

# %% 
label = [0] * len(ids)
dataset_list = []

for a in dataset:
    dataset_list.append(a)    

for a in range(len(mutations_dataset)):
    row = mutations_dataset.iloc[a]
    index = ids.index(row["gene"])
    seq = dataset[index]
    seq[(int(row["mutation_from"])), (int(row["location"])-1)] = False
    seq[(int(row["mutation_to"])), (int(row["location"])-1)] = True
    dataset_list.append(seq)
    ids.append(row["gene"])
    label.append(1)

print(len(ids), lend(label), len(dataset_list))
# %%

row = mutations_dataset.iloc[1]
int(row["mutation_from"])
type(int(row["location"]))

# %%

# %%
for a in ads:
    index = ids.index(a)
    seq = dataset[a]
a = 
print(a)
copy = dataset[a]
for a in range(5):
    print(copy[a, 2]) 

