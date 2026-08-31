# %% import necessary packages
import pandas as pd
from Bio import SeqIO

# %% Initial data exploration
# in: clinvar_results.txt (txt file)
# out: column names and data types
data = pd.read_csv("data/raw/clinvar_result.txt", sep="\t")

for col in data.columns:
    print(f"\n{col}")
    print(data[col].unique())

# %% Filter for germline classification and variant type
# in: data table
# out: data_filtered_pathogenic and data_filtered_benign
def filter_data(germline_classification):
    data_filtered = data[data["Germline classification"] == germline_classification]
    data_filtered = data_filtered[data_filtered["Variant type"] == "single nucleotide variant"]
    return(data_filtered)

data_filtered_pathogenic = filter_data("Pathogenic")
data_filtered_benign = filter_data("Benign")
print(data_filtered_pathogenic.info())
print(data_filtered_benign.info())

# %% Check for na values
# in: data_filtered_pathogenic
# out: na details
print(data_filtered_pathogenic[data_filtered_pathogenic["Gene(s)"].isna()])

# %% manage na values
# in: data_filtered_pathogenic with na values
# out: data_filtered_pathogenic without na values
data_filtered_pathogenic = data_filtered_pathogenic.dropna(subset=["Gene(s)"]).reset_index(drop=True)

# %% Export tables as csv files
# in: data tables in code
# out: data tables as csv files
def export_data_as_files(dataframe):
    dataframe.to_csv("data/processed/hcm_pathogenic_and_snp.csv", index=False)
export_data_as_files(data_filtered_pathogenic)
export_data_as_files(data_filtered_benign)

# %% Store unique values for gene names
# in: data_filtered_pathogenic and data_filtered_benign
# out: unique_genes_pathogenic and unique_genes_benign (lists)
def store_gene_names(dataset):
    genes = []

    for a in dataset["Name"]:
        end = a.find("(")
        genes.append(a[0:end])

    genes = list(set(genes))
    return(genes)
unique_genes_pathogenic = store_gene_names(data_filtered_pathogenic)
unique_genes_benign = store_gene_names(data_filtered_benign)
print(f"len pathogenic: {len(unique_genes_pathogenic)} len benign: {len(unique_genes_benign)}")

# %% Download gene sequences to fasta files
# in: unique_genes_pathogenic and unique_genes_benign (lists)
# out: gene_sequences.fasta (fasta file)

sequences_list = unique_genes_pathogenic + unique_genes_benign
sequences_list = list(set(sequences_list))

def sequence_search(accessions):
    from Bio import Entrez

    Entrez.email = "j.camilo.ariasospina@gmail.com"

    filename = f"results/gene_sequences.fasta"

    with open(filename, "w") as outfile:

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
sequence_search(sequences_list)

# %% Check fasta file result
# in: gene_sequences.fasta (fasta files)
# out: visualization

records = list(SeqIO.parse("results/gene_sequences.fasta", "fasta"))
print(f"Number of sequences: {len(records)}\n")
for i, record in enumerate(records, start=1):
    print(f"{i:>4} | {record.id} | length = {len(record.seq)}")

# %% Encode sequences with One Hot Encoding
# in: gene_sequences (fasta file)
# out: ids_ohe (list), sequences_ohe (numpy array)

from one_hot_encoding import create_dataset_from_fasta

ids_ohe, sequences_ohe = create_dataset_from_fasta("results/gene_sequences_pathogenic.fasta")

# %% Create a dataset with labels and sequences
# in: data_filtered_pathogenic, data_filtered_benign (datasets)
# out: mutations_dataset (dataset)

nucleotide_numbers = {"A":0, "C":1, "G":2, "T":3, "U":4}
label_1 = [1] * (len(data_filtered_pathogenic))
label_2 = [0] * (len(data_filtered_benign))
label = label_1 + label_2
gene = []
location = []
mutation_from = []
mutation_to = []

data_filtered_merge = pd.concat([data_filtered_pathogenic, data_filtered_benign], ignore_index=True)

for a in data_filtered_merge["Name"]:
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

# %% induce mutations to sequences
# in: ids_ohe (list), sequences_ohe (numpy array), mutations_dataset (dataset)
# out: ids_complete (list), sequences_complete (numpy array), label (list)
label = [0] * len(ids_ohe)
ids_complete = ids_ohe
sequences_complete = sequences_ohe   

for a in range(len(mutations_dataset)):
    row = mutations_dataset.iloc[a]
    index = ids_ohe.index(row["gene"])
    seq = sequences_ohe[index]
    seq[(int(row["mutation_from"])), (int(row["location"])-1)] = False
    seq[(int(row["mutation_to"])), (int(row["location"])-1)] = True
    label.append(row["label"])
    ids_complete.append(row["gene"])
    sequences_complete.append(seq)

print(len(label), len(ids_complete), len(sequences_complete))
# %% Center sequences and filter size

# %% Run ML/DL

# %%


