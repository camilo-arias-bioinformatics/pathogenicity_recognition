# %% import necessary packages
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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

# %% Store unique values for gene names
# in: data_filtered_pathogenic and data_filtered_benign
# out: unique_genes_pathogenic, unique_genes_benign, gene_chromosome (lists)
def store_gene_names(dataset):
    genes = []

    for a in dataset["Name"]:
        start = a.find("(")
        end = a.find(")")
        genes.append(a[start+1:end])

    genes = list(set(genes))
    return(genes)
unique_genes_pathogenic = store_gene_names(data_filtered_pathogenic)
unique_genes_benign = store_gene_names(data_filtered_benign)

print(f"len pathogenic: {len(unique_genes_pathogenic)} len benign: {len(unique_genes_benign)}")

# %% Download gene sequences to fasta files and store information from downloaded sequences
# in: unique_genes_pathogenic and unique_genes_benign (lists)
# out: gene_sequences.fasta (fasta file), gene_positions (dataset)

sequences_list = unique_genes_pathogenic + unique_genes_benign
sequences_list = list(set(sequences_list))

email = "j.camilo.ariasospina@gmail.com"

def sequence_search(gene_symbols, file_name, organism="Homo sapiens"):
    from Bio import Entrez
    import re
    import time

    Entrez.email = email

    filename = file_name
    failed = []
    records = []

    with open(filename, "w") as outfile:

        for gene in gene_symbols:

            try:
                search_handle = Entrez.esearch(
                    db="gene",
                    term=f"{gene}[sym] AND {organism}[orgn]"
                )
                search_result = Entrez.read(search_handle)
                search_handle.close()

                if not search_result["IdList"]:
                    failed.append(gene)
                    continue

                gene_id = search_result["IdList"][0]

                gene_handle = Entrez.efetch(
                    db="gene",
                    id=gene_id,
                    retmode="xml"
                )
                gene_record = Entrez.read(gene_handle)
                gene_handle.close()

                loc = gene_record[0]["Entrezgene_locus"][0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]

                chr_acc_gi = loc["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
                start = int(loc["Seq-interval_from"])
                stop = int(loc["Seq-interval_to"])
                strand = loc.get("Seq-interval_strand", {}).get("Na-strand", {}).attributes.get("value", "plus")

                chromosome = None
                biosource_subtypes = gene_record[0]["Entrezgene_source"]["BioSource"].get("BioSource_subtype", [])
                for subtype in biosource_subtypes:
                    if subtype.attributes.get("value") == "chromosome":
                        chromosome = str(subtype)
                        break

                summary_handle = Entrez.esummary(db="nucleotide", id=chr_acc_gi)
                summary_result = Entrez.read(summary_handle)
                summary_handle.close()

                accession = summary_result[0].get("AccessionVersion", "")
                title = summary_result[0].get("Title", "")
                assembly_match = re.search(r"(GRCh\d+(\.p\d+)?|GRCm\d+(\.p\d+)?)", title)
                assembly = assembly_match.group(0) if assembly_match else title

                fetch = Entrez.efetch(
                    db="nucleotide",
                    id=chr_acc_gi,
                    rettype="fasta",
                    retmode="text",
                    seq_start=start + 1,
                    seq_stop=stop + 1,
                    strand=2 if strand == "minus" else 1
                )

                outfile.write(fetch.read())
                outfile.write("\n")
                fetch.close()

                records.append({
                    "gene": gene,
                    "start": start + 1,
                    "end": stop + 1,
                    "chromosome": chromosome if chromosome else accession,
                    "assembly": assembly
                })

                print(f"Downloaded {gene}")
                time.sleep(0.4)

            except Exception as e:
                failed.append(gene)
                print(f"Failed {gene}: {e}")

    if failed:
        with open("results/failed_genes.txt", "w") as failfile:
            for gene in failed:
                failfile.write(gene + "\n")

    gene_positions = pd.DataFrame(records, columns=["gene", "start", "end", "chromosome", "assembly"])

    print(f"Failed genes: {failed}")
    return gene_positions, failed

gene_positions, failed = sequence_search(sequences_list, "results/gene_sequences.fasta")

# %% Double check failed genes and handle
# in: failed genes manually checked
# out: gene_sequences.fasta (fasta file)

sequences_manually_verified = ["ALPK3", "SOS1", "ACADVL"]

# All failed genes were already downloaded from other sequences, no need for further handling

# %% Analyze fasta file results
# in: gene_sequences.fasta (fasta file), gene_positions, data_filtered_merge (datasets)
# out: Analysis

for col in gene_positions.columns:
    print(f"\n{col}")
    print(gene_positions[col].unique())

# Elliminate sequences that come from 'Homo sapiens mitochondrion, complete genome'

compare_data = pd.merge(gene_positions, data_filtered_merge, on="Name", how="inner")

# %% ****Elliminate inadequate sequences
# in: 
# out:



# %% Check for outliers
# in: gene_sequences.fasta (fasta files)
# out: visualization

lengths = [len(record.seq) for record in SeqIO.parse("results/gene_sequences.fasta", "fasta")]

plt.figure(figsize=(6, 8))
sns.boxplot(y=lengths)
plt.ylabel("Sequence length (bp)")
plt.title("Distribution of gene sequence lengths")
plt.tight_layout()
plt.savefig("results/sequence_length_violin.png", dpi=300)
plt.show()

# %% Elliminate outliers
# in: gene_sequences.fasta (fasta file)
# out: gene_sequences_filtered.fasta (fasta file)

records = list(SeqIO.parse("results/gene_sequences.fasta", "fasta"))
lengths = np.array([len(r.seq) for r in records])

q1 = np.percentile(lengths, 25)
q3 = np.percentile(lengths, 75)

filtered = [r for r in records if q1 <= len(r.seq) <= q3]

SeqIO.write(filtered, "results/gene_sequences_filtered.fasta", "fasta")

print(f"Kept {len(filtered)} of {len(records)} sequences")

# %% Encode sequences with One Hot Encoding
# in: gene_sequences (fasta file)
# out: ids_ohe (list), sequences_ohe (numpy array)

from one_hot_encoding import create_dataset_from_fasta

ids_ohe, sequences_ohe = create_dataset_from_fasta("results/gene_sequences_filtered.fasta")

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

for a in data_filtered_merge["GRCh38Location"]
    start = a.index(":")
    end = a.index(a[-4])
    location.append(a[(start+1):(end)])

    mutation_from.append(a[-3])
    mutation_to.append(a[-1])

mutations_dataset = pd.DataFrame({"gene": gene, "location": location, "mutation_from": mutation_from, "mutation_to": mutation_to})

print(mutations_dataset.head())

# %% induce mutations to sequences and center
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

# %% ****Run ML/DL
# in: 
# out: 

# %%
# print(gene_positions.columns)
print(data_filtered_merge.head())

