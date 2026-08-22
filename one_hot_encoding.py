"""
fasta_to_onehot.py

Convierte secuencias de nucleótidos en formato FASTA a una representación
One-Hot Encoding (vectores de unos y ceros).

Basado en la lógica de codificación utilizada en Inpactor2.py
(funciones: fasta2one_hot, one_hot2fasta, create_dataset_master).

Alfabeto usado: A, C, G, T, N  (5 canales)
Cada secuencia se representa como una matriz de forma (5, total_win_len):
    - Fila 0 -> A
    - Fila 1 -> C
    - Fila 2 -> G
    - Fila 3 -> T
    - Fila 4 -> N
"""

import sys
import numpy as np
from Bio import SeqIO

# Alfabeto de nucleótidos (orden fijo, igual que en Inpactor2.py)
LANGU = ['A', 'C', 'G', 'T', 'N']


def fasta2one_hot(sequence: str, total_win_len: int) -> np.ndarray:
    """
    Convierte una secuencia de nucleótidos (str) en una matriz One-Hot.

    Parámetros:
        sequence      : secuencia de nucleótidos (A, C, G, T, N)
        total_win_len : longitud fija de la ventana (padding con ceros si la
                         secuencia es más corta)

    Retorna:
        np.ndarray booleano de forma (5, total_win_len)
    """
    rep2d = np.zeros((5, total_win_len), dtype=bool)

    for pos_nucl, nucl in enumerate(sequence):
        if pos_nucl >= total_win_len:
            break  # evita desbordamiento si la secuencia es más larga que la ventana
        pos_lang = LANGU.index(nucl.upper())
        rep2d[pos_lang][pos_nucl] = 1

    return rep2d


def one_hot2fasta(dataset: np.ndarray) -> str:
    """
    Convierte una matriz One-Hot (5, L) de vuelta a su secuencia de nucleótidos.
    """
    fasta_seq = ""
    for j in range(dataset.shape[1]):
        if dataset[:, j].sum() > 0:
            pos = np.argmax(dataset[:, j])
            fasta_seq += LANGU[pos]
    return fasta_seq


def create_dataset_from_fasta(fasta_file: str, total_win_len: int = None):
    """
    Lee un archivo FASTA y convierte todas sus secuencias en un dataset
    One-Hot de forma (n_secuencias, 5, total_win_len).

    Si total_win_len es None, se usa automáticamente la longitud de la
    secuencia más larga del archivo.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        print(f"ERROR: no se encontraron secuencias en {fasta_file}")
        sys.exit(1)

    seqs = [str(r.seq) for r in records]
    ids = [str(r.id) for r in records]

    if total_win_len is None:
        total_win_len = max(len(s) for s in seqs)

    n = len(seqs)
    dataset = np.zeros((n, 5, total_win_len), dtype=bool)

    for i, seq in enumerate(seqs):
        dataset[i, :, :] = fasta2one_hot(seq, total_win_len)

    return ids, dataset


# if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python fasta_to_onehot.py <archivo.fasta> [longitud_ventana]")
        sys.exit(1)

    fasta_path = sys.argv[1]
    win_len = int(sys.argv[2]) if len(sys.argv) > 2 else None

    ids, one_hot_dataset = create_dataset_from_fasta(fasta_path, win_len)

    print(f"Secuencias procesadas: {len(ids)}")
    print(f"Forma del dataset One-Hot: {one_hot_dataset.shape}")  # (n, 5, L)

    # Guardar el resultado en un archivo .npy para uso posterior
    out_path = fasta_path.rsplit(".", 1)[0] + "_onehot.npy"
    np.save(out_path, one_hot_dataset)
    print(f"Dataset guardado en: {out_path}")

    # Ejemplo: mostrar la codificación de la primera secuencia
    print("\nEjemplo (primera secuencia, primeros 10 nucleótidos):")
    print(f"ID: {ids[0]}")
    print(one_hot_dataset[0, :, :10].astype(int))
# %%
