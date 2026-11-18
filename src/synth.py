"""Deterministic synthetic sequence helpers (used ONLY for the sanity panel and smoke test).

Proteins are random (seeded); codons are chosen from synonymous codons with a mild
GC-biased usage so Prodigal's metagenomic models can call the ORFs.
These are NOT biological sequences.
"""

import random

from Bio.Data import CodonTable

_TABLE = CodonTable.unambiguous_dna_by_id[11]
_AA_TO_CODONS = {}
for _codon, _aa in _TABLE.forward_table.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)

_AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_W = [8.3, 1.4, 5.4, 6.7, 3.9, 7.1, 2.3, 5.9, 5.8, 9.7, 2.4, 4.0, 4.7, 3.9, 5.5, 6.6, 5.4, 6.9, 1.1, 2.9]


def random_protein(seed, n_aa=280):
    rng = random.Random(seed)
    return "M" + "".join(rng.choices(_AA, weights=_AA_W, k=n_aa - 1))


def encode_orf(protein, seed=0):
    rng = random.Random(seed)
    parts = []
    for aa in protein:
        codons = _AA_TO_CODONS[aa]
        weights = [1.0 + 0.6 * sum(c in "GC" for c in cod) for cod in codons]
        parts.append(rng.choices(codons, weights=weights, k=1)[0])
    parts[0] = "ATG"
    return "".join(parts) + "TAA"


def spacer(seed, n):
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", weights=[0.3, 0.2, 0.2, 0.3], k=n))
