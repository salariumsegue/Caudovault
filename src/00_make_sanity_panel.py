#!/usr/bin/env python3
"""
src/00_make_sanity_panel.py
Builds a SYNTHETIC positive/negative control panel and matching control HMM profiles so the
whole pipeline can be verified before any real data is staged. Uses config.sanity.yaml.
Nothing here is biological data.
"""

import logging
import os
from pathlib import Path

from pyhmmer.easel import Alphabet, TextMSA, TextSequence
from pyhmmer.plan7 import Background, Builder

from common import load_cfg, setup_logging
from synth import encode_orf, random_protein, spacer

# profile_id -> (description in "category | annotation" form, seed)
FAMILIES = {
    "phrog_101": ("head and packaging | terminase large subunit", 101),
    "phrog_102": ("head and packaging | major capsid protein", 102),
    "phrog_103": ("head and packaging | portal protein", 103),
    "phrog_104": ("tail | tail tape measure protein", 104),
    "phrog_999": ("cellular_housekeeping | DNA-directed RNA polymerase subunit beta rpoB", 999),
}
MIN_PRIORITY_LEN = 10500


def gene(family_id):
    return encode_orf(random_protein(FAMILIES[family_id][1], 280), seed=FAMILIES[family_id][1])


def filler_gene(i):
    return encode_orf(random_protein(2000 + i, 260), seed=2000 + i)


def assemble(core_genes, seed_base):
    """Join genes with AT-rich spacers; pad with unannotated filler genes up to >=10.5 kb."""
    parts, i = [], 0
    for g in core_genes:
        parts += [spacer(seed_base + len(parts), 90), g]
    while sum(len(p) for p in parts) < MIN_PRIORITY_LEN:
        parts += [spacer(seed_base + 500 + i, 90), filler_gene(seed_base + i)]
        i += 1
    return "".join(parts) + spacer(seed_base + 999, 90)


def main():
    cfg = load_cfg(os.environ.get("PHAGE_CONFIG", "config.sanity.yaml"))
    setup_logging("00_make_sanity_panel", cfg)

    dtr = spacer(77, 150)  # 150 bp direct terminal repeat
    phage = dtr + assemble([gene("phrog_101"), gene("phrog_102"), gene("phrog_103"),
                            gene("phrog_104")], 100) + dtr
    cell = assemble([gene("phrog_999")] * 4, 300)
    short = spacer(5, 2500)

    ctrl = Path(cfg["paths"]["controls_fasta"])
    ctrl.parent.mkdir(parents=True, exist_ok=True)
    with open(ctrl, "w") as fh:
        fh.write(f">CTRL_POS_TAILED_PHAGE_1\n{phage}\n")
        fh.write(f">CTRL_NEG_CELLULAR_HOUSEKEEPING_1\n{cell}\n")
        fh.write(f">CTRL_SHORT_CONTIG_LT5KB\n{short}\n")
    logging.info(f"Wrote synthetic control panel -> {ctrl} "
                 f"(pos {len(phage)} bp, neg {len(cell)} bp, short {len(short)} bp)")

    hall = Path(cfg["paths"]["phrogs_hallmark_hmm"])
    allp = Path(cfg["paths"]["phrogs_all_hmm"])
    hall.parent.mkdir(parents=True, exist_ok=True)
    alpha = Alphabet.amino()
    builder, bg = Builder(alpha), Background(alpha)
    with open(hall, "wb") as fh_h, open(allp, "wb") as fh_a:
        for pid, (desc, seed) in FAMILIES.items():
            prot = random_protein(seed, 280)
            msa = TextMSA(name=pid, accession=pid, description=desc,
                          sequences=[TextSequence(name="s1", sequence=prot),
                                     TextSequence(name="s2", sequence=prot)]).digitize(alpha)
            hmm, _, _ = builder.build_msa(msa, bg)
            hmm.write(fh_h)
            hmm.write(fh_a)
    logging.info(f"Wrote synthetic control HMMs -> {hall}")


if __name__ == "__main__":
    main()
