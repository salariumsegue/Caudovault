#!/usr/bin/env python3
"""
src/00_preflight.py
Hardware audit, package lock, optional-binary check, input file audit, and a synthetic
end-to-end smoke test (FASTA -> Pyrodigal -> PyHMMER).
"""

import hashlib
import importlib
import importlib.metadata as md
import logging
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from common import load_cfg, setup_logging

PACKAGES = [
    ("biopython", "Bio"), ("pandas", "pandas"), ("numpy", "numpy"),
    ("pyrodigal", "pyrodigal"), ("pyhmmer", "pyhmmer"), ("sourmash", "sourmash"),
    ("matplotlib", "matplotlib"), ("seaborn", "seaborn"), ("networkx", "networkx"),
    ("pyyaml", "yaml"),
]

OPTIONAL_BINARIES = [
    ("megahit", ["--version"]), ("genomad", ["--version"]), ("checkv", ["-h"]),
    ("blastn", ["-version"]), ("makeblastdb", ["-version"]), ("diamond", ["version"]),
    ("mmseqs", ["version"]), ("hmmsearch", ["-h"]), ("minimap2", ["--version"]),
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_hardware():
    logging.info("=== 1. HARDWARE & OS ===")
    u = platform.uname()
    logging.info(f"OS / kernel   : {u.system} {u.release}")
    logging.info(f"CPU arch      : {u.machine}")
    logging.info(f"Python        : {sys.version.split()[0]} ({sys.executable})")
    logging.info(f"Logical cores : {os.cpu_count()}")
    total, _, free = shutil.disk_usage(".")
    logging.info(f"Disk (cwd)    : {free / 1024**3:.2f} GB free / {total / 1024**3:.2f} GB total")
    if Path("/proc/meminfo").exists():
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith(("MemTotal:", "MemAvailable:")):
                k, v = line.split(":")[:2]
                logging.info(f"{k:<14}: {int(v.split()[0]) / 1024**2:.2f} GB")
    elif u.system == "Darwin":
        try:
            mem = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True))
            logging.info(f"Total RAM     : {mem / 1024**3:.2f} GB")
        except Exception as e:  # noqa: BLE001
            logging.warning(f"Could not read RAM via sysctl: {e}")


def audit_packages():
    logging.info("=== 2. PYTHON PACKAGES ===")
    locked, ok = [], True
    for pip_name, mod_name in PACKAGES:
        try:
            importlib.import_module(mod_name)
            ver = md.version(pip_name)
            logging.info(f"  [OK]      {pip_name:<12} {ver}")
            locked.append(f"{pip_name}=={ver}")
        except Exception as e:  # noqa: BLE001
            logging.error(f"  [MISSING] {pip_name:<12} {e}")
            ok = False
    if ok:
        Path("requirements.lock.txt").write_text("\n".join(locked) + "\n")
        logging.info("  -> wrote requirements.lock.txt")
    return ok


def audit_binaries():
    logging.info("=== 3. OPTIONAL EXTERNAL BINARIES ===")
    for name, args in OPTIONAL_BINARIES:
        exe = shutil.which(name)
        if not exe:
            logging.info(f"  [ABSENT]  {name:<12} (Python fallback used where one exists)")
            continue
        try:
            res = subprocess.run([exe] + args, capture_output=True, text=True, timeout=10)
            out = (res.stdout or res.stderr).splitlines()
            logging.info(f"  [READY]   {name:<12} {(out[0] if out else 'ok')[:70]}")
        except Exception as e:  # noqa: BLE001
            logging.warning(f"  [WARN]    {name:<12} {e}")


def audit_files(cfg):
    logging.info("=== 4. INPUT & REFERENCE FILES ===")
    for key, rel in cfg["paths"].items():
        if key.endswith("_dir") and key in ("results_dir", "logs_dir", "report_dir"):
            continue
        p = Path(rel)
        if p.is_file():
            logging.info(f"  [PRESENT] {key:<22} {rel} ({p.stat().st_size / 1024**2:.2f} MB, sha256 {sha256_file(p)[:16]}...)")
        elif p.is_dir():
            logging.info(f"  [DIR]     {key:<22} {rel} ({len(list(p.iterdir()))} items)")
        else:
            logging.info(f"  [PENDING] {key:<22} {rel}")


def smoke_test():
    logging.info("=== 5. SMOKE TEST (synthetic FASTA -> Pyrodigal -> PyHMMER) ===")
    try:
        import pyhmmer
        import pyrodigal
        from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextMSA, TextSequence
        from pyhmmer.plan7 import Background, Builder

        from synth import encode_orf, random_protein, spacer

        prot = random_protein(seed=11, n_aa=300)
        contig = spacer(1, 400) + encode_orf(prot, seed=3) + spacer(2, 400)
        genes = pyrodigal.GeneFinder(meta=True).find_genes(contig.encode())
        assert len(genes) >= 1, "Pyrodigal found 0 genes on the synthetic ORF"

        alpha = Alphabet.amino()
        seqs = [TextSequence(name=f"g{i}", sequence=g.translate().rstrip("*")).digitize(alpha)
                for i, g in enumerate(genes, 1)]
        block = DigitalSequenceBlock(alpha, seqs)
        msa = TextMSA(name="smoke", description="smoke | terminase large subunit",
                      sequences=[TextSequence(name="a", sequence=prot),
                                 TextSequence(name="b", sequence=prot)]).digitize(alpha)
        hmm, _, _ = Builder(alpha).build_msa(msa, Background(alpha))
        hits = list(pyhmmer.hmmsearch([hmm], block, Z=100000))[0]
        sig = [h for h in hits if h.evalue <= 1e-5]
        assert sig, "PyHMMER did not recover the synthetic hit"
        logging.info(f"  [PASS] {len(contig)} bp -> {len(genes)} ORF(s) -> HMM hit E={sig[0].evalue:.2e}")
        return True
    except Exception as e:  # noqa: BLE001
        logging.error(f"  [FAIL] {type(e).__name__}: {e}")
        return False


if __name__ == "__main__":
    cfg = load_cfg()
    setup_logging("00_preflight", cfg)
    logging.info(f"Preflight {datetime.now(timezone.utc).isoformat()} | config={cfg['_config_path']} | mode={cfg['run']['mode']}")
    audit_hardware()
    pk = audit_packages()
    audit_binaries()
    audit_files(cfg)
    sm = smoke_test()
    sys.exit(0 if (pk and sm) else 1)
