"""Shared helpers for all pipeline stages."""

import logging
import os
import re
import sys
from pathlib import Path

import yaml

# Keywords that must match as whole tokens (short or ambiguous substrings)
STRICT_KEYWORDS = {"terl"}


def load_cfg(path=None):
    path = path or os.environ.get("PHAGE_CONFIG", "config.yaml")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["_config_path"] = str(path)
    return cfg


def setup_logging(stage, cfg):
    logs = Path(cfg["paths"]["logs_dir"])
    logs.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(logs / f"{stage}.log", mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )


def rdir(cfg, sub):
    """Return (and create) results_dir/<sub>."""
    p = Path(cfg["paths"]["results_dir"]) / sub
    p.mkdir(parents=True, exist_ok=True)
    return p


def contig_of(gene_id):
    """CTG_000001_12 -> CTG_000001"""
    return gene_id.rsplit("_", 1)[0]


def text(x):
    """pyhmmer/pyrodigal names may be bytes or str depending on version."""
    if x is None:
        return ""
    if isinstance(x, (bytes, bytearray)):
        return x.decode()
    return str(x)


def kw_match(keyword, haystack):
    k, h = keyword.lower(), haystack.lower()
    if k in STRICT_KEYWORDS:
        return re.search(r"(?<![a-z0-9])" + re.escape(k) + r"(?![a-z0-9])", h) is not None
    return k in h
