#!/usr/bin/env python3
"""Write SHA-256 checksums for the deposited alignment, model, and data files."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
OUTPUT = REPOSITORY / "provenance" / "checksums_sha256.csv"
DEPOSIT_DIRECTORIES = ("alignments", "covariance_models", "data")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


files = sorted(
    path
    for directory in DEPOSIT_DIRECTORIES
    for path in (REPOSITORY / directory).rglob("*")
    if path.is_file()
)

with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=("path", "sha256", "bytes"))
    writer.writeheader()
    for path in files:
        writer.writerow(
            {
                "path": path.relative_to(REPOSITORY).as_posix(),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
        )

print(f"Wrote {len(files)} checksums to {OUTPUT.relative_to(REPOSITORY)}")
