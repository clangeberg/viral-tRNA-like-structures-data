# Viral tRNA-like structure data

This repository contains the Stockholm alignments, calibrated Infernal covariance models, and derived data tables used to analyze viral tRNA-like structures (TLSs) and related 3′ RNA elements.

## Contents

- `alignments/`: eight curated Stockholm alignments.
- `covariance_models/`: the existing calibrated covariance models supplied with the study, preserved unchanged.
- `data/sequence_counts.csv`: record and alignment-column counts for every alignment.
- `data/entropy/`: Henikoff-weighted normalized Shannon entropy tables used for the TLSVal, TLSHis, and Bromovirus TLSTyr structure-mapped panels.
- `data/insertions/`: insertion-boundary frequencies, residue-level colors, and color-scale metadata used for the structure-mapped panels.
- `data/mappings/`: reference-sequence-to-PDB residue maps and the corresponding entropy values.

## Alignment counts

The current alignments contain 201 canonical TLS subtype records: 116 TLSVal, 33 TLSHis, 46 Bromoviridae TLSTyr, and 6 hordeivirus TLSTyr. The 52 TLSTyr subtype records reduce to 49 nonduplicate sequences in the shared TLSTyr core alignment because three subtype records are exact sequence duplicates. The related-element alignments contain 12 tobravirus 3′ pseudoknots, 50 type 1/2 TLEs, and 16 type 3 TLEs.

## Data definitions

Normalized Shannon entropy is calculated as `H/log(4)` from Henikoff-weighted, non-gap A/C/G/U frequencies. The scale therefore ranges from 0 for an invariant alignment position to 1 for a maximally diverse position. Insertion frequency is the fraction of evaluable alignment records containing at least one nucleotide within a reference-gap block bounded by two retained reference nucleotides. Frequencies are mapped to both experimentally resolved flanking residues.

The deposited models are the existing study models from `Working_SI/CM`; they were built and calibrated with Infernal 1.1.2 and were copied into this repository without recalculation. Updated covariance models are intended for submission to Rfam following bioRxiv deposition.
