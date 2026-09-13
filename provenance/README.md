# Provenance

The alignments in this repository are the curated September 3, 2026 Stockholm files. The update combined searches of `Viral_RefSeq_08262026.fasta` and a taxonomy-filtered GenBank viral nucleotide data set excluding SARS-CoV-2 sequences. The archived `cmsearch` reports identify Infernal 1.1.2; the inclusion E-value used for candidate collection was 0.01. Candidate regions were projected onto the original seed consensus coordinates and manually filtered for complete, terminal, class-compatible structures.

Consensus structures were evaluated with the online CaCoFold implementation of R-scape v2.6.16 using default parameters and an E-value threshold of 0.05. The covariance models deposited here are the existing calibrated study files from `Working_SI/CM`. Their embedded command records show that they were built with Infernal 1.1.2 using automatic `cmbuild` architecture selection and calibrated with `cmcalibrate --cpu 8`; the files were copied into this repository unchanged and were not recalculated for the deposit.

The three entropy and insertion analyses use the alignments corresponding to the experimental structures shown in the manuscript: TLSVal, TLSHis, and the Bromoviridae E/B3 TLSTyr subclass.
