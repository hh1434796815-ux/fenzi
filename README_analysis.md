# Xylocopa mitochondrial genome annotation and comparative analysis

This repository now contains a complete, reproducible Python 3 workflow for the five mounted mitochondrial FASTA assemblies (`F1`, `F2`, `F4`, `F5`, and `F6`). The workflow uses a curated insect/Xylocopa-like mitochondrial gene order as the closest-reference background, produces 37-gene annotation tables, nucleotide-composition/skew summaries, codon usage/RSCU summaries, phylogenetic output, and publication-style vector figures.

## Re-run command

```bash
python3 analyze_mitogenomes.py
```

No third-party Python packages are required for the production workflow. If users later add package installation or network downloads, the first shell line should be:

```bash
export https_proxy=http://127.0.0.1:7897 http_proxy=http://127.0.0.1:7897 all_proxy=socks5://127.0.0.1:7897
```

## Main outputs

- `analysis_results/tables/all_gene_annotation.csv`: combined 37-locus annotations for all five genomes.
- `analysis_results/tables/F*_gene_annotation.csv`: per-genome annotation tables.
- `analysis_results/tables/nucleotide_composition_skew.csv`: T/C/A/G, AT/GC, GC skew, and AT skew for 13 PCGs plus 2 rRNAs.
- `analysis_results/tables/codon_usage_RSCU.csv`: codon counts and RSCU values for 13 PCGs.
- `analysis_results/tables/rearrangement_scores_RS.csv`: rearrangement score table.
- `analysis_results/tables/rearrangement_frequency_RF.csv`: per-gene rearrangement-frequency table.
- `analysis_results/phylogeny_newick.nwk`: Newick tree from the concatenated 13 PCGs + rRNAs matrix.
- `analysis_results/figures/`: Figures 1–7 as SVG vector files and same-named `.pdf` companions.

## Notes for interpretation

The five input files contain complete mitochondrial-length sequences but no external feature table. Therefore, the script performs reference-guided structural annotation using a conserved Xylocopa/Hymenoptera mitochondrial template rather than remote NCBI annotation. The generated figures and CSVs are appropriate as a reproducible first-pass analysis and should be cross-validated with specialist tools such as MITOS/tRNAscan-SE/IQ-TREE before journal submission if external execution is available.
