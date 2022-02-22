# Perform secondary structure alignment of MSTRG.5949.8 and AfuSINE2-1a
conda activate rna
RNAfold -in m5949_comparison.fa | RNAforester -l
