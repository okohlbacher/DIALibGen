#!/bin/bash
# dn.sh THREADS OUT -- DIA-NN 2.0 predicted library for the full proteome (owner's command line)
T=$1; O=$2; TMP=${O%.parquet}_tmp; mkdir -p $TMP
exec /ceph/ibmi/abi/oliver/opt/diann/diann-2.0/diann-linux --fasta /scratch/kohlbach/bid10-instrument/in/proteome.fasta --fasta-search --predictor --gen-spec-lib --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --min-pr-charge 2 --max-pr-charge 4 --min-pr-mz 350 --max-pr-mz 1200 --min-fr-mz 200 --max-fr-mz 1800 --unimod4 --met-excision --var-mods 0 --pg-level 1 --threads $T --temp $TMP --out-lib $O
