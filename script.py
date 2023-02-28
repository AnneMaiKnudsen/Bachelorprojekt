import sys
import os
import gzip
import pandas as pd
from collections import defaultdict
from contextlib import redirect_stdout
from Bio import SeqIO
import re 
from ete3 import Tree

excluded_species = ["tupBel1", "mm10", "canFam3"]

def write_phylip(seqs, output_file):
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    
    with open(output_file, "w") as dest, redirect_stdout(dest):
        lengths=[len(x) for x in seqs.values()]
        assert(all(lengths[0] == l for l in lengths))
        seq_length = lengths[0]

        for name in seqs:
            seqs[name]=seqs[name].replace("-", "?")
        
        print(f"{len(seqs)} {seq_length}")
        for name in sorted(seqs.keys(), key=lambda x: x!="hg38"):
            print(f"{name:<10}{seqs[name]}")
            