#!/usr/bin/env python

import click
from pyfaidx import Fasta
from tqdm import tqdm



#@click.command()
#@click.option('--EPIC-manifest','-e', help='EPIC manifest file.')
#@click.option('--fasta','-f', help='Reference FASTA file.')
#@click.option('--size','-s', default=8, help='size of motif [on both sizes].')
#def main(EPIC_manifest, fasta, size):
EPIC_manifest = "/mnt/neuro-genomic-1-ro/catnon/Methylation - EPIC arrays/EPIC.hg38.manifest.tsv.gz"
fasta = "/data/bio/hg38/fasta/hg38.fa"

size=64

fh_fasta = Fasta(fasta)

with open(EPIC_manifest, "r") as fh_manifest:
    
    for line in tqdm(fh_manifest):
        line = line.split("\t")
        if line[0] != "CpG_chrm":
            #print(line)
            #chr1	10524	10526	-	cg14817997	
            
            if line[0] in fh_fasta:
                pre = fh_fasta[line[0]][int(line[1]) - size:int(line[1])]
                subseq = fh_fasta[line[0]][int(line[1]):int(line[2])]
                post = fh_fasta[line[0]][int(line[2]): int(line[2]) + size]
                
                if line[4][0:2] == "cg": # there are odd non-cg probes
                    print(line[4] + "\t" + line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t" + str(pre) + "|" + str(subseq) + "|" + str(post))


#main()

