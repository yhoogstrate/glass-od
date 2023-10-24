#!/usr/bin/env python


#' horrible to implement neatly w/ GRanges


from tqdm import tqdm
import pandas as pdd
import HTSeq
import statistics
import numpy as np
import glob

export = {}
repli_scores = sorted([_.replace("data/wgEncodeUwRepliSeq/","").replace(".bigWig.hg38.out.sorted.bgr","") for _ in glob.glob("data/wgEncodeUwRepliSeq/*.bgr")])

for repli_score in tqdm(repli_scores):

    # build index 
    repli = pdd.read_table("data/wgEncodeUwRepliSeq/"+repli_score+".bigWig.hg38.out.sorted.bgr", header=None, names = ["chr","start","end","repli"])

    idx = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for i in tqdm(range(len(repli['chr']))):
        _chr = repli['chr'][i]
        _pos_s = repli['start'][i]
        _pos_e = repli['end'][i]
        _repli_score = repli['repli'][i]
            
        idx[HTSeq.GenomicInterval(_chr, _pos_s, _pos_e)] += _repli_score


    # get all probes + chr & pos & query index
    dff1 = pdd.read_table("data/Improved DNA Methylation Array Probe Annotation/EPIC/EPIC.hg38.manifest.tsv")
    chrs = dff1['CpG_chrm'].fillna('')
    for i in tqdm(range(len(chrs))):
        _chr_p = chrs[i]
        
        if _chr_p != "":
            _probe_id = dff1['Probe_ID'][i]
            _pos_s_p = dff1['CpG_beg'][i]
            _pos_e_p = dff1['CpG_end'][i]
            
            try:
                out_l = list(idx[HTSeq.GenomicPosition(_chr_p, int(_pos_s_p))])
            except:
                print("probe_id")
                print(_probe_id)
                print("")
                
                print('chr_p')
                print(_chr_p)
                print("")
                
                print("pos_s_p")
                print(_pos_s_p)
                print("")
            
            if len(out_l) == 0:
                out = np.nan
            elif len(out_l) == 1:
                out = out_l[0]
            else:
                out = statistics.mean(out_l)
        else:
            out = np.nan
        
        if _probe_id not in export:
            export[_probe_id] = {}
        
        export[_probe_id][repli_score] = out


with open("output/tables/repli-seq_data.txt", "w") as fh:
    fh.write("probe_id\t"+("\t".join(repli_scores)) + "\n")
    
    for probe_id in tqdm(sorted(export.keys())):
        fh.write(probe_id)
        
        for repli_score in repli_scores:
            fh.write("\t" + str(export[probe_id][repli_score]))
        
        fh.write("\n")



