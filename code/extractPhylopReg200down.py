import numpy as np
import pyBigWig
bw = pyBigWig.open("../data/PhyloP/hg38.phyloP100way.bw")
region=open("../data/PhyloP/PAS_200downpregions.bed","r")
fout=open("../data/PhyloP/PAS_phyloP_200downstream.txt","w")
for i, line in enumerate(region):
    if i >1:
        #print(line)
        cols = line.strip().split()
        vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
        valsNp=np.array(vals)
        meanVal=np.mean(valsNp)
        fout.write("%s\t%s\t%s\t%s\t%s\n"%(cols[0], cols[1], cols[2], cols[3], meanVal))
fout.close()
bw.close()
