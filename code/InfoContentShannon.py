
def main(species, outFile):
    inPAS=open("../data/PAS_doubleFilter/PAS_5perc_either_HumanCoord_BothUsage_meta_doubleFilter.txt", "r")
    fout=open(outFile, "w")
    fout.write("gene\tnumPAS\tbase2\tbasee\n")
    geneDic={}
    #PAS disc gene loc chr start end Chimp Human strandFix
    for i, ln in enumerate(inPAS):
        if i > 0:
           if species== "Human":
               usage=ln.split()[8]
           else:
               usage=ln.split()[7]
           gene = ln.split()[2]
           if gene not in geneDic.keys():
               geneDic[gene]=[float(usage)]
           else:
               geneDic[gene].append(float(usage))
    # loop through and write out
    #print(geneDic)
    for key, value in geneDic.items():
        nPAS=len(value)
        base2= entropy(value, base=2)
        basee= entropy(value, base=e)
        fout.write("%s\t%s\t%s\t%s\n"%(key, nPAS, base2, basee))

if __name__ == "__main__":
    import sys
    from math import log, e
    from scipy.stats import entropy
    import numpy as np
    np.seterr(divide='ignore', invalid='ignore')
    species =sys.argv[1]
    if species == "Human":
      outFile= "../data/InfoContent/Human_InfoContent.txt"
    else:
      outFile= "../data/InfoContent/Chimp_InfoContent.txt"
    main(species, outFile)
