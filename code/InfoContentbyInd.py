def main(inFile, outFile):
    #ind dic
    fout=open(outFile,"w")
    indDic={}
    fin=open(inFile,"r")
    for i, ln in enumerate(fin):
        if i==0:
            header= ln.split()
            numCol=len(header)
            #print(numCol)
            for num, name in enumerate(header):
                if num > 1:
                    indDic[num]=name
    fin.close()
    #print(indDic)
    for ind in list(range(2, numCol)):
        fin=open(inFile,"r")
        geneDic={}
        for i, ln in enumerate(fin):
            if i > 0:
                gene=ln.split()[0]
                #print(ind)
                usage=ln.split()[ind]
                if gene not in geneDic.keys():
                    geneDic[gene]=[float(usage)]
                else:
                    geneDic[gene].append(float(usage))
        fin.close()
        #print(geneDic)
        for key, value in geneDic.items():
            #print(key)
            #print(value)
            if len(value) > 1:
                sh=shannon(value)
                eq=equit(value)
                simp=simpson(value)
                #print(ind)
                indiv=indDic[ind]
                fout.write("%s\t%s\t%s\t%s\t%s\n"%(key, indiv, sh, eq, simp))
    fout.close()

def shannon(usage):
    base2= entropy(usage, base=2)
    return(base2)

def equit(usage):
    base2= entropy(usage, base=2)
    nPAS=len(usage)
    eq=base2/log2(nPAS)
    return(eq)

def simpson(usage):
    nPAS=len(usage)
    simp=0
    for i in usage:
      simp += i**2
      sipOpp=1-simp
    return(sipOpp)

if __name__ == "__main__":
    import sys
    from math import log, e, log2
    from scipy.stats import entropy
    import numpy as np
    np.seterr(divide='ignore', invalid='ignore')
    infile =sys.argv[1]
    outFile =sys.argv[2]
    main(infile, outFile)
