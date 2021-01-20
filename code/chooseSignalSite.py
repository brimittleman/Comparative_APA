
def main(PASF, outF):
    signalDIC={'AATAAA':1, 'ATTAAA':2, 'AAAAAG':3, 'AAAAAA':4, 'TATAAA':5, 'AATATA':6, 'AGTAAA':7, 'AATACA':8, 'GATAAA':9, 'AATAGA':10, 'CATAAA':11, 'ACTAAA':12}
    PASdoc=open(PASF, "r")
    fout=open(outF, "w")
    PAS_dic={}
    for ln in PASdoc:
        PAS, spe, signal, count, ident = ln.split()
        if ident=="Y":
            if PAS not in PAS_dic.keys():
                 PAS_dic[PAS]=[signalDIC[signal]]
            else:
                 PAS_dic[PAS].append(signalDIC[signal])
    print(PAS_dic)
    for key, value in PAS_dic.items():
        choose=min(value)
        fout.write("%s\t%d\n"%(key, choose))
    fout.close()

if __name__ == "__main__":
  import numpy as np
  import sys
  PASfile =sys.argv[1]
  outFile= sys.argv[2]
  main(PASfile, outFile)
