#python prepareCleanLiftedFC_5perc4LC.py fc bed out


def main(fc,anno, output):
    pas=open(anno, "r")
    pas_dic={}
    for i, ln in enumerate(pas):
        if i>0:
            okPAS=ln.split()[3]
            pas_dic[okPAS]=""
    fout=open(output, "w")
    inFC=open(fc, "r")
    for i, ln in enumerate(inFC):
        if i == 1:
            print(ln)
            lines=ln.split()[6:]
            fout.write(" ".join(lines)+'\n')
        #Human:human31:chr1:778450:778650:+:LOC100288069_utr3
        if i >1:
            fcPAS=ln.split()[0].split(":")[1]
            if fcPAS in pas_dic.keys():
                chrom=ln.split()[0].split(":")[2]
                start=ln.split()[0].split(":")[3]
                end=ln.split()[0].split(":")[4]
                gene=ln.split()[0].split(":")[6].split("_")[0]
                new_ID="%s:%s:%s:%s"%(chrom, start, end, gene)
                pheno=ln.split()[6:]
                pheno.insert(0, new_ID)
                fout.write(" ".join(pheno)+'\n')
    fout.close()



if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    anno=sys.argv[2]
    outFile = sys.argv[3]
    main(inFile,anno, outFile)
