#fixClusterformat
#change chr21:9127547:9129560:clu_1:? to chr2:3548720:3549058:clu_5_NA


def main(infile, outfile):
    fin=open(infile, "r")
    fout=open(outfile,"w")
    for i, ln in enumerate(fin):
        if i == 0:
            fout.write(ln)
        if i > 0:
            lineList=ln.split()
            cluster=lineList[0]
            newCluster=rreplace(cluster, ":", "_", 3)
            print(newCluster)
            lineList[0]=newCluster
            newLine=" ".join(lineList)
            fout.write("%s\n"%(newLine))
    fout.close()




def rreplace(s, old, new, occurrence):
    li = s.rsplit(old)
    li_new=old.join(li[:occurrence]), new.join(li[occurrence:])
    print(old.join(li_new))
    return old.join(li_new)


if __name__ == "__main__":
  import sys
  inFile =sys.argv[1]
  outFile= sys.argv[2]
  main(inFile, outFile)
