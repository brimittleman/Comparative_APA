def main(inFile,letter, libraryName):
    Fin=open(inFile, "r")
    numNs=0
    numReads=0
    #print(letter)
    for i, ln in enumerate(Fin):
        if i >0:
            numReads +=1
            seq=ln.split()[15]
            seqU=seq.upper()
            #print(seqU)
            newN=seqU.count(letter)
            #print(newN)
            numNs += newN
            #print(numNs)
    print("%s\t%s\t%d\t%d\n"%(libraryName,letter, numNs, numReads))


if __name__ == "__main__":
  import sys
  inFile =sys.argv[1]
  Letter=sys.argv[2]
  libraryName=inFile.split("/")[-1].split("-")[0]
  main(inFile, Letter, libraryName)
