#SAF to 15 up bed
#usage python bed215upbed_gen.py  input outfile
def main(inputF, output):
    fin=open(inputF,"r")
    fout=open(output,"w")
    for ln in fin:
      chr, start, end, name, score, strand = ln.split()
      if strand=="+":
          start_new=int(start)-15
          end_new=int(start)
      else:
          start_new=int(end)
          end_new=int(end)+15
      score="."
      fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chr,start_new,end_new, name, score, strand))


if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    main(inFile, outFile)
