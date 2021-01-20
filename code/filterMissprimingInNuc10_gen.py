def main(Fin, Fout):
  outBed=open(Fout, "w")
  inBed=open(Fin, "r")
  for ind, ln in enumerate(inBed):
    if ind >=1:
      chrom,start, end, name, score, strand, pctAT, pctGC, A, C, G, T, N, Other, Length, Sequence = ln.split()
      sequenceU=Sequence.upper()
      sequenceFixed=sequenceU.replace("N", "T")
      nT=sequenceFixed.count("T")
      Tperc= float(nT) / float(Length)
      if Tperc < .7:
          if "TTTTTT" not in sequenceU:
              start_new=int(start)
              end_new=int(end)
              outBed.write("%s\t%d\t%d\t%s\t%s\t%s\n"%(chrom, start_new, end_new , name, score, strand))
  outBed.close()

if __name__ == "__main__":
    import sys
    inFile = sys.argv[1]
    outFile=sys.argv[2]
    main(inFile, outFile)
