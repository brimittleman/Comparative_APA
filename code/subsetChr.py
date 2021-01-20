def main(fin, fout):
  small=open(fout, "w")
  chrkeep=["chr1", "chr2", 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21','chr22']
  for ln in open(fin,"r"):
    chrom=ln.split()[0]
    if chrom in chrkeep:
        small.write(ln)
  small.close()


if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    Outfile =sys.argv[2]
    main(infile, Outfile)
