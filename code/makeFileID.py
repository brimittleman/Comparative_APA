def main(infile,outfile):
  fout=open(outfile, "w")
  inf=open(infile,"r")
  for line, i in enumerate(inf):
      if line == 1:
          print(line)
          i_list=i.split()
          files= i_list[6:]
          for each in files:
              print(each)
              full = each.split("/")[4]
              samp= full.split("-")[0].split("_")[2:]
              lim="_"
              samp_st=lim.join(samp)
              outLine= full + "\t" + samp_st
              fout.write(outLine + "\n")
  fout.close()


if __name__ == "__main__":
    import numpy as np
    from misc_helper import *
    import sys
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
