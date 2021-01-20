def main(inputFile, fileID, outputF):
  #PYTHON 3

  dic_IND = {}
  dic_BAM = {}

  for ln in open(fileID, "r"):
      bam, IND = ln.split("\t")
      IND = IND.strip()
      dic_IND[bam] = IND
      if IND not in dic_BAM:
          dic_BAM[IND] = []
      dic_BAM[IND].append(bam)


  #now I have ind dic with keys as the bam and ind as the values
  #I also have a bam dic with ind as the keys and bam as the values

  inds=list(dic_BAM.keys()) #list of ind libraries


  count_file=open(inputFile, "r")
  genes=[]
  for line , i in enumerate(count_file):
      if line > 1:
          i_list=i.split()
          id=i_list[0]
          id_list=id.split(":")
          gene=id_list[6].split("_")[0]
          if gene not in genes:
              genes.append(gene)

  #make the ind and gene dic
  dic_dub={}
  for g in genes:
      dic_dub[g]={}
      for i in inds:
          dic_dub[g][i]=0


  #populate the dictionary
  count_file=open(inputFile, "r")
  for line, i in enumerate(count_file):
      if line > 1:
          i_list=i.split()
          id=i_list[0]
          id_list=id.split(":")
          g= id_list[6].split("_")[0]
          #print(g)
          values=list(i_list[6:])
          list_list=[]
          for ind,val in zip(inds, values):
              list_list.append([ind, val])
          for num, name in enumerate(list_list):
              dic_dub[g][list_list[num][0]] += int(list_list[num][1])


  #write the file by acessing the dictionary and putting values in the table ver the value in the dic


  fout=open(outputF,"w")
  peak=["chrom"]
  inds_noL=[]
  for each in inds:
      indsNA= "NA" + each
      inds_noL.append(indsNA)
  fout.write(" ".join(peak + inds_noL) + '\n' )


  count_file=open(inputFile, "r")
  for line , i in enumerate(count_file):
      if line > 1:
          i_list=i.split()
          id=i_list[0]
          id_list=id.split(":")
          genefull=id_list[6]
          gene=id_list[6].split("_")[0]
          strand=id_list[5]
          loc=id_list[6].split("_")[1]
          pasid=id_list[0]+ "-" + id_list[1]
          #start=dic_geneS[id_list[5]]
          start=int(id_list[3])
          #end=dic_geneE[id_list[5]]
          end=int(id_list[4])
          buff=[]
          buff.append("%s:%d:%d:%s_%s_%s-%s"%(id_list[2], start, end, gene, strand, loc, pasid))
          for x,y in zip(i_list[6:], inds):
              b=int(dic_dub[gene][y])
              t=int(x)
              buff.append("%d/%d"%(t,b))
          fout.write(" ".join(buff)+ '\n')

  fout.close()


if __name__ == "__main__":
  import numpy as np
  from misc_helper import *
  import sys
  inFile =sys.argv[1]
  fileID=sys.argv[2]
  outFile= sys.argv[3]
  main(inFile, fileID,outFile)
