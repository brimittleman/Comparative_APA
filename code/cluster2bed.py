#chr2:219001:224864:clu_1_NA 0 8 4 12 24 4

#cluster 2 bed file


def main(clus, bed):
    fout=open(bed, "w")
    fin=open(clus, "r")
    for i, ln in enumerate(fin):
        if i > 0:
            cluster=ln.split()[0]
            #print(cluster)
            chr, start, end, site= cluster.split(":")
            score="."
            strand="."
            fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chr, start, end, site, score, strand))
    fout.close()

if __name__ == "__main__":
  from misc_helper import *
  import sys
  cluster =sys.argv[1]
  bed=sys.argv[2]
  main(cluster,bed)
