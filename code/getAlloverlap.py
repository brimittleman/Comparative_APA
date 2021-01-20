infile=open("../data/OverlappingPAS/FileMoreThan1Overlapping.txt", "r")
outfile=open("../data/OverlappingPAS/ListOverlappingPAS.txt","w")

for ln in infile:
    PASlist=ln.split()[3].split(",")
    for pas in PASlist:
        outfile.write("%s\n"%(pas))
outfile.close()
