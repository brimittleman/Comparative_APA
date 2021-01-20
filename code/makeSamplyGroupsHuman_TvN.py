outfile=open("../Human/data/DiffIso_Human/sample_groups.txt", "w")
infile=open("../Human/data/CleanLiftedPeaks_FC/ALLPAS_postLift_LocParsed_Human_fixed.fc", "r")

for line, i in enumerate(infile):
    if line == 1:
        i_list=i.split()
        libraries=[]
        for sample in i_list[6:]:
            libraries.append(sample)
        for l in libraries:
            if l[-1] == "T":
                outfile.write("%s\tTotal\n"%(l))
            else:
                outfile.write("%s\tNuclear\n"%(l))
    else:
          next
