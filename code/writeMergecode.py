files=open("../data/HC_filenames.txt", "r")
fout=open("comands2Mege.sh","w")

fout.write("#!/bin/bash\n\n\n#SBATCH --job-name=mergethreeprime\n#SBATCH --mail-type=END\n#SBATCH --error=merge.err\n#SBATCH --partition=broadwl\n#SBATCH --mem=36G\n\n\n")
for i, ln in enumerate(files):
    if i>0:
        id ,lane1, lane2, lane3, final = ln.split()
        fout.write("cat %s %s %s > %s\n"%(lane1,lane2,lane3, final))

fout.close()
