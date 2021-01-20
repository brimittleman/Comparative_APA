def main(inFile, outFile, target):
    ifile=open(inFile, "r")
    ofile=open(outFile, "w")
    target=int(target)
    for num, ln in enumerate(ifile):
        if num == 0:
            ofile.write(ln)
        else:
            ID=ln.split()[0]
            chrom=ID.split(":")[0][3:]
            #print(chrom)
            if chrom == str(target):
                ofile.write(ln)

if __name__ == "__main__":
    import sys

    target = sys.argv[1]
    inFile = "../data/NuclearHvC_DF//ALLPAS_postLift_LocParsed_HvC_Nuclear_fixed4LC_Doublefilt.fc"
    outFile = "../data/DiffIso_Nuclear_DF/ALLPAS_postLift_LocParsed_HvC_Nuclear_fixed4LC_Doublefilt_chr%s.fc"%(target)
    main(inFile, outFile, target)
