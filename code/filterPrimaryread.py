def main(inFile, outFile):
    samfile = pysam.AlignmentFile(inFile, "rb")
    outfile = pysam.AlignmentFile(outFile, "wb", template=samfile)
    for read in samfile.fetch():
        if read.is_secondary:
            continue
        else:
            outfile.write(read)
    samfile.close()
    outfile.close()


if __name__ == "__main__":
    import sys
    import pysam
    inFile =sys.argv[1]
    outFile= sys.argv[2]
    main(inFile, outFile)
