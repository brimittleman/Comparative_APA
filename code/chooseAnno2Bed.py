def main(inFile, outF):
    outFile=open(outF, "w")
    for ln in open(inFile, "r"):
        chrom, start, end, peak, cov, strand, anno = ln.split()
        if anno==".":
            continue
        anno_lst=anno.split(",")
        if len(anno_lst)==1:
            gene=anno_lst[0].split(":")[1]
            loc=anno_lst[0].split(":")[0]
            #print("1 gene")
            start_i=int(start)
            end_i=int(end)
            ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
            outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
        else:
            type_dic={}
            for each in anno_lst:
                type_dic[each.split(":")[0]]=each.split(":")[1]
            if "utr3" in type_dic.keys():
                gene=type_dic["utr3"]
                loc="utr3"
                #peak_i=int(peak)
                start_i=int(start)
                end_i=int(end)
                ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
                #continue
            elif "end" in type_dic.keys():
                gene=type_dic["end"]
                loc="end"
                #peak_i=int(peak)
                start_i=int(start)
                end_i=int(end)
                ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
                #continue
            elif "cds" in type_dic.keys():
                gene=type_dic["cds"]
                loc="cds"
                #peak_i=int(peak)
                start_i=int(start)
                end_i=int(end)
                ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
                #continue
            elif "utr5" in type_dic.keys():
                gene=type_dic["utr5"]
                loc="utr5"
                #peak_i=int(peak)
                start_i=int(start)
                end_i=int(end)
                ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
                #continue
            elif "intron" in type_dic.keys():
                gene=type_dic["intron"]
                loc="intron"
                #peak_i=int(peak)
                start_i=int(start)
                end_i=int(end)
                ID="%s:%s:%d:%d:%s:%s_%s"%(peak, chrom, start_i, end_i, strand, gene, loc)
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, start_i, end_i,ID, cov, strand))
                #continue
            #else:
                #continue
    outFile.close()


if __name__ == "__main__":
    import numpy as np
    from misc_helper import *
    import sys
    infile =sys.argv[1]
    outfile= sys.argv[2]
    main(infile, outfile)
