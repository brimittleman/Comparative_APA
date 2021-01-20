library(dplyr)
library(optparse)

option_list = list(
    make_option(c("-J", "--Junction"), action="store", default=NA, type='character', help="Original Junction"),
    make_option(c("-L", "--lifted"), action="store", default=NA, type='character', help="Reverse lifted"),
    make_option(c("-O", "--output"), action="store", default=NA, type='character', help="Output file")
)


opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


InJunc= read.table(opt$Junction, col.names=c("chr","start", "end", "JuncName", "score", "strand", "thickS", "thickE", "RGB", "count", "blockSize", "blockStart"),stringsAsFactors = F)

InLift= read.table(opt$lifted, col.names=c("chr","start", "end", "JuncName", "score", "strand", "thickS", "thickE", "RGB", "count", "blockSize", "blockStart"),stringsAsFactors = F)

Final=InLift %>% semi_join(InJunc, by=c("chr", "start","end","JuncName"))


write.table(Final, opt$output, quote=F, sep="\t", row.name=F, col.names=F)
