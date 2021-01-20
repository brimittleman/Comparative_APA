library(ggplot2)
library(dplyr)
library(optparse)

option_list = list(
    make_option(c("-g", "--gene"), action="store", default=NA, type='character', help="Gene Input"),
    make_option(c("-x", "--exta"), action="store", default="extra", type='character', help="extra")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

geneUse= opt$gene

usage=read.table("../data/files4viz_nuclear/NuclearPASUsage.txt", header=T) %>% filter(Gene==geneUse)

usage$start=as.integer(usage$start)

plot=ggplot(usage,aes(x=(reorder(PAS, start)), y=numUsage, by=species, fill=species))+ geom_boxplot() + scale_fill_brewer(palette = "Dark2")+theme(axis.text.x = element_text(angle = 90)) + labs(y="Usage", title=paste("Usage by species for ", geneUse))


outputfile=paste("../data/files4viz_nuclear/NuclearUsageBoxplot", geneUse, ".pdf",sep="_")


ggsave(plot, filename=outputfile, height=5, width=10)
