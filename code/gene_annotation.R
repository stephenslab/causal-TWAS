# credit to https://gist.github.com/jwaageSnippets/5133941
# refGene format: http://genome.ucsc.edu/cgi-bin/hgTables
# Get the longest transcript start and end.

## not working for hg17, genomedb only supports hg18 and up, see https://bioconductor.org/packages/devel/bioc/manuals/GenomeInfoDb/man/GenomeInfoDb.pdf

# library("rtracklayer")
# session <- browserSession("UCSC")
# genome(session)<-"hg17"
# query <- ucscTableQuery(session, "refGene")
# tableName(query) <- "refGene"
# getTable(query) -> refseq

#-- download refGene from UCSC first and then read in
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two arguments must be supplied (refGene file, flanking size).n", call.=FALSE)
} else {
  refgfile <- args[1]
  fsize<- as.numeric(args[2]) # will add fsize upstream of tx start and minus fsize downstream of tx end
}
#--
refseq <- read.table(refgfile,header=F)
refseq <- refseq[,c(2,3,4,5,6,13)] #name,chr,strand,txStart,txEnd,name2
refseq$"width" <- abs(refseq$V5-refseq$V6)

refseq$V2 <- as.character(refseq$V2)
refseq$V13 <- as.character(refseq$V13)

split(refseq, refseq$V13) -> refseqSplit

getGene <- function(x)
{
  x[which.max(x$"width"),]
}
outdf <- do.call(rbind, lapply(refseqSplit, FUN=getGene))
outdf.NM <- outdf[grepl("NM",outdf$V2),] # refseq name start with NM means RNA, see here: https://en.wikipedia.org/wiki/RefSeq

outdf.NM$V5 <-  outdf.NM$V5 - fsize
outdf.NM$V6 <-  outdf.NM$V6 + fsize
outdf.NM$width <- NULL

write.table(outdf.NM , file=paste0("refGene_processed_", fsize,".txt"), row.names=F, col.names=F, sep="\t", quote = F)



