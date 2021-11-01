
#coverage.files<-list.files("~/coverage_plotting", full.names = TRUE, pattern = ".txt")
#coverage.names<-list.files("~/coverage_plotting", full.names = F, pattern=".txt")
args<- commandArgs(trailingOnly = TRUE)
coverage.file <-args[1]
#setwd("/Volumes/Georgia's Hard drive/temp_work")
pdf.file <- gsub("txt","pdf", coverage.file)
plot.colors <- c("red","blue","green","yellow","purple")
coverage <- read.delim(coverage.file)
pdf(pdf.file, width = 5, height= 4)
plot(-100,-100, xlim=c(0,250), ylim=c(1,1e6), xlab="Coverage", ylab="Number of basepairs", log="y")
colnames(coverage) <- c("contig", "position", "coverage")
contigs <- unique(coverage[,1])
for(j in 1:length(contigs)) {
    contig.cov <- subset(coverage,contig==contigs[j])
    cov.hist <- hist(contig.cov$coverage, plot=F, breaks=50)
    points(cov.hist$mids, cov.hist$counts, ty="o", col=plot.colors[j], pch=19, cex=0.5)
  }
  dev.off()






