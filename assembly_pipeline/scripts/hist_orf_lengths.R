data.file <- commandArgs(trailingOnly=T)[1]
pdf.file <- commandArgs(trailingOnly=T)[2]
pdf(file=pdf.file, height=3.5, width=4)

short.length <- 0.90
par(las=1)
nf <- layout(c(1,2), widths=c(1,1), heights=c(0.3,0.7) )
par(lwd=0.3)

d <- read.table(data.file)

## plot full hist
par(mar=c(0.05,4,1,1))
h <- hist(d[,1]/d[,2], ylim=c(100,5000), xlim=c(0,1.2), ylab="", xlab="", yaxt="n", main="", xaxt="n", col="light blue")
box(which="plot", lty=1)
axis(2,at=c(0,2000,4000), labels=c("0", "2K", "4K"), lwd=0.3)

## plot zoomed hist
par(mar=c(4,4,0.25,1))
h <- hist(d[,1]/d[,2], ylim=c(0,165), xlim=c(0,1.2), xlab="Ratio of query to hit length", main="", col="light blue", lwd=0.2)
box(which="plot", lty=1)

### write text of short reading frames
small <- length(which(d[,1]/d[,2]<short.length))
text(0.0, 80, labels=paste("short: ", small, "  (", round(small/dim(d)[1]*100,1), "% )", sep=""), cex=1., pos=4, offset=0.2)

dev.off()
