#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: plot_split_qc.R DIAGNOSTICS_TSV OUTPUT_PREFIX")
d <- read.delim(args[[1]], check.names=FALSE)
ordered <- d[order(d$F1, d$RefOG), ]
draw <- function() {
  old <- par(mfrow=c(1,2), mar=c(7,4,2.5,1), mgp=c(2.3,0.7,0), tcl=-0.25,
             cex.axis=0.75, cex.lab=0.9, cex.main=0.75, las=1)
  on.exit(par(old))
  x <- seq_len(nrow(ordered))
  plot(x, ordered$F1, type="o", pch=16, cex=0.45, lwd=1,
       col="#2563A6", ylim=c(0,1.03), xaxt="n", xlab="",
       ylab="Best-group F1", main="Best-group recovery varies across 70 RefOGs")
  ticks <- seq(1,nrow(ordered),by=7)
  axis(1, at=ticks, labels=ordered$RefOG[ticks], las=2, cex.axis=0.6)
  mtext("RefOG (ordered by F1)", side=1, line=5.7, cex=0.9)
  mtext("higher = better", side=3, line=-1.2, adj=0.02, col="#2563A6", cex=0.7)
  plot(d$true_size, d$split_count, pch=21, bg="#D97706", col="white", cex=0.8,
       xlab="True RefOG size", ylab="Split count",
       main="Fragmentation increases for some RefOGs")
}
png(paste0(args[[2]], ".png"), width=2160, height=900, res=300); draw(); dev.off()
pdf(paste0(args[[2]], ".pdf"), width=7.2, height=3); draw(); dev.off()
