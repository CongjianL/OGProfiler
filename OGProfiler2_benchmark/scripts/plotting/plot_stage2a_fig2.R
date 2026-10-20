#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: plot_stage2a_fig2.R FIGURE_DATA_DIR OUTPUT_DIR")
data_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(patchwork)
library(ggrepel)

official <- read.delim(file.path(data_dir, "official_precision_recall.tsv"), check.names = FALSE)
refog <- read.delim(file.path(data_dir, "refog_F1_long.tsv"), check.names = FALSE)
split_contam <- read.delim(file.path(data_dir, "split_contamination.tsv"), check.names = FALSE)
fp <- read.delim(file.path(data_dir, "fp_concentration_long.tsv"), check.names = FALSE)
methods <- c("OGProfiler2", "OGProfiler2_default", "OGProfiler2_soft42", "OGProfiler1First", "OrthoFinder3", "FastOMA", "SonicParanoid2", "Proteinortho6")
palette <- c(OGProfiler2="#888888", OGProfiler2_default="#C44E52", OGProfiler2_soft42="#008C95", OGProfiler1First="#E39C37", OrthoFinder3="#4C78A8", FastOMA="#59A14F", SonicParanoid2="#B07AA1", Proteinortho6="#9C755F")
display <- c(OGProfiler2="V2 legacy", OGProfiler2_default="V2 default", OGProfiler2_soft42="V2 soft42", OGProfiler1First="V1 first", OrthoFinder3="OrthoFinder3", FastOMA="FastOMA", SonicParanoid2="SonicParanoid2", Proteinortho6="Proteinortho6")
shapes <- c(OGProfiler2=16, OGProfiler2_default=17, OGProfiler2_soft42=15, OGProfiler1First=18, OrthoFinder3=16, FastOMA=17, SonicParanoid2=15, Proteinortho6=18)
set.seed(42)
for (x in list(official, refog, split_contam, fp)) if (!all(unique(x$method) %in% methods)) stop("unknown method")
official$method <- factor(official$method, methods)
refog$method <- factor(refog$method, methods)
split_contam$method <- factor(split_contam$method, methods)
fp$method <- factor(fp$method, methods)

theme_set(theme_classic(base_size=8, base_family="sans") + theme(
  axis.line=element_line(linewidth=0.35), axis.ticks=element_line(linewidth=0.35),
  legend.position="none", plot.title=element_text(size=8, face="plain", hjust=0),
  axis.title=element_text(size=8), axis.text=element_text(size=6.5),
  plot.tag=element_text(size=10, face="bold"), plot.tag.position=c(0,1)
))

pA <- ggplot(official, aes(x=recall, y=precision, colour=method, shape=method)) +
  geom_point(size=2.2) + geom_text_repel(aes(label=display[as.character(method)]), size=2.1, seed=42, box.padding=0.25, point.padding=0.15, max.overlaps=Inf, show.legend=FALSE) +
  scale_colour_manual(values=palette) + scale_shape_manual(values=shapes) + coord_cartesian(xlim=c(0,100), ylim=c(0,103), clip="off") +
  labs(title="Official precision-recall separates method trade-offs", x="Recall (%)", y="Precision (%)", tag="A")

pB <- ggplot(refog, aes(x=method, y=F1, fill=method, colour=method)) +
  geom_violin(width=0.78, alpha=0.28, linewidth=0.35, trim=TRUE) +
  geom_jitter(width=0.13, height=0, size=0.65, alpha=0.55) +
  stat_summary(fun=median, geom="crossbar", width=0.5, linewidth=0.45, colour="black") +
  scale_fill_manual(values=palette) + scale_colour_manual(values=palette) +
  scale_x_discrete(labels=c(OGProfiler2="V2\nlegacy", OGProfiler2_default="V2\ndefault", OGProfiler2_soft42="V2\nsoft42", OGProfiler1First="V1\nfirst", OrthoFinder3="Ortho-\nFinder 3", FastOMA="Fast-\nOMA", SonicParanoid2="Sonic-\nParanoid 2", Proteinortho6="Protein-\northo 6")) +
  coord_cartesian(ylim=c(-0.03,1.03)) + theme(axis.text.x=element_text(size=5.8)) +
  labs(title="Paired RefOG performance reveals accuracy heterogeneity", x=NULL, y="Best-group F1 (70 RefOGs)", tag="B")

pC <- ggplot(split_contam, aes(x=median_split, y=median_contamination, colour=method, shape=method)) +
  geom_point(size=2.2) + geom_text_repel(aes(label=display[as.character(method)]), size=2.1, seed=42, direction="both", nudge_y=0.004, ylim=c(0.002,0.011), box.padding=0.3, point.padding=0.2, max.overlaps=Inf, show.legend=FALSE) +
  scale_colour_manual(values=palette) + scale_shape_manual(values=shapes) + coord_cartesian(ylim=c(0,0.0125), clip="off") + expand_limits(x=0) +
  labs(title="Split and contamination profiles", x="Median split count", y="Median contamination", tag="C")

pD <- ggplot(fp, aes(x=top_n, y=cumulative_FP_fraction, colour=method, group=method, shape=method)) +
  geom_line(linewidth=0.7) + geom_point(size=1.5) +
  scale_colour_manual(values=palette) + scale_shape_manual(values=shapes) + scale_x_continuous(breaks=c(1,5,10)) +
  scale_y_continuous(labels=scales::label_percent(accuracy=1), limits=c(0,1.03)) +
  geom_text_repel(data=subset(fp, top_n==10), aes(label=display[as.character(method)]), size=2.0, seed=42, direction="y", nudge_x=0.6, hjust=0, box.padding=0.2, point.padding=0.15, max.overlaps=Inf, show.legend=FALSE) +
  coord_cartesian(xlim=c(1,13), clip="off") +
  labs(title="False-positive concentration across families", x="Top-ranked predicted families", y="Cumulative pairwise FP burden", tag="D")

figure <- (pA | pC) / (pB | pD) + plot_layout(widths=c(1.08,1), heights=c(1,1.1))
width <- 183 / 25.4
height <- 160 / 25.4
ggsave(file.path(out_dir, "Fig2_candidate.pdf"), figure, width=width, height=height, units="in", device=grDevices::pdf)
ragg::agg_png(file.path(out_dir, "Fig2_candidate.png"), width=width, height=height, units="in", res=600, background="white")
print(figure)
dev.off()
