#!/usr/bin/env Rscript
################# define functions 
get.args <- function () {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) == 5) {
    if (!file.exists(args[1])) {
      stop("input rank value not specified.n", call.=FALSE)
    }
  } else {
    stop("optimize_tangent.R has no input", call.=FALSE)
  }
  return(args)
}
load.packages <- function() {
	package.list = c(
		'ggplot2'
	)
	suppressMessages(suppressWarnings(sapply(package.list, library, character.only = TRUE, quietly = TRUE)))

	return()
}
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}
set.custom.ggplot.theme <- function () {
	custom.ggplot.theme = list(
		theme_bw() +
			theme(
				plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
				legend.position="right",
				legend.margin=margin(0,0,0,0),
				legend.text=element_text(size=8),
				legend.title=element_text(size=8),
				axis.title.y=element_text(size=8), 
				axis.title.x=element_text(size=8), 
				axis.text.y=element_text(lineheight=1.2, size=8, angle=0, hjust = 1), 
				axis.text.x=element_text(lineheight=1.2, size=8, angle=0, hjust = 1), 
				strip.text.x=element_text(size = 8), 
				strip.text.y=element_text(size = 8, angle=0, hjust = 0),
				strip.background=element_rect(color="white", fill="white", size=1.5, linetype="solid"),
				plot.title = element_text(size=12, hjust = 0.5)
			)
	)
	return(custom.ggplot.theme)
}
plot.elbow <- function (input.rank.df, custom.ggplot.theme, y_cutoff, plot.title, xaxis.title, yaxis.title, output.pdf.path) {
	
	elbow.plot <- ggplot(input.rank.df, aes(x=rank, y=score)) +
		geom_point(alpha = 0.25) +
		geom_line() +
		geom_hline(yintercept = y_cutoff) +
		annotate("text",x=10,y=y_cutoff,label=paste("cutoff=",y_cutoff), hjust= 0, vjust = -1) +
		ggtitle(plot.title) +
		xlab(xaxis.title) +
		ylab(yaxis.title) +
		custom.ggplot.theme
	ggsave(output.pdf.path, elbow.plot, width = 5, height = 5)
	
	return()
	
}
################# run
null <- load.packages()
args <- get.args()
#args <- c(
#	'/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/test_run/annotate_super_enhancer/demo/bed/tmp/distal_locus_count_rank.txt',
#	'/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/test_run/annotate_super_enhancer/demo/bed/tmp/distal_locus_count_rank.elbow.pdf', 
#	'hyperactive_distal_locus', 
#	'signal_count',
#	'distal_locus_rank'
#)
input.rank.path <- args[1]
output.pdf.path <- args[2]
plot.title <- args[3]
xaxis.title <- args[4]
yaxis.title <- args[5]
input.rank.df <- read.table(input.rank.path, header = FALSE)

colnames(input.rank.df) <- c('score')
input.rank.df$rank <- NA
input.rank.df$rank[order(input.rank.df$score)] <- 1:nrow(input.rank.df)
input.rank.vtr <- sort(input.rank.df$score)

#This is the slope of the line we want to slide. This is the diagonal.
slope <- (max(input.rank.vtr)-min(input.rank.vtr))/length(input.rank.vtr)

#Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
opt_res <- optimize(numPts_below_line,lower=1,upper=length(input.rank.vtr), myVector=input.rank.vtr,slope=slope)
x_point <- floor(opt_res$minimum)

#The y-value at this x point. This is our cutoff.
y_cutoff <- input.rank.vtr[x_point]

#plot elbow
custom.ggplot.theme <- set.custom.ggplot.theme()
null <- plot.elbow(input.rank.df, custom.ggplot.theme, y_cutoff, plot.title, xaxis.title, yaxis.title, output.pdf.path)
cat(paste(y_cutoff,"\n"))
