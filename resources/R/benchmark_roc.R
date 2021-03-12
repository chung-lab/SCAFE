#!/usr/bin/env Rscript
################# set parameters 
lenient_specificity <- 0.5
################# define functions 
get.args <- function () {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) > 0) {
    if (!file.exists(args[1])) {
      stop("input file not specified.n", call.=FALSE)
    } else if (length(args)==1) {
      # default output dir
      args[2] = "./output/"
    }
  } else {
	stop("benchmark_roc.R has no input", call.=FALSE)
	args <- c(
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/v1.0/resources/R/00_test/00_reg_roc/roc.input.tsv',      
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/v1.0/resources/R/00_test/00_reg_roc'      
    )
  }
  return(args)
}
load.packages <- function() {
  package.list = c(
    'ROCR','PRROC', 'ggplot2'
  )
  suppressMessages(suppressWarnings(sapply(package.list, library, character.only = TRUE, quietly = TRUE)))
  return()
}
run.benchmark <- function() {
  data_test.df <- read.table(file=input.table,  sep="\t", header=TRUE)
  data_test.df <- data_test.df[data_test.df$binary > -1,]
  prob_cutoff.list <- run.ROC(data_test.df$prob, data_test.df$binary, output.dir)
}
run.ROC <- function(predict, binary, ROC.dir) {
  #predict <- data_test.df$prob
  #binary <- data_test.df$test_binary
  
  dir.create(ROC.dir,recursive = TRUE, showWarnings = FALSE)
  ROCR_pred <- prediction(predict, binary)
  ROCR_f_score	<- performance(ROCR_pred, "f")
  ROCR_fpr_tpr	<- performance(ROCR_pred, measure = "tpr", x.measure = "fpr")
  ROCR_prec_acc	<- performance(ROCR_pred, measure = "acc", x.measure = "prec")
  ROCR_rec_prec	<- performance(ROCR_pred, measure = "prec", x.measure = "rec")
  ROCR_optcut 	<- function(ROCR_fpr_tpr, ROCR_pred){cut.ind = mapply(FUN=function(x, y, p){d = (x - 0)^2 + (y-1)^2; ind = which(d == min(d)); c(opt_specificity = 1-x[[ind]], opt_sensitivity = y[[ind]], opt_cutoff = p[[ind]])}, ROCR_fpr_tpr@x.values, ROCR_fpr_tpr@y.values, ROCR_pred@cutoffs)}
  PRROC_pr		<- pr.curve(scores.class0 = predict,  weights.class0 = binary)
  PRROC_roc		<- roc.curve(scores.class0 = predict,  weights.class0 = binary)
  f_score_data_table <- cbind('f_score'=ROCR_f_score@y.values[[1]],'cutoff'=ROCR_f_score@x.values[[1]])
  ROC_data_table 	<- cbind('specificity'=1-ROCR_fpr_tpr@x.values[[1]],'sensitivity'=ROCR_fpr_tpr@y.values[[1]],'cutoff'=ROCR_fpr_tpr@alpha.values[[1]])
  FDR_data_table 	<- cbind('false_discovery_rate'=1-ROCR_prec_acc@x.values[[1]],'accuracy'=ROCR_prec_acc@y.values[[1]],'cutoff'=ROCR_prec_acc@alpha.values[[1]])
  PR_data_table		<- cbind('recall'=ROCR_rec_prec@x.values[[1]],'precision'=ROCR_rec_prec@y.values[[1]],'cutoff'=ROCR_rec_prec@alpha.values[[1]])
  
  trim.ROC_data_table.df <- as.data.frame(ROC_data_table)
  trim.ROC_data_table.df <- trim.ROC_data_table.df[trim.ROC_data_table.df$specificity <= lenient_specificity,]
  trim.ROC_data_table.df <- trim.ROC_data_table.df[order(trim.ROC_data_table.df$cutoff, decreasing = TRUE),]
  robust_prob_cutoff.list <- ROCR_optcut(ROCR_fpr_tpr, ROCR_pred)
  prob_cutoff.list <- list(
    'lenient' = trim.ROC_data_table.df$cutoff[1],
    'robust' = robust_prob_cutoff.list[3]
  )
  cat(c(prob_cutoff.list$lenient, prob_cutoff.list$robust), file=paste0(ROC.dir,"/lenient.robust.cutoff.txt"))
  cat(PRROC_pr$auc.integral, file=paste0(ROC.dir,"/PR_AUC.txt"))
  cat(PRROC_roc$auc, file=paste0(ROC.dir,"/ROC_AUC.txt"))
  write.table(f_score_data_table,file=paste0(ROC.dir,"/f_score.out_data.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  write.table(ROC_data_table,file=paste0(ROC.dir,"/ROC_curve.out_data.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  write.table(FDR_data_table,file=paste0(ROC.dir,"/FDR_accuracy.out_data.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  write.table(PR_data_table,file=paste0(ROC.dir,"/PR_curve.out_data.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  pdf(paste0(ROC.dir,"/ROC_curve.pdf"))
  plot(ROCR_fpr_tpr,colorize=TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7))
  dev.off()
  pdf(paste0(ROC.dir,"/f_score.pdf"))
  plot(ROCR_f_score)
  dev.off()
  pdf(paste0(ROC.dir,"/PR_curve.pdf"))
  plot(ROCR_rec_prec)
  dev.off()
  
  return(prob_cutoff.list)
}
################# run
load.packages()
args <- get.args()
input.table <- args[1]
output.dir <- args[2]
dir.create(output.dir,recursive = TRUE, showWarnings = FALSE)
run.benchmark()
