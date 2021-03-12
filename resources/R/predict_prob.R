#!/usr/bin/env Rscript
################# set parameters 
num_digit <- 3
lenient_specificity <- 0.5
################# define functions 
get.args <- function () {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) == 4) {
    if (!file.exists(args[1])) {
      stop("input file not specified.n", call.=FALSE)
    }
    if (!file.exists(args[2])) {
      stop("input model not specified.n", call.=FALSE)
    } 
  } else {
    args <- c(
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/1.0/demo/output/sc.solo/filter/demo/glm/tssCluster.glm.predictors.tsv',      
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/1.0/demo/output/sc.solo/filter/demo/glm/build_glm/combined.predictors.glm.model.RDS',
      '0.5',
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/1.0/demo/output/sc.solo/filter/demo/glm/predict_prob'      
    )
  }
  return(args)
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
load.packages <- function() {
  package.list = c(
    'ROCR','PRROC', 'caret', 'e1071', 'ggplot2', 'scales'
  )
  suppressMessages(suppressWarnings(sapply(package.list, library, character.only = TRUE, quietly = TRUE)))
  return()
}
read.tssCluster.info <- function() {
  all.tssCluster.df <- read.table(file=input.table,  sep="\t", header=TRUE)
  all.tssCluster.df$ung_pct <- round(rescale(all.tssCluster.df$ung_pct, to=c(0,1)),num_digit)
  all.tssCluster.df$bkgd_rltv_expr <- round(rescale(all.tssCluster.df$bkgd_rltv_expr, to=c(0,1)), num_digit)
  all.tssCluster.df$summit_count <- round(rescale(log2(all.tssCluster.df$summit_count), to=c(0,1)),num_digit)
  all.tssCluster.df$flank_count <- round(rescale(log2(all.tssCluster.df$flank_count), to=c(0,1)),num_digit)
  all.tssCluster.df$cluster_count <- round(rescale(log2(all.tssCluster.df$cluster_count), to=c(0,1)),num_digit)
  
  return(all.tssCluster.df)  
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
predict.prob <- function() {
  
  glm.model <- readRDS(input.model)
  writeLines(capture.output(summary(glm.model)),con=paste0(output.dir,"/glm.model.summary.txt"))
  writeLines(capture.output(glm.model$results),con=paste0(output.dir,"/glm.model.result.txt"))
  predictions <- predict(glm.model, all.tssCluster.df, type = 'prob')
  all.tssCluster.df$prob <- predictions$"1"
  
  data_test.df <- all.tssCluster.df[all.tssCluster.df$test_binary > -1,]
  data_test.df <- subset(data_test.df, select = c(test_binary, prob))
  prob_cutoff.list <- run.ROC(data_test.df$prob, data_test.df$test_binary, output.dir)
  #vline_robust_cutoff <- list(geom_vline(xintercept = prob_cutoff.list$robust, linetype="longdash", color = "#d95f02", size=0.5))
  #vline_lenient_cutoff <- list(geom_vline(xintercept = prob_cutoff.list$lenient, linetype="longdash", color = "#7570b3", size=0.5))
  #vline_default_cutoff <- list(geom_vline(xintercept = default.cutoff, linetype="longdash", color = "#1b9e77", size=0.5))
  vline_robust_cutoff <- list(geom_vline(aes(xintercept = prob_cutoff.list$robust, color = "robust"), linetype="longdash", size=0.5))
  vline_lenient_cutoff <- list(geom_vline(aes(xintercept = prob_cutoff.list$lenient, color = "lenient"), linetype="longdash", size=0.5))
  vline_default_cutoff <- list(geom_vline(aes(xintercept = default.cutoff, color = "default"), linetype="longdash", size=0.5))
  geom_text_default_cutoff <- list(annotate(geom = 'text', label = paste0("default=",round(default.cutoff,digits=2)), x = default.cutoff, y = Inf, angle=270, hjust="inward", vjust=1.5, color = "#1b9e77"))
  geom_text_robust_cutoff <- list(annotate(geom = 'text', label = paste0("robust=",round(prob_cutoff.list$robust,digits=2)), x = prob_cutoff.list$robust, y = Inf, angle=270, hjust="inward", vjust=1.5, color = "#d95f02"))
  geom_text_lenient_cutoff <- list(annotate(geom = 'text', label = paste0("lenient=",round(prob_cutoff.list$lenient,digits=2)), x = prob_cutoff.list$lenient, y = Inf, angle=270, hjust="inward", vjust=-1.5, color = "#7570b3"))
  
  all.tssCluster.df <- all.tssCluster.df[order(all.tssCluster.df$prob, decreasing = TRUE),]
  write.table(data.frame(tssClusterID=all.tssCluster.df$tssClusterID,prob=all.tssCluster.df$prob),file=paste0(output.dir,"/combined.predictors.glm.prob.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  
  all.cum.plot <- ggplot(all.tssCluster.df,aes(prob)) + 
    stat_ecdf(geom = "step") +
    ggtitle("Cumulative distribution of tssCluster\nagainst logistic model probablity") +
    vline_robust_cutoff + vline_lenient_cutoff + vline_default_cutoff +
    #geom_text_default_cutoff + geom_text_robust_cutoff + geom_text_lenient_cutoff +
    xlab('Logistic model probablity') +
    ylab('Fraction of tssClusters') +
    scale_color_manual(name = "cutoff", values = c('robust'='#d95f02', 'lenient'='#7570b3', 'default'='#1b9e77'), labels= c('robust'='robust', 'lenient'='lenient', 'default'='default')) +
    custom.ggplot.theme
  
  all.hist.plot <- ggplot(all.tssCluster.df, aes(x=prob)) +
  	 geom_histogram(binwidth=0.005, alpha=0.5, size = 0.1, position="identity") +
  	 ggtitle("Percentage of tssCluster with\nvarious logistic model probablity") +
    vline_robust_cutoff + vline_lenient_cutoff + vline_default_cutoff +
    #geom_text_default_cutoff + geom_text_robust_cutoff + geom_text_lenient_cutoff +
    xlab('Logistic model probablity') +
    ylab('Percentage of tssClusters') +
    scale_color_manual(name = "cutoff", values = c('robust'='#d95f02', 'lenient'='#7570b3', 'default'='#1b9e77'), labels= c('robust'='robust', 'lenient'='lenient', 'default'='default')) +
    custom.ggplot.theme
  
  plot.data.df <- subset(all.tssCluster.df, select = c(test_binary, prob, anno_region))
  plot.data.df$test_binary <- factor(plot.data.df$test_binary)
  plot.data.df$anno_region <- factor(plot.data.df$anno_region)
  benchmark.hist.plot <- ggplot(plot.data.df, aes(x=prob, fill=test_binary, group=test_binary)) +
    geom_histogram(aes(y = stat(width*density)), bins=100, alpha=0.3, size = 0.1, position="identity") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    vline_robust_cutoff + vline_lenient_cutoff + vline_default_cutoff +
    #geom_text_default_cutoff + geom_text_robust_cutoff + geom_text_lenient_cutoff +
    ggtitle("Percentage of tssCluster with\nvarious logistic model probablity") +
    xlab('Logisitic model probability') +
    ylab('% of tssClusters') +
    scale_fill_manual(name = "tssCluster testing set", values = c('-1'='#636363', '0'='#e41a1c', '1'='#377eb8'), labels= c('-1'='Others', '0'='Negative\ntesting region', '1'='Positive\ntesting region')) +
    scale_color_manual(name = "cutoff", values = c('robust'='#d95f02', 'lenient'='#7570b3', 'default'='#1b9e77'), labels= c('robust'='robust', 'lenient'='lenient', 'default'='default')) +
    custom.ggplot.theme
  
  anno_reg.pct.hist.plot <- ggplot(plot.data.df, aes(x=prob)) +
    geom_histogram(aes(y = stat(width*density)), bins=100, alpha=0.8, size = 0.1, position="identity") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
    vline_robust_cutoff + vline_lenient_cutoff + vline_default_cutoff +
    #geom_text_default_cutoff + geom_text_robust_cutoff + geom_text_lenient_cutoff +
    ggtitle("Percentage of tssCluster with\nvarious logistic model probablity") +
    xlab('Logisitic model probability') +
    ylab('% of tssClusters') +
    facet_grid(
      anno_region~., 
      scales = 'free_y'
    ) +
    scale_color_manual(name = "cutoff", values = c('robust'='#d95f02', 'lenient'='#7570b3', 'default'='#1b9e77'), labels= c('robust'='robust', 'lenient'='lenient', 'default'='default')) +
    custom.ggplot.theme
  
  anno_reg.count.hist.plot <- ggplot(plot.data.df, aes(x=prob)) +
    geom_histogram(bins=100, alpha=0.8, size = 0.1, position="identity") +
    vline_robust_cutoff + vline_lenient_cutoff + vline_default_cutoff +
    #geom_text_default_cutoff + geom_text_robust_cutoff + geom_text_lenient_cutoff +
    ggtitle("Percentage of tssCluster with\nvarious logistic model probablity") +
    xlab('Logisitic model probability') +
    ylab('Number of tssClusters') +
    facet_grid(
      anno_region~.
    ) +
    scale_color_manual(name = "cutoff", values = c('robust'='#d95f02', 'lenient'='#7570b3', 'default'='#1b9e77'), labels= c('robust'='robust', 'lenient'='lenient', 'default'='default')) +
    custom.ggplot.theme
  
  ggsave(paste0(output.dir, "/prob.cum_distr.all.tssCluster.pdf"), all.cum.plot, width = 5, height = 5)
  ggsave(paste0(output.dir, "/prob.hist.all.tssCluster.pdf"), all.hist.plot, width = 5, height = 5)
  ggsave(paste0(output.dir, "/prob.hist.benchmark.tssCluster.pdf"), benchmark.hist.plot, width = 8, height = 4)
  ggsave(paste0(output.dir, "/prob.hist.anno_reg.pct.tssCluster.pdf"), anno_reg.pct.hist.plot, width = 5, height = 8)
  ggsave(paste0(output.dir, "/prob.hist.anno_reg.count.tssCluster.pdf"), anno_reg.count.hist.plot, width = 5, height = 8)
  
  return()
}
################# run
load.packages()
args <- get.args()
input.table <- args[1]
input.model <- args[2]
default.cutoff <- as.numeric(args[3])
output.dir <- args[4]
dir.create(output.dir,recursive = TRUE, showWarnings = FALSE)
custom.ggplot.theme <- set.custom.ggplot.theme()
all.tssCluster.df <- read.tssCluster.info()
predict.prob()
