#!/usr/bin/env Rscript
################# set parameters 
max_train_sample_num <- 10000
num_digit <- 3
validation_fold <- 10
predictor.comb.operator <- "*"
################# define functions 
get.args <- function () {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) > 0) {
    if (length(args)==0) {
      stop("At least one argument must be supplied (input file).n", call.=FALSE)
    } else if (length(args)==1) {
      # default output dir
      args[2] = "./testing"
    }
  } else {
    args <- c(
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/1.0/demo/output/sc.solo/filter/demo/glm/tssCluster.glm.predictors.tsv',      
      '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/1.0/demo/output/sc.solo/filter/demo/glm/build_glm'      
    )
  }
  return(args)
}
load.packages <- function() {
  package.list = c(
    'ROCR','PRROC', 'scales', 'caret', 'e1071', 'ggplot2', 'reshape2'
  )
  suppressMessages(suppressWarnings(sapply(package.list, library, character.only = TRUE, quietly = TRUE)))
  return()
}
read.tssCluster.info <- function() {
  all.tssCluster.df <- read.table(file=input.file,  sep="\t", header=TRUE)
  all.tssCluster.df$ung_pct <- round(rescale(all.tssCluster.df$ung_pct, to=c(0,1)),num_digit)
  all.tssCluster.df$bkgd_rltv_expr <- round(rescale(all.tssCluster.df$bkgd_rltv_expr, to=c(0,1)), num_digit)
  all.tssCluster.df$summit_count <- round(rescale(log2(all.tssCluster.df$summit_count), to=c(0,1)),num_digit)
  all.tssCluster.df$flank_count <- round(rescale(log2(all.tssCluster.df$flank_count), to=c(0,1)),num_digit)
  all.tssCluster.df$cluster_count <- round(rescale(log2(all.tssCluster.df$cluster_count), to=c(0,1)),num_digit)
  all.tssCluster.df$size <- round(rescale(all.tssCluster.df$size, to=c(0,1)),num_digit)
  all.tssCluster.df$peakness <- round(rescale(all.tssCluster.df$peakness, to=c(0,1)),num_digit)
  all.tssCluster.df$stability <- round(rescale(all.tssCluster.df$stability, to=c(0,1)), num_digit)

  write.table(all.tssCluster.df,file=paste0(output.dir,"/scaled.predictor.values.tsv"), row.names = FALSE, quote = FALSE, sep="\t")
  
  return(all.tssCluster.df)  
}
run.ROC <- function(predict, binary, ROC.dir) {
  
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
  cat(ROCR_optcut(ROCR_fpr_tpr, ROCR_pred), file=paste0(ROC.dir,"/optimal_cutoff.txt"))
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
}
train.glm.model <- function() {
  
  predictor.comb.list <- list(
    "combined.predictors" = c("ung_pct", "bkgd_rltv_expr", "flank_count"),
    #00"combined.predictors" = c("ung_pct", "bkgd_rltv_expr", "summit_count", "flank_count", "cluster_count"),
    #02"combined.predictors" = c("ung_pct", "bkgd_rltv_expr", "cluster_count"),
    #01"combined.predictors" = c("ung_pct", "bkgd_rltv_expr"),
    "solo.stability" = c("stability"),
    "solo.size" = c("size"),
    "solo.peakness" = c("peakness"),
    "solo.unencoded_G" = c("ung_pct"),
    "solo.bkgd_score" = c("bkgd_rltv_expr"),
    "solo.summit_count" = c("summit_count"),
    "solo.flank_count" = c("flank_count"),
    "solo.cluster_count" = c("cluster_count")
  )
  
  golden_train.tssCluster.df <- all.tssCluster.df[all.tssCluster.df$train_binary > -1,]
  train_sample_num <- max_train_sample_num
  golden_sample_num <- nrow(golden_train.tssCluster.df)
  if (max_train_sample_num > golden_sample_num) {
    train_sample_num <- golden_sample_num 
  }
  train.data.df <- golden_train.tssCluster.df[sample(1:train_sample_num, ), ]
  train.data.df$train_binary <- factor(train.data.df$train_binary)
  train.data.df <- subset(train.data.df, select = c(train_binary, ung_pct, bkgd_rltv_expr, summit_count, flank_count, cluster_count, size, peakness, stability))
  
  glm.model.list <- list()
  for (predictor.comb in names(predictor.comb.list)) {
    training.formula <- as.formula(paste("train_binary~", paste(predictor.comb.list[[predictor.comb]], collapse=predictor.comb.operator)))
    set.seed(123) 
    #https://daviddalpiaz.github.io/r4sl/the-caret-package.html
    glm.model <- train(
      training.formula, 
      data = train.data.df, 
      trControl = trainControl(method = "cv", number = validation_fold),
      method = "glm",
      family = "binomial"
    )
    glm.model.list[[predictor.comb]] <- glm.model
  }
  
  return(glm.model.list)
}
evaluate.glm.model <- function() {
  
  plot.title.list <- list(
    "combined.predictors" = "Logistic model with combined predictors",
    "solo.stability" = "Logistic model with paraclu clustering stability",
    "solo.size" = "Logistic model with span of the cluster ",
    "solo.peakness" = "Logistic model with peakness",
    "solo.unencoded_G" = "Logistic model with % UMI with unencoded G",
    "solo.bkgd_score" = "Logistic model with background expression normalized score",
    "solo.flank_count" = "Logistic model with flanking UMI count",
    "solo.cluster_count" = "Logistic model with cluster UMI count",
    "solo.summit_count" = "Logistic model with summit UMI count"
  )
  
  predictors.list <- c(
    'prob' = 'Logistic model\nprobability',
    'ung_pct' = '% UMI with\nunencoded G',
    'bkgd_rltv_expr' = 'Background\nexpression\nnormalized score',
    'flank_count' = 'Flanking\nUMI count',
    'cluster_count' = 'Cluster\nUMI count',
    'summit_count' = 'Summit\nUMI count',
    'size' = 'Span',
    'peakness' = 'Peakness',
    'stability' = 'Stability'
  )
  
  for (predictor.comb in names(glm.model.list)) {
    predictor.comb.dir <- paste0(output.dir,"/",predictor.comb)
    dir.create(predictor.comb.dir,recursive = TRUE, showWarnings = FALSE)
    glm.model <- glm.model.list[[predictor.comb]]
    writeLines(capture.output(summary(glm.model)),con=paste0(predictor.comb.dir,"/glm.model.summary.txt"))
    writeLines(capture.output(glm.model$results),con=paste0(predictor.comb.dir,"/glm.model.result.txt"))
    data_test.df <- all.tssCluster.df[all.tssCluster.df$train_binary > -1,]
    data_test.df <- subset(data_test.df, select = c(train_binary, ung_pct, bkgd_rltv_expr, summit_count, flank_count, cluster_count, size, peakness, stability))
    predictions <- predict(glm.model, data_test.df, type = 'prob')
    data_test.df$prob <- predictions$"1"
    run.ROC(data_test.df$prob, data_test.df$train_binary, predictor.comb.dir)
    anova.res <- anova(glm.model$finalModel,glm.model.list$combined.predictors$finalModel,test="Chisq")
    writeLines(capture.output(anova.res),con=paste0(predictor.comb.dir,"/anova.vs.combined.predictors.txt"))
    plot.data.df <- melt(data_test.df, id.vars=c("train_binary"))
    plot.data.df$train_binary <- factor(plot.data.df$train_binary)
    plot.data.df$predictors <- factor(plot.data.df$variable)
    plot <- ggplot(plot.data.df, aes(x=value, color=train_binary, fill=train_binary, group=train_binary)) +
      #geom_density(alpha = 0.4)+
      geom_histogram(binwidth=0.005, alpha=0.5, size = 0.1, position="identity") +
      facet_grid(
        predictors~., 
        scales = 'free_y',
        labeller = labeller(predictors = predictors.list)
      ) +
      ggtitle(plot.title.list[[predictor.comb]]) +
      xlab('scaled predictor values') +
      ylab('% of tssClusters') +
      scale_fill_manual(name = "tssCluster testing set", values = c('0'='#e41a1c', '1'='#377eb8'), labels= c('0'='Low ATAC signal', '1'='High ATAC signal')) +
      scale_color_manual(name = "tssCluster testing set", values = c('0'='#e41a1c', '1'='#377eb8'), labels= c('0'='Low ATAC signal', '1'='High ATAC signal')) +
      theme_bw() +
      theme(
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="top",
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
    
    ggsave(paste0(predictor.comb.dir, "/predictor.distribution.pdf"), plot, width = 5, height = 8)
  }
  saveRDS(glm.model.list$combined.predictors, paste0(output.dir,"/combined.predictors.glm.model.RDS"))
}
################# run
load.packages()
args <- get.args()
input.file <- args[1]
output.dir <- args[2]
dir.create(output.dir,recursive = TRUE, showWarnings = FALSE)
all.tssCluster.df <- read.tssCluster.info()
glm.model.list <- train.glm.model()
evaluate.glm.model()
