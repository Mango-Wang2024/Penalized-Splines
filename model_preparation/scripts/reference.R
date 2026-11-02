#!/usr/bin/env Rscript
run <- function() {
  suppressPackageStartupMessages(library(mgcv))
  script <- gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)),fixed=TRUE)
  project <- dirname(dirname(dirname(normalizePath(script))))
  source(file.path(project,"data","shared.R"))
  args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
  reference <- args[1]
  cases <- unique(read.csv(file.path(reference,"source_checksums.csv"))$case)
  for(name in cases) {
    path <- file.path(reference,name);terms <- read.delim(file.path(path,"terms.tsv"))
    X0 <- list()
    for(key in terms$predictor) {
      dst <- file.path(path,key);sm <- readRDS(file.path(dst,"objects.rds"))$prepared
      re <- mgcv::smooth2random(sm,c("y",terms$predictor),type=2)
      stopifnot(length(re$rand)==1L)
      U <- if(is.null(re$trans.U))diag(ncol(sm$X)) else re$trans.U
      D <- if(is.null(re$trans.D))rep(1,ncol(sm$X)) else re$trans.D
      C <- sweep(U,2,D,"*")
      stopifnot(qr(C)$rank==ncol(sm$X))
      for(nm in c("X0","Z","C"))save_matrix(switch(nm,X0=re$Xf,Z=as.matrix(re$rand[[1]]),C=C),nm,dst)
      X0[[key]] <- re$Xf
    }
    save_matrix(cbind(1,do.call(cbind,X0)),"X_fixed",path)
    cat("R preparation:",name,"\n")
  }
}
run()
