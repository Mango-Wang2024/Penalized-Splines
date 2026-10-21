#!/usr/bin/env Rscript
# Convert the R-constructed smooths for the combined-workflow reference.
run <- function() {
  if(!requireNamespace("mgcv",quietly=TRUE)) stop("Install mgcv before running this script")
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=1L) stop("Supply the R construction reference directory")
  reference <- args[1]
  cases <- read.csv(file.path(reference,"source_checksums.csv"),stringsAsFactors=FALSE)$case
  exports <- c("X_fixed","X0","Z","C")
  # Check every destination before writing any conversion output.
  for(name in cases) for(nm in exports) {
    if(file.exists(file.path(reference,name,paste0(nm,".csv"))))
      stop("Converted reference already exists; choose a new run directory")
  }
  for(name in cases) {
    dst <- file.path(reference,name)
    sm <- readRDS(file.path(dst,"objects.rds"))$prepared
    re <- mgcv::smooth2random(sm,c("x","y","truth"),type=2)
    if(length(re$rand)!=1L) stop("Expected one random block")
    U <- if(is.null(re$trans.U))diag(ncol(sm$X)) else re$trans.U
    D <- if(is.null(re$trans.D))rep(1,ncol(sm$X)) else re$trans.D
    C <- sweep(U,2,D,"*")
    if(qr(C)$rank!=ncol(sm$X)) stop("R conversion lost coefficient directions")
    matrices <- list(X_fixed=cbind(1,re$Xf),X0=re$Xf,Z=as.matrix(re$rand[[1]]),C=C)
    for(nm in names(matrices)) {
      M <- as.matrix(matrices[[nm]])
      strings <- matrix(sprintf("%.17g",as.numeric(M)),nrow(M),ncol(M))
      write.table(strings,file.path(dst,paste0(nm,".csv")),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
    }
    cat("R construction and conversion:",name,"\n")
  }
}
run()
