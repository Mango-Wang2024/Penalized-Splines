#!/usr/bin/env Rscript
# Reference fitting uses the existing mgcv -> nlme implementation.
# Sources: mgcv::gamm and nlme::pdMatrix official manuals (linked in VALIDATION.md).
run <- function() {
  suppressPackageStartupMessages(library(mgcv))
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Expected shared input directory and new reference output directory")
  inputs <- args[1]; out <- args[2]
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
    stop("Reference output directory is not empty")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  options(digits=17)
  cases <- read.delim(file.path(inputs,"cases.tsv"),stringsAsFactors=FALSE)
  rows <- list(); fingerprints <- list()
  for(i in seq_len(nrow(cases))) {
    name <- cases$case[i]; k <- cases$k[i]
    src <- file.path(inputs,name); dst <- file.path(out,name)
    dir.create(dst)
    dat <- read.csv(file.path(src,"data.csv"))
    knots <- scan(file.path(src,"knots.csv"),quiet=TRUE)
    supplied <- if(cases$knot_mode[i]=="default") NULL else list(x=knots)
    # Fail on a convergence warning rather than accepting a partial lme fit.
    fit <- withCallingHandlers(
      gamm(y~s(x,bs="bs",m=c(3,2),k=k),data=dat,knots=supplied,
           method="REML",control=list(niterEM=0,msMaxIter=200,
             msTol=1e-9,tolerance=1e-8,returnObject=FALSE)),
      warning=function(w) stop(conditionMessage(w),call.=FALSE))
    stopifnot(fit$lme$method=="REML",length(fit$gam$sp)==1L,
              length(fit$lme$modelStruct$reStruct)==1L)
    relative_cov <- as.matrix(fit$lme$modelStruct$reStruct[[1]])
    ratio <- mean(diag(relative_cov))
    shared_error <- max(abs(relative_cov-diag(ratio,nrow(relative_cov))))
    sigma2 <- fit$lme$sigma^2
    r2 <- sigma2*ratio
    r <- sqrt(r2)
    lambda <- unname(fit$gam$sp[1])
    fitted <- as.numeric(predict(fit$gam,newdata=dat,type="response"))
    stopifnot(all(is.finite(c(lambda,sigma2,r2,fitted))),sigma2>0,
              r2>0,nrow(relative_cov)==k-2L,shared_error<1e-12,
              abs(lambda-sigma2/r2)<1e-10*max(1,lambda),
              max(abs(fitted-as.numeric(stats::fitted(fit$lme))))<1e-9)
    write.csv(data.frame(x=dat$x,y=dat$y,truth=dat$truth,fitted=fitted),
              file.path(dst,"fitted.csv"),row.names=FALSE)
    write.table(matrix(fit$gam$smooth[[1]]$knots,ncol=1),file.path(dst,"knots.csv"),
                sep=",",row.names=FALSE,col.names=FALSE)
    saveRDS(fit,file.path(dst,"model.rds"))
    rows[[i]] <- data.frame(case=name,lambda=lambda,sigma2=sigma2,r2=r2,r=r,
      random_directions=nrow(relative_cov),variance_parameters=1L,
      shared_covariance_error=shared_error,converged=TRUE,method="REML",
      status="lme returned without warnings")
    fingerprints[[i]] <- data.frame(case=name,
      data_md5=unname(tools::md5sum(file.path(src,"data.csv"))),
      knots_md5=unname(tools::md5sum(file.path(src,"knots.csv"))))
    cat("R REML:",name,"| lambda =",lambda,"\n")
  }
  write.table(do.call(rbind,rows),file.path(out,"fit_summary.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,fingerprints),file.path(out,"source_checksums.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  write.table(data.frame(key=c("R","mgcv","nlme","method"),
    value=c(R.version.string,as.character(packageVersion("mgcv")),
      as.character(packageVersion("nlme")),"REML")),file.path(out,"versions.tsv"),
    sep="\t",quote=FALSE,row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
}
run()
