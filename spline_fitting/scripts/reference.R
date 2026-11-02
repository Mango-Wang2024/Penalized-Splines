#!/usr/bin/env Rscript
run <- function() {
  suppressPackageStartupMessages(library(mgcv))
  script <- gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)),fixed=TRUE)
  project <- dirname(dirname(dirname(normalizePath(script))))
  source(file.path(project,"data","shared.R"))
  args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
  inputs <- args[1];out <- args[2]
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty")
  dir.create(out,recursive=TRUE)
  rows <- list();estimates <- list();hashes <- list();index <- 0L;hi <- 0L
  for(name in read.delim(file.path(inputs,"cases.tsv"))$case) {
    case <- load_case(inputs,name);dst <- file.path(out,name);dir.create(dst)
    fit <- withCallingHandlers(gamm(case$formula,data=case$data,knots=case$supplied,method="REML",
      control=list(niterEM=0,msMaxIter=200,msMaxEval=1000,msTol=1e-13,tolerance=1e-10,
                   rel.tol=1e-14,x.tol=1e-12,sing.tol=1e-25,returnObject=FALSE)),
      warning=function(w)stop(conditionMessage(w),call.=FALSE))
    m <- nrow(case$terms);sigma2 <- fit$lme$sigma^2
    stopifnot(fit$lme$method=="REML",length(fit$gam$sp)==m,length(fit$lme$modelStruct$reStruct)==m)
    curves <- predict(fit$gam,newdata=case$data,type="terms")
    colnames(curves) <- case$terms$predictor
    fitted <- as.numeric(predict(fit$gam,newdata=case$data,type="response"))
    stopifnot(max(abs(fitted-as.numeric(stats::fitted(fit$lme))))<1e-8)
    write.csv(cbind(case$data[,case$terms$predictor,drop=FALSE],y=case$data$y,truth=case$data$truth,fitted=fitted),file.path(dst,"fitted.csv"),row.names=FALSE)
    write.csv(curves,file.path(dst,"components.csv"),row.names=FALSE)
    grid <- grid_inputs(case)
    write.csv(cbind(grid,fitted=as.numeric(predict(fit$gam,newdata=grid,type="response"))),file.path(dst,"grid.csv"),row.names=FALSE)
    covariance <- lapply(fit$lme$modelStruct$reStruct,as.matrix)
    unused <- seq_along(covariance)
    for(j in seq_len(m)) {
      key <- case$terms$predictor[j];lambda <- unname(fit$gam$sp[j]);r2 <- sigma2/lambda
      q <- case$terms$k[j]-2L
      errors <- vapply(unused,function(i) {
        A <- covariance[[i]]
        if(nrow(A)!=q)return(Inf)
        max(abs(A-diag(1/lambda,q)))/max(1,1/lambda)
      },0.)
      selected <- which.min(errors);stopifnot(errors[selected]<1e-8)
      unused <- unused[-selected]
      index <- index+1L
      estimates[[index]] <- data.frame(case=name,component=key,lambda=lambda,r2=r2,r=sqrt(r2),theta=sqrt(1/lambda),random_directions=q)
      path <- file.path(dst,key);dir.create(path)
      save_matrix(fit$gam$smooth[[j]]$knots,"knots",path)
    }
    rows[[name]] <- data.frame(case=name,components=m,sigma2=sigma2,fixed_columns=length(fit$lme$coefficients$fixed),
      random_directions=sum(case$terms$k-2L),converged=TRUE,method="REML",status="lme returned without warnings")
    for(file in c("data.csv","terms.tsv","grid.csv",paste0("knots_",case$terms$predictor,".csv"))) {
      hi <- hi+1L;hashes[[hi]] <- data.frame(case=name,file=file,md5=unname(tools::md5sum(file.path(case$path,file))))
    }
    saveRDS(fit,file.path(dst,"model.rds"))
    cat("R joint REML:",name,"|",m,"components\n")
  }
  write.table(do.call(rbind,rows),file.path(out,"fit_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,estimates),file.path(out,"component_estimates.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,hashes),file.path(out,"source_checksums.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(data.frame(key=c("R","mgcv","nlme","method"),value=c(R.version.string,as.character(packageVersion("mgcv")),as.character(packageVersion("nlme")),"REML")),file.path(out,"versions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}
run()
