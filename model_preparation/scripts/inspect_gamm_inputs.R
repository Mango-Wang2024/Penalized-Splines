#!/usr/bin/env Rscript
# Capture the actual prepared arguments at the mgcv -> nlme boundary.
run <- function() {
  suppressPackageStartupMessages(library(mgcv))
  script <- gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)),fixed=TRUE)
  project <- dirname(dirname(dirname(normalizePath(script))))
  source(file.path(project,"data","shared.R"))
  args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L)
  out <- args[1];case <- load_case(args[2],args[3]);m <- nrow(case$terms)
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty")
  dir.create(out,recursive=TRUE)
  captured <- new.env(parent=emptyenv());wrapped <- mgcv::gamm
  env <- new.env(parent=environment(wrapped))
  env$lme <- function(fixed,random,data,correlation=NULL,control=list(),weights=NULL,method=c("REML","ML"),...) {
    captured$inputs <- list(fixed=fixed,random=random,data=data,correlation=correlation,weights=weights,method=method,control=control)
    nlme::lme(fixed=fixed,random=random,data=data,correlation=correlation,weights=weights,method=method,control=control,...)
  }
  environment(wrapped) <- env
  fit <- wrapped(case$formula,data=case$data,knots=case$supplied,method="REML")
  baseline <- mgcv::gamm(case$formula,data=case$data,knots=case$supplied,method="REML")
  input <- captured$inputs
  stopifnot(!is.null(input),length(input$random)==m)
  rel <- function(a,b)sqrt(sum((a-b)^2))/max(1,sqrt(sum(a*a)),sqrt(sum(b*b)))
  projector <- function(A) {v <- svd(A);r <- sum(v$d>max(v$d)*1e-10);tcrossprod(v$u[,seq_len(r),drop=FALSE])}
  records <- list()
  check <- function(component,label,error) records[[length(records)+1L]] <<- data.frame(component=component,check=label,relative_error=error,passed=is.finite(error)&&error<1e-8)
  check("all","wrapped vs normal predictions",rel(predict(fit$gam),predict(baseline$gam)))
  check("all","wrapped vs normal smoothing",rel(fit$gam$sp,baseline$gam$sp))
  mf <- model.frame(input$fixed,data=input$data)
  check("all","original response passed to lme",rel(model.response(mf),case$data$y))
  Xactual <- model.matrix(input$fixed,data=mf)
  Zactual <- lapply(input$random,function(term)model.matrix(formula(term),data=input$data))
  unused <- seq_along(Zactual);X0 <- list()
  for(j in seq_len(m)) {
    key <- case$terms$predictor[j];sm <- fit$gam$smooth[[j]]
    B <- PredictMat(sm,case$data);S <- sm$S[[1]]
    sm$X <- B;sm$fixed <- FALSE
    re <- smooth2random(sm,names(case$data),type=2)
    U <- if(is.null(re$trans.U))diag(ncol(B)) else re$trans.U
    D <- if(is.null(re$trans.D))rep(1,ncol(B)) else re$trans.D
    C <- sweep(U,2,D,"*");Z <- as.matrix(re$rand[[1]]);X0[[key]] <- re$Xf
    target <- diag(c(rep(1,ncol(Z)),rep(0,ncol(re$Xf))))
    check(key,"curve identity",rel(B%*%C,cbind(Z,re$Xf)))
    check(key,"normalized penalty",rel(t(C)%*%S%*%C,target))
    errors <- vapply(unused,function(i)rel(tcrossprod(Z),tcrossprod(Zactual[[i]])),0.)
    best <- which.min(errors)
    check(key,"random covariance passed to lme",errors[best]);unused <- unused[-best]
    dst <- file.path(out,key);dir.create(dst)
    for(nm in c("B","S","X0","Z","C"))save_matrix(switch(nm,B=B,S=S,X0=re$Xf,Z=Z,C=C),nm,dst)
  }
  Xfixed <- cbind(1,do.call(cbind,X0))
  check("all","combined fixed space passed to lme",rel(projector(Xfixed),projector(Xactual)))
  checks <- do.call(rbind,records);stopifnot(all(checks$passed))
  save_matrix(Xactual,"X_fixed_lme",out)
  write.csv(checks,file.path(out,"reference_checks.csv"),row.names=FALSE)
  metadata <- data.frame(key=c("R_version","mgcv_version","nlme_version","observations","components"),
    value=c(R.version.string,as.character(packageVersion("mgcv")),as.character(packageVersion("nlme")),nrow(case$data),m))
  write.table(metadata,file.path(out,"metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  report <- c("WHAT DOES MGCV PASS TO LME?","",
    paste("Data case:",case$name,"| observations:",nrow(case$data),"| spline components:",m),
    paste("Fixed formula:",paste(deparse(input$fixed),collapse=" ")),
    paste("Fixed design:",paste(dim(Xactual),collapse=" x ")),
    paste("Random designs:",paste(vapply(Zactual,function(Z)paste(dim(Z),collapse=" x "),""),collapse="; ")),
    "Each spline block has its own shared coefficient variance; residuals share sigma squared.",
    "The response is unchanged at the lme input boundary. Rotations occur internally during fitting.",
    paste("Checks:",sum(checks$passed),"/",nrow(checks),"passed."),"",
    "ACTUAL RANDOM FORMULAS",vapply(input$random,function(term)paste(deparse(formula(term)),collapse=" "),""),"",
    "R controls only this captured reference; Julia independently constructs its own inputs.")
  writeLines(report,file.path(out,"INPUTS_EXPLAINED.txt"))
  saveRDS(input,file.path(out,"lme_inputs.rds"))
  cat(paste(report,collapse="\n"),"\n")
}
run()
