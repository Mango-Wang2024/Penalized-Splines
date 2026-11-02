#!/usr/bin/env Rscript
run <- function() {
  suppressPackageStartupMessages(library(mgcv))
  file_arg <- grep("^--file=",commandArgs(),value=TRUE)
  script <- gsub("~+~"," ",sub("^--file=","",file_arg),fixed=TRUE)
  root <- dirname(dirname(normalizePath(script)));project <- dirname(root)
  source(file.path(project,"data","shared.R"))
  args <- commandArgs(trailingOnly=TRUE)
  inputs <- if(length(args)>=1L)args[1] else file.path(project,"data","inputs")
  out <- if(length(args)>=2L)args[2] else file.path(root,"outputs","shared","r_reference")
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Output directory is not empty")
  dir.create(out,recursive=TRUE)
  checksums <- list();index <- 0L
  for(name in read.delim(file.path(inputs,"cases.tsv"))$case) {
    case <- load_case(inputs,name)
    dst <- file.path(out,name);dir.create(dst)
    write.table(case$terms,file.path(dst,"terms.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    for(j in seq_len(nrow(case$terms))) {
      key <- case$terms$predictor[j];k <- case$terms$k[j]
      path <- file.path(dst,key);dir.create(path)
      spec <- eval(parse(text=sprintf("s(%s,k=%d,bs='bs',m=c(3,2))",key,k)))
      supplied <- if(case$terms$knot_mode[j]=="default")NULL else case$knots[key]
      raw <- smoothCon(spec,case$data,knots=supplied,absorb.cons=FALSE,scale.penalty=FALSE,diagonal.penalty=FALSE,sparse.cons=0)[[1]]
      scaled <- smoothCon(spec,case$data,knots=supplied,absorb.cons=FALSE,scale.penalty=TRUE,diagonal.penalty=FALSE,sparse.cons=0)[[1]]
      prepared <- smoothCon(spec,case$data,knots=supplied,absorb.cons=TRUE,scale.penalty=TRUE,diagonal.penalty=FALSE,sparse.cons=0)[[1]]
      qrc <- attr(prepared,"qrc")
      stopifnot(inherits(qrc,"qr"),attr(prepared,"nCons")==1L)
      H <- qr.Q(qrc,complete=TRUE)[,-1,drop=FALSE]
      xe <- scan(file.path(case$path,paste0("evaluation_",key,".csv")),quiet=TRUE)
      new <- setNames(data.frame(xe),key)
      matrices <- list(knots=matrix(raw$knots,ncol=1),D2=splines::splineDesign(raw$knots,case$data[[key]],ord=4,derivs=2),
        B_raw=raw$X,S_raw=raw$S[[1]],constraint=t(raw$C),penalty_scale=matrix(scaled$S.scale[1],1,1),
        S_scaled=scaled$S[[1]],H=H,B=prepared$X,S=prepared$S[[1]],B_eval_raw=PredictMat(raw,new),B_eval=PredictMat(prepared,new))
      for(nm in names(matrices))save_matrix(matrices[[nm]],nm,path)
      saveRDS(list(raw=raw,scaled=scaled,prepared=prepared),file.path(path,"objects.rds"))
      index <- index+1L
      checksums[[index]] <- data.frame(case=name,component=key,data_md5=unname(tools::md5sum(file.path(case$path,"data.csv"))),
        knots_md5=unname(tools::md5sum(file.path(case$path,paste0("knots_",key,".csv")))))
      cat("R basis:",name,key,"\n")
    }
  }
  write.csv(do.call(rbind,checksums),file.path(out,"source_checksums.csv"),row.names=FALSE)
  write.table(data.frame(key=c("R","mgcv"),value=c(R.version.string,as.character(packageVersion("mgcv")))),file.path(out,"versions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
}
run()
