#!/usr/bin/env Rscript
# Independently construct reference matrices from Julia's saved predictor values.
run <- function() {
  if (!requireNamespace("mgcv",quietly=TRUE)) stop("Install mgcv before running this script")
  suppressPackageStartupMessages(library(mgcv))
  file_arg <- grep("^--file=",commandArgs(),value=TRUE)
  script_path <- gsub("~+~"," ",sub("^--file=","",file_arg),fixed=TRUE)
  root <- dirname(dirname(normalizePath(script_path)))
  args <- commandArgs(trailingOnly=TRUE)
  inputs <- if(length(args)>=1L)args[1] else file.path(dirname(root),"data","inputs")
  out <- if(length(args)>=2L)args[2] else file.path(root,"outputs","shared","r_reference")
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Reference directory is not empty; choose a new run directory")
  cases <- read.delim(file.path(inputs,"cases.tsv"),stringsAsFactors=FALSE)
  source_checksums <- data.frame(case=character(),data_md5=character(),knots_md5=character())
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  save_matrix <- function(M,name,dir) {
    M <- as.matrix(M)
    # Preserve enough decimal digits for a Float64 round trip.
    strings <- matrix(sprintf("%.17g",as.numeric(M)),nrow(M),ncol(M))
    write.table(strings,file.path(dir,paste0(name,".csv")),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
  }
  for(i in seq_len(nrow(cases))) {
    name <- cases$case[i];k <- cases$k[i]
    p <- file.path(inputs,name);dst <- file.path(out,name)
    dir.create(dst,showWarnings=FALSE)
    dat <- read.csv(file.path(p,"data.csv"))
    source_checksums <- rbind(source_checksums,data.frame(case=name,
      data_md5=unname(tools::md5sum(file.path(p,"data.csv"))),
      knots_md5=unname(tools::md5sum(file.path(p,"knots.csv")))))
    t <- scan(file.path(p,"knots.csv"),quiet=TRUE)
    xe <- scan(file.path(p,"evaluation_x.csv"),quiet=TRUE)
    supplied <- if(cases$knot_mode[i]=="default")NULL else list(x=t)
    spec <- s(x,k=k,bs="cr")
    raw <- smoothCon(spec,dat,knots=supplied,absorb.cons=FALSE,scale.penalty=FALSE,
                     diagonal.penalty=FALSE,sparse.cons=0)[[1]]
    scaled <- smoothCon(spec,dat,knots=supplied,absorb.cons=FALSE,scale.penalty=TRUE,
                        diagonal.penalty=FALSE,sparse.cons=0)[[1]]
    prepared <- smoothCon(spec,dat,knots=supplied,absorb.cons=TRUE,scale.penalty=TRUE,
                          diagonal.penalty=FALSE,sparse.cons=0)[[1]]
    qrc <- attr(prepared,"qrc")
    if(!inherits(qrc,"qr")||nrow(qrc$qr)!=k||attr(prepared,"nCons")!=1L)
      stop("This example requires one full centering constraint")
    H <- qr.Q(qrc,complete=TRUE)[,-1,drop=FALSE]
    matrices <- list(knots=matrix(raw$xp,ncol=1),F=t(matrix(raw$F,k,k)),
      B_raw=raw$X,S_raw=raw$S[[1]],constraint=t(raw$C),
      penalty_scale=matrix(scaled$S.scale[1],1,1),S_scaled=scaled$S[[1]],
      H=H,B=prepared$X,S=prepared$S[[1]],
      B_eval_raw=PredictMat(raw,data.frame(x=xe)),
      B_eval=PredictMat(prepared,data.frame(x=xe)))
    for(nm in names(matrices))save_matrix(matrices[[nm]],nm,dst)
    saveRDS(list(raw=raw,scaled=scaled,prepared=prepared),file.path(dst,"objects.rds"))
    cat("R reference:",name,"|",nrow(dat),"observations |",k,"knots\n")
  }
  metadata <- data.frame(key=c("R_version","mgcv_version","spline_type","constraint","penalty_scaling"),
    value=c(R.version.string,as.character(packageVersion("mgcv")),"cr","one mean-zero constraint","one-norm(S_raw) / infinity-norm(B_raw)^2"))
  write.table(metadata,file.path(out,"metadata.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.csv(source_checksums,file.path(out,"source_checksums.csv"),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
}
run()
