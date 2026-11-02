#!/usr/bin/env Rscript
run <- function() {
  args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==4L)
  inputs <- normalizePath(args[1]);basis <- normalizePath(args[2]);workflow <- normalizePath(args[3]);report <- args[4]
  if(file.exists(report))stop("Report already exists")
  source(file.path(dirname(inputs),"shared.R"))
  tabs <- list("Basis construction only"=read.delim(file.path(basis,"julia","checks.tsv")),
               "Combined additive preparation"=read.delim(file.path(workflow,"julia","checks.tsv")))
  for(tab in tabs)stopifnot(all(as.character(tab$passed)%in%c("true","TRUE")),all(is.finite(tab$scaled_error)),all(tab$scaled_error<1e-8))
  captured <- file.path(workflow,"r_reference","lme_inputs")
  rc <- read.csv(file.path(captured,"reference_checks.csv"));stopifnot(all(rc$passed))
  for(path in c(basis,workflow)) {
    hashes <- read.csv(file.path(path,"r_reference","source_checksums.csv"))
    for(i in seq_len(nrow(hashes))) {
      src <- file.path(inputs,hashes$case[i])
      stopifnot(hashes$data_md5[i]==unname(tools::md5sum(file.path(src,"data.csv"))),
        hashes$knots_md5[i]==unname(tools::md5sum(file.path(src,paste0("knots_",hashes$component[i],".csv")))))
    }
  }
  cross <- c()
  main <- load_case(inputs,"uniform_k8")
  for(key in main$terms$predictor)for(nm in c("B","S")) {
    a <- as.matrix(read.csv(file.path(basis,"r_reference","uniform_k8",key,paste0(nm,".csv")),header=FALSE))
    b <- as.matrix(read.csv(file.path(captured,key,paste0(nm,".csv")),header=FALSE))
    cross <- c(cross,sqrt(sum((a-b)^2))/max(1,sqrt(sum(a*a)),sqrt(sum(b*b))))
  }
  stopifnot(max(cross)<1e-8)
  relative <- function(path) {
    a <- strsplit(normalizePath(dirname(report)),"/",fixed=TRUE)[[1]];b <- strsplit(normalizePath(path),"/",fixed=TRUE)[[1]]
    n <- 0L;while(n<min(length(a),length(b))&&a[n+1L]==b[n+1L])n <- n+1L
    paste(c(rep("..",length(a)-n),b[seq_along(b)>n]),collapse="/")
  }
  lines <- c("# Stage 1 results","","R and Julia independently construct each component and assemble the additive model from the same saved inputs.","",
    "## Overall results","","| Check | Passed | Largest scaled error |","| --- | ---: | ---: |")
  for(name in names(tabs)) {
    tab <- tabs[[name]]
    lines <- c(lines,sprintf("| %s | %d/%d | %.3g |",name,nrow(tab),nrow(tab),max(tab$scaled_error)))
  }
  lines <- c(lines,"",sprintf("Captured R fitting-input checks: **%d/%d passed**. All preparation errors are below $10^{-8}$.",nrow(rc),nrow(rc)),
    sprintf("The basis references and captured fitting inputs agree; their largest scaled difference is %.3g.",max(cross)),"",
    "## Results by case","","| Case | Components | Basis checks | Preparation checks | Largest scaled error |","| --- | ---: | ---: | ---: | ---: |")
  for(name in read.delim(file.path(inputs,"cases.tsv"))$case) {
    case <- load_case(inputs,name);a <- tabs[[1]][tabs[[1]]$case==name,];b <- tabs[[2]][tabs[[2]]$case==name,]
    stopifnot(nrow(a)==13L*nrow(case$terms),nrow(b)==1L+6L*nrow(case$terms))
    lines <- c(lines,sprintf("| %s | %d | %d/%d | %d/%d | %.3g |",name,nrow(case$terms),nrow(a),nrow(a),nrow(b),nrow(b),max(a$scaled_error,b$scaled_error)))
  }
  lines <- c(lines,"","## Prepared model","",
    "The main two-predictor case produces $X_{\\mathrm{fixed}}$ with dimensions 200 × 3 and two penalized designs $Z_1,Z_2$, each 200 × 6. Each component retains its own basis, penalty, knots and coefficient transformation.","",
    "## Detailed evidence","",
    paste0("- [Basis comparisons](",relative(file.path(basis,"julia","checks.tsv")),")"),
    paste0("- [Additive preparation comparisons](",relative(file.path(workflow,"julia","checks.tsv")),")"),
    paste0("- [Captured R inputs](",relative(file.path(captured,"INPUTS_EXPLAINED.txt")),")"),
    paste0("- [Captured-input checks](",relative(file.path(captured,"reference_checks.csv")),")"),"",
    "Sources: [mgcv construction](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smoothCon.html) and [mixed-model preparation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html).")
  dir.create(dirname(report),recursive=TRUE,showWarnings=FALSE);writeLines(lines,report)
  cat("Stage 1 reports generated; all comparisons passed.\n")
}
run()
