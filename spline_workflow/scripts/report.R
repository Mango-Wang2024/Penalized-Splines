#!/usr/bin/env Rscript
# Summarize the three completed checks and verify their shared input fingerprints.
run <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=5L) stop("Supply inputs, basis output, conversion output, workflow output and report path")
  inputs <- args[1];basis <- args[2];conversion <- args[3];workflow <- args[4];report <- args[5]
  if(file.exists(report)) stop("Summary already exists; choose a new run directory")
  checks <- list(
    "Basis construction only"=read.delim(file.path(basis,"julia","checks.tsv")),
    "Conversion only"=read.delim(file.path(conversion,"julia_conversion","validation.tsv")),
    "Combined workflow"=read.delim(file.path(workflow,"julia","checks.tsv")))
  for(tab in checks) {
    if(!nrow(tab)||any(!as.character(tab$passed)%in%c("true","TRUE"))) stop("A Stage 1 comparison failed")
  }
  rchecks <- read.csv(file.path(conversion,"r_reference","reference_checks.csv"))
  if(!nrow(rchecks)||any(!rchecks$passed)) stop("R reference checks failed")
  cases <- read.delim(file.path(inputs,"cases.tsv"),stringsAsFactors=FALSE)$case
  hashes <- lapply(c(basis,workflow),function(out)read.csv(file.path(out,"r_reference","source_checksums.csv"),stringsAsFactors=FALSE))
  for(h in hashes) {
    if(!setequal(h$case,cases)) stop("Reference case lists differ")
    for(name in cases) {
      row <- h[h$case==name,,drop=FALSE]
      if(nrow(row)!=1L||row$data_md5!=unname(tools::md5sum(file.path(inputs,name,"data.csv")))||
         row$knots_md5!=unname(tools::md5sum(file.path(inputs,name,"knots.csv"))))
        stop("Reference input fingerprints do not match the current shared inputs")
    }
  }
  md <- read.delim(file.path(conversion,"r_reference","metadata.tsv"),stringsAsFactors=FALSE)
  md <- setNames(md$value,md$key)
  main <- file.path(inputs,"uniform_k8")
  data_hash <- unname(tools::md5sum(file.path(main,"data.csv")))
  knots_hash <- unname(tools::md5sum(file.path(main,"knots.csv")))
  if(md[["source_data_md5"]]!=data_hash||md[["source_knots_md5"]]!=knots_hash)
    stop("Conversion reference does not use the shared main inputs")
  read_matrix <- function(path)as.matrix(read.table(path,sep=","))
  err <- function(a,b)sqrt(sum((a-b)^2))/max(1,sqrt(sum(a^2)),sqrt(sum(b^2)))
  cross <- vapply(c("B","S"),function(nm) {
    a <- read_matrix(file.path(basis,"r_reference","uniform_k8",paste0(nm,".csv")))
    b <- read_matrix(file.path(conversion,"r_reference",paste0(nm,".csv")))
    if(!identical(dim(a),dim(b))) stop("Reference dimensions differ")
    err(a,b)
  },numeric(1))
  if(any(!is.finite(cross))||any(cross>=1e-8)) stop("The basis and captured fitter references disagree")
  dir.create(dirname(report),recursive=TRUE,showWarnings=FALSE)
  relative <- function(path) {
    from <- strsplit(normalizePath(dirname(report)),"/",fixed=TRUE)[[1]]
    to <- strsplit(normalizePath(path),"/",fixed=TRUE)[[1]]
    i <- 0L
    while(i<min(length(from),length(to))&&from[i+1L]==to[i+1L]) i <- i+1L
    tail <- if(i<length(to))to[seq.int(i+1L,length(to))] else character()
    paste(c(rep("..",length(from)-i),tail),collapse="/")
  }
  lines <- c("# Stage 1 validation","",
    "Three separate checks passed for the tested natural cubic spline model. Inputs were generated once in Julia.","",
    "| Check | Passed | Largest scaled error |","| --- | --- | --- |")
  for(name in names(checks)) {
    tab <- checks[[name]]
    errors <- if("scaled_error"%in%names(tab))tab$scaled_error else tab$relative_error
    if(any(!is.finite(errors))||any(errors>=1e-8)) stop("Invalid numerical comparison result")
    lines <- c(lines,sprintf("| %s | %d/%d | %.3g |",name,nrow(tab),nrow(tab),max(errors)))
  }
  lines <- c(lines,"",sprintf("The mathematical tests completed before these comparisons. R input-capture checks passed: %d/%d. Comparison threshold: 1e-8.",nrow(rchecks),nrow(rchecks)),"",
    "## What each check establishes","",
    "1. **Basis:** R and Julia independently construct B and S from the same data and settings. No conversion is run.",
    "2. **Conversion:** both converters use the exact same R-generated B and S. Julia is compared with mgcv and the actual lme inputs.",
    "3. **Workflow:** each language constructs and converts its own spline from the shared raw inputs. Final model representations are compared.","",
    "The basis and workflow checks cover five cases. The captured lme conversion check uses the main eight-knot case.","",
    "## Shared inputs","",
    "All basis and workflow input fingerprints match the shared files. The conversion check uses the same main data and knots.",
    paste0("Main data MD5: `",data_hash,"`. Main knots MD5: `",knots_hash,"`."),
    sprintf("Cross-reference B and S scaled errors: %.3g and %.3g.",cross[1],cross[2]),"",
    "## Reports","",
    sprintf("- [Basis construction](<%s>)",relative(file.path(basis,"RESULTS.md"))),
    sprintf("- [Conversion](<%s>)",relative(file.path(conversion,"julia_conversion","RESULTS_EXPLAINED.txt"))),
    sprintf("- [Combined workflow](<%s>)",relative(file.path(workflow,"RESULTS.md"))),"",
    "## Stage 2 inputs","",
    sprintf("The [main Julia model inputs](<%s>) contain y, X_fixed and Z. Matching B, S, C, knots, centering H, penalty scale and variance assumptions are retained.",relative(file.path(workflow,"julia","uniform_k8"))),"",
    "Stage 1 prepares the model. Julia coefficient fitting, smoothing estimation and fitted predictions remain later stages.",
    "The R input-capture example performs a reference fit. Its fitted parameters are not Julia estimates.",
    paste0("Reference versions: ",md[["R_version"]],", mgcv ",md[["mgcv_version"]],", nlme ",md[["nlme_version"]],". The historical source revision is not automatically reproduced."))
  writeLines(lines,report)
  cat("All three Stage 1 checks passed. Summary:",normalizePath(report),"\n")
}
run()
