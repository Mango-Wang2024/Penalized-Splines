#!/usr/bin/env Rscript
# Summarize construction and the combined workflow, checking shared inputs.
run <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=4L) stop("Supply inputs, basis output, workflow output and report path")
  inputs <- args[1];basis <- args[2];workflow <- args[3];report <- args[4]
  if(file.exists(report)) stop("Summary already exists; choose a new run directory")
  checks <- list(
    "Basis construction only"=read.delim(file.path(basis,"julia","checks.tsv")),
    "Combined workflow"=read.delim(file.path(workflow,"julia","checks.tsv")))
  for(tab in checks) {
    if(!nrow(tab)||any(!as.character(tab$passed)%in%c("true","TRUE"))) stop("A Stage 1 comparison failed")
    if(any(!is.finite(tab$scaled_error))||any(tab$scaled_error>=1e-8)) stop("Invalid numerical comparison result")
  }
  captured <- file.path(workflow,"r_reference","lme_inputs")
  rchecks <- read.csv(file.path(captured,"reference_checks.csv"))
  if(!nrow(rchecks)||any(!rchecks$passed)||any(!is.finite(rchecks$relative_error))||
     any(rchecks$relative_error>=1e-8)) stop("R fitting-input checks failed")
  settings <- read.delim(file.path(inputs,"cases.tsv"),stringsAsFactors=FALSE)
  cases <- settings$case
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
  md <- read.delim(file.path(captured,"metadata.tsv"),stringsAsFactors=FALSE)
  md <- setNames(md$value,md$key)
  main <- file.path(inputs,"uniform_k8")
  data_hash <- unname(tools::md5sum(file.path(main,"data.csv")))
  knots_hash <- unname(tools::md5sum(file.path(main,"knots.csv")))
  if(md[["source_data_md5"]]!=data_hash||md[["source_knots_md5"]]!=knots_hash)
    stop("Captured fitter reference does not use the shared main inputs")
  read_matrix <- function(path)as.matrix(read.table(path,sep=","))
  err <- function(a,b)sqrt(sum((a-b)^2))/max(1,sqrt(sum(a^2)),sqrt(sum(b^2)))
  cross <- vapply(c("B","S"),function(nm) {
    a <- read_matrix(file.path(basis,"r_reference","uniform_k8",paste0(nm,".csv")))
    b <- read_matrix(file.path(captured,paste0(nm,".csv")))
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
  lines <- c("# Stage 1 results","",
    "Saved comparisons for cubic B-splines with one curvature penalty: mgcv bs=\"bs\", m=c(3,2). All comparisons use shared Julia-generated inputs.","",
    "Stage 1 checks basis construction and the combined workflow. The required fixed/random preparation is part of that workflow.","",
    "## Overall results","",
    "| Check | Passed | Largest scaled error |","| --- | --- | --- |")
  for(name in names(checks)) {
    tab <- checks[[name]]
    lines <- c(lines,sprintf("| %s | %d/%d | %.3g |",name,nrow(tab),nrow(tab),max(tab$scaled_error)))
  }
  lines <- c(lines,"",sprintf("Mathematical tests completed before these comparisons. R fitting-input checks passed: %d/%d. Comparison threshold: 1e-8.",nrow(rchecks),nrow(rchecks)))
  for(name in names(checks)) {
    tab <- checks[[name]]
    if(!setequal(unique(tab$case),cases)) stop("Comparison cases do not match the shared settings")
    count_per_case <- if(name=="Basis construction only")13L else 7L
    lines <- c(lines,"",paste0("## ",name,": results by case"),"",
      sprintf("Five test cases, with %d comparisons per case.",count_per_case),"",
      "| Case | Raw basis size k | Passed | Largest scaled error |",
      "| --- | --- | --- | --- |")
    for(case in cases) {
      rows <- tab[tab$case==case,,drop=FALSE]
      if(nrow(rows)!=count_per_case) stop("Unexpected number of comparisons for a case")
      k <- settings$k[settings$case==case]
      lines <- c(lines,sprintf("| %s | %d | %d/%d | %.3g |",case,k,nrow(rows),nrow(rows),max(rows$scaled_error)))
    }
  }
  lines <- c(lines,"","## Shared inputs","",
    "Data and knot fingerprints match the shared files for all five cases. The captured R fitting inputs use the same main observations and knots.",
    paste0("Main data MD5: `",data_hash,"`. Main knots MD5: `",knots_hash,"`."),
    sprintf("Cross-reference B and S scaled errors: %.3g and %.3g.",cross[1],cross[2]),"",
    "## Detailed evidence","",
    sprintf("- [Basis construction checks](<%s>)",relative(file.path(basis,"julia","checks.tsv"))),
    sprintf("- [Combined workflow checks](<%s>)",relative(file.path(workflow,"julia","checks.tsv"))),
    sprintf("- [R fitting-input checks](<%s>)",relative(file.path(captured,"reference_checks.csv"))),"",
    "## Prepared model","",
    sprintf("The [main Julia model inputs](<%s>) contain y, X_fixed and Z, together with B, S, the coefficient map C, training settings and variance assumptions.",relative(file.path(workflow,"julia","uniform_k8"))),"",
    paste0("R reference versions: ",md[["R_version"]],", mgcv ",md[["mgcv_version"]],", nlme ",md[["nlme_version"]],"."),"",
    "## Sources","",
    "- [mgcv cubic B-splines and curvature penalty](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth.construct.bs.smooth.spec.html)",
    "- [Basis dimensions](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/splineDesign.html)",
    "- [Fixed/random preparation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/smooth2random.html)")
  writeLines(lines,report)
  cat("Both Stage 1 comparisons passed. Summary:",normalizePath(report),"\n")
}
run()
