#!/usr/bin/env Rscript
run <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=4L) stop("Expected inputs, fitting outputs, report path and unit-test log")
  inputs <- normalizePath(args[1]); out <- normalizePath(args[2]); report <- args[3]
  if(file.exists(report)) stop("The report already exists; choose a new report path")
  if(!file.exists(args[4])) stop("The Julia unit-test log is required")
  rdir <- file.path(out,"r_reference"); jdir <- file.path(out,"julia")
  r <- read.delim(file.path(rdir,"fit_summary.tsv"),check.names=FALSE)
  j <- read.delim(file.path(jdir,"fit_summary.tsv"),check.names=FALSE)
  j$converged <- as.logical(j$converged)
  stopifnot(!anyNA(j$converged))
  cases <- read.delim(file.path(inputs,"cases.tsv"),stringsAsFactors=FALSE)
  stopifnot(identical(r$case,cases$case),identical(j$case,cases$case))
  checks <- list(); index <- 0L
  scaled <- function(a,b) sqrt(sum((a-b)^2))/max(1,sqrt(sum(a*a)),sqrt(sum(b*b)))
  relative <- function(a,b) {
    if(is.infinite(a)&&is.infinite(b)&&sign(a)==sign(b)) return(0)
    if(!is.finite(a)||!is.finite(b)) return(Inf)
    abs(a-b)/max(abs(a),abs(b),.Machine$double.eps)
  }
  check <- function(case,label,error,tolerance) {
    index <<- index+1L
    checks[[index]] <<- data.frame(case=case,check=label,error=error,
      tolerance=tolerance,passed=is.finite(error)&&error<=tolerance)
  }
  source <- read.delim(file.path(rdir,"source_checksums.tsv"))
  for(i in seq_len(nrow(cases))) {
    name <- cases$case[i]; src <- file.path(inputs,name)
    # Fingerprints verify that source files did not change during fitting.
    stopifnot(source$data_md5[i]==unname(tools::md5sum(file.path(src,"data.csv"))),
              source$knots_md5[i]==unname(tools::md5sum(file.path(src,"knots.csv"))))
    rf <- read.csv(file.path(rdir,name,"fitted.csv"))
    jf <- read.csv(file.path(jdir,name,"fitted.csv"))
    raw <- read.csv(file.path(src,"data.csv"))
    stopifnot(nrow(rf)==nrow(raw),nrow(jf)==nrow(raw))
    observations_error <- max(scaled(as.matrix(rf[,c("x","y")]),as.matrix(raw[,c("x","y")])),
                              scaled(as.matrix(jf[,c("x","y")]),as.matrix(raw[,c("x","y")])))
    rk <- scan(file.path(rdir,name,"knots.csv"),quiet=TRUE)
    jk <- scan(file.path(jdir,name,"knots.csv"),quiet=TRUE)
    stopifnot(length(rk)==length(jk))
    check(name,"same predictor and response values",observations_error,1e-12)
    check(name,"same constructed knots",scaled(rk,jk),1e-12)
    check(name,"R REML converged",as.numeric(!r$converged[i]||r$method[i]!="REML"),0)
    check(name,"Julia REML converged",as.numeric(!j$converged[i]||j$method[i]!="REML"),0)
    one_variance <- r$variance_parameters[i]==1L && j$variance_parameters[i]==1L &&
      r$random_directions[i]==cases$k[i]-2L && j$random_directions[i]==cases$k[i]-2L &&
      r$shared_covariance_error[i]<1e-12
    check(name,"one shared spline variance",as.numeric(!one_variance),0)
    check(name,"R variance-to-smoothing relationship",relative(r$lambda[i],r$sigma2[i]/r$r2[i]),1e-10)
    check(name,"estimated lambda",relative(r$lambda[i],j$lambda[i]),1e-5)
    check(name,"estimated residual variance",relative(r$sigma2[i],j$sigma2[i]),1e-5)
    check(name,"estimated spline variance",relative(r$r2[i],j$r2[i]),1e-5)
    check(name,"fitted responses",scaled(rf$fitted,jf$fitted),1e-6)
  }
  checks <- do.call(rbind,checks)
  dir.create(file.path(out,"comparison"),showWarnings=FALSE)
  write.table(checks,file.path(out,"comparison","checks.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  comparison <- data.frame(case=cases$case,R_lambda=r$lambda,Julia_lambda=j$lambda,
    lambda_relative_error=mapply(relative,r$lambda,j$lambda),
    R_sigma2=r$sigma2,Julia_sigma2=j$sigma2,
    sigma2_relative_error=mapply(relative,r$sigma2,j$sigma2),
    R_r2=r$r2,Julia_r2=j$r2,
    r2_relative_error=mapply(relative,r$r2,j$r2),
    fitted_scaled_error=checks$error[checks$check=="fitted responses"])
  write.table(comparison,file.path(out,"comparison","estimates.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

  # Scientific comparison in the original observation coordinates.
  rf <- read.csv(file.path(rdir,"uniform_k8","fitted.csv"))
  jf <- read.csv(file.path(jdir,"uniform_k8","fitted.csv"))
  png(file.path(out,"comparison","fitted_curves.png"),width=1400,height=720,res=150)
  par(mar=c(4.3,4.5,2.8,1.1),cex=1)
  plot(rf$x,rf$y,pch=16,col="#B1BBC780",xlab="Predictor x",ylab="Response y",
       main="REML fit: R and Julia",cex.main=1.15)
  lines(rf$x,rf$truth,col="#757575",lty=3,lwd=2)
  lines(rf$x,rf$fitted,col="#1675AE",lwd=3)
  lines(jf$x,jf$fitted,col="#D75E1D",lty=2,lwd=2)
  legend("topleft",c("Observations","Known curve","R: mgcv / nlme","Julia: MixedModels"),
    col=c("#B1BBC7","#757575","#1675AE","#D75E1D"),
    pch=c(16,NA,NA,NA),lty=c(NA,3,1,2),lwd=c(NA,2,3,2),bty="n",cex=.85)
  dev.off()

  versions <- rbind(read.delim(file.path(rdir,"versions.tsv")),
                    read.delim(file.path(jdir,"versions.tsv")))
  version <- function(key) paste(unique(versions$value[versions$key==key]),collapse=", ")
  # Relative links make both the published and custom-run reports navigable.
  relative_path <- function(target) {
    base <- strsplit(normalizePath(dirname(report)),"/",fixed=TRUE)[[1]]
    dest <- strsplit(normalizePath(target),"/",fixed=TRUE)[[1]]
    n <- 0L
    while(n<min(length(base),length(dest)) && base[n+1L]==dest[n+1L]) n <- n+1L
    paste(c(rep("..",length(base)-n),dest[seq_along(dest)>n]),collapse="/")
  }
  link <- function(label,path) paste0("[",label,"](",relative_path(path),")")
  fmt <- function(x) sprintf("%.9g",x)
  lines <- c("# Stage 2 results","",
    "Both branches repeat Stage 1 preparation from the same saved data, then fit with **REML**. R uses `mgcv::gamm` with `nlme::lme`; Julia uses prepared `FeTerm` and `ReMat` inputs with `MixedModels.LinearMixedModel`.","",
    sprintf("**%d/%d reference comparisons passed across five cases.** No smoothing parameter is fixed before fitting.",sum(checks$passed),nrow(checks)),"",
    "## Estimated smoothing and fitted curves","",
    "Each fit estimates $\\sigma^2$ and $r^2$; $\\lambda=\\sigma^2/r^2$ is derived from the same fit. Parameter symbols follow the derivation and denote fitted estimates here.","",
    "The table compares smoothing strength $\\lambda$ and all 200 fitted responses, $\\hat y_R$ and $\\hat y_{\\mathrm{Julia}}$, in their original observation order. Matching $\\lambda$ checks smoothing; matching $\\hat y$ checks the resulting curve. Acceptance thresholds are $10^{-5}$ for the relative error in $\\lambda$ and $10^{-6}$ for the scaled fitted-response error, as defined in VALIDATION.md.","",
    "| Case | R $\\lambda$ | Julia $\\lambda$ | $\\lambda$ relative error | Fitted-response scaled error |",
    "| --- | ---: | ---: | ---: | ---: |")
  for(i in seq_len(nrow(comparison))) lines <- c(lines,sprintf("| `%s` | %s | %s | %.2e | %.2e |",
    comparison$case[i],fmt(comparison$R_lambda[i]),fmt(comparison$Julia_lambda[i]),
    comparison$lambda_relative_error[i],comparison$fitted_scaled_error[i]))
  lines <- c(lines,"",
    paste0("![R and Julia fitted curves](",relative_path(file.path(out,"comparison","fitted_curves.png")),")"),"",
    "Gray points are the observed $y$ values supplied to both fitters; the dotted curve is the known simulation curve. The blue R and orange Julia curves show their fitted responses $\\hat y$ and nearly coincide.","",
    "## Supporting checks","",
    paste0("The ",link("detailed estimates",file.path(out,"comparison","estimates.tsv")),
      " retain $\\sigma^2$ and $r^2$ for every case. These check the variance scale, which is not determined by their ratio alone."),
    sprintf("Their largest relative error is **%.2e**, below the $10^{-5}$ threshold.",
      max(comparison$sigma2_relative_error,comparison$r2_relative_error)),"",
    "Independent Julia tests check that the input adapter preserves the covariance model, REML objective and fitted coefficients. They also check invariance to rotated spline directions, the zero-variance boundary and invalid inputs.","",
    paste0(link("Detailed comparisons",file.path(out,"comparison","checks.tsv")),"; ",
      link("Julia test log",args[4]),". ",
      link("R fitting status",file.path(rdir,"fit_summary.tsv"))," and ",
      link("Julia fitting status",file.path(jdir,"fit_summary.tsv"))," record the individual fits."),"",
    "## Reproducibility","",
    paste0("Versions: ",version("R"),"; mgcv ",version("mgcv"),"; nlme ",version("nlme"),
      "; Julia ",version("Julia"),"; MixedModels ",version("MixedModels"),
      "; StatsModels ",version("StatsModels"),"; BSplineKit ",version("BSplineKit"),"."),"",
    paste0("The same saved observations and five knot settings are reused from Stage 1. ",
      link("Input fingerprints",file.path(rdir,"source_checksums.tsv")),
      " confirm the R source files; both languages' saved predictor, response and knot values are compared directly. Both stages share the dependency versions pinned in the root `Manifest.toml`."))
  writeLines(lines,report)
  cat(sprintf("%d/%d Stage 2 comparisons passed\n",sum(checks$passed),nrow(checks)))
  if(!all(checks$passed)) stop("Stage 2 comparisons failed; see checks.tsv")
}
run()
