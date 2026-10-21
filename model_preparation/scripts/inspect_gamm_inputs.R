#!/usr/bin/env Rscript
# Stage 1: capture the actual lme inputs produced by a simple mgcv::gamm fit.
# Uses the installed mgcv; it does NOT silently install the historical revision.
# The temporary local wrapper does not modify installed packages or global bindings.

run <- function() {
  for (pkg in c("mgcv","nlme")) {
    if (!requireNamespace(pkg,quietly=TRUE))
      stop(sprintf('Missing required R package: %s',pkg))
  }
  suppressPackageStartupMessages(library(mgcv))
  file_arg <- grep("^--file=",commandArgs(),value=TRUE)
  if (length(file_arg)!=1L) stop("Run this script with Rscript.")
  script_path <- gsub("~+~"," ",sub("^--file=","",file_arg),fixed=TRUE)
  root <- dirname(dirname(normalizePath(script_path)))
  project <- dirname(root)
  args <- commandArgs(trailingOnly=TRUE)
  out <- if (length(args)>=1L) args[1] else file.path(root,"outputs","shared","r_reference","lme_inputs")
  data_csv <- if(length(args)>=2L)args[2] else file.path(project,"data","inputs","uniform_k8","data.csv")
  knots_csv <- if(length(args)>=3L)args[3] else file.path(project,"data","inputs","uniform_k8","knots.csv")
  if(!file.exists(data_csv)||!file.exists(knots_csv)) stop("Shared data are missing. Run data/generate.jl from the project root.")
  dat <- read.csv(data_csv)
  if(!all(c("x","y")%in%names(dat))||!is.numeric(dat$x)||!is.numeric(dat$y)||any(!is.finite(dat$x))||any(!is.finite(dat$y)))
    stop("Data must contain finite numeric x and y columns")
  knots <- scan(knots_csv,quiet=TRUE)
  n <- nrow(dat);k <- length(knots)-4
  if(n<20L||k<5L||k>=length(unique(dat$x))||any(!is.finite(knots))||any(diff(knots)<=0)||min(dat$x)<knots[4]||max(dat$x)>knots[length(knots)-3])
    stop("Require at least 20 observations and 5<=k<unique(x), with increasing finite knots")
  # Refuse to overwrite a previous reference silently.
  if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
    stop("Output directory is not empty. Supply a new directory for a new run.")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  f <- as.formula(sprintf("y ~ s(x, bs='bs', m=c(3,2), k=%d)",k),env=environment())

  # Clone gamm's closure and shadow only its unqualified lme lookup. This is
  # the lookup used in the reference source. All other helpers resolve
  # through the original mgcv namespace. A normal gamm fit cross-checks the wrapper.
  captured <- new.env(parent=emptyenv())
  reference_gamm <- mgcv::gamm
  wrapper_env <- new.env(parent=environment(reference_gamm))
  wrapper_env$lme <- function(fixed,random,data,correlation=NULL,control=list(),
                             weights=NULL,method=c("REML","ML"),...) {
    captured$inputs <- list(fixed=fixed,random=random,data=data,
                           correlation=correlation,control=control,weights=weights,
                           method=method,extra_arguments=list(...))
    nlme::lme(fixed=fixed,random=random,data=data,correlation=correlation,
              control=control,weights=weights,method=method,...)
  }
  environment(reference_gamm) <- wrapper_env
  fit <- reference_gamm(f,data=dat,knots=list(x=knots),method="REML",correlation=NULL,weights=NULL)
  if (is.null(captured$inputs))
    stop("Installed gamm did not call the intercepted lme lookup. Inspect its source before adapting the capture.")
  baseline <- mgcv::gamm(f,data=dat,knots=list(x=knots),method="REML",correlation=NULL,weights=NULL)
  rel <- function(a,b) sqrt(sum((a-b)^2))/max(1,sqrt(sum(a^2)),sqrt(sum(b^2)))
  wrapper_prediction_error <- rel(predict(fit$gam),predict(baseline$gam))
  wrapper_sp_error <- rel(fit$gam$sp,baseline$gam$sp)
  if (wrapper_prediction_error>1e-8 || wrapper_sp_error>1e-8)
    stop("Captured and normal gamm fits disagree. Do not use this reference.")

  input <- captured$inputs
  if (length(input$random)!=1L) stop("This script supports exactly one spline random block.")
  sm <- fit$gam$smooth[[1]]
  if (length(sm$S)!=1L) stop("This script supports exactly one penalty.")
  B <- mgcv::PredictMat(sm,dat)
  S <- sm$S[[1]]
  if (!identical(dim(S),c(ncol(B),ncol(B)))) stop("Basis/penalty dimensions disagree.")
  # Give smooth2random the original basis, as gamm.setup does before conversion.
  sm_for_conversion <- sm
  sm_for_conversion$X <- B
  sm_for_conversion$fixed <- FALSE
  re <- mgcv::smooth2random(sm_for_conversion,names(dat),type=2)
  if (length(re$rand)!=1L) stop("Expected one converted random-effect matrix.")
  X0 <- re$Xf
  Z <- as.matrix(re$rand[[1]])
  U <- if (is.null(re$trans.U)) diag(ncol(B)) else re$trans.U
  D <- if (is.null(re$trans.D)) rep(1,ncol(B)) else re$trans.D
  C <- sweep(U,2,D,"*") # z = C %*% c(r, alpha), for this single-penalty smooth.
  r <- ncol(Z)
  Xordinary <- matrix(1,n,1,dimnames=list(NULL,"intercept"))
  Xfixed <- cbind(Xordinary,X0)
  mf <- model.frame(input$fixed,data=input$data)
  y_lme <- model.response(mf)
  X_lme <- model.matrix(input$fixed,data=mf)
  Z_lme <- model.matrix(formula(input$random[[1]]),data=input$data)
  group_name <- names(input$random)[1]
  group <- input$data[[group_name]]
  if (is.null(group)||length(unique(group))!=1L) stop("Expected one dummy spline group.")

  projector <- function(A) {
    v <- svd(A)
    q <- sum(v$d>max(v$d)*1e-10)
    if (!q) return(matrix(0,nrow(A),nrow(A)))
    tcrossprod(v$u[,seq_len(q),drop=FALSE])
  }
  expected_penalty <- diag(c(rep(1,r),rep(0,ncol(X0))))
  checks <- data.frame(
    check=c("wrapper predictions vs normal gamm","wrapper smoothing vs normal gamm",
            "response passed to lme","curve identity","normalized penalty",
            "fixed space: smooth2random vs actual lme","random covariance: smooth2random vs actual lme"),
    relative_error=c(wrapper_prediction_error,wrapper_sp_error,rel(y_lme,dat$y),
                     rel(B%*%C,cbind(Z,X0)),rel(t(C)%*%S%*%C,expected_penalty),
                     rel(projector(Xfixed),projector(X_lme)),
                     rel(tcrossprod(Z),tcrossprod(Z_lme))))
  checks$passed <- checks$relative_error<1e-8
  if (any(!checks$passed)) { print(checks); stop("Reference construction checks failed.") }
  if (qr(C)$rank!=ncol(B)) stop("Coefficient transformation lost directions.")

  matrices <- list(y=matrix(dat$y,ncol=1),B=B,S=S,X_ordinary=Xordinary,
                   X0_r=X0,Z_r=Z,C_r=C,X_fixed_r=Xfixed,
                   X_fixed_lme=X_lme,Z_lme=Z_lme)
  registry <- data.frame(file=character(),rows=integer(),columns=integer(),meaning=character())
  meanings <- c(y="Observed responses",B="Constrained spline basis in original coordinates",
    S="Matching scaled spline penalty",X_ordinary="Ordinary fixed predictor: intercept",
    X0_r="Unpenalized spline columns",Z_r="Penalized and rescaled spline columns",
    C_r="Map [r; alpha] back to original spline coefficients",
    X_fixed_r="Intercept plus unpenalized spline components",
    X_fixed_lme="Fixed matrix from the actual captured lme formula/data",
    Z_lme="Random matrix from the actual captured lme random formula/data")
  for (nm in names(matrices)) {
    M <- matrices[[nm]]
    write.table(M,file.path(out,paste0(nm,".csv")),sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)
    registry <- rbind(registry,data.frame(file=paste0(nm,".csv"),rows=nrow(M),columns=ncol(M),meaning=meanings[[nm]]))
  }
  write.csv(registry,file.path(out,"matrix_dictionary.csv"),row.names=FALSE)
  write.csv(dat,file.path(out,"data.csv"),row.names=FALSE)
  write.csv(checks,file.path(out,"reference_checks.csv"),row.names=FALSE)
  metadata <- c(R_version=R.version.string,mgcv_version=as.character(packageVersion("mgcv")),
    nlme_version=as.character(packageVersion("nlme")),n=n,k_requested=k,
    constrained_basis_columns=ncol(B),penalty_rank=r,null_rank=ncol(X0),
    spline_type="bs",degree=3,penalty_derivative_order=2,full_knot_entries=length(knots),method="REML",data_source="shared CSV; no R simulation",
    source_data_case=basename(dirname(data_csv)),
    source_data_md5=unname(tools::md5sum(data_csv)),
    source_knots_md5=unname(tools::md5sum(knots_csv)),
    reference_scope="installed version; historical revision not automatically reproduced")
  write.table(data.frame(key=names(metadata),value=unname(metadata)),file.path(out,"metadata.tsv"),
              sep="\t",row.names=FALSE,quote=FALSE)
  # Preserve exact R objects as well as portable numerical exports.
  saveRDS(input,file.path(out,"lme_inputs.rds"))
  saveRDS(list(smooth=sm,conversion=re,fit=fit,metadata=metadata),file.path(out,"reference_objects.rds"))
  writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

  report <- c("STAGE 1: WHAT DOES MGCV PASS TO LME?", "",
    "This example uses one cubic B-spline, one penalty, an intercept,",
    "Gaussian noise, one dummy spline group and REML. It is not an EEG overlap model.","",
    "Data source: saved Julia-generated CSV; R does not regenerate data.",
    paste("Source data case:",basename(dirname(data_csv))),
    paste("Source data MD5:",unname(tools::md5sum(data_csv))),
    paste("Source knots MD5:",unname(tools::md5sum(knots_csv))),"",
    paste("Requested basis size:",k),paste("Full knot entries:",length(knots)),paste("Constrained basis columns:",ncol(B)),
    paste("Penalized columns:",r),paste("Unpenalized spline columns:",ncol(X0)),"",
    "ACTUAL INPUTS (captured before lme fitting)",
    paste("fixed.formula:",paste(deparse(input$fixed),collapse=" ")),
    paste("random block class:",paste(class(input$random[[1]]),collapse=", ")),
    paste("random formula:",paste(deparse(formula(input$random[[1]])),collapse=" ")),
    paste("dummy group:",group_name,"; distinct levels:",length(unique(group))),
    paste("method:",paste(input$method,collapse=", ")),
    paste("correlation:",if(is.null(input$correlation))"NULL (independent residuals)" else "see RDS"),
    paste("weights:",if(is.null(input$weights))"NULL (one residual variance)" else "see RDS"),
    paste("data fields:",paste(names(input$data),collapse=", ")),
    "control settings:",capture.output(str(input$control)),"",
    "MATRIX DICTIONARY",capture.output(print(registry,row.names=FALSE)),"",
    "REFERENCE CHECKS",capture.output(print(checks,row.names=FALSE)),"",
    "ESTIMATED OUTPUTS (not inputs)",paste("lambda:",format(fit$gam$sp[1],digits=12)),
    paste("noise variance:",format(fit$lme$sigma^2,digits=12)),"",
    "Important: lambda is learned by fitting; k was chosen before fitting.",
    "The final Julia implementation need not depend on R.",
    "This run records the installed versions. It does not claim historical-version compatibility.")
  writeLines(report,file.path(out,"INPUTS_EXPLAINED.txt"))
  plot_file <- file.path(out,"reference_curve.pdf")
  grDevices::pdf(plot_file,width=8,height=5)
  plot(dat$x,dat$y,col="grey60",pch=16,xlab="Predictor x",ylab="Response y",
       main="R reference: one penalized spline")
  if("truth"%in%names(dat))lines(dat$x,dat$truth,col="black",lwd=2,lty=2)
  lines(dat$x,as.numeric(predict(fit$gam)),col="steelblue4",lwd=2)
  if("truth"%in%names(dat))legend("topleft",c("Known simulated relationship","mgcv fit"),col=c("black","steelblue4"),lty=c(2,1),lwd=2,bty="n")
  grDevices::dev.off()
  cat(paste(report,collapse="\n"),"\n\nSaved outputs to:",normalizePath(out),"\n")
}
run()
