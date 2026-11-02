#!/usr/bin/env Rscript
run <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  stopifnot(length(args)==4L,file.exists(args[4]))
  inputs <- normalizePath(args[1]);out <- normalizePath(args[2]);report <- args[3]
  script <- gsub("~+~"," ",sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)),fixed=TRUE)
  project <- dirname(dirname(dirname(normalizePath(script))))
  if(file.exists(report))stop("The report already exists; choose a new report path")
  rdir <- file.path(out,"r_reference");jdir <- file.path(out,"julia")
  r <- read.delim(file.path(rdir,"fit_summary.tsv"));j <- read.delim(file.path(jdir,"fit_summary.tsv"))
  re <- read.delim(file.path(rdir,"component_estimates.tsv"));je <- read.delim(file.path(jdir,"component_estimates.tsv"))
  cases <- read.delim(file.path(inputs,"cases.tsv"))$case
  stopifnot(identical(r$case,cases),identical(j$case,cases))
  hashes <- read.delim(file.path(rdir,"source_checksums.tsv"))
  for(i in seq_len(nrow(hashes)))stopifnot(hashes$md5[i]==unname(tools::md5sum(file.path(inputs,hashes$case[i],hashes$file[i]))))
  scaled <- function(a,b)sqrt(sum((a-b)^2))/max(1,sqrt(sum(a*a)),sqrt(sum(b*b)))
  relative <- function(a,b) {
    if(is.infinite(a)&&is.infinite(b)&&sign(a)==sign(b))return(0)
    if(!is.finite(a)||!is.finite(b))return(Inf)
    abs(a-b)/max(abs(a),abs(b),.Machine$double.eps)
  }
  checks <- list();estimates <- list();ci <- 0L;ei <- 0L
  check <- function(case,component,label,error,tolerance) {
    ci <<- ci+1L;checks[[ci]] <<- data.frame(case,component,check=label,error,tolerance,passed=is.finite(error)&&error<=tolerance)
  }
  for(i in seq_along(cases)) {
    name <- cases[i];src <- file.path(inputs,name);terms <- read.delim(file.path(src,"terms.tsv"));keys <- terms$predictor
    raw <- read.csv(file.path(src,"data.csv"));rf <- read.csv(file.path(rdir,name,"fitted.csv"));jf <- read.csv(file.path(jdir,name,"fitted.csv"))
    rc <- read.csv(file.path(rdir,name,"components.csv"));jc <- read.csv(file.path(jdir,name,"components.csv"))
    rg <- read.csv(file.path(rdir,name,"grid.csv"));jg <- read.csv(file.path(jdir,name,"grid.csv"))
    check(name,"all","same observations",max(scaled(as.matrix(rf[,c(keys,"y")]),as.matrix(raw[,c(keys,"y")])),scaled(as.matrix(jf[,c(keys,"y")]),as.matrix(raw[,c(keys,"y")]))),1e-12)
    check(name,"all","R REML converged",as.numeric(!as.logical(r$converged[i])||r$method[i]!="REML"),0)
    check(name,"all","Julia REML converged",as.numeric(!as.logical(j$converged[i])||j$method[i]!="REML"),0)
    check(name,"all","one variance per component",as.numeric(r$components[i]!=length(keys)||j$variance_parameters[i]!=length(keys)||r$random_directions[i]!=sum(terms$k-2L)||j$random_directions[i]!=sum(terms$k-2L)),0)
    check(name,"all","residual variance",relative(r$sigma2[i],j$sigma2[i]),1e-5)
    check(name,"all","fitted responses",scaled(rf$fitted,jf$fitted),1e-6)
    check(name,"all","same grid predictors",scaled(as.matrix(rg[,keys,drop=FALSE]),as.matrix(jg[,keys,drop=FALSE])),1e-12)
    check(name,"all","grid predictions",scaled(rg$fitted,jg$fitted),1e-6)
    for(key in keys) {
      a <- re[re$case==name&re$component==key,];b <- je[je$case==name&je$component==key,];stopifnot(nrow(a)==1L,nrow(b)==1L)
      rk <- scan(file.path(rdir,name,key,"knots.csv"),quiet=TRUE);jk <- scan(file.path(jdir,name,key,"knots.csv"),quiet=TRUE)
      check(name,key,"same knots",scaled(rk,jk),1e-12)
      check(name,key,"smoothing lambda",relative(a$lambda,b$lambda),1e-5)
      check(name,key,"spline variance",relative(a$r2,b$r2),1e-5)
      check(name,key,"component curve",scaled(rc[[key]],jc[[key]]),1e-6)
      check(name,key,"variance-to-smoothing identity",max(relative(a$lambda,r$sigma2[i]/a$r2),relative(b$lambda,j$sigma2[i]/b$r2)),1e-10)
      ei <- ei+1L;estimates[[ei]] <- data.frame(case=name,component=key,R_lambda=a$lambda,Julia_lambda=b$lambda,
        lambda_relative_error=relative(a$lambda,b$lambda),R_r2=a$r2,Julia_r2=b$r2,r2_relative_error=relative(a$r2,b$r2),
        R_sigma2=r$sigma2[i],Julia_sigma2=j$sigma2[i],sigma2_relative_error=relative(r$sigma2[i],j$sigma2[i]))
    }
  }
  checks <- do.call(rbind,checks);estimates <- do.call(rbind,estimates)
  dst <- file.path(out,"comparison");dir.create(dst,showWarnings=FALSE)
  for(pair in list(list("checks.tsv",checks),list("estimates.tsv",estimates)))write.table(pair[[2]],file.path(dst,pair[[1]]),sep="\t",quote=FALSE,row.names=FALSE)
  # Show the joint surface for two predictors; component plots also work for any count.
  main <- cases[1];raw <- read.csv(file.path(inputs,main,"data.csv"));keys <- read.delim(file.path(inputs,main,"terms.tsv"))$predictor
  rg <- read.csv(file.path(rdir,main,"grid.csv"));jg <- read.csv(file.path(jdir,main,"grid.csv"));grid <- read.csv(file.path(inputs,main,"grid.csv"))
  if(length(keys)==2L) {
    source(file.path(project,"assets","plot_surfaces.R"))
    surface_stats <- plot_response_surfaces(grid,rg,jg,keys,file.path(dst,"fitted_surfaces.png"))
  }
  rc <- read.csv(file.path(rdir,main,"components.csv"));jc <- read.csv(file.path(jdir,main,"components.csv"))
  png(file.path(dst,"component_curves.png"),width=1800,height=700*ceiling(length(keys)/2),res=200,type="quartz")
  par(mfrow=c(ceiling(length(keys)/2),min(2,length(keys))),mar=c(4,4,2.3,1),cex=.9)
  for(key in keys) {
    ix <- order(raw[[key]]);truth <- raw[[paste0("truth_",key)]];truth <- truth-mean(truth)
    plot(raw[[key]][ix],truth[ix],type="l",lty=3,lwd=2,col="#666666",xlab=key,ylab="Centered component",main=paste("Spline component:",key),ylim=range(truth,rc[[key]],jc[[key]]))
    lines(raw[[key]][ix],rc[[key]][ix],lwd=2,col="#1675AE");lines(raw[[key]][ix],jc[[key]][ix],lwd=2,lty=2,col="#D75E1D")
    legend("topleft",c("Known component","R","Julia"),col=c("#666666","#1675AE","#D75E1D"),lty=c(3,1,2),lwd=2,bty="n",cex=.8)
  }
  dev.off()
  relative_path <- function(path) {
    a <- strsplit(normalizePath(dirname(report)),"/",fixed=TRUE)[[1]];b <- strsplit(normalizePath(path),"/",fixed=TRUE)[[1]];n <- 0L
    while(n<min(length(a),length(b))&&a[n+1L]==b[n+1L])n <- n+1L
    paste(c(rep("..",length(a)-n),b[seq_along(b)>n]),collapse="/")
  }
  link <- function(label,path)paste0("[",label,"](",relative_path(path),")")
  image_tag <- function(file,alt)paste0('<a href="',relative_path(file),'"><img src="',relative_path(file),'" alt="',alt,'" width="800"></a>')
  lines <- c("# Stage 2 results","",
    "R and Julia independently prepare and fit the same additive model. R uses `mgcv::gamm` / `nlme::lme`; Julia uses `MixedModels.jl`. REML jointly estimates a separate spline variance for each predictor and one common noise variance.","",
    sprintf("**%d/%d comparisons passed across %d cases.** Spline-component counts: %s.",sum(checks$passed),nrow(checks),length(cases),paste(unique(r$components),collapse=", ")),"",
    "## Smoothing estimates","",
    "Each component has its own $\\lambda_j=\\sigma^2/r_j^2$. The table compares these smoothing strengths; the variance values and their errors are retained in the detailed estimates.","",
    "| Case | Component | R λ<sub>j</sub> | Julia λ<sub>j</sub> | Relative error |","| --- | --- | ---: | ---: | ---: |")
  for(i in seq_len(nrow(estimates)))lines <- c(lines,sprintf("| `%s` | $%s$ | %.8g | %.8g | %.2e |",estimates$case[i],gsub("x","x_",estimates$component[i]),estimates$R_lambda[i],estimates$Julia_lambda[i],estimates$lambda_relative_error[i]))
  lines <- c(lines,"",sprintf("Largest scaled difference in training fitted responses: **%.2e**; on the shared evaluation grid: **%.2e**.",max(checks$error[checks$check=="fitted responses"]),max(checks$error[checks$check=="grid predictions"])),"",
    "## Joint response and individual components","")
  if(length(keys)==2L)lines <- c(lines,image_tag(file.path(dst,"fitted_surfaces.png"),"Known, R and Julia response surfaces and their difference"),"",
    sprintf("Panels (a–c) share the same viewing angle and response limits. Panel (d) shows Julia minus R with a symmetric colour scale centred at zero; the multiplier above its colourbar gives the scale. Its maximum absolute difference is **%.2e**. The known surface is a simulation reference; only the noisy observations are used for fitting.",surface_stats$maximum_absolute_difference),"")
  lines <- c(lines,image_tag(file.path(dst,"component_curves.png"),"Centered additive spline components"),"",
    "Each panel isolates one spline contribution. Known components are centered over the training observations to match the fitted smooths; their removed average levels belong to the overall intercept.","",
    "## Supporting checks","",
    paste0(link("Individual comparisons",file.path(dst,"checks.tsv"))," and ",link("detailed estimates",file.path(dst,"estimates.tsv"))," record all cases. The acceptance criteria are defined in ",link("VALIDATION.md",file.path(project,"VALIDATION.md")),"."),"",
    paste0(link("Julia unit tests",args[4])," independently check covariance, the REML objective, coefficients and restoration of observation coordinates for 1, 2, 3 and 10 components, including unequal block sizes and the zero-variance boundary."),"",
    paste0(link("R fitting status",file.path(rdir,"fit_summary.tsv")),"; ",link("Julia fitting status",file.path(jdir,"fit_summary.tsv")),"; ",link("source fingerprints",file.path(rdir,"source_checksums.tsv")),"."),"",
    "Both languages use the same saved observations, per-predictor settings and evaluation grid. Package versions are recorded in each branch's `versions.tsv`; Julia dependencies are pinned by the root `Manifest.toml`.")
  writeLines(lines,report)
  cat(sprintf("%d/%d Stage 2 comparisons passed\n",sum(checks$passed),nrow(checks)))
  if(!all(checks$passed))stop("Stage 2 comparisons failed; see checks.tsv")
  if(normalizePath(dirname(report))==normalizePath(project)) {
    source(file.path(project,"assets","plot_matrices.R"))
    plot_fitting_matrices(project,out)
  }
}
run()
