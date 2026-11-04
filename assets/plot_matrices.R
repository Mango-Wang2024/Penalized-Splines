# Numerical heatmaps from the saved main-case Julia matrices.
plot_fitting_matrices <- function(project, fitting_outputs) {
  out <- file.path(project,"assets");case <- "uniform_k8"
  keys <- read.delim(file.path(project,"data","inputs",case,"terms.tsv"))$predictor
  read_matrix <- function(path)as.matrix(read.table(path,sep=","))
  prep <- file.path(project,"model_preparation","outputs","shared","julia",case)
  fit <- file.path(fitting_outputs,"julia",case)
  parts <- lapply(keys,function(key)read_matrix(file.path(prep,key,"B_raw.csv")))
  raw <- do.call(cbind,parts)
  centered <- do.call(cbind,lapply(keys,function(key)read_matrix(file.path(prep,key,"B.csv"))))
  Zf <- read_matrix(file.path(fit,"Zf.csv"))
  Xf <- read_matrix(file.path(fit,"Xf.csv"))
  loading <- read_matrix(file.path(fit,"relative_loading.csv"))
  stopifnot(all(is.finite(c(raw,centered,Zf,Xf,loading))))
  widths <- list(vapply(parts,ncol,0L),vapply(parts,ncol,0L)-1L,vapply(parts,ncol,0L)-2L)
  diverging <- colorRampPalette(c("#2166AC","#FFFFFF","#B2182B"))(255)
  positive <- colorRampPalette(c("#FFFFFF","#E6A4A2","#B2182B"))(255)
  panel <- function(M,title,tag="",positive_only=FALSE,boundaries=numeric(),
                    row_label="Observation",column_label="Basis column",labels=NULL) {
    nr <- nrow(M);nc <- ncol(M);limit <- max(abs(M))
    limits <- if(positive_only)c(0,limit) else c(-limit,limit)
    palette <- if(positive_only)positive else diverging
    par(mar=c(3.5,3.8,2.1,4.8),mgp=c(2.1,.5,0),tcl=-.2,cex.axis=.82,cex.lab=.92,lwd=.6,xaxs="i",yaxs="i")
    image(seq(.5,nc+.5,1),seq(.5,nr+.5,1),t(M[nr:1,,drop=FALSE]),
      zlim=limits,col=palette,axes=FALSE,xlab="",ylab="",useRaster=TRUE)
    ticks <- if(nc<=16)seq_len(nc) else unique(round(seq(1,nc,length.out=8)))
    axis(1,at=ticks,labels=if(is.null(labels))ticks else labels[ticks])
    rows <- if(nr<=20)unique(c(1,seq(4,nr,4),nr)) else unique(c(1,seq(50,nr,50),nr))
    axis(2,at=nr+1-rows,labels=rows,las=1)
    if(length(boundaries))abline(v=boundaries+.5,col="#555555",lty=2,lwd=.65)
    box(col="#555555");title(main=title,font.main=1,cex.main=1.1,line=.6)
    if(nzchar(tag))mtext(tag,3,.6,adj=0,font=2,cex=.9)
    mtext(column_label,1,2.15,cex=.87);mtext(row_label,2,2.75,cex=.87)
    left <- nc+.5+.075*nc;right <- left+.04*nc
    par(xpd=NA)
    rasterImage(as.raster(matrix(rev(palette),ncol=1)),left,.5,right,nr+.5,interpolate=FALSE)
    ct <- if(positive_only)c(0,limit/2,limit) else c(-limit,-limit/2,0,limit/2,limit)
    positions <- .5+(ct-limits[1])/diff(limits)*nr
    text(right+.025*nc,positions,format(signif(ct,3),trim=TRUE),adj=c(0,.5),cex=.75)
    text((left+right)/2,nr+.5+.045*nr,"Value",cex=.75)
    par(xpd=FALSE)
  }
  draw <- function(name,code,width=5.5,height=4.55) {
    png(file.path(out,name),width=width,height=height,units="in",res=400,pointsize=9.5,type="quartz",bg="white")
    force(code);dev.off()
  }
  boundary <- function(w)head(cumsum(w),-1L)
  fixed_labels <- c("X",paste0("X0,",seq_along(keys)))
  draw("fitting_matrices.png",{
    par(mfrow=c(2,2),oma=c(.15,.15,.1,.3))
    panel(raw,expression(B[plain(raw)]),"(a)",TRUE,boundary(widths[[1]]))
    panel(centered,expression(B),"(b)",boundaries=boundary(widths[[2]]))
    panel(Zf,expression(Z[f]),"(c)",boundaries=boundary(widths[[3]]))
    panel(Xf,expression(X[f]),"(d)",labels=fixed_labels,column_label="Fixed model part")
  },8.1,6.5)
  draw("b_raw_heatmap.png",panel(raw,expression(B[plain(raw)]),positive_only=TRUE,boundaries=boundary(widths[[1]])))
  draw("b_heatmap.png",panel(centered,expression(B),boundaries=boundary(widths[[2]])))
  draw("zf_heatmap.png",panel(Zf,expression(Z[f]),boundaries=boundary(widths[[3]]),column_label="Penalized basis column"))
  draw("xf_heatmap.png",panel(Xf,expression(X[f]),labels=fixed_labels,column_label="Fixed model part"))
  draw("loading_heatmap.png",panel(loading,expression(Lambda),positive_only=TRUE,boundaries=boundary(widths[[3]]),
    row_label="Penalized coefficient row",column_label="Penalized coefficient column"))
}
