# Plot the recorded benchmark; this script does not fit or change any models.
# Arguments: project directory, figure output directory.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
project <- normalizePath(args[1])
output <- args[2]
results <- read.delim(file.path(project,"spline_fitting/outputs/shared/comparison/qr_benchmark_summary.tsv"))
labels <- c(uniform_k8="Uniform knots",default_k8="Default knots",uneven_k10="Uneven knots",
            duplicates_k6="Repeated values",small_k5="Small basis")
cases <- unname(labels[results$case])
rows <- rev(seq_len(nrow(results)))
colours <- c(direct="#007A83",QR="#B55D32")
symbols <- c(direct=16L,QR=17L)
dir.create(output,recursive=TRUE,showWarnings=FALSE)

draw <- function(filename,title,x_label,direct,qr,direct_interval=NULL,qr_interval=NULL) {
  png(file.path(output,filename),width=7.6,height=4.15,units="in",res=300,
      pointsize=11,type="quartz",bg="white")
  on.exit(dev.off())
  par(mar=c(4.3,8.2,4.1,1.0),mgp=c(2.5,.55,0),tcl=-.22,
      family="sans",fg="#333333",col.axis="#333333",col.lab="#222222",lwd=.7)
  maximum <- max(direct,qr,direct_interval,qr_interval)
  ticks <- pretty(c(0,maximum*1.08),n=5)
  limits <- c(0,max(ticks))
  plot.new()
  plot.window(xlim=limits,ylim=c(.5,length(rows)+.5),xaxs="i",yaxs="i")
  abline(v=ticks,col="#E5E8EB",lwd=.65)
  abline(h=rows,col="#F0F1F3",lwd=.6)
  axis(1,at=ticks,labels=format(ticks,trim=TRUE),col="#777777",col.ticks="#777777",cex.axis=.9)
  axis(2,at=rows,labels=cases,las=1,tick=FALSE,cex.axis=.92)
  mtext(x_label,side=1,line=2.5,cex=.95)
  mtext(title,side=3,line=2.65,adj=0,font=1,cex=1.12)
  legend(limits[1],length(rows)+1.05,legend=c("Direct fitting","With QR rotation"),
         col=colours,pch=symbols,horiz=TRUE,bty="n",xpd=NA,
         xjust=0,yjust=.5,cex=.9,pt.cex=1.05,x.intersp=.65)
  # Offset the methods vertically so overlapping intervals remain visible.
  direct_y <- rows+.13
  qr_y <- rows-.13
  if(!is.null(direct_interval)) {
    arrows(direct_interval[,1],direct_y,direct_interval[,2],direct_y,
           angle=90,code=3,length=.035,col=colours["direct"],lwd=1.1)
    arrows(qr_interval[,1],qr_y,qr_interval[,2],qr_y,
           angle=90,code=3,length=.035,col=colours["QR"],lwd=1.1)
  }
  points(direct,direct_y,col=colours["direct"],pch=symbols["direct"],cex=1.1)
  points(qr,qr_y,col=colours["QR"],pch=symbols["QR"],cex=1.15)
}

draw("qr_runtime_comparison.png","Runtime comparison","Median time per fit (ms)",
     results$direct_median_ms,results$qr_median_ms,
     cbind(results$direct_q25_ms,results$direct_q75_ms),
     cbind(results$qr_q25_ms,results$qr_q75_ms))

draw("qr_allocation_comparison.png","Memory allocation comparison","Median allocated memory per fit (KiB)",
     results$direct_median_allocated_bytes/1024,results$qr_median_allocated_bytes/1024)
