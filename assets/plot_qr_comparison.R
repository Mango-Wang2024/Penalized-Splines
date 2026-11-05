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

panel <- function(title,x_label,direct,qr,limits,ticks,show_cases,
                  direct_interval=NULL,qr_interval=NULL) {
  par(mar=c(4.1,if(show_cases)7.9 else 1.4,2.2,1.0))
  plot.new()
  plot.window(xlim=limits,ylim=c(.5,length(rows)+.5),xaxs="i",yaxs="i")
  abline(v=ticks,col="#E5E8EB",lwd=.65)
  abline(h=rows,col="#F0F1F3",lwd=.6)
  axis(1,at=ticks,labels=format(ticks,trim=TRUE),col="#777777",col.ticks="#777777",cex.axis=.9)
  if(show_cases)axis(2,at=rows,labels=cases,las=1,tick=FALSE,cex.axis=.92)
  mtext(x_label,side=1,line=2.5,cex=.95)
  mtext(title,side=3,line=.75,adj=0,font=1,cex=1.08)
  # Keep both methods on matching rows in the two panels.
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

png(file.path(output,"qr_comparison.png"),width=10,height=4.2,units="in",res=300,
    pointsize=11,type="quartz",bg="white")
layout(matrix(1:2,nrow=1),widths=c(1.24,1))
par(oma=c(2.0,0,.1,0),mgp=c(2.5,.55,0),tcl=-.22,
    family="sans",fg="#333333",col.axis="#333333",col.lab="#222222",lwd=.7)
panel("(a) Runtime","Median time per fit (ms)",
      results$direct_median_ms,results$qr_median_ms,
      c(.4,.95),seq(.4,.9,.1),TRUE,
      cbind(results$direct_q25_ms,results$direct_q75_ms),
      cbind(results$qr_q25_ms,results$qr_q75_ms))
panel("(b) Memory allocation","Median allocation per fit (KiB)",
      results$direct_median_allocated_bytes/1024,results$qr_median_allocated_bytes/1024,
      c(300,700),seq(300,700,100),FALSE)

# One shared legend sits below both panels.
par(fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE)
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
legend(.5,.025,legend=c("Direct fitting","With QR rotation"),
       col=colours,pch=symbols,horiz=TRUE,bty="n",xjust=.5,yjust=0,
       cex=.9,pt.cex=1.05,x.intersp=.65)
invisible(dev.off())

# Numerical agreement uses the project's existing response and smoothing tolerances.
agreement <- read.delim(file.path(project,"spline_fitting/outputs/shared/comparison/qr_benchmark_agreement.tsv"))
stopifnot(identical(agreement$case,results$case))
ratios <- cbind(agreement$qr_direct_training_scaled_difference/1e-6,
                agreement$qr_direct_grid_scaled_difference/1e-6,
                agreement$qr_direct_lambda_relative_difference/1e-5)
metric_colours <- c("#007A83","#536EA5","#8B5A83")
metric_symbols <- c(16L,17L,15L)

png(file.path(output,"qr_agreement.png"),width=6.7,height=4.6,units="in",res=300,
    pointsize=11,type="quartz",bg="white")
par(mar=c(3.8,7.0,2.0,1.3),oma=c(2.0,0,.1,0),mgp=c(2.5,.55,0),tcl=-.22,
    family="sans",fg="#333333",col.axis="#333333",col.lab="#222222",lwd=.7)
plot.new()
plot.window(xlim=c(1e-5,1.8),ylim=c(.5,length(rows)+.5),log="x",xaxs="i",yaxs="i")
ticks <- 10^(-5:0)
abline(v=ticks,col="#E5E8EB",lwd=.65)
abline(h=rows,col="#F0F1F3",lwd=.6)
abline(v=1,col="#555555",lty=2,lwd=1.0)
axis(1,at=ticks,labels=c("0.00001","0.0001","0.001","0.01","0.1","1"),
     col="#777777",col.ticks="#777777",cex.axis=.9)
axis(2,at=rows,labels=cases,las=1,tick=FALSE,cex.axis=.92)
mtext("Direct–QR difference / numerical tolerance (log scale)",side=1,line=2.5,cex=.95)
mtext("(c) Agreement of fitted results",side=3,line=.75,adj=0,font=1,cex=1.08)
mtext("Tolerance",side=3,line=.4,at=1,cex=.85,col="#555555")
for(j in seq_len(ncol(ratios))) {
  points(ratios[,j],rows+c(.2,0,-.2)[j],col=metric_colours[j],pch=metric_symbols[j],cex=1.1)
}
par(fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE)
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
legend(.5,.025,legend=c("Fitted responses","Grid predictions","Smoothing strengths"),
       col=metric_colours,pch=metric_symbols,horiz=TRUE,bty="n",xjust=.5,yjust=0,
       cex=.9,pt.cex=1.05,x.intersp=.65)
invisible(dev.off())
