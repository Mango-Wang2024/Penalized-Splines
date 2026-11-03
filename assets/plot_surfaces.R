# Six panels from the saved evaluation grid: three surfaces and three differences.
plot_response_surfaces <- function(grid, r_grid, julia_grid, keys, output, resolution=400, observations=NULL) {
  stopifnot(length(keys)==2L, nrow(grid)==nrow(r_grid), nrow(grid)==nrow(julia_grid))
  x <- sort(unique(grid[[keys[1]]])); y <- sort(unique(grid[[keys[2]]]))
  stopifnot(!anyDuplicated(grid[keys]))
  stopifnot(nrow(grid)==length(x)*length(y))
  for(key in keys)stopifnot(isTRUE(all.equal(grid[[key]],r_grid[[key]],tolerance=1e-12)),
                          isTRUE(all.equal(grid[[key]],julia_grid[[key]],tolerance=1e-12)))
  # Reorder the prediction rows together into the Cartesian plotting grid.
  order <- order(grid[[keys[2]]],grid[[keys[1]]])
  grid <- grid[order,];r_grid <- r_grid[order,];julia_grid <- julia_grid[order,]
  labels <- lapply(keys,function(key) {
    if(grepl("^x[0-9]+$",key))bquote(x[.(as.integer(sub("^x","",key)))]) else key
  })
  if(!is.null(observations)) {
    stopifnot(nrow(observations)>0L,all(c(keys,"y")%in%names(observations)),
              all(is.finite(as.matrix(observations[,c(keys,"y")]))))
    stopifnot(all(observations[[keys[1]]]>=min(x)&observations[[keys[1]]]<=max(x)),
              all(observations[[keys[2]]]>=min(y)&observations[[keys[2]]]<=max(y)))
  }
  # All three surfaces and all observed responses share one vertical range.
  z_limits <- extendrange(range(grid$truth,r_grid$fitted,julia_grid$fitted,observations$y),f=.025)
  point_colour <- adjustcolor("#35414D",alpha.f=.58)
  difference <- julia_grid$fitted-r_grid$fitted
  largest <- max(abs(difference))
  exponent <- if(largest>0)floor(log10(largest)) else 0L
  magnitude <- 10^exponent
  bound <- max(1,ceiling(largest/magnitude))
  truth_errors <- list(r_grid$fitted-grid$truth,julia_grid$fitted-grid$truth)
  truth_largest <- max(abs(unlist(truth_errors)))
  truth_magnitude <- if(truth_largest>0)10^floor(log10(truth_largest)) else 1
  truth_bound <- max(1,ceiling(truth_largest/truth_magnitude))*truth_magnitude
  truth_stats <- data.frame(method=c("R","Julia"),grid_points=nrow(grid),
    rmse=vapply(truth_errors,function(e)sqrt(mean(e^2)),0.),
    maximum_absolute_error=vapply(truth_errors,function(e)max(abs(e)),0.))
  colours <- colorRampPalette(c("#326795","#FFFFFF","#A65048"))(201)
  png(output,width=11.0,height=7.45,units="in",res=resolution,pointsize=11,type="quartz",bg="white")
  on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,3),oma=c(0,0,if(is.null(observations))0 else 1.6,0),family="sans",cex=1,col.axis="#333333",col.lab="#222222",
      fg="#555555",font.main=1,cex.main=1.0,lwd=.6)
  surface <- function(values,title,colour) {
    par(mar=c(1.5,1.2,2.5,1.2),mgp=c(1.6,.4,0),cex.axis=.75,cex.lab=.9)
    projection <- persp(x,y,matrix(values,length(x),length(y)),theta=-40,phi=25,expand=.65,
      col=adjustcolor(colour,alpha.f=if(is.null(observations))1 else .42),
      border=if(is.null(observations))"#FFFFFFA0" else "#FFFFFF70",shade=.2,lwd=.3,
      xlab="",ylab="",zlab="Response",
      zlim=z_limits,ticktype="detailed",nticks=3)
    if(!is.null(observations)) {
      positions <- trans3d(observations[[keys[1]]],observations[[keys[2]]],observations$y,projection)
      points(positions,pch=16,cex=.42,col=point_colour)
    }
    # Draw plotmath axis labels separately; persp treats expressions as literal text.
    usr <- par("usr");dx <- diff(usr[1:2]);dy <- diff(usr[3:4])
    px <- trans3d(mean(x),min(y),z_limits[1],projection)
    py <- trans3d(min(x),mean(y),z_limits[1],projection)
    text(px$x+.08*dx,px$y-.11*dy,labels=labels[[1]],srt=35,xpd=NA,cex=.9)
    text(py$x-.11*dx,py$y-.07*dy,labels=labels[[2]],srt=-45,xpd=NA,cex=.9)
    title(main=title,line=.8,adj=0,font.main=1,cex.main=1.05)
  }
  surface(grid$truth,"(a) Known surface","#A9B7C4")
  surface(r_grid$fitted,"(b) R: mgcv / nlme","#709BB5")
  surface(julia_grid$fitted,"(c) Julia: MixedModels","#D09B78")
  heatmap <- function(values,title,limit,power=0L,legend_label=expression(Delta*hat(y))) {
    par(mar=c(4.1,3.9,2.5,5.0),mgp=c(2,.4,0),tcl=-.2,cex.axis=.82,cex.lab=.9,xaxs="i",yaxs="i")
    image(x,y,matrix(values,length(x),length(y)),col=colours,zlim=c(-limit,limit),
      axes=FALSE,xlab=as.expression(labels[[1]]),ylab=as.expression(labels[[2]]),useRaster=TRUE,
      xlim=range(x),ylim=range(y))
    axis(1,at=pretty(range(x),n=3));axis(2,at=pretty(range(y),n=3),las=1)
    box(col="#555555",lwd=.6)
    title(main=title,line=.8,adj=0,font.main=1,cex.main=1.05)
    span <- diff(range(x));yrange <- range(y);height <- diff(yrange)
    left <- max(x)+.08*span;right <- left+.04*span
    par(xpd=NA)
    rasterImage(as.raster(matrix(rev(colours),ncol=1)),left,min(y),right,max(y),interpolate=FALSE)
    rect(left,min(y),right,max(y),border="#777777",lwd=.4)
    ticks <- pretty(c(-limit,limit),n=5);ticks <- ticks[abs(ticks)<=limit]
    positions <- min(y)+(ticks+limit)/(2*limit)*height
    segments(right,positions,right+.018*span,positions,col="#555555",lwd=.5)
    text(right+.035*span,positions,format(signif(ticks,3),trim=TRUE),adj=c(0,.5),cex=.76)
    if(power!=0L)text((left+right)/2,max(y)+.075*height,bquote("×"~10^.(power)),cex=.78)
    text(right+.25*span,mean(yrange),legend_label,srt=90,cex=.82)
    par(xpd=FALSE)
  }
  heatmap(difference/magnitude,"(d) Julia − R",bound,exponent)
  heatmap(truth_errors[[1]],"(e) R − known surface",truth_bound,legend_label="Prediction − truth")
  heatmap(truth_errors[[2]],"(f) Julia − known surface",truth_bound,legend_label="Prediction − truth")
  if(!is.null(observations)) {
    # One legend for the same observed data shown in all three surface panels.
    par(fig=c(0,1,0,1),mar=rep(0,4),oma=rep(0,4),new=TRUE)
    plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
    legend(.5,.995,legend=sprintf("Training observations (n = %d)",nrow(observations)),
      pch=16,col=point_colour,pt.cex=.65,bty="n",cex=.82,xjust=.5,yjust=1)
  }
  invisible(list(maximum_absolute_difference=largest,colour_exponent=exponent,
                 truth_comparison=truth_stats,truth_colour_limit=truth_bound))
}
