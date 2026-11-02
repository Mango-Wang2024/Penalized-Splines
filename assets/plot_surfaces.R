# Draw the saved two-predictor evaluation grid without refitting either model.
plot_response_surfaces <- function(grid, r_grid, julia_grid, keys, output, resolution=400) {
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
  z_limits <- range(grid$truth,r_grid$fitted,julia_grid$fitted)
  difference <- julia_grid$fitted-r_grid$fitted
  largest <- max(abs(difference))
  exponent <- if(largest>0)floor(log10(largest)) else 0L
  magnitude <- 10^exponent
  bound <- max(1,ceiling(largest/magnitude))
  colours <- colorRampPalette(c("#326795","#FFFFFF","#A65048"))(201)
  png(output,width=7.6,height=6.4,units="in",res=resolution,pointsize=9.5,type="quartz",bg="white")
  on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,2),family="sans",cex=1,col.axis="#333333",col.lab="#222222",
      fg="#555555",font.main=1,cex.main=1.0,lwd=.6)
  surface <- function(values,title,colour) {
    par(mar=c(1.5,1.2,2.5,1.2),mgp=c(1.6,.4,0),cex.axis=.75,cex.lab=.9)
    projection <- persp(x,y,matrix(values,length(x),length(y)),theta=-40,phi=25,expand=.65,
      col=colour,border="#FFFFFFA0",shade=.2,lwd=.45,
      xlab="",ylab="",zlab="Response",
      zlim=z_limits,ticktype="detailed",nticks=3)
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
  par(mar=c(4.2,4,2.5,5.2),mgp=c(2,.4,0),tcl=-.2,cex.axis=.85,cex.lab=.95,xaxs="i",yaxs="i")
  image(x,y,matrix(difference/magnitude,length(x),length(y)),col=colours,zlim=c(-bound,bound),
    axes=FALSE,xlab=as.expression(labels[[1]]),ylab=as.expression(labels[[2]]),useRaster=TRUE,
    xlim=range(x),ylim=range(y))
  axis(1,at=pretty(range(x),n=3));axis(2,at=pretty(range(y),n=3),las=1)
  box(col="#555555",lwd=.6)
  title(main="(d) Difference: Julia − R",line=.8,adj=0,font.main=1,cex.main=1.05)
  # A separate, centered colour scale replaces contour labels and in-panel text.
  span <- diff(range(x));yrange <- range(y);height <- diff(yrange)
  left <- max(x)+.08*span;right <- left+.035*span
  par(xpd=NA)
  rasterImage(as.raster(matrix(rev(colours),ncol=1)),left,min(y),right,max(y),interpolate=FALSE)
  rect(left,min(y),right,max(y),border="#777777",lwd=.4)
  ticks <- pretty(c(-bound,bound),n=4);ticks <- ticks[abs(ticks)<=bound]
  positions <- min(y)+(ticks+bound)/(2*bound)*height
  segments(right,positions,right+.018*span,positions,col="#555555",lwd=.5)
  text(right+.035*span,positions,format(ticks,trim=TRUE),adj=c(0,.5),cex=.78)
  text((left+right)/2,max(y)+.075*height,bquote(times~10^.(exponent)),cex=.8)
  text(right+.14*span,mean(yrange),expression(Delta*hat(y)),srt=90,cex=.9)
  par(xpd=FALSE)
  invisible(list(maximum_absolute_difference=largest,colour_exponent=exponent))
}
