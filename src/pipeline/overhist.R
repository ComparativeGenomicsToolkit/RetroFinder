 overhist <- function(x1, x2, breaks="Sturges",
 	col1="red", col2="green",
 	coli=rgb( (col2rgb(col1) + col2rgb(col2))[1]/2,
  		  (col2rgb(col1) + col2rgb(col2))[2]/2,
     		  (col2rgb(col1) + col2rgb(col2))[3]/2,
 		  max=255
 		),
 	...) {

 	h12 <- hist(c(x1,x2), breaks, freq=FALSE, plot=FALSE)

 	h1 <- hist(x1, breaks=h12$breaks, freq=FALSE, plot=FALSE)
 	h2 <- hist(x2, breaks=h12$breaks, freq=FALSE, plot=FALSE)

 	h0 <- h1
 	h0$intensities <- pmin(h1$intensities, h2$intensities)
 	h0$density     <- pmin(h1$density,     h2$density)
 	h0$counts      <- pmin(h1$counts,      h2$counts)

 	xlim <- range(c(x1,x2))
 	ylim <- range(c(h1$intensities, h2$intensities))

 	xlab <- list(...)$xlab
 	main <- list(...)$main

 	main0 <- paste("Histogram of", deparse(substitute(x1)), "and", deparse(substitute(x2)))
 	xlab0 <- paste(deparse(substitute(x1)), "and", deparse(substitute(x2)))

 	if(is.null(xlab) && is.null(main)) {
 		plot(h1, freq=FALSE, xlim=xlim, ylim=ylim, col=col1, border=col1, xlab=xlab0, main=main0, ...)
 		plot(h2, freq=FALSE, xlim=xlim, ylim=ylim, col=col2, border=col2, add=TRUE, ...)
 		plot(h0, freq=FALSE, xlim=xlim, ylim=ylim, col=coli, border=coli, add=TRUE, ...)
 	} else if(is.null(xlab) && !is.null(main)) {
 		plot(h1, freq=FALSE, xlim=xlim, ylim=ylim, col=col1, border=col1, xlab=xlab0, ...)
 		plot(h2, freq=FALSE, xlim=xlim, ylim=ylim, col=col2, border=col2, add=TRUE, ...)
 		plot(h0, freq=FALSE, xlim=xlim, ylim=ylim, col=coli, border=coli, add=TRUE, ...)
 	} else if(!is.null(xlab) && is.null(main)) {
 		plot(h1, freq=FALSE, xlim=xlim, ylim=ylim, col=col1, border=col1, main=main0, ...)
 		plot(h2, freq=FALSE, xlim=xlim, ylim=ylim, col=col2, border=col2, add=TRUE, ...)
 		plot(h0, freq=FALSE, xlim=xlim, ylim=ylim, col=coli, border=coli, add=TRUE, ...)
 	} else if(!is.null(xlab) && !is.null(main)) {
 		plot(h1, freq=FALSE, xlim=xlim, ylim=ylim, col=col1, border=col1, ...)
 		plot(h2, freq=FALSE, xlim=xlim, ylim=ylim, col=col2, border=col2, add=TRUE, ...)
 		plot(h0, freq=FALSE, xlim=xlim, ylim=ylim, col=coli, border=coli, add=TRUE, ...)
 	}
}
