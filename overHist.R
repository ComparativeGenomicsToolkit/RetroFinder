overhist <- function(x1, x2, breaks="Sturges",
 	col1="blue", col2="green", xlab0 = paste(deparse(substitute(x1)), "and", deparse(substitute(x2))),
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

 	main0 <- ""

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
fPlot <- function(feature , brk, title, p, neg)
{ 
overhist(feature[p$label==neg], feature[p$label==1], breaks=brk, xlab0 = paste(title))
}
doModel <-function(vec, feature,title,cutoff=0.5)
{
model<-glm(vec$label~feature,family=binomial)
summary(model)
yfit<-exp(model$fitted.values)/(1+exp(model$fitted.values))
pred<-model$fitted.values>cutoff
y<- c(table(vec$label[pred==1]), table(vec$label[pred==0]))
sens<-1-(y[1]/(y[1]+y[2]))
spec<-y[3]/(y[3]+y[4])
result<-c(sens,spec)
#fPlot(yfit,100,title,vec,0)
result
}
