plot.ssf=function(coxph,habitat,add=F,type=NULL,ylim=NULL,lty=NULL,col=NULL,lwd=NULL,lty.error=NULL,par=1,int=0,int_sd=0){
####################################################################
########GRAPHS TO COXPH OBJETC - SSF CHANGING IN TIME###############
####################################################################
#ARGUMENTS##########################################################
#coxph = a object coxph returned from the function survival::coxph()##
#habitat = the name of the habitat type to be ploted################
#add = logical (T or F) to plot on previous graph###################
#ylim = vector indicating the min and max of the y scale############
####################################################################
require(survival)
modcoef=coef(coxph)
which.coef=names(modcoef)[grep(habitat,names(modcoef))]
modcoef=coef(coxph)[which.coef]
sevals=(summary(coxph)[["coefficients"]])[which.coef,3]
n=10000
nplot=200
cycle.len=24
betalist=list(b1=modcoef[1],b1.se=sevals[1],
              c1=modcoef[2],c1.se=sevals[2],
              c2=modcoef[3],c2.se=sevals[3],
              s1=modcoef[4],s1.se=sevals[4],
              s2=modcoef[5],s2.se=sevals[5])

#####UNPACK
unpack.list <- function(object) {
  for(.x in names(object)){
    assign(value = object[[.x]], x=.x, envir = parent.frame())
  }
}

#####CURVE
curvefun.sim<-function(t,betalist,se.mod=1.96,cycle.len=24,n=1){
  ##set defaul values to zero
  b1=b1.se=b1.2=b1.2.se=d1=c1=c1.se=c2=c2.se=s1=s1.se=s2=s2.se=0
  ##extract coefficients to plot
  unpack.list(betalist)
  #calculate pointwise CI
  y= rnorm(n,int,int_sd)+par*rnorm(n,b1,b1.se)+
    par*rnorm(n,b1.2,b1.2.se)*d1+
    par*rnorm(n,c1,c1.se)*cos(t*2*pi/cycle.len)+
    par*rnorm(n,c2,c2.se)*cos(t*4*pi/cycle.len)+
    par*rnorm(n,s1,s1.se)*sin(t*2*pi/cycle.len)+
    par*rnorm(n,s2,s2.se)*sin(t*4*pi/cycle.len)  
  return(y)
}

##Simulate point-wise confidence intervals
cihi=rep(NA,nplot)
cilo=rep(NA,nplot)
expected=rep(NA,nplot)
t=seq(0,cycle.len,length.out=nplot)
for(i in 1:nplot){
  y=curvefun.sim(t[i],betalist,cycle.len=cycle.len,n=n)
  tiles=quantile(y,c(0.025,0.975))
  cihi[i]=tiles[2]
  cilo[i]=tiles[1]
  expected[i]=mean(y)
  i=i+1
}

curvedat<-data.frame(cihi=cihi,cilo=cilo,mean=expected)
	if(is.null(ylim)){ylim=c(min(-0.000001,cilo),max(cihi))}
	if(is.null(lty)){lty=1}
	if(is.null(col)){col=2}
	if(is.null(lwd)){lwd=1}
	if(is.null(lty.error)){lty.error=1}
	if(add==F){
	plot(t,curvedat$cihi,type="l",col=col,lty=lty.error,xlab="Daytime",
	ylab="Relative Selection Strength",ylim=ylim,xaxp=c(0,24,12),lwd=lwd,xlim=c(0,24),las=1)
	lines(t,curvedat$mean,col=col,lwd=lwd,lty=lty)
	lines(t,curvedat$cilo,lty=lty.error,col=col,lwd=lwd)
	abline(h=0,lty=3,col=1)
	}
	if(add==T){
	points(t,curvedat$cihi,type="l",col=col,lty=lty.error,xlab="Daytime",
	ylab="Relative Selection Strength",ylim=ylim,lwd=lwd,xaxp=c(0,24,12),las=2,xlim=c(0,24))
	lines(t,curvedat$mean,col=col,lwd=lwd,lty=lty)
	lines(t,curvedat$cilo,lty=lty.error,col=col,lwd=lwd)
	abline(h=0,lty=3,col=1)
	}
	#if(!is.null(type) & type=="circular"){
	#circle.control(type="l")
	#plot.circular(circular(t,unit="hours"),units="hours",shrink=2)
	#lines.circular(circular(t,unit="hours"),curvedat$mean,offset=0.2)
	#}
return(expected)
}

