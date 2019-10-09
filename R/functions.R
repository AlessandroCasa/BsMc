dnorm.deriv <- function (x, mu = 0, sigma = 1, deriv.order = 0)
{
  r <- deriv.order
  phi <- dnorm(x, mean = mu, sd = sigma)
  x <- (x - mu)
  arg <- x/sigma
  hmold0 <- 1
  hmold1 <- arg
  hmnew <- 1
  if (r == 1)
    hmnew <- hmold1
  if (r >= 2)
    for (i in (2:r)) {
      hmnew <- arg * hmold1 - (i - 1) * hmold0
      hmold0 <- hmold1
      hmold1 <- hmnew
    }
  derivt <- (-1)^r * phi * hmnew/sigma^r
  return(derivt)
}

dnorm.deriv.mixt <- function (x, mus = 0, sigmas = 1, props = 1, deriv.order = 0)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("Proportions don't sum to one")
  if (identical(all.equal(props[1], 1), TRUE))
    dens <- dnorm.deriv(x, mu = mus[1], sigma = sigmas[1],
                        deriv.order = deriv.order)
  else {
    k <- length(props)
    dens <- 0
    for (i in 1:k) dens <- dens + props[i] * dnorm.deriv(x = x,
                                                         mu = mus[i], sigma = sigmas[i], deriv.order = deriv.order)
  }
  return(dens)
}

pnorm.mixt <- function (x, mus = 0, sigmas = 1, props = 1)
{
  if (!(identical(all.equal(sum(props), 1), TRUE)))
    stop("proportions don't sum to one")
  if (identical(all.equal(props[1], 1), TRUE))
    cdf <- pnorm(x, mean = mus[1], sd = sigmas[1])
  else {
    k <- length(props)
    cdf <- 0
    for (i in 1:k) cdf <- cdf + props[i] * pnorm(x, mean = mus[i],
                                                 sd = sigmas[i])
  }
  return(cdf)
}

psi<-function(mu,sigma){
  r<-mu/sigma
  value<-2*sigma*dnorm(r)+mu*(1-2*pnorm(-r))
  return(value)
}


extremes.kde <- function (x, h, plot=FALSE)
{
  eps=10^(-8)
  fderiv.hat<-function(eval.points){
    val<-kdde(x=x,h=h,eval.points=eval.points,deriv.order=1)$estimate
    return(val)
  }
  exts.hat<-uniroot.all(f=fderiv.hat,lower=min(x)-h,upper=max(x)+h,tol=eps)
  exts.hat <- exts.hat[order(exts.hat)]
  maxs.hat <- exts.hat[which((1:length(exts.hat))%%2!=0)]
  nmodes<-length(maxs.hat)
  mins.hat<-numeric(nmodes-1)
  if(nmodes>1){ mins.hat <- exts.hat[which((1:length(exts.hat))%%2==0)] }

  if(plot){
    xs<-seq(min(x)-3*h,max(x)+3*h,length=500)
    plot(xs,kde(x=x,h=h,eval.points=xs)$estimate,t="l",lwd=2)
    if(length(mins.hat)>0){
      points(mins.hat,kde(x=x,h=h,eval.points=mins.hat)$estimate,pch=19,col=4,lwd=2)
    }
    points(maxs.hat,kde(x=x,h=h,eval.points=maxs.hat)$estimate,pch=19,col=2,lwd=2)
  }
  return(list(maxim=maxs.hat,minim=mins.hat))
}

extremes.mixt <- function (mus = 0, sigmas = 1, props = 1, plot=FALSE)
{
  eps = 10^(-8)
  fderiv.hat<-function(eval.points){
    #val<-ks:::dnorm.deriv.mixt(x=eval.points,mus=mus,sigmas=sigmas,props=props,deriv.order=1)
    val<-dnorm.deriv.mixt(x=eval.points,mus=mus,sigmas=sigmas,props=props,deriv.order=1)
    return(val)
  }
  exts.hat<-uniroot.all(f=fderiv.hat,
                        lower=min(mus)-3*sigmas[which.min(mus)],upper=max(mus)+3*sigmas[which.max(mus)],
                        tol=eps)
  exts.hat <- exts.hat[order(exts.hat)]
  maxs.hat <- exts.hat[which((1:length(exts.hat))%%2!=0)]
  nmodes<-length(maxs.hat)
  mins.hat<-numeric(nmodes-1)
  if(nmodes>1){ mins.hat <- exts.hat[which((1:length(exts.hat))%%2==0)] }

  if(plot){
    xs<-seq(min(mus)-3*sigmas[which.min(mus)],max(mus)+3*sigmas[which.max(mus)],length=1000)
    plot(xs,dnorm.mixt(x=xs,mus=mus,sigmas=sigmas,props=props),t="l",lwd=2)
    if(length(mins.hat)>0){
      points(mins.hat,dnorm.mixt(x=mins.hat,mus=mus,sigmas=sigmas,props=props),pch=19,col=4,lwd=2)
    }
    points(maxs.hat,dnorm.mixt(x=maxs.hat,mus=mus,sigmas=sigmas,props=props),pch=19,col=2,lwd=2)
  }
  return(list(maxim=maxs.hat,minim=mins.hat))
}

#' Compute the Distance in Measure between the clustering induced by a kernel density estimator (based on a sample x and a bandwidth h) and the population clustering defined by a K-component normal mixture density.
#'
#' @param x (vector) the data to be partitioned.
#' @param h the bandwidth to be used to estimate the density via KDE.
#' @param mus vector of means of the mixture components.
#' @param sigmas vector of standard deviations of the mixture components.
#' @param props vector of mixing proportions of the mixture components.
#' @param plot if true, the true density and the estimated one are displayed.
#' @return the value of the Distance in Measure.
#' @references Chacón, J.E. (2015). A population background for nonparametric density-based clustering. Statistical Science 30(4): 518-532.

dist.meas <- function (x, h, mus=0, sigmas=1, props=1, plot=FALSE)
{
  n<-length(x)
  mh<-extremes.kde(x=x,h=h)$minim
  mf<-extremes.mixt(mus=mus,sigmas=sigmas,props=props)$minim
  mh<-c(-Inf,mh,Inf)
  mf<-c(-Inf,mf,Inf)
  Ih<-length(mh)-1
  If<-length(mf)-1
  prob.sym.diff<-matrix(0,nrow=Ih,ncol=If)
  # A matrix containing the probabilities of all the symmetric differences between clusters of f and clusters of f_h
  for(i in 1:Ih){for(j in 1:If){
    pts<-sort(c(mh[i],mh[i+1],mf[j],mf[j+1]))
    # The symmetric difference between two intervals I and J is the union of two intervals (a,b) and (c,d)
    # where (a,b,c,d) is the set of all the four extremes of I and J, sorted in increasing order
    prob.sym.diff[i,j]<-pnorm.mixt(pts[2],mus=mus,sigmas=sigmas,props=props)-pnorm.mixt(pts[1],mus=mus,sigmas=sigmas,props=props)+
      pnorm.mixt(pts[4],mus=mus,sigmas=sigmas,props=props)-pnorm.mixt(pts[3],mus=mus,sigmas=sigmas,props=props)
  }}
  if(Ih!=If){
    if(Ih<If){
      prob.clus.f<-diff(pnorm.mixt(x=sort(mf),mus=mus,sigmas=sigmas,props=props)) # The probabilities of each cluster of f
      prob.sym.diff<-rbind(prob.sym.diff,matrix(rep(prob.clus.f,If-Ih),nrow=If-Ih,ncol=If,byrow=TRUE))
      # Those probabilities are appended to the previous matrix, as the result of extending the clusters of f_h with If-Ih empty sets
    }
    if(If<Ih){
      prob.sym.diff<-t(prob.sym.diff)
      prob.clus.fh<-diff(pnorm.mixt(x=sort(mh),mus=mus,sigmas=sigmas,props=props)) # The probabilities of each cluster of f_h
      prob.sym.diff<-rbind(prob.sym.diff,matrix(rep(prob.clus.fh,Ih-If),nrow=Ih-If,ncol=Ih,byrow=TRUE))
      # Those probabilities are appended to the previous matrix, as the result of extending the clusters of f with Ih-If empty sets
      prob.sym.diff<-t(prob.sym.diff)
    }
  }
  lsap<-solve_LSAP(x=prob.sym.diff) # Solve the linear sum assignment problem
  distance<-sum(prob.sym.diff[cbind(seq_along(lsap), lsap)])/2 # Compute the optimal value
  if(plot){
    xs<-seq(min(c(x,mus))-3*sigmas[which.min(mus)],max(c(x,mus))+3*sigmas[which.max(mus)],length=1000)
    plot(xs,dnorm.mixt(xs,mus=mus,sigmas=sigmas,props=props),t="l",xlab="",ylab="",lwd=2)
    lines(xs,rep(0,1000),col="grey")
    lines(density(x,bw=h))
  }
  return(distance)
}

#' Compute the Expected Distance in Measure (EDM) between a kernel estimator-induced partition and the population one defined by a K-component normal mixture density over a grid of bandwidths.
#'
#' @param n sample size of the simulated Monte Carlo samples.
#' @param hmin lower value for the grid of bandwidths.
#' @param hmax upper value for the grid of bandwidths.
#' @param byh increment of the grid sequence
#' @param mus vector of means of the mixture components.
#' @param sigmas  vector of standard deviations of the mixture components.
#' @param props vector of mixing proportions of the mixture components.
#' @param plot if TRUE the curve of the EDM as a function of bandwidth values in the grid is displayed along with the single DM for each Monte Carlo sample.
#' @param B the number of Monte Carlo samples.
#' @param verbose if TRUE, the computational progression is displayed.
#' @return the value of the EDM for each value in the grid.
#' @references Chacón, J.E. (2015). A population background for nonparametric density-based clustering. Statistical Science 30(4): 518-532.

Edist.meas <- function (n=100, hmin=0.05, hmax=1, byh=hmin, mus=0, sigmas=1,
                        props=1, plot=TRUE, B=100, verbose=TRUE)
{
  hs<-seq(hmin,hmax,by=byh)
  nhs<-length(hs)
  dmh<-matrix(0,nrow=B,ncol=nhs)
  for(i in 1:B){
    if(verbose){cat(i," ")}
    set.seed(i)
    x<-rnorm.mixt(n=n,mus=mus,sigmas=sigmas,props=props)
    for(j in 1:nhs){
      dmh[i,j]<-dist.meas(h=hs[j],x=x,mus=mus,sigmas=sigmas,props=props)
    }
  }
  cat("\n")
  Edmh<-colMeans(dmh)
  if(plot){
    plot(hs,hs,col=0,ylim=c(0,max(dmh)),xlab="h",ylab="Expected distance in measure")
    for(i in 1:B){
      lines(hs,dmh[i,],t="l",col="grey")
    }
    lines(hs,Edmh,t="l",lwd=2,col=1)
  }
  return(list(hs=hs,Edmh=Edmh))
}

#' Compute the Expected Distance in Measure (EDM) between a kernel estimator-induced partition and the population one defined by a K-component normal mixture density for a single bandwidth value.
#'
#' @param n sample size of the simulated Monte Carlo samples.
#' @param h value of the bandwidth for which the EDM has to be computed.
#' @param mus vector of means of the mixture components.
#' @param sigmas  vector of standard deviations of the mixture components.
#' @param props vector of mixing proportions of the mixture components.
#' @param B the number of Monte Carlo samples.
#' @return the value of the EDM for the given value \verb{h}
#' @references Chacón, J.E. (2015). A population background for nonparametric density-based clustering. Statistical Science 30(4): 518-532.

Edist.meas_singleh <- function (n=100, h, mus=0, sigmas=1, props=1, B=100)
{
  dmh<-matrix(0,nrow=B,ncol=1)
  for(i in 1:B){
    set.seed(i)
    x<-rnorm.mixt(n=n,mus=mus,sigmas=sigmas,props=props)
    dmh[i,1]<-dist.meas(h=h,x=x,mus=mus,sigmas=sigmas,props=props)
  }
  Edmh<-colMeans(dmh)
  return(list(hs=h,Edmh=Edmh,dmh=dmh))
}


#' Asymptotic Expected Distance in Measure (AEDM)
#' @description Compute the Asymptotic Expected Distance in Measure (AEDM) between a kernel estimator-induced partition and the population one defined by a K-component normal mixture density over a grid of bandwidths.
#' @param n sample size of the simulated Monte Carlo samples.
#' @param hmin lower value for the grid of bandwidths.
#' @param hmax upper value for the grid of bandwidths.
#' @param byh increment of the grid sequence
#' @param mus vector of means of the mixture components.
#' @param sigmas  vector of standard deviations of the mixture components.
#' @param props vector of mixing proportions of the mixture components.
#' @param plot if TRUE the curve of the AEDM as a function of bandwidth values in the grid is displayed along with the single DM for each Monte Carlo sample.
#' @return the value of the EDM for each value in the grid.
#' @references Casa A., Chacón, J.E. and Menardi, G. (2019). Modal clustering asymptotics with applications to bandwidth selection (https://arxiv.org/pdf/1901.07300.pdf).

AEdist.meas <- function (n=100, hmin=0.05, hmax=1, byh=hmin, mus=0,
                         sigmas=1, props=1, plot=TRUE)
{
  hs<-seq(hmin,hmax,by=byh)
  nhs<-length(hs)
  betas<-hs*n^(1/7)

  mins<-extremes.mixt(mus=mus,sigmas=sigmas,props=props)$minim
  fmins<-dnorm.mixt(mins,mus=mus,sigmas=sigmas,props=props)
  #f2mins<-ks:::dnorm.deriv.mixt(mins,mus=mus,sigmas=sigmas,props=props,deriv.order=2)
  f2mins<-dnorm.deriv.mixt(mins,mus=mus,sigmas=sigmas,props=props,deriv.order=2)
  #f3mins<-ks:::dnorm.deriv.mixt(mins,mus=mus,sigmas=sigmas,props=props,deriv.order=3)
  f3mins<-dnorm.deriv.mixt(mins,mus=mus,sigmas=sigmas,props=props,deriv.order=3)
  mu2K<-1
  RK1<-1/(4*sqrt(pi))

  AEhs<-rep(0,length(hs))
  for(i in 1:length(mins)){
    AEhs<-AEhs+n^(-2/7)*(fmins[i]/f2mins[i])*psi(mu=betas^2*f3mins[i]*mu2K/2,
                                                 sigma=sqrt(RK1*fmins[i]/betas^3))
  }

  if(plot){
    plot(hs,hs,col=0,ylim=c(0,max(AEhs)),xlab="h",ylab="AEDM")
    lines(hs,AEhs,t="l",lwd=2,col=1)
  }
  return(list(hs=hs,aedm=AEhs))
}


#' Bandwidth \verb{h_AB1} (3.6)
#' @description Compute the optimal bandwidth \verb{h_AB1} (3.6) in Casa et al. (2019).
#' @param data data
#' @return a value for the modal clustering-oriented optimal bandwidth
#' @references Casa A., Chacón, J.E. and Menardi, G. (2019). Modal clustering asymptotics with applications to bandwidth selection (https://arxiv.org/pdf/1901.07300.pdf).

hab1 <- function(data) {
  n <- length(data)
  mu2K<-1
  RK1<-1/(4*sqrt(pi))

  bw <- hpi(data,deriv.order = 1)
  minima <- extremes.kde(data,h=bw)$minim

  if (length(minima)==0) {
    out <- bw.crit(data)
  } else {

    min_obj <- list(minima=minima,
                    fmins=kde(data, h=hpi(data,deriv.order = 0), eval.points = minima)$estimate,
                    f2mins=kdde(data, h=hpi(data,deriv.order = 2), deriv.order=2,eval.points=minima)$estimate,
                    f3mins=kdde(data, h=hpi(data,deriv.order = 3), deriv.order=3,eval.points=minima)$estimate)

    mins <- min_obj$minima
    nmin <- length(mins)
    fmins <- min_obj$fmins
    if (any(fmins<0)) {fmins[fmins<0] <- 0}
    f2mins <- min_obj$f2mins
    f3mins <- min_obj$f3mins

    f_ratio <- sum(sapply(1:nmin,function(y) (fmins[y]^(3/2))/f2mins[y]))
    f_ratio2 <- sum(sapply(1:nmin,function(y) (fmins[y]*abs(f3mins[y]))/f2mins[y]))

    out <- (n^(-1/7))*((3*sqrt(RK1)*f_ratio)/(mu2K*sqrt(2*pi)*f_ratio2))^(2/7)
  }
  return(out)
}

#' Bandwidth \verb{h_AB2} (3.7)
#' @description Compute the optimal bandwidth \verb{h_AB2} (3.7) in Casa et al. (2019).
#' @param data data
#' @return a value for the modal clustering-oriented optimal bandwidth
#' @references Casa A., Chacón, J.E. and Menardi, G. (2019). Modal clustering asymptotics with applications to bandwidth selection (https://arxiv.org/pdf/1901.07300.pdf).

hab2 <- function(data) {
  n <- length(data)
  mu2K<-1
  RK1<-1/(4*sqrt(pi))

  bw <- hpi(data,deriv.order = 1)
  minima <- extremes.kde(data,h=bw)$minim

  if (length(minima)==0) {
    out <- bw.crit(data)
  } else {
    min_obj <- list(minima=minima,
                    fmins=kde(data, h=hpi(data,deriv.order = 0), eval.points = minima)$estimate,
                    f2mins=kdde(data, h=hpi(data,deriv.order = 2), deriv.order=2,eval.points=minima)$estimate,
                    f3mins=kdde(data, h=hpi(data,deriv.order = 3), deriv.order=3,eval.points=minima)$estimate)

    mins <- min_obj$minima
    nmin <- length(mins)
    fmins <- min_obj$fmins
    if (any(fmins<0)) {fmins[fmins<0] <- 0}
    f2mins <- min_obj$f2mins
    f3mins <- min_obj$f3mins

    f_ratio <- sum(sapply(1:nmin,function(y) (fmins[y]^(3/2))/f2mins[y]))
    f_ratio2 <- sum(sapply(1:nmin,function(y) (sqrt(fmins[y])*f3mins[y]^2)/f2mins[y]))

    out <- (n^(-1/7))*((24*RK1*f_ratio)/(11*f_ratio2*mu2K^2))^(1/7)
  }
  return(out)
}





