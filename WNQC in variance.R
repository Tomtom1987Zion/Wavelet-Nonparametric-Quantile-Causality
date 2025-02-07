
## Wavelet Nonparametric Quantile Causality

## by Tomiwa Sunday Adebayo & Oktay ?zkan

## https://doi.org/10.1016/j.jclepro.2023.140321


rm(list=ls(all=TRUE))

#Set working directory
setwd("C:/Users/Admin/Desktop/BRAINBOX/NEW CODESS")

# Load packages

library(tidyquant) #for stock data
library(tidyverse)  # for data analysis
library(crypto2)    #for cryptocurrency download
library(bootstrap)  #for bootstrapping
library(waveslim)   #for wavelet
library(wranglR) #install from Github
library(lattice)   #for plotting
library(viridisLite) #for plotting
library(RColorBrewer) #for plotting
library(QCSIS) #for quantile estimation
library("readxl")
library(quantreg)
library(wesanderson)


# quantile causality test
# x = cause (independent) variable
# y = dependent variable
# q = vector of quatiles, default = 0.01,....,0.99
# moment = postive integer, moment to test for, defaults to 1
# h = bandwidth (postive real) or bandwidth claculation method (string)
#     available methods are:
#       CV = cross validation
#       SJ = Venables and Ripley (2002)
#       YJ = Yu and Jones 1998
#       Silverman = Silverman's plug-in bandwidth
#
#      *** SJ works well usually
#
# code only considers first lags of x and y
#
# Returns:
#    stat = vector t statistics for causality at each quantile
#    q = vector of quantiles
#
# Mehmet Balcilar, 2016-8-19

lrq.causality.test <- function(x,y, hm=c("SJ","CV","YJ","Silverman"), q=NULL,moment=1) {
  if(is.null(q)) {
    qvec <-seq(0.01, 0.99, by = 0.01)
  }
  else {
    qvec <- q
  }
  
  if (moment<0) stop("Negative moment")
  if (trunc(moment) != moment) stop("Noninteger moment")
  
  nq <- length(qvec)
  
  tstatvec <- vector(length=nq, mode="numeric") # initilize the tstat vector
  
  tn <- length(y)-1
  
  yall <- embed(y,2)
  yt1 <- yall[,-1]   # y(t-1), ... y(t-p)
  yt <- yall[,1] # y(t)
  xall <- embed(x,2)
  xt1 <- xall[,-1]
  xt <- xall[,1]
  
  if (moment >= 1) {
    yt1m <- yt1^moment
    ytm <- yt^moment
    xt1m <- xt1^moment
    xtm <- xt^moment
  }
  else{
    stop("moment must be a nonnegative integer")
  }
  
  data.m1 <- data.frame(yt,yt1,xt,xt1)
  data.mhigh <- data.frame(ytm,yt1m,xtm,xt1m)
  
  
  if(hm=="CV") {
    h <- npregbw(ytm~yt1m,regtype="ll")$bw
  }
  else if (hm=="Silverman") {
    h = 1.06*(min(sd(ytm),IQR(ytm)/1.34)/length(ytm)^(1/5))
  }
  else if (hm=="SJ") {
    h = density(ytm,bw="SJ")$bw
  }
  else if (hm=="YJ") {
    h <- dpill(yt1m, ytm, gridsize = tn)
  }
  else{
    h <- hm
  }
  cat("\nbandwith = ", h,"\n")
  
  for  ( jj in 1:nq) {
    qj <- qvec[jj]
    qrh <- h*((qj*(1-qj)/(dnorm(qnorm(p=qj))^2))^(1/5))
    fit <- lprq2(data_m_1=data.m1, data_m_high=data.mhigh, h=qrh, tau=qj, moment=moment)
    iftemp <- (ytm <= fit$fv) - qj
    ifvector <- data.matrix(iftemp)
    kk <- matrix(data = 0, nrow = tn, ncol = tn)
    ymatrix = kronecker(yt1, t(vector(length= tn, mode="numeric")+1))-t(kronecker(yt1, t(vector(length= tn, mode="numeric")+1)))
    wmatrix = kronecker(xt1, t(vector(length= tn, mode="numeric")+1))-t(kronecker(xt1, t(vector(length= tn, mode="numeric")+1)))
    kk=dnorm(ymatrix/qrh)*dnorm(wmatrix/(qrh/sd(yt1)*sd(xt1)))
    tstat <-  t(ifvector)%*%kk%*%ifvector*sqrt(tn/2/qj/(1-qj)/(tn-1)/sum(kk^2)) # Theorem 3.1, Song et al. (2012)
    tstatvec[jj] <- tstat
  }
  return(list(stat=tstatvec,q=qvec))
}


## estimates the quantile of y given x
"lprq2.p" <- function(x, y, h, tau, moment) 
{       
  
  x <- as.matrix(x)
  xx <- as.matrix(x0)
  
  nx <- NROW(xx)
  
  fv <- numeric(nx)
  dv <- numeric(nx)
  
  for(i in 1:nx) {
    z <- x - xx[i,]
    wx <- dnorm(z/h)
    r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
    fv[i] <- r$coef[1.]
    dv[i] <- r$coef[2.]
  }
  list(xx = xx, fv = fv, dv = dv)
}

## estimates the quantile of y given x
"lprq2" <- function(data_m_1, data_m_high, h, tau, moment) # modified from lprq, s.t. we can specify where to estimate quantiles
{       
  
  yt1 <- data_m_1$yt1
  yt <- data_m_1$yt
  yt1m <- data_m_high$yt1m
  ytm <- data_m_high$ytm
  
  x <- yt1
  xx <- yt1  # y(t-1)
  fv <- xx
  dv <- xx
  y <- ytm
  
  for(i in 1:length(xx)) {
    
    z <- x - xx[i]
    wx <- dnorm(z/h)
    z2 <- z
    r <- rq(y~z2, weights=wx, tau=tau, ci=FALSE)
    fv[i] <- r$coef[1.]
    dv[i] <- r$coef[2.]
  }
  list(xx = xx, fv = fv, dv = dv)
}

do.causality.figure <- function(obj1,obj2=NULL,title="") {
  
  #browser()
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  lab_title <- title
  
  nq <- length(obj1$q)
  
  if (is.null(obj2))  {
    pdata <-  data.frame(cbind(q=obj1$q,stats=obj1$stat),CV=rep(1.645,length(obj1$q)))
    p0 <- ggplot(data = pdata, aes(x = q)) +
      geom_line(aes(y=stats,colour="Statistic",
                    linetype="Statistic"), size=2) +
      geom_line(aes(y=CV,colour="CV (5%)",
                    linetype="CV (5%)"), size=1) +
      ggtitle(lab_title)
    #scale_x_continuous(breaks=obj1$q) 
  }
  else {
    pdata <-  data.frame(cbind(q=obj1$q,statsm=obj1$stat,statsv=obj2$stat),CV=rep(1.645,length(obj1$q)))
    p0 <- ggplot(data = pdata, aes(x = q)) +
      geom_line(aes(y=statsm, colour="Mean",
                    linetype="Mean"), size=2) +
      geom_line(aes(y=statsv,colour="Variance",
                    linetype="Variance"), size=2) +
      geom_line(aes(y=CV, colour="CV (5%)",
                    linetype="CV (5%)"), size=1) +
      ggtitle(lab_title)
    #scale_x_continuous(breaks=obj1$q)
  }
  
  p1 <- p0 +
    xlab("Kantiller") +
    ylab("Test  statisti i") +
    theme_bw() +
    scale_size(range=c(1, 2, 2)) +
    theme(legend.position = "none") + 
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(angle = 0, size=12)) +
    theme(axis.text.y = element_text(angle = 0, size=12)) +
    theme(plot.margin = unit(c(5, 5, 0, 0), "mm")) + 
    theme(legend.spacing = unit(0, "cm")) +
    scale_linetype_manual(values=c("solid", "solid", "solid")) +
    scale_color_manual(values=cbPalette[c(1,2,3)])
  
  return(p1)  
}



# Loading DATA, where it is located and name of file.
Data<- read_excel("DATA.xlsx")
attach(Data)

DEP <- DEP
IND <- IND


## To install WranglR
## install.packages("remotes")
## remotes::install_github("phively/wranglR")
# all other packages can be installed from CRAN normally
#estimation
# you can change wavelet parameters in d1 and d2 as you require

## moment=1 (causality in mean) and moment=2 (causality in variance)
## Note: moment=2 works well with % series not first difference series.


## if u get singularity matrix error (this is because of the Balcilar et al.'s (2016) method error) try to use DEP <- DEP*10 IND <- IND*10 or DEP*100 IND <- IND*100 or DEP*1000 IND <- IND*1000. These do not cause a big change in the estimates. That is to get rid of the singularity matrix error.

tau <- seq(0.05, 0.95, 0.05)
d1=as.data.frame(mra(DEP, wf = "la8", J = 5, method = "modwt", boundary = "periodic"))
d2=as.data.frame(mra(IND, wf = "la8", J = 5, method = "modwt", boundary = "periodic"))
cols <- intersect(colnames(d1), colnames(d2))
res <- lapply(cols, function(x) {
  model <-lrq.causality.test(x=d2[, x], y=d1[, x], moment=2, hm="Silverman",q=tau)
  result <- list(tau = tau, rho = model$stat)
  return(result)
})
names(res) <- cols
op=cbind.data.frame(tau,ListExtract(res,"rho"))

## Plotting
plot=op[,c(-1,-10)] 
plot=as.matrix(plot)
row.names(plot)=tau   

a <- plot[,c(1,2,3)]

## you can use any results you want to use. for example if you want to use results of 1,2,3 use this: a <- plot[,c(1,2,3)] and arrange the following colnames as colnames(a)=c("D1", "D2", "D3"). you can also change colnames as colnames(a)=c("Short", "Medium", "Long")

colnames(a)=c("Short", "Medium", "Long")

pallete <- colorRampPalette(c("darkseagreen1","orange","red"))(20)

levelplot(a,
          col.regions = pallete,
          
          panel = function(...) {
            panel.levelplot(...)
            for (i in 1:nrow(a)) {
              for (j in 1:ncol(a)) {
                if (a[i, j] >= 1.96) {
                  panel.text(i, j, "**", col = "black", cex = 1.5)
                } else if (a[i, j] >= 1.645 & a[i, j] < 1.96) {
                  panel.text(i, j, "*", col = "black", cex = 1.5) 
                }
              }
            }
          },
          xlab = "Quantiles",
          ylab = "Periods",
          main = ""
)