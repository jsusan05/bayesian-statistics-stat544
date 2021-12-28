#######################################
# STAT 544: Final Project 
# jsusan@iastate.edu
#######################################
library(reshape2)
library(plyr)
library(ggplot2)
library(xtable)
library(rstan)
library(igraph)
source("multiplot.R")

set.seed(1)
rstan_options(auto_write = TRUE)

d = read.csv("./data/data_epa.csv",header=TRUE)
d$County=factor(d$County)
d$Time=factor(d$Time)
d$wind 
d$tmax
d$pm25
dd <- read.csv("./data/county_graph.csv")
gg <- graph.data.frame(dd, directed=FALSE)
g=simplify(gg)
plot(g,)
A=as_adjacency_matrix(g)

############################
# Data plots
###########################
df = rbind(data.frame(county = factor(d$County),         
                      time = d$Time, 
                      wind = d$wind,
                      tmax = d$tmax,
                      sites = d$ToxicSites,
                      pm25 = d$pm25)) # Each group has its own mean

#explanatory variables
p1<-ggplot(df, aes(x=wind, y=pm25, color=time)) + 
  geom_point() + scale_colour_discrete(guide = FALSE)+
  facet_wrap(~county,ncol=3) +
  theme_bw()
p2<-ggplot(df, aes(x=tmax, y=pm25, color=time)) + 
  geom_point() + scale_colour_discrete(guide = FALSE)+
  facet_wrap(~county,ncol=3) +
  theme_bw()
p3<-ggplot(df, aes(x=wind, y=pm25, color=county)) + 
  geom_point() + scale_colour_discrete(guide = FALSE)+
  facet_wrap(~time,ncol=4) +
  theme_bw()
p4<-ggplot(df, aes(x=tmax, y=pm25, color=county)) + 
  geom_point() + scale_colour_discrete(guide = FALSE)+ 
  facet_wrap(~time,ncol=4) +
  theme_bw()
multiplot(p1, p2, p3, p4, cols=2)
#Space data
p1<-ggplot(df, aes(x=county, y=pm25)) + 
  geom_point() + 
  facet_wrap(~time,ncol=2) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2<-ggplot(df, aes(x=county, y=wind)) + 
  geom_point() + 
  facet_wrap(~time,ncol=2) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3<-ggplot(df, aes(x=county, y=tmax)) + 
  geom_point() + 
  facet_wrap(~time,ncol=2) +
  theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
multiplot(p1, p2, p3, cols=3)
#Time data
p1<-ggplot(df, aes(x=time, y=pm25, group = 1)) + 
  geom_point() + geom_line()+ 
  facet_wrap(~county,ncol=2) +
  theme_bw()
p2<-ggplot(df, aes(x=time, y=wind, group = 1)) + 
  geom_point() + geom_line()+
  facet_wrap(~county,ncol=2) +
  theme_bw()
p3<-ggplot(df, aes(x=time, y=tmax, group = 1)) + 
  geom_point() + geom_line()+ 
  facet_wrap(~county,ncol=2) +
  theme_bw()
multiplot(p1, p2, p3, cols=3)

#m = stan_model(model_code = ar_spatial_mod)
M=model.matrix(pm25~tmax+wind,data=d)
X=matrix(M,nrow = nrow(M),ncol=ncol(M),byrow = FALSE)
Zm = model.matrix(~0+Time,data=d)
Z1=matrix(Zm,nrow = nrow(Zm),ncol=ncol(Zm),byrow = FALSE)
Zm = model.matrix(~0+County,data=d)
Z2=matrix(Zm,nrow = nrow(Zm),ncol=ncol(Zm),byrow = FALSE)

W=matrix(A[levels(d$County),levels(d$County),drop=FALSE],nrow=nrow(A),ncol = ncol(A),byrow = TRUE)
m_diag = c()
for (i in 1:nrow(W))
  m_diag[i]=sum(W[i,])
   #m_diag[i]=unique(d$ToxicSites[row.names(A[levels(d$County),levels(d$County),drop=FALSE])[i]==d$County])
m_diag
# D=diag(m_diag)

#D_m05=diag(sqrt(1/m_diag))
#r_b = 1/min(eigen(D_m05%*%W%*%D_m05)$values) #not <1

dat = list(n = nrow(X),               #Number of samples
           p = ncol(X),               #Number of fixed parameters
           t = nlevels(d$Time),       #Number of time points
           c = nlevels(d$County),     #Number of counties
           X = X,                     #Model Matrix for fixed component
           Z1=Z1,                     #Model Matrix for time component
           Z2=Z2,                     #Model Matrix for county component
           W=W,                       #Neighbor Matrix from county graph
           y = d$pm25                #Average of observed pm2.5 matter measured
)                

#Linear Regression Model
full_fit <- stan('stan/epa_lm.stan', data = dat)
print(full_fit, pars = c('beta', 'sigma', 'lp__'))
to_plot <- c('beta', 'mu', 'sigma')
traceplot(full_fit, pars = to_plot)

#Linear Regression Temporal Model
full_fit <- stan('stan/epa_lm_temporal.stan', data = dat)
print(full_fit, pars = c('beta', 'sigma', 'omega', 'phi', 'sigma_w', 'lp__'))
to_plot <- c('beta', 'sigma', 'omega', 'phi', 'sigma_w')
traceplot(full_fit, pars = to_plot)

#Linear Regression Spatial Model
full_fit <- stan('stan/epa_lm_spatial.stan', data = dat)
print(full_fit, pars = c('beta', 'sigma', 'rho','mu','tau', 'theta', 'lp__'))
to_plot <- c('beta', 'sigma', 'rho', 'tau', 'theta')
traceplot(full_fit, pars = to_plot)

#Linear Regression Temporal Spatial Model
full_fit <- stan('stan/epa_lm_temporal_spatial.stan', data = dat)
print(full_fit, pars = c('beta', 'sigma', 'omega', 'theta', 'phi', 'sigma_w', 'tau','rho','mu','lp__'))
to_plot <- c('beta', 'sigma', 'omega', 'theta', 'phi', 'sigma_w', 'tau','rho')
traceplot(full_fit, pars = to_plot)

############################
# Prior Distribution plots
###########################
library("dplyr")
library("tidyr")
library("ggplot2")
source("multiplot.R")
par(mfrow=c(1,1))
#beta
prec=2
sd_b = sqrt(1/prec)
beta_sample = rstan::extract(full_fit)$beta
sim = data.frame(x = as.numeric(beta_sample))
p1<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dnorm, args=list(mean=0, sd=sd_b)) +
  theme_bw()+ggtitle("beta")
#sigma
sigma_sample = rstan::extract(full_fit)$sigma
sim = data.frame(x = as.numeric(sigma_sample))
p2<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dcauchy, args=list(location=0, scale=5)) +
  theme_bw()+ggtitle("sigma")
#phi
phi_sample = rstan::extract(full_fit)$phi
sim = data.frame(x = as.numeric(phi_sample))
p3<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dnorm, args=list(mean=0, sd=1)) +
  theme_bw()+ggtitle("phi")
#sigma_w
sigmaw_sample = rstan::extract(full_fit)$sigma_w
sim = data.frame(x = as.numeric(sigmaw_sample))
p4<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dcauchy, args=list(location=0, scale=5)) +
  theme_bw()+ggtitle("sigma_w")
#tau
tau_sample = rstan::extract(full_fit)$tau
sim = data.frame(x = as.numeric(tau_sample))
p5<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dgamma, args=list(shape=2, rate=2)) +
  theme_bw()+ggtitle("tau")
#rho
rho_sample = rstan::extract(full_fit)$rho
sim = data.frame(x = as.numeric(rho_sample))
p6<-ggplot(sim, aes(x=x)) + 
  geom_histogram(aes(y=..density..), binwidth=.01) + 
  stat_function(color = 'red', fun=dunif, args=list(min=0, max=1)) +
  theme_bw()+ggtitle("rho")

multiplot(p1, p2, p3, p4, p5, p6, cols=3)

################################
# Posterior distribution plots
################################
plot(full_fit, pars=c('theta'))
plot(full_fit, pars=c('omega'))
plot(full_fit, pars=c('mu'))

theta_df=as.data.frame(summary(full_fit, pars = c("theta"), probs = c(0.5))$summary)
theta_df[,"50%"]
theta_df[,"sd"]

omega_df=as.data.frame(summary(full_fit, pars = c("omega"), probs = c(0.5))$summary)
omega_df[,"50%"]
omega_df[,"sd"]

mu_df=as.data.frame(summary(full_fit, pars = c("mu"), probs = c(0.5))$summary)
mu_df[,"50%"]
mu_df[,"sd"]

omega_sample = as.matrix(rstan::extract(full_fit)$omega)
res <- melt(omega_sample)
ddf = data.frame(time=res$Var2,omega=res$value)
ggplot(ddf, aes(x=as.factor(time), y=omega)) + geom_boxplot() +
  xlab("time")+
  theme_bw()

theta_sample = as.matrix(rstan::extract(full_fit)$theta)
res <- melt(theta_sample)
ddf = data.frame(time=res$Var2,omega=res$value)
ggplot(ddf, aes(x=as.factor(time), y=omega)) + geom_boxplot() +
  xlab("county")+scale_x_discrete(labels=levels(d$County))+
  theme_bw()

library(choroplethr)
library(choroplethrMaps)
data(df_pop_county) 
head(df_pop_county) 
dfdmed=data.frame(region=c(19045,19055,19103,19113,19139,19163),value=theta_df[,"50%"])
dfdsd=data.frame(region=c(19045,19055,19103,19113,19139,19163),value=theta_df[,"sd"])
dfd_mumed=data.frame(county=levels(d$County)[rep(c(1:6),each=8)],time=rep(c(1:8),times=6),median=exp(mu_df[,"50%"]),sd=mu_df[,"sd"])
county_choropleth(dfdmed, state_zoom = "iowa")
county_choropleth(dfdsd, state_zoom = "iowa")

p1<-ggplot(dfd_mumed, aes(x=time, y=median, group = 1)) + 
  geom_point() + geom_line()+ scale_x_discrete(limits = c(1:8))+
  facet_wrap(~county,ncol=2) +ylab("mu_median")+
  theme_bw()
p2<-ggplot(dfd_mumed, aes(x=county, y=median, group = 1)) + 
  geom_point() + geom_point()+
  facet_wrap(~time,ncol=2) +ylab("mu_median")+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
multiplot(p1, p2, cols=2)

mu_sample = as.matrix(rstan::extract(full_fit)$mu)[1:20,]
ddf=rbind(data.frame(county=rep(levels(d$County)[1],nrow(mu_sample)),mu=mu_sample[,1:8]),
          data.frame(county=rep(levels(d$County)[2],nrow(mu_sample)),mu=mu_sample[,9:16]),
          data.frame(county=rep(levels(d$County)[3],nrow(mu_sample)),mu=mu_sample[,17:24]),
          data.frame(county=rep(levels(d$County)[4],nrow(mu_sample)),mu=mu_sample[,25:32]),
          data.frame(county=rep(levels(d$County)[5],nrow(mu_sample)),mu=mu_sample[,33:40]),
          data.frame(county=rep(levels(d$County)[6],nrow(mu_sample)),mu=mu_sample[,41:48]))
res <- melt(mu_sample)
ddf = data.frame(time=res$Var2,mu=res$value)
ggplot(ddf, aes(x=as.factor(time), y=mu)) + geom_boxplot() +
  xlab("county")+scale_x_discrete(labels=levels(d$County))+
  theme_bw()

#Model Checking
####################################
# Problem2: Model Checking
####################################
library(data.table)
nreps=19

df = rbind(data.frame(county = factor(d$County),         
                      y = d$pm25,
                      t = d$Time,
                      group=as.numeric(factor(d$County)),
                      nsites = d$ToxicSites))

#Sample mu posterior values
mu_sample = rstan::extract(full_fit)$mu
sigma_sample = rstan::extract(full_fit)$sigma
#########################
# Visually observe data
#########################
rep = rdply(nreps, { 
  tmp=df
  lvec=c()
  n=nrow(df)
  id = sample(nreps,1)
  tmp$y = exp(rnorm(n,mu_sample[id, df$group],sigma_sample[id]))
  tmp 
})

id = sample(20,1) # Choose a random spot to insert observed data set
if (id==20) {
  df$.n = 20
  rep = rbind(rep,df)
} else {
  rep$.n[rep$.n==id] = 20
  df$.n = id
  rep = rbind(rep,df) 
}

ggplot(rep, aes(x=y, fill=factor(1))) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~.n) +
  theme_bw() +
  theme(legend.position="none")

ggplot(rep, aes(x=y, fill=.n==id)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~.n) +
  theme_bw() +
  theme(legend.position="none")



