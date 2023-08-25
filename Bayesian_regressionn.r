library(rethinking)


R2 <- function( quap_fit, dlist, xvar ) {
s <- rethinking::link( quap_fit, data=dlist )
r <- apply(s,2,mean) - xvar
1 - var2(r)/var2(xvar)
}
backScale<-function(x,y){
return((x*sd(y))+mean(y))

}
data<-readShapePoly("C:/Users/seboc/Box/Parasite_S_gradients_NA/Ecoregions_Sep27_wgs84.shp")
data<-data@data[,c(4,38:42)]

er<-readShapePoly("C:/Users/seboc/Box/Parasite_S_gradients_NA/Ecoregions_Sep27_wgs84.shp")
er1<-er@data
er1<-er1[which(er1$SampleCove>0.7 & er1$H_S_Cov < 0.6),]
er1$Rar90SP<-er1$estimator_
er1$Rar90SP_Lo<-er1$asym_Lower
er1$Rar90SP_Up<-er1$asym_Uper
#er1<-merge(er1[,-c(38:43)],data,by="ECO_NAME")

dlist <- list( 
Pa_obs =er1$Rar90SP,
D_sd = ((er1$Rar90SP - er1$Rar90SP_Lo)/1.96 ),
Lat = scale(er1$Lat),
Lat2=scale(er1$Lat)^2,
Host = scale(er1$HostS),
Host2 = scale(er1$HostS) ^ 2,
Temp = scale(er1$Temperatur),
TemSeason = scale(er1$Temperatur.1),
Prec = scale(er1$Precipitat),
PercSeason=scale(er1$Precipitat.1),
Area =scale(er1$EcoRegionA),
N = nrow(er1)
)



#### Latitude
Lat <- ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Lat ,
a ~ dnorm(100,30),
bL ~ dnorm( 0 , 10),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Lat)
precis(Lat , depth=2 )



setwd("C:/Users/seboc/Box/Parasite_S_gradients_NA/New_Analyses/Bayesian_model_Figs")
tiff(filename = "Rodent_Asym_Lat.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
plot(er1$Lat, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Latitude" , ylab="Species richness-Asypmtotic est",cex=3,cex.lab=2.8,cex.axis=2.5 )

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Lat[i]
lines( c(x,x) , ci,lwd=2 )
}

p <- rethinking::link( Lat , data=data.frame( Lat=seq(-4,4, by=0.1) ) )
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-4,4, by=0.1),er1$Lat), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-4,4, by=0.1),er1$Lat), col=col.alpha(rangi2,0.3) )

dev.off()
R2(Lat,dlist, er1$Rar90SP)

#### Latitude quadr

Lat2 <- ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Lat + bL2*Lat2 ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 10),
bL2~ dnorm( 0 , 10),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Lat2)
precis(Lat2 , depth=2 )



setwd("C:/Users/seboc/Box/Parasite_S_gradients_NA/New_Analyses/Bayesian_model_Figs")
tiff(filename = "Rodent_Asym_LatQ.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$Lat, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Latitude" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Lat[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Lat2 , data=data.frame( Lat=seq(-4,4, by=0.1), Lat2=seq(-4,4, by=0.1)^2 ) )
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-4,4, by=0.1),er1$Lat), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-4,4, by=0.1),er1$Lat), col=col.alpha(rangi2,0.3) )

dev.off()
R2(Lat2,dlist, er1$Rar90SP)



##### Host
Host<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Host ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Host)
precis(Host , depth=2 )

tiff(filename = "Rodent_Asym_Host.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$HostS, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Host richness" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$HostS[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Host , data=data.frame( Host=seq(-4,4, by=0.1), Lat2=seq(-4,4, by=0.1)^2 ) )
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-4,4, by=0.1),er1$HostS), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-4,4, by=0.1),er1$HostS), col=col.alpha(rangi2,0.3) )

dev.off()


R2(Host,dlist, er1$Rar90SP)


##### Host Quadratic
Host2<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Host + bL2*Host2 ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 10),
bL2~ dnorm( 0 , 10),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Host2)
precis(Host2 , depth=2 )

tiff(filename = "Rodent_Asym_HostQua.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$HostS, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Host richness" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$HostS[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Host2 , data=data.frame( Host=seq(-4,4, by=0.1), Host2=seq(-4,4, by=0.1)^2 ) )
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-4,4, by=0.1),er1$HostS), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-4,4, by=0.1),er1$HostS), col=col.alpha(rangi2,0.3) )

dev.off()

R2(Host2,dlist, er1$Rar90SP)


######### Temperature

Temp<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Temp ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Temp)
precis(Temp , depth=2 )

tiff(filename = "Rodent_Asym_Temp.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$Temperatur, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Annual Mean Temperature" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Temperatur[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Temp , data=data.frame( Temp=seq(-5,5, by=0.1) ))
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-5,5, by=0.1),er1$Temperatur), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-5,5, by=0.1),er1$Temperatur), col=col.alpha(rangi2,0.3) )

dev.off()


R2(Temp,dlist, er1$Rar90SP)



######### Temperature seasonality

TempS<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*TemSeason ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(TempS)
precis(TempS , depth=2 )

tiff(filename = "Rodent_Asym_TempSeasonality.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$Temperatur.1, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Temperature Seasonality" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Temperatur.1[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( TempS , data=data.frame( TemSeason=seq(-5,5, by=0.1) ))
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-5,5, by=0.1),er1$Temperatur.1), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-5,5, by=0.1),er1$Temperatur.1), col=col.alpha(rangi2,0.3) )

dev.off()

R2(TempS,dlist, er1$Rar90SP)

################Precipitation
Pre<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Prec ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Pre)
precis(Pre , depth=2 )

tiff(filename = "Rodent_Asym_Precipitation.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$Precipitat, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Annual Precipitation" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Precipitat[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Pre , data=data.frame( Prec=seq(-5,5, by=0.1) ))
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-5,5, by=0.1),er1$Precipitat), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-5,5, by=0.1),er1$Precipitat), col=col.alpha(rangi2,0.3) )

dev.off()

R2(Pre,dlist, er1$Rar90SP)



################Precipitation Seasonality
PreS<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*PercSeason ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(PreS)
precis(PreS , depth=2 )

tiff(filename = "Rodent_Asym_PrecipitationSeasonality.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$Precipitat.1, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Precipitation Seasonality" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$Precipitat.1[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( PreS , data=data.frame( PercSeason=seq(-5,5, by=0.1) ))
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-5,5, by=0.1),er1$Precipitat.1), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-5,5, by=0.1),er1$Precipitat.1), col=col.alpha(rangi2,0.3) )

dev.off()
DIC(PreS) 
R2(PreS,dlist, er1$Rar90SP)


################Area
Area<-ulam(
alist(
Pa_obs ~ dnorm( Pa_real , D_sd ),
vector[N]:Pa_real ~ dnorm(mu,sigma),
mu <- a + bL*Area ,
a ~ dnorm(0,50),
bL ~ dnorm( 0 , 20),
sigma ~ dexp(1)
) , data=dlist , chains=4 , cores=4 , log_lik=TRUE)


trankplot(Area)
precis(Area , depth=2 )

tiff(filename = "Rodent_Asym_Area.tiff",width = 1200, height = 1200)
par(mar=c(8, 7, 4, 2)+2,mgp=c(6, 2.5, 0)) 
	 
plot(er1$EcoRegionA, dlist$Pa_obs, pch=16 , col=rangi2 ,
xlab="Ecoregion Area (km2)" , ylab="Species richness-Asypmtotic est" ,cex=3,cex.lab=2.8,cex.axis=2.5)

for ( i in 1:length(er1$Rar90SP) ) {
ci <- c(er1$Rar90SP_Lo[i], er1$Rar90SP_Up[i])
x <- er1$EcoRegionA[i]
lines( c(x,x) , ci, lwd=2 )
}

p <- rethinking::link( Area , data=data.frame( Area=seq(-5,5, by=0.1) ))
p_mean <- apply( p , 2 , mean )
p_ci <- apply( p , 2 , PI , prob=0.97 )
lines( backScale(seq(-5,5, by=0.1),er1$EcoRegionA), p_mean , lwd=3 )
shade( p_ci ,  backScale(seq(-5,5, by=0.1),er1$EcoRegionA), col=col.alpha(rangi2,0.3) )

dev.off()
DIC(PreS) 
R2(Area,dlist, er1$Rar90SP)