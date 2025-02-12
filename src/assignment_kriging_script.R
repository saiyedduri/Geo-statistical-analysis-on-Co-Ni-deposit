objects()
library(sp)
library(gstat)
#clean up with 
rm(list=ls())
#inking
cf = function(x) {
  crange = range(x)
  hsv(0.7*(crange[2]-x)/(crange[2]-crange[1]),0.7,0.7)
}

##1. load data
testdata<-read.table("testData123.csv",header=T,sep=",")

sapply(testdata,class)
testdata

require("MASS")
##2. Explorative Data analysis
eqscplot(testdata$x,testdata$y,col=cf(testdata$Co),pch=19,cex=1.7,xlab="x [m]",ylab="y [m]")
legend("right",title="Co",legend=round(seq(min(testdata$Co),max(testdata$Co),length.out=10),digits=2),fill=cf(round(seq(min(testdata$Co),max(testdata$Co),length.out=10),digits=2)))

eqscplot(testdata$x,testdata$y,col=cf(testdata$Ni),pch=19,cex=1.7,xlab="x [m]",ylab="y [m]")
legend("right",title="Ni",legend=round(seq(min(testdata$Ni),max(testdata$Ni),length.out=10),digits=2),fill=cf(round(seq(min(testdata$Ni),max(testdata$Ni),length.out=10),digits=2)))

summary(testdata[,-c(1:3)])
par(mfrow=c(2,3))		# tile display
hist(testdata$Co)
qqnorm(testdata$Co)
boxplot(testdata$Co)
hist(testdata$Ni)
qqnorm(testdata$Ni)
boxplot(testdata$Ni)
par(mfrow=c(1,1))
cor(testdata[,4:5])
### convert to spatial dataset (ds) 
ds = SpatialPointsDataFrame(testdata[,c("x","y")],testdata)

##3. Empirical Variogram
?variogram
plot(variogram((Co)~1, data=ds, cloud=TRUE),main="Variogram Cloud -Co",pch=".")
plot(variogram((Ni)~1, data=ds, cloud=TRUE),main="Variogram Cloud -Ni",pch=".")

###unidirectional Variogram
vCo <- variogram((Co)~1, ds)
vCo
plot(vCo,main="experimental Semivariogram - Co")

vNi <- variogram((Ni)~1, ds)
vNi
plot(vNi,main="experimental Semivariogram - Ni")

###Check for Anisothrophy
vCo_map <- variogram((Co)~1,ds,cutoff=450,width=60,map=TRUE)
vNi_map <- variogram((Ni)~1,ds,cutoff=450,width=60,map=TRUE)
plot(vCo_map,main="Variogram Map - Co")
plot(vNi_map,main="Variogram Map - Ni")

##4. Variogrammodel

mCo_sph<-fit.variogram(vCo, vgm(18, "Sph", 1250,2)) #vgm(sill,model,range,nugget)
mCo_sph
plot(vCo,mCo_sph,main="Spherical Variogrammodell - Co")

mCo_exp<-fit.variogram(vCo, vgm(18, "Exp", 1250,2)) #vgm(sill,model,range,nugget)
mCo_exp
plot(vCo,mCo_exp,main="Exponential Variogrammodell - Co")

plot(vCo$dist,vCo$gamma,main="comparison of unidirectional variogrammmodels - Cu",xlab="distance",ylab="gamma")
mCo_exp_L<-variogramLine(mCo_exp,maxdist=1500,n=20)
lines(mCo_exp_L$dist,mCo_exp_L$gamma,col="blue")
mCo_sph_L<-variogramLine(mCo_sph,maxdist=1500,n=20)
lines(mCo_sph_L$dist,mCo_sph_L$gamma,col="green")
legend("bottomright",lty=1,col=c("blue","green"),legend=c("Exponentielles Mod.: Nug=3.077, sill=15.37, range=633.83","Sphaerisches Mod.: Nug=5.96, sill=22.36, range=1566.909"))

# Spherical variogram for Ni
mNi_sph<-fit.variogram(vNi, vgm(30, "Sph", 1400,2)) #vgm(sill,model,range,nugget)
mNi_sph
plot(vNi,mNi_sph,main="SphericalVariogrammodell - Ni")

mNi_exp<-fit.variogram(vNi, vgm(30, "Exp", 1400,2)) #vgm(sill,model,range,nugget)
mNi_exp
plot(vNi,mNi_exp,main="Exponential Variogrammodell - Co")

plot(vNi$dist,vNi$gamma,main="comparison of unidirectional variogrammmodels - Ni",xlab="distance",ylab="gamma")
mNi_exp_L<-variogramLine(mNi_exp,maxdist=1500,n=20)
lines(mNi_exp_L$dist,mNi_exp_L$gamma,col="blue")
mNi_sph_L<-variogramLine(mNi_sph,maxdist=1500,n=20)
lines(mNi_sph_L$dist,mNi_sph_L$gamma,col="green")
legend("bottomright",lty=1,col=c("blue","green"),legend=c("Exponentielles Mod.: Nug=4.32, sill=28.05, range=829.92","Sphaerisches Mod.: Nug=5.96, sill=22.36, range=1566.90"))


###directional Variogram
###Co
vCo_dir <- variogram((Co)~1, ds,alpha=seq(0,180,by=45),tol.hor=22.5,cutoff=1500,width=100)
mCo_dir<-fit.variogram(vCo_dir, vgm(18, "Exp", 1250,1,anis=c(45,0.5))) #vgm(sill,model,range,nugget)
plot(vCo_dir,mCo_dir,main="directional variogram-Co")

###Ni
vNi_dir <- variogram((Ni)~1, ds,alpha=seq(0,180,by=45),tol.hor=22.5,cutoff=1500,width=100)
mNi_dir<-fit.variogram(vNi_dir, vgm(22,"Sph",1500,0.2,anis=c(90,1))) #vgm(sill,model,range,nugget)
plot(vNi_dir,mNi_dir,main="directional variogram-Ni")

# Kriging
###Create a grid for the kriging estimation
min(testdata$x)
max(testdata$x)
x.grid<-seq(min(testdata$x)-50,max(testdata$x)+50,by=10)
x.grid
y.grid<-seq(min(testdata$x)-50,max(testdata$y)+50,by=10)
y.grid
testdata.grid<-expand.grid(x=x.grid,y=y.grid)
testdata.grid
gridded(testdata.grid) = ~x+y

# Since no strength trend on x,y variables and no global mean has been recognized, ordinary kriging has been chosen
# as the kriging method. ordinary kriging estimates the unknown point by considering the 
# neighbood means.

# Using ordinary kriging for kriging estimation
?krige
tdk_Co <- krige(formula = Co~1,ds, testdata.grid, model = mCo_dir,nmax=10)
plot(tdk_Co)
plot(tdk_Co["var1.var"],main="Kriging-variance for Co")
spplot(tdk_Co, main = "Comparision of predicted Co concentrations with variance")
summary(tdk_Co)

tdk_Ni <- krige(formula = Ni~1,ds, testdata.grid, model = mNi_dir,nmax=10)
plot(tdk_Ni)
plot(tdk_Ni["var1.var"],main="Kriging-variance for Ni")
spplot(tdk_Ni, main = "Comparision of predicted Ni concentrations with variance")
summary(tdk_Ni)

# Cross-validation through leave-one-out cross validation.
#cross validation (leave-one-out) for Co:
?krige.cv
mcv_co <- krige.cv(Co~1, ds, mCo_dir, nmax = 40)
mcv_co
sd_co <- sd(mcv_co$residual)
sd_co
summary(mcv_co)
bubble(mcv_co, "residual", main = "Co-leave-one-out cross validation")
plot(mcv_co$var1.pred,mcv_co$observed,main="Co: leave-one-out cross validation",sub="estimated vs observed values")
abline(0,1)

# Cross-validation for Ni:
?krige.cv
mcv_ni <- krige.cv(Ni~1, ds, mNi_dir, nmax = 40)
mcv_ni
sd_ni<- sd(mcv_ni$residual)
sd_ni
summary(mcv_ni)
bubble(mcv_ni, "residual", main = "Ni-leave-one-out cross validation")
plot(mcv_ni$var1.pred,mcv_ni$observed,main="Ni:leave-one-out cross validation",sub="estimated vs observed values")
abline(0,1)
