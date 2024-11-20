###############################################################
# CARBayes version 6.1.1: An R Package for Spatial            #
# Areal Unit Modelling with Conditional                       #
# Autoregressive Priors                                       #
###############################################################
# Author: Duncan Lee (University of Glasgow)                  #
###############################################################
# Example 1 - property prices in Greater Glasgow              #
###############################################################
# Source:
# https://cran.r-project.org/web/packages/CARBayes/vignettes/CARBayes.pdf
# Installing R Packages
# install.packages("GGally")
# install.packages("CARBayes")
# install.packages("CARBayesdata")
# install.packages("sf")
# Loading R Packages 
###############################################################
library(CARBayesdata)
library(CARBayes)
library(sf)
library(GGally)
library(mapview)
library(RColorBrewer)
library(spdep)
data(pricedata)
data(GGHB.IZ)
head(pricedata)
head(GGHB.IZ)
mode(GGHB.IZ)
class(GGHB.IZ)
# Property prices are positive and skewed to the right, and an initial linear regression
# model including all the covariates showed the residuals from this model were non-normal and
# skewed to the right. Therefore we model property price on the natural log scale, and the log
# transformed variable can be added to the data set using the code below.
library(dplyr)
pricedata <- pricedata %>% mutate(logprice = log(pricedata$price))
head(pricedata)
# Then the relationships between these variables can be visualised using the GGally package
# as follows
GGally::ggpairs(data = pricedata, columns = c(8, 3:7))
# The next step is to produce a map of property prices across the Greater Glasgow and Clyde 
# region, but to do this the data must be merged with the sf object GGHB.IB. This merging can 
# be done as follows:
pricedata.sf <- merge(x=GGHB.IZ, y=pricedata, by="IZ", all.x=FALSE)
head(pricedata.sf)
# The pricedata.sf object needs to have its coordinate reference system changed to longitude 
# and latitude as this is what the mapview package requires.
# The st_transform function transform or convert coordinates of simple feature
# Arguments:
# x: object of class sf, sfc or sfg
# crs: target coordinate reference system
# The World Geodetic System 1984 (WGS 84) is a 3-dimensional coordinate reference frame for 
# establishing latitude, longitude and heights for navigation, positioning and targeting
pricedata.sf <- sf::st_transform(x=pricedata.sf,crs="+proj=longlat +datum=WGS84 +no_defs")
# Then a map of price can be drawn using the following code:
mapview::mapview(pricedata.sf, zcol = "price", wd=0.5, layer.name="Price")
# The at argument can be used to define your own breakpoints.
mapview::mapview(pricedata.sf, zcol = "price", wd=0.5, layer.name="Price",at=c(50,100,150,200,250,300,350,400))
# Non-spatial modelling
form  <- logprice~crime+rooms+sales+factor(type) + driveshop
model <- lm(formula=form, data=pricedata.sf)
summary(model)
W.nb   <- spdep::poly2nb(pricedata.sf, row.names = pricedata.sf$IZ)
W.list <- spdep::nb2listw(W.nb, style="B")
# The Moran’s I test has a p-value much less than 0.05, which suggests that the 
# residuals contain substantial positive spatial autocorrelation
spdep::moran.mc(x=residuals(model), listw=W.list, nsim=1000)
# Spatial modelling with CARBayes
# First we need to create the neighbourhood matrix W
W <- nb2mat(W.nb, style="B")
# Then the model can be run with 3 parallel Markov chains on 3 processors of the same computer
# as follows, which is much computationally quicker than running the 3 chains one after the other
# S.CARleroux: Fit a spatial generalised linear mixed model to data, where the random effects have
# a Leroux conditional autoregressive prior
# Arguments:
# W is the non-negative K by K neighbourhood matrix
# n.chains is the number of MCMC chains to run when fitting the model. Defaults to 1.
# n.cores is the number of computer cores to run the MCMC chains on. Must be less than or equal to 
# n.chains. Defaults to 1.
chain <- CARBayes::S.CARleroux(formula=form, data=pricedata.sf, family="gaussian", W=W,
            burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
print(chain)
summary(chain)	
summary(chain$samples)
summary.beta <- summary(chain$samples$beta, quantiles=c(0.025, 0.975))
summary.beta
# Convergence plots of the Markov chains
plot(chain$samples$beta[ ,2:4])
# Then from this summary output a table of posterior means and 95% credible intervals can
# be computed as follows:
beta.mean <- summary.beta$statistics[ ,"Mean"]
beta.ci <- summary.beta$quantiles
beta.results <- cbind(beta.mean, beta.ci)
rownames(beta.results) <- colnames(chain$X)
round(beta.results, 5)
