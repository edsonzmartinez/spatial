####################################################
# Códigos R para o curso                           #
# VISUALIZAÇÃO DE MAPAS E INFORMAÇÕES              #
# ESPACIAIS EM SAÚDE UTILIZANDO R PARA             #
# VIGILÂNCIA EPIDEMIOLÓGICA                        #
####################################################
# Ministrado em novembro de 2024 no EPI            #
####################################################

library(sf)
library(spdep)
library(mapview)
library(CARBayes)
library(glmtoolbox)
library(sandwich)

shp <- sf::st_read("/vsizip//vsicurl/https://raw.githubusercontent.com/edsonzmartinez/spatial/main/BR_UF_2021.zip")

mapview::mapview(shp)
class(shp)
mode(shp)
names(shp)
# Rondônia
mapaRO <- data.frame(shp$geometry[[1]][1])
head(mapaRO,25)
plot(mapaRO$X1,mapaRO$X2,pch=16,cex=0.4)
plot(mapaRO$X1,mapaRO$X2,type="l",xlab="Longitude",ylab="Latitude")
points(-63.903889, -8.761944,pch=19)
points(-61.455960,-11.437710,pch=19)
points(-62.717149,-11.696318,pch=19)
text(-63.903889, -8.761944,"Porto\nVelho",pos=4)
text(-61.455960,-11.437710,"Cacoal",pos=1)
text(-62.717149,-11.696318,"São Miguel\ndo Guaporé",pos=1)


urlfile="https://raw.githubusercontent.com/edsonzmartinez/basesdedados/main/IDH_UF_Brasil.csv"
w  <- read.csv2(urlfile)
w2 <- data.frame(SIGLA=w$UF, IDH2021=w$IDH2021, E2021=w$E2021, R2021=w$R2021, L2021=w$L2021)
shp.sf <- merge(x=shp, y=w2, by="SIGLA", all.x=FALSE)
mapview::mapview(shp.sf, zcol = "IDH2021",layer.name="IDH 2021")

# AIDS
urlfile="https://raw.githubusercontent.com/edsonzmartinez/basesdedados/main/AIDS2020_2023.csv"
aids <- read.csv2(urlfile)
aids$incd <- 10000*aids$AIDS2022/aids$POP2022
aids2 <- data.frame(SIGLA=aids$SIGLA, incd=aids$incd, AIDS2022=aids$AIDS2022, pop2022=aids$POP2022)
shp.sf <- merge(x=shp.sf, y=aids2, by="SIGLA", all.x=FALSE)
mapview::mapview(shp.sf, zcol = "incd",layer.name="Incidência AIDS 2022")
plot(shp.sf$IDH2021,shp.sf$incd)

# Vizinhos
cr <- rep("#e5f5f9",27)
cr[3] <- "#de2d26"
cr[c(1,2,4,5,25)] <- "#fc9272"
plot(shp[1],col=cr,main="")

# Criando uma matriz de vizinhanças
W.nb   <- spdep::poly2nb(shp.sf, row.names = shp.sf$SIGLA)
W.list <- spdep::nb2listw(W.nb, style="B")
W      <- spdep::nb2mat(W.nb, style="B")
print(matrix(W,ncol=27,dimnames=list(shp.sf$SIGLA,shp.sf$SIGLA)))

plot(shp.sf$IDH2021,shp.sf$incd,ylim=c(0.5,4),pch=19)
text(shp.sf$IDH2021,shp.sf$incd,shp.sf$SIGLA,pos=1,cex=0.8,col="red")

# Modelo de regressão
model <- lm(incd ~ IDH2021, data=shp.sf)
summary(model)
abline(model)
# Diagnósticos de resíduos
par(mar = c(4, 4, 2, 2), mfrow = c(2, 2))
plot(model)
glmtoolbox::envelope(model, rep=5000, col="red", type="internal")
# Moran
spdep::moran.mc(x=residuals(model), listw=W.list, nsim=1000)

# Modelo de Poisson não espacial
poisson.model <- glm(AIDS2022 ~ IDH2021 + offset(log(pop2022)), data=shp.sf, family = poisson(link = "log"))
summary(poisson.model)
# Erros padrão robustos (Cameron and Trivedi, 2009)
# Cameron, A. C. and Trivedi, P. K. 2009. Microeconometrics Using Stata. College Station, TX: Stata Press.
# Cameron, A. C. and Trivedi, P. K. 1998. Regression Analysis of Count Data. New York: Cambridge Press.
cov.m1 <- vcovHC(poisson.model, type="HC0")
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(Estimate= coef(poisson.model), "Robust SE" = std.err, "Pr(>|z|)" = 2 * pnorm(abs(coef(poisson.model)/std.err), lower.tail=FALSE),
LL = coef(poisson.model) - 1.96 * std.err,
UL = coef(poisson.model) + 1.96 * std.err)
r.est
# Pearson's goodness-of-fit
Pearson <- sum((shp.sf$AIDS2022 - poisson.model$fitted.values)^2 / poisson.model$fitted.values)
message("Pearson's goodness-of-fit = ",round(Pearson,3),", p-value ",1 - pchisq(Pearson, df = poisson.model$df.residual))
message("Estimated dispersion parameter = ",round(Pearson / poisson.model$df.residual,4))
# Gráfico
a <- as.numeric(poisson.model$coefficients[1])
b <- as.numeric(poisson.model$coefficients[2])
x <- seq(0.5,1,0.001)
curve <- function(x,a,b,pop) pop*exp(a+b*x)
plot(shp.sf$IDH2021,shp.sf$incd,ylim=c(0.5,4),pch=19)
text(shp.sf$IDH2021,shp.sf$incd,shp.sf$SIGLA,pos=1,cex=0.8,col="blue")
points(x,curve(x,a,b,10000),type="l",col="red")
spdep::moran.mc(x=residuals(poisson.model), listw=W.list, nsim=100000)

# Modelo Quasi Poisson não espacial
qpoisson.model <- glm(AIDS2022 ~ IDH2021 + offset(log(pop2022)), data=shp.sf, family = quasipoisson(link = "log"))
summary(qpoisson.model)
spdep::moran.mc(x=residuals(qpoisson.model), listw=W.list, nsim=100000)

# Modelo de Poisson espacial
form <- AIDS2022 ~ IDH2021 + offset(log(pop2022))
chain  <- CARBayes::S.CARleroux(formula=form, data=shp.sf, family="poisson", W=W,
   burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
summary.beta <- summary(chain$samples$beta, quantiles=c(0.025, 0.975))
a <- summary.beta$statistics[1,1]
b <- summary.beta$statistics[2,1]
plot(shp.sf$IDH2021,shp.sf$incd,ylim=c(0.5,4),pch=19)
text(shp.sf$IDH2021,shp.sf$incd,shp.sf$SIGLA,pos=1,cex=0.8,col="blue")
points(x,curve(x,a,b,10000),type="l",col="red")
# Convergência
plot(chain$samples$beta)
spdep::moran.mc(x=residuals(chain), listw=W.list, nsim=100000)

# Modelo de Poisson espacial
form <- AIDS2022 ~ IDH2021 + offset(log(pop2022))
chain  <- CARBayes::S.CARbym(formula=form, data=shp.sf, family="poisson", W=W,
   burnin=100000, n.sample=300000, thin=100, n.chains=3, n.cores=3)
summary.beta <- summary(chain$samples$beta, quantiles=c(0.025, 0.975))
summary.beta <- summary(chain$samples$beta, quantiles=c(0.025, 0.975))
a <- summary.beta$statistics[1,1]
b <- summary.beta$statistics[2,1]
plot(shp.sf$IDH2021,shp.sf$incd,ylim=c(0.5,4),pch=19)
text(shp.sf$IDH2021,shp.sf$incd,shp.sf$SIGLA,pos=1,cex=0.8,col="blue")
points(x,curve(x,a,b,10000),type="l",col="red")
# Convergência
plot(chain$samples$beta)
spdep::moran.mc(x=residuals(chain), listw=W.list, nsim=100000)
