#########################################################################################################
#########################################################################################################
#                                                                                                       #
# Aplicação do método Integrated Nested Laplace Approximation (INLA)                                    #
# Exemplo: Low Birth Weigth in Georgia (US)                                                             #
# Blangiardo et al (2012)                                                                               #
#                                                                                                       #
# Este arquivo acrescenta comentários e adaptações ao script introduzido por                            #
# Blangiardo et al (2012) disponível em                                                                 #
# https://inla.r-inla-download.org/r-inla.org/case-studies/Blangiardo-et-al-2012/Code/LowBirthWeight.R  #
#                                                                                                       #
# Ver também: https://www.paulamoraga.com/book-geospatial/sec-inla.html                                 #
#                                                                                                       #
# Referência para o método INLA:                                                                        #
# Rue, H. and Martino, S. and Chopin, N. (2009) Approximate Bayesian Inference for latent Gaussian      # 
# models using Integrated Nested Laplace Approximations, JRSS-series B (with discussion), vol 71,       #
# no 2, pp 319-392.                                                                                     #
#                                                                                                       #
#########################################################################################################
#########################################################################################################

# Importando os dados - Low birth weight in Georgia
urlfile="https://raw.githubusercontent.com/edsonzmartinez/spatial/main/data.final.csv"
data <- read.csv(urlfile)
# Dados introduzidos por 
# Lawson, A., 2009. Bayesian Disease Mapping. Hierarchical Modeling in Spatial Epidemiology. CRC Press.
# Contagens de crianças com baixo peso ao nascer (menos de 2500 gramas) em 159 cidades do estado da
# Georgia, de 2000 a 2010
head(data)
# As seis primeiras linhas da base de dados é mostrada a seguir, em que obs2000 a obs2010 são as contagens 
# observadas e Exp2000 são valores esperados, para o cálculo das SMR
#
#         ID obs2000 obs2001 obs2002 obs2003 obs2004 obs2005 obs2006 obs2007 obs2008 obs2009 obs2010   Exp2000
# 1  Appling      20      24      25      31      24      40      29      35      26      25      35 25.815340
# 2 Atkinson      10      12      11      15      15       9      14       8      12      12       5 14.114862
# 3    Bacon      16       8       9      19      18      19      20      18      11      21      16 14.114862
# 4    Baker       4       8       5       1       2       1       4       7       3       3       2  3.342994
# 5  Baldwin      56      64      44      58      63      56      70      65      53      46      54 51.166375
# 6    Banks      22      17      12      17      19      23      22      17      14      23      24 19.222214
#
#     Exp2001  Exp2002  Exp2003   Exp2004   Exp2005   Exp2006   Exp2007   Exp2008   Exp2009   Exp2010
# 1 25.072453 22.75093 25.16531 24.886731 26.651088 26.001062 27.208254 25.072453 24.515287 24.329565
# 2 13.093392 15.69350 14.20772 12.907670 15.693498 14.207723 13.557697 15.229193 11.886200 12.536226
# 3 13.929140 13.18625 14.57917 15.136332 15.693498 15.879220 17.179273 14.393445 15.136332 13.929140
# 4  4.457325  3.99302  2.97155  2.785828  2.507245  3.621576  3.064411  2.785828  3.064411  3.157272
# 5 48.659130 49.12343 50.98065 52.466428 53.023594 53.766482 53.487899 51.816402 50.980653 45.780441
# 6 20.429406 19.96510 18.38647 22.565207 19.779379 22.658068 18.850770 19.407935 19.593657 15.972081
#
dim(data)
# [1] 159  23
#############################################################################################################

# Pacotes
library(sp)
library(spdep)
library(lattice)
library(sf)
library(mapview)
library(INLA)

## default pc-priors
hyper <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

## Uma lista das distribuições a priori é obtida de:
inla.list.models("prior")
## Entre as possibilidades, temos:
# pc.prec: PC prior for log(precision)
## Para obter mais informações sobre uma distribuição específica:
inla.doc("pc.prec") 

# Mapa
georgia <- sf::st_read("/vsizip//vsicurl/https://raw.githubusercontent.com/edsonzmartinez/spatial/main/co13_d00.zip")
# Obsrvação de Blangiardo et al:
## Need to drop extra polygons 98 (Macon county polygon), 100, 105 (Taylor county polygons) + 137 (Lee county polygon)
## These are very small and always adjacent to "main" polygon, so we can base neighborhood structure on "main" polygon:
rmIdx <- c(98, 100, 105, 137)
georgia <- georgia[-rmIdx,]
# Visualizando o mapa
mapview::mapview(shp)

#######################################################################
#      Criando a matriz de adjacências a partir de um shapefile       #           
#######################################################################
zzz <- poly2nb(georgia)
summary(zzz)
# Se usarmos:
head(zzz)
# Teremos a lista dos polígonos que são vizinhos a cada um dos 
# polígonos do mapa
# A função nb2INLA() converte estas vizinhanças para o formato 
# usado pelo pacote INLA
nb2INLA("Georgia.graph", zzz)
# O código original de Blangiardo et al (2012) utiliza esta linha
# para criar um objeto Georgia.adj que será posteriormente usado
# nos modelos INLA:
# Georgia.adj <- paste(getwd(),"/Georgia.graph",sep="")
# Alternativamente, podemos usar:
Georgia.adj <- inla.read.graph(filename = "Georgia.graph")

# A base de dados 'data' é reordenada de acordo com a sequência dos
# polígonos no shapefile
order <- match(georgia$NAME,data[,1])
data  <- data[order,]

#--Transform the data to be in the right format for INLA--#
low.vector <- as.vector(as.matrix(data[,2:12]))  #by column
E.vector   <- as.vector(as.matrix(data[,13:23])) #by column
year <- numeric(0)
for(i in 1:11){
  year <- append(year,rep(i,dim(data)[1]))
}
county <- as.factor(rep(data[,1],11))

data<- data.frame(y= low.vector, E= E.vector, ID.area=as.numeric(county), ID.area1=as.numeric(county), year=year,
                  ID.year = year, ID.year1=year, ID.area.year = seq(1,length(county)))
dim(data)
# [1] 1749    8
head(data)
#    y          E ID.area ID.area1 year ID.year ID.year1 ID.area.year
# 1 13  16.436386     119      119    1       1        1            1
# 2  2   8.264623     139      139    1       1        1            2
# 3 19  21.172293      55       55    1       1        1            3
# 4 37  53.023594     105      105    1       1        1            4
# 5 97 162.042332     155      155    1       1        1            5
# 6 65  65.745542      23       23    1       1        1            6
tail(data)
#        y         E ID.area ID.area1 year ID.year ID.year1 ID.area.year
# 1744  29  19.50080      14       14   11      11       11         1744
# 1745  63  55.34512     136      136   11      11       11         1745
# 1746  17  12.90767      24       24   11      11       11         1746
# 1747  15  10.40042     125      125   11      11       11         1747
# 1748 216 161.29944      92       92   11      11       11         1748
# 1749   6   5.94310      50       50   11      11       11         1749

# O motivo de criar duas variáveis idênticas, ID.area e ID.area1, é explicado
# por Blangiardo et al.:
# Note that each function f(.) can only be assigned to one covariate
# in R-INLA, so in this case we need to create a new variable ID.area1 which
# is a duplicate of ID.area.


#######################################################################
#                     Especificando o modelo                          #
#######################################################################
# formula <- y ~ 1 + x1 + x2 + f(z1, ...)
# y, x1, x2 e z1 são as variáveis
# f(.) especifica a estrutura de uma função dada por
# f(z1, model = "...", ...)
# o default é model="iid"
# uma lista de possibilidades pode ser visualizada usando
names(inla.models()$latent)
# Exemplo:
formula <- y ~ 1 + f(ID, model="bym", graph=mapa.adj)
# ID identifica a área no shapefile
# bym é o Besag-York-Mollie (BYM) (Besag et al., 1991),
# que assume uma intrinsic conditional autoregressive structure (iCAR)
# em graph=mapa.adj, mapa.adj é o conjunto de adjacências
# No presente exemplo, usaremos graph=Georgia.adj
# Outras especificações em model= são:
# rw1: Random walk of order 1                  
# rw2: Random walk of order 2       
# hyper é usado para definir as hiperprioris:
prior.prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
formula    <- y ~ f(ID, model = "iid", hyper = prior.prec)
# A interpretação de 1 e 0.01 em param = c(1, 0.01) pode ser
# encontrada em:
inla.doc("pc.prec")
# O argumento scale.model = TRUE é usado para tornar o parâmetro de
# precisão de diferentes modelos com prioris CAR comparáveis:
# formula <- y ~ 1 + f(ID, ..., scale.model = TRUE)
#######################################################################
# Modelo baseado na distribuição de Poisson:
m1 <- inla(formula, family="poisson", data=data, E=E, 
                    control.predictor=list(compute=TRUE), 
					control.compute=list(dic=TRUE,cpo=TRUE))
# E=E é usado para especificar um componente conhecido na média da
# distribuição de Poisson, definida por E exp(eta), em que eta é
# o preditor linear. Se não especificado, E é um vetor [1,1,1,...,1]'
# control.predictor=list(compute=TRUE) calcula as posterioris das
# predições
# control.compute=list(dic=TRUE,cpo=TRUE) 
# Observar que o default é dic=FALSE e cpo=FALSE, temos que incluir
# estes argumentos para obtermos estas medidas
# Para outras opções, usar
?control.compute
#######################################################################


#######################################################################
# Parametric model alpha + csi_i + (delta_i + beta)*year              #
#######################################################################
# alpha é o intercepto (efeito fixo)
# beta é um efeito geral do tempo (efeito fixo)
# delta_i é um efeito linear do tempo em cada área, delta_i ~ N(0,taud)
# ID.area e ID.area1 são variáveis idênticas e correspontem à
# identificação de cada área no shapefile

formula.ST1 <- y ~ 1 + f(ID.area, model="bym2", graph=Georgia.adj) +
                       f(ID.area1,year,model="rw1", hyper = hyper, scale.model = TRUE) + 
					   (year-mean(year))
model.inla.ST1 <- inla(formula.ST1, family="poisson", data=data, E=E, 
                       control.predictor=list(compute=TRUE), 
					   control.compute=list(dic=TRUE,cpo=TRUE,return.marginals.predictor = TRUE))
model.inla.ST1
summary(model.inla.ST1)

# Gráficos dos resultados:
plot(model.inla.ST1)
# Gráficos dos resultados, comparando com as distribuições a priori:
plot(model.inla.ST1,plot.prior = TRUE)
#
# Deviance Information Criterion (DIC):
model.inla.ST1$dic$dic
# Resumos para os efeitos fixos:
model.inla.ST1$summary.fixed
# Observação:
# The column kld represents the symmetric Kullback-Leibler divergence 
# (Kullback and Leibler 1951) that describes the difference between the 
# Gaussian and the simplified or full Laplace approximations for each posterior
# Resumos para os efeitos aleatórios:
model.inla.ST1$summary.random
# Resumos para os hiperparâmetros:
model.inla.ST1$summary.hyperpar
# Preditores lineares:
model.inla.ST1$summary.linear.predictor
# Gráficos das marginais (spline smoothing):
par(mfrow = c(1, 2))
plot(inla.smarginal(model.inla.ST1$marginals.fixed[[1]]),type="l",xlab="Intercept")
plot(inla.smarginal(model.inla.ST1$marginals.fixed[[2]]),type="l",xlab="year")
# Valores ajustados pelo modelo:
head(model.inla.ST1$summary.fitted.values)
# Observados versus preditos:
a <- model.inla.ST1$summary.fitted.values$mean
b <- data$y/data$E
plot(b,a,las=1,xlim=c(0,3),ylim=c(0,3))
lines(c(0,3),c(0,3),col="blue")


#Non Parametric model alpha + csii + gammaj + phij #No space time interaction yet!
#csi_i and are modelled through BYM2
#gammaj are modelled as RW1
#phij are modelled as exchangeable
formula.ST2    <- y ~ 1 + f(ID.area,model="bym2",graph=Georgia.adj) +
                  f(ID.year,model="rw1", hyper = hyper, scale.model = TRUE) + 
				  f(ID.year1,model="iid", hyper = hyper)
model.inla.ST2 <- inla(formula.ST2, family="poisson", data=data, E=E, 
                  control.predictor=list(compute=TRUE), 
				  control.compute=list(dic=TRUE,cpo=TRUE))

#Non Parametric model alpha + csii + gammaj + phij + deltaij
#csii are modelled through BYM2
#gammaj are modelled as RW1
#phij are modelled as exchangeable
#Interaction (deltaij) is modelled as exchangeable
formula.ST3<- y ~ 1 + f(ID.area,model="bym2",graph=Georgia.adj) +
  f(ID.year,model="rw1", hyper = hyper, scale.model = TRUE) + f(ID.year1,model="iid", hyper = hyper) + f(ID.area.year,model="iid", hyper = hyper)

#To obtain the marginal of phij + gammaj we need to create the corresponding linear combinations and include these in the model
lcs = inla.make.lincombs(ID.year = diag(11),  ID.year1 = diag(11))

model.inla.ST3 <- inla(formula.ST3,family="poisson",data=data,E=E,
                       control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE),
                       lincomb=lcs)

#Put the temporal effect  (gammaj+phij) on the natural scale
temporal<-lapply(model.inla.ST3$marginals.lincomb.derived, function(X){
  marg <- inla.tmarginal(function(x) exp(x), X)
  inla.emarginal(mean, marg)
})
# Observados versus preditos:
a <- model.inla.ST3$summary.fitted.values$mean
b <- data$y/data$E
plot(b,a,las=1,xlim=c(0,3),ylim=c(0,3))
lines(c(0,3),c(0,3),col="blue")

#############################################################
# Computethe DIC as a tool for model choice
model.inla.ST1$dic$dic
model.inla.ST2$dic$dic
model.inla.ST3$dic$dic

# DIC components: Effective number of parameter (pd)
model.inla.ST1$dic$p.eff
model.inla.ST2$dic$p.eff
model.inla.ST3$dic$p.eff
# DIC components: mean.deviance
model.inla.ST1$dic$mean.deviance
model.inla.ST2$dic$mean.deviance
model.inla.ST3$dic$mean.deviance

######################################################
# The last model (with interaction) shows the best fit. 
# Obtain zeta_i exponentiating csi_i
m <- model.inla.ST3$marginals.random[[1]][seq_len(nrow(georgia))]
zeta <- unlist(lapply(m,function(x)inla.emarginal(exp,x)))

#Probability that theta>1
a=0
inlaprob <- lapply(model.inla.ST3$marginals.random[[1]][seq_len(nrow(georgia))], function(X){
  1-inla.pmarginal(a, X)
})

Spatial.results<- data.frame(NAME=georgia$NAME,zeta=unlist(zeta), pp=unlist(inlaprob))

# Maps
# Create classes of SMRs
zeta.cutoff<- c(0.6, 0.9, 1.0, 1.1,1.8)
pp.cutoff <- c(0,0.2,0.8,1)
zeta=cut(Spatial.results$zeta,breaks=zeta.cutoff,include.lowest=TRUE)
pp=cut(Spatial.results$pp,breaks=pp.cutoff,include.lowest=TRUE)

maps.factors <- data.frame(NAME=georgia$NAME, zeta=zeta,pp=pp)
georgia <- cbind(georgia, maps.factors)

georgia_sp <- sf::as_Spatial(georgia)

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol= "zeta", col.regions=gray(3.5:0.5/4),main="",par.settings=list(fontsize=list(text=17)))

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol= "pp", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))

#Plot the National temporal trend
plot(seq(1,11),seq(0.8,1.2,length=11),type="n",xlab="year",ylab=expression(exp(gamma[t]+phi[t])))
lines(unlist(temporal))
abline(h=1,lty=2)

########################
# Space-Time Interaction
delta <- data.frame(delta=model.inla.ST3$summary.random$ID.area.year[,2],year=data$ID.year,ID.area=data$ID.area)
delta.matrix <- matrix(delta[,1], nrow(georgia),11,byrow=FALSE)
rownames(delta.matrix)<- delta[seq_len(nrow(georgia)),3]

# Space time probability>1
a=0
inlaprob.delta<-lapply(model.inla.ST3$marginals.random[[4]], function(X){
  1-inla.pmarginal(a, X)
})
pp.delta<-unlist(inlaprob.delta)

pp.cutoff.interaction <- c(0,0.2,0.8,1)
pp.delta.matrix <- matrix(pp.delta, nrow(georgia),11,byrow=FALSE)
pp.delta.factor <- data.frame(NAME=georgia$NAME)
for(i in 1:11){
  pp.delta.factor.temp <- cut(pp.delta.matrix[,i],breaks=pp.cutoff.interaction,include.lowest=TRUE)
  pp.delta.factor <- cbind(pp.delta.factor,pp.delta.factor.temp)
}
colnames(pp.delta.factor)<- c("NAME",seq(2000,2010))

#Maps
georgia <- cbind(georgia, pp.delta.factor)
georgia_sp <- sf::as_Spatial(georgia)

par(mfrow = c(2, 2))

trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol="X2001", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol="X2004", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol="X2007", col.regions=gray(2.5:0.5/3),main="",par.settings=list(fontsize=list(text=17)))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=georgia_sp, zcol="X2010", col.regions=gray(2:0/2),main="",par.settings=list(fontsize=list(text=17)))

#############################################################
