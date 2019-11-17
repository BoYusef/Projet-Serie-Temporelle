##################################################################################
# TITRE      : Pratiques sur des données Temporelle                              #
# ANNEE      : 2018-2019                                                         #
# DATE       : Mai - Juin                                                        #
# ECOLE      : ENSAE - SENEGAL                                                   #
# CLASSE     : ITS - 3                                                           #
# AUTEURS    : Youssouf BANCE                                                    #                                                    #
# PROFESSEUR : M. Souleymane FOFANA, Chef de filière ITS                         #                      *
##################################################################################


#================================================================================#
#                 PARTIE : Importation, prélimaire et packages                   #
#================================================================================#


######################################
##### LES PACKAGES DU PROJETS

library(ggfortify)
# Définition de l'espace de travail
setwd(choose.dir())

# importation de la base dInflation

bInflation<-read.csv2("tsRAPPORT.csv",header=TRUE,sep=";",dec=",")


############### ACP sur variable temporelle ###################################

###### Imoportation de la base

bTemporelle<-read.csv2("ACPtemporelle.csv",header=TRUE,sep=";",dec=",",row.names = 1)
  

install.packages("corrplot")
library(corrplot)

coor=cor(bTemporelle)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white","cyan", "#007FFF", "blue", "#00007F"))
corrplot(coor,is.corr=TRUE, method = "circle")  #

corrplot(coor, method = "circle")
corrplot.mixed(coor,lower.col = "black", number.cex = .7)

# Plot avec ggcorrplot
install.packages("ggcorrplot")

library(ggcorrplot)
ggcorrplot(coor, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Corrélation entre IHPC et les variables susceptibles de l'expliquée", 
           ggtheme=theme_bw)

## Mise en oeuvre de l'ACP
library(FactoMineR)
library(factoextra)


res<-PCA(bTemporelle)

summary(res)

## description des dimmension (dismdesc: de facteur miner permet de voir les variables plus significatif sur un axe)
res.desc<-dimdesc(res,axes=1:2,proba=0.05)
res.desc

# descripton sur l'axe 1
res.desc$Dim.1

# description sur l'axe 2
res.desc$Dim.2

# Contribution à la formation des axes
var<-get_pca_var(res)
corrplot(var$contrib,is.corr = FALSE)

# Répresentation du cercle de correlation selon la valeur du  cos2 des variables

fviz_pca_var(res,col.var = "cos2",col.quanti.sup = "black",gradient.cols=c("white","blue","red"),repel = TRUE)  #avoid text overlapping

# Répresentation des individus 
fviz_pca_ind(res,col.ind = "cos2",gradient.cols=c("yellow","blue","red"),repel=TRUE)

# Répresentation simultannée des années et des variables
fviz_pca_biplot(res,col.var = "red",col.ind = "blue",repel=TRUE)

#==============================================================================#
#                 PARTIE : Serie temporelle descriptive                        #
#==============================================================================#


########################### Méthode descriptive  #############################

#####mise en forme de la serie temporelle

IHPCciv <- ts(bInflation[,"IPC._GLOBAL"],start=2009,freq=12)
IHPCciv <- window(IHPCciv,end=c(2018,12))
IHPCciv

attach(bInflation)

######description de la série IHPCciv
# dimension
length(IHPCciv)


# résumé
summary(IHPCciv)
sd(IHPCciv)

# Graphe de la série
plot.ts(IHPCciv,xlab="Période",ylab="IHPC",col="red",main="Evolution de l'IHPC de la Côte d'ivoire")
monthplot(IHPCciv,ylab="IHPC",main="",cex.main=8)

# Boxplot
boxplot(bInflation[,3]~bInflation[,2],col="lightblue",pch=20,cex=0.5)

# Nature additive ou multiplicateur de la série ; Shéma de décomposition
### Avec l'allure de la série, elle semble être multiplicateur

### Courbe de profils (contruction soi même)

colors()
c<-c("red","blue","yellow","cyan")
c<-matrix(c,nrow=1,ncol=4,byrow=F)
c[1,2]

cp<-ts(IHPCciv,start=c(1,1),end=c(1,12)) 
v=as.numeric(cp)  ##creation de vecteur de donnée
plot(c(1:12),v,"l",xlim=c(1,12),ylim=c(95,120),col=c[1,1])## c=matrice contenant des noms de couleurs

# les autres droites de profils
X<-matrix(IHPCciv,nrow=10,ncol=12,byrow=T)   ###Creation de matrice avec la serie
x=nrow(X)
for ( i in 2:x) { 
  cp<-ts(X[i,],start=c(1,1),end=c(1,12))
  v=as.numeric(cp)
  lines(v)
}

ggseasonplot(IHPCciv, year.labels=TRUE, year.labels.left=TRUE) +
  ylab("IHPC") +
  ggtitle("LES DROITES DE PROFILS : IHPC COTE D'IVOIRE")

ggseasonplot(IHPCciv, polar=TRUE) + ylab("IHPC")

## On voit bien que les profils ne sont pas parallèle donc modèle multiplicative

### On peut passer par la significativité d'une regression sur le temps
### donc P<0.05 donc le modèle est multiplicateur

###########décomposition de la série temporelle

####### Methode de decomposition 1 : Decomposition classique
ihpcComp <- decompose(IHPCciv,type="multiplicative")
autoplot(ihpcComp)

####tendance seulement

tendance=ihpcComp$trend
tendance
plot(tendance,col="yellow")

####random seulement


irregulier=ihpcComp$random
plot(irregulier,col="green")

####Saisonnalité seulement

saison=ihpcComp$seasonal
plot(saison,col="red")

#####Estimation des coefficients du saisonnier
coef=IHPCciv/tendance
s=tapply(coef,cycle(coef),mean,na.rm=T)
s
sd=s-mean(s)


#### Vérification que la somme des coefficient saisonier est null
(sum(sd=FALSE)==0)

###### Serie Cvs( corrige de ses variations saisonnier)
IHPCI_CVS <- IHPCciv/saison
plot(IHPCI_CVS)
ts.plot(IHPCI_CVS,IHPCciv)

########## Décomposition par moyenne mobile

ma(IHPC,13)

ma13 <- ma(IHPCciv, order=12, centre=TRUE)

Saison=IHPCciv/ma4

autoplot(IHPCciv, series="IHPC") +
  autolayer(ma(IHPCciv,13), series="Tendance 13-MA") +
  xlab("Période") + ylab("IHPC") +
  ggtitle("Indice Harmonisé des Prix à la Consommation") +
  scale_colour_manual(values=c("IHPC"="grey50","Tendance 13-MA"="red"),
                      breaks=c("IHPC","Tendance 13-MA"))

########## Décomposition par l'algorithme x11()
install.packages("seasonal")
library(seasonal)

decomp.x11 = seas(IHPCciv2018,transform.function="none",x11="")
summary(decomp.x11)

autoplot(decomp.x11$est$,facet=T) + xlab("Year") + ggtitle("X11 Decomposition de l'IHPC")
bas=data.frame(decomp.x11$trend)
# tendance
Tt = trendcycle(decomp.x11)
 
# Saisonnalité
St= seasonal(decomp.x11)

# le bruit
Rt = remainder(decomp.x11)

# Saisonnalité ajusté
decomp.seasadj = seasadj(decomp.x11)


Composante = cbind(Tt,St,Rt,decomp.seasadj)
autoplot(Composante,facet=T) + xlab("Période") + ggtitle("IHPC : Décomposition par l'algorithme X11()")

#### Calibrage des données
# Série IHPCciv1
sum=0
j=0
i=0
for(i in 0:9){
  if(i!=2){
    j=j+1
  sum=  IHPCciv[4+i*12]+sum
}}

round(sum/j,2)

IHPCciv1=IHPCciv
IHPCciv1[28]=round(sum/j,2)

# Série IHPCciv1
sum=0
j=0
i=0
for(i in 1:12){
  if(i!=4){
    j=j+1
    sum=  IHPCciv[24+i]+sum
  }}

round(sum/j,2)

IHPCciv2=IHPCciv
IHPCciv2[28]=round(sum/j,2)

# Graphique 
par(mfrow=c(2,1))
plot.ts(IHPCciv1)
plot.ts(IHPCciv2)

########################### PREVISION AVEC LES METHODES CLASSIQUES #############

# Sugestion pour les coefficients 
sugCOEF=hw=HoltWinters(log(IHPCciv),beta=FALSE,gamma=FALSE)
sugCOEF

# Lissage exponentiel double (LED)
# Sur la série IHPCciv
IHPCciv.LED<- HoltWinters(IHPCciv, gamma=0.9)
p <- predict(IHPCciv.LED, n.ahead=4)
round(p,2)

# Sur la série IHPCciv1
IHPCciv.LED1<- HoltWinters(IHPCciv1, gamma=0.9)
p <- predict(IHPCciv.LED1, n.ahead=4)
round(p,2)

# Sur la série IHPCciv12
IHPCciv.LED2<- HoltWinters(IHPCciv2, gamma=0.9)
p <- predict(IHPCciv.LED2, n.ahead=4)
round(p,2)

par(mfrow=c(3,1))
plot(IHPCciv.LED)
plot(IHPCciv.LED1)
plot(IHPCciv.LED2)



###Holt winters multiplicative
ihpclisse1 <- HoltWinters(IHPCciv, seasonal="multiplicative",alpha=0.7,beta=0.7,gamma=0.8)

ihpclisse2 <- HoltWinters(IHPCciv1, seasonal="multiplicative",alpha=0.7,beta=0.7,gamma=0.8)

ihpclisse3 <- HoltWinters(IHPCciv2, seasonal="multiplicative",alpha=0.7,beta=0.7,gamma=0.8)

h <-4
p1 <- predict(ihpclisse1, n.ahead=h)
p1

p <- predict(ihpclisse2, n.ahead=h)
p

p <- predict(ihpclisse3, n.ahead=h)
p


#########représentation des prévision
x_p <- ts(c(IHPCciv,p), start=2009, frequency=12)
######################créer une chronologie obtenue en ralongeant la serie de donnée
plot(x_p, type="o",col="green")
lines(p, col="purple")


#==============================================================================#
#        PARTIE : Modélisation univariée de l'IHPC                             #                                         #
#==============================================================================#


###Importation de la base
bInflation<-read.csv2("tsRAPPORT.csv",header=TRUE,sep=";",dec=",")

###CONFIGURATION DE L'INDICE TEMPORELLE'

IHPCciv<-ts(bInflation[,3], start=c(2009,1), frequency=12)

IHPCciv2018<-window(IHPCciv,end=c(2018,12))
IHPCciv2018

### Répresentation de la série, PACF, ACF

# Série
autoplot(IHPCciv2018) +
  xlab("Période") + ylab("IHPC-CIV") +
  ggtitle("l'Indice Harmonisé des prix à la consommation de la Côte d'ivoire")

install.packages("fpp2")
library(fpp2)
# ACF : fonction d'autocorrelation
ggAcf(IHPCciv2018, lag.max=50, plot=T)

# PACF : fonction d'autocorrelation partielle
ggPacf(IHPCciv2018, plot=T)

ban=ggPacf(IHPCciv2018, plot=T)

# Les trois d'un coup
ggtsdisplay(IHPCciv2018)

# Densité spectrale

periodogram<-function(z){
  nfreq <- n%/%2 + 1 
  freq<-seq(1:nfreq)/n
  
  periodogram<-((Mod(fft(z))^2/(2 * pi * n))/var(z))[1:nfreq]
  plot(freq,log(periodogram),type="l")
  return (freq)
}

x11()
par(mfrow=c(2,2))
plot(IHPCciv2018)
acf(IHPCciv2018)
pacf(IHPCciv2018)
frequence=periodogram(IHPCciv2018)
dev.off()

x11()
par(mfrow=c(2,1))
frequence=periodogram(IHPCciv2018)
spectrum(IHPCciv2018)

## densité spectral par estimation multitaper
install.packages("multitaper")
library(multitaper)
x11()
mt.spec <- spec.mtm(IHPCciv2018, nw = 16, k = 2 * 16 - 1, jackknife = TRUE, dtUnits = "month")

### Stationnarisation

##Centrer la série
IHPCciv2018c<-IHPCciv2018-mean(IHPCciv2018)

autoplot(IHPCciv2018c) +
  xlab("Période") + ylab("IHPC-CIV") +
  ggtitle("l'Indice Harmonisé des prix à la consommation de la Côte d'ivoire")


# "détendancialisé"
IHPCciv2018d<-diff(IHPCciv2018c, lag=1)
ts.plot(IHPCciv2018d)

ggtsdisplay(IHPCciv2018d)


##Tests de Stationnarité 
adf.test(IHPCciv2018d, k=0) # DF Simple
adf.test(IHPCciv2018d, k=1) # DF Augmanté
pp.test(IHPCciv2018d)       # Phillips-Perron
kpss.test(IHPCciv2018d)     # KPSS  

## Bien vrai que la série soit stationnaire, on observe une repercution
## d'un choc sur la serie, notamment en 2011
 

######## SPECIFICATION ou IDENTIFICATION DU MODELE ########

# Identification de l'orde de la partie MA(q)
ggAcf(IHPCciv2018d, lag.max=50, plot=T)

# q=2 ou q=2
# Identification de l'orde de la partie AR(p)
ggPacf(IHPCciv2018d, plot=T)



#ARIMA(1,1,0) # a cause de la nature du PACF et l'ACF de la série brute
#ARIMA(1,1,2)
#ARIMA(1,1,1)
#ARIMA(2,1,1)
#ARIMA(2,1,2)
#ARIMA(6,1,6)

######## ESTIMATION DES PARAMETRES ########
### Cas d'un ARMA par MLE

IHPCmodele1<-arima(IHPCciv2018c, order = c(1,1,0))
IHPCmodele2<-arima(IHPCciv2018c, order = c(1,1,1))
IHPCmodele3<-arima(IHPCciv2018c, order = c(1,1,2))
IHPCmodele4<-arima(IHPCciv2018c, order = c(2,1,1))
IHPCmodele5<-arima(IHPCciv2018c, order = c(2,1,2))
IHPCmodele6<-arima(IHPCciv2018c, order = c(6,1,6))

### Résumé des estimations 
IHPCmodele1
IHPCmodele2
IHPCmodele3
IHPCmodele4
IHPCmodele5
IHPCmodele6

### Variance des coefficients des modèles estimées
IHPCmodele1$var.coef
IHPCmodele2$var.coef
IHPCmodele3$var.coef
IHPCmodele4$var.coef
IHPCmodele5$var.coef
IHPCmodele6$var.coef


### Code de convergence des modèles estimées
IHPCmodele1$code
IHPCmodele2$code
IHPCmodele3$code
IHPCmodele4$code
IHPCmodele5$code
IHPCmodele6$code

### Le modele 6 ne converge pas donc on continue avec les autres modele

### Finalement le modèle 1 estimé est : (I-B)(1-B-B^2)(X_t-mean(IHPCc))=(1-B-B^2)epsilon_t
### Finalement le modèle 2 estimé est : (I-B)(1-B-B^2)(X_t-mean(IHPC))=(1-B-B^2-B^3-B^4-B^5-B^6)epsilon_t
### Finalement le modèle 3 estimé est : (1-B-B^2-B^3-B^4-B^5-B^6)(X_t-mean(IHPC))=(1-B-B^2)epsilon_t
### Finalement le modèle 4 estimé est :
### Finalement le modèle 6 estimé est :

######## VALIDATION DU MODELE ########

##Test de significativité des paramètres
abs(IHPCmodele1$coef)>1.96*sqrt(diag(IHPCmodele1$var.coef))
abs(IHPCmodele2$coef)>1.96*sqrt(diag(IHPCmodele2$var.coef))
abs(IHPCmodele3$coef)>1.96*sqrt(diag(IHPCmodele3$var.coef))
abs(IHPCmodele4$coef)>1.96*sqrt(diag(IHPCmodele4$var.coef))
abs(IHPCmodele5$coef)>1.96*sqrt(diag(IHPCmodele5$var.coef))

# ou
library(lmtest)
coeftest(IHPCmodele1)
coeftest(IHPCmodele2)
coeftest(IHPCmodele3)
coeftest(IHPCmodele4)
coeftest(IHPCmodele5)

# SARIMA 
sarima_11=Arima(IHPCciv2018c,order=c(1,1,0),seasonal = list(order=c(0,1,0), period=12)) 
sarima_22=Arima(IHPCciv2018c,order=c(1,1,1),seasonal = list(order=c(0,1,0), period=12)) 
sarima_33=Arima(IHPCciv2018c,order=c(1,1,2),seasonal = list(order=c(0,1,0), period=12)) 
sarima_44=Arima(IHPCciv2018c,order=c(2,1,1),seasonal = list(order=c(0,1,0), period=12)) 
sarima_55=Arima(IHPCciv2018c,order=c(2,1,2),seasonal = list(order=c(0,1,0), period=12)) 

coeftest(sarima_11) 
coeftest(sarima_22)
coeftest(sarima_33) 
coeftest(sarima_44)
coeftest(sarima_55)

# SARIMA
sarima_111=Arima(IHPCciv2018c,order=c(1,1,0),seasonal = list(order=c(1,1,0), period=12)) 
sarima_221=Arima(IHPCciv2018c,order=c(1,1,1),seasonal = list(order=c(1,1,0), period=12)) 
sarima_331=Arima(IHPCciv2018c,order=c(1,1,2),seasonal = list(order=c(1,1,0), period=12)) 
sarima_441=Arima(IHPCciv2018c,order=c(2,1,1),seasonal = list(order=c(0,1,0), period=12)) 
sarima_551=Arima(IHPCciv2018c,order=c(2,1,2),seasonal = list(order=c(0,1,0), period=12)) 

coeftest(sarima_111) 
coeftest(sarima_221)
coeftest(sarima_331) 
coeftest(sarima_441)
coeftest(sarima_551)

############################### ARIMA AVEC INTERVENTION ########################


library("strucchange")

break_point <- breakpoints(IHPCciv2018 ~ 1) 

plot(IHPCciv2018)
lines(fitted(break_point, breaks = 2), col = 4)
lines(confint(break_point, breaks = 1))

summary(break_point)

# BIC minimal entre pour m=4

plot(IHPCciv2018)
fitted.ts <- fitted(break_point, breaks = 4)
lines(fitted.ts, col = 4)
lines(confint(break_point, breaks = 4))

unique(as.integer(fitted.ts))

breakdates(break_point, breaks = 4)

fitted.ts <- fitted(break_point, breaks = 4)
autoplot(fitted.ts)


checkresiduals(abhutondot_xreg)

# Mise en oeuvre de l'intervention
arima_Crise1=Arima(IHPCciv2018,order=c(1,1,0),xreg = fitted.ts, include.mean = TRUE) 
arima_Crise2=Arima(IHPCciv2018,order=c(1,1,1),xreg = fitted.ts, include.mean = TRUE) 
arima_Crise3=Arima(IHPCciv2018,order=c(1,1,2),xreg = fitted.ts, include.mean = TRUE) 
arima_Crise4=Arima(IHPCciv2018,order=c(2,1,1),xreg = fitted.ts, include.mean = TRUE) 
arima_Crise5=Arima(IHPCciv2018,order=c(2,1,2),xreg = fitted.ts, include.mean = TRUE) 

summary(arima_Crise1) 
summary(arima_Crise2)
summary(arima_Crise3) 
summary(arima_Crise4)
summary(arima_Crise5)

coeftest(arima_Crise1) 
coeftest(arima_Crise2)
coeftest(arima_Crise3) 
coeftest(arima_Crise4)
coeftest(arima_Crise5)

# SARIMA avec intervention
sarima_Crise11=Arima(IHPCciv2018c,order=c(1,1,0),seasonal = list(order=c(0,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise22=Arima(IHPCciv2018c,order=c(1,1,1),seasonal = list(order=c(0,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise33=Arima(IHPCciv2018c,order=c(1,1,2),seasonal = list(order=c(0,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise44=Arima(IHPCciv2018c,order=c(2,1,1),seasonal = list(order=c(0,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise55=Arima(IHPCciv2018c,order=c(2,1,2),seasonal = list(order=c(0,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 

summary(sarima_Crise11) 
summary(sarima_Crise22)
summary(sarima_Crise33) 
summary(sarima_Crise44)
summary(sarima_Crise55)

coeftest(sarima_Crise11) 
coeftest(sarima_Crise22)
coeftest(sarima_Crise33) 
coeftest(sarima_Crise44)
coeftest(sarima_Crise55)

#
sarima_Crise111=Arima(IHPCciv2018c,order=c(1,1,0),seasonal = list(order=c(1,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise221=Arima(IHPCciv2018c,order=c(1,1,1),seasonal = list(order=c(1,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise331=Arima(IHPCciv2018c,order=c(1,1,2),seasonal = list(order=c(1,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise441=Arima(IHPCciv2018c,order=c(2,1,1),seasonal = list(order=c(1,1,0), period=12),xreg = fitted.ts, include.mean = TRUE) 
sarima_Crise551=Arima(IHPCciv2018c,order=c(2,1,2),seasonal = list(order=c(1,1,0), period=12),xreg = fitted.ts, include.mean = TRUE)

coeftest(sarima_Crise111) 
coeftest(sarima_Crise221)
coeftest(sarima_Crise331) 
coeftest(sarima_Crise441)
coeftest(sarima_Crise551)

(model_1 <- auto.arima(IHPCciv2018c, stepwise = FALSE, trace = TRUE))

# On retient donc les modèles
# SARIMA (1,1,1)(0,1,0) :  sarima_22
# SARIMA (2,1,2)(0,1,0) :  sarima_55
# SARIMA (1,1,1)(1,1,0) :  sarima_221
# SARIMA (2,1,2)(1,1,0) :  sarima_551
# ARIMA-Intervention(1,1,1)  : arima_Crise2
# SARIMA-Intervention(1,1,1)(0,1,0) : sarima_Crise22
# SARIMA-Intervention(1,1,1)(1,1,0) : sarima_Crise221



############### ANALYSE DES RESIDUS #####################################

# Récuperation des residus
sarima_22.resid<-residuals(sarima_22)
sarima_55.resid<-residuals(sarima_55)
sarima_221.resid<-residuals(sarima_221)
sarima_551.resid<-residuals(sarima_551)
arima_Crise2.resid<-residuals(arima_Crise2)
sarima_Crise22.resid<-residuals(sarima_Crise22)
sarima_Crise221.resid<-residuals(sarima_Crise221)

# Nulluté de la moyenne des résidus

abs(mean(sarima_22.resid))<1.96*sqrt(var(sarima_22.resid))
abs(mean(sarima_55.resid))<1.96*sqrt(var(sarima_55.resid))
abs(mean(sarima_221.resid))<1.96*sqrt(var(sarima_221.resid))
abs(mean(sarima_551.resid))<1.96*sqrt(var(sarima_551.resid))
abs(mean(arima_Crise2.resid))<1.96*sqrt(var(arima_Crise2.resid))
abs(mean(sarima_Crise22.resid))<1.96*sqrt(var(sarima_Crise22.resid))
abs(mean(sarima_Crise221.resid))<1.96*sqrt(var(sarima_Crise221.resid))

# Test de non autocorrelation

Box.test(sarima_22.resid, lag = length(sarima_22.resid)/3, type = c("Box-Pierce", "Ljung-Box"))
Box.test(sarima_55.resid, lag = length(sarima_55.resid)/3, type = c("Box-Pierce", "Ljung-Box"))
Box.test(sarima_551.resid, lag = length(sarima_551.resid)/3, type = c("Box-Pierce", "Ljung-Box"))
Box.test(arima_Crise2.resid, lag = length(arima_Crise2.resid)/3, type = c("Box-Pierce", "Ljung-Box"))
Box.test(sarima_Crise22.resid, lag = length(sarima_Crise22.resid)/3, type = c("Box-Pierce", "Ljung-Box"))
Box.test(sarima_Crise221.resid, lag = length(sarima_Crise221.resid)/3, type = c("Box-Pierce", "Ljung-Box"))


BP=function(h) Box.test(sarima_22.resid,lag=h,type="Box-Pierce")$p.value
LB=function(h) Box.test(sarima_22.resid,lag=h,type="Ljung-Box")$p.value
plot(1:24,Vectorize(LB)(1:24),ylim=c(0,1),type="b",col="blue")
points(1:24,Vectorize(BP)(1:24),ylim=c(0,1),type="b",col="red",pch=2)
abline(h=.05,lty=2)
legend(20,.4,c("Box-Pierce", "Ljung-Box"),col=c("blue","red"),lty=1,pch=c(1,2))

BP=function(h) Box.test(arima_Crise2.resid,lag=h,type="Box-Pierce")$p.value
LB=function(h) Box.test(arima_Crise2.resid,lag=h,type="Ljung-Box")$p.value
plot(1:24,Vectorize(LB)(1:24),ylim=c(0,1),type="b",col="blue")
points(1:24,Vectorize(BP)(1:24),ylim=c(0,1),type="b",col="red",pch=2)
abline(h=.05,lty=2)
legend(20,.4,c("Box-Pierce", "Ljung-Box"),col=c("blue","red"),lty=1,pch=c(1,2))


library(car)
?durbin.watson

library(lmtest)
?dwtest

dwtest(sarima_55.resid ~ 1)
dwtest(sarima_22.resid ~ 1)
dwtest(arima_Crise1.resid ~ 1)
dwtest(sarima_Crise55.resid ~ 1)


# Normalité des erreurs

qqnorm(sarima_22.resid)
qqnorm(sarima_55.resid)
qqnorm(sarima_221.resid)
qqnorm(sarima_551.resid)
qqnorm(arima_Crise2.resid)
qqnorm(sarima_Crise22.resid)
qqnorm(sarima_Crise221.resid)

shapiro.test(sarima_22.resid)
shapiro.test(sarima_55.resid)
shapiro.test(sarima_221.resid)
shapiro.test(sarima_551.resid)
shapiro.test(arima_Crise2.resid)
shapiro.test(sarima_Crise22.resid)
shapiro.test(sarima_Crise221.resid)

jarque.bera.test(sarima_22.resid)
jarque.bera.test(sarima_55.resid)
jarque.bera.test(sarima_221.resid)
jarque.bera.test(sarima_551.resid)
jarque.bera.test(arima_Crise2.resid)
jarque.bera.test(sarima_Crise22.resid)
jarque.bera.test(sarima_Crise221.resid)


# Graphe global
checkresiduals(sarima_22)
checkresiduals(sarima_Crise221)


# Analyse global des résidus
tsdiag(sarima_22, plot=T)
tsdiag(sarima_55, plot=T)
tsdiag(sarima_Crise22, plot=T)
tsdiag(arima_Crise2, plot=T)

############### Donc seule sarima_Crise55 est normale et arima_Crise1 semble tende vers la normale
##### arima_Crise1.resid

#########Normalité des résidus############## 
plot(sarima_Crise55.resid, col="light blue")
abline(a=0, b=0) 

hist(sarima_Crise55.resid, breaks=40, col="pink")
plot(density(sarima_Crise55.resid))

hist(sarima_Crise55.resid, breaks=30, col="pink", freq=F)
#plot(density(ihpcresid))
lines(density(sarima_Crise55.resid))
Z=seq(min(sarima_Crise55.resid), 100, 0.01)
lines(Z, dnorm(Z,mean(sarima_Crise55.resid), sqrt(var(sarima_Crise55.resid))), col="red")

qqnorm(sarima_Crise55.resid) 
qqline(sarima_Crise55.resid)


### Donc nous allons continuer avec sarima_Crise55.resid et 
#Critère de sélection


RMSE1=sqrt(var(sarima_22.resid))
RMSE2=sqrt(var(sarima_55.resid))
RMSE3=sqrt(var(sarima_221.resid))
RMSE4=sqrt(var(sarima_551.resid))
RMSE5=sqrt(var(arima_Crise2.resid))
RMSE6=sqrt(var(sarima_Crise22.resid))
RMSE7=sqrt(var(sarima_Crise221.resid))
RMSE1
RMSE2
RMSE3
RMSE4
RMSE5
RMSE6
RMSE7

AIC1=log(sarima_Crise55$sigma2)+2*(2+2)/length(sarima_Crise55.resid)
AIC2=log(arima_Crise1$sigma2)+2*(1+1)/length(arima_Crise1.resid)
AIC3=log(sarima_Crise55$sigma2)+2*(2+2)/length(sarima_Crise55.resid)
AIC4=log(arima_Crise1$sigma2)+2*(1+1)/length(arima_Crise1.resid)
AIC5=log(sarima_Crise55$sigma2)+2*(2+2)/length(sarima_Crise55.resid)
AIC6=log(arima_Crise1$sigma2)+2*(1+1)/length(arima_Crise1.resid)
AIC7=log(sarima_Crise55$sigma2)+2*(2+2)/length(sarima_Crise55.resid)
AIC1
AIC2
AIC3
AIC4
AIC5
AIC6
AIC7

BIC1=log(sarima_Crise55$sigma2)+(2+2)*log(length(sarima_Crise55.resid))/length(sarima_Crise55.resid)
BIC2=log(arima_Crise1$sigma2)+(1+1)*log(length(arima_Crise1.resid))/length(arima_Crise1.resid)
BIC3=log(sarima_Crise55$sigma2)+(2+2)*log(length(sarima_Crise55.resid))/length(sarima_Crise55.resid)
BIC4=log(arima_Crise1$sigma2)+(1+1)*log(length(arima_Crise1.resid))/length(arima_Crise1.resid)
BIC5=log(sarima_Crise55$sigma2)+(2+2)*log(length(sarima_Crise55.resid))/length(sarima_Crise55.resid)
BIC6=log(arima_Crise1$sigma2)+(1+1)*log(length(arima_Crise1.resid))/length(arima_Crise1.resid)
BIC7=log(sarima_Crise55$sigma2)+(2+2)*log(length(sarima_Crise55.resid))/length(sarima_Crise55.resid)

BIC1
BIC2
BIC3
BIC4
BIC5
BIC6
BIC7


# BIC ET AIC des modèles

sarima_22$aic
sarima_55$aic
sarima_221$aic
sarima_551$aic
arima_Crise2$aic
sarima_Crise22$aic
sarima_Crise221$aic

sarima_22$bic
sarima_55$bic
sarima_221$bic
sarima_551$bic
arima_Crise2$bic
sarima_Crise22$bic
sarima_Crise221$bic

### Donc le modèle SARIMA_Crise221 se demarque des autres selon les critères d'information
######################### PREVISION  avec sarima_Crise55 ##################
###### predict(busbanc, n.ahead=h) ######## 
#ARMA
h_fut <- 4
plot(forecast(IHPCmodele5, h = h_fut))

Prevision=predict(IHPCmodele5,4)

##
IHPChprev<-Prevision$pred+mean(IHPCciv2018)
IHPChprevSE<-Prevision$se 

autoplot(IHPChprev)

IHPCprev<-window(IHPChprev,start=c(2019, 1), units="months", frequency=12)

borneinf<-IHPChprev-1.96*IHPChprevSE
bornesup<-IHPChprev+1.96*IHPChprevSE


ts.plot(window(IHPCciv, start=c(2019,1)), IHPChprev, bornesup, borneinf,lty=c(1,3,2,2), main="Prévision à l'horizon h=4")

lines(IHPCprev, col='blue', lty=3)
lines(borneinf, col='red', lty=3) 
lines(bornesup, col='red', lty=3) 

# SARIMA
h_fut <- 4
plot(forecast(sarima_221, h = h_fut))

Prevision=predict(sarima_221,4)

##
IHPChprev<-Prevision$pred+mean(IHPCciv2018)
IHPChprevSE<-Prevision$se 

autoplot(IHPChprev)

#==============================================================================#
#        PARTIE : Modélisation Multivarié de l'IHPC                            #
#==============================================================================#

##########################  Modelisation MCE  #############################
#### Packages necessaires
library(timeSeries)
library(stats)
library(graphics)
library(tseries)
library(foreign)
library(forecast)
library(urca)
library(vars)

install.packages("fUnitRoots")
library(fUnitRoots)


#### Importation de la base

baseMCE <- read.csv2("tsRAPPORT.csv", head=T, sep=";", dec=",")

#### Recuperations des differents series
IHPCg <- baseMCE$IPC._GLOBAL
IHPCpa <-baseMCE$IPC._Produits_alimentaires
IHPCbt <-baseMCE$IPC_Boissons_taba 
IHPCha <- baseMCE$IPC_habillement 
IHPClog <-baseMCE$IPC_logement
IHPCae <-baseMCE$IPC_Ameublemen_equip
IHPCls <- baseMCE$IPC_Loisirs_spectacles
IHPChr <- baseMCE$IPC._Hotels_restaurants
IHPCe <- baseMCE$IPC_Enseignement
IHPCnca <- baseMCE$IPC._Autres.biens_services


IHPCs <- baseMCE$IPC._Sante
IHPCc<- baseMCE$IPC_Communication
IHPCt <- baseMCE$IPC_Transports


# Mise en forme des différents séries

IHPCg <- ts(IHPCg, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCpa <- ts(IHPCpa, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCbt <- ts(IHPCbt, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCha <- ts(IHPCha, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPClog <- ts(IHPClog, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCae <- ts(IHPCae, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCls <- ts(IHPCls, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPChr <- ts(IHPChr, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCe <- ts(IHPCe, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCnca <- ts(IHPCnca, start = c(2009, 1), frequency = 12, end=c(2018,12))

IHPCs <- ts(IHPCs, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCc <- ts(IHPCc, start = c(2009, 1), frequency = 12, end=c(2018,12))
IHPCt <- ts(IHPCt, start = c(2009, 1), frequency = 12, end=c(2018,12))

# Passage au logarithme des series 
IHPCg <- log(IHPCg)
IHPCpa <- log(IHPCpa)
IHPCbt <- log(IHPCbt)
IHPCha <- log(IHPCha)
IHPClog <- log(IHPClog)
IHPCae <- log(IHPCae)
IHPCls <- log(IHPCls)
IHPChr <- log(IHPChr)
IHPCe <- log(IHPCe)
IHPCnca <- log(IHPCnca)



IHPCg <- log(IHPCg)
IHPCpa <- log(IHPCpa)
IHPCbt <- log(IHPCbt)
IHPCha <- log(IHPCha)
IHPClog <- log(IHPClog)
IHPCae <- log(IHPCae)
IHPCls <- log(IHPCls)
IHPChr <- log(IHPChr)
IHPCe <- log(IHPCe)
IHPCnca <- log(IHPCnca)

data=data.frame(IHPCg,IHPCpa,IHPCbt,IHPCha,IHPClog,IHPCae,IHPCls,IHPChr,IHPCe,IHPCnca)
data$obs <- 1:120

trend1 <- as.numeric(data$obs==28)
trend2 <- data$obs

model <- lm(IHPCg~trend1+trend2)
model1 <- lm(IHPCpa~trend1+trend2)
model2 <- lm(IHPCbt~trend1+trend2)
model3 <- lm(IHPCha~trend1+trend2)
model4 <- lm(IHPClog~trend1+trend2)
model5 <- lm(IHPCae~trend1+trend2)
model6 <- lm(IHPCls~trend1+trend2)
model7 <- lm(IHPChr~trend1+trend2)
model8 <- lm(IHPCe~trend1+trend2)
model9 <- lm(IHPCnca~trend1+trend2)

summary(model)

bIHPCg <- residuals(model)
bIHPCpa <- residuals(model1)
bIHPCbt <- residuals(model2)
bIHPCha <- residuals(model3)
bIHPClog <- residuals(model4)
bIHPCae <- residuals(model5)
bIHPCls <- residuals(model6)
bIHPChr <- residuals(model7)
bIHPCe <- residuals(model8)
bIHPCnca <- residuals(model9)

bIHPCg <- log(bIHPCg)
bIHPCpa <- log(bIHPCpa)
bIHPCbt <- log(bIHPCbt)
bIHPCha <- log(bIHPCha)
bIHPClog <- log(bIHPClog)
bIHPCae <- log(bIHPCae)
bIHPCls <- log(bIHPCls)
bIHPChr <- log(bIHPChr)
bIHPCe <- log(bIHPCe)
bIHPCnca <- log(bIHPCnca)

####### Test de stationnrité sur les séries loguées

## Test de dickey-fuller
adfTest(IHPCg,lags=1)
adfTest(IHPCpa,lags=1)
adfTest(IHPCbt,lags=1)
adfTest(IHPCha,lags=1)
adfTest(IHPClog,lags=1)
adfTest(IHPCae,lags=1)
adfTest(IHPCls,lags=1)
adfTest(IHPChr,lags=1)
adfTest(IHPCe,lags=1)
adfTest(IHPCnca,lags=1)

adfTest(IHPCs,lags=1) 
adfTest(IHPCc,lags=1)
adfTest(IHPCt,lags=1)

# PP
library(urca)
urppTest(IHPCg)
urppTest(IHPCpa)
urppTest(IHPCbt)
urppTest(IHPCha)
urppTest(IHPClog)
urppTest(IHPCae)
urppTest(IHPCls)
urppTest(IHPChr)
urppTest(IHPCe)
urppTest(IHPCnca)

## Differenciation des series
bIHPCgd=diff(bIHPCg,lag=1)
bIHPCpad=diff(bIHPCpa,lag=1)
bIHPCbtd=diff(bIHPCbt,lag=1)
bIHPChad=diff(bIHPCha,lag=1)
bIHPClogd=diff(bIHPClog,lag=1)
bIHPCaed=diff(bIHPCae,lag=1)
bIHPClsd=diff(bIHPCls,lag=1)
bIHPChrd=diff(bIHPChr,lag=1)
bIHPCed=diff(bIHPCe,lag=1)
bIHPCncad=diff(bIHPCnca,lag=1)

IHPCcd=diff(IHPCc,lags=1)
IHPCtd=diff(IHPCt,lags=1)

## Tets de stationnarité sur les séries différenciés
adfTest(bIHPCgd)
adfTest(IHPCpad,lags=1)
adfTest(IHPCbtd,lags=1)
adfTest(IHPChad,lags=1)
adfTest(IHPClogd,lags=1)
adfTest(IHPCaed,lags=1)
adfTest(IHPClsd,lags=1)
adfTest(IHPChrd,lags=1)
adfTest(IHPCed,lags=1)
adfTest(IHPCncad,lags=1)

adfTest(IHPCcd)
adfTest(IHPCtd)

##### Tous les séries différencies sont donc stationnaire et les series
##### loguées sont integré d'ordre 1

######################### Test de cointégration entre les variables#######

## Test de cointégration d'Engle et Granger

### Entre IHPCgd et IHPCpad

# Etape 1
reg1=lm(IHPCgd~IHPCpad)
reg1.resid=reg1$residuals

# Etape 2
library(urca)
y=ur.df(reg1.resid,type="none",selectlags = "AIC")
summary(y)
y@teststat
y@cval

### Test de johasen

# Selection du lag approprié
library(vars)

baseMCEt<-cbind(bIHPCg,bIHPCpa,bIHPCbt,bIHPCha,bIHPClog,bIHPCae,bIHPCls,bIHPChr,bIHPCe,bIHPCnca)
VARselect(baseMCEt)

mod=lm(bIHPCg~bIHPCpa+bIHPCbt+bIHPCha+bIHPClog+bIHPCae+bIHPCls+bIHPChr+bIHPCe+bIHPCnca)

summary(mod)
baseMCE2<-cbind(IHPCg,IHPCpa,IHPCae,IHPCls,IHPCc,IHPCt)
# Test de cointegration
xy.cit <- ca.jo(baseMCEt,ecdet="const",type="eigen", K=2, spec="transitory" )

xy.cit
xy.cit@teststat[2]  # test statistics H0 r=0
xy.cit@teststat[1]  # test statistiques H0 r=1

# Les valeurs critiques
xy.cit@cval

## Estimation des paramètres
vecm<-cajorls(xy.cit)
summary(vecm$rlm)

# estraction des termes d'erreurs
vecm$rlm$coefficient[1,]

# Relation de long terme
vecm$beta

### Test 2
vecm2=ca.jo(baseMCEt, type="trace", K=2, ecdet="trend", spec="longrun")
summary(vecm2)
vecm2@teststat[2]  # test statistics H0 r=0
vecm2@teststat[1]  # test statistiques H0 r=1

# Statistique de trace LRtrace
vecm2@teststat  

## Estimation des paramètres
vecm<-cajorls(vecm2)
summary(vecm$rlm)

# estraction des termes d'erreurs
vecm$rlm$coefficient[1,]

# Relation de long terme
vecm$beta

### Transformation en var
ban=vec2var(vecm2)


### Prediction
predictin=predict(ban, n.ahead = 12, ci = 0.95)
x11()
plot(predictin)

### Impulsion
svec.irf <- irf(ban, response = "bIHPCg", n.ahead = 4, boot = TRUE)
plot(svec.irf)


### decomposition de la variance de prevision
fevd <- fevd(ban, n.ahead = 4)$bIHPCg
plot(fevd)





############################################################## renommage des variables
IHPC <- basets$IPC._GLOBAL
IHPPA <-  IPC._Produits_alimentaires
IHPBT <-  IPC_Boissons_taba 
IHPHabil <- IPC_habillement 
IHPLOGE <-IPC_logement
IHPAmEq <- IPC_Ameublemen_equip
IHPSant <- IPC._Sante
IHPCom <- IPC_Communication
IHPTrsp <- IPC_Transports

IHPLS <- IPC_Loisirs_spectacles
IHPHR <- IPC._Hotels_restaurants
IHPEnseig <- IPC_Enseignement
IHPABS <- IPC._Autres.biens_services
summary(IHPC)
####################################################################
IHPC = log(IHPC)
IHPPA <-  log(IHPPA)
IHPBT <-  log(IHPBT)
IHPHabil <- log(IHPHabil)
IHPLOGE <- log(IHPLOGE)
IHPAmEq <- log(IHPAmEq)
IHPSant <- log(IHPSant)
IHPCom <- log(IHPCom)
IHPTrsp <- log(IHPTrsp)
IHPLS <- log(IHPLS)
IHPHR <- log(IHPHR)
IHPEnseig <- log(IHPEnseig)
IHPABS <- log(IHPABS)



####################### MODELISATION VAR ###############

##Installation des packages
library(stats)
library(tseries)
library(vars)
library(stats)
library(grDevices)
library(graphics)
library(forecast)
library(timeSeries)
library(urca)
library(foreign)


#************************* Pr?sentation des donn?es************************#
summary(basets)
str(basets)
#Declaration des s?ries comme temporelles de 2009  Ã  2019 avec une frÃ©quence mensuelle

IHPC <- ts(IHPC, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPPA <- ts(IHPPA, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPBT <- ts(IHPBT, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPHabil <- ts(IHPHabil, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPLOGE <- ts(IHPLOGE, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPAmEq <- ts(IHPAmEq, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPLS <- ts(IHPLS, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPHR <- ts(IHPHR, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPEnseig <- ts(IHPEnseig, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPABS <- ts(IHPABS, start = c(2009, 1), frequency = 12, end=c(2019,4))

IHPSant <- ts(IHPSant, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPCom <- ts(IHPCom, start = c(2009, 1), frequency = 12, end=c(2019,4))
IHPTrsp <- ts(IHPTrsp, start = c(2009, 1), frequency = 12, end=c(2019,4))



#**************************Etape 1 : Etudes prÃ©liminaires******************#
par(mfrow=c(5,2))
x11()
plot(IHPC)
x11()
plot(IHPPA)
ts.plot(IHPBT)
ts.plot(IHPHabil)
ts.plot(IHPLOGE)
ts.plot(IHPAmEq)
ts.plot(IHPLS)
ts.plot(IHPHR)
ts.plot(IHPEnseig)
ts.plot(IHPABS)

#********************** Centrage des sÃ©ries ******************************
IHPCc <- IHPC-mean(IHPC)
IHPPAc <- IHPPA-mean(IHPPA)
IHPBTc <- IHPBT-mean(IHPBT)
IHPHabilc <- IHPHabil-mean(IHPHabil)
IHPLOGEc <- IHPLOGE-mean(IHPLOGE)
IHPAmEqc <- IHPAmEq-mean(IHPAmEq)
IHPLSc <- IHPLS-mean(IHPLS)
IHPHRc <- IHPHR-mean(IHPHR)
IHPEnseigc <- IHPEnseig-mean(IHPEnseig)
IHPABSc <- IHPABS-mean(IHPABS)

IHPSantc <- IHPSant-mean(IHPSant)
IHPComc <- IHPCom-mean(IHPCom)
IHPTrspc <- IHPTrsp-mean(IHPTrsp)

#************************TEST DE STATIONNARITE********************#

#*********** Variable : Indice HarmonisÃ© des Prix Ã  la Consommation
ts.plot(IHPCc)
kpss.test(IHPCc)
#le test de kpss donne une p-value de 0.01 donc on rejette la stationnaritÃ©

#*********** Variable :Indice de Production Alimentaire
ts.plot(IHPPAc)
kpss.test(IHPPAc)
# p-value=0.01 donc la sÃ©rie est non stationnaire

#*********** Variable :Indice de Boissons et Tabac
ts.plot(IHPBTc)
kpss.test(IHPBTc)
# p-value=0.01 donc la sÃ©rie est non stationnaire

#*********** Variable :Indice d'Habillement
ts.plot(IHPHabilc)
kpss.test(IHPHabilc)
# p-value=0.01 donc la sÃ©rie n'est pas  stationnaire

#*********** Variable : Indice de Logement
ts.plot(IHPLOGEc)
kpss.test(IHPLOGEc)
# p-value=0.01 donc la sÃ©rie n'est pas stationnaire

#*********** Variable : Indice d'Ameublement et Ã©quipement
ts.plot(IHPAmEqc)
kpss.test(IHPAmEqc)
# p-value=0.01 donc la sÃ©rie est  stationnaire


#*********** Variable :Indice de Loisirs et Spectacles
ts.plot(IHPLSc)
kpss.test(IHPLSc)
# p-value=0.01 donc la sÃ©rie n'est pas  stationnaire

#*********** Variable : Indice d'HÃ´tel et Restaurant
ts.plot(IHPHRc)
kpss.test(IHPHRc)
# p-value=0.01 donc la sÃ©rie n'est pas stationnaire

#*********** Variable : Indice d'Enseignement
ts.plot(IHPEnseigc)
kpss.test(IHPEnseigc)
# p-value=0.01 donc la sÃ©rie est  non stationnaire

#*********** Variable : Indice d'autres Biens et Services
ts.plot(IHPABSc)
kpss.test(IHPABSc)
# p-value=0.01 donc la sÃ©rie n'est pas stationnaire



#*********** Variable : Indice de Sante
ts.plot(IHPSantc)
kpss.test(IHPSantc)
# p-value=0.1 donc la sÃ©rie est stationnaire

#*********** Variable : Indice de communication
ts.plot(IHPComc)
kpss.test(IHPComc)
# p-value=0.01 donc la sÃ©rie n'est pas non stationnaire

#*********** Variable : Indice de transport
ts.plot(IHPTrspc)
kpss.test(IHPTrspc)
# p-value=0.02572 donc la sÃ©rie n'est pas stationnaire


#******************Stationnarisation de la sÃ©rie Indice HarmonisÃ© des Prix Ã  la consommation
IHPCcd<-diff(IHPCc)
ts.plot(IHPCcd)
kpss.test(IHPCcd)
# On constate qu'aprÃ¨s diffÃ©renciation la sÃ©rie devient stationnaire, intÃ©grÃ© d'ordre 1

IHPPAcd<-diff(IHPPAc)
ts.plot(IHPPAcd)
kpss.test(IHPPAcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPBTcd <- diff(IHPBTc)
x11()
ts.plot(IHPBTcd)
kpss.test(IHPBTcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1



IHPHabilcd <- diff(IHPHabilc)
kpss.test(IHPHabilcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPLOGEcd <- diff(IHPLOGEc)
ts.plot(IHPLOGEcd)
kpss.test(IHPLOGEcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPAmEqcd <- diff(IHPAmEqc)
kpss.test(IHPAmEqcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1



IHPLScd <- diff(IHPLSc)
kpss.test(IHPLScd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPHRcd <- diff(IHPHRc)
kpss.test(IHPHRcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPEnseigcd <- diff(IHPEnseigc)
kpss.test(IHPEnseigcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPABScd <- diff(IHPABSc)
kpss.test(IHPABScd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1



IHPComcd <- diff(IHPComc)
kpss.test(IHPComcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1


IHPTrspcd <- diff(IHPTrspc)
kpss.test(IHPTrspcd)
# Ici Ã©galement la sÃ©rie devient stationnaire aprÃ¨s diffÃ©renciation, intÃ©grÃ© d'ordre 1





#**************************Etapes Identification ******************#


## Matrice de corrÃ©lation
cor(cbind(as.vector(IHPCcd),as.vector(IHPPAcd)))

library(corrplot)
M=cor(cbind(as.vector(IHPCcd),as.vector(IHPPAcd),as.vector(IHPBTcd),as.vector(IHPComcd),as.vector(IHPTrspcd),as.vector(IHPHabilcd),as.vector(IHPLOGEcd),as.vector(IHPAmEqcd),as.vector(IHPLScd),as.vector(IHPHRcd),as.vector(IHPEnseigcd),as.vector(IHPABScd)))
x11()
corrplot(M, method = "circle")
#* on va prendre Indice de production alimentaire, Boissons& Tabac et logement

#***********Selection de l'ordre p du VAR
# Pour des soucis de corrÃ©lation avec la variable indice global et indice de la production alimentaire, on a enlevÃ© les variables communication, santÃ© et transport
base_var1<-cbind(as.vector(IHPCcd),as.vector(IHPPAcd),as.vector(IHPBTcd),as.vector(IHPLOGEcd))
base_var2<-cbind(IHPCcd,IHPPAcd,IHPBTcd,IHPLOGEcd,IHPComcd)
x11()
plot(base_var1)
plot(base_var2)

#Detection de l'ordre optimal
library(vars)
VARselect(base_var2,lag=12)$selection

#Candidat 
#  p=1;3;12

##Estimation du VAR 

######## Validation du mod?le
model1<-vars::VAR(base_var2,p=12, type="both")
model2<-vars::VAR(base_var2,p=3, type="both")

summary(model1)

# Les variables communication, logement, Production Alimentaire, Boissons et Tabac

#ReprÃ©sentation de l'inverse des racines
roots(model1)
roots(model2)

#####Analyse des r?sidus

#Normalit?

normality.test(model1)

## Les erreurs sont normales pour le modÃ¨le VAR(12)


# Tests de non autocorrelation des erreurs
serial.test(model1, lags.pt = 5, type ="PT.adjusted")
serial.test(model2, lags.pt = 5, type ="PT.adjusted")
serial.test(model3, lags.pt = 5, type ="PT.adjusted")# le PT.ajustÃ© envoie une erreur
serial.test(model4, lags.pt = 5, type ="PT.adjusted")# le PT.ajustÃ© envoie une erreur
serial.test(model5, lags.pt = 5, type ="PT.adjusted")# le PT.ajustÃ© envoie une erreur
serial.test(model6, lags.pt = 5, type ="PT.adjusted")# le PT.ajustÃ© envoie une erreur

serial.test(model1, type="BG")


serial.test(model1, lags.pt = 16, type = "PT.asymptotic")
serial.test(model6, lags.pt = 16, type = "PT.asymptotic")
################################################################# on choisi le mod?le 1; VAR(12)

## Seul le modÃ¨le VAR(12) a des rÃ©sidus non autocorrÃ©lÃ©s. On utilisera donc p=12 comme p optimal. Puisque -----
## d_max =1 (car aucune sÃ©rie n'est I(2)) on estime donc par la mÃ©thodologie de Toda et Yamamoto, le modÃ¨le ----
## VAR(13) ----

## ModÃ©lisation selon Toda-Yamamoto ---- boisson et tabac, logement, habillement, hÃ´tel & restaurant
var.IHPC.mod.est <-  vars::VAR(base_var2,p=13, type="both")
normality.test(var.IHPC.mod.est) # RÃ©sidus normaux 
serial.test(var.IHPC.mod.est, lags.pt = 5, type ="PT.asymptotic") # RÃ©sidus non autocorrelÃ©s

summary(var.IHPC.mod.est)


#Causalit? et Stabilit?
causality(model1,cause ="IHPCcd") 
causality (model1,cause = "IHPLOGEcd")
causality(model1,cause = "IHPComcd")
causality(model1,cause = "IHPPAcd")
causality(model1,cause = "IHPBTcd")
x11()
plot(model1, names = "IHPCcd")
chartSeries(IHPCcd)
plot(model1)



#Fonctions de r?ponses impulsionnelles

impul1 <- irf(model1, reponse="IHPCcd" , n.ahead = 30, boot = TRUE)
plot(impul1)

#DÃ©composition de la variance
fevd.IHPCcd <- fevd(model1, n.ahead = 20)$IHPCcd
fevd.IHPCcd

#PREVISION

prevision <- predict(model1, n.ahead=71)
par(mfrow=c(1,1))
ts.plot(prevision)

########################################"

#Repr?sentation graphique
prevision2 <- predict(model1, n.ahead=150)
x11()
plot(prevision2)$IHPCcd

plot(prevision)$pib_tdf


fanchart(prevision)
#################################################### coint?gration
mod1=lm(IHPCcd~IHPLOGEcd)
mod2=lm(IHPCcd~IHPComcd)
mod3=lm(IHPCcd~IHPPAcd)
mod4=lm(IHPCcd~IHPBTcd)
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)
res1 <- mod1$residuals
res2 <- mod2$residuals
res3 <- mod3$residuals
res4 <- mod4$residuals
x11()
acf(res1)
adf.test(res1)

# DÃ©composition de la variance
fevd.Freq <- fevd(var.IHPC.mod.est, n.ahead = 18)
fevd.Freq

# PrÃ©visions sur une annÃ©e
prevision <- predict(var.IHPC.mod.est, n.ahead=12)
x11()
plot(prevision)

