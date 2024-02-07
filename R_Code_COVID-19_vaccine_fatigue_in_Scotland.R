#this file only contains the code for fitting the models for which the results are presented in Table 1.
#if you have other questions, please contact the corresponding author at robin.muegge@glasgow.ac.uk
library(readxl)
library(ggplot2)
library(INLA) #INLA needs to be installed manually (see https://www.r-inla.org/).
library(spdep)
library(spatialreg)
library(sf)
library(readxl)
library(gridExtra)
library(cowplot)
library(grid)
library(pkgbuild)

#set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load the data. The data are publicly available at https://www.opendata.nhs.scot/dataset/covid-19-vaccination-in-scotland
vacc <- read_excel("vaccdata_paper.xlsx",sheet = "Data")
areas <- as.data.frame(unique(vacc[,"CA"]))
areas$idarea <- 1:nrow(areas)
colnames(areas) <- c("CA","idarea")
vacc <- merge(vacc,areas)

ages <- as.data.frame(unique(vacc[,"AgeGroup"]))
ages$idage <- 1:nrow(ages)
colnames(ages) <- c("AgeGroup","idage")
vacc <- merge(vacc,ages)

sex_dose <- as.data.frame(unique(vacc[,c("Sex","d")]))
sex_dose$idsexdose <- c(1,4,3,2)
colnames(sex_dose) <- c("Sex","d","idsexdose")
vacc <- merge(vacc,sex_dose)
vacc$idsexdose2 <- vacc$idsexdose

vacc$iddose <- vacc$d-1
vacc$iddose2 <- vacc$iddose
vacc$idsex <- vacc$Sex
vacc$idsex[which(vacc$idsex=="Female")] <- 1
vacc$idsex[which(vacc$idsex=="Male")] <- 2
vacc$idsex <- as.integer(vacc$idsex)
vacc$idsex2 <- vacc$idsex

vacc <- vacc[order(vacc$CA,vacc$Sex,vacc$AgeGroup),]

###############
#Model fitting#
###############

#the shapefile can be obtained from https://geoportal.statistics.gov.uk/datasets/c370f21b4b4649a5b6813bf48469836f_0/explore?location=55.215415%2C-3.313872%2C6.14
Bdry <- st_read("Local_Authority_Districts__May_2020__Boundaries_UK_BFC.shp")
Bdry2 <- st_transform(x=Bdry, crs='+proj=longlat +datum=WGS84 +no_defs')
Vacc.Bdry <- merge(x=Bdry2,y=areas,by.x="LAD20CD",by.y="CA",all.x=FALSE)
sf_use_s2(FALSE) 
W.nb <- poly2nb(Vacc.Bdry,row.names=Vacc.Bdry$LAD20CD)
summary(W.nb)

#add links between islands and nearest council areas:
W.nb[[6]] <- as.integer(8) #make Highlands neighbour of western isles Na h-Eileanan Siar
W.nb[[13]] <- as.integer(c(8,15)) #make Highlands and Shetland islands neighbour of Orkney islands
W.nb[[15]] <- as.integer(13) #make Orkney islands neighbour of Shetland islands
W.nb[[8]] <- as.integer(c(sort(c(W.nb[[8]],6,13)))) #add Na h-Eileanan Siar and Orkney islands to neighbours of the Highlands
nb2INLA("map.adj",W.nb)
W <- inla.read.graph(filename="map.adj")

#create structure matrix for age groups
W_age <- matrix(c(0,1,rep(0,8),
                  1,0,1,rep(0,7),
                  0,1,0,1,rep(0,6),
                  0,0,1,0,1,rep(0,5),
                  rep(0,3),1,0,1,rep(0,4),
                  rep(0,4),1,0,1,rep(0,3),
                  rep(0,5),1,0,1,rep(0,2),
                  rep(0,6),1,0,1,0,
                  rep(0,7),1,0,1,
                  rep(0,8),1,0),nrow=10,byrow=T) 
nb_age <- mat2listw(W_age)
nb2INLA("age.adj",nb_age$neighbours)
W_age <- inla.read.graph(filename="age.adj")

#change initial value from -3 to 3 (needed for some models that don't converge otherwise)
prior_new <- list(
  prec = list(
    prior = "pc.prec",
    param = c(1,0.01)),
  phi = list(
    prior = "pc",
    param = c(0.5,0.5),
    initial=3)
)

#1st row of Table 1
res_binom_sex <- inla(NumVacc_d ~ factor(Sex) + factor(d),
                      family="binomial",
                      Ntrials=NumVacc_d_minus_1,
                      data=vacc,
                      control.inla=list(strategy="laplace"),
                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                      control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex)

#2nd row of Table 1
res_binom_trn <- inla(NumVacc_d ~ factor(d),
                       family="binomial",
                       Ntrials=NumVacc_d_minus_1,
                       data=vacc,
                       control.inla=list(strategy="laplace"),
                       control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                       control.fixed = list(prec.intercept=0.001))
summary(res_binom_trn)

#3rd row of Table 1
res_binom_age <- inla(NumVacc_d ~ f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new),
                      family="binomial",
                      Ntrials=NumVacc_d_minus_1,
                      data=vacc,
                      control.inla=list(strategy="laplace"),
                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                      control.fixed = list(prec.intercept=0.001))
summary(res_binom_age)

#4th row of Table 1
res_binom_area <- inla(NumVacc_d ~ f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                       family="binomial",
                       Ntrials=NumVacc_d_minus_1,
                       data=vacc,
                       control.inla=list(strategy="laplace"),
                       control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                       control.fixed = list(prec.intercept=0.001))
summary(res_binom_area)

#5th row of Table 1
res_binom_age_sex <- inla(NumVacc_d ~ factor(Sex) + f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new),
                          family="binomial",
                          Ntrials=NumVacc_d_minus_1,
                          data=vacc,
                          control.inla=list(strategy="laplace"),
                          control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                          control.fixed = list(prec.intercept=0.001))
summary(res_binom_age_sex)

#6th row of Table 1
res_binom_age_trn <- inla(NumVacc_d ~ factor(d) + f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new),
                            family="binomial",
                            Ntrials=NumVacc_d_minus_1,
                            data=vacc,
                            control.inla=list(strategy="laplace"),
                            control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                            control.fixed = list(prec.intercept=0.001))
summary(res_binom_age_trn)

#7th row of Table 1
res_binom_age_area <- inla(NumVacc_d ~ f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new) + 
                             f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                           family="binomial",
                           Ntrials=NumVacc_d_minus_1,
                           data=vacc,
                           control.inla=list(strategy="laplace"),
                           control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                           control.fixed = list(prec.intercept=0.001))
summary(res_binom_age_area)

#8th row of Table 1
res_binom_age_trn_sex <- inla(NumVacc_d ~ factor(Sex) + factor(d) +
                                  f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new),
                                family="binomial",
                                Ntrials=NumVacc_d_minus_1,
                                data=vacc,
                                control.inla=list(strategy="laplace"),
                                control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                control.fixed = list(prec.intercept=0.001))
summary(res_binom_age_trn_sex)

#9th row of Table 1
res_binom_age_trn_area <- inla(NumVacc_d ~ factor(d) +
                                   f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE,hyper=prior_new) +
                                   f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                                 family="binomial",
                                 Ntrials=NumVacc_d_minus_1,
                                 data=vacc,
                                 control.inla=list(strategy="laplace"),
                                 control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                 control.fixed = list(prec.intercept=0.001))
summary(res_binom_age_trn_area)

#10th row of Table 1
res_binom_sex_trn_age_area <- inla(NumVacc_d ~ factor(Sex) + factor(d) +
                                      f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE) +
                                      f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                                    family="binomial",
                                    Ntrials=NumVacc_d_minus_1,
                                    data=vacc,
                                    control.inla=list(strategy="laplace"),
                                    control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                    control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_age_area)

#11th row of Table 1
res_binom_sex_trn_ageS_area <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                       f(idage,model="bym2",graph=W_age,replicate=idsex,scale.model=TRUE,constr=TRUE) +
                                       f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                                     family="binomial",
                                     Ntrials=NumVacc_d_minus_1,
                                     data=vacc,
                                     control.inla=list(strategy="laplace"),
                                     control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                     control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageS_area)

#12th row of Table 1
res_binom_sex_trn_ageD_area <- inla(NumVacc_d ~ factor(Sex) + factor(d) +
                                       f(idage,model="bym2",graph=W_age,replicate=iddose,scale.model=TRUE,constr=TRUE) +
                                       f(idarea,model="bym2",graph=W,scale.model=TRUE,constr=TRUE),
                                     family="binomial",
                                     Ntrials=NumVacc_d_minus_1,
                                     data=vacc,
                                     control.inla=list(strategy="laplace"),
                                     control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                     control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageD_area)

#13th row of Table 1
res_binom_sex_trn_age_areaS <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                       f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE) +
                                       f(idarea,model="bym2",graph=W,replicate=idsex,scale.model=TRUE,constr=TRUE),
                                     family="binomial",
                                     Ntrials=NumVacc_d_minus_1,
                                     data=vacc,
                                     control.inla=list(strategy="laplace"),
                                     control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                     control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_age_areaS)

#14th row of Table 1
res_binom_sex_trn_age_areaD <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                       f(idage,model="bym2",graph=W_age,scale.model=TRUE,constr=TRUE) +
                                       f(idarea,model="bym2",graph=W,replicate=iddose,scale.model=TRUE,constr=TRUE),
                                     family="binomial",
                                     Ntrials=NumVacc_d_minus_1,
                                     data=vacc,
                                     control.inla=list(strategy="laplace"),
                                     control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                     control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_age_areaD)

#15th row of Table 1
res_binom_sex_trn_ageS_areaD <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                        f(idage,model="bym2",graph=W_age,replicate=idsex,scale.model=TRUE,constr=TRUE) +
                                        f(idarea,model="bym2",graph=W,replicate=iddose,scale.model=TRUE,constr=TRUE),
                                      family="binomial",
                                      Ntrials=NumVacc_d_minus_1,
                                      data=vacc,
                                      control.inla=list(strategy="laplace"),
                                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                      control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageS_areaD)

#16th row of Table 1
res_binom_sex_trn_ageD_areaD <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                        f(idage,model="bym2",graph=W_age,replicate=iddose,scale.model=TRUE,constr=TRUE) +
                                        f(idarea,model="bym2",graph=W,replicate=iddose2,scale.model=TRUE,constr=TRUE),
                                      family="binomial",
                                      Ntrials=NumVacc_d_minus_1,
                                      data=vacc,
                                      control.inla=list(strategy="laplace"),
                                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                      control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageD_areaD)

#Model 1 (17th row of Table 1)
res_binom_sex_trn_ageSD_areaD <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                         f(idage,model="bym2",graph=W_age,replicate=idsexdose,scale.model=TRUE,constr=TRUE) +
                                         f(idarea,model="bym2",graph=W,replicate=iddose,scale.model=TRUE,constr=TRUE),
                                       family="binomial",
                                       Ntrials=NumVacc_d_minus_1,
                                       data=vacc,
                                       control.inla=list(strategy="laplace"),
                                       control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                       control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageSD_areaD)

#with IMD covariate
res_binom_sex_trn_ageSD_areaD <- inla(NumVacc_d ~ factor(Sex)+factor(d)+log(IMD_CA_rank)+
                                        f(idage,model="bym2",graph=W_age,replicate=idsexdose,scale.model=TRUE,constr=TRUE) +
                                        f(idarea,model="bym2",graph=W,replicate=iddose,scale.model=TRUE,constr=TRUE),
                                      family="binomial",
                                      Ntrials=NumVacc_d_minus_1,
                                      data=vacc,
                                      control.inla=list(strategy="laplace"),
                                      control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                      control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageSD_areaD)

#Model 2 (18th row of Table 1)
res_binom_sex_trn_ageSD_areaSD <- inla(NumVacc_d ~ factor(Sex)+factor(d) +
                                          f(idage,model="bym2",graph=W_age,replicate=idsexdose,scale.model=TRUE,constr=TRUE) +
                                          f(idarea,model="bym2",graph=W,replicate=idsexdose2,scale.model=TRUE,constr=TRUE),
                                        family="binomial",
                                        Ntrials=NumVacc_d_minus_1,
                                        data=vacc,
                                        control.inla=list(strategy="laplace"),
                                        control.compute=list(dic=TRUE,waic=TRUE,cpo=TRUE),
                                        control.fixed = list(prec.intercept=0.001))
summary(res_binom_sex_trn_ageSD_areaSD)

