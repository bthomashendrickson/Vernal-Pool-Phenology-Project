#Install Dependencies
install.packages(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools"))
lapply(c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools"), require, character.only = TRUE)
list.of.packages <- c("nortest","readr","sjPlot","ggplot2","lmtest","AICcmodavg","DescTools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#Load Dataframe
HOME="~/"
VP_DIR=paste(HOME,"Documents/Experiments/Vernal_Pool /",sep="")
Pheno_data <- read_csv(paste(VP_DIR,"VP_Phenology_data.csv",sep=""))
CIMIS <- read_csv(paste(VP_DIR,"CIMIS_SITES.csv",sep=""))
climate.data <- read_csv(paste(VP_DIR,"MERCED_CLIMATE_2000_2022.csv",sep=""))

#General: Climate trends over study period & 1999-2022, Climate comparison between CIMIS cites,
# and plotting pool characteristics.
water_year <- na.omit(climate.data[(climate.data$MONTH>=10 | climate.data$MONTH<4),])
GDH<-function(x,t,GDH){
    if(any(class(x)=="data.frame")==TRUE){
        output <- character (nrow(x))
        for (i in 1:nrow(x)){
            if(x[i,"TEMP"] > t){
                output[i]<-1
            } else {
                output[i]<-0
            }
        }
    } else {
        print("Not a data.frame")
    }
    x[[GDH]]<-output
    return(x)
}

water_year<-GDH(water_year,10,GDH="GDH_MERCED")
EARLY_climate <- subset(water_year,water_year$MONTH>=10)
EARLY_gdh<-split(as.numeric(EARLY_climate$GDH_MERCED),EARLY_climate$YEAR)
EARLY_precip <- split(as.numeric(EARLY_climate$PRECIP),EARLY_climate$YEAR)
LATE_climate <- subset(water_year,water_year$MONTH<=3)
LATE_gdh<-split(as.numeric(LATE_climate$GDH_MERCED),LATE_climate$YEAR)
LATE_precip <- split(as.numeric(LATE_climate$PRECIP),LATE_climate$YEAR)
tab_sum <- data.frame("Accumulated_GDH"=NULL,"Accumulated_PRECIP"=NULL,"SEASON"=NULL)
for(i in 1:length(EARLY_gdh)){
    sum_gdh <- sum(EARLY_gdh[[i]])
    sum_precip <- sum(EARLY_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_GDH"]<-sum_gdh
    tab_sum[nrow(tab_sum),"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"EARLY"
    tab_sum[nrow(tab_sum),"YEAR"]<-as.numeric(names(EARLY_gdh[i]))
    sum_gdh <- sum(LATE_gdh[[i]])
    sum_precip <- sum(LATE_precip[[i]])
    tab_sum[nrow(tab_sum)+1,"Accumulated_GDH"]<-sum_gdh
    tab_sum[nrow(tab_sum),"Accumulated_PRECIP"]<-sum_precip
    tab_sum[nrow(tab_sum),"SEASON"]<-"LATE"
    tab_sum[nrow(tab_sum),"YEAR"]<-as.numeric(names(EARLY_gdh[i]))
}

#Regress climate variables by year for Past 20 Years.
SEASON_LIST<-list("LATE","EARLY")
CLIM_var_LIST<-list("Accumulated_GDH","Accumulated_PRECIP")
lm_list_total<-list()
lm_list_study<-list()
for (q in CLIM_var_LIST){
    for (i in SEASON_LIST){
        lm <- lm(tab_sum[tab_sum$SEASON==i & tab_sum$YEAR >=2016,][[q]]~tab_sum[tab_sum$SEASON==i & tab_sum$YEAR>=2016,][["YEAR"]])
        name <- paste(i,q,sep="")
        lm_list_study[[name]]<-lm
        }
    }

for (q in CLIM_var_LIST){
    for (i in SEASON_LIST){
        lm <- lm(tab_sum[tab_sum$SEASON==i,][[q]]~tab_sum[tab_sum$SEASON==i,][["YEAR"]])
        name <- paste(i,q,sep="")
        lm_list_total[[name]]<-lm
    }
}

#Create Plots of Climate Varibles by Year


#Split dataframe into early and late winter seasons
EARLY_CIMIS <- na.omit(CIMIS[CIMIS$MONTH >= 10,])
LATE_CIMIS <- na.omit(CIMIS[CIMIS$MONTH <= 3,])
#Run ANOVA of Merced to four other CIMIS stations.
CIMIS_LIST <- list(EARLY_CIMIS,LATE_CIMIS)
aov.cimis_list <- list()
for (i in seq(1,2)){
    if(any(CIMIS_LIST[[i]][,"MONTH"]>=10)){
        name = "Early_CIMIS"
    } else {
        name = "Late_CIMIS"
    }
    lm_CIMIS<-lm(TEMP ~ SITE,data=CIMIS_LIST[[i]])
    anova <- aov(lm_CIMIS)
    name1 <- paste(name,"_TEMP",sep="")
    aov.cimis_list[[name1]]<-anova
    lm_CIMIS <- lm(PRECIP ~ SITE,data=CIMIS_LIST[[i]])
    anova <- aov(lm_CIMIS)
    name1 <- paste(name,"_PRECIP",sep="")
    aov.cimis_list[[name1]]<-anova
}
#Perform Dunnett test with Merced CIMIS station as the control.
EARLY_TEMP <- DunnettTest(EARLY_CIMIS$TEMP,EARLY_CIMIS$GROUP)
EARLY_PRECIP <- DunnettTest(EARLY_CIMIS$PRECIP,EARLY_CIMIS$GROUP)
LATE_TEMP <- DunnettTest(LATE_CIMIS$PRECIP,LATE_CIMIS$GROUP)
LATE_PRECIP <- DunnettTest(LATE_CIMIS$PRECIP,LATE_CIMIS$GROUP)
DUNLIST <- list(EARLY_TEMP,EARLY_PRECIP,LATE_PRECIP,LATE_TEMP)
for (i in DUNLIST){
    print(i)
}
#Create tables

#Create Plots

##########
###DONE###
##########

#Questions 1: Is flowering time of meadowfoam and whitetip clover significantly different between year and pool?

#Create Histograms & QQ PLOT
destination='Histograms_Phenology.pdf'
pdf(file=destination)
par(mfrow=c(2,4))
for (i in names(var_pheno)){
    name <- paste(i,"Histogram",sep="")
    plotNormalHistogram(var_pheno[,i],main=name)
}
destination='QQPlot_Phenology.pdf'
pdf(file=destination)
par(mfrow=c(2,4))
for (i in colnames(var_pheno)){
    name <- paste(i,"QQPlot",sep="_")
    qqnorm(var_pheno[,i],pch=1,frame=FALSE,main=name)
    qqline(var_pheno[,i],col="red",lwd=2)
}
#Determine if variables are normally distributed, and BOXCOX transform
#variables that are not normally distributed
var_pheno <- as.matrix(Pheno_data[,3:10],header=T)
AD_list <- list()
var_pheno1<-as.data.frame(var_pheno)
for (i in colnames(var_pheno)){
    name = paste(i,"_normality_test",sep="")
    ad_test <- ad.test(var_pheno[,i])
    if(ad_test$p.value < 0.05){
        BOXCOX<-transformTukey(var_pheno[,i])
        var_pheno1[,ncol(var_pheno1)+1]<-BOXCOX
        colnames(var_pheno1)[ncol(var_pheno1)]<- paste("BOX_COX",i,sep=".")
    }
    AD_list[[name]]<-ad_test
}
#Of the variables that are non-normally distributed, do transformed variables produce a lower AIC score for polynomial regressions?
log_group<-grepl("BOX_COX.",colnames(var_pheno1))
log_data <- var_pheno1[log_group]
len <- length(log_data)
lm_list <- list()
lm_list_Untransformed <- list()
x <- sub("BOX_COX.","",colnames(var_pheno1[log_group]))
for (i in seq(1,len)){
    for (q in seq(1,4)){
    name <- paste(colnames(log_data[i]),q,sep="_")
    lm <- lm(log_data[,i]~poly(Pheno_data$Year,q,raw=TRUE))
    lm_list[[name]] <- lm
    }
}
AIC_Scores <- data.frame(BOX_COX=NULL)
for (i in lm_list){
    aic <- AIC(i)
    AIC_Scores[nrow(AIC_Scores)+1,"BOX_COX"] <- aic
}
for (i in seq(1:len)){
    for (q in seq(1,4)){
        name <- paste(x[i],q,sep=".Poly_")
        lm <- lm(as.matrix(var_pheno1[x[i]])~poly(Pheno_data$Year,q,raw=TRUE))
        lm_list_Untransformed[[name]]<-lm
    }
}
AIC_Scores_UnTransformed<-data.frame(UnTransformed=NULL)
for (i in lm_list_Untransformed){
    aic <- AIC(i)
    AIC_Scores_UnTransformed[nrow(AIC_Scores_UnTransformed)+1,"UnTransformed"]<-aic
}
AIC_Scores <- cbind("UnTransformed"=AIC_Scores_UnTransformed$UnTransformed,AIC_Scores)

for(i in 1:nrow(AIC_Scores)){
    if(AIC_Scores$UnTransformed[i]<AIC_Scores$BOX_COX[i]){
        AIC_Scores$Winner <- print("UnTransformed")
    } else {
        AIC_Scores$Winner <- print("Transformed")
    }
}
row.names(AIC_Scores)<-names(lm_list_Untransformed)

#Which polynomial regression best fits the data?
list_Unshuff<-list()
for (i in names(var_pheno1)){
    Var_SHUFF<- sample(nrow(var_pheno1[i]))
    var_pheno1[,ncol(var_pheno1)+1]<-Var_SHUFF
    colnames(var_pheno1)[ncol(var_pheno1)]<-paste(i,"Shuffled",sep="_")
}
for (i in names(var_pheno1)){
    a<-as.list(var_pheno1[i])
    a <- as.data.frame(a)
    a<-cbind(Year=Pheno_data$Year,a)
    list_Unshuff[[i]]<-a
}
v_length <- length(var_pheno1)
h <- v_length/2
list_Shuff <-list_Unshuff[(h+1):v_length]
list_Unshuff<-list_Unshuff[0:h]
degree = 4
K = 10
k_list=list()
regress_list=list()
Breusch_Pagan_list=list()
BP_table<-data.frame(Heteroscedastic=NULL)
mse=matrix(data=NA,nrow=K,ncol=degree)
for (x in list_Shuff){
    for(i in 1:K){
    folds <- cut(seq(1,nrow(x)),breaks=K,labels=FALSE)
    index <- which(folds==i,arr.ind=TRUE)
    test <- x[index,]
    train <- x[-index,]

    for(j in 1:degree){
        fittrain = lm(train[,2] ~ poly(train[,1],j,raw=TRUE),train)
        fittest = predict(fittrain, new = test)
        mse[i,j]=mean((fittest-as.numeric(test[,2]))^2)
     }
    }
    name <- names(x[2])
    mse_5 <- na.omit(mse)
    min_mse <- min(colMeans(mse_5))
    v <- colMeans(mse_5)
    k = which(v==min_mse)
    k_list[[name]]<-k
    #Run Polynomial regression using resulting best fit polynomial for
    #each zone
    for (y in list_Unshuff){
        q = names(y[2])
        d = names(x[2])
        d <- sub("_Shuffled","",d)
        if(q==d){
            y<-na.omit(y)
            regress <- lm(y[,2]~poly(y[,1],k,raw=TRUE),data = y)
            #Perform Breusch-Pagan test to determine if model residuals exhibit heteroscedasticity
            bpmodel <- bptest(regress)
            #Plot the Regression
            ggplot(y,aes(x=y[,1],y=y[,2]))+
            geom_point()+
            stat_smooth(method="lm", formula = y ~ poly(x,k,raw=TRUE),linewidth=1)+
            xlab('Year')+
            ylab('Phenological_Trait')
            #Save Regressions as pdf
            name=paste(d,"Polynomial_regression_VP_Phenology.pdf",sep="_")
            ggsave(
                filename=name,
                plot=last_plot(),
                device=pdf,
                dpi=300,
                limitsize=TRUE
            )
            }
        regress_list[[d]]<-regress
        Breusch_Pagan_list[[d]]<- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table[nrow(BP_table)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table[nrow(BP_table)+1,"Heteroscedastic"]<-print("False")
    }
}
row.names(BP_table)<-names(Breusch_Pagan_list)

##########
###DONE###
##########

#Question 2: What climate variables best explains meadowfoam and whitetip clover floral phenology?

#Run Stepwise regression of each phenological variable
lm_list <- list()
BP_table_PHENO<-data.frame("Heteroscedastic"=NULL)
Breusch-Pagan_list_PHENO<-list()
for(W in names(Pheno_data[,3:10])){
    Climate_List<-list(Pheno_data$GDDEARLY,Pheno_data$GDDLATE,Pheno_data$PRECIP_EARLY,Pheno_data$PRECIP_LATE)
    names(Climate_List)<-c("GDDEARLY","GDDLATE","PRECIP_EARLY","PRECIP_LATE")
    for(pp in seq(1,5)){
        a <- length(Climate_List)
        b <- a - 1
        c <- b - 1
        d <- c - 1
        ax <- names(Climate_List[a])
        bx <- names(Climate_List[b])
        cx <- names(Climate_List[c])
        dx <- names(Climate_List[d])
        if(a>=1 && b>=1 && c>=1 && d>= 1 ){
            lm <- lm(Pheno_data_matrix[,W] ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]]+Climate_List[[d]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx,dx)
            fail <- "NO"
        }else if(a>=1 && b>=1 && c>=1){
            lm <- lm(Pheno_data_matrix[,W] ~ Climate_List[[a]]+Climate_List[[b]]+Climate_List[[c]])
            names(lm$coefficients)<-c("Intercept",ax,bx,cx)
            fail <- "NO"
        }else if(a>=1 && b>=1){
            lm <- lm(Pheno_data_matrix[,W] ~ Climate_List[[a]]+Climate_List[[b]])
            names(lm$coefficients)<-c("Intercept",ax,bx)
            fail <- "NO"
        }else if(a>=1){
            lm <- lm(Pheno_data_matrix[,W] ~ Climate_List[[a]])
            names(lm$coefficients)<-c("Intercept",ax)
            fail <- "NO"
        }else {
            fail <- "YES"
        }
        temp <- data.frame(Variable=NULL,Value=NULL)
        for(i in 1:(length(Climate_List))+1){
            if(summary(lm)$coefficients[i,4]>0.05){
            Var<-(names(summary(lm)$coefficients[,4][i]))
            Val<-summary(lm)$coefficients[Var,4]
            temp[nrow(temp)+1,"Variable"]<-Var
            temp[nrow(temp),"Value"]<-Val
            } 
    }
        if(any(temp$Value>0.05)==TRUE){
        x <- temp[which.max(temp$Value),]
        x <- x[,1]
        Climate_List<-Climate_List[names(Climate_List) != x]
        }   
    }
    if(fail == "NO"){
    lm_list[[W]] <- lm
    bpmodel <- bptest(lm)
    Breusch-Pagan_list_PHENO[[W]] <- bpmodel
    }
    if(bpmodel$p.value<0.05){
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("True")
        }else   {
    BP_table_PHENO[nrow(BP_table_PHENO)+1,"Heteroscedastic"]<-print("False")
    }
}

##########
###DONE###
##########

#Question 3: Does distance from center of vernal pool influence the floral phenology and
#growth patterns of meadowfoam and whtietip clover?









