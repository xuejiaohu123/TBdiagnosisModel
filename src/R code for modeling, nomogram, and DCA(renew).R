
### Part 1: Select the best feature subsets of logistic regression model in the Selection Cohort using the bestglm function (AIC/BIC) and 10-fold cross-validation.
#The bestglm() function begins with a data frame containing explanatory variables and response variables. The response variable should be in the last column. 

# Read the data of the Selection Cohort ("Selection.csv")
Selection <- read.csv("Selection.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
Vars <- colnames(Selection )[c(2:10,16)]
Catvars <- colnames(Selection )[c(1,11:15)]
Selection[Vars] <- lapply(Selection[Vars], as.numeric)
Selection[Catvars] <- lapply(Selection[Catvars], as.factor)

library(bestglm)
set.seed(123)
# (1) for EHR+lncRNA model
Combined.for.bestglm <- within(Selection, {
  y   <- Group_value
  Group <- NULL
  Group_value  <- NULL        # Delete the columns of  Group and Group_value
})

library(bestglm)
res.bestglm.AIC <- bestglm(Xy = Combined.for.bestglm, family = binomial(link='logit'), IC = "AIC", TopModels = 10, method = "exhaustive")
res.bestglm.AIC$BestModels # Show top 10 models
res.bestglm.BIC <- bestglm(Xy = Combined.for.bestglm, family = binomial(link='logit'), IC = "BIC", TopModels = 10, method = "exhaustive")
res.bestglm.BIC$BestModels 

# (2) for EHR model
EHR.for.bestglm <- within(Selection[,5:16], {
  y   <- Group_value
  Group_value  <- NULL        # Delete the columns of  Group and Group_value
})
res.bestglm.AIC <- bestglm(Xy = EHR.for.bestglm, family = binomial(link='logit'), IC = "AIC", TopModels = 10, method = "exhaustive")
res.bestglm.AIC$BestModels 
res.bestglm.BIC <- bestglm(Xy = EHR.for.bestglm, family = binomial(link='logit'), IC = "BIC", TopModels = 10, method = "exhaustive")
res.bestglm.BIC$BestModels

# (3) for LncRNA model
Lnc.for.bestglm <- within(Selection[,c(2:4,16)], {
  y   <- Group_value
  Group_value  <- NULL  
})
res.bestglm.AIC <- bestglm(Xy = Lnc.for.bestglm, family = binomial(link='logit'), IC = "AIC", TopModels = 10, method = "exhaustive")
res.bestglm.AIC$BestModels 
res.bestglm.BIC <- bestglm(Xy = Combined.for.bestglm, family = binomial(link='logit'), IC = "BIC", TopModels = 10, method = "exhaustive")
res.bestglm.BIC$BestModels 


# 10-fold cross-validation to further identify a prediction model，taking "EHR+lncRNA" modeling as an example
library(caret)
set.seed(1)
ctrl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)  #method = "repeatedcv"
grid <- expand.grid(.fL=c(0), .usekernel=c(FALSE))
model <- train(Group_value~ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,
               data=Selection, method="glm", trControl = ctrl)
model$results 

# Plot 10-fold AUC, mena AUC and AUC of the whole Selection Cohort
df <- model$pred
df <- as.list(df)
list(df)
split(df, 1:300)
predlist <- split(df$pred, df$Resample)
obslist <- split(df$obs, df$Resample)
library(cvAUC)
out <- cvAUC(predlist, obslist, folds = 10)
out

##Plot CV AUC
plot(out$perf, col="grey82", lty=1)

# plot(out$perf, col="red", avg="vertical", add=TRUE)

#Plot AUC using the entire Selection Cohort
Combined_model <- glm(Group~ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,
                      family=binomial(link="logit"),data=Selection)
Pre_com <- predict.glm(Combined_model,type='response')
Pred<-prediction(Pre_com,Combined_model$y)
perf<-performance(Pred,"tpr","fpr")
plot(perf, col="blue", add = TRUE)
abline(a = 0, b = 1, col="gray",lty=2)
legend(0.4,0.15,cex=0.7,legend=c("The entire training set", "10-fold cross-validation"),
       col=c("blue", "grey82"), bty='n', adj = c(0, 0.5), lwd = 1.5, seg.len = 0.8, text.width = 1)

# (1)Evaluate the models and select the best feature subsets based on the criteria of the minimum of AIC, and the maximum of performance (accuracy and AUC). 
# (2)The best EHR+lncRNA model is Combined_model <- glm(Group ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link="logit"),data=Selection);  
#The best EHR model is EHR_model <- glm(Group~Age+HB+M+Low_grade_fever+Weight_loss+CT_calcification+CT_bronchus_sign+TB_IGRA,family=binomial(link="logit"),data=Selection);
#The best LncRNA model is LncRNA_model <- glm(Group~ENST00000497872+n333737+n335265,family=binomial(link="logit"),data=Selection)


##################################################################################################################################
### Part2: Model evaluation and comparison: EHR+lncRNA model, EHR model, and LncRNA model

## Model evaluation 1: performance
Combined_model <- glm(Group~ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,
                      family=binomial(link="logit"),data=Selection)
EHR_model <- glm(Group~Age+HB+M+Low_grade_fever+Weight_loss+CT_calcification+CT_bronchus_sign+TB_IGRA,
                      family=binomial(link="logit"),data=Selection)
LncRNA_model <- glm(Group~ENST00000497872+n333737+n335265,family=binomial(link="logit"),data=Selection)

Pre_com <- predict.glm(Combined_model,type='response')
Pre_ehr <- predict.glm(EHR_model,type='response')
Pre_lnc <- predict.glm(LncRNA_model,type='response')

# calculate AUC
library(pROC)
roc_com <- roc(Selection$Group_value, Pre_com)
roc_com$auc 
ci(roc_com)
roc_ehr <- roc(Selection$Group_value, Pre_ehr)
roc_ehr$auc
ci(roc_ehr)
roc_lnc <- roc(Selection$Group_value, Pre_lnc)
roc_lnc$auc
ci(roc_lnc)
roc_com

# draw ROC curves
roc1<-roc(Selection$Group_value, Pre_com,
          plot=TRUE, 
          ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
          main="Statistical comparison", percent=TRUE, col="#1c61b6",
          show.thres=TRUE)
plot(roc_com,axes=F,legacy.axes=F,asp = F,lwd = 0.85,lty = 1, identity.lty=2, show.thresholds=TRUE,
     col='blue', main='Selection Cohort')
plot(roc_ehr,axes=F,legacy.axes=F,asp = F,lwd = 0.85,lty = 1, identity.lty=2,
     col='orange', add = TRUE)
plot(roc_lnc,axes=F,legacy.axes=F,asp = F,lwd = 0.85,lty = 1, identity.lty=2,
     col='green4', add = TRUE)
axis(side=1,at=seq(0,1,0.2),lwd=1,tcl=-0.2, mgp=c(4, 1.6, 1.5))
axis(side=2,at=seq(0,1,0.2),lwd=1,tcl=-0.2,mgp=c(4, 0.8, 0.4))
legend(0.7,0.3,cex=0.7,legend=c("EHR+lncRNA  AUC = 0.92 (0.89, 0.95)", "EHR only  AUC = 0.87 (0.83, 0.91)","lncRNA only  AUC = 0.82 (0.77, 0.86)"),
       col=c("blue", "orange","green4"), bty='n', adj = c(0, 0.5), lwd = 1.5, seg.len = 0.8)

# get an “optimal” cut point (A cutoff of each model was determined by combining the Youden’s index and the sensitivity for the samples in the training dataset equal to or greater than 0.85)
library(dplyr)
roc_com %>%
  coords(transpose=FALSE) %>%
  filter(sensitivity > 0.85)

# AUC comparison for EHR+lncRNA model, EHR model, and lncRNA model  p-value =0；05/3=0.016
roc.test(roc_com, roc_ehr) # default Delong method
roc.test(roc_com, roc_lnc) 
roc.test(roc_ehr, roc_lnc) 

## Model evaluation 2: model features(multicollinearity and coefficients), Hosmer-Lemeshow goodness of fit test, Likelihood ratio test,Nagelkerke R2, and McFadden R2 
summary(Combined_model)
varImp(Combined_model) # Calculating feature importance
sort(vif(Combined_model),decreasing = TRUE) # VIF > 10 means multicollinearity exists among features, here VIF ranged from 1.04-1.29
pscl::pR2(Combined_model) 
ResourceSelection::hoslem.test (Combined_model$y, fitted(Combined_model), g=10)  # The result is X-squared = 2.6079, df = 8, p-value = 0.9565,indicating no evidence of poor fit.


##################################################################################################################################
###Part 3: Nomogram visulization for the optimal model (EHR+lncRNA model)
library(rms)
dd <- datadist(Selection)
options(datadist = 'dd')
Combined <- lrm(Group ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA, x = TRUE,y = TRUE, data = Selection)
print(Combined)
nom <- nomogram(Combined, fun= plogis,
                fun.at = c(.01,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95),
                lp=F, funlabel="Risk of NPTB")
plot(nom,xfrac=.35)

#Calibration curve of the nomogram in the Selection cohort
cal<-calibrate(Combined,method="boot",B=500)
plot(cal)


##################################################################################################################################
###Part 4: Validate the Nomogram (EHR+lncRNA model) in the independent Validaiotn cohort

# (1) Performance of the Validation Cohort (read "Validation.csv")
Validation <- read.csv("Validation.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
Vars <- colnames(Validation)[c(2:10,16)]
Catvars <- colnames(Validation)[c(1,11:15)]
Validation[Vars] <- lapply(Validation[Vars], as.numeric)
Validation[Catvars] <- lapply(Validation[Catvars], as.factor)
str(Validation)

Pre_comv <- predict.glm(Combined_model,type='response',newdata=Validation)
Pre_ehrv <- predict.glm(EHR_model,type='response',newdata=Validation)
Pre_lncv <- predict.glm(LncRNA_model,type='response',newdata=Validation)

# calculate sensitivity,specificity, accuracy, PPV, NPV
# The cutoff probability in the Selection Cohort was 0.37 for "EHR+lncRNA" model, 0.26 for "EHR only" model, and 0.32 for "lncRNA" model, respectively. 
Pre2 =ifelse(Pre_comv>0.37,1,0)
Validation$Pre2 = Pre2
true_value2=Validation$Group_value
predict_value2=Validation$Pre2
cft2 <- table(predict_value2,true_value2)
confusionMatrix(cft2, positive = "1")

Pre3 =ifelse(Pre_ehrv>0.26,1,0)
Validation$Pre3 = Pre3
true_value3=Validation$Group_value
predict_value3=Validation$Pre3
cft3 <- table(predict_value3,true_value3)
confusionMatrix(cft3, positive = "1")

Pre4 =ifelse(Pre_lncv>0.32,1,0)
Validation$Pre4 = Pre4
true_value4=Validation$Group_value
predict_value4=Validation$Pre4
cft4 <- table(predict_value4,true_value4)
confusionMatrix(cft4, positive = "1")

# calculate AUC
roc_comv <- roc(Validation$Group_value, Pre_comv)
roc_comv$auc 
ci(roc_comv)
roc_ehrv <- roc(Validation$Group_value, Pre_ehrv)
roc_ehrv$auc
ci(roc_ehrv)
roc_lncv <- roc(Validation$Group_value, Pre_lncv)
roc_lncv$auc
ci(roc_lncv)

# AUC comparison for EHR+lncRNA model, EHR model, and lncRNA model  p-value =0；05/3=0.016
roc.test(roc_comv, roc_ehrv)
roc.test(roc_comv, roc_lncv) 
roc.test(roc_ehrv, roc_lncv) 

# draw ROC curves for Validation Cohort
plot(roc_comv,axes=F,legacy.axes=F,asp = F,lwd = 0.85, lty = 1, identity.lty=2,
     col='blue', main='Validation Cohort')
plot(roc_ehrv,axes=F,legacy.axes=F,asp = F,lwd = 0.85, lty = 1, identity.lty=2,
     col='orange', add = TRUE)
plot(roc_lncv,axes=F,legacy.axes=F,asp = F,lwd = 0.85, lty = 1, identity.lty=2,
     col='green4', add = TRUE)
axis(side=1,at=seq(0,1,0.2),lwd=1,tcl=-0.2, mgp=c(4, 1.6, 1.5))
axis(side=2,at=seq(0,1,0.2),lwd=1,tcl=-0.2,mgp=c(4, 0.8, 0.4))
legend(0.7,0.3,cex=0.7,legend=c("EHR+lncRNA  AUC = 0.89 (0.84, 0.93)", "EHR only  AUC = 0.83 (0.78, 0.89)","lncRNA only  AUC = 0.80 (0.74, 0.85)"),
       col=c("blue", "orange","green4"), bty='n', adj = c(0, 0.5), lwd=1.5, seg.len = 0.8)


# (2) calibration characteristics (calibration plot and Hosmer-Lemeshow test)
Nomo.validation <- lrm(Group ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA, x = TRUE,y = TRUE, data=Validation)
plot(calibrate(Nomo.validation,method="boot",B=500))

Combined.validation <- glm(Group ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link="logit"), data=Validation)
ResourceSelection::hoslem.test (Combined.validation$y, fitted(Combined.validation), g=10) 


##################################################################################################################################
#Part 5: Nomogram DCA analysis http://mdbrown.github.io/rmda/ for the nomogram vs the coventional EHR model

#Read the data of the Allpatients cohort ("Allpatients.csv")
Allpatients <- read.csv("Allpatients.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
Vars <- colnames(Allpatients)[c(2:10,16)]
Catvars <- colnames(Allpatients)[c(1,11:15)]
Allpatients[Vars] <- lapply(Allpatients[Vars], as.numeric)
Allpatients[Catvars] <- lapply(Allpatients[Catvars], as.factor)
str(Allpatients)

#DCA http://mdbrown.github.io/rmda/
library(rmda)
set.seed(123) 
EHR <- decision_curve(Group_value ~ Age+HB+M+Low_grade_fever+Weight_loss+CT_calcification+CT_bronchus_sign+TB_IGRA,family=binomial(link = "logit"),data = Allpatients, policy = "opt-in", thresholds= seq(0,1, by = 0.01), confidence.intervals =0.95, bootstraps = 500)
Nomogram <- decision_curve(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link = "logit"),data = Allpatients, policy = "opt-in", thresholds= seq(0,1, by = 0.01), confidence.intervals =0.95, bootstraps = 500)
plot_decision_curve( list(EHR, Nomogram),curve.names = c("EHR only model", "Nomogram"),xlab=c("Threshold Probability"),standardize = F,col = c("red", "blue"), confidence.intervals =FALSE, xlim = c(0, 1), lty = c(2,1),cost.benefit.axis = FALSE,legend.position = "none")
legend("topright", cex=0.7, legend=c("EHR only model", "Nomogram","All","None"),col=c('blue','red','grey','black'), lty = c(2,1,1,1),lwd=c(2, 2, 1, 1))
library(gridGraphics)
grid.echo()
grid.ls()
grid.remove("abline", grep=TRUE, global=TRUE)
