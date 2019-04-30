
#Part 1: Select the best feature subsets of logistic regression model in the Selection Cohort using the bestglm function (AIC/BIC).
#The bestglm() function begins with a data frame containing explanatory variables and response variables. The response variable should be in the last column. 

# Read the data of the Selection Cohort ("Selection.csv")
Selection <- read.csv("Selection.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
catVars <- c("Low_grade_fever","Weight_loss","CT_calcification","CT_bronchus_sign","TB_IGRA")
Selection[catVars] <- lapply(Selection[catVars], as.factor)
unique(Selection)
str(Selection)

library(bestglm)
set.seed(123)
# (1) for EHR+lncRNA model
Combined.for.bestglm <- within(Selection, {
  y   <- Group_value
  Group <- NULL
  Group_value  <- NULL        # Delete the columns of  Group and Group_value
})
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



# Part2: Model evaluation and selection: the best EHR+lncRNA model, EHR model, and LncRNA model
# (1)Evaluating the models and select the best model with the minimum of AIC, and the maximum of performance (the most important C-statistics, accuracy, sensitivity,specificity, PPV, NPV). 
# (2)The best EHR+lncRNA model is Combined_model <- glm(Group ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link="logit"),data=Selection); The best EHR model is EHR_model <- glm(Group~Age+HB+M+Low_grade_fever+Weight_loss+CT_calcification+CT_bronchus_sign+TB_IGRA,family=binomial(link="logit"),data=Selection); The best LncRNA model is LncRNA_model <- glm(Group~ENST00000497872+n333737+n335265,family=binomial(link="logit"),data=Selection)

# Take "EHR+lncRNA" modeling as an example
# Model evaluation 1: performance
Combined_model <- glm(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link="logit"),data=Selection)
Pre <- predict.glm(Combined_model,type='response')
Pre =ifelse(Pre>0.5,1,0)
Selection$Predict1 = Pre
true_value=Selection$Group_value
predict_value=Selection$Predict1
cft1 <- table(true_value,predict_value)
library(caret)
confusionMatrix(cft1, positive = "1")

# Model evaluation 2: model features, coefficients of features, C-statistics,bias_corrected_c_statistic, Likelihood ratio test,Nagelkerke R2, and McFadden R2  
library(rms)
dd <- datadist(Selection)
options(datadist = 'dd')
Combined <- lrm(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA, x = TRUE,y = TRUE, data = Selection)
print(Combined)

pscl::pR2(Combined) 

v <- validate(Combined, dxy=TRUE, B=500)
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
orig_c_statistic <- abs(orig_Dxy)/2+0.5
bias_corrected_c_statistic  <- abs(Dxy)/2+0.5
print(bias_corrected_c_statistic) 
print(orig_c_statistic)

# Model evaluation 3: Hosmer-Lemeshow goodness of fit test 
library(ResourceSelection)
hoslem.test (Selection$Group_value, fitted(Combined_model), g=10)  # The result is X-squared = 2.6079, df = 8, p-value = 0.9565,indicating no evidence of poor fit.

# Model evaluation 4: Multicollinearity among features
library(car)  
sort(vif(Combined_model),decreasing = TRUE) # VIF > 10 means multicollinearity exists among features, here VIF ranged from 1.04-1.29



#Part 3: Nomogram visulization for the optimal model (EHR+lncRNA model)
nom <- nomogram(Combined, fun= plogis,
                fun.at = c(.01,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95),
                lp=F, funlabel="Risk of PTB without pathogenic evidence")
plot(nom,xfrac=.35)
#Calibration curve of the nomogram in the selection cohort
cal<-calibrate(Combined,method="boot",B=500)
plot(cal)
# Calculating feature importance
varImp(Combined_model)



#Part 4: Validate the Nomogram (EHR+lncRNA model) in the independent Validaiotn cohort

# Read the data of the Validation Cohort ("Validation.csv")
Validation <- read.csv("Validation.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
Validation [catVars] <- lapply(Validation[catVars], as.factor)
str(Validation)

# Take "EHR+lncRNA" modeling as an example
# (1) performance
Pre2 <- predict.glm(Combined_model,type='response',newdata=Validation)
Pre2 =ifelse(Pre2>0.5,1,0)
Validation$Pre2 = Pre2
true_value2=Validation$Group_value
predict_value2=Validation$Pre2
cft2 <- table(true_value2,predict_value2)
confusionMatrix(cft2, positive = "1")
# (2) calibration characteristics (calibration plot and Hosmer-Lemeshow test)
Nomo.validation <- lrm(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA, x = TRUE,y = TRUE, data=Validation)
plot(calibrate(Nomo.validation,method="boot",B=500))

Combined.validation <- glm(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link="logit"), data=Validation)
hoslem.test (Validation$Group_value, fitted(Combined.validation), g=10) 

# (3) bias_corrected c_statistic
v <- validate(Nomo.validation, dxy=TRUE, B=500)
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"]
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"]
orig_c_statistic <- abs(orig_Dxy)/2+0.5
bias_corrected_c_statistic  <- abs(Dxy)/2+0.5
print(bias_corrected_c_statistic) 
print(orig_c_statistic)



#Part 5: Nomogram DCA analysis http://mdbrown.github.io/rmda/ for the nomogram vs the coventional EHR model

# Read the data of the Allpatients cohort ("Allpatients.csv")
Allpatients <- read.csv("Allpatients.csv",as.is = T,sep = ",",header = TRUE, na.strings = c(""))
Allpatients[catVars] <- lapply(Allpatients[catVars], as.factor)
str(Allpatients)

library(rmda)
set.seed(123) 
EHR <- decision_curve(Group_value ~ Age+HB+M+Low_grade_fever+Weight_loss+CT_calcification+CT_bronchus_sign+TB_IGRA,family=binomial(link = "logit"),data = Allpatients, policy = "opt-in", thresholds= seq(0,1, by = 0.01), confidence.intervals =0.95, bootstraps = 500)
Nomogram <- decision_curve(Group_value ~ ENST00000497872+n333737+n335265+Age+HB+Low_grade_fever+Weight_loss+CT_calcification+TB_IGRA,family=binomial(link = "logit"),data = Allpatients, policy = "opt-in", thresholds= seq(0,1, by = 0.01), confidence.intervals =0.95, bootstraps = 500)
plot_decision_curve( list(EHR, Nomogram),curve.names = c("EHR only model", "Nomogram"),xlab=c("Threshold Probability"),standardize = F,col = c("#F8766D","#619CFF"), confidence.intervals =FALSE, xlim = c(0, 1), lty = c(2,1),cost.benefit.axis = FALSE,legend.position = "none")
legend("topright", cex=0.7, legend=c("EHR only model", "Nomogram","All","None"),col=c("#F8766D","#619CFF",'grey','black'), lty = c(2,1,1,1),lwd=c(2, 2, 1, 1))
library(gridGraphics)
grid.echo()
grid.ls()
grid.remove("abline", grep=TRUE, global=TRUE)

