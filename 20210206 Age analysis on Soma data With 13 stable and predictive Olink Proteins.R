####################################################################
# 20210206 Age analysis on SomaLogic data With 13 stable and predictive Olink Proteins
# Date:      10-01-2021
# Author:    Adama J.Berndt
####################################################################
# Purpose: To repeat final iteration of Olink analysis using 
#          only the set of overlapping Soma/OLink proteins that were predictive in the olink dataset.
####################################################################

# Load settings
options(stringsAsFactors = FALSE)

# Load packages
library(dplyr)
library(gplots)
library(GGally)
library(ggplot2)
library(dendextend)
library(glmnet)
library(ggpubr)
library(caret)
library(leaps)
library(stringr)



# Directory
getwd()
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/SomaLogic")

####################################################################
# Data formatting
####################################################################

# Load data
soma.dat <- read.csv("Soma compiled data with PC1,2 5-28-20.csv")
penn.age <- read.csv("Soma_Penn_AgeWithDecimal.csv")
names(penn.age)[1] <- "PatientID"
penn.age$PatientID <- as.character(penn.age$PatientID)

# For Penn patients, replace age with appropriate values (i.e. with
# 1 decimal place to match PDBP format and improve accuracy).
for (i in 1:nrow(Olink.age)) {
  soma.dat$Age[soma.dat$PatientID == Olink.age$PatientID[i]] <- 
    Olink.age$AGE..YRS[i]
}

# Rename 1st column
names(soma.dat)[1] <- "INDDID"

# Round to 1st decimal place
soma.dat$Age <- round(soma.dat$Age, digits = 1)

# Load OLink data
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/OLink")
olink.dat <- read.csv("OLink data 20210120.csv")

####################################################################
# Understand overlapping proteins
####################################################################

# Load Somalogic identifiers
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/SomaLogic")
soma.ids <- read.csv("IDs for shared (UPenn, PDBP) somalogic proteins.csv")

# Load list of OLink Proteins
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/OLink")
olink.ids <- read.csv("IDs for OLink proteins.csv")
names(olink.ids)[1] <- "Panels"
summary(as.factor(olink.ids$Panels))

# Remove panels used in only 20 samples
panels.to.remove <- c("Olink Target 96 Neuro Exploratory(v.3911)",
                      "Olink Target 96 Oncology II(v.7004)", 
                      "Olink Target 96 Immuno-Oncology(v.3111)", 
                      "Olink Target 96 Immune Response(v.3203)")
panel.idx <- rep(TRUE,1288)
for (i in 1:length(panels.to.remove)) {
  panel.idx[olink.ids$Panels == panels.to.remove[i]] <- FALSE
}
sum(panel.idx)
# Remove
olink.ids <- olink.ids[panel.idx,]

# Distinct OLink Assays
dim(olink.ids %>% select(Assay) %>% distinct())   # 900 (20 repeated)

# Generate soma list with OLink matches
soma.match.list <- soma.ids
soma.match.list$Match <- rep(FALSE, nrow(soma.match.list))
soma.match.list$OLink.Assay <- rep('',nrow(soma.match.list))
soma.match.list$OLink.Panel <- rep('',nrow(soma.match.list))
soma.match.list$OLink.Uniprot <- rep('',nrow(soma.match.list))
soma.match.list$OLink.idx <- rep('',nrow(soma.match.list))

# Run loop to find matches
for (i in 1:nrow(soma.ids)) {
  for (j in 1:nrow(olink.ids)) {
    # Any match?
    if (soma.olink.match(soma.ids$Uniprot[i],olink.ids$Uniprot.ID[j])) {
      # Assign TRUE
      soma.match.list$Match[i] <- soma.olink.match(soma.ids$Uniprot[i],
                                                   olink.ids$Uniprot.ID[j])
      # If fields are empty (no prior matches), just fill in.
      if (soma.match.list$OLink.Assay[i] == '') {
        soma.match.list$OLink.Assay[i] <- olink.ids$Assay[j]
        soma.match.list$OLink.Panel[i] <- olink.ids$Panels[j]
        soma.match.list$OLink.Uniprot[i] <- olink.ids$Uniprot.ID[j]
        soma.match.list$OLink.idx[i] <- as.character(j)
      } else {
        # Otherwise, do some gymnastics
        soma.match.list$OLink.Assay[i] <- paste(soma.match.list$OLink.Assay[i],
                                                olink.ids$Assay[j], sep = ',')
        soma.match.list$OLink.Panel[i] <- paste(soma.match.list$OLink.Panel[i],
                                                olink.ids$Panels[j], sep = ',')
        soma.match.list$OLink.Uniprot[i] <- paste(soma.match.list$OLink.Uniprot[i],
                                                  olink.ids$Uniprot.ID[j], sep = ',')
        soma.match.list$OLink.idx[i] <- paste(soma.match.list$OLink.idx[i],
                                                  as.character(j), sep = ',')
      }
    }
  }
}

# Summary
sum(soma.match.list$Match) # 213/940 Somalogic assays that have an OLink counterpart

# Save
setwd("C:/Users/Maria/Desktop/CP Lab/2021/Aging Project")
write.csv(soma.match.list[soma.match.list$Match,],
          file = "Somamers with matching OLink assays (2-6-21).csv",
          row.names = TRUE)

# Which somamers have >1 OLink assay?
sum(grepl(',', soma.match.list$OLink.idx)) # 12 have more than 1 assay
View(soma.match.list[grepl(',', soma.match.list$OLink.idx),]) # 2 assays/somamer (not >2)
soma.match.list <- soma.match.list[soma.match.list$Match,]

# Therefore, at most, you will have 357 + 12 = 369 correlations

####################################################################
# Look at repeated proteins
####################################################################

# Repeated proteins
list <- olink.ids %>% select(Assay) %>% distinct()
list$Sum <- numeric(nrow(list))
for (i in 1:nrow(list)) {
  list$Sum[i] <- sum(olink.ids$Assay == list$Assay[i])
}
list <- list[list$Sum > 1,]
list$cor <- numeric(nrow(list))
list$Assay <- gsub(pattern = '-', replacement = '.', list$Assay)
olink.dat1 <- olink.dat[,12:1299]
olink.dat1 <- olink.dat1[,panel.idx]
for (i in 1:nrow(list)) {
  idx1 <- c(which(grepl(list$Assay[i], names(olink.dat1))))[1]
  idx2 <- c(which(grepl(list$Assay[i], names(olink.dat1))))[2]
  list$cor[i] <- cor.test(olink.dat1[,idx1], olink.dat1[,idx2],
                          method = 'spearman', exact = FALSE)$estimate
}
View(list) # This is alarming

# Save
write.csv(list, file = "Repeated protein correlations 2-6-2021.csv",
          row.names = FALSE)

# NOTE: Repeated proteins DO NOT have a correlation of 1
hist(list$cor, main = 'Spearman correlations for
     duplicate OLink proteins', xlab = 'rho')
summary(list$cor)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1585  0.7764  0.9035  0.7892  0.9420  0.9852 

# beta-NGF, CXCL1, ENTPD6 are problematic

####################################################################
# Correlations between overlapping proteins
####################################################################

# Overlapping samples
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/SomaLogic")
overlapping.samples <- read.csv("Soma-Olink-Overlapping-Samples.csv")
names(overlapping.samples)[1] <- "INDDID"

# Generate dataframe of interest
corr.df <- data.frame(Soma.Target = character(369),
                      Soma.Uniprot = character(369),
                      OLink.Assay = character(369),
                      OLink.Panel = character(369),
                      OLink.Uniprot = character(369),
                      rho = numeric(369),
                      pval = numeric(369),
                      FDR = numeric(369))

# Lets generate data1, data2 for Soma & olink with matching INDDIDs for ease
data1 <- soma.dat[soma.dat$INDDID %in% overlapping.samples$INDDID,]
data2 <- olink.dat[olink.dat$INDDID %in% overlapping.samples$INDDID,]
data1 <- data1 %>% arrange(INDDID)
data2 <- data2 %>% arrange(INDDID)
all(data1$INDDID == data2$INDDID) # TRUE, proceed.
data2 <- data2[,12:ncol(data2)][,panel.idx]

# Generate loop
# Start count that wil be modified if >1 match
j <- 0
for (i in 1:nrow(soma.match.list)) {
  # Execute regardless
  # Fill table
  corr.df$Soma.Target[i+j] <- soma.match.list$Target[i]
  corr.df$Soma.Uniprot[i+j] <- soma.match.list$Uniprot[i]
  corr.df$OLink.Assay[i+j] <- strsplit(soma.match.list$OLink.Assay[i], split = ',')[[1]][1]
  corr.df$OLink.Panel[i+j] <- strsplit(soma.match.list$OLink.Panel[i], split = ',')[[1]][1]
  corr.df$OLink.Uniprot[i+j] <- strsplit(soma.match.list$OLink.Uniprot[i], split = ',')[[1]][1]
  # Indeces for correlation
  # idx1 --> Where names data1 match Target
  idx1 <- which(names(data1) == soma.match.list$Target[i])
  # idx2 --> Saved in soma.match.list
  idx2 <- as.numeric(strsplit(soma.match.list$OLink.idx[i], split = ',')[[1]][1])
  # Correlation Test
  corr.df$rho[i+j] <- cor.test(data1[,idx1], data2[,idx2], 
                               method = 'spearman', exact = FALSE)$estimate
  corr.df$pval[i+j] <- cor.test(data1[,idx1], data2[,idx2], 
                               method = 'spearman', exact = FALSE)$p.value
  # If there is MORE than 1 match...
  if (grepl(',', soma.match.list$OLink.idx[i])) {
    # Fill table
    j <- j+1
    corr.df$Soma.Target[i+j] <- soma.match.list$Target[i]
    corr.df$Soma.Uniprot[i+j] <- soma.match.list$Uniprot[i]
    corr.df$OLink.Assay[i+j] <- strsplit(soma.match.list$OLink.Assay[i], split = ',')[[1]][2]
    corr.df$OLink.Panel[i+j] <- strsplit(soma.match.list$OLink.Panel[i], split = ',')[[1]][2]
    corr.df$OLink.Uniprot[i+j] <- strsplit(soma.match.list$OLink.Uniprot[i], split = ',')[[1]][2]
    # Indeces for correlation
    # idx1 --> Where names data1 match Target
    idx1 <- which(names(data1) == soma.match.list$Target[i])
    # idx2 --> Saved in soma.match.list
    idx2 <- as.numeric(strsplit(soma.match.list$OLink.idx[i], split = ',')[[1]][2])
    # Correlation Test
    corr.df$rho[i+j] <- cor.test(data1[,idx1], data2[,idx2], 
                                 method = 'spearman', exact = FALSE)$estimate
    corr.df$pval[i+j] <- cor.test(data1[,idx1], data2[,idx2], 
                                  method = 'spearman', exact = FALSE)$p.value
  }
}
corr.df$FDR <- p.adjust(corr.df$pval, method = 'BH')

# Summary
summary(corr.df$rho)     
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.5247  0.1922  0.5234  0.4533  0.7481  0.9857 

# All proteins
sum(corr.df$pval < 0.05) # 217 nominally correlated proteins
sum(corr.df$FDR < 0.05)  # 206 correlated at FDR < 0.05
hist(corr.df$rho, main = "Spearman's Correlation for
     Overlapping Proteins",
     xlab = 'Rho')

# Distinct
dim(corr.df[corr.df$pval < 0.05,] %>% select(Soma.Target) %>% distinct()) # 213
dim(corr.df[corr.df$FDR < 0.05,] %>% select(Soma.Target) %>% distinct())  # 203

# Save corr.df
getwd()
setwd("C:/Users/Maria/Desktop/CP Lab/2021/Aging Project")
write.csv(corr.df, file = "Soma-Olink Correlations 2-6-21.csv",
          row.names = FALSE)

####################################################################
# Figure for 9 duplicate assays
####################################################################

# Make data frame
df1 <- data.frame(Assay = repeated,
                  Rho1 = numeric(9),
                  Rho2 = numeric(9))2


corr.df1 <- corr.df[corr.df$OLink.Assay %in% repeated,]
corr.df1 <- corr.df1 %>% arrange(OLink.Assay)
df1 <- df1 %>% arrange(Assay)
df1$Rho1 <- corr.df1$rho[seq(from = 1, to = 17, by = 2)]
df1$Rho2 <- corr.df1$rho[seq(from = 2, to = 18, by = 2)]
rm(corr.df1)

# Plot graph
ggpaired(df1, cond1 = 'Rho1', cond2 = 'Rho2', fill = 'condition',
         ylab = "Spearman's Rho", )

####################################################################
# Using Top15 predictive proteins from Olink to train on Soma-mer data
####################################################################

# Select targets

Top13OlinkID <- c("CCDC80", "EDA2R", "IL.17D", "FUT3.FUT5", "IL.27", "LHB", "CST5", "ANGPTL3", "EGFR", "CDON", "ACAN", "CLEC11A.1", "CD200R1")

#Used the following command to check all the Olink identifiers for the same proteins that do not match the somalogic dataset, and then searched through soma.match,list to find the correct proteins
#print(Top13 %in% soma.ids$Target)
#[1] FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
Top15SomaID <- c("URB", "XEDAR", "IL.17D", "Fucosyltransferase.3","FUT5", "IL.27", "Luteinizing.hormone", "CYTD", "ANGL3", "ERBB1", "CDON", "Aggrecan", "SCGF.beta","SCGF.alpha", "MO2R1") #Fut3 and Fut5 have independent rows in somalogic dataset but only one rowname in Olink, top13 proteins becomes a top14, same for CLEC11A which becomes SCGF.beta and SCGF.alpha



# Select data
ctl.dat <- cbind(soma.dat[soma.dat$Group == 'Control',c(1:6,947:949)], 
                 soma.dat[soma.dat$Group == 'Control',Top15SomaID])
ND.dat <- cbind(soma.dat[soma.dat$Group != 'Control',c(1:6,947:949)], 
                soma.dat[soma.dat$Group != 'Control',Top15SomaID])

####################################################################

# 1) IDENTIFY AGE-ASSOCIATED PROTEINS:

# Generate data frame
modelOlink15.df <- data.frame(Protein = character(14), Uniprot = character(14),
                        Age.coef = numeric(14), Age.pval = numeric(14),
                        Age.FDR = numeric(14),
                        Male.coef = numeric(14), Male.pval = numeric(14),
                        Male.FDR = numeric(14),
                        PC1.coef = numeric(14), PC1.pval = numeric(14),
                        PC1.FDR = numeric(14)
                        )

# Run loop and apply linear model to all proteins
for (i in 1:14) {
  modelOlink15.df$Protein[i] <- names(ctl.dat)[i+9]
  formula_i <- paste(modelOlink15.df$Protein[i], "~Age+Sex+PC1", sep = "")
  model_i <- lm(formula_i, ctl.dat)
  modelOlink15.df$Age.coef[i] <- summary(model_i)$coefficients[2,1]
  modelOlink15.df$Age.pval[i] <- summary(model_i)$coefficients[2,4]
  modelOlink15.df$Male.coef[i] <- summary(model_i)$coefficients[3,1]
  modelOlink15.df$Male.pval[i] <- summary(model_i)$coefficients[3,4]
}

# Adjust p-values
modelOlink15.df$Age.FDR <- p.adjust(modelOlink15.df$Age.pval, method = 'BH')
modelOlink15.df$Male.FDR <- p.adjust(modelOlink15.df$Male.pval, method = 'BH')
modelOlink15.df$PC1.FDR <- p.adjust(modelOlink15.df$PC1.pval, method = 'BH')

# Summary
sum(modelOlink15.df$Age.FDR < 0.05) # 100 age-associated proteins
modelOlink15.df <- modelOlink15.df %>% arrange(Age.FDR)
View(modelOlink15.df)
setwd("C:/Users/nezne/Desktop/Aging Project/Aging Project/SomaLogic/Analysis with Top15 proteins in Olink dataset")
write.csv(modelOlink15.df, 
          file = "Age Associated Linear model results for top 15 Olink proteins with Somalogic values 10-01-2021.csv",
          row.names = FALSE)


sig.proteins <- modelOlink15.df$Protein[modelOlink15.df$Age.FDR < 0.05]
sum(sig.proteins %in% modelOlink15.df$Protein)
####################################################################

# 2) MODEL SELECTION WITH ONLY 6 TOP PROTEINS:

# Will run elastic net with 10-fold CV repeated times with different seeds

# First, using 100% samples:
idx.sig <- match(sig.proteins, names(ctl.dat))
Y <- ctl.dat$Age
X <- model.matrix(Age~., data = ctl.dat[,c(4,5,7,idx.sig)])[, -1]  # proteins,sex,PC1
count.matrix <- data.frame(Param = colnames(X), Non.Zero = numeric(8))
for (i in 1:1000) {
  set.seed(i)
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix$Param)
  count.matrix$Non.Zero[idx_i] <- count.matrix$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (100% samples, 6 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")
# Summary
sum(count.matrix$Non.Zero == 1000) # 7
sum(count.matrix$Non.Zero >= 900)  # 8

# Second, using 90% samples:
count.matrix_0.9 <- data.frame(Param = colnames(X), Non.Zero = numeric(8))
for (i in 1:1000) {
  set.seed(i)
  rand_0.9 <- sample(1:147, 132, replace = FALSE)
  Y <- ctl.dat$Age[rand_0.9]
  X <- model.matrix(Age~., data = ctl.dat[rand_0.9,c(4,5,7,idx.sig)])[, -1]
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.9$Param)
  count.matrix_0.9$Non.Zero[idx_i] <- count.matrix_0.9$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.9) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (90% samples, 6 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")

# Lastly, using 70% samples:
count.matrix_0.7 <- data.frame(Param = colnames(X), Non.Zero = numeric(8))
for (i in 1:1000) {
  set.seed(i)
  rand_0.7 <- sample(1:147, 103, replace = FALSE)
  Y <- ctl.dat$Age[rand_0.7]
  X <- model.matrix(Age~., data = ctl.dat[rand_0.7,c(4,5,7,idx.sig)])[, -1]
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.9$Param)
  count.matrix_0.7$Non.Zero[idx_i] <- count.matrix_0.7$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.7) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (70% samples, 6 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")
# Save results
write.csv(count.matrix %>% arrange(1000-Non.Zero), 
          "Stability selection (100-pct data, 6proteins).csv",
          row.names = FALSE)
write.csv(count.matrix_0.9 %>% arrange(1000-Non.Zero), 
          "Stability selection (90-pct data, 6proteins).csv",
          row.names = FALSE)
write.csv(count.matrix_0.7 %>% arrange(1000-Non.Zero), 
          "Stability selection (70-pct data, 6proteins).csv",
          row.names = FALSE)

# Understand results
sum(count.matrix_0.7$Non.Zero > 500) # 8
sum(count.matrix_0.7$Param[count.matrix_0.7$Non.Zero > 500] %in%
      count.matrix$Param[count.matrix$Non.Zero == 1000])   #7 out 8 parameters were in 1000 iterations on
count.matrix$Param[count.matrix$Non.Zero == 1000][!(count.matrix$Param[
  count.matrix$Non.Zero == 1000] %in% 
                       count.matrix_0.7$Param[count.matrix_0.7$Non.Zero > 500])]
# None Left

####################################################################

# 3) MODEL SELECTION WITH EXHAUSTIVE SEARCH:

# Now, let's perform exhaustive search for the top protein/parameter models up to 7
ss_7 <- count.matrix$Param[count.matrix$Non.Zero == 1000]
ss_7[1] <- "Sex"
idx3 <- match(ss_7, names(ctl.dat))
ctl.dat3 <- ctl.dat[,c(4,idx3)]

# Exhaustive search
fit.exh <- regsubsets(Age~., ctl.dat3, nvmax = 7,
                      method = 'exhaustive', really.big = TRUE)
sum.exh <- summary(fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(sum.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
plot(sum.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
plot(sum.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)

# Results of exh selection
all.models <- data.frame(Var.Num = 1:7,
                         Formula = character(7),
                         R.sqrd = numeric(7),
                         Cp = numeric(7),
                         BIC = numeric(7),
                         AIC = numeric(7),
                         Full.Cor = numeric(7))
for (i in 1:7){
  # Save the essentials
  all.models$R.sqrd[i] <- sum.exh$rsq[i]
  all.models$Cp[i] <- sum.exh$cp[i]
  all.models$BIC[i] <- sum.exh$bic[i]
  # Get the names
  names_i <- names(ctl.dat3[sum.exh$which[i,]])[-1]
  formula_i <- paste("Age~", paste(names_i, collapse = "+"), sep = "")
  # Save the formula
  all.models$Formula[i] <- formula_i
  # Fit model
  model_i <- lm(formula_i, ctl.dat)
  # AIC
  all.models$AIC[i] <- AIC(model_i)
  # Get correlation between actual and predicted
  all.models$Full.Cor[i] <- cor(ctl.dat$Age, model_i$fitted.values,
                                method = 'pearson')
}
all.models.exh <- all.models
View(all.models.exh)

# Let's generate the following loop

# Create empty model list
# Create empty data frames from RMSE, R.sqr, MAE
# For each model in fit.exh
# Perform 10-fold CV repeated 100x and save in list
# Save resample data to data.tables

model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = 7)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = 7)
MAE.exh <- matrix(0, nrow = 1000, ncol = 7)
for (i in 1:7){
  model_i <- train(formula(all.models.exh$Formula[i]), 
                   data = ctl.dat, method = "lm",
                   trControl = trainControl(method = "repeatedcv",
                                            number = 10,
                                            repeats = 100))
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$resample$RMSE
  R.sqr.exh[,i] <- model.list[[i]]$resample$Rsquared
  MAE.exh[,i] <- model.list[[i]]$resample$MAE
}

which.min(apply(RMSE.exh, 2, mean))  # 7
which.max(apply(R.sqr.exh, 2, mean)) # 7
which.min(apply(MAE.exh, 2, mean))   # 7

boxplot(RMSE.exh, main = 'RMSE for 100x 10-fold CV', 
        ylab = 'RMSE', xlab = 'Parameter #')
boxplot(R.sqr.exh, main = 'Rsquared for 100x 10-fold CV', 
        ylab = 'Rsqrd', xlab = 'Parameter #')
boxplot(MAE.exh, main = 'MAE for 100x 10-fold CV',
        ylab = 'MAE', xlab = 'Parameter #')

# Save results
all.models.exh <- as.data.frame(all.models.exh)
write.csv(all.models.exh, "All models exhaustive search 10-01-21.csv",
          row.names = FALSE)

# My repeatCV
my.repeat.cv <- function(formula, data, kfold, repeats) {
  # Generate output data frame
  output.df <- data.frame(Iter = character(kfold*repeats), 
                          Rsquared.train = numeric(kfold*repeats),
                          RMSE.train = numeric(kfold*repeats),
                          Train.cor = numeric(kfold*repeats),
                          Test.cor = numeric(kfold*repeats))
  # Y variable
  yvar <- strsplit(formula, split = '~')[[1]][1]
  # Generate kfold.vector
  kfold.vector <- numeric(nrow(data))
  # Calculate AVERAGE size per fold (will vary if kfold is not a factor of nrow(data))
  size.per.fold <- nrow(data)/kfold
  # For k folds 1:kfold
  for (k in 1:kfold) {
    # Determine size of fold
    kfold.size <- length((round(size.per.fold*(k-1))+1):round(size.per.fold*k))
    # Replace vector entries with values indicating fold
    kfold.vector[(round(size.per.fold*(k-1))+1):round(size.per.fold*k)] <- rep(k, kfold.size)
    # Ex: kfold = 10, nrow = 147, size.per.fold = 14.7
    # k = 1, size.per.fold*k = 14.7, 1:15 , n = 15
    # k = 2, size.per.fold*k = 29.4, 16:29, n = 14
    # k = 3, size.per.fold*k = 44.1, 30:44, n = 15
  }
  for (i in 1:repeats){
    # Generate random fold vector
    set.seed(i)
    random.folds.i <- sample(kfold.vector, nrow(data), replace = FALSE)
    for (j in 1:kfold) {
      # Define Iter
      output.df$Iter[((i-1)*kfold)+j] <- paste('Repeat', i, 'Fold', j)
      # Define train/test sets 1:kfold
      test.idx <- c(which(random.folds.i == j))
      test.set <- data[test.idx,]
      train.set <- data[-test.idx,]
      # Train model
      model.i <- lm(formula, train.set)
      # Calculate rsquared for this fold
      output.df$Rsquared.train[((i-1)*kfold)+j] <- summary(model.i)$r.squared
      # Calculate RMSE for this fold
      output.df$RMSE.train[((i-1)*kfold)+j] <- RMSE(model.i$fitted.values, train.set[,yvar])
      # Train correlation
      output.df$Train.cor[((i-1)*kfold)+j] <- cor(train.set[,yvar], model.i$fitted.values)
      # Test correlation
      agepred <- predict(model.i, test.set)
      output.df$Test.cor[((i-1)*kfold)+j] <- cor(test.set[,yvar], agepred)
    }
    # if (i %in% seq(from = 0, to = 1000, by = 100)) {
    #   print(paste("Just completed ", i, "th iteration"))
    # }
  }
  return(output.df)
}

# Apply to 7 models
model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = 7)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = 7)
Train.cor.exh <- matrix(0, nrow = 1000, ncol = 7)
Test.cor.exh <- matrix(0, nrow = 1000, ncol = 7)
for (i in 1:7){
  model_i <- my.repeat.cv(formula = all.models.exh$Formula[i],
                          data = ctl.dat,
                          kfold = 10, 
                          repeats = 100)
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$RMSE.train
  R.sqr.exh[,i] <- model.list[[i]]$Rsquared.train
  Train.cor.exh[,i] <- model.list[[i]]$Train.cor
  Test.cor.exh[,i] <- model.list[[i]]$Test.cor
  print(paste("Just completed ", i, "th iteration"))
}
which.min(apply(RMSE.exh, 2, mean))       # 7
which.max(apply(R.sqr.exh, 2, mean))      # 7
which.max(apply(Train.cor.exh, 2, mean))  # 7
which.max(apply(Test.cor.exh, 2, mean))   # 7 
apply(Test.cor.exh, 2, mean)              # Max Average cor = 0.6914832

# Save
write.csv(RMSE.exh, "RMSE for exh search top7 Olink on Soma models 10-01-21.csv",
          row.names = FALSE)
write.csv(R.sqr.exh, "R.sqr for exh search top7 Olink on Soma models 10-01-21.csv",
          row.names = FALSE)
write.csv(Train.cor.exh, "Train set correlations top7 Olink on Soma 10-01-21.csv",
          row.names = FALSE)
write.csv(Test.cor.exh, "Test set correlations top7 Olink on Soma 10-01-21.csv",
          row.names = FALSE)

# Plot
boxplot(Test.cor.exh, main = 'Test correlation for 100x 10-fold CV',
        ylab = 'Pearson Cor', xlab = 'Parameter #')

####################################################################

# 4) Apply to NC and ND's

# Predict
final_model <- lm(all.models.exh$Formula[7], ctl.dat)
ctl.dat$AgePred <- predict(final_model, ctl.dat)
ctl.dat$DeltaAge <- ctl.dat$AgePred - ctl.dat$Age
ND.dat$AgePred <- predict(final_model, ND.dat)
ND.dat$DeltaAge <- ND.dat$AgePred - ND.dat$Age
data1 <- rbind(ctl.dat, ND.dat)

# Diagnostic plots
plot(final_model,1)
plot(final_model,2)
plot(final_model,3)

top7_coef <- coef(final_model)
print(top7_coef)
write.csv(top7_coef, "Top7_coefficients_10-01-21.csv",
          row.names = TRUE)

# Generate residual plots
ggplot(data1, aes(x = Group, y = DeltaAge)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Site)) +
  labs(title = "Delta Age by Disease Category\n(7-parameter model)", 
       y = "Delta Age", x = "Disease Group") +
  theme(legend.position = "right")

# Stats
anova(lm(DeltaAge ~ Group, data1))
#Maria's results
# Response: DeltaAge
#            Df  Sum Sq Mean Sq F value    Pr(>F)    
# Group       3  3491.9 1163.98   23.82 1.733e-14 ***
# Residuals 538 26289.9   48.87 

#TopOlink On Soma 7 paramter model
#Response: DeltaAge
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Group       3    912 304.002  6.5158 0.0002464 ***
#  Residuals 538  25101  46.656                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Maria's results
# Is Delta Age different than 0? Answer: Use 1-sample t-test
t.test(data1$DeltaAge[data1$Group == 'Control'])$p.value # 1
t.test(data1$DeltaAge[data1$Group == 'PD'])$p.value      # 0.1369226
t.test(data1$DeltaAge[data1$Group == 'AD'])$p.value      # 0.07302485   (n.s.)
t.test(data1$DeltaAge[data1$Group == 'ALS'])$p.value     # 7.846516e-08 ****

# Top Olink 7 on Soma Is Delta Age different than 0? Answer: Use 1-sample t-test
t.test(data1$DeltaAge[data1$Group == 'Control'])$p.value # 1
t.test(data1$DeltaAge[data1$Group == 'PD'])$p.value      # 0.01143263
t.test(data1$DeltaAge[data1$Group == 'AD'])$p.value      # 0.001001572   
t.test(data1 $DeltaAge[data1$Group == 'ALS'])$p.value     # 0.01565371 ****

# Sensitivity test
AD.pvals <- numeric(35)
for (i in 1:35) {
  final_model <- lm(all.models.exh$Formula[i], ctl.dat)
  ctl.dat$AgePred <- predict(final_model, ctl.dat)
  ctl.dat$DeltaAge <- ctl.dat$AgePred - ctl.dat$Age
  ND.dat$AgePred <- predict(final_model, ND.dat)
  ND.dat$DeltaAge <- ND.dat$AgePred - ND.dat$Age
  data1 <- rbind(ctl.dat, ND.dat)
  AD.pvals[i] <- t.test(data1$DeltaAge[data1$Group == 'AD'])$p.value
}

# NOTE: Did some sensitivity tests* and the more parameters you include after 20, 
#       the smaller your p-value gets. Why is this? Is this because the 
# * Tried model size = c(1:35 and saved pvals)

# Scatterplots

# Figure out axis limits
min(data1$Age[data1$Group == 'Control'])     # 39.4
max(data1$Age[data1$Group == 'Control'])     # 91
min(data1$AgePred[data1$Group == 'Control']) # 54.74849
max(data1$AgePred[data1$Group == 'Control']) # 95.2525


# Train data
plot(x = data1$Age[data1$Group == 'Control'],
     y = data1$AgePred[data1$Group == 'Control'],
     xlim = c(40,92), ylim = c(40,92),
     type = 'p', col = 'blue', pch = 16,
     main = 'NC Full Set',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')
text(x = 80, y = 50, labels = "cor = 0.723331")

cor.test(data1$Age[data1$Group == 'Control'],data1$AgePred[data1$Group == 'Control'],
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

# AD: Axis limits
min(data1$Age[data1$Group == 'AD'])     # 56
max(data1$Age[data1$Group == 'AD'])     # 84
min(data1$AgePred[data1$Group == 'AD']) # 62.9038
max(data1$AgePred[data1$Group == 'AD']) # 92.00934

# AD: scatter
plot(x = data1$Age[data1$Group == 'AD'],
     y = data1$AgePred[data1$Group == 'AD'],
     xlim = c(56,90), ylim = c(56,90),
     type = 'p', col = 'blue', pch = 16,
     main = 'AD',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')
text(x = 80, y = 50, labels = "cor = 0.5284426")

cor.test(data1$Age[data1$Group == 'AD'],data1$AgePred[data1$Group == 'AD'],
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

# PD: Axis limits
min(data1$Age[data1$Group == 'PD'])     # 46.3
max(data1$Age[data1$Group == 'PD'])     # 86.7
min(data1$AgePred[data1$Group == 'PD']) # 39.6
max(data1$AgePred[data1$Group == 'PD']) # 103.1

# PD: scatter
plot(x = data1$Age[data1$Group == 'PD'],
     y = data1$AgePred[data1$Group == 'PD'],
     xlim = c(46,87), ylim = c(40,102),
     type = 'p', col = 'blue', pch = 16,
     main = 'PD',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')
text(x = 80, y = 50, labels = "cor = 0.7015891")

cor.test(data1$Age[data1$Group == 'PD'],data1$AgePred[data1$Group == 'PD'],
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)

# ALS: Axis limits
min(data1$Age[data1$Group == 'ALS'])     # 39
max(data1$Age[data1$Group == 'ALS'])     # 85
min(data1$AgePred[data1$Group == 'ALS']) # 52
max(data1$AgePred[data1$Group == 'ALS']) # 89

# ALS: scatter
plot(x = data1$Age[data1$Group == 'ALS'],
     y = data1$AgePred[data1$Group == 'ALS'],
     xlim = c(39,86), ylim = c(52,89),
     type = 'p', col = 'blue', pch = 16,
     main = 'ALS',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')
text(x = 80, y = 55, labels = "cor = 0.5642727")

cor.test(data1$Age[data1$Group == 'ALS'],data1$AgePred[data1$Group == 'ALS'],
         method = c("pearson"),
         exact = NULL, conf.level = 0.95, continuity = FALSE)
          
####################################################################
# NOW, REPEAT ALL OF THE ABOVE WITH SMALLER SUBSET
####################################################################

# Select targets
Targets <- corr.df[corr.df$pval < 0.05,] %>% select(Soma.Target) %>% distinct()
Targets <- as.character(Targets$Soma.Target)

# Select data
ctl.dat <- cbind(soma.dat[soma.dat$Group == 'Control',c(1:6,947:949)], 
                 soma.dat[soma.dat$Group == 'Control',Targets])
ND.dat <- cbind(soma.dat[soma.dat$Group != 'Control',c(1:6,947:949)], 
                soma.dat[soma.dat$Group != 'Control',Targets])

####################################################################

# 1) IDENTIFY AGE-ASSOCIATED PROTEINS:

# Generate data frame
model213.df <- data.frame(Protein = character(213), 
                          Age.coef = numeric(213), Age.pval = numeric(213),
                          Age.FDR = numeric(213),
                          Male.coef = numeric(213), Male.pval = numeric(213),
                          Male.FDR = numeric(213),
                          PC1.coef = numeric(213), PC1.pval = numeric(213),
                          PC1.FDR = numeric(213))

# Run loop and apply linear model to all proteins
for (i in 1:213) {
  model213.df$Protein[i] <- names(ctl.dat)[i+9]
  formula_i <- paste(model213.df$Protein[i], "~Age+Sex+PC1", sep = "")
  model_i <- lm(formula_i, ctl.dat)
  model213.df$Age.coef[i] <- summary(model_i)$coefficients[2,1]
  model213.df$Age.pval[i] <- summary(model_i)$coefficients[2,4]
  model213.df$Male.coef[i] <- summary(model_i)$coefficients[3,1]
  model213.df$Male.pval[i] <- summary(model_i)$coefficients[3,4]
  model213.df$PC1.coef[i] <- summary(model_i)$coefficients[4,1]
  model213.df$PC1.pval[i] <- summary(model_i)$coefficients[4,4]
}

# Adjust p-values
model213.df$Age.FDR <- p.adjust(model213.df$Age.pval, method = 'BH')
model213.df$Male.FDR <- p.adjust(model213.df$Male.pval, method = 'BH')
model213.df$PC1.FDR <- p.adjust(model213.df$PC1.pval, method = 'BH')
View(model213.df %>% arrange(Age.pval))

# Summary
sum(model213.df$Age.FDR < 0.05) # 58 age-associated proteins
model213.df <- model213.df %>% arrange(Age.FDR)
View(model213.df)
setwd("C:/Users/Maria/Desktop/CP Lab/2021/Aging Project/SomaLogic/Analysis with p = 213")
write.csv(model213.df, 
          file = "Linear model results (p = 213) 2-9-2021.csv",
          row.names = FALSE)
sig.proteins <- model213.df$Protein[model213.df$Age.FDR < 0.05]

#############################################################################

# 2) MODEL SELECTION WITH ONLY 58 TOP PROTEINS:

# Will run elastic net with 10-fold CV repeated times with different seeds

# First, using 100% samples:
idx.sig <- match(sig.proteins, names(ctl.dat))
Y <- ctl.dat$Age
X <- model.matrix(Age~., data = ctl.dat[,c(4,5,7,idx.sig)])[, -1]  # proteins,sex,PC1
count.matrix <- data.frame(Param = colnames(X), Non.Zero = numeric(60))
for (i in 1:1000) {
  set.seed(i)
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix$Param)
  count.matrix$Non.Zero[idx_i] <- count.matrix$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (100% samples, 58 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")
# Summary
sum(count.matrix$Non.Zero == 1000) # 25
sum(count.matrix$Non.Zero >= 900)  # 36

# Second, using 90% samples:
count.matrix_0.9 <- data.frame(Param = colnames(X), Non.Zero = numeric(60))
for (i in 1:1000) {
  set.seed(i)
  rand_0.9 <- sample(1:147, 132, replace = FALSE)
  Y <- ctl.dat$Age[rand_0.9]
  X <- model.matrix(Age~., data = ctl.dat[rand_0.9,c(4,5,7,idx.sig)])[, -1]
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.9$Param)
  count.matrix_0.9$Non.Zero[idx_i] <- count.matrix_0.9$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.9) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (90% samples, 58 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")

# Lastly, using 70% samples:
count.matrix_0.7 <- data.frame(Param = colnames(X), Non.Zero = numeric(60))
for (i in 1:1000) {
  set.seed(i)
  rand_0.7 <- sample(1:147, 103, replace = FALSE)
  Y <- ctl.dat$Age[rand_0.7]
  X <- model.matrix(Age~., data = ctl.dat[rand_0.7,c(4,5,7,idx.sig)])[, -1]
  fit_i <- cv.glmnet(X, Y, alpha = 0.5, nfolds = 10)
  coef_i <- coef(fit_i, s = "lambda.min")
  coef_i <- coef_i[which(coef_i !=0),]   # get the non=zero coefficients
  names_i <- names(coef_i)[-1]
  idx_i <- match(names_i, count.matrix_0.9$Param)
  count.matrix_0.7$Non.Zero[idx_i] <- count.matrix_0.7$Non.Zero[idx_i] + 1
  if (i %in% seq(from = 0, to = 1000, by = 100)) {
    print(paste("Just completed ", i, "th iteration"))
  }
}
ggplot(count.matrix_0.7) +
  geom_histogram(aes(x = Non.Zero), bins = 20) +
  labs(title = "Repeated Elastic Net with 10-fold CV
       (70% samples, 58 proteins+Sex,PC1)",
       x = "Times with non-zero coefficient", y = "Parameters")
# Save results
write.csv(count.matrix %>% arrange(1000-Non.Zero), 
          "Stability selection (100-pct data, 58proteins).csv",
          row.names = FALSE)
write.csv(count.matrix_0.9 %>% arrange(1000-Non.Zero), 
          "Stability selection (90-pct data, 58proteins).csv",
          row.names = FALSE)
write.csv(count.matrix_0.7 %>% arrange(1000-Non.Zero), 
          "Stability selection (70-pct data, 58proteins).csv",
          row.names = FALSE)

# Understand results
sum(count.matrix_0.7$Non.Zero > 500) # 27
sum(count.matrix_0.7$Param[count.matrix_0.7$Non.Zero > 500] %in%
      count.matrix$Param[count.matrix$Non.Zero == 1000])   # 25 of these were included using 100% data 
# Which 2 are not?
count.matrix_0.7$Param[count.matrix_0.7$Non.Zero > 500][!(count.matrix_0.7$Param[
  count.matrix_0.7$Non.Zero > 500] %in% count.matrix$Param[count.matrix$Non.Zero == 1000])]
# suPAR, TNF.sR.I

############################################################################

# 3) MODEL SELECTION WITH EXHAUSTIVE SEARCH:

# Now, let's perform exhaustive search for the top 25 proteins
ss_25 <- count.matrix$Param[count.matrix$Non.Zero == 1000]
idx3 <- match(ss_25, names(ctl.dat))
ctl.dat3 <- ctl.dat[,c(4,idx3)]

# Exhaustive search
fit.exh <- regsubsets(Age~., ctl.dat3, nvmax = 35,
                      method = 'exhaustive', really.big = TRUE)
sum.exh <- summary(fit.exh)

# Graph Rsq, Cp, BIC using training data
par(mfrow=c(3,1), mar=c(2.5,4,0.5,1), mgp=c(1.5,0.5,0))
plot(sum.exh$rsq, xlab="Number of predictors",
     ylab="R-squared", col="green", type="p", pch=16)
plot(sum.exh$cp, xlab="Number of predictors",
     ylab="Cp", col="red", type="p", pch=16)
plot(sum.exh$bic, xlab="Number of predictors",
     ylab="BIC", col="blue", type="p", pch=16)

# Results of exh selection
all.models <- data.frame(Var.Num = 1:25,
                         Formula = character(25),
                         R.sqrd = numeric(25),
                         Cp = numeric(25),
                         BIC = numeric(25),
                         AIC = numeric(25),
                         Full.Cor = numeric(25))
for (i in 1:25){
  # Save the essentials
  all.models$R.sqrd[i] <- sum.exh$rsq[i]
  all.models$Cp[i] <- sum.exh$cp[i]
  all.models$BIC[i] <- sum.exh$bic[i]
  # Get the names
  names_i <- names(ctl.dat3[sum.exh$which[i,]])[-1]
  formula_i <- paste("Age~", paste(names_i, collapse = "+"), sep = "")
  # Save the formula
  all.models$Formula[i] <- formula_i
  # Fit model
  model_i <- lm(formula_i, ctl.dat)
  # AIC
  all.models$AIC[i] <- AIC(model_i)
  # Get correlation between actual and predicted
  all.models$Full.Cor[i] <- cor(ctl.dat$Age, model_i$fitted.values,
                                method = 'pearson')
}
all.models.exh <- all.models
View(all.models.exh)

# Let's generate the following loop

# Create empty model list
# Create empty data frames from RMSE, R.sqr, MAE
# For each model in fit.exh
# Perform 10-fold CV repeated 100x and save in list
# Save resample data to data.tables

model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = 35)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = 35)
MAE.exh <- matrix(0, nrow = 1000, ncol = 35)
for (i in 1:35){
  model_i <- train(formula(all.models.exh$Formula[i]), 
                   data = ctl.dat, method = "lm",
                   trControl = trainControl(method = "repeatedcv",
                                            number = 10,
                                            repeats = 100))
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$resample$RMSE
  R.sqr.exh[,i] <- model.list[[i]]$resample$Rsquared
  MAE.exh[,i] <- model.list[[i]]$resample$MAE
}

which.min(apply(RMSE.exh, 2, mean))  # 23
which.max(apply(R.sqr.exh, 2, mean)) # 23
which.min(apply(MAE.exh, 2, mean))   # 28

boxplot(RMSE.exh, main = 'RMSE for 100x 10-fold CV', 
        ylab = 'RMSE', xlab = 'Parameter #')
boxplot(R.sqr.exh, main = 'Rsquared for 100x 10-fold CV', 
        ylab = 'Rsqrd', xlab = 'Parameter #')
boxplot(MAE.exh, main = 'MAE for 100x 10-fold CV',
        ylab = 'MAE', xlab = 'Parameter #')

# Save results
all.models.exh <- as.data.frame(all.models.exh)
write.csv(all.models.exh, "All models exhaustive search (p = 213) 2-9-21.csv",
          row.names = FALSE)

# My repeatCV
my.repeat.cv <- function(formula, data, kfold, repeats) {
  # Generate output data frame
  output.df <- data.frame(Iter = character(kfold*repeats), 
                          Rsquared.train = numeric(kfold*repeats),
                          RMSE.train = numeric(kfold*repeats),
                          Train.cor = numeric(kfold*repeats),
                          Test.cor = numeric(kfold*repeats))
  # Y variable
  yvar <- strsplit(formula, split = '~')[[1]][1]
  # Generate kfold.vector
  kfold.vector <- numeric(nrow(data))
  # Calculate AVERAGE size per fold (will vary if kfold is not a factor of nrow(data))
  size.per.fold <- nrow(data)/kfold
  # For k folds 1:kfold
  for (k in 1:kfold) {
    # Determine size of fold
    kfold.size <- length((round(size.per.fold*(k-1))+1):round(size.per.fold*k))
    # Replace vector entries with values indicating fold
    kfold.vector[(round(size.per.fold*(k-1))+1):round(size.per.fold*k)] <- rep(k, kfold.size)
    # Ex: kfold = 10, nrow = 147, size.per.fold = 14.7
    # k = 1, size.per.fold*k = 14.7, 1:15 , n = 15
    # k = 2, size.per.fold*k = 29.4, 16:29, n = 14
    # k = 3, size.per.fold*k = 44.1, 30:44, n = 15
  }
  for (i in 1:repeats){
    # Generate random fold vector
    set.seed(i)
    random.folds.i <- sample(kfold.vector, nrow(data), replace = FALSE)
    for (j in 1:kfold) {
      # Define Iter
      output.df$Iter[((i-1)*kfold)+j] <- paste('Repeat', i, 'Fold', j)
      # Define train/test sets 1:kfold
      test.idx <- c(which(random.folds.i == j))
      test.set <- data[test.idx,]
      train.set <- data[-test.idx,]
      # Train model
      model.i <- lm(formula, train.set)
      # Calculate rsquared for this fold
      output.df$Rsquared.train[((i-1)*kfold)+j] <- summary(model.i)$r.squared
      # Calculate RMSE for this fold
      output.df$RMSE.train[((i-1)*kfold)+j] <- RMSE(model.i$fitted.values, train.set[,yvar])
      # Train correlation
      output.df$Train.cor[((i-1)*kfold)+j] <- cor(train.set[,yvar], model.i$fitted.values)
      # Test correlation
      agepred <- predict(model.i, test.set)
      output.df$Test.cor[((i-1)*kfold)+j] <- cor(test.set[,yvar], agepred)
    }
    # if (i %in% seq(from = 0, to = 1000, by = 100)) {
    #   print(paste("Just completed ", i, "th iteration"))
    # }
  }
  return(output.df)
}

# Apply to 25 models
model.list <- list()
RMSE.exh <- matrix(0, nrow = 1000, ncol = 25)
R.sqr.exh <- matrix(0, nrow = 1000, ncol = 25)
Train.cor.exh <- matrix(0, nrow = 1000, ncol = 25)
Test.cor.exh <- matrix(0, nrow = 1000, ncol = 25)
for (i in 1:25){
  model_i <- my.repeat.cv(formula = all.models.exh$Formula[i],
                          data = ctl.dat,
                          kfold = 10, 
                          repeats = 100)
  model.list[[i]] <- model_i
  RMSE.exh[,i] <- model.list[[i]]$RMSE.train
  R.sqr.exh[,i] <- model.list[[i]]$Rsquared.train
  Train.cor.exh[,i] <- model.list[[i]]$Train.cor
  Test.cor.exh[,i] <- model.list[[i]]$Test.cor
  print(paste("Just completed ", i, "th iteration"))
}
which.min(apply(RMSE.exh, 2, mean))       # 25
which.max(apply(R.sqr.exh, 2, mean))      # 25
which.max(apply(Train.cor.exh, 2, mean))  # 25
which.max(apply(Test.cor.exh, 2, mean))   # 14 
apply(Test.cor.exh, 2, mean)              # Average cor = 0.8064776

# Save
write.csv(RMSE.exh, "RMSE for exh search models 2-9-21.csv",
          row.names = FALSE)
write.csv(R.sqr.exh, "R.sqr for exh search models 2-9-21.csv",
          row.names = FALSE)
write.csv(Train.cor.exh, "Train set correlations 2-9-21.csv",
          row.names = FALSE)
write.csv(Test.cor.exh, "Test set correlations 2-9-21.csv",
          row.names = FALSE)

# Plot
boxplot(Test.cor.exh, main = 'Test correlation for 100x 10-fold CV',
        ylab = 'Pearson Cor', xlab = 'Parameter #')

############################################################################

# 4) Apply to NC and ND's

# Predict
final_model <- lm(all.models.exh$Formula[15], ctl.dat)
ctl.dat$AgePred <- predict(final_model, ctl.dat)
ctl.dat$DeltaAge <- ctl.dat$AgePred - ctl.dat$Age
ND.dat$AgePred <- predict(final_model, ND.dat)
ND.dat$DeltaAge <- ND.dat$AgePred - ND.dat$Age
data1 <- rbind(ctl.dat, ND.dat)

# Diagnostic plots
plot(final_model,1)
plot(final_model,2)
plot(final_model,3)

# Generate residual plots
ggplot(data1, aes(x = Group, y = DeltaAge)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size = 1, aes(color = Site)) +
  labs(title = "Delta Age by Disease Category\n(21-parameter model)", 
       y = "Delta Age", x = "Disease Group") +
  theme(legend.position = "right")

# Stats
anova(lm(DeltaAge ~ Group, data1))
# Response: DeltaAge
#            Df  Sum Sq Mean Sq F value    Pr(>F)    
# Group       3  3491.9 1163.98   23.82 1.733e-14 ***
# Residuals 538 26289.9   48.87 

# Is Delta Age different than 0? Answer: Use 1-sample t-test
t.test(data1$DeltaAge[data1$Group == 'Control'])$p.value # 1
t.test(data1$DeltaAge[data1$Group == 'PD'])$p.value      # 0.1369226
t.test(data1$DeltaAge[data1$Group == 'AD'])$p.value      # 0.07302485   (n.s.)
t.test(data1$DeltaAge[data1$Group == 'ALS'])$p.value     # 7.846516e-08 ****

# Sensitivity test
AD.pvals <- numeric(35)
for (i in 1:35) {
  final_model <- lm(all.models.exh$Formula[i], ctl.dat)
  ctl.dat$AgePred <- predict(final_model, ctl.dat)
  ctl.dat$DeltaAge <- ctl.dat$AgePred - ctl.dat$Age
  ND.dat$AgePred <- predict(final_model, ND.dat)
  ND.dat$DeltaAge <- ND.dat$AgePred - ND.dat$Age
  data1 <- rbind(ctl.dat, ND.dat)
  AD.pvals[i] <- t.test(data1$DeltaAge[data1$Group == 'AD'])$p.value
}

# NOTE: Did some sensitivity tests* and the more parameters you include after 20, 
#       the smaller your p-value gets. Why is this? Is this because the 
# * Tried model size = c(1:35 and saved pvals)

# Scatterplots

# Figure out axis limits
min(data1$Age[data1$Group == 'Control'])     # 39.4
max(data1$Age[data1$Group == 'Control'])     # 91
min(data1$AgePred[data1$Group == 'Control']) # 43.7
max(data1$AgePred[data1$Group == 'Control']) # 92.6

# Train data
plot(x = data1$Age[data1$Group == 'Control'],
     y = data1$AgePred[data1$Group == 'Control'],
     xlim = c(40,92), ylim = c(40,92),
     type = 'p', col = 'blue', pch = 16,
     main = 'NC Full Set',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')
text(x = 80, y = 50, labels = "cor = 0.899")

# AD: Axis limits
min(data1$Age[data1$Group == 'AD'])     # 56.3
max(data1$Age[data1$Group == 'AD'])     # 84.3
min(data1$AgePred[data1$Group == 'AD']) # 59.21
max(data1$AgePred[data1$Group == 'AD']) # 90.15

# AD: scatter
plot(x = data1$Age[data1$Group == 'AD'],
     y = data1$AgePred[data1$Group == 'AD'],
     xlim = c(56,90), ylim = c(56,90),
     type = 'p', col = 'blue', pch = 16,
     main = 'AD',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')

# PD: Axis limits
min(data1$Age[data1$Group == 'PD'])     # 46.3
max(data1$Age[data1$Group == 'PD'])     # 86.7
min(data1$AgePred[data1$Group == 'PD']) # 39.6
max(data1$AgePred[data1$Group == 'PD']) # 103.1

# PD: scatter
plot(x = data1$Age[data1$Group == 'PD'],
     y = data1$AgePred[data1$Group == 'PD'],
     xlim = c(46,87), ylim = c(40,102),
     type = 'p', col = 'blue', pch = 16,
     main = 'PD',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')

# ALS: Axis limits
min(data1$Age[data1$Group == 'ALS'])     # 38.7
max(data1$Age[data1$Group == 'ALS'])     # 85.2
min(data1$AgePred[data1$Group == 'ALS']) # 53.4
max(data1$AgePred[data1$Group == 'ALS']) # 85.0

# ALS: scatter
plot(x = data1$Age[data1$Group == 'ALS'],
     y = data1$AgePred[data1$Group == 'ALS'],
     xlim = c(38,86), ylim = c(50,85),
     type = 'p', col = 'blue', pch = 16,
     main = 'ALS',
     xlab = 'Chronological Age',
     ylab = 'Proteomic Age')
abline(a = 1, b = 1, col = 'blue')




