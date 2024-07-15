library(vegan)
# An R script to analyze a panel of metabolites using linear models
library(tidyverse)
library(sjstats)
library(performance)
library(forestmodel)
library(emmeans)
library(pwr)
library(ggplot2)
#### Chapter 2 Walleye RDA ####
# Setting the document directory
# loading set of data
setwd("C:/Users/user/Desktop/Chapter 2 Analysis") 
Data_2018 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/2018_MET_FM.csv", quote="")
   View(Data_2018) 
Factors18 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Factors18.csv", quote="")
View(Factors18)
dim(Factors18)
Factors18T <- decostand(Factors18[,2:7], "normalize")
# factors only FL, , W,SEX, Age
FactorsT4 <- Factors18T[,c(2,3,5,6)]
Metabolites18 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Metabolites18.csv", quote="")
View(Metabolites18) 
Metabolites18T <- decostand(Metabolites18[,2:164], "normalize")
# factors FL SEX AGE SITE
Factors <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Factors_site.csv", quote="")
# rename
Factors4 <- Factors
Factors4 <- Factors[,c(1,2,4,5)]
head(Factors4)

# normalize
Factors4_transformed <- decostand(Factors4[,c(2,4)], "normalize")

dim(Factors18)
library(vegan)


#### Biogenic Amines Linear Models ####
#### Spermine, Spermidine & Putrescine ####
# linear model versus age
bio_lm <- glm(Spine.Age ~ Spermine + Spermidine + Putrescine, 
              data = Data_2018, family = "gaussian")
summary(bio_lm)
# Call:
#   glm(formula = Spine.Age ~ Spermine + Spermidine + Putrescine, 
#       family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.932816   1.726525   4.595 2.96e-05 ***
#   Spermine    -0.005353   0.046018  -0.116   0.9079    
# Spermidine  -0.191022   0.095837  -1.993   0.0517 .  
# Putrescine   0.140257   0.074318   1.887   0.0649 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 8.950572)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 447.53  on 50  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 277.44
# 
# Number of Fisher Scoring iterations: 2
check_model(bio_lm)              

forest_model(bio_lm)

#### Linear model only spermidine ####

# linear model versus age
spermidine_lm <- glm(Spine.Age ~ Spermidine, 
              data = Data_2018, family = "gaussian")
summary(spermidine_lm)
# Call:
#   glm(formula = Spine.Age ~ Spermidine, family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  8.76183    0.76258  11.490 6.93e-16 ***
#   Spermidine  -0.06856    0.04152  -1.651    0.105    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 9.232341)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 480.08  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 277.23
# 
# Number of Fisher Scoring iterations: 2

#### Linear model age vs spermidine & putrescine ####
# linear model versus age
spermidine_putrescine_lm <- glm(Spine.Age ~ Spermidine + Putrescine, 
              data = Data_2018, family = "gaussian")
summary(spermidine_putrescine_lm)

# Call:
#   glm(formula = Spine.Age ~ Spermidine + Putrescine, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.76258    0.90725   8.556 1.98e-11 ***
#   Spermidine  -0.19733    0.07827  -2.521   0.0149 *  
#   Putrescine   0.14098    0.07334   1.922   0.0602 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 8.777445)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 447.65  on 51  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 275.46
# 
# Number of Fisher Scoring iterations: 2

check_model(spermidine_putrescine_lm)              

forest_model(spermidine_putrescine_lm)


#### Amino acids linear models ####

#### Essential ####

#### branched chained iso, leu, val ####

# linear model versus age
iso_leu_val_lm <- glm(Spine.Age ~ Valine_Average 
                      + Leucine + Isoleucine_Average, 
                     data = Data_2018, family = "gaussian")
summary(iso_leu_val_lm)
# Call:
#   glm(formula = Spine.Age ~ Valine_Average + Leucine + Isoleucine_Average, 
#       family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         7.183180   0.763117   9.413 1.18e-12 ***
#   Valine_Average     -0.011262   0.011541  -0.976   0.3339    
# Leucine             0.015010   0.006059   2.477   0.0167 *  
#   Isoleucine_Average  0.005768   0.017935   0.322   0.7491    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 8.998823)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 449.94  on 50  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 277.73
# 
# Number of Fisher Scoring iterations: 2
check_model(iso_leu_val_lm)              

forest_model(iso_leu_val_lm)

#### linear model only leucine ####

leu_lm <- glm(Spine.Age ~ Leucine, 
              data = Data_2018, family = "gaussian")
summary(leu_lm)
# Call:
#   glm(formula = Spine.Age ~ Leucine, family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 6.669588   0.699064   9.541 5.14e-13 ***
#   Leucine     0.007999   0.004373   1.829   0.0731 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 9.12915)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 474.72  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 276.63
# 
# Number of Fisher Scoring iterations: 2

cor.test(Data_2018$Spine.Age, Data_2018$Leucine)
# Pearson's product-moment correlation
# 
# data:  Data_2018$Spine.Age and Data_2018$Leucine
# t = 1.8291, df = 52, p-value = 0.07312
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.02343531  0.48190313
# sample estimates:
#       cor 
# 0.2458681 

#### rest of essentials ####
#### linear models his, lys, met, phe, thr, trp ####

# linear model versus age
his_lys_met_phe_thr_trp_lm <- glm(Spine.Age ~ Histidine_Average
                                  + Lysine + Methionine + Phenylalanine_Average
                                  + Threonine_Average + Tryptophan, 
                     data = Data_2018, family = "gaussian")
summary(his_lys_met_phe_thr_trp_lm)
# Call:
#   glm(formula = Spine.Age ~ Histidine_Average + Lysine + Methionine + 
#         Phenylalanine_Average + Threonine_Average + Tryptophan, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            5.595475   1.018053   5.496 1.54e-06 ***
#   Histidine_Average     -0.043939   0.028275  -1.554   0.1269    
# Lysine                 0.027667   0.017368   1.593   0.1179    
# Methionine             0.010111   0.044589   0.227   0.8216    
# Phenylalanine_Average  0.054817   0.030246   1.812   0.0763 .  
# Threonine_Average      0.011868   0.009642   1.231   0.2245    
# Tryptophan            -0.094207   0.068921  -1.367   0.1782    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 6.389124)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 300.29  on 47  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 261.9
# 
# Number of Fisher Scoring iterations: 2

check_model(his_lys_met_phe_thr_trp_lm)              

forest_model(his_lys_met_phe_thr_trp_lm)

#### linear model only phenylalanine ####

phe_lm <- glm(Spine.Age ~ Phenylalanine_Average, 
              data = Data_2018, family = "gaussian")
summary(phe_lm)
# Call:
#   glm(formula = Spine.Age ~ Phenylalanine_Average, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            4.59794    0.94895   4.845 1.18e-05 ***
#   Phenylalanine_Average  0.04667    0.01306   3.572 0.000774 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 7.802209)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 405.71  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 268.15
# 
# Number of Fisher Scoring iterations: 2

check_model(phe_lm)              

forest_model(phe_lm)

plot(Spine.Age ~ Phenylalanine_Average, data= Data_2018)
plot(Phenylalanine_Average ~Spine.Age, data= Data_2018)
cor.test(Data_2018$Spine.Age, Data_2018$Phenylalanine_Average)
# Pearson's product-moment correlation
# 
# data:  Data_2018$Spine.Age and Data_2018$Phenylalanine_Average
# t = 3.5719, df = 52, p-value = 0.0007739
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1998579 0.6360336

#### Conditional AAs ####
#### linear model arg glut gly pro tyr ####
# linear model versus age
arg_glut_gly_pro_tyr_lm <- glm(Spine.Age ~ Arginine_Average
                                  + Glutamine + Glycine_Average + Proline
                                  + Tyrosine_Average, 
                                  data = Data_2018, family = "gaussian")
summary(arg_glut_gly_pro_tyr_lm)
# Call:
#   glm(formula = Spine.Age ~ Arginine_Average + Glutamine + Glycine_Average + 
#         Proline + Tyrosine_Average, family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       4.479556   1.013367   4.420 5.61e-05 ***
#   Arginine_Average  0.016706   0.011397   1.466 0.149239    
# Glutamine         0.029786   0.029950   0.995 0.324946    
# Glycine_Average   0.006039   0.003098   1.949 0.057099 .  
# Proline          -0.028225   0.007683  -3.674 0.000601 ***
#   Tyrosine_Average  0.012209   0.006953   1.756 0.085493 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 5.158598)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 247.61  on 48  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 249.48
# 
# Number of Fisher Scoring iterations: 2

check_model(arg_glut_gly_pro_tyr_lm)  

#### Glycine linear model ####

# linear model versus age
gly_lm <- glm(Spine.Age ~ Glycine_Average, 
                               data = Data_2018, family = "gaussian")
summary(gly_lm)
# Call:
#   glm(formula = Spine.Age ~ Glycine_Average, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     3.553631   1.100835   3.228 0.002160 ** 
#   Glycine_Average 0.008942   0.002233   4.004 0.000199 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 7.426803)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 386.19  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 265.48
# 
# Number of Fisher Scoring iterations: 2
check_model(gly_lm)


#### Proline linear model ####

# linear model versus age
pro_lm <- glm(Spine.Age ~ Proline, 
              data = Data_2018, family = "gaussian")
summary(pro_lm)
#Call:
#   glm(formula = Spine.Age ~ Proline, family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  8.492115   0.680693  12.476   <2e-16 ***
#   Proline     -0.011276   0.007709  -1.463     0.15    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 9.332527)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 485.29  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 277.82
# 
# Number of Fisher Scoring iterations: 2
check_model(pro_lm)

# linear model versus age natural log
proln_lm <- glm(Spine.Age ~ log(Proline), 
              data = Data_2018, family = "gaussian") # not significant
summary(proln_lm)
check_model(proln_lm)

#### Tyrosine linear model ####

# linear model versus age
tyr_lm <- glm(Spine.Age ~ Tyrosine_Average, 
              data = Data_2018, family = "gaussian")
summary(tyr_lm)

check_model(tyr_lm)

# Call:
#   glm(formula = Spine.Age ~ Tyrosine_Average, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      6.490472   0.583977  11.114 2.39e-15 ***
#   Tyrosine_Average 0.020481   0.007261   2.821  0.00676 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 8.426961)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 438.20  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 272.3
# 
# Number of Fisher Scoring iterations: 2

# linear model versus  lon
tyrln_lm <- glm(Spine.Age ~ log(Tyrosine_Average), 
              data = Data_2018, family = "gaussian")
summary(tyrln_lm)
# Call:
#   glm(formula = Spine.Age ~ log(Tyrosine_Average), family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            -0.9254     2.5086  -0.369  0.71372   
# log(Tyrosine_Average)   2.2298     0.6407   3.480  0.00102 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 7.880742)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 409.80  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 268.69
# 
# Number of Fisher Scoring iterations: 2
check_model(tyrln_lm)
 # better


#### non essential amino acids ####
#### linear model ala aspn asp glut ser ####
# linear model versus  lon
ala_apn_asp_glu_lm <- glm(Spine.Age ~ Alanine + Asparagine + Aspartate + Glutamate_Average, 
                data = Data_2018, family = "gaussian")
summary(ala_apn_asp_glu_lm)
# Call:
#   glm(formula = Spine.Age ~ Alanine + Asparagine + Aspartate + 
#         Glutamate_Average, family = "gaussian", data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        4.1337595  0.9909565   4.171 0.000123 ***
#   Alanine           -0.0002525  0.0028899  -0.087 0.930741    
# Asparagine        -0.0748388  0.0201091  -3.722 0.000511 ***
#   Aspartate          0.0374811  0.0529044   0.708 0.482011    
# Glutamate_Average  0.0455461  0.0092741   4.911 1.05e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 6.147836)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 301.24  on 49  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 258.07
# 
# Number of Fisher Scoring iterations: 2

check_model(ala_apn_asp_glu_lm)
forest_model(ala_apn_asp_glu_lm)


# linear model versus  lon
glutm_lm <- glm(Spine.Age ~ Glutamate_Average, 
                          data = Data_2018, family = "gaussian")
summary(glutm_lm)
# Call:
#   glm(formula = Spine.Age ~ Glutamate_Average, family = "gaussian", 
#       data = Data_2018)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       4.760954   0.948968   5.017 6.49e-06 ***
#   Glutamate_Average 0.031816   0.009383   3.391  0.00134 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 7.957147)
# 
# Null deviance: 505.26  on 53  degrees of freedom
# Residual deviance: 413.77  on 52  degrees of freedom
# (6 observations deleted due to missingness)
# AIC: 269.21
# 
# Number of Fisher Scoring iterations: 2

#### Cancor factominer amino acids ####
# reading data of amino acids
Amino_Acids18 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Amino_Acids18.csv", quote="")
   View(Amino_Acids18)
str(Amino_Acids18)
dim(Amino_Acids18)
num_amino<- Amino_Acids18[,9:27]
factor_amino<- Amino_Acids18[,1:8]
library(FactoMineR)
library(factoextra)
AA.ca <- CA(Amino_Acids18, 
             row.sup = 15:18,  # Supplementary rows
             col.sup = 6:8,    # Supplementary columns
             graph = FALSE)
# another try
# Install and load the FactoMineR package
install.packages("FactoMineR")
library(FactoMineR)

# Generate sample data
set.seed(123)
data <- data.frame(
  X1 = rnorm(100),
  X2 = rnorm(100),
  X3 = rnorm(100),
  Y1 = rnorm(100),
  Y2 = rnorm(100)
)

# Perform canonical correlation analysis
cancor_result <- cancor(num_amino, factor_amino$Spine.Age)

# Summary of canonical correlations
summary(cancor_result)
# Length Class  Mode   
# cor       1    -none- numeric
# xcoef   361    -none- numeric
# ycoef     1    -none- numeric
# xcenter  19    -none- numeric
# ycenter   1    -none- numeric
cancor_result$cor
#[1] 0.8191136
# Canonical loadings for Y variables
cancor_result$cor.y

# Canonical coefficients for X variables
cancor_result$xcoef
#[,1]          [,2]          [,3]          [,4]
# Serine                -8.855798e-04  6.836656e-04  4.034046e-03 -1.236359e-03
# Asparagine             4.327983e-04  5.757311e-03  5.775242e-03 -5.526846e-03
# Glutamine              1.391110e-03 -2.593404e-04 -1.693521e-02 -8.064533e-04
# Tryptophan             7.247519e-03 -1.351133e-03 -7.116755e-05  2.167737e-02
# Lysine                 2.173915e-04 -4.052763e-05 -2.134692e-06  1.039134e-05
# Arginine_Average      -7.832878e-04  1.460260e-04  7.691553e-06 -3.744124e-05
# Glutamate_Average     -3.249601e-04  6.058132e-05  3.190970e-06 -1.553313e-05
# Glycine_Average       -1.864387e-04  3.475719e-05  1.830748e-06 -8.911788e-06
# Histidine_Average      3.819441e-04 -7.120467e-05 -3.750528e-06  1.825697e-05
# Isoleucine_Average     2.392521e-04 -4.460304e-05 -2.349354e-06  1.143628e-05
# Leucine                5.606616e-04 -1.045224e-04 -5.505458e-06  2.679968e-05
# Methionine             4.878769e-04 -9.095340e-05 -4.790743e-06  2.332057e-05
# Phenylalanine_Average -2.994478e-03  5.582515e-04  2.940450e-05 -1.431364e-04
# Proline                2.717484e-03 -5.066123e-04 -2.668454e-05  1.298960e-04
# Threonine_Average     -6.991158e-05  1.303340e-05  6.865020e-07 -3.341781e-06
# Tyrosine_Average      -3.656333e-04  6.816390e-05  3.590363e-06 -1.747731e-05
# Valine_Average        -1.202252e-03  2.241321e-04  1.180560e-05 -5.746776e-05
# Alanine               -3.854077e-05  7.185038e-06  3.784539e-07 -1.842253e-06
# Aspartate             -4.479217e-03  8.350469e-04  4.398401e-05 -2.141071e-04
# [,5]          [,6]          [,7]          [,8]
# Serine                 3.266578e-03 -2.050983e-03 -1.486415e-03 -1.710626e-03
# Asparagine            -2.316787e-03  3.141740e-03 -1.499510e-03  7.995848e-04
# Glutamine              8.849942e-04 -1.825454e-03 -3.972293e-03 -1.589338e-03
# Tryptophan            -3.106497e-03 -4.242382e-03 -1.854748e-03 -8.058796e-03
# Lysine                -7.786173e-03 -4.949694e-03 -6.166278e-04 -1.048330e-03
# Arginine_Average       3.975578e-06  7.829218e-03 -4.286372e-04 -4.130123e-04
# Glutamate_Average      1.649335e-06 -3.665537e-07  5.080807e-03 -1.138227e-03
# Glycine_Average        9.462696e-07 -2.103021e-07 -2.579529e-05  1.483229e-03
# Histidine_Average     -1.938557e-06  4.308313e-07  5.284504e-05  1.748945e-05
# Isoleucine_Average    -1.214324e-06  2.698754e-07  3.310246e-05  1.095550e-05
# Leucine               -2.845638e-06  6.324240e-07  7.757206e-05  2.567303e-05
# Methionine            -2.476219e-06  5.503231e-07  6.750170e-05  2.234017e-05
# Phenylalanine_Average  1.519848e-05 -3.377759e-06 -4.143103e-04 -1.371190e-04
# Proline               -1.379259e-05  3.065311e-06  3.759858e-04  1.244352e-04
# Threonine_Average      3.548363e-07 -7.885998e-08 -9.672832e-06 -3.201293e-06
# Tyrosine_Average       1.855772e-06 -4.124328e-07 -5.058832e-05 -1.674257e-05
# Valine_Average         6.102029e-06 -1.356135e-06 -1.663412e-04 -5.505182e-05
# Alanine                1.956137e-07 -4.347383e-08 -5.332426e-06 -1.764805e-06
# Aspartate              2.273427e-05 -5.052539e-06 -6.197359e-04 -2.051060e-04
# [,9]         [,10]         [,11]         [,12]
# Serine                 1.673839e-03 -1.413092e-03 -2.992344e-03 -2.229212e-03
# Asparagine            -8.426117e-04 -2.746607e-03  5.408042e-04 -3.764451e-03
# Glutamine             -7.460644e-03  2.654346e-03 -2.160178e-03 -2.284634e-03
# Tryptophan            -7.301792e-03 -2.271235e-02  4.242711e-03 -4.856691e-03
# Lysine                 1.162858e-03 -2.416227e-03  3.555793e-03  3.403477e-03
# Arginine_Average      -3.981160e-03  2.514252e-03  4.132554e-03 -4.898538e-03
# Glutamate_Average      6.920706e-04 -4.912166e-04 -1.389039e-04 -3.879231e-04
# Glycine_Average        4.687146e-05  1.030655e-04  3.297643e-04  8.697380e-05
# Histidine_Average      1.100303e-02 -4.278899e-04 -6.501986e-04 -6.996397e-04
# Isoleucine_Average    -1.511248e-05  4.184381e-03  2.292784e-03  7.740413e-04
# Leucine               -3.541447e-05 -4.838530e-05 -2.867933e-03 -9.686200e-04
# Methionine            -3.081699e-05 -4.210395e-05 -1.390935e-05  2.092174e-02
# Phenylalanine_Average  1.891477e-04  2.584246e-04  8.537247e-05 -1.310754e-04
# Proline               -1.716512e-04 -2.345199e-04 -7.747537e-05  1.189507e-04
# Threonine_Average      4.416000e-06  6.033395e-06  1.993177e-06 -3.060195e-06
# Tyrosine_Average       2.309541e-05  3.155429e-05  1.042419e-05 -1.600463e-05
# Valine_Average         7.594084e-05  1.037548e-04  3.427615e-05 -5.262539e-05
# Alanine                2.434447e-06  3.326082e-06  1.098796e-06 -1.687020e-06
# Aspartate              2.829320e-04  3.865581e-04  1.277023e-04 -1.960659e-04
# [,13]         [,14]         [,15]         [,16]
# Serine                -1.225608e-03  1.080900e-03  1.900226e-04 -6.305231e-04
# Asparagine             7.879609e-05 -8.893880e-03  2.986244e-03  1.768159e-05
# Glutamine             -9.371157e-04  1.500495e-03  9.164046e-05 -4.413041e-03
# Tryptophan            -1.833819e-02  6.106442e-03 -1.971187e-03 -1.160469e-02
# Lysine                -2.279155e-03  1.834744e-03 -1.230858e-03  5.329640e-03
# Arginine_Average      -3.106956e-04 -1.470064e-03  2.942717e-03  9.166776e-04
# Glutamate_Average     -6.898692e-04  3.945312e-04 -1.490468e-03 -2.223314e-03
# Glycine_Average       -2.276400e-04 -1.255171e-04  6.201637e-04  1.408555e-04
# Histidine_Average      1.339902e-03 -1.998221e-04  5.230628e-03  4.714980e-03
# Isoleucine_Average    -3.964285e-05 -2.472764e-03  9.802233e-04  7.837330e-04
# Leucine                4.171810e-04 -1.388216e-04 -8.466674e-04 -3.682521e-03
# Methionine            -7.925249e-03  3.823364e-06  3.372413e-03 -1.310801e-03
# Phenylalanine_Average  1.312750e-02 -2.072291e-03 -2.017470e-04 -2.818342e-03
# Proline                1.727151e-04  5.872360e-03 -1.165860e-03  3.526472e-03
# Threonine_Average     -4.443370e-06  1.332036e-05 -4.453320e-03 -3.630550e-04
# Tyrosine_Average      -2.323855e-05  6.966467e-05  2.260402e-05  6.749811e-03
# Valine_Average        -7.641150e-05  2.290669e-04  7.432509e-05 -6.324218e-05
# Alanine               -2.449535e-06  7.343231e-06  2.382651e-06 -2.027364e-06
# Aspartate             -2.846856e-04  8.534321e-04  2.769122e-04 -2.356208e-04
# [,17]         [,18]         [,19]
# Serine                -5.469888e-04  1.968351e-03 -0.0018072646
# Asparagine            -3.479302e-03 -5.713533e-03 -0.0018816678
# Glutamine              1.685583e-03  4.599950e-03  0.0030360909
# Tryptophan            -2.694638e-02 -3.963453e-03 -0.0064016069
# Lysine                -1.453147e-03 -5.104154e-04  0.0008621998
# Arginine_Average       4.317362e-04 -2.012257e-03  0.0054594379
# Glutamate_Average     -1.227959e-03  2.731515e-04 -0.0023684743
# Glycine_Average        4.396125e-05 -1.621769e-04  0.0005770208
# Histidine_Average      5.571720e-04 -8.625690e-06  0.0043488062
# Isoleucine_Average    -5.135382e-03 -1.872980e-03  0.0029145364
# Leucine               -1.142421e-03 -2.644456e-04 -0.0014823883
# Methionine             5.057628e-03  8.694562e-04 -0.0096564585
# Phenylalanine_Average  3.394693e-03 -1.432014e-03  0.0001589290
# Proline               -2.513044e-03  4.412996e-03 -0.0008710108
# Threonine_Average      1.221224e-04  4.104107e-04 -0.0024340800
# Tyrosine_Average      -2.912286e-03  1.346933e-03  0.0022072747
# Valine_Average         7.205996e-03  1.415790e-05 -0.0006191544
# Alanine               -3.829468e-06 -1.395110e-03 -0.0004449174
# Aspartate             -4.450617e-04  1.940390e-04  0.0288149552
# Canonical coefficients for Y variables
cancor_result$ycoef
# [,1]
# [1,] -0.044488
plot(cancor_result, choix = "ind", habillage = 1)

#### correlation analysis ####
# Generate example data
set.seed(123)
continuous_var <- rnorm(100)  # One continuous variable
other_continuous_vars <- matrix(rnorm(100*20), ncol = 20)  # 20 other continuous variables

#### Compute correlations ####
correlations <- cor(Amino_Acids18$Spine.Age, num_amino)

# Print correlations
print(correlations)
#      Serine     Asparagine Glutamine Tryptophan    Lysine Arginine_Average
# [1,] 0.6193454 -0.2720322  0.325698 0.06018503 0.5335207        0.5717329
#       Glutamate_Average Glycine_Average Histidine_Average Isoleucine_Average
# [1,]          0.425524       0.4854403           0.16368         0.01377828
#       Leucine Methionine Phenylalanine_Average    Proline Threonine_Average
# [1,] 0.2458681  0.4148477             0.4438654 -0.1987963         0.3349239
#      Tyrosine_Average Valine_Average   Alanine   Aspartate
# [1,]        0.3643056      0.0133911 0.1418183 -0.08054285

# Load the car package
library(car)

#### Amino Acids Perform MANCOVA ####
manova_result <- manova(cbind(Amino_Acids18$Serine, Amino_Acids18$Asparagine, 
                              Amino_Acids18$Glutamine, Amino_Acids18$Tryptophan,
                              Amino_Acids18$Lysine, Amino_Acids18$Arginine_Average,
                              Amino_Acids18$Glutamate_Average, Amino_Acids18$Glycine_Average,
                              Amino_Acids18$Histidine_Average, Amino_Acids18$Isoleucine_Average,
                              Amino_Acids18$Leucine, Amino_Acids18$Methionine,
                              Amino_Acids18$Phenylalanine_Average, Amino_Acids18$Proline,
                              Amino_Acids18$Threonine_Average, Amino_Acids18$Tyrosine_Average,
                              Amino_Acids18$Alanine, Amino_Acids18$Aspartate)  ~ Amino_Acids18$Spine.Age + Amino_Acids18$Site)
# Summarize MANOVA results
summary(manova_result)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Amino_Acids18$Spine.Age  1 0.73598   4.9557     18     32 4.219e-05 ***
#   Amino_Acids18$Site       3 1.95045   3.5102     54    102 2.371e-08 ***
#   Residuals               49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract coefficients
coefficients <- coef(manova_result)
coefficients
#                                [,1]        [,2]       [,3]        [,4]
# (Intercept)                  29.841298  43.9359791  13.676394  25.3730246
# Amino_Acids18$Spine.Age       9.900848  -0.4852388   2.0 how c h57201   0.7253119
# Amino_Acids18$SiteMatheson  -17.001300 -31.5653672 -16.329708 -17.0001449
# Amino_Acids18$SiteRed River   9.963954 -39.2907025  -9.461043 -14.4296647
# Amino_Acids18$SiteSandy Bar   6.619717 -37.8449426 -15.382178 -17.7484373
#                                 [,5]      [,6]       [,7]       [,8]       [,9]
# (Intercept)                 -21.576714  2.893556  67.579769  319.90031  49.087087
# Amino_Acids18$Spine.Age       4.883039  6.832202   6.864115   29.30521   1.185197
# Amino_Acids18$SiteMatheson    8.540807  9.455624 -46.897477 -158.26133  -7.347302
# Amino_Acids18$SiteRed River  13.343075 12.666121 -30.570170  -74.67669  -8.972630
# Amino_Acids18$SiteSandy Bar  32.089120 18.554266 -38.679427 -101.17952 -23.227843
#                                [,10]     [,11]      [,12]      [,13]       [,14]
# (Intercept)                  134.185205  83.96989  18.103429  50.702541 129.1906485
# Amino_Acids18$Spine.Age        5.007468  12.00304   2.674064   5.195651   0.3650292
# Amino_Acids18$SiteMatheson  -104.339784 -72.10467 -16.254398 -42.409999 -80.9847114
# Amino_Acids18$SiteRed River -108.531396 -83.48307 -10.819171 -25.173584 -92.4239267
# Amino_Acids18$SiteSandy Bar  -92.300267 -23.67351 -12.758299 -32.412439 -83.4908431
#                                 [,15]      [,16]      [,17]       [,18]
# (Intercept)                  79.969730  33.379121 244.012647  19.7420823
# Amino_Acids18$Spine.Age       7.450964   7.827507   5.082214   0.2301975
# Amino_Acids18$SiteMatheson   -2.748978 -55.682178 -13.992728  -8.2863029
# Amino_Acids18$SiteRed River -34.027250 -37.381356   4.670250 -13.0809655
# Amino_Acids18$SiteSandy Bar -49.796988 -51.555836 -51.116293 -15.1656214
# Plot the relationship between X and each dependent variable while controlling for C
for (i in 1:18) {
  # Extract coefficients for the ith dependent variable
  coeffs <- coefficients[, i]
  
  # Plot the relationship
  ggplot(data = data.frame(X, Y = coeffs[1] + coeffs[2] * as.numeric(X) 
                           + coeffs[3] * C),
         aes(x = as.numeric(X), y = Y)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "X", y = paste0("Y", i)) +
    ggtitle(paste0("Effect of X on Y", i))
}

# trying to plot
# Load necessary packages
library(tidyr)
library(dplyr)

# Reshape your data to long format
Amino_Acids18_long <- Amino_Acids18 %>%
  pivot_longer(cols = c(Serine, Asparagine, Glutamine, Tryptophan, Lysine, Arginine_Average,
                        Glutamate_Average, Glycine_Average, Histidine_Average, Isoleucine_Average,
                        Leucine, Methionine, Phenylalanine_Average, Proline, Threonine_Average,
                        Tyrosine_Average, Alanine, Aspartate),
               names_to = "Amino_Acid",
               values_to = "Concentration")

# Plot using ggplot2
library(ggplot2)

ggplot(data = Amino_Acids18_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Amino_Acid, scales = "free_y") +
  labs(title = "Relationship between Spine Age, Site, and Amino Acid Concentration")
#### trying plot MANOVA ####

# Load necessary libraries
library(ggplot2)
library(GGally)

# Select only the dependent variables from your data
dep_vars <- Amino_Acids18[, c("Serine", "Asparagine", "Glutamine", "Tryptophan",
                              "Lysine", "Arginine_Average", "Glutamate_Average", 
                              "Glycine_Average", "Histidine_Average", "Isoleucine_Average",
                              "Leucine", "Methionine", "Phenylalanine_Average", 
                              "Proline", "Threonine_Average", "Tyrosine_Average",
                              "Alanine", "Aspartate")]

# Create a scatterplot matrix
ggpairs(dep_vars)

# Add a factor (e.g., 'Site') as a different color
ggpairs(dep_vars, mapping = ggplot2::aes(color = Amino_Acids18$Site))

# Create a scatterplot matrix with smaller text
ggpairs(dep_vars, mapping = ggplot2::aes(color = Amino_Acids18$Site)) +
  theme(text = element_text(size = 8))  # Adjust size as needed

# separate graphs versus age
# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Create a list to store plots
plots_list <- list()

# Loop through each amino acid
for (amino_acid in colnames(dep_vars)) {
  # Create a scatterplot for each amino acid vs. Spine.Age
  p <- ggplot(Amino_Acids18, aes_string(x = "Spine.Age", y = amino_acid)) +
    geom_point() +
    labs(x = "Spine Age", y = amino_acid) +
    theme(text = element_text(size = 8))  # Adjust text size as needed
  
  # Add the plot to the list
  plots_list[[amino_acid]] <- p
}

# Arrange the plots in a grid
grid.arrange(grobs = plots_list, ncol = 3)

# In this example, aes_string() is used 
#to specify the x and y variables 
#in the ggplot() call. 
#This allows you to loop through 
#each amino acid and create 
#a scatterplot for each one.
# 
# The grid.arrange() function 
#from the gridExtra library 
#is then used to arrange all the scatterplots 
#in a grid, with 3 plots per row. 
#You can adjust the ncol argument 
#as needed to change the number of plots 
#per row

#### PCA Amino Acids ####      
#### load packages ####

library(tidyverse)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(sjstats)
library(performance)
library(forestmodel)
library(emmeans)
library(pwr)
# as data.frame
AA_data.frame <- data.frame(Amino_Acids18)
names(AA_data.frame)
dim(AA_data.frame)
# 54 27

# dataframe with just met data columns
AA_data_pca <- AA_data.frame[, c(9:27)] %>%
  print()

dim(AA_data_pca) 
#54 20
#### Correlations plot AAs ####

AA_cor <- cor(AA_data_pca)

summary(AA_cor)

corrplot(AA_cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
#dev.off()
library(corrplot)

#### PCA AAs ####

# PCA with FactoMineR package

AA_pca_FM <- PCA(AA_data_pca, graph = FALSE)


#### Eigenvalues / Variances AAs   ####

# extract eigenvalues/varianes
get_eig(AA_pca_FM)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  9.78439685       48.9219843                    48.92198
# Dim.2  3.87964687       19.3982344                    68.32022
# Dim.3  1.55198625        7.7599313                    76.08015
# Dim.4  0.97695765        4.8847883                    80.96494
# Dim.5  0.70445330        3.5222665                    84.48720
# Dim.6  0.63483362        3.1741681                    87.66137
# Dim.7  0.48732103        2.4366052                    90.09798
# Dim.8  0.40485900        2.0242950                    92.12227
# Dim.9  0.29732012        1.4866006                    93.60887
# Dim.10 0.27108970        1.3554485                    94.96432
# Dim.11 0.23655148        1.1827574                    96.14708
# Dim.12 0.18331540        0.9165770                    97.06366
# Dim.13 0.14662609        0.7331304                    97.79679
# Dim.14 0.12998310        0.6499155                    98.44670
# Dim.15 0.08666878        0.4333439                    98.88005
# Dim.16 0.07311227        0.3655614                    99.24561
# Dim.17 0.06886092        0.3443046                    99.58991
# Dim.18 0.03768843        0.1884421                    99.77835
# Dim.19 0.02784165        0.1392083                    99.91756
# Dim.20 0.01648746        0.0824373                   100.00000

corrplot(AA_cor, type = "upper", order = "hclust", 
         +          tl.col = "black", tl.srt = 60, 
         +          tl.cex = 1, cl.cex = 1)
AA_pca_FM <- PCA(AA_data_pca, graph = FALSE)
# extract eigenvalues/varianes
get_eig(AA_pca_FM)

fviz_screeplot(AA_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
fviz_screeplot(AA_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = 16
)

# extract the results for variables
var <- get_pca_var(AA_pca_FM)
var
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
var$coord  # coordinates of variables
# Dim.1       Dim.2        Dim.3        Dim.4        Dim.5
# Serine                0.7230922  0.55670020  0.145408849  0.061686837 -0.183945501
# Asparagine            0.5311752 -0.73223312  0.147563364 -0.117507759 -0.117460358
# Glutamine             0.8433966  0.09679767  0.282074923 -0.066800589 -0.227214762
# Tryptophan            0.7969608 -0.43581061 -0.198918811 -0.113309981 -0.008525766
# Lysine                0.4295974  0.77855332 -0.113669025  0.079232490 -0.058799775
# Arginine_Average      0.5685023  0.76207187 -0.007751513 -0.036789496  0.132302017
# Glutamate_Average     0.7723778  0.08686899  0.087922019  0.207271652 -0.357570572
# Glycine_Average       0.8083054  0.29095979 -0.014021850  0.021789435 -0.172926795
# Histidine_Average     0.7366175  0.09415033  0.157220141 -0.466144847  0.191315193
# Isoleucine_Average    0.7725156 -0.44012375 -0.257686120 -0.020714158 -0.127952426
# Leucine               0.7413415  0.08442662 -0.488223559  0.172688560  0.313575909
# Methionine            0.8861882  0.20490839  0.086145124 -0.019674310  0.176571046
# Phenylalanine_Average 0.9167304  0.11707764 -0.017697043  0.005519357 -0.088586941
# Proline               0.7104739 -0.58360391  0.098844334  0.029149447 -0.048057782
# Threonine_Average     0.7321295  0.18329453  0.278803287 -0.413608354  0.210898849
# Tyrosine_Average      0.7040285 -0.02383442 -0.461267372  0.294025477  0.208793775
# Valine_Average        0.7828489 -0.50032832 -0.307660821  0.058826111 -0.025263560
# Alanine               0.3475975  0.02174957  0.632091737  0.546340051  0.176637605
# Aspartate             0.3910232 -0.55601386  0.486883738  0.144488056  0.258654448
var$contrib # contributions of variables to the PCs
# Dim.1       Dim.2        Dim.3        Dim.4       Dim.5
# Serine                5.420025  8.79293876  1.362619133  0.415813041  5.14610601
# Asparagine            2.924755 15.21213032  1.403297986  1.508850414  2.09837463
# Glutamine             7.373567  0.26584047  5.127690066  0.487611203  7.85187680
# Tryptophan            6.583969  5.38873586  2.550025470  1.402973375  0.01105521
# Lysine                1.913098 17.19760679  0.832677974  0.685993037  0.52583770
# Arginine_Average      3.350259 16.47719022  0.003872271  0.147897482  2.66215487
# Glutamate_Average     6.184056  0.21410205  0.498182095  4.694538628 19.44571552
# Glycine_Average       6.772747  2.40191277  0.012670783  0.051880605  4.54804733
# Histidine_Average     5.624681  0.25149824  1.592975779 23.744054077  5.56671911
# Isoleucine_Average    6.186263  5.49592625  4.279320956  0.046886485  2.48998906
# Leucine               5.697056  0.20223210 15.361381443  3.258668893 14.95497766
# Methionine            8.140776  1.19127162  0.478249150  0.042297247  4.74175782
# Phenylalanine_Average 8.711585  0.38890116  0.020183396  0.003328818  1.19354774
# Proline               5.232511  9.66334954  0.629645877  0.092848251  0.35125930
# Threonine_Average     5.556352  0.95321380  5.009433201 18.693546453  6.76470418
# Tyrosine_Average      5.138003  0.01611761 13.711920309  9.446763073  6.63033513
# Valine_Average        6.352867  7.10234006  6.100103443  0.378140651  0.09707100
# Alanine               1.252469  0.01342123 25.748565303 32.616636670  4.74533331
# Aspartate             1.584962  8.77127115 15.277185365  2.281271595 10.17513763
var$cos2  # Cos2: quality on the factore map
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# Serine                0.5228624 0.3099151107 2.114373e-02 0.0038052659 3.383595e-02
# Asparagine            0.2821471 0.5361653459 2.177495e-02 0.0138080735 1.379694e-02
# Glutamine             0.7113179 0.0093697886 7.956626e-02 0.0044623186 5.162655e-02
# Tryptophan            0.6351465 0.1899308884 3.956869e-02 0.0128391517 7.268868e-05
# Lysine                0.1845539 0.6061452669 1.292065e-02 0.0062777875 3.457414e-03
# Arginine_Average      0.3231949 0.5807535310 6.008595e-05 0.0013534670 1.750382e-02
# Glutamate_Average     0.5965674 0.0075462213 7.730281e-03 0.0429615379 1.278567e-01
# Glycine_Average       0.6533576 0.0846575967 1.966123e-04 0.0004747795 2.990368e-02
# Histidine_Average     0.5426053 0.0088642839 2.471817e-02 0.2172910184 3.660150e-02
# Isoleucine_Average    0.5967803 0.1937089110 6.640214e-02 0.0004290764 1.637182e-02
# Leucine               0.5495872 0.0071278539 2.383622e-01 0.0298213389 9.832985e-02
# Methionine            0.7853295 0.0419874499 7.420982e-03 0.0003870785 3.117733e-02
# Phenylalanine_Average 0.8403947 0.0137071747 3.131853e-04 0.0000304633 7.847646e-03
# Proline               0.5047731 0.3405935291 9.770202e-03 0.0008496902 2.309550e-03
# Threonine_Average     0.5360137 0.0335968860 7.773127e-02 0.1710718706 4.447832e-02
# Tyrosine_Average      0.4956561 0.0005680797 2.127676e-01 0.0864509810 4.359484e-02
# Valine_Average        0.6128523 0.2503284244 9.465518e-02 0.0034605113 6.382475e-04
# Alanine               0.1208240 0.0004730436 3.995400e-01 0.2984874518 3.120084e-02
# Aspartate             0.1528992 0.3091514161 2.370558e-01 0.0208767983 6.690212e-02

#### graph of variables default plot ####
fviz_pca_var(AA_pca_FM, col.var = "black")
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(AA_pca_FM, col.var = "black", axes = c(2,3))
fviz_pca_var(AA_pca_FM, col.var = "black", axes = c(3,4))
fviz_pca_var(AA_pca_FM, col.var = "black", axes = c(4,5))

fviz_pca_var(AA_pca_FM, col.var = "black", axes = c(1,3))

#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(AA_pca_FM, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(AA_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Color by cos2 values: quality on the factor map
fviz_pca_var(AA_pca_FM, col.var = "cos2", axes = 2:3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

variables_PCA<- fviz_pca_var(AA_pca_FM, col.var = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                             repel = TRUE # Avoid text overlapping
)

# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(AA_pca_FM, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#### Color by a custom continuous variable ####
# Age

fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

fviz_pca_var(AA_pca_FM, col.var = AA_data.frame$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

indyear1_2<- fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Spine.Age,
                          geom.ind = "point", # show points only (nbut not "text"),
                          pointsize = 2.5,
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          legend.title = "Age") 
ggsave("indyear1_2.pdf")

fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Year, axes = c(2,3),
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Year")  
indyear2_3<- fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Year, axes = c(2,3),
                          geom.ind = "point", # show points only (nbut not "text"),
                          pointsize = 2.5,
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          legend.title = "Year") 
ggsave("indyear2_3.pdf")

# Load necessary packages
library(FactoMineR)
library(ggplot2)

# Example data
# Assuming your data is stored in a data frame called "my_data"
# It contains the PCA scores in columns PC1, PC2, etc., the variable you want to plot against (e.g., "variable"), and the factor (e.g., "factor")



# Extract the principal component scores for the axis you're interested in (e.g., PC1)
pc_scores <- AA_pca_FM$ind$coord[, "Dim.1"]

# Combine the principal component scores with the variable and factor
pca_AA_data <- data.frame(PC1 = AA_pca_FM, 
                       variable = factor_amino$Spine.Age,
                       factor = factor_amino$Site)

# Plot using ggplot2
ggplot(pca_data, aes(x = PC1, y = variable, color = factor)) +
  geom_point() +
  labs(title = "PCA Axis 1 vs Variable, Colored by Factor")
length(AA_data.frame)
dim(AA_data.frame)
length(AA_pca_FM)
length(factor_amino$Spine.Age)
dim(factor_amino)
length(factor_amino$Site)

#### trying heatmap on aas ####
# Define a distance function based on euclidean norm
# calculated between PCA values of the i-th and j-th items
dst <- Vectorize(function(i,j,dtset) sqrt(sum((dtset[i,2:3]-dtset[j,2:3])^2)), vectorize.args=c("i","j"))

# Calculate the distance matrix
nr <- nrow(AA_pca_FM)
mtx <- outer(1:nr, 1:nr, "dst", dtset=AA_pca_FM)
colnames(mtx) <- rownames(mtx) <- AA_pca_FM[,1]

# Plot the heatmap using ggplot2
library(reshape2)
library(ggplot2)
mtx.long <- melt(mtx)
ggplot(mtx.long, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+xlab("")+ylab("")

# Default plot
library(stats)
summary(AA_data_pca)
dim(AA_data_pca)
str(AA_data_pca)
# Assuming your data frame is named 'df'
AA_data_pca[,1:19] <- lapply(AA_data_pca[,1:19], as.numeric)
sum(is.na(AA_data_pca))

heatmap(AA_data_pca, scale = "column")

#### heatmap on correlations try 2 ####
view(Amino_Acids18)
dim(Amino_Acids18)
corAA_df <- Amino_Acids18[,8:27]
cor_df <- round(cor(corAA_df), 2)

library(reshape2)
melted_cor <- melt(cor_df)

ggplot(data = melted_cor, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  #geom_text(aes(Var2, Var1, label = value), size = 5) +
  scale_fill_gradient2(low = "blue", high = "red", 
                       limit = c(-1, 1), 
                       name = "Correlation") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.background = element_blank())

ggplot(data = melted_cor, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", limit = c(-1, 1), 
                       name = "Correlation") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
    panel.background = element_blank()
  )

# Example: Rename the columns
library(dplyr)
new_amino_acid <- Amino_Acids18 %>%
  rename(Valine_Average = Valine, 
         Tyrosine_Average = Tyrosine,
         Threonine_Average = Threonine,
         Phenylalanine_Average = Phenylalanine,
         Isoleucine_Average = Isoleucine,
         Histidine_Average = Histidine,
         Glycine_Average = Glycine,
         Glutamate_Average = Glutamate,
         Arginine_Average = Arginine,
         Spine.Age = Age)

names(Amino_Acids18)

names(Amino_Acids18) <- c( "ID", "Site", "TL", "FL", "W", "K", "Sex", 
                           "Age", "Serine", "Asparagine", "Glutamine", 
                           "Tryptophan", "Lysine", "Arginine", "Glutamate", 
                           "Glycine", "Histidine",   "Isoleucine",   "Leucine", 
                           "Methionine", "Phenylalanine", "Proline", "Threonine",  
                           "Tyrosine", "Valine", "Alanine", "Aspartate")

#### Correspondence analysis trying the aas ####
# normalize
aa_transformed <- decostand(num_amino, "normalize")
aa.ca2 <- cca(aa_transformed)
aa.ca2
# Call: cca(X = aa_transformed)
# 
# Inertia Rank
# Total          0.1424     
# Unconstrained  0.1424   18
# Inertia is scaled Chi-square 
# 
# Eigenvalues for unconstrained axes:
#   CA1     CA2     CA3     CA4     CA5     CA6     CA7     CA8 
# 0.05080 0.03212 0.01600 0.01168 0.00662 0.00547 0.00469 0.00357 
# (Showing 8 of 18 unconstrained eigenvalues)

summary(aa.ca2)

# 
# Call:
#   cca(X = aa_transformed) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          0.1424          1
# Unconstrained  0.1424          1
# 
# Eigenvalues, and their contribution to the scaled Chi-square 
# 
# Importance of components:
#   CA1     CA2    CA3     CA4      CA5      CA6      CA7      CA8
# Eigenvalue            0.0508 0.03212 0.0160 0.01168 0.006616 0.005475 0.004691 0.003575
# Proportion Explained  0.3568 0.22555 0.1124 0.08203 0.046460 0.038448 0.032940 0.025103
# Cumulative Proportion 0.3568 0.58231 0.6947 0.77674 0.823196 0.861644 0.894584 0.919687
# CA9     CA10     CA11     CA12     CA13      CA14      CA15
# Eigenvalue            0.002667 0.002259 0.001622 0.001143 0.001011 0.0008896 0.0008061
# Proportion Explained  0.018729 0.015865 0.011393 0.008026 0.007098 0.0062474 0.0056608
# Cumulative Proportion 0.938416 0.954281 0.965673 0.973700 0.980798 0.9870451 0.9927059
# CA16      CA17      CA18
# Eigenvalue            0.0005285 0.0003265 0.0001836
# Proportion Explained  0.0037115 0.0022930 0.0012896
# Cumulative Proportion 0.9964174 0.9987104 1.0000000
# 
# Scaling 2 for species and site scores
# * Species are scaled proportional to eigenvalues
# * Sites are unscaled: weighted dispersion equal on all dimensions
# 
# 
# Species scores
# 
# CA1      CA2      CA3       CA4       CA5        CA6
# Serine                -0.15813 -0.13080 -0.02801  0.002757 -0.061198  7.820e-02
# Asparagine             1.03500  0.51531 -0.06238  0.414600 -0.212082  4.218e-02
# Glutamine              0.07311 -0.05426  0.04336  0.241040 -0.193027  2.436e-01
# Tryptophan             0.15834  0.07162 -0.09123 -0.040480  0.071440  5.175e-02
# Lysine                -0.21255 -0.65570  0.34301  0.071419 -0.339573  1.718e-02
# Arginine_Average      -0.16327 -0.29997  0.13041  0.094714 -0.014900  1.443e-02
# Glutamate_Average     -0.06178  0.01631 -0.10368 -0.058289 -0.077306  1.370e-01
# Glycine_Average       -0.09140 -0.04421 -0.15146 -0.047790 -0.012561 -6.020e-02
# Histidine_Average     -0.01804 -0.05208 -0.04403  0.138350  0.082313  8.128e-03
# Isoleucine_Average     0.33402  0.06031  0.01094 -0.054136 -0.029327 -2.047e-05
# Leucine                0.17973 -0.22483  0.24079 -0.085006  0.066348 -1.086e-01
# Methionine             0.04397 -0.11207  0.07525  0.089595  0.044120  9.406e-02
# Phenylalanine_Average  0.03430 -0.02456 -0.04178 -0.026220 -0.005377  8.198e-02
# Proline                0.40937  0.21575  0.05356  0.093198 -0.121565 -9.191e-02
# Threonine_Average      0.00768 -0.11554 -0.06249  0.280436  0.142935 -7.069e-03
# Tyrosine_Average       0.14338 -0.08286  0.07599 -0.228429  0.203104  1.966e-01
# Valine_Average         0.33154  0.08218  0.01096 -0.066720  0.000233  1.471e-02
# Alanine               -0.30609  0.28619  0.13651  0.001510  0.005152 -1.036e-02
# Aspartate              0.05900  0.25460  0.02510  0.209074  0.144307  7.433e-02
# 
# 
# Site scores (weighted averages of species scores)
# 
# CA1      CA2      CA3       CA4      CA5      CA6
# sit1  -0.49474 -0.99745 -0.12938  0.067368  0.71143  2.46094
# sit2  -1.19381  0.61241 -1.30353 -0.876900  0.08496  0.02415
# sit3  -1.12990  0.99021 -0.55846 -0.688799  0.39121  1.64263
# sit4  -0.87235 -0.19172  1.09752 -0.574233 -0.97416  0.90373
# sit5   0.14032 -1.42380  0.57917 -0.633424  0.16424  0.13994
# sit6  -0.80181 -0.54000  0.72876  0.309401 -0.15571  0.32609
# sit7  -0.52901 -0.77253 -0.52880 -0.085293  1.37814  0.33778
# sit8  -0.58332 -1.10079  0.08926  0.783440 -0.22201  1.33873
# sit9  -0.94763 -1.05149 -0.23834 -0.596581 -2.58411  1.41502
# sit10 -0.50178  0.39790  0.43686  0.168156  0.99831  0.53441
# sit11 -0.97529 -0.48146 -0.47235 -0.318289 -0.51383 -0.47184
# sit12 -1.63160  1.44630 -0.69020 -0.892002  0.18441 -0.54783
# sit13 -0.60929 -0.47190 -2.02562 -0.992104  0.44430 -1.20063
# sit14 -0.50542  0.03739  0.18235 -1.664519  0.43234  0.98937
# sit15 -1.12506  0.45758 -1.52567 -0.716844 -0.96717  0.77498
# sit16 -0.27398 -1.77925  0.64146 -0.773905 -2.36669  0.53033
# sit17 -0.60086  0.57792 -0.99619 -1.077692  0.60409  0.10124
# sit18 -2.30854  2.58450  1.32716 -0.539561  0.28946 -1.15084
# sit19 -0.34696  0.26603 -0.07531 -1.067625  0.17617 -0.80911
# sit20  0.07266 -0.44037 -3.07457 -0.666778 -0.14293 -0.09056
# sit21 -0.22504 -1.90148  0.00863 -0.503435 -2.56262 -0.96673
# sit22 -0.19751 -0.82598  0.95043 -0.355975 -1.03416 -0.81639
# sit23  1.12449 -1.04835  1.26462  0.058050  0.41186 -2.58219
# sit24 -0.71095 -0.31017  1.63401 -0.584319 -0.68485 -1.43034
# sit25 -1.27062  3.26167  2.58770  1.063926 -0.41080 -1.13480
# sit26 -0.36745 -0.64609  0.82520  1.474230  0.65214 -0.09029
# sit27  0.42922 -0.99844 -0.41774  1.613355  1.37172  0.68988
# sit28  0.17281 -0.78528  0.82779  1.140220 -0.80432  0.07893
# sit29  0.08876 -0.50977 -0.49963  0.339553 -0.16200 -1.39792
# sit30 -0.26803 -0.68678 -0.21135  1.478094  0.67243 -0.74592
# sit31  0.24377 -0.36626 -0.03696  1.202431  0.89901 -0.15902
# sit32  0.59548 -0.84623  0.64670  2.160433  0.72358 -0.12809
# sit33  0.34564 -1.16101  0.41088 -0.076596  1.00881 -0.33665
# sit34 -0.40358  0.24949  1.45462  0.014248  0.35284  0.41643
# sit35 -1.02083 -0.06386 -2.89872  0.273915  1.34563 -2.61325
# sit36 -0.89427 -0.15812  0.46603  0.637658  0.17110 -0.03891
# sit37 -1.07168 -0.39448 -0.92714  0.286255  0.23054  0.38283
# sit38 -0.57358  0.01451  0.21788 -0.156968  0.34197 -0.20637
# sit39 -1.04086  0.91284 -0.07083  0.005662  0.67781 -0.42624
# sit40  1.10831  0.26669 -0.11316  1.395040  0.40778  0.85849
# sit41  0.54205  0.17452 -0.85491  1.416742 -0.35596  0.24664
# sit42  2.32895  0.79885 -0.26787 -0.460618 -1.08841  0.34131
# sit43  1.05561  0.56690 -0.63159 -0.473762  0.33830 -0.44328
# sit44  0.79431  1.16875 -0.26928  0.075259 -0.82656 -0.28182
# sit45 -1.17485  1.73172  1.38066  0.101783  0.42103  1.10889
# sit46  1.29720  1.00034  0.36802  0.461472  0.35125 -0.87192
# sit47  0.40031  0.09785 -0.66421  0.246122  1.36286  1.13991
# sit48  1.44073 -0.33255  0.35560 -0.963413 -0.14761 -0.69625
# sit49  0.97644  1.25231 -1.54333  1.170056 -0.36134  0.02520
# sit50  1.89245  0.90984 -0.12505 -1.264293 -0.69799 -0.85333
# sit51  1.76191  1.34786 -0.69963 -0.708463 -0.90334  0.75547
# sit52  0.56020 -0.01453 -0.09843 -0.267693 -0.31589 -1.77905
# sit53  2.19000 -0.41610  1.71722 -3.339693  2.67044  1.30090
# sit54  1.42696  1.78011 -0.15122  1.606117 -1.98205  1.50670
plot(aa.ca2, display = "species")


plot(aa.ca2, type="n")
text(aa.ca2,display="species")
points(aa.ca2,col="blue")


#### Carnitine book ####
# read data
Carnitine_book <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Carnitine_book.csv", quote="")
   View(Carnitine_book)
dim(Carnitine_book)
Carnitine_only <- Carnitine_book[,10:49]
dim(Carnitine_only)
names(Carnitine_only)

#### Compute correlations ####
correlations2 <- cor(Carnitine_book$Spine.Age, Carnitine_only)

# Print correlations
print(correlations2)
# C2      C3_1         C3       C4_1         C4       C3OH
# [1,] 0.3402395 -0.126548 -0.2052054 -0.2383346 -0.1764019 -0.1193076
# C5_1         C5      C4OH       C6_1         C6      C5OH
# [1,] -0.1920974 -0.2457064 -0.163636 -0.4467556 -0.0382808 0.0850139
# C5_1DC       C5DC        C8    C5MDC        C9      C7DC     C10_2
# [1,] -0.2111345 0.05670399 0.1578173 0.159258 0.2244565 0.1220914 0.2661331
# C10_1       C10     C12_1        C12       C14_2     C14_1        C14
# [1,] 0.1897864 0.2484103 0.1068175 0.07959819 -0.03558301 0.0935537 0.05424027
# C12DC   C14_2OH    C14_1OH     C16_2     C16_1        C16    C16_2OH
# [1,] 0.2535054 0.1887444 0.03678913 0.2784666 0.1826252 0.03926562 0.02541022
# C16_1OH      C16OH     C18_2     C18_1        C18    C18_1OH
# [1,] -0.001632518 0.05298148 0.2056525 0.1416085 0.04720983 0.08718539
# C0_Average
# [1,]  0.2796395

#### Perform MANCOVA ####
manova_result2 <- manova(cbind(Carnitine_book$C2, Carnitine_book$C3_1,
                              Carnitine_book$C3OH, Carnitine_book$C4,
                              Carnitine_book$C4_1, Carnitine_book$C4OH,
                              Carnitine_book$C5, Carnitine_book$C5_1,
                              Carnitine_book$C5_1DC, Carnitine_book$C5DC,
                              Carnitine_book$C5MDC, Carnitine_book$C5OH,
                              Carnitine_book$C6, Carnitine_book$C6_1,
                              Carnitine_book$C7DC, Carnitine_book$C8,
                              Carnitine_book$C9, Carnitine_book$C10,
                              Carnitine_book$C10_1, Carnitine_book$C10_2,
                              Carnitine_book$C12, Carnitine_book$C12_1,
                              Carnitine_book$C12DC, Carnitine_book$C14,
                              Carnitine_book$C14_1, Carnitine_book$C14_1OH,
                              Carnitine_book$C14_2, Carnitine_book$C14_2OH,
                              Carnitine_book$C16, Carnitine_book$C16_1,
                              Carnitine_book$C16OH, Carnitine_book$C16_1OH,
                              Carnitine_book$C16_2, Carnitine_book$C16_2OH,
                              Carnitine_book$C18, Carnitine_book$C18_1,
                              Carnitine_book$C18_1OH, Carnitine_book$C18_2) 
                        ~ Carnitine_book$Spine.Age + Carnitine_book$Site)
# Summarize MANOVA results
summary(manova_result2)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Carnitine_book$Spine.Age  1 0.91685   3.4822     38     12 0.0119917 *  
#   Carnitine_book$Site       3 2.63393   2.6508    114     42 0.0002722 ***
#   Residuals                49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract coefficients
coefficients2 <- coef(manova_result2)
coefficients2
# [,1]          [,2]          [,3]
# (Intercept)                  10.0888393  0.0608727063  0.0273300629
# Carnitine_book$Spine.Age      0.7463222 -0.0005988484 -0.0001479047
# Carnitine_book$SiteMatheson   5.4325259  0.0006915675 -0.0029932128
# Carnitine_book$SiteRed River  0.2574222 -0.0105836852 -0.0019333717
# Carnitine_book$SiteSandy Bar -5.3931319 -0.0237026061 -0.0021884762
# [,4]          [,5]         [,6]
# (Intercept)                   0.395223715  0.0880296534  0.106061723
# Carnitine_book$Spine.Age     -0.002993218 -0.0003920505 -0.003739219
# Carnitine_book$SiteMatheson  -0.175397065 -0.0143921892  0.011673470
# Carnitine_book$SiteRed River -0.197155035 -0.0318654383 -0.027521164
# Carnitine_book$SiteSandy Bar -0.268846705 -0.0265382054 -0.049300260
# [,7]          [,8]         [,9]
# (Intercept)                   0.539734314  2.830463e-02  0.146530218
# Carnitine_book$Spine.Age     -0.005231734  3.400512e-05 -0.001195298
# Carnitine_book$SiteMatheson  -0.308111341 -3.332141e-03 -0.038920544
# Carnitine_book$SiteRed River -0.383418452 -1.004882e-02 -0.042108381
# Carnitine_book$SiteSandy Bar -0.407556089 -7.353557e-03 -0.043994901
# [,10]         [,11]        [,12]
# (Intercept)                   0.202168780  0.0285583312  0.168616733
# Carnitine_book$Spine.Age      0.004705982  0.0002528951  0.002235253
# Carnitine_book$SiteMatheson   0.010184161 -0.0014744391  0.007192612
# Carnitine_book$SiteRed River -0.079745809  0.0002440691 -0.015759475
# Carnitine_book$SiteSandy Bar -0.055586438  0.0004445905 -0.016562862
# [,13]        [,14]       [,15]        [,16]
# (Intercept)                   0.85891658  0.239688251 0.018677431  0.211810813
# Carnitine_book$Spine.Age     -0.00855104 -0.004976040 0.001324616  0.003771977
# Carnitine_book$SiteMatheson   0.03372493 -0.003350073 0.016255589 -0.019563698
# Carnitine_book$SiteRed River  0.05780999 -0.025650006 0.012774562  0.032756340
# Carnitine_book$SiteSandy Bar -0.11463188 -0.011292431 0.014347350 -0.023912881
# [,17]         [,18]        [,19]
# (Intercept)                   0.0217333547  0.1327246055  0.186974544
# Carnitine_book$Spine.Age      0.0005210493  0.0007855886  0.002194546
# Carnitine_book$SiteMatheson  -0.0007656089 -0.0022912359 -0.009424970
# Carnitine_book$SiteRed River  0.0014859757  0.0029061164 -0.023215880
# Carnitine_book$SiteSandy Bar -0.0029070164 -0.0004063073 -0.023449293
# [,20]         [,21]        [,22]
# (Intercept)                  0.056891697  3.460815e-02  0.057315488
# Carnitine_book$Spine.Age     0.002084469  7.133717e-06  0.001087028
# Carnitine_book$SiteMatheson  0.015217148 -1.457605e-03 -0.008420110
# Carnitine_book$SiteRed River 0.030783033  7.537172e-03 -0.011738490
# Carnitine_book$SiteSandy Bar 0.025422955 -2.477933e-03 -0.016595676
# [,23]         [,24]        [,25]
# (Intercept)                  0.0137913596  0.0603630700  0.097310695
# Carnitine_book$Spine.Age     0.0003918906  0.0004700416  0.001756206
# Carnitine_book$SiteMatheson  0.0008791011 -0.0146607304 -0.022166183
# Carnitine_book$SiteRed River 0.0064726168  0.0119958669  0.001880141
# Carnitine_book$SiteSandy Bar 0.0035293698 -0.0167189028 -0.031818735
# [,26]         [,27]        [,28]
# (Intercept)                   4.070768e-02  0.0192690339 0.0166418618
# Carnitine_book$Spine.Age     -2.226488e-05 -0.0002635317 0.0006976008
# Carnitine_book$SiteMatheson  -8.492527e-03 -0.0010625848 0.0004903455
# Carnitine_book$SiteRed River  9.057914e-03  0.0029899680 0.0102343442
# Carnitine_book$SiteSandy Bar -7.659245e-03 -0.0016566006 0.0034519109
# [,29]        [,30]         [,31]
# (Intercept)                   0.1269388889  0.100201295  0.0122957624
# Carnitine_book$Spine.Age     -0.0008166667  0.005117690  0.0000775112
# Carnitine_book$SiteMatheson   0.0080055556 -0.013747125 -0.0014121838
# Carnitine_book$SiteRed River  0.0508000000  0.088190058  0.0004252975
# Carnitine_book$SiteSandy Bar  0.0081833333  0.005058548 -0.0027680593
# [,32]         [,33]         [,34]
# (Intercept)                   0.0460929004  0.0110594391  0.0261771646
# Carnitine_book$Spine.Age     -0.0006873001  0.0004864044 -0.0003384997
# Carnitine_book$SiteMatheson  -0.0026794732 -0.0009324867 -0.0005409672
# Carnitine_book$SiteRed River  0.0129193602  0.0050168394  0.0103498656
# Carnitine_book$SiteSandy Bar  0.0012557667 -0.0002576903  0.0010906110
# [,35]       [,36]         [,37]
# (Intercept)                   0.0406804702 0.135340671  1.893164e-02
# Carnitine_book$Spine.Age     -0.0001758637 0.003601999  1.097409e-04
# Carnitine_book$SiteMatheson   0.0011563244 0.003776934  3.839731e-05
# Carnitine_book$SiteRed River  0.0184427639 0.107806935  5.208829e-03
# Carnitine_book$SiteSandy Bar  0.0123919546 0.015101556 -4.632470e-04
# [,38]
# (Intercept)                   0.018590335
# Carnitine_book$Spine.Age      0.001120473
# Carnitine_book$SiteMatheson  -0.002824726
# Carnitine_book$SiteRed River  0.011821152
# Carnitine_book$SiteSandy Bar -0.002549047

# Reshape your data to long format
Carnitines_long <- Carnitine_book %>%
  pivot_longer(cols = c(C2, C3_1,C3OH, C4, C4_1, C4OH, C5, C5_1,
                        C5_1DC, C5DC,
                        C5MDC, C5OH,
                        C6, C6_1,
                        C7DC, C8,
                        C9, C10,
                        C10_1, C10_2,
                        C12, C12_1,
                        C12DC, C14,
                        C14_1, C14_1OH,
                        C14_2, C14_2OH,
                        C16, C16_1,
                        C16OH, C16_1OH,
                        C16_2, C16_2OH,
                        C18, C18_1,
                        C18_1OH, C18_2),
               names_to = "Carnitines",
               values_to = "Concentration")

# Plot using ggplot2
library(ggplot2)

ggplot(data = Carnitines_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Carnitines, scales = "free_y") +
  labs(title = "Relationship between Spine Age, Site, and Carnitine Concentration")

# as data.frame
carn_data.frame <- data.frame(Carnitine_only)
names(carn_data.frame)
dim(carn_data.frame)
# 54 27

# dataframe with just met data columns
carn_data_pca <- carn_data.frame %>%
  print()

#### Correlations plot ####

carn_cor <- cor(carn_data_pca)

summary(carn_cor)
# C2               C3_1                C3                C4_1         
# Min.   :0.07077   Min.   :-0.02872   Min.   :-0.11996   Min.   :-0.06577  
# 1st Qu.:0.30974   1st Qu.: 0.17057   1st Qu.: 0.07586   1st Qu.: 0.12227  
# Median :0.47357   Median : 0.23583   Median : 0.16928   Median : 0.18810  
# Mean   :0.41602   Mean   : 0.27177   Mean   : 0.22381   Mean   : 0.24810  
# 3rd Qu.:0.50554   3rd Qu.: 0.34246   3rd Qu.: 0.29330   3rd Qu.: 0.39628  
# Max.   :1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C4                C3OH              C5_1                C5          
# Min.   :-0.05129   Min.   :0.01279   Min.   :-0.12031   Min.   :-0.12575  
# 1st Qu.: 0.42655   1st Qu.:0.12016   1st Qu.: 0.04513   1st Qu.: 0.01529  
# Median : 0.54400   Median :0.16824   Median : 0.13816   Median : 0.11870  
# Mean   : 0.50516   Mean   :0.20467   Mean   : 0.20579   Mean   : 0.19751  
# 3rd Qu.: 0.59690   3rd Qu.:0.24150   3rd Qu.: 0.28357   3rd Qu.: 0.25488  
# Max.   : 1.00000   Max.   :1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C4OH              C6_1                 C6               C5OH          
# Min.   :0.05115   Min.   :-0.008159   Min.   :0.06763   Min.   :-0.008159  
# 1st Qu.:0.24440   1st Qu.: 0.191465   1st Qu.:0.42617   1st Qu.: 0.154102  
# Median :0.46544   Median : 0.327378   Median :0.68091   Median : 0.194533  
# Mean   :0.41612   Mean   : 0.310329   Mean   :0.58940   Mean   : 0.227887  
# 3rd Qu.:0.53232   3rd Qu.: 0.397840   3rd Qu.:0.76166   3rd Qu.: 0.244638  
# Max.   :1.00000   Max.   : 1.000000   Max.   :1.00000   Max.   : 1.000000  
# C5_1DC              C5DC                C8              C5MDC         
# Min.   :-0.01482   Min.   :-0.03117   Min.   :0.06763   Min.   :-0.06456  
# 1st Qu.: 0.39760   1st Qu.: 0.35097   1st Qu.:0.43750   1st Qu.: 0.22732  
# Median : 0.49772   Median : 0.40787   Median :0.78210   Median : 0.51076  
# Mean   : 0.47185   Mean   : 0.42403   Mean   :0.62790   Mean   : 0.40623  
# 3rd Qu.: 0.54461   3rd Qu.: 0.47986   3rd Qu.:0.84125   3rd Qu.: 0.55982  
# Max.   : 1.00000   Max.   : 1.00000   Max.   :1.00000   Max.   : 1.00000  
# C9              C7DC              C10_2             C10_1        
# Min.   :0.1303   Min.   :-0.12031   Min.   :-0.1258   Min.   :0.05851  
# 1st Qu.:0.4213   1st Qu.: 0.02493   1st Qu.: 0.1785   1st Qu.:0.22036  
# Median :0.7844   Median : 0.07774   Median : 0.5778   Median :0.24995  
# Mean   :0.6250   Mean   : 0.08944   Mean   : 0.4431   Mean   :0.28233  
# 3rd Qu.:0.8645   3rd Qu.: 0.13184   3rd Qu.: 0.6659   3rd Qu.:0.31510  
# Max.   :1.0000   Max.   : 1.00000   Max.   : 1.0000   Max.   :1.00000  
# C10             C12_1              C12              C14_2       
# Min.   :0.0154   Min.   :0.02987   Min.   :0.03221   Min.   :0.0760  
# 1st Qu.:0.4162   1st Qu.:0.47967   1st Qu.:0.39176   1st Qu.:0.3963  
# Median :0.7687   Median :0.74671   Median :0.82469   Median :0.7559  
# Mean   :0.5994   Mean   :0.62575   Mean   :0.63929   Mean   :0.6313  
# 3rd Qu.:0.8439   3rd Qu.:0.80904   3rd Qu.:0.92166   3rd Qu.:0.8957  
# Max.   :1.0000   Max.   :1.00000   Max.   :1.00000   Max.   :1.0000  
# C14_1              C14              C12DC             C14_2OH        
# Min.   :0.06973   Min.   :0.06925   Min.   :-0.09428   Min.   :-0.04236  
# 1st Qu.:0.43001   1st Qu.:0.40004   1st Qu.: 0.21826   1st Qu.: 0.32875  
# Median :0.83503   Median :0.81753   Median : 0.71101   Median : 0.80264  
# Mean   :0.67107   Mean   :0.65911   Mean   : 0.53285   Mean   : 0.63199  
# 3rd Qu.:0.93226   3rd Qu.:0.94996   3rd Qu.: 0.80224   3rd Qu.: 0.93938  
# Max.   :1.00000   Max.   :1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C14_1OH            C16_2             C16_1                C16         
# Min.   :0.08976   Min.   :-0.0557   Min.   :-0.007021   Min.   :0.02844  
# 1st Qu.:0.39864   1st Qu.: 0.2844   1st Qu.: 0.355819   1st Qu.:0.39675  
# Median :0.81100   Median : 0.7766   Median : 0.828562   Median :0.80244  
# Mean   :0.65700   Mean   : 0.6030   Mean   : 0.640844   Mean   :0.64116  
# 3rd Qu.:0.95033   3rd Qu.: 0.8912   3rd Qu.: 0.944138   3rd Qu.:0.93600  
# Max.   :1.00000   Max.   : 1.0000   Max.   : 1.000000   Max.   :1.00000  
# C16_2OH           C16_1OH            C16OH             C18_2        
# Min.   :0.03783   Min.   :0.06144   Min.   :0.02326   Min.   :0.01496  
# 1st Qu.:0.41092   1st Qu.:0.40387   1st Qu.:0.41105   1st Qu.:0.34269  
# Median :0.80773   Median :0.79002   Median :0.77672   Median :0.83064  
# Mean   :0.65632   Mean   :0.64961   Mean   :0.64307   Mean   :0.64859  
# 3rd Qu.:0.95248   3rd Qu.:0.93506   3rd Qu.:0.91748   3rd Qu.:0.95009  
# Max.   :1.00000   Max.   :1.00000   Max.   :1.00000   Max.   :1.00000  
# C18_1                C18              C18_1OH          C0_Average     
# Min.   :-0.009674   Min.   :-0.01104   Min.   :0.05228   Min.   :0.01117  
# 1st Qu.: 0.359742   1st Qu.: 0.38682   1st Qu.:0.39725   1st Qu.:0.29251  
# Median : 0.829898   Median : 0.77814   Median :0.80721   Median :0.48189  
# Mean   : 0.649463   Mean   : 0.62485   Mean   :0.65245   Mean   :0.42099  
# 3rd Qu.: 0.964931   3rd Qu.: 0.92585   3rd Qu.:0.94192   3rd Qu.:0.54321  
# Max.   : 1.000000   Max.   : 1.00000   Max.   :1.00000   Max.   :1.00000 
corrplot(carn_cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
dev.off()


#### PCA ####

# PCA with FactoMineR package

carn_pca_FM <- PCA(carn_data_pca, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(carn_pca_FM)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  2.263695e+01     56.592369213                    56.59237
# Dim.2  4.988967e+00     12.472417553                    69.06479
# Dim.3  1.895176e+00      4.737939943                    73.80273
# Dim.4  1.483528e+00      3.708820650                    77.51155
# Dim.5  1.287992e+00      3.219979673                    80.73153
# Dim.6  1.207845e+00      3.019613320                    83.75114
# Dim.7  8.469295e-01      2.117323852                    85.86846
# Dim.8  7.494791e-01      1.873697774                    87.74216
# Dim.9  6.919937e-01      1.729984183                    89.47215
# Dim.10 6.194931e-01      1.548732868                    91.02088
# Dim.11 5.284162e-01      1.321040520                    92.34192
# Dim.12 4.628633e-01      1.157158234                    93.49908
# Dim.13 4.015538e-01      1.003884431                    94.50296
# Dim.14 3.468046e-01      0.867011493                    95.36997
# Dim.15 3.316164e-01      0.829040909                    96.19901
# Dim.16 2.617825e-01      0.654456268                    96.85347
# Dim.17 2.373848e-01      0.593461953                    97.44693
# Dim.18 1.761680e-01      0.440420084                    97.88735
# Dim.19 1.539630e-01      0.384907382                    98.27226
# Dim.20 1.389589e-01      0.347397155                    98.61966
# Dim.21 9.638396e-02      0.240959906                    98.86062
# Dim.22 8.261812e-02      0.206545309                    99.06716
# Dim.23 6.802450e-02      0.170061238                    99.23722
# Dim.24 5.979856e-02      0.149496410                    99.38672
# Dim.25 5.559813e-02      0.138995335                    99.52572
# Dim.26 3.804372e-02      0.095109303                    99.62082
# Dim.27 3.178670e-02      0.079466746                    99.70029
# Dim.28 2.186383e-02      0.054659576                    99.75495
# Dim.29 1.792189e-02      0.044804717                    99.79976
# Dim.30 1.744082e-02      0.043602041                    99.84336
# Dim.31 1.572161e-02      0.039304036                    99.88266
# Dim.32 1.413258e-02      0.035331438                    99.91799
# Dim.33 8.468063e-03      0.021170157                    99.93916
# Dim.34 7.231397e-03      0.018078492                    99.95724
# Dim.35 4.890095e-03      0.012225238                    99.96947
# Dim.36 4.467680e-03      0.011169201                    99.98064
# Dim.37 2.755202e-03      0.006888004                    99.98752
# Dim.38 2.090729e-03      0.005226824                    99.99275
# Dim.39 1.933174e-03      0.004832936                    99.99758
# Dim.40 9.662544e-04      0.002415636                   100.00000
corrplot(carn_cor, type = "upper", order = "hclust", 
         +          tl.col = "black", tl.srt = 60, 
         +          tl.cex = 1, cl.cex = 1)
carn_pca_FM <- PCA(carn_data_pca, graph = FALSE)
# extract eigenvalues/varianes
get_eig(carn_pca_FM)

fviz_screeplot(carn_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
fviz_screeplot(carn_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = 16
)

# extract the results for variables
var3 <- get_pca_var(carn_pca_FM)
var3
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" 

var3$coord
# Dim.1       Dim.2       Dim.3        Dim.4        Dim.5
# C2         0.5705240  0.08878494  0.37582068 -0.053118510  0.033320711
# C3_1       0.3009570  0.57064471  0.08180654  0.032796640  0.011333119
# C3         0.2085785  0.73753659  0.07593569  0.232794714 -0.335729789
# C4_1       0.2379540  0.72997203  0.04912729  0.051624563  0.199114372
# C4         0.6373766  0.67197554 -0.19422267  0.050476018 -0.178784582
# C3OH       0.2179220  0.39171028  0.02956060  0.481756315  0.155236229
# C5_1       0.1775164  0.76344802  0.15266421 -0.100130877 -0.178147175
# C5         0.1625366  0.84783400 -0.19346253 -0.003250965 -0.230632304
# C4OH       0.5777071  0.21648618 -0.33068590  0.040779246  0.077351378
# C6_1       0.3934419  0.33736197 -0.33074540  0.023584426  0.555348344
# C6         0.8280668  0.14723359  0.10554197 -0.073230113  0.071957539
# C5OH       0.2543685  0.27775588  0.65795984  0.348783928 -0.012408795
# C5_1DC     0.5962648  0.61372650 -0.27836925 -0.075621009 -0.025804516
# C5DC       0.5291722  0.52478323  0.16481548 -0.415295790  0.117541594
# C8         0.8995066  0.05107445  0.12367972 -0.083121051 -0.072907328
# C5MDC      0.5990714 -0.08901571  0.04947379 -0.384805808 -0.035312872
# C9         0.9054296 -0.03169274  0.05601607 -0.021645116  0.001930235
# C7DC       0.1043711 -0.10517994  0.18944539  0.508499007  0.518179154
# C10_2      0.6732497 -0.28087183  0.33073596 -0.237729164 -0.005091980
# C10_1      0.3439148  0.32063008  0.31624985 -0.401809978  0.507216434
# C10        0.8870782 -0.13078922  0.09400225 -0.202488264  0.032614871
# C12_1      0.8816240  0.16800489  0.02127086 -0.246639375  0.150279670
# C12        0.9510426 -0.13135164 -0.06339824 -0.057992182 -0.068766026
# C14_2      0.9198261 -0.04024529 -0.10715097  0.174345223  0.024097144
# C14_1      0.9747282 -0.01164996 -0.02779031  0.034778870 -0.060256142
# C14        0.9726928 -0.07656936 -0.14065989  0.054703143 -0.058203822
# C12DC      0.8218828 -0.35523670  0.07553913 -0.010348409 -0.151062469
# C14_2OH    0.9518170 -0.25094079  0.02372481  0.042478929 -0.060289682
# C14_1OH    0.9704861 -0.08744574 -0.15424386  0.068004566 -0.014093342
# C16_2      0.9127627 -0.27736346  0.08784088 -0.001487240 -0.094833720
# C16_1      0.9606953 -0.20421987 -0.01361691  0.007395946 -0.041432932
# C16        0.9545146 -0.14582993 -0.14602655  0.072354595  0.044071892
# C16_2OH    0.9719150 -0.13192509 -0.08918920  0.073165236  0.043061782
# C16_1OH    0.9581388 -0.09600861 -0.14895232  0.140394266  0.042450711
# C16OH      0.9411718 -0.02094882 -0.13753202  0.039574013 -0.099310948
# C18_2      0.9659669 -0.17719586  0.02226547  0.066011625 -0.083412648
# C18_1      0.9702809 -0.19445164 -0.03290288  0.074561155 -0.005260420
# C18        0.9340043 -0.19364421 -0.11190754  0.131770336  0.082080666
# C18_1OH    0.9660829 -0.12202434 -0.09292027  0.044028801  0.029087566
# C0_Average 0.5781520  0.02313829  0.65550447  0.176680857 -0.153147148
var3$cor
# Dim.1       Dim.2       Dim.3        Dim.4        Dim.5
# C2         0.5705240  0.08878494  0.37582068 -0.053118510  0.033320711
# C3_1       0.3009570  0.57064471  0.08180654  0.032796640  0.011333119
# C3         0.2085785  0.73753659  0.07593569  0.232794714 -0.335729789
# C4_1       0.2379540  0.72997203  0.04912729  0.051624563  0.199114372
# C4         0.6373766  0.67197554 -0.19422267  0.050476018 -0.178784582
# C3OH       0.2179220  0.39171028  0.02956060  0.481756315  0.155236229
# C5_1       0.1775164  0.76344802  0.15266421 -0.100130877 -0.178147175
# C5         0.1625366  0.84783400 -0.19346253 -0.003250965 -0.230632304
# C4OH       0.5777071  0.21648618 -0.33068590  0.040779246  0.077351378
# C6_1       0.3934419  0.33736197 -0.33074540  0.023584426  0.555348344
# C6         0.8280668  0.14723359  0.10554197 -0.073230113  0.071957539
# C5OH       0.2543685  0.27775588  0.65795984  0.348783928 -0.012408795
# C5_1DC     0.5962648  0.61372650 -0.27836925 -0.075621009 -0.025804516
# C5DC       0.5291722  0.52478323  0.16481548 -0.415295790  0.117541594
# C8         0.8995066  0.05107445  0.12367972 -0.083121051 -0.072907328
# C5MDC      0.5990714 -0.08901571  0.04947379 -0.384805808 -0.035312872
# C9         0.9054296 -0.03169274  0.05601607 -0.021645116  0.001930235
# C7DC       0.1043711 -0.10517994  0.18944539  0.508499007  0.518179154
# C10_2      0.6732497 -0.28087183  0.33073596 -0.237729164 -0.005091980
# C10_1      0.3439148  0.32063008  0.31624985 -0.401809978  0.507216434
# C10        0.8870782 -0.13078922  0.09400225 -0.202488264  0.032614871
# C12_1      0.8816240  0.16800489  0.02127086 -0.246639375  0.150279670
# C12        0.9510426 -0.13135164 -0.06339824 -0.057992182 -0.068766026
# C14_2      0.9198261 -0.04024529 -0.10715097  0.174345223  0.024097144
# C14_1      0.9747282 -0.01164996 -0.02779031  0.034778870 -0.060256142
# C14        0.9726928 -0.07656936 -0.14065989  0.054703143 -0.058203822
# C12DC      0.8218828 -0.35523670  0.07553913 -0.010348409 -0.151062469
# C14_2OH    0.9518170 -0.25094079  0.02372481  0.042478929 -0.060289682
# C14_1OH    0.9704861 -0.08744574 -0.15424386  0.068004566 -0.014093342
# C16_2      0.9127627 -0.27736346  0.08784088 -0.001487240 -0.094833720
# C16_1      0.9606953 -0.20421987 -0.01361691  0.007395946 -0.041432932
# C16        0.9545146 -0.14582993 -0.14602655  0.072354595  0.044071892
# C16_2OH    0.9719150 -0.13192509 -0.08918920  0.073165236  0.043061782
# C16_1OH    0.9581388 -0.09600861 -0.14895232  0.140394266  0.042450711
# C16OH      0.9411718 -0.02094882 -0.13753202  0.039574013 -0.099310948
# C18_2      0.9659669 -0.17719586  0.02226547  0.066011625 -0.083412648
# C18_1      0.9702809 -0.19445164 -0.03290288  0.074561155 -0.005260420
# C18        0.9340043 -0.19364421 -0.11190754  0.131770336  0.082080666
# C18_1OH    0.9660829 -0.12202434 -0.09292027  0.044028801  0.029087566
# C0_Average 0.5781520  0.02313829  0.65550447  0.176680857 -0.153147148
var3$cos2
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# C2         0.32549769 0.0078827652 0.1412411799 2.821576e-03 1.110270e-03
# C3_1       0.09057512 0.3256353838 0.0066923095 1.075620e-03 1.284396e-04
# C3         0.04350497 0.5439602184 0.0057662286 5.419338e-02 1.127145e-01
# C4_1       0.05662211 0.5328591615 0.0024134909 2.665096e-03 3.964653e-02
# C4         0.40624888 0.4515511304 0.0377224469 2.547828e-03 3.196393e-02
# C3OH       0.04749001 0.1534369464 0.0008738289 2.320891e-01 2.409829e-02
# C5_1       0.03151207 0.5828528780 0.0233063615 1.002619e-02 3.173642e-02
# C5         0.02641816 0.7188224868 0.0374277489 1.056878e-05 5.319126e-02
# C4OH       0.33374549 0.0468662649 0.1093531664 1.662947e-03 5.983236e-03
# C6_1       0.15479656 0.1138130976 0.1093925229 5.562252e-04 3.084118e-01
# C6         0.68569454 0.0216777307 0.0111391079 5.362649e-03 5.177887e-03
# C5OH       0.06470334 0.0771483307 0.4329111487 1.216502e-01 1.539782e-04
# C5_1DC     0.35553176 0.3766602128 0.0774894393 5.718537e-03 6.658730e-04
# C5DC       0.28002318 0.2753974341 0.0271641411 1.724706e-01 1.381603e-02
# C8         0.80911204 0.0026085994 0.0152966736 6.909109e-03 5.315478e-03
# C5MDC      0.35888654 0.0079237973 0.0024476562 1.480755e-01 1.246999e-03
# C9         0.81980282 0.0010044295 0.0031377997 4.685110e-04 3.725806e-06
# C7DC       0.01089333 0.0110628206 0.0358895540 2.585712e-01 2.685096e-01
# C10_2      0.45326515 0.0788889847 0.1093862763 5.651516e-02 2.592826e-05
# C10_1      0.11827740 0.1028036460 0.1000139705 1.614513e-01 2.572685e-01
# C10        0.78690770 0.0171058198 0.0088364223 4.100150e-02 1.063730e-03
# C12_1      0.77726089 0.0282256432 0.0004524493 6.083098e-02 2.258398e-02
# C12        0.90448208 0.0172532543 0.0040193370 3.363093e-03 4.728766e-03
# C14_2      0.84608009 0.0016196836 0.0114813309 3.039626e-02 5.806723e-04
# C14_1      0.95009502 0.0001357217 0.0007723013 1.209570e-03 3.630803e-03
# C14        0.94613122 0.0058628666 0.0197852042 2.992434e-03 3.387685e-03
# C12DC      0.67549130 0.1261931100 0.0057061598 1.070896e-04 2.281987e-02
# C14_2OH    0.90595552 0.0629712817 0.0005628668 1.804459e-03 3.634846e-03
# C14_1OH    0.94184332 0.0076467574 0.0237911668 4.624621e-03 1.986223e-04
# C16_2      0.83313573 0.0769304884 0.0077160195 2.211883e-06 8.993434e-03
# C16_1      0.92293553 0.0417057553 0.0001854202 5.470002e-05 1.716688e-03
# C16        0.91109804 0.0212663686 0.0213237531 5.235187e-03 1.942332e-03
# C16_2OH    0.94461884 0.0174042298 0.0079547134 5.353152e-03 1.854317e-03
# C16_1OH    0.91802999 0.0092176531 0.0221867941 1.971055e-02 1.802063e-03
# C16OH      0.88580437 0.0004388531 0.0189150572 1.566103e-03 9.862664e-03
# C18_2      0.93309200 0.0313983731 0.0004957509 4.357535e-03 6.957670e-03
# C18_1      0.94144496 0.0378114408 0.0010825998 5.559366e-03 2.767202e-05
# C18        0.87236395 0.0374980804 0.0125232974 1.736342e-02 6.737236e-03
# C18_1OH    0.93331625 0.0148899406 0.0086341765 1.938535e-03 8.460865e-04
# C0_Average 0.33425974 0.0005353806 0.4296861050 3.121613e-02 2.345405e-02
var3$contrib
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# C2         1.43790450  0.158003956  7.452668333 1.901936e-01 8.620161e-02
# C3_1       0.40012071  6.527110371  0.353123380 7.250415e-02 9.972080e-03
# C3         0.19218569 10.903263464  0.304258215 3.653006e+00 8.751180e+00
# C4_1       0.25013137 10.680751330  0.127349171 1.796458e-01 3.078166e+00
# C4         1.79462748  9.050994494  1.990445602 1.717411e-01 2.481687e+00
# C3OH       0.20978983  3.075525368  0.046108063 1.564440e+01 1.870997e+00
# C5_1       0.13920637 11.682836858  1.229772948 6.758343e-01 2.464023e+00
# C5         0.11670373 14.408242903  1.974895701 7.124081e-04 4.129782e+00
# C4OH       1.47433965  0.939398170  5.770079807 1.120941e-01 4.645399e-01
# C6_1       0.68382258  2.281295850  5.772156476 3.749340e-02 2.394517e+01
# C6         3.02909452  0.434513410  0.587761138 3.614794e-01 4.020124e-01
# C5OH       0.28583066  1.546378848 22.842794225 8.200061e+00 1.195490e-02
# C5_1DC     1.57058172  7.549863754  4.088772771 3.854687e-01 5.169854e-02
# C5DC       1.23701829  5.520129377  1.433330806 1.162570e+01 1.072680e+00
# C8         3.57429830  0.052287365  0.807137376 4.657214e-01 4.126950e-01
# C5MDC      1.58540163  0.158826412  0.129151920 9.981307e+00 9.681730e-02
# C9         3.62152542  0.020133016  0.165567721 3.158086e-02 2.892724e-04
# C7DC       0.04812189  0.221745714  1.893732004 1.742948e+01 2.084715e+01
# C10_2      2.00232449  1.581268915  5.771826872 3.809510e+00 2.013076e-03
# C10_1      0.52249714  2.060619875  5.277292012 1.088292e+01 1.997439e+01
# C10        3.47620938  0.342872979  0.466258671 2.763783e+00 8.258824e-02
# C12_1      3.43359404  0.565761271  0.023873735 4.100426e+00 1.753426e+00
# C12        3.99560087  0.345828190  0.212082519 2.266956e-01 3.671426e-01
# C14_2      3.73760677  0.032465310  0.605818722 2.048917e+00 4.508354e-02
# C14_1      4.19709861  0.002720436  0.040750900 8.153332e-02 2.818964e-01
# C14        4.17958830  0.117516643  1.043977154 2.017106e-01 2.630207e-01
# C12DC      2.98402111  2.529443659  0.301088651 7.218573e-03 1.771740e+00
# C14_2OH    4.00210986  1.262210824  0.029699976 1.216330e-01 2.822103e-01
# C14_1OH    4.16064627  0.153273361  1.255353967 3.117312e-01 1.542108e-02
# C16_2      3.68042431  1.542012366  0.407140003 1.490961e-04 6.982524e-01
# C16_1      4.07712005  0.835959732  0.009783798 3.687157e-03 1.332841e-01
# C16        4.02482724  0.426267974  1.125159531 3.528876e-01 1.508031e-01
# C16_2OH    4.17290729  0.348854376  0.419734817 3.608392e-01 1.439696e-01
# C16_1OH    4.05544953  0.184760754  1.170698361 1.328627e+00 1.399126e-01
# C16OH      3.91309103  0.008796472  0.998063368 1.055661e-01 7.657396e-01
# C18_2      4.12198682  0.629356197  0.026158571 2.937278e-01 5.401952e-01
# C18_1      4.15888651  0.757901197  0.057123972 3.747395e-01 2.148463e-03
# C18        3.85371723  0.751620130  0.660798658 1.170414e+00 5.230806e-01
# C18_1OH    4.12297745  0.298457387  0.455587061 1.306706e-01 6.569036e-02
# C0_Average 1.47661137  0.010731292 22.672623026 2.104181e+00 1.820978e+00
#### graph of variables default plot ####
fviz_pca_var(carn_pca_FM, col.var = "black")
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(carn_pca_FM, col.var = "black", axes = c(2,3))
fviz_pca_var(carn_pca_FM, col.var = "black", axes = c(3,4))
fviz_pca_var(carn_pca_FM, col.var = "black", axes = c(4,5))

fviz_pca_var(carn_pca_FM, col.var = "black", axes = c(1,3))

#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var3$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(carn_pca_FM, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(carn_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Color by cos2 values: quality on the factor map
fviz_pca_var(carn_pca_FM, col.var = "cos2", axes = 2:3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

variables_PCA<- fviz_pca_var(carn_pca_FM, col.var = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                             repel = TRUE # Avoid text overlapping
)

# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(carn_pca_FM, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#### Color by a custom continuous variable ####
# Age

fviz_pca_ind(carn_pca_FM, col.ind = AA_data.frame$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

fviz_pca_var(carn_pca_FM, col.var = AA_data.frame$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

indyear1_2<- fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Spine.Age,
                          geom.ind = "point", # show points only (nbut not "text"),
                          pointsize = 2.5,
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          legend.title = "Age") 
ggsave("indyear1_2.pdf")

fviz_pca_ind(AA_pca_FM, col.ind = AA_data.frame$Year, axes = c(2,3),
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Year")  
indyear2_3<- fviz_pca_ind(AA_pca_FM, col.ind = DM_data.frame$Year, axes = c(2,3),
                          geom.ind = "point", # show points only (nbut not "text"),
                          pointsize = 2.5,
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          legend.title = "Year") 
ggsave("indyear2_3.pdf")

#### phospholipids correlation analysis ####


# Load the car package
library(car)
# Example data: three dependent variables (Y1, Y2, Y3), one independent variable (X), and one covariate (C)
set.seed(123)
Y1 <- rnorm(100)
Y2 <- rnorm(100)
Y3 <- rnorm(100)
X <- factor(rep(1:2, each = 50))
C <- rnorm(100)

#### load phospholipid data ####
Phosphatidylcholines <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Phosphatidylcholines.csv", quote="")
  View(Phosphatidylcholines)
dim(Phosphatidylcholines)
head(Phosphatidylcholines)
names(Phosphatidylcholines)
# [1] "ID"              "Site"            "TL"             
# [4] "FL"              "W"               "K"              
# [7] "Sex"             "Scale.Age"       "Spine.Age"      
# [10] "x14_1SMOH"       "x16_1SM"         "x16_0SM"        
# [13] "x16_1SMOH"       "x18_1SM"         "PC32_2AA"       
# [16] "x18_0SM"         "x20_2SM"         "PC36_0AE"       
# [19] "PC36_6AA"        "PC36_0AA"        "x22_2SMOH"      
# [22] "x22_1SMOH"       "PC38_6AA"        "PC38_0AA"       
# [25] "PC40_6AE"        "x24_1SMOH"       "PC40_6AA"       
# [28] "PC40_2AA"        "PC401AA"         "Choline_Average"
Phospha<- Phosphatidylcholines[1:54,]
dim(Phospha)
PTC <- Phosphatidylcholines[1:54,10:30]
dim(PTC)
#### Perform MANCOVA ####
manova_result3 <- manova(cbind(Phospha$x14_1SMOH, Phospha$x16_1SM,Phospha$x16_0SM,        
                               Phospha$x16_1SMOH, Phospha$x18_1SM, Phospha$PC32_2AA,       
                               Phospha$x18_0SM, Phospha$x20_2SM, Phospha$PC36_0AE,       
                               Phospha$PC36_6AA, Phospha$PC36_0AA, Phospha$x22_2SMOH,     
                               Phospha$x22_1SMOH, Phospha$PC38_6AA, Phospha$PC38_0AA,       
                               Phospha$PC40_6AE, Phospha$x24_1SMOH, Phospha$PC40_6AA,       
                               Phospha$PC40_2AA, Phospha$PC401AA, 
                               Phospha$Choline_Average) 
                         ~ Phospha$Spine.Age + Phospha$Site)
# Summarize MANOVA results
summary(manova_result3)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Phosphatidylcholines$Spine.Age  1 0.78549   5.0568     21     29 3.984e-05 ***
#   Phosphatidylcholines$Site       3 1.92453   2.6416     63     93 1.019e-05 ***
#   Residuals                      49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Extract coefficients
coefficients3 <- coef(manova_result3)
coefficients3
# [,1]        [,2]       [,3]       [,4]        [,5]      [,6]
# (Intercept)                         5.4837320 12.35677323 11.2918985 12.4495915  2.93101182 19.657064
# Phosphatidylcholines$Spine.Age     -0.1027777 -0.04076104 -0.1393313 -0.4154313  0.06574445  1.203095
# Phosphatidylcholines$SiteMatheson  -0.1145542 -4.20741974 -3.3278453 -3.3084787 -0.55113482  2.573721
# Phosphatidylcholines$SiteRed River -0.9239314 -5.88848468 -4.6434063 -3.8686466 -0.65041557  8.936929
# Phosphatidylcholines$SiteSandy Bar -0.7972585 -6.05906632 -4.9342340 -3.4732607 -0.46747482 11.258468
# [,7]        [,8]       [,9]       [,10]      [,11]      [,12]
# (Intercept)                        7.3994572  7.04966392 22.2075396  96.2826203 181.461143 15.1199485
# Phosphatidylcholines$Spine.Age     0.6061731  0.28428464 -0.4138631  -0.8892074  -5.439511 -0.4258161
# Phosphatidylcholines$SiteMatheson  0.1238613 -0.01928312 -6.8497899 -11.2735886 -67.442943 -3.0675887
# Phosphatidylcholines$SiteRed River 1.8799661  0.39688914 -6.9056914  -6.5633696 -70.528485 -4.3389286
# Phosphatidylcholines$SiteSandy Bar 3.6214668  0.05566067 -7.6388634  -9.2829486 -76.236547 -4.7498746
# [,13]       [,14]      [,15]     [,16]       [,17]      [,18]
# (Intercept)                         5.92141108 762.9657963  88.174284 81.300762 11.18775566  75.833360
# Phosphatidylcholines$Spine.Age     -0.08376491   0.1664258  -1.676043 -0.569254 -0.03926668   7.547012
# Phosphatidylcholines$SiteMatheson  -1.00924048 -86.4044130 -21.394855 -4.481130  0.56235558 -11.232418
# Phosphatidylcholines$SiteRed River -1.24924563 -85.6682760 -27.725050 -9.039347  0.92071338 -20.825211
# Phosphatidylcholines$SiteSandy Bar -1.72374503 -19.2165908 -28.657366 -5.111246  0.20537778 -23.951459
# [,19]       [,20]     [,21]
# (Intercept)                         3.067508559  2.55299741 20.619449
# Phosphatidylcholines$Spine.Age     -0.023131878  0.01628672  2.558771
# Phosphatidylcholines$SiteMatheson   0.007555546  0.12238702  7.314051
# Phosphatidylcholines$SiteRed River -0.181644658  0.11671582 17.707766
# Phosphatidylcholines$SiteSandy Bar -0.276429374 -0.21018669 32.926243
# Extract coefficients for the ith dependent variable
coeffs3 <- coefficients3[, i]

# Plot the relationship
ggplot(data = data.frame(X, Y = coeffs3[1] + coeffs3[2] * as.numeric(X) 
                         + coeffs[3] * C),
       aes(x = as.numeric(X), y = Y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "X", y = paste0("Y", i)) +
  ggtitle(paste0("Effect of X on Y", i))


# trying to plot
# Load necessary packages
library(tidyr)
library(dplyr)

# Reshape your data to long format
PC18_long <- Phospha %>%
  pivot_longer(cols = c(x14_1SMOH,x16_1SM,x16_0SM,       
                        x16_1SMOH,x18_1SM,PC32_2AA,      
                        x18_0SM,x20_2SM,PC36_0AE,      
                        PC36_6AA,PC36_0AA,x22_2SMOH,     
                        x22_1SMOH,PC38_6AA,PC38_0AA,      
                        PC40_6AE,x24_1SMOH,PC40_6AA,      
                        PC40_2AA,PC401AA,Choline_Average),
               names_to = "Phosphatidylcholines",
               values_to = "Concentration")


# trying to plot
# Load necessary packages
library(tidyr)
library(dplyr)

# Reshape your data to long format
Amino_Acids18_long <- Amino_Acids18 %>%
  pivot_longer(cols = c(Serine, Asparagine, Glutamine, Tryptophan, Lysine, Arginine_Average,
                        Glutamate_Average, Glycine_Average, Histidine_Average, Isoleucine_Average,
                        Leucine, Methionine, Phenylalanine_Average, Proline, Threonine_Average,
                        Tyrosine_Average, Alanine, Aspartate),
               names_to = "Amino_Acid",
               values_to = "Concentration")
# Plot using ggplot2
library(ggplot2)

ggplot(data = PC18_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~"Phosphatidylcholines", scales = "free_y") +
  labs(title = "Relationship between Spine Age, Site, and Choline Metabolite Concentration")

library(ggplot2)
library(tidyr)

# Assuming Amino_Acids18 is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets
ggplot(PC18_long, aes(x = Spine.Age, y = Concentration)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Phosphocholines, scales = "free") +
  theme_minimal() +
  labs(x = "Phosphocholines/Sphingomyelins", y = "Concentration")

library(ggplot2)
library(tidyr)

# Assuming Amino_Acids18_long is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets
ggplot(PC18_long, aes(x = Spine.Age, y = Concentration)) +
  geom_point() +  # Change to the appropriate geom for your data
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "Age", y = "Concentration")

library(ggplot2)
library(tidyr)

# Assuming Amino_Acids18_long is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets and color by location
ggplot(PC18_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +  # Change to the appropriate geom for your data
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "Age", y = "Concentration", color = "Location")

library(ggplot2)
library(tidyr)

# Assuming Amino_Acids18_long is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets, color by location, and regression lines
PC18PLOT <-ggplot(PC18_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +  # Change to the appropriate geom for your data
  geom_smooth(method = "lm", se = FALSE, aes(group = Site)) +  # Add regression lines
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "Age", y = "Concentration", color = "Site")
ggsave("PC18PLOT.pdf")
ggsave("PC18PLOT.png")

#### PCA ####      
#### load packages ####

library(tidyverse)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(sjstats)
library(performance)
library(forestmodel)
library(emmeans)
library(pwr)
# as data.frame
PC18_data.frame <- data.frame(PTC)
names(PC18_data.frame)
dim(PC18_data.frame)
# 54 21

# dataframe with just met data columns
PC_data_pca <- PC18_data.frame %>%
  print()

dim(PC_data_pca) 
#60 21
#### Correlations plot ####

PC_cor <- cor(PC_data_pca)

# summary(PC_cor)
# x14_1SMOH           x16_1SM           x16_0SM          x16_1SMOH           x18_1SM      
# Min.   :-0.02821   Min.   :-0.2374   Min.   :-0.2564   Min.   :-0.26024   Min.   :0.2225  
# 1st Qu.: 0.43118   1st Qu.: 0.3419   1st Qu.: 0.3185   1st Qu.: 0.06417   1st Qu.:0.4772  
# Median : 0.64722   Median : 0.4924   Median : 0.5139   Median : 0.41364   Median :0.6466  
# Mean   : 0.57109   Mean   : 0.4953   Mean   : 0.5087   Mean   : 0.39529   Mean   :0.5860  
# 3rd Qu.: 0.75354   3rd Qu.: 0.7753   3rd Qu.: 0.7944   3rd Qu.: 0.72001   3rd Qu.:0.6889  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000   Max.   :1.0000  
# PC32_2AA           x18_0SM            x20_2SM            PC36_0AE          PC36_6AA     
# Min.   :-0.19180   Min.   :-0.16756   Min.   :-0.01752   Min.   :-0.1639   Min.   :0.1702  
# 1st Qu.: 0.03204   1st Qu.: 0.06581   1st Qu.: 0.29942   1st Qu.: 0.3803   1st Qu.:0.5869  
# Median : 0.46180   Median : 0.47748   Median : 0.66578   Median : 0.6194   Median :0.6472  
# Mean   : 0.38087   Mean   : 0.39397   Mean   : 0.55773   Mean   : 0.5535   Mean   :0.6605  
# 3rd Qu.: 0.62346   3rd Qu.: 0.58687   3rd Qu.: 0.79799   3rd Qu.: 0.8219   3rd Qu.:0.7960  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000   Max.   :1.0000  
# PC36_0AA         x22_2SMOH         x22_1SMOH           PC38_6AA         PC38_0AA       
# Min.   :-0.2715   Min.   :-0.1789   Min.   :-0.01552   Min.   :0.1866   Min.   :-0.09387  
# 1st Qu.: 0.1599   1st Qu.: 0.4523   1st Qu.: 0.62032   1st Qu.:0.5720   1st Qu.: 0.59868  
# Median : 0.4403   Median : 0.6739   Median : 0.73279   Median :0.6647   Median : 0.77411  
# Mean   : 0.4228   Mean   : 0.5764   Mean   : 0.67731   Mean   :0.6541   Mean   : 0.65483  
# 3rd Qu.: 0.7377   3rd Qu.: 0.8430   3rd Qu.: 0.82133   3rd Qu.:0.7589   3rd Qu.: 0.80415  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000   Max.   :1.0000   Max.   : 1.00000  
# PC40_6AE        x24_1SMOH         PC40_6AA          PC40_2AA          PC401AA       
# Min.   :0.1820   Min.   :0.2011   Min.   :-0.1342   Min.   :0.08242   Min.   :0.06417  
# 1st Qu.:0.5305   1st Qu.:0.3803   1st Qu.: 0.2188   1st Qu.:0.49240   1st Qu.:0.38400  
# Median :0.6871   Median :0.5534   Median : 0.5097   Median :0.66578   Median :0.60305  
# Mean   :0.6662   Mean   :0.5521   Mean   : 0.4612   Mean   :0.63221   Mean   :0.58240  
# 3rd Qu.:0.8213   3rd Qu.:0.7247   3rd Qu.: 0.6963   3rd Qu.:0.81303   3rd Qu.:0.76328  
# Max.   :1.0000   Max.   :1.0000   Max.   : 1.0000   Max.   :1.00000   Max.   :1.00000  
# Choline_Average  
# Min.   :-0.2715  
# 1st Qu.:-0.1639  
# Median : 0.1702  
# Mean   : 0.1211  
# 3rd Qu.: 0.2833  
# Max.   : 1.0000   
Phoscorplot<- corrplot(PC_cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
ggsave("Phoscorplot.pdf")
ggsave("Phoscorplot.png")
dev.off()


#### PCA ####

# PCA with FactoMineR package

PC_pca_FM <- PCA(PC_data_pca, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(PC_pca_FM)
corrplot(PC_cor, type = "upper", order = "hclust", 
         +          tl.col = "black", tl.srt = 60, 
         +          tl.cex = 1, cl.cex = 1)
PC_pca_FM <- PCA(PC_data_pca, graph = FALSE)
# extract eigenvalues/varianes
get_eig(PC_pca_FM)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  11.834600261      56.35523934                    56.35524
# Dim.2   4.953314535      23.58721207                    79.94245
# Dim.3   0.917724297       4.37011570                    84.31257
# Dim.4   0.829004864       3.94764221                    88.26021
# Dim.5   0.648498030       3.08808586                    91.34830
# Dim.6   0.481235641       2.29159829                    93.63989
# Dim.7   0.313536480       1.49303086                    95.13292
# Dim.8   0.267859835       1.27552302                    96.40845
# Dim.9   0.172115848       0.81959927                    97.22805
# Dim.10  0.152244578       0.72497418                    97.95302
# Dim.11  0.121370608       0.57795527                    98.53098
# Dim.12  0.076112573       0.36244082                    98.89342
# Dim.13  0.055434190       0.26397233                    99.15739
# Dim.14  0.045478911       0.21656625                    99.37396
# Dim.15  0.040382007       0.19229527                    99.56625
# Dim.16  0.026976523       0.12845963                    99.69471
# Dim.17  0.022830143       0.10871497                    99.80343
# Dim.18  0.015123797       0.07201808                    99.87544
# Dim.19  0.013020203       0.06200097                    99.93744
# Dim.20  0.006978134       0.03322921                    99.97067
# Dim.21  0.006158542       0.02932639                   100.00000
fviz_screeplot(PC_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)

phos_screeplot<- fviz_screeplot(PC_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
ggsave("phos_screeplot.pdf")
ggsave("phos_screeplot.png")
fviz_screeplot(PC_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = 16
)

# extract the results for variables
var4 <- get_pca_var(PC_pca_FM)
var4
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
#### graph of variables default plot ####
var4$coord #  "Coordinates for the variables" 
# Dim.1      Dim.2         Dim.3        Dim.4        Dim.5
# x14_1SMOH       0.80480192 -0.2773314 -0.2247628019  0.224439133  0.191277426
# x16_1SM         0.71012983 -0.5375898  0.1126643933 -0.156664923  0.339808226
# x16_0SM         0.72931266 -0.5627289  0.1253872123 -0.139302207  0.268793506
# x16_1SMOH       0.57211852 -0.7034312  0.2052793167  0.230029953 -0.142350417
# x18_1SM         0.78067799  0.3293272  0.3563923349 -0.187924973  0.147326246
# PC32_2AA        0.47338240  0.7794562  0.2031160335  0.006980527 -0.157983950
# x18_0SM         0.48657574  0.7551973  0.3006474230  0.014221677  0.005236863
# x20_2SM         0.73092541  0.6233690  0.0870793884 -0.044199777  0.061621787
# PC36_0AE        0.78536108 -0.4901975  0.1970924582  0.011624554 -0.177737543
# PC36_6AA        0.90274258  0.2096062  0.0118858387  0.071107326 -0.253371482
# PC36_0AA        0.61164484 -0.6597036  0.2525549594  0.064248227 -0.256515885
# x22_2SMOH       0.82378223 -0.4729131 -0.0466967348  0.025920917  0.125184977
# x22_1SMOH       0.94756447 -0.1776770 -0.0502597293 -0.037770869  0.107579906
# PC38_6AA        0.89279132  0.2071670 -0.0477323986  0.085892702 -0.230196383
# PC38_0AA        0.92365314 -0.2985325 -0.0005335868  0.010872492 -0.080937789
# PC40_6AE        0.91604057  0.1815660 -0.1960517971  0.090834068 -0.068501003
# x24_1SMOH       0.75080341  0.2931032 -0.3663020872  0.235142990  0.076104593
# PC40_6AA        0.60485780  0.5838661  0.0945977471 -0.364016901  0.057386558
# PC40_2AA        0.87933530  0.1183945 -0.2963979714 -0.058184385 -0.083461708
# PC401AA         0.79722939  0.3969879 -0.3099576897 -0.203271310  0.025005543
# Choline_Average 0.09375146  0.6024240  0.2348702929  0.620144528  0.289303182
var4$cor    # "Correlations between variables and dimensions"
# Dim.1      Dim.2         Dim.3        Dim.4        Dim.5
# x14_1SMOH       0.80480192 -0.2773314 -0.2247628019  0.224439133  0.191277426
# x16_1SM         0.71012983 -0.5375898  0.1126643933 -0.156664923  0.339808226
# x16_0SM         0.72931266 -0.5627289  0.1253872123 -0.139302207  0.268793506
# x16_1SMOH       0.57211852 -0.7034312  0.2052793167  0.230029953 -0.142350417
# x18_1SM         0.78067799  0.3293272  0.3563923349 -0.187924973  0.147326246
# PC32_2AA        0.47338240  0.7794562  0.2031160335  0.006980527 -0.157983950
# x18_0SM         0.48657574  0.7551973  0.3006474230  0.014221677  0.005236863
# x20_2SM         0.73092541  0.6233690  0.0870793884 -0.044199777  0.061621787
# PC36_0AE        0.78536108 -0.4901975  0.1970924582  0.011624554 -0.177737543
# PC36_6AA        0.90274258  0.2096062  0.0118858387  0.071107326 -0.253371482
# PC36_0AA        0.61164484 -0.6597036  0.2525549594  0.064248227 -0.256515885
# x22_2SMOH       0.82378223 -0.4729131 -0.0466967348  0.025920917  0.125184977
# x22_1SMOH       0.94756447 -0.1776770 -0.0502597293 -0.037770869  0.107579906
# PC38_6AA        0.89279132  0.2071670 -0.0477323986  0.085892702 -0.230196383
# PC38_0AA        0.92365314 -0.2985325 -0.0005335868  0.010872492 -0.080937789
# PC40_6AE        0.91604057  0.1815660 -0.1960517971  0.090834068 -0.068501003
# x24_1SMOH       0.75080341  0.2931032 -0.3663020872  0.235142990  0.076104593
# PC40_6AA        0.60485780  0.5838661  0.0945977471 -0.364016901  0.057386558
# PC40_2AA        0.87933530  0.1183945 -0.2963979714 -0.058184385 -0.083461708
# PC401AA         0.79722939  0.3969879 -0.3099576897 -0.203271310  0.025005543
# Choline_Average 0.09375146  0.6024240  0.2348702929  0.620144528  0.289303182
var4$cos2   # "Cos2 for the variables" 
# Dim.1      Dim.2        Dim.3        Dim.4        Dim.5
# x14_1SMOH       0.647706125 0.07691270 5.051832e-02 5.037292e-02 3.658705e-02
# x16_1SM         0.504284369 0.28900275 1.269327e-02 2.454390e-02 1.154696e-01
# x16_0SM         0.531896951 0.31666384 1.572195e-02 1.940510e-02 7.224995e-02
# x16_1SMOH       0.327319598 0.49481539 4.213960e-02 5.291378e-02 2.026364e-02
# x18_1SM         0.609458120 0.10845644 1.270155e-01 3.531580e-02 2.170502e-02
# PC32_2AA        0.224090893 0.60755204 4.125612e-02 4.872776e-05 2.495893e-02
# x18_0SM         0.236755953 0.57032295 9.038887e-02 2.022561e-04 2.742474e-05
# x20_2SM         0.534251959 0.38858892 7.582820e-03 1.953620e-03 3.797245e-03
# PC36_0AE        0.616792022 0.24029361 3.884544e-02 1.351303e-04 3.159063e-02
# PC36_6AA        0.814944163 0.04393476 1.412732e-04 5.056252e-03 6.419711e-02
# PC36_0AA        0.374109405 0.43520883 6.378401e-02 4.127835e-03 6.580040e-02
# x22_2SMOH       0.678617164 0.22364684 2.180585e-03 6.718939e-04 1.567128e-02
# x22_1SMOH       0.897878423 0.03156910 2.526040e-03 1.426639e-03 1.157344e-02
# PC38_6AA        0.797076343 0.04291816 2.278382e-03 7.377556e-03 5.299037e-02
# PC38_0AA        0.853135131 0.08912164 2.847149e-07 1.182111e-04 6.550926e-03
# PC40_6AE        0.839130317 0.03296621 3.843631e-02 8.250828e-03 4.692387e-03
# x24_1SMOH       0.563705765 0.08590947 1.341772e-01 5.529223e-02 5.791909e-03
# PC40_6AA        0.365852954 0.34089963 8.948734e-03 1.325083e-01 3.293217e-03
# PC40_2AA        0.773230566 0.01401726 8.785176e-02 3.385423e-03 6.965857e-03
# PC401AA         0.635574704 0.15759937 9.607377e-02 4.131923e-02 6.252772e-04
# Choline_Average 0.008789336 0.36291465 5.516405e-02 3.845792e-01 8.369633e-02
var4$contrib #"contributions of the variables"
# Dim.1      Dim.2        Dim.3        Dim.4        Dim.5
# x14_1SMOH       5.47298693  1.5527523 5.504738e+00  6.076312269  5.641814147
# x16_1SM         4.26110183  5.8345326 1.383124e+00  2.960645861 17.805702616
# x16_0SM         4.49442262  6.3929685 1.713146e+00  2.340770933 11.141120780
# x16_1SMOH       2.76578499  9.9895815 4.591749e+00  6.382806846  3.124703601
# x18_1SM         5.14979895  2.1895730 1.384027e+01  4.260022694  3.346968180
# PC32_2AA        1.89352312 12.2655655 4.495481e+00  0.005877862  3.848728498
# x18_0SM         2.00054035 11.5139659 9.849240e+00  0.024397454  0.004228963
# x20_2SM         4.51432196  7.8450281 8.262634e-01  0.235658487  0.585544508
# PC36_0AE        5.21176895  4.8511680 4.232800e+00  0.016300297  4.871353925
# PC36_6AA        6.88611482  0.8869769 1.539386e-02  0.609918239  9.899352779
# PC36_0AA        3.16114948  8.7862142 6.950236e+00  0.497926470 10.146584269
# x22_2SMOH       5.73417901  4.5150946 2.376079e-01  0.081048249  2.416549903
# x22_1SMOH       7.58689270  0.6373328 2.752505e-01  0.172090490  1.784652490
# PC38_6AA        6.73513533  0.8664534 2.482643e-01  0.889929190  8.171246834
# PC38_0AA        7.20882085  1.7992325 3.102401e-05  0.014259395  1.010168946
# PC40_6AE        7.09048298  0.6655384 4.188219e+00  0.995268942  0.723577749
# x24_1SMOH       4.76320072  1.7343835 1.462065e+01  6.669710648  0.893126705
# PC40_6AA        3.09138412  6.8822527 9.751005e-01 15.984020136  0.507822211
# PC40_2AA        6.53364329  0.2829875 9.572783e+00  0.408371863  1.074152341
# PC401AA         5.37047885  3.1816952 1.046870e+01  4.984195776  0.096419286
# Choline_Average 0.07426813  7.3267031 6.010962e+00 46.390467900 12.906181268

#### graph of variables default plot ####

fviz_pca_var(PC_pca_FM, col.var = "black")
phosplot1<- fviz_pca_var(PC_pca_FM, col.var = "black")

ggsave("phosplot1.png")

# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(2,3))
phosplot2<- fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(2,3))
ggsave("phosplot2.png")

fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(3,4))
phosplot3<- fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(3,4))
ggsave("phosplot3.png")
fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(4,5))
phosplot4<-fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(4,5))
ggsave("phosplot4.png")

fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(1,3))
phosplot5<-fviz_pca_var(PC_pca_FM, col.var = "black", axes = c(1,3))
ggsave("phosplot5.png")

#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var4$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(PC_pca_FM, choice = "var", axes = 1:2)
phosbargraph<- fviz_cos2(PC_pca_FM, choice = "var", axes = 1:2)
ggsave("phosbargraoh.png")
# Color by cos2 values: quality on the factor map
fviz_pca_var(PC_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
phos_var_pca_1_2<- fviz_pca_var(PC_pca_FM, col.var = "cos2",
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                repel = TRUE # Avoid text overlapping
)
ggsave("phos_var_pca_1_2.png")
# Color by cos2 values: quality on the factor map
fviz_pca_var(PC_pca_FM, col.var = "cos2", axes = 2:3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
phos_var_pca_2_3<- fviz_pca_var(PC_pca_FM, col.var = "cos2", axes = 2:3,
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                repel = TRUE # Avoid text overlapping
)
ggsave("phos_var_pca_2_3.png")

phos_variables_PCA<- fviz_pca_var(PC_pca_FM, col.var = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                             repel = TRUE # Avoid text overlapping
)
ggsave("phos_variables_PCA.png")
# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(PC_pca_FM, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
phos_var_pva_transparency_cos2 <-fviz_pca_var(PC_pca_FM, alpha.var = "cos2",
                                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                              repel = TRUE # Avoid text overlapping
)
ggsave("phos_var_pva_transparency_cos2.png")
#### Color by a custom continuous variable ####
# Age

fviz_pca_ind(PC_pca_FM, col.ind = Phospha$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

fviz_pca_var(PC_pca_FM, col.var = Phospha$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

phos_indyear1_2PC<- fviz_pca_ind(PC_pca_FM, col.ind = Phospha$Spine.Age,
                          geom.ind = "point", # show points only (nbut not "text"),
                          pointsize = 2.5,
                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                          legend.title = "Age") 
ggsave("phos_indyear1_2PC.pdf")
indyear1_2PC


#### Correspondence analysis trying the aas ####
# normalize
PC_transformed <- decostand(PTC, "normalize")
PC.ca2 <- cca(PC_transformed)
PC.ca2
# Call: cca(X = PC_transformed)
# 
# Inertia Rank
# Total          0.0513     
# Unconstrained  0.0513   20
# Inertia is scaled Chi-square 
# 
# Eigenvalues for unconstrained axes:
#   CA1      CA2      CA3      CA4      CA5      CA6      CA7      CA8 
# 0.027014 0.011162 0.006012 0.002264 0.001925 0.001158 0.000559 0.000438 
# (Showing 8 of 20 unconstrained eigenvalues)
summary(PC.ca2)
# Call:
#   cca(X = PC_transformed) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          0.0513          1
# Unconstrained  0.0513          1
# 
# Eigenvalues, and their contribution to the scaled Chi-square 
# 
# Importance of components:
#   CA1     CA2      CA3      CA4      CA5      CA6
# Eigenvalue            0.02701 0.01116 0.006012 0.002264 0.001925 0.001158
# Proportion Explained  0.52664 0.21760 0.117213 0.044145 0.037534 0.022583
# Cumulative Proportion 0.52664 0.74424 0.861451 0.905596 0.943130 0.965713
# CA7       CA8       CA9      CA10      CA11
# Eigenvalue            0.0005587 0.0004384 0.0002518 0.0001388 0.0001056
# Proportion Explained  0.0108919 0.0085467 0.0049097 0.0027055 0.0020592
# Cumulative Proportion 0.9766054 0.9851520 0.9900618 0.9927673 0.9948265
# CA12      CA13      CA14      CA15      CA16
# Eigenvalue            7.867e-05 5.391e-05 3.739e-05 3.345e-05 2.461e-05
# Proportion Explained  1.534e-03 1.051e-03 7.289e-04 6.521e-04 4.798e-04
# Cumulative Proportion 9.964e-01 9.974e-01 9.981e-01 9.988e-01 9.993e-01
# CA17      CA18      CA19      CA20
# Eigenvalue            0.0000141 1.308e-05 5.642e-06 4.525e-06
# Proportion Explained  0.0002749 2.550e-04 1.100e-04 8.822e-05
# Cumulative Proportion 0.9995468 9.998e-01 9.999e-01 1.000e+00
# 
# Scaling 2 for species and site scores
# * Species are scaled proportional to eigenvalues
# * Sites are unscaled: weighted dispersion equal on all dimensions
# 
# 
# Species scores
# 
# CA1       CA2      CA3       CA4        CA5       CA6
# x14_1SMOH        0.071664  0.102116  0.04594 -0.058191 -0.0533139 -0.144387
# x16_1SM          0.180816  0.042619 -0.09672 -0.057986 -0.0416232 -0.188483
# x16_0SM          0.196964  0.046753 -0.08915 -0.026077 -0.0416429 -0.170801
# x16_1SMOH        0.329187  0.150929  0.02101  0.060744  0.0280037 -0.112222
# x18_1SM         -0.043305  0.019154 -0.07326 -0.058457  0.1101516 -0.073078
# PC32_2AA        -0.223353  0.006765  0.04813  0.012548  0.1693643 -0.057729
# x18_0SM         -0.206078  0.032715  0.01934  0.009909  0.1592140 -0.098136
# x20_2SM         -0.167069 -0.009353 -0.05579 -0.095352  0.0697235  0.032324
# PC36_0AE         0.230446  0.061871 -0.08235  0.007316  0.0184732 -0.028655
# PC36_6AA         0.008643 -0.011048  0.04881 -0.085203  0.0777691  0.057800
# PC36_0AA         0.439339  0.092764 -0.11112  0.045048  0.0341059  0.025996
# x22_2SMOH        0.198550  0.069743 -0.02130 -0.114477 -0.0351949 -0.111241
# x22_1SMOH        0.082236  0.028574 -0.04997 -0.113790 -0.0006123 -0.072045
# PC38_6AA        -0.005068 -0.026829  0.04720  0.024825 -0.0132885  0.001120
# PC38_0AA         0.183891 -0.006933 -0.03812 -0.109784 -0.0223591 -0.008160
# PC40_6AE        -0.009046 -0.022032  0.01895 -0.087208 -0.0557341  0.009785
# x24_1SMOH       -0.056362  0.017228  0.03331 -0.116740 -0.0452422 -0.030394
# PC40_6AA        -0.217400 -0.132892 -0.19193  0.012937 -0.0060288  0.007045
# PC40_2AA        -0.014194  0.027669 -0.01373 -0.045114 -0.0207313 -0.072133
# PC401AA         -0.087566 -0.008553 -0.03433 -0.065646 -0.0187449 -0.062635
# Choline_Average -0.316220  0.439517 -0.04397  0.005437 -0.0293489  0.019041
# 
# 
# Site scores (weighted averages of species scores)
# 
# CA1      CA2       CA3       CA4      CA5      CA6
# 1  -1.91792  1.60337 -0.729933 -0.076640  0.41975  1.22801
# 2   0.18885 -0.08485  1.748834  0.838881 -0.65623  1.46757
# 3   0.58414  0.48991  1.281660 -0.873884 -0.46778  0.19983
# 4  -1.61122 -1.13141 -1.182524  0.494167  1.83332 -0.26800
# 5  -1.34631 -0.71734 -0.610394  0.127517  0.22571  0.95237
# 6  -0.99739 -0.73226  1.153497  1.172254  1.40549  0.11159
# 7  -0.92267 -1.22472  0.421006  1.333493 -0.98208 -0.08463
# 8  -1.06632 -0.97863 -0.420225 -0.484010  0.98795 -0.07247
# 9  -1.51035 -0.49480 -0.553810 -0.245591  1.70086 -0.32893
# 10 -1.30995 -0.31500 -0.970779  0.164258 -1.11704 -1.27270
# 11 -1.59326  1.17217 -0.077747  0.082742  0.36981  0.42680
# 12 -0.12026  0.77545  0.751735 -1.633412 -0.24230  1.61555
# 13  0.88959  0.47566  1.049687 -1.177828  2.30854  0.89175
# 14  1.13533  0.66997  0.657335  1.346870  0.26267 -0.63157
# 15  0.67247 -0.08026  1.487998  0.748283 -0.64968 -0.16448
# 16 -1.32515 -0.54018 -0.001314  0.773847 -0.92789  0.78418
# 17 -0.11479  1.21013  2.026822  0.054537  2.05372 -2.10568
# 18  0.17506  1.03201  1.508137  1.267813 -1.55759  1.17026
# 19  0.84516  1.15885  0.934906  1.468378  0.81315 -1.50304
# 20  0.15671  0.12957  1.723771 -0.253737 -0.81766  1.06006
# 21 -1.17284  0.37046  0.760354  0.340198  1.62958 -0.06956
# 22 -1.59736  1.26175  0.045263  0.316907 -1.39970  1.47115
# 23  0.31333 -0.48525  1.411547 -1.081391  0.74763 -0.33083
# 24 -0.73387 -0.09915  0.426342  1.601874  1.33721  0.12816
# 25  0.12234 -1.22361 -0.024280 -0.217250  0.95718  0.83984
# 26  0.44873 -0.60589 -0.084060  2.488632 -0.42041 -0.23512
# 27  0.12092 -1.38115  0.361511 -0.618581  0.72176  1.71647
# 28  0.30006  0.05597  0.181719 -1.680984 -0.30096  0.51811
# 29 -1.03245  3.79935 -1.431838  0.073276 -0.42319 -0.34492
# 30  0.03863  2.40359  0.058751 -1.360808 -1.01688 -0.28536
# 31  0.13981  0.29075  0.735506 -2.403437  0.03979 -1.13257
# 32  0.30674 -0.42675  0.863904 -1.444961 -0.75545  0.75363
# 33 -1.08130  0.04593 -1.189254 -0.583800 -0.60250 -0.47171
# 34 -0.90513 -1.26069  0.206123  0.337337 -0.01185 -0.96522
# 35 -0.07795 -0.38214  1.369374  0.608294 -1.34084 -2.53894
# 36 -0.77795 -1.87099  0.145704 -0.369521 -1.35001  0.52554
# 37 -1.03448 -0.32102 -0.928437 -0.631524 -0.03787 -2.09731
# 38  0.38375 -0.29592 -0.975838  1.927638 -1.77223 -0.75967
# 39 -0.72593 -0.62029 -0.008393 -0.654831 -1.33087 -0.76466
# 40  1.53044  0.65449 -1.163683  1.406087  0.23106  0.80772
# 41  1.00336 -0.42420 -0.348725 -0.008408 -0.02825  0.24551
# 42  1.07463 -0.32622 -0.010119  0.266037  1.11139  0.38852
# 43  0.75552 -1.04922  0.411378 -1.885496 -1.32879 -1.65193
# 44  1.08267 -0.42322 -0.372132  0.914876  0.02516  0.13150
# 45  1.29811 -0.03953  0.519168  0.378027  0.11706  0.70033
# 46  1.87279  0.54226 -0.728340  0.600780 -0.29215 -1.10361
# 47  0.50676 -0.32469  0.460439 -0.052043 -0.88815 -0.03583
# 48  0.35598 -0.57278 -1.625693  0.068056 -0.50720  0.42222
# 49  1.18120  0.70320 -0.334454  0.622360  0.19011  0.34250
# 50  0.49909 -0.70928 -0.855429 -1.064422 -0.94935  0.71138
# 51  0.53015 -0.74536 -0.864333 -0.835740 -0.06000 -1.67582
# 52  1.15817 -0.74380 -1.356282  0.013413  0.17545  1.74715
# 53  1.52727  0.39408 -2.607032 -0.506871  0.34955  0.29398
# 54  1.47034  0.16709 -1.338548 -0.520343  1.64071 -0.64285
plot(PC.ca2, display = "species")
phos_ca_species<-plot(PC.ca2, display = "species")
ggsave("phos_ca_species.png")
plot(PC.ca2, type="n")
text(PC.ca2,display="species")
points(PC.ca2,col="blue")

#### lysophospholipids analysis ####
Lysophospholipids <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Lysophospholipids.csv", quote="")
View(Lysophospholipids)
names(Lysophospholipids)
column_names <- names(Lysophospholipids)

# Print the column names
print(column_names)
Lysophospholipids$LYSOC26_1
#### Perform MANCOVA lyso ####
manova_result4 <- manova(cbind(Lysophospholipids$LYSOC14_0,
                               Lysophospholipids$LYSOC16_1, 
                               Lysophospholipids$LYSOC16_0,
                               Lysophospholipids$LYSOC17_0, 
                               Lysophospholipids$LYSOC18_2, 
                               Lysophospholipids$LYSOC18_1, 
                               Lysophospholipids$LYSOC18_0, 
                               Lysophospholipids$LYSOC20_4,
                               Lysophospholipids$LYSOC20_3,
                               Lysophospholipids$LYSOC24_0,
                               Lysophospholipids$LYSOC26_1,
                               Lysophospholipids$LYSOC26_0, 
                               Lysophospholipids$LYSOC28_1,
                               Lysophospholipids$LYSOC28_0) 
                         ~ Lysophospholipids$Spine.Age + Lysophospholipids$Site)
# Summarize MANOVA results
summary(manova_result4)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Lysophospholipids$Spine.Age  1 0.60591   3.9536     14     36 0.0004352 ***
#   Lysophospholipids$Site       3 1.10864   1.5910     42    114 0.0279532 *  
#   Residuals                   49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Extract coefficients
coefficients4 <- coef(manova_result4)
coefficients4

# Reshape your data to long format
lyso18_long <- Lysophospholipids %>%
  pivot_longer(cols = c(LYSOC14_0,
                        LYSOC16_1, 
                        LYSOC16_0,
                        LYSOC17_0, 
                        LYSOC18_2, 
                        LYSOC18_1, 
                        LYSOC18_0, 
                        LYSOC20_4,
                        LYSOC20_3,
                        LYSOC24_0,
                        LYSOC26_1,
                        LYSOC26_0, 
                        LYSOC28_1,
                        LYSOC28_0),
               names_to = "Lysophospholipids",
               values_to = "Concentration")
# Plot using ggplot2
library(ggplot2)


# Example plot using ggplot2 with facets, color by location, and regression lines
lys18PLOT <-ggplot(lyso18_long, aes(x = Spine.Age, y = Concentration, color = Site)) +
  geom_point() +  # Change to the appropriate geom for your data
  geom_smooth(method = "lm", se = FALSE, aes(group = Site)) +  # Add regression lines
  facet_wrap(~ Lysophospholipids, scales = "free") +
  theme_minimal() +
  labs(x = "Age", y = "Concentration", color = "Site")
ggsave("lys18PLOT.pdf")
ggsave("lys18PLOT.png")

lys18PLOT


#### PCA lyso ####
dim(Lysophospholipids)
lyso18_data.frame <- data.frame(Lysophospholipids[,10:23])
names(lyso18_data.frame)
dim(lyso18_data.frame)
# 54 21

# dataframe with just met data columns
lyso_data_pca <- lyso18_data.frame %>%
  print()

dim(lyso_data_pca) 
#60 21
#### Correlations plot ####

lyso_cor <- cor(lyso_data_pca)

lyso_cor
# LYSOC14_0  LYSOC16_1  LYSOC16_0 LYSOC17_0 LYSOC18_2 LYSOC18_1 LYSOC18_0
# LYSOC14_0 1.0000000 0.83199900 0.80631287 0.8539243 0.6511887 0.6886374 0.4827952
# LYSOC16_1 0.8319990 1.00000000 0.85340667 0.8763429 0.6031389 0.5958808 0.2958889
# LYSOC16_0 0.8063129 0.85340667 1.00000000 0.8526515 0.4740136 0.5026367 0.3109810
# LYSOC17_0 0.8539243 0.87634292 0.85265146 1.0000000 0.5761475 0.5642660 0.3137862
# LYSOC18_2 0.6511887 0.60313886 0.47401365 0.5761475 1.0000000 0.9111146 0.7328276
# LYSOC18_1 0.6886374 0.59588075 0.50263670 0.5642660 0.9111146 1.0000000 0.8505661
# LYSOC18_0 0.4827952 0.29588890 0.31098098 0.3137862 0.7328276 0.8505661 1.0000000
# LYSOC20_4 0.3211490 0.27218028 0.06593617 0.1778867 0.6111966 0.5401914 0.4401687
# LYSOC20_3 0.4742875 0.20743674 0.20349122 0.4304893 0.3645590 0.3773631 0.4308913
# LYSOC24_0 0.3961212 0.17124488 0.27149075 0.4432091 0.2353449 0.2593819 0.2761937
# LYSOC26_1 0.4591675 0.24143734 0.17144091 0.4331945 0.4932677 0.4696037 0.4124363
# LYSOC26_0 0.4966495 0.25480990 0.33155135 0.5247856 0.3524790 0.3553251 0.3021140
# LYSOC28_1 0.5370578 0.29137335 0.22118310 0.5065735 0.4596140 0.4218488 0.3678261
# LYSOC28_0 0.3632407 0.08897578 0.10560702 0.3699150 0.3205842 0.2771307 0.2575217
# LYSOC20_4 LYSOC20_3  LYSOC24_0 LYSOC26_1  LYSOC26_0 LYSOC28_1  LYSOC28_0
# LYSOC14_0  0.32114902 0.4742875  0.3961212 0.4591675 0.49664948 0.5370578 0.36324072
# LYSOC16_1  0.27218028 0.2074367  0.1712449 0.2414373 0.25480990 0.2913733 0.08897578
# LYSOC16_0  0.06593617 0.2034912  0.2714908 0.1714409 0.33155135 0.2211831 0.10560702
# LYSOC17_0  0.17788665 0.4304893  0.4432091 0.4331945 0.52478564 0.5065735 0.36991502
# LYSOC18_2  0.61119660 0.3645590  0.2353449 0.4932677 0.35247895 0.4596140 0.32058421
# LYSOC18_1  0.54019142 0.3773631  0.2593819 0.4696037 0.35532511 0.4218488 0.27713067
# LYSOC18_0  0.44016874 0.4308913  0.2761937 0.4124363 0.30211396 0.3678261 0.25752173
# LYSOC20_4  1.00000000 0.1535794 -0.0575427 0.2441072 0.01166202 0.4055096 0.12209730
# LYSOC20_3  0.15357941 1.0000000  0.3812340 0.4266790 0.36543936 0.5292829 0.38967041
# LYSOC24_0 -0.05754270 0.3812340  1.0000000 0.8420170 0.92678182 0.7112483 0.89075105
# LYSOC26_1  0.24410716 0.4266790  0.8420170 1.0000000 0.85608245 0.8306397 0.89395293
# LYSOC26_0  0.01166202 0.3654394  0.9267818 0.8560825 1.00000000 0.7331524 0.90321136
# LYSOC28_1  0.40550960 0.5292829  0.7112483 0.8306397 0.73315240 1.0000000 0.81146116
# LYSOC28_0  0.12209730 0.3896704  0.8907511 0.8939529 0.90321136 0.8114612 1.00000000

 
corrplot(lyso_cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
lyscorplot<- corrplot(lyso_cor, type = "upper", order = "hclust", 
                      tl.col = "black", tl.srt = 60, 
                      tl.cex = 1, cl.cex = 1)

save("lyscorplot.png")
#### PCA ####
  
# PCA with FactoMineR package

lyso_pca_FM <- PCA(lyso_data_pca, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(lyso_pca_FM)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  7.14889322       51.0635230                    51.06352
# Dim.2  2.77840440       19.8457457                    70.90927
# Dim.3  1.70889494       12.2063924                    83.11566
# Dim.4  0.78597547        5.6141105                    88.72977
# Dim.5  0.69595573        4.9711124                    93.70088
# Dim.6  0.23592607        1.6851862                    95.38607
# Dim.7  0.15118231        1.0798736                    96.46594
# Dim.8  0.13017787        0.9298420                    97.39579
# Dim.9  0.11075976        0.7911411                    98.18693
# Dim.10 0.07618688        0.5441920                    98.73112
# Dim.11 0.05547233        0.3962310                    99.12735
# Dim.12 0.04925216        0.3518011                    99.47915
# Dim.13 0.04319047        0.3085033                    99.78765
# Dim.14 0.02972840        0.2123457                   100.00000

fviz_screeplot(lyso_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)

lyso_scree-plot<- fviz_screeplot(lyso_pca_FM, addlabels = TRUE, 
                                 ggtheme = theme_classic(),
                                 main = "",
                                 font.x = c(14, "bold"), font.y = c(14, "bold"),
                                 font.tickslab = 12,
                                 barfill = "#99d8c9", barcolor = "#66c2a4",
                                 font.submain = 16,
                                 ncp = Inf
)
save("lyso_scree-plot.png")
fviz_screeplot(lyso_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = 16
)
lyso_scree_plot2<- fviz_screeplot(lyso_pca_FM, addlabels = TRUE, 
                                  ggtheme = theme_classic(),
                                  main = "",
                                  font.x = c(14, "bold"), font.y = c(14, "bold"),
                                  font.tickslab = 12,
                                  barfill = "#99d8c9", barcolor = "#66c2a4",
                                  font.submain = 16,
                                  ncp = 16
)
save("lyso_scree_plot2.png")
# extract the results for variables
var5 <- get_pca_var(lyso_pca_FM)
var5
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
#### graph of variables default plot ####
var5$coord #

# Dim.1       Dim.2       Dim.3        Dim.4       Dim.5
# LYSOC14_0 0.8468572  0.33935160 -0.23685025  0.043887106  0.10455004
# LYSOC16_1 0.6668464  0.57029604 -0.37847133 -0.150257131  0.08803224
# LYSOC16_0 0.6297442  0.48886419 -0.53052475 -0.013131915 -0.10089706
# LYSOC17_0 0.8079102  0.28794630 -0.43900719  0.005935978  0.12659597
# LYSOC18_2 0.7753506  0.37967360  0.36101647 -0.086873548 -0.11332538
# LYSOC18_1 0.7794901  0.40329775  0.35271867  0.023622357 -0.25394375
# LYSOC18_0 0.6399814  0.24008570  0.50927166  0.251208963 -0.38307216
# LYSOC20_4 0.4021662  0.31050571  0.63918177 -0.358247237  0.38930991
# LYSOC20_3 0.5707935 -0.09987329  0.11578878  0.721242533  0.34196295
# LYSOC24_0 0.6999679 -0.62414801 -0.18850482 -0.010542943 -0.16393503
# LYSOC26_1 0.7958002 -0.50836234  0.12025422 -0.121300355 -0.03430450
# LYSOC26_0 0.7682127 -0.54172618 -0.18036153 -0.082075831 -0.15953503
# LYSOC28_1 0.7926125 -0.39299807  0.12866765 -0.083449665  0.32401905
# LYSOC28_0 0.7001532 -0.67034884  0.01517816 -0.114001903  0.00145308
var5$cor
# Dim.1       Dim.2       Dim.3        Dim.4       Dim.5
# LYSOC14_0 0.8468572  0.33935160 -0.23685025  0.043887106  0.10455004
# LYSOC16_1 0.6668464  0.57029604 -0.37847133 -0.150257131  0.08803224
# LYSOC16_0 0.6297442  0.48886419 -0.53052475 -0.013131915 -0.10089706
# LYSOC17_0 0.8079102  0.28794630 -0.43900719  0.005935978  0.12659597
# LYSOC18_2 0.7753506  0.37967360  0.36101647 -0.086873548 -0.11332538
# LYSOC18_1 0.7794901  0.40329775  0.35271867  0.023622357 -0.25394375
# LYSOC18_0 0.6399814  0.24008570  0.50927166  0.251208963 -0.38307216
# LYSOC20_4 0.4021662  0.31050571  0.63918177 -0.358247237  0.38930991
# LYSOC20_3 0.5707935 -0.09987329  0.11578878  0.721242533  0.34196295
# LYSOC24_0 0.6999679 -0.62414801 -0.18850482 -0.010542943 -0.16393503
# LYSOC26_1 0.7958002 -0.50836234  0.12025422 -0.121300355 -0.03430450
# LYSOC26_0 0.7682127 -0.54172618 -0.18036153 -0.082075831 -0.15953503
# LYSOC28_1 0.7926125 -0.39299807  0.12866765 -0.083449665  0.32401905
# LYSOC28_0 0.7001532 -0.67034884  0.01517816 -0.114001903  0.00145308
var5$cos2
# Dim.1       Dim.2        Dim.3        Dim.4        Dim.5
# LYSOC14_0 0.7171672 0.115159510 0.0560980424 1.926078e-03 1.093071e-02
# LYSOC16_1 0.4446841 0.325237578 0.1432405459 2.257721e-02 7.749676e-03
# LYSOC16_0 0.3965778 0.238988199 0.2814565065 1.724472e-04 1.018022e-02
# LYSOC17_0 0.6527189 0.082913071 0.1927273128 3.523584e-05 1.602654e-02
# LYSOC18_2 0.6011686 0.144152040 0.1303328901 7.547013e-03 1.284264e-02
# LYSOC18_1 0.6076049 0.162649078 0.1244104632 5.580157e-04 6.448743e-02
# LYSOC18_0 0.4095762 0.057641146 0.2593576281 6.310594e-02 1.467443e-01
# LYSOC20_4 0.1617376 0.096413798 0.4085533394 1.283411e-01 1.515622e-01
# LYSOC20_3 0.3258052 0.009974674 0.0134070415 5.201908e-01 1.169387e-01
# LYSOC24_0 0.4899550 0.389560737 0.0355340685 1.111536e-04 2.687469e-02
# LYSOC26_1 0.6332980 0.258432265 0.0144610781 1.471378e-02 1.176799e-03
# LYSOC26_0 0.5901507 0.293467255 0.0325302820 6.736442e-03 2.545143e-02
# LYSOC28_1 0.6282346 0.154447481 0.0165553635 6.963847e-03 1.049883e-01
# LYSOC28_0 0.4902145 0.449367571 0.0002303765 1.299643e-02 2.111440e-06
var5$contrib
# Dim.1      Dim.2       Dim.3        Dim.4        Dim.5
# LYSOC14_0 10.031863  4.1448074  3.28270868  0.245055746 1.570604e+00
# LYSOC16_1  6.220321 11.7059121  8.38205689  2.872507655 1.113530e+00
# LYSOC16_0  5.547401  8.6016348 16.47008837  0.021940530 1.462768e+00
# LYSOC17_0  9.130349  2.9841974 11.27789126  0.004483071 2.302810e+00
# LYSOC18_2  8.409254  5.1883031  7.62673510  0.960209790 1.845324e+00
# LYSOC18_1  8.499286  5.8540462  7.28017038  0.070996584 9.266024e+00
# LYSOC18_0  5.729225  2.0746132 15.17692061  8.028996585 2.108529e+01
# LYSOC20_4  2.262415  3.4701139 23.90745798 16.328891734 2.177756e+01
# LYSOC20_3  4.557422  0.3590073  0.78454451 66.184100442 1.680260e+01
# LYSOC24_0  6.853579 14.0210236  2.07935945  0.014142126 3.861552e+00
# LYSOC26_1  8.858686  9.3014633  0.84622394  1.872040129 1.690911e-01
# LYSOC26_0  8.255134 10.5624384  1.90358583  0.857080455 3.657047e+00
# LYSOC28_1  8.787858  5.5588553  0.96877597  0.886013217 1.508549e+01
# LYSOC28_0  6.857208 16.1735840  0.01348102  1.653541937 3.033871e-04

#### graph of variables default plot lyso ####

fviz_pca_var(lyso_pca_FM, col.var = "black")
lyso1<- fviz_pca_var(lyso_pca_FM, col.var = "black")
save("lyso1.png")
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(lyso_pca_FM, col.var = "black", axes = c(2,3))
lyso2<- fviz_pca_var(lyso_pca_FM, col.var = "black",axes = c(2,3))
save("lyso2.png")

fviz_pca_var(lyso_pca_FM, col.var = "black", axes = c(3,4))
lyso3<- fviz_pca_var(lyso_pca_FM, col.var = "black",axes = c(3,4))
save("lyso3.png")

fviz_pca_var(lyso_pca_FM, col.var = "black", axes = c(4,5))
lyso4<- fviz_pca_var(lyso_pca_FM, col.var = "black",axes = c(4,5))
save("lyso4.png")

fviz_pca_var(lyso_pca_FM, col.var = "black", axes = c(1,3))

#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.
# contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(lyso_pca_FM, choice = "var", axes = 1:2)

# Color by cos2 values: quality on the factor map
fviz_pca_var(lyso_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Color by cos2 values: quality on the factor map
fviz_pca_var(lyso_pca_FM, col.var = "cos2", axes = 2:3,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

variables5_PCA<- fviz_pca_var(lyso_pca_FM, col.var = "cos2",
                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                             repel = TRUE # Avoid text overlapping
)
variables5_PCA
# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(lyso_pca_FM, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#### Color by a custom continuous variable ####
# Age

fviz_pca_ind(lyso_pca_FM, col.ind = Phospha$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

#### correlation analysis FL ####

#### Compute correlations FL ####
correlationsFL <- cor(Amino_Acids18$FL, num_amino)

# Print correlations
print(correlationsFL)
# Serine Asparagine Glutamine Tryptophan    Lysine Arginine_Average
# [1,] 0.6009596 -0.2010011 0.2919843 0.05031348 0.4947476        0.5490666
# Glutamate_Average Glycine_Average Histidine_Average Isoleucine_Average   Leucine
# [1,]         0.4020331       0.5075164        0.09505531         0.07558305 0.2867163
# Methionine Phenylalanine_Average    Proline Threonine_Average Tyrosine_Average
# [1,]  0.4448718             0.3942967 -0.1049786         0.3475866        0.3590063
# Valine_Average    Alanine  Aspartate
# [1,]     0.09183021 0.09910765 0.02349854
library(car)

#### Amino Acids Perform MANCOVA ####
manova_resultFL<- manova(cbind(Amino_Acids18$Serine, Amino_Acids18$Asparagine, 
                               Amino_Acids18$Glutamine, Amino_Acids18$Tryptophan,
                               Amino_Acids18$Lysine, Amino_Acids18$Arginine_Average,
                               Amino_Acids18$Glutamate_Average, Amino_Acids18$Glycine_Average,
                               Amino_Acids18$Histidine_Average, Amino_Acids18$Isoleucine_Average,
                               Amino_Acids18$Leucine, Amino_Acids18$Methionine,
                               Amino_Acids18$Phenylalanine_Average, Amino_Acids18$Proline,
                               Amino_Acids18$Threonine_Average, Amino_Acids18$Tyrosine_Average,
                               Amino_Acids18$Alanine, Amino_Acids18$Aspartate)  ~ Amino_Acids18$FL + Amino_Acids18$Site)
# Summarize MANOVA results
summary(manova_resultFL)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Amino_Acids18$FL    1 0.76359   5.7420     18     32 9.431e-06 ***
#   Amino_Acids18$Site  3 2.04177   4.0248     54    102 7.185e-10 ***
#   Residuals          49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       
# Extract coefficients
coefficientsFL <- coef(manova_resultFL)
coefficientsFL
# [,1]        [,2]         [,3]         [,4]
# (Intercept)                 -43.6241934  52.4979511   3.51804929  24.75283465
# Amino_Acids18$FL              0.2862336  -0.0244573   0.04874009   0.01095957
# Amino_Acids18$SiteMatheson  -11.2847907 -31.5006787 -14.78700554 -16.25039640
# Amino_Acids18$SiteRed River  27.2777404 -39.6157104  -5.32475237 -12.65883700
# Amino_Acids18$SiteSandy Bar  11.9427174 -38.1795205 -14.35201450 -17.42921905
# [,5]        [,6]        [,7]         [,8]
# (Intercept)                 -55.6601814 -46.9533185  30.0342686  102.9881285
# Amino_Acids18$FL              0.1366509   0.1957346   0.1703018    0.8460872
# Amino_Acids18$SiteMatheson   11.5095440  13.4593751 -42.0038120 -141.3039492
# Amino_Acids18$SiteRed River  22.1089255  24.7032885 -17.1541491  -23.3735761
# Amino_Acids18$SiteSandy Bar  34.6824663  22.2148580 -35.1879230  -85.4320951
# [,9]        [,10]       [,11]        [,12]
# (Intercept)                  49.60517012  111.8914692   7.0642843  -1.02568914
# Amino_Acids18$FL              0.01468924    0.1135251   0.3214508   0.07580897
# Amino_Acids18$SiteMatheson   -6.01572435 -100.4155649 -64.3293080 -14.66091673
# Amino_Acids18$SiteRed River  -5.91740007  -98.2064625 -61.2101658  -6.06777715
# Amino_Acids18$SiteSandy Bar -22.72897378  -89.8288667 -17.4009140 -11.33122770
# [,13]         [,14]       [,15]       [,16]
# (Intercept)                  26.0222869 132.286394856  34.0703429 -12.7998695
# Amino_Acids18$FL              0.1210469  -0.001647761   0.1956744   0.2012749
# Amino_Acids18$SiteMatheson  -38.4459509 -80.370514035   2.2055339 -50.3354917
# Amino_Acids18$SiteRed River -14.6240545 -91.173115729 -20.0070204 -22.4373358
# Amino_Acids18$SiteSandy Bar -29.8251574 -83.380810845 -45.9305677 -47.5243237
# [,17]         [,18]
# (Intercept)                 229.4432636  17.815022246
# Amino_Acids18$FL              0.0982835   0.007115284
# Amino_Acids18$SiteMatheson   -9.4499077  -8.168612056
# Amino_Acids18$SiteRed River  15.9995017 -12.701520588
# Amino_Acids18$SiteSandy Bar -48.7276855 -15.038607549
# Reshape your data to long format
Amino_Acids18_long <- Amino_Acids18 %>%
pivot_longer(cols = c(Serine, Asparagine, Glutamine, Tryptophan, Lysine, Arginine_Average,
                                Glutamate_Average, Glycine_Average, Histidine_Average, Isoleucine_Average,
                                Leucine, Methionine, Phenylalanine_Average, Proline, Threonine_Average,
                                Tyrosine_Average, Alanine, Aspartate),
                       names_to = "Amino_Acid",
                       values_to = "Concentration")
        
# Plot using ggplot2
        
        
ggplot(data = Amino_Acids18_long, aes(x = FL, y = Concentration, color = Site)) +
          geom_point() +
          geom_smooth(method = "lm", se = FALSE) +
          facet_wrap(~Amino_Acid, scales = "free_y") +
          labs(title = "Relationship between Fork Length, Site, and Amino Acid Concentration")
        
        
        

#### FL Compute correlations Carnitines ####
correlationsCFL <- cor(Carnitine_book$FL, Carnitine_only)

# Print correlations
print(correlationsCFL)
#C2        C3_1         C3       C4_1         C4       C3OH
# [1,] 0.3799869 -0.06969024 -0.1331719 -0.2279919 -0.1450002 -0.1455589
# C5_1         C5       C4OH       C6_1          C6      C5OH
# [1,] -0.0530221 -0.2077195 -0.1364356 -0.3601744 -0.05203432 0.1071554
# C5_1DC       C5DC        C8    C5MDC        C9       C7DC     C10_2
# [1,] -0.1517193 0.07146401 0.1213146 0.222465 0.2195445 0.05760186 0.2874442
# C10_1       C10     C12_1       C12       C14_2     C14_1        C14
# [1,] 0.2004984 0.2764728 0.1374464 0.1040009 -0.07105483 0.1261562 0.07742114
# C12DC   C14_2OH    C14_1OH     C16_2     C16_1        C16    C16_2OH
# [1,] 0.3121473 0.2123942 0.04542568 0.2802805 0.1602683 0.04062809 0.03341544
# C16_1OH      C16OH     C18_2     C18_1        C18    C18_1OH
# [1,] -0.007551621 0.08857781 0.2000688 0.1181948 0.03008532 0.08353118
# C0_Average
# [1,]   0.376303
#### FL Perform MANCOVA ####
manova_resultCFL <- manova(cbind(Carnitine_book$C2, Carnitine_book$C3_1,
                                 Carnitine_book$C3OH, Carnitine_book$C4,
                                 Carnitine_book$C4_1, Carnitine_book$C4OH,
                                 Carnitine_book$C5, Carnitine_book$C5_1,
                                 Carnitine_book$C5_1DC, Carnitine_book$C5DC,
                                 Carnitine_book$C5MDC, Carnitine_book$C5OH,
                                 Carnitine_book$C6, Carnitine_book$C6_1,
                                 Carnitine_book$C7DC, Carnitine_book$C8,
                                 Carnitine_book$C9, Carnitine_book$C10,
                                 Carnitine_book$C10_1, Carnitine_book$C10_2,
                                 Carnitine_book$C12, Carnitine_book$C12_1,
                                 Carnitine_book$C12DC, Carnitine_book$C14,
                                 Carnitine_book$C14_1, Carnitine_book$C14_1OH,
                                 Carnitine_book$C14_2, Carnitine_book$C14_2OH,
                                 Carnitine_book$C16, Carnitine_book$C16_1,
                                 Carnitine_book$C16OH, Carnitine_book$C16_1OH,
                                 Carnitine_book$C16_2, Carnitine_book$C16_2OH,
                                 Carnitine_book$C18, Carnitine_book$C18_1,
                                 Carnitine_book$C18_1OH, Carnitine_book$C18_2) 
                           ~ Carnitine_book$FL + Carnitine_book$Site)
# Summarize MANOVA results
summary(manova_resultCFL)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Carnitine_book$FL    1 0.90854   3.1371     38     12   0.01858 *  
#   Carnitine_book$Site  3 2.66908   2.9715    114     42 6.509e-05 ***
#   Residuals           49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Extract coefficients
coefficientsCFL<- coef(manova_resultCFL)
coefficientsCFL
# (Intercept)                   3.11409836  6.395491e-02  2.968882e-02
# Carnitine_book$FL             0.02459665 -1.445118e-05 -6.927176e-06
# Carnitine_book$SiteMatheson   5.76355730  2.511856e-04 -2.990941e-03
# Carnitine_book$SiteRed River  1.41090129 -1.177455e-02 -2.058922e-03
# Carnitine_book$SiteSandy Bar -4.97054142 -2.400434e-02 -2.286730e-03
# [,4]          [,5]          [,6]
# (Intercept)                   0.4526491117  1.014895e-01  0.1335967710
# Carnitine_book$FL            -0.0001605573 -3.351205e-05 -0.0001076586
# Carnitine_book$SiteMatheson  -0.1746775730 -1.388520e-02  0.0094999101
# Carnitine_book$SiteRed River -0.1986733584 -3.143770e-02 -0.0340822064
# Carnitine_book$SiteSandy Bar -0.2709790491 -2.690571e-02 -0.0513074539
# [,7]          [,8]          [,9]
# (Intercept)                   0.6580953025  2.761189e-02  1.567515e-01
# Carnitine_book$FL            -0.0003184458  1.908853e-06 -3.739802e-05
# Carnitine_book$SiteMatheson  -0.3056033908 -3.343119e-03 -3.951671e-02
# Carnitine_book$SiteRed River -0.3841740187 -1.003582e-02 -4.405595e-02
# Carnitine_book$SiteSandy Bar -0.4115503506 -7.328733e-03 -4.465761e-02
# [,10]         [,11]         [,12]
# (Intercept)                   1.952118e-01  2.500649e-02  1.521940e-01
# Carnitine_book$FL             7.727344e-05  1.083277e-05  6.427841e-05
# Carnitine_book$SiteMatheson   1.484482e-02 -1.444870e-03  8.494527e-03
# Carnitine_book$SiteRed River -6.856579e-02  5.095280e-04 -1.183344e-02
# Carnitine_book$SiteSandy Bar -5.347171e-02  6.054405e-04 -1.536354e-02
# [,13]         [,14]        [,15]
# (Intercept)                   0.9106808689  0.2687023511 1.977883e-02
# Carnitine_book$FL            -0.0002226476 -0.0001272331 1.531916e-05
# Carnitine_book$SiteMatheson   0.0279755462 -0.0067728261 1.778011e-02
# Carnitine_book$SiteRed River  0.0416235745 -0.0351862327 1.624431e-02
# Carnitine_book$SiteSandy Bar -0.1190555984 -0.0138502249 1.489714e-02
# [,16]         [,17]         [,18]
# (Intercept)                   0.1877167502  1.728387e-02  1.238959e-01
# Carnitine_book$FL             0.0001008617  1.628951e-05  2.901642e-05
# Carnitine_book$SiteMatheson  -0.0171151586 -5.053064e-04 -2.046143e-03
# Carnitine_book$SiteRed River  0.0397634117  2.335600e-03  3.963376e-03
# Carnitine_book$SiteSandy Bar -0.0219427999 -2.618221e-03  6.060490e-05
# [,19]        [,20]        [,21]
# (Intercept)                   1.781624e-01 2.634177e-02  0.028931212
# Carnitine_book$FL             4.773887e-05 9.196658e-05  0.000012028
# Carnitine_book$SiteMatheson  -7.638565e-03 1.537231e-02 -0.001844393
# Carnitine_book$SiteRed River -1.858982e-02 3.283661e-02  0.006956194
# Carnitine_book$SiteSandy Bar -2.238042e-02 2.676767e-02 -0.002390558
# [,22]        [,23]         [,24]
# (Intercept)                   5.001968e-02 5.927875e-03  4.105702e-02
# Carnitine_book$FL             2.980728e-05 2.174634e-05  4.683923e-05
# Carnitine_book$SiteMatheson  -7.738961e-03 7.609211e-04 -1.548882e-02
# Carnitine_book$SiteRed River -9.756326e-03 6.635000e-03  1.114867e-02
# Carnitine_book$SiteSandy Bar -1.602270e-02 3.813674e-03 -1.623123e-02
# [,25]         [,26]         [,27]
# (Intercept)                   6.823689e-02  3.389946e-02  2.087134e-02
# Carnitine_book$FL             8.449364e-05  1.401459e-05 -6.876409e-06
# Carnitine_book$SiteMatheson  -2.226726e-02 -8.990083e-03 -1.239287e-03
# Carnitine_book$SiteRed River  3.258419e-03  8.283134e-03  2.491862e-03
# Carnitine_book$SiteSandy Bar -3.063624e-02 -7.567630e-03 -1.793038e-03
# [,28]        [,29]        [,30]
# (Intercept)                  4.980083e-03 0.1106775124  0.041090337
# Carnitine_book$FL            3.380027e-05 0.0000233096  0.000192383
# Carnitine_book$SiteMatheson  4.423376e-04 0.0059825625 -0.012261464
# Carnitine_book$SiteRed River 1.076989e-02 0.0470165247  0.094909041
# Carnitine_book$SiteSandy Bar 3.923300e-03 0.0080758323  0.008123951
# [,31]         [,32]         [,33]
# (Intercept)                   9.845719e-03  0.0465083241  5.132890e-03
# Carnitine_book$FL             6.181926e-06 -0.0000100231  1.893311e-05
# Carnitine_book$SiteMatheson  -1.497749e-03 -0.0034019027 -8.127214e-04
# Carnitine_book$SiteRed River  3.630007e-04  0.0112231598  5.622891e-03
# Carnitine_book$SiteSandy Bar -2.698537e-03  0.0009558367  3.823840e-05
# [,34]        [,35]        [,36]
# (Intercept)                   2.317941e-02 3.781172e-02 0.0854542712
# Carnitine_book$FL             1.794953e-06 3.688939e-06 0.0001528147
# Carnitine_book$SiteMatheson  -1.119353e-03 7.646858e-04 0.0042469262
# Carnitine_book$SiteRed River  9.176560e-03 1.769482e-02 0.1116620338
# Carnitine_book$SiteSandy Bar  9.904621e-04 1.235940e-02 0.0173821130
# [,37]         [,38]
# (Intercept)                   1.607737e-02  5.587607e-03
# Carnitine_book$FL             7.460680e-06  4.224858e-05
# Carnitine_book$SiteMatheson  -4.003315e-05 -2.503687e-03
# Carnitine_book$SiteRed River  5.185474e-03  1.328579e-02
# Carnitine_book$SiteSandy Bar -3.739445e-04 -1.876999e-03
# Reshape your data to long format
Carnitines_long <- Carnitine_book %>%
  pivot_longer(cols = c(C2, C3_1,C3OH, C4, C4_1, C4OH, C5, C5_1,
                        C5_1DC, C5DC,
                        C5MDC, C5OH,
                        C6, C6_1,
                        C7DC, C8,
                        C9, C10,
                        C10_1, C10_2,
                        C12, C12_1,
                        C12DC, C14,
                        C14_1, C14_1OH,
                        C14_2, C14_2OH,
                        C16, C16_1,
                        C16OH, C16_1OH,
                        C16_2, C16_2OH,
                        C18, C18_1,
                        C18_1OH, C18_2),
               names_to = "Carnitines",
               values_to = "Concentration")

# Plot using ggplot2


ggplot(data = Carnitines_long, aes(x = FL, y = Concentration, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Carnitines, scales = "free_y") +
  labs(title = "Relationship between FL, Site, and Carnitine Concentration")

#### phospholipids correlation analysis ####

#### load phospholipid and lyso data ####
Phosphatidylcholines <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Phosphatidylcholines.csv", quote="")
View(Phosphatidylcholines)
dim(Phosphatidylcholines)
head(Phosphatidylcholines)
names(Phosphatidylcholines)
# [1] "ID"              "Site"            "TL"             
# [4] "FL"              "W"               "K"              
# [7] "Sex"             "Scale.Age"       "Spine.Age"      
# [10] "x14_1SMOH"       "x16_1SM"         "x16_0SM"        
# [13] "x16_1SMOH"       "x18_1SM"         "PC32_2AA"       
# [16] "x18_0SM"         "x20_2SM"         "PC36_0AE"       
# [19] "PC36_6AA"        "PC36_0AA"        "x22_2SMOH"      
# [22] "x22_1SMOH"       "PC38_6AA"        "PC38_0AA"       
# [25] "PC40_6AE"        "x24_1SMOH"       "PC40_6AA"       
# [28] "PC40_2AA"        "PC401AA"         "Choline_Average"
Phospha<- Phosphatidylcholines[1:54,]
dim(Phospha)
PTC <- Phosphatidylcholines[1:54,10:30]
dim(PTC)
#### Perform MANCOVA ####
manova_resultPHFL <- manova(cbind(Phospha$x14_1SMOH, Phospha$x16_1SM,Phospha$x16_0SM,        
                               Phospha$x16_1SMOH, Phospha$x18_1SM, Phospha$PC32_2AA,       
                               Phospha$x18_0SM, Phospha$x20_2SM, Phospha$PC36_0AE,       
                               Phospha$PC36_6AA, Phospha$PC36_0AA, Phospha$x22_2SMOH,     
                               Phospha$x22_1SMOH, Phospha$PC38_6AA, Phospha$PC38_0AA,       
                               Phospha$PC40_6AE, Phospha$x24_1SMOH, Phospha$PC40_6AA,       
                               Phospha$PC40_2AA, Phospha$PC401AA, 
                               Phospha$Choline_Average) 
                         ~ Phospha$FL + Phospha$Site)
# Summarize MANOVA results
summary(manova_resultPHFL)
# Df  Pillai approx F num Df den Df    Pr(>F)    
# Phospha$FL    1 0.86112   8.5624     21     29 1.620e-07 ***
#   Phospha$Site  3 1.95867   2.7766     63     93 3.793e-06 ***
#   Residuals    49                                             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Extract coefficients
coefficientsPHFL <- coef(manova_resultPHFL)
coefficientsPHFL
# [,1]         [,2]         [,3]        [,4]
# (Intercept)            6.450334562 12.933422945 12.643404654 17.06090865
# Phospha$FL            -0.003400073 -0.001754771 -0.004695778 -0.01522362
# Phospha$SiteMatheson  -0.159717598 -4.211895580 -3.386212950 -3.44207907
# Phospha$SiteRed River -1.082136353 -5.930830508 -4.853538625 -4.43380110
# Phospha$SiteSandy Bar -0.855544958 -6.085053713 -5.013861274 -3.71931800
# [,5]        [,6]       [,7]        [,8]
# (Intercept)            2.593759836 14.40530713 4.22892485 4.769383998
# Phospha$FL             0.001584151  0.02705583 0.01473436 0.008577808
# Phospha$SiteMatheson  -0.502709264  3.52382060 0.56611034 0.132980472
# Phospha$SiteRed River -0.519557720 11.42863074 3.08005489 0.875994023
# Phospha$SiteSandy Bar -0.434365332 11.85069452 3.92764729 0.211038734
# [,9]        [,10]       [,11]      [,12]
# (Intercept)           27.33758296 108.37402971 244.4326312 19.0157292
# Phospha$FL            -0.01629311  -0.03725418  -0.2047822 -0.0138578
# Phospha$SiteMatheson  -6.94562111 -11.40516837 -69.0120615 -3.2622755
# Phospha$SiteRed River -7.41213913  -7.53867338 -77.6548529 -5.0058786
# Phospha$SiteSandy Bar -7.89195577  -9.84261400 -79.4968454 -4.9897418
# [,13]        [,14]        [,15]        [,16]
# (Intercept)            6.759643164 789.49733183 104.24919896  84.37465110
# Phospha$FL            -0.002877123  -0.05355417  -0.05610255  -0.01403972
# Phospha$SiteMatheson  -1.042543143 -84.37836868 -22.10966240  -4.88973981
# Phospha$SiteRed River -1.372861775 -82.44729386 -30.27203874 -10.15616588
# Phospha$SiteSandy Bar -1.771998334 -19.53956504 -29.61250467  -5.40021182
# [,17]      [,18]         [,19]         [,20]
# (Intercept)           10.653495374  10.172983  3.3164675570  2.3288134746
# Phospha$FL             0.000600276   0.238491 -0.0008312659  0.0006880602
# Phospha$SiteMatheson   0.482297541  -7.546436 -0.0004261414  0.1246081414
# Phospha$SiteRed River  0.764926146  -8.647021 -0.2139371192  0.1342927092
# Phospha$SiteSandy Bar  0.196530839 -19.750453 -0.2900142789 -0.1998954854
# [,21]
# (Intercept)           -0.7028417
# Phospha$FL             0.0788842
# Phospha$SiteMatheson   8.6290626
# Phospha$SiteRed River 21.9358467
# Phospha$SiteSandy Bar 34.3366150

# Reshape your data to long format
PC18_long <- Phospha %>%
  pivot_longer(cols = c(x14_1SMOH,x16_1SM,x16_0SM,       
                        x16_1SMOH,x18_1SM,PC32_2AA,      
                        x18_0SM,x20_2SM,PC36_0AE,      
                        PC36_6AA,PC36_0AA,x22_2SMOH,     
                        x22_1SMOH,PC38_6AA,PC38_0AA,      
                        PC40_6AE,x24_1SMOH,PC40_6AA,      
                        PC40_2AA,PC401AA,Choline_Average),
               names_to = "Phosphatidylcholines",
               values_to = "Concentration")


# trying to plot


ggplot(data = PC18_long, aes(x = FL, y = Concentration, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~"Phosphatidylcholines", scales = "free_y") +
  labs(title = "Relationship between Fork Length, Site, and Choline Metabolite Concentration")


# Assuming is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets
ggplot(PC18_long, aes(x = FL, y = Concentration)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "Phosphocholines/Sphingomyelins", y = "Concentration")

# Assuming Amino_Acids18_long is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets
ggplot(PC18_long, aes(x = FL, y = Concentration)) +
  geom_point() +  # Change to the appropriate geom for your data
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "FL", y = "Concentration")



# Assuming Amino_Acids18_long is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets and color by location
ggplot(PC18_long, aes(x = FL, y = Concentration, color = Site)) +
  geom_point() +  # Change to the appropriate geom for your data
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "FL", y = "Concentration", color = "Location")


# Assuming is your dataset after reshaping to long format
# Check if your dataset is indeed in long format
# If not, reshape it using pivot_longer as you've described in your code

# Example plot using ggplot2 with facets, color by location, and regression lines
PC18PLOTFL <-ggplot(PC18_long, aes(x = FL, y = Concentration, color = Site)) +
  geom_point() +  # Change to the appropriate geom for your data
  geom_smooth(method = "lm", se = FALSE, aes(group = Site)) +  # Add regression lines
  facet_wrap(~ Phosphatidylcholines, scales = "free") +
  theme_minimal() +
  labs(x = "FL", y = "Concentration", color = "Site")
ggsave("PC18PLOTFL.pdf")
ggsave("PC18PLOT.png")

#### Walleye PCR ####
#### Mike's Data ####
Walleye_PCR_17 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Walleye_PCR_17.csv", quote="")
Walleye_PCR_18 <- read.csv("C:/Users/user/Desktop/Chapter 2 Analysis/Walleye_PCR_18.csv", quote="")

summary(Walleye_PCR_17)
#    ID            Site            Total_Length    Fork_Length   
# Min.   :  1.0   Length:67          Min.   :482.0   Min.   :452.0  
# 1st Qu.: 19.5   Class :character   1st Qu.:556.5   1st Qu.:522.5  
# Median :161.0   Mode  :character   Median :627.0   Median :594.0  
# Mean   :127.1                      Mean   :620.4   Mean   :587.8  
# 3rd Qu.:187.5                      3rd Qu.:673.5   3rd Qu.:639.0  
# Max.   :205.0                      Max.   :755.0   Max.   :713.0  
# Weight_Kg           K              Tel_Ct          Ef1_Ct     
# Min.   :1.000   Min.   :0.7413   Min.   :12.06   Min.   :18.25  
# 1st Qu.:1.700   1st Qu.:0.8797   1st Qu.:12.85   1st Qu.:18.66  
# Median :2.310   Median :0.9904   Median :13.16   Median :18.84  
# Mean   :2.506   Mean   :0.9956   Mean   :13.13   Mean   :18.81  
# 3rd Qu.:3.225   3rd Qu.:1.1018   3rd Qu.:13.44   3rd Qu.:18.94  
# Max.   :4.750   Max.   :1.3233   Max.   :13.94   Max.   :19.45  
# PanX2_Ct     EF1a_TS_Ratio    Panx2_TS_Ratio   Mean_TS_Ratio   
# Min.   :22.06   Min.   :0.5778   Min.   :0.5497   Min.   :0.5638  
# 1st Qu.:22.42   1st Qu.:0.8646   1st Qu.:0.8512   1st Qu.:0.8639  
# Median :22.58   Median :1.0674   Median :1.0682   Median :1.0442  
# Mean   :22.56   Mean   :1.1916   Mean   :1.1549   Mean   :1.1732  
# 3rd Qu.:22.70   3rd Qu.:1.3762   3rd Qu.:1.3274   3rd Qu.:1.3267  
# Max.   :23.14   Max.   :2.6925   Max.   :2.5920   Max.   :2.6169  
# CV          
# Min.   : 0.02726  
# 1st Qu.: 2.12217  
# Median : 4.23317  
# Mean   : 5.18851  
# 3rd Qu.: 7.65196  
# Max.   :18.21896  

one.way1 <- aov(Mean_TS_Ratio ~ Fork_Length, data = Walleye_PCR_17)

summary(one.way1)
#Df Sum Sq Mean Sq F value Pr(>F)
# Fork_Length  1  0.168  0.1683   0.825  0.367
# Residuals   65 13.261  0.2040 


attach(Walleye_PCR_17)
# Create a scatterplot to visualize the relationship
plot(Fork_Length, Mean_TS_Ratio, pch=16)
plot(Walleye_PCR_17$Fork_Length, Walleye_PCR_17$Mean_TS_Ratio, 
     col=as.factor (Walleye_PCR_17$Site)) 
library(ggpubr)
# Scatter plot by group
ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
stat_cor(method = "pearson", label.x = 650, label.y = 2.5) +
  theme_classic()

ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = 650, label.y = 2.5) +
  scale_color_brewer(palette = "Set1") +  # Change the color palette
  theme_classic()

ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
  stat_cor(method = "pearson", label.x = 650, label.y = 2.5) +
  theme_classic() +
  theme(plot.background = element_rect(fill = "white"))  # Change the plot background to white

# Calculate the correlation coefficient
correlation_coefficient17 <- cor(Walleye_PCR_17$Fork_Length, Walleye_PCR_17$Mean_TS_Ratio)

# Create the scatterplot
ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
  geom_text(x = 650, y = 2.5, label = paste("R =", round(correlation_coefficient17, 3))) +
  labs(y = "Mean T/S Ratio", x = "Fork Length (mm)",
       title = "Relationship between Mean T/S Ratio and Fork-Length Walleye 2017") +
  theme_classic()


# another code

library(ggpubr)  # Load the ggpubr package

ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x = 630) +  # Add correlation coefficient
  labs(y = "Mean T/S Ratio", x = "Fork-Length (mm)",
       title = "Mean T/S Ratio vs. Fork_Length Walleye 2017") +
  #scale_color_brewer(palette = "Set1") +  # Change the color palette
  theme_classic()

ggplot(Walleye_PCR_17, aes(x = Fork_Length, y = Mean_TS_Ratio, color = Site)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x = 630) + # Add correlation coefficient
  scale_color_brewer(palette = "Set1") +
  labs(y = "Mean T/S Ratio", x = "Fork-Length (mm)",
       title = "Mean T/S Ratio vs. Fork_Length Walleye 2017") +
  #scale_color_brewer(palette = "Set1") +  # Change the color palette
  theme_classic()


ggplot(Walleye_PCR_17, aes(x = Site, y = Mean_TS_Ratio, fill= Site)) +
  geom_bar(stat = "identity")
# Perform the correlation test
cor.test(Fork_Length, Mean_TS_Ratio)
# Pearson's product-moment correlation
# 
# data:  Fork_Length and Mean_TS_Ratio
# t = -0.90832, df = 65, p-value = 0.3671
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3429415  0.1317975
# sample estimates:
#        cor 
# -0.1119555 
# Print the results
print(correlation_result)

one.way2 <- aov(Mean_TS_Ratio ~ Site, data = Walleye_PCR_17)

summary(one.way2)
# Df Sum Sq Mean Sq F value Pr(>F)
# Site         4  0.499  0.1247   0.598  0.665
# Residuals   62 12.930  0.2086  

ggplot(Walleye_PCR_17, aes(x=Site, y=Mean_TS_Ratio)) + 
  geom_boxplot()

p <- ggplot(Walleye_PCR_17, aes(x=Site, y=Mean_TS_Ratio, fill = Site)) + 
  geom_boxplot() + theme_classic()
p

#### 2018 Jen data ####
summary(Walleye_PCR_18)
# ID            Site               Sex                 Age        
# Min.   :206.0   Length:80          Length:80          Min.   : 3.000  
# 1st Qu.:237.8   Class :character   Class :character   1st Qu.: 5.000  
# Median :279.5   Mode  :character   Mode  :character   Median : 7.000  
# Mean   :274.4                                         Mean   : 7.414  
# 3rd Qu.:320.2                                         3rd Qu.: 8.750  
# Max.   :355.0                                         Max.   :15.000  
# NA's   :10      
#   Total_length    Fork_length        Weight      Condition_factor
#  Min.   :373.0   Min.   :350.0   Min.   :0.500   Min.   :0.7391  
#  1st Qu.:432.2   1st Qu.:404.5   1st Qu.:0.750   1st Qu.:0.9261  
#  Median :480.0   Median :450.0   Median :1.100   Median :1.0150  
#  Mean   :518.9   Mean   :490.9   Mean   :1.691   Mean   :1.0217  
#  3rd Qu.:617.8   3rd Qu.:592.5   3rd Qu.:2.513   3rd Qu.:1.0798  
#  Max.   :740.0   Max.   :703.0   Max.   :4.850   Max.   :1.4816  
#                                                                  
#  EF1a_TS_ratio    PanX2_TS_ratio   Mean_TS_ratio   
#  Min.   :0.6004   Min.   :0.6490   Min.   :0.6369  
#  1st Qu.:0.8670   1st Qu.:0.8606   1st Qu.:0.8836  
#  Median :1.0316   Median :1.0602   Median :1.0624  
#  Mean   :1.2103   Mean   :1.1877   Mean   :1.1990  
#  3rd Qu.:1.2933   3rd Qu.:1.3134   3rd Qu.:1.3364  
#  Max.   :9.8965   Max.   :8.6905   Max.   :9.2935  

one.way3 <- aov(Mean_TS_ratio ~ Fork_length, data = Walleye_PCR_18)

summary(one.way3)
# Df Sum Sq Mean Sq F value Pr(>F)
# Fork_length  1  0.134  0.1340   1.802  0.183
# Residuals   77  5.729  0.0744               
# 1 observation deleted due to missingness


attach(Walleye_PCR_18)
# Create a scatterplot to visualize the relationship
plot(Fork_length, Mean_TS_ratio, pch=16)
plot(Walleye_PCR_18$Fork_length, Walleye_PCR_18$Mean_TS_ratio, 
     col=as.factor (Walleye_PCR_18$Site)) 

# Scatter plot by group
ggplot(Walleye_PCR_18, aes(x = Fork_length, y = Mean_TS_ratio, color = Site)) +
  geom_point() +theme_classic()


 ggplot(Walleye_PCR_18, aes(x = Site, y = Mean_TS_ratio, fill= Site)) +
  geom_bar(stat = "identity") + theme_classic() 

# Perform the correlation test
cor.test(Fork_length, Mean_TS_ratio)
# Pearson's product-moment correlation
# 
# data:  Fork_length and Mean_TS_ratio
# t = -1.3422, df = 77, p-value = 0.1835
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.36026554  0.07232903
# sample estimates:
#        cor 
# -0.1511995 
print(correlation_result)

one.way4 <- aov(Mean_TS_ratio ~ Site, data = Walleye_PCR_18)

summary(one.way4)
# Df Sum Sq Mean Sq F value Pr(>F)
# Site         3  0.235 0.07845   1.046  0.377
# Residuals   75  5.627 0.07503               
# 1 observation deleted due to missingness

ggplot(Walleye_PCR_18, aes(x=Site, y=Mean_TS_ratio)) + 
  geom_boxplot()

p2 <- ggplot(Walleye_PCR_18, aes(x=Site, y=Mean_TS_ratio, fill = Site)) + 
  geom_boxplot() + theme_classic()
p2

cor.test(Age, Mean_TS_ratio)
# Pearson's product-moment correlation
# 
# data:  Age and Mean_TS_ratio
# t = -1.658, df = 67, p-value = 0.102
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.41566998  0.04004223
# sample estimates:
#        cor 
# -0.1985199 

# Scatter plot by group
ggplot(Walleye_PCR_18, aes(x = Age, y = Mean_TS_ratio, color = Site)) +
  geom_point() + theme_classic()
ggplot(Walleye_PCR_18, aes(x = Fork_length, y = Mean_TS_ratio, color = Site)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x = 630) + # Add correlation coefficient
  scale_color_brewer(palette = "Set1") +
  labs(y = "Mean T/S Ratio", x = "Fork-Length (mm)",
       title = "Mean T/S Ratio vs. Fork Length Walleye 2018") +
  #scale_color_brewer(palette = "Set1") +  # Change the color palette
  theme_classic()

# bar plot
ggplot(Walleye_PCR_18, aes(x = Site, y = Mean_TS_ratio, fill= Site)) +
  geom_bar(stat = "identity")

ggplot(Walleye_PCR_17, aes(x=Site, y=Mean_TS_Ratio)) + 
  geom_boxplot()

p3 <- ggplot(Walleye_PCR_18, aes(x=Site, y=Mean_TS_ratio, fill = Site)) + 
  geom_boxplot() +  theme_classic()
p3

# age graph

# Scatter plot by group
ggplot(Walleye_PCR_18, aes(x = Age, y = Mean_TS_ratio, color = Site)) +
  geom_point() + theme_classic()


ggplot(Walleye_PCR_18, aes(x = Age, y = Mean_TS_ratio, color = Site)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x = 10) + # Add correlation coefficient
  scale_color_brewer(palette = "Set1") +
  labs(y = "Mean T/S Ratio", x = "Age (yr)",
       title = "Mean T/S Ratio vs. Fork Length Walleye 2018") +
  #scale_color_brewer(palette = "Set1") +  # Change the color palette
  theme_classic()
library(ggpubr)
