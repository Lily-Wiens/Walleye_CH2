#### Chapter 2 rerun of PCA ####
#no asparagine no DR outlier
#------------------------------------------#
#              Walleye ProjectII           #
#------------------------------------------#


#### 2018 Metabolite Only with no DR outlier  PCA ####
#### set up ####
getwd()
setwd("C:/Users/user/Desktop/Chapter 2 Analysis")
# for using at Jay's PC
setwd("C:/Users/Treberg's Lab PC/Desktop/Chapter 2 Analysis")


options(stringsAsFactors = FALSE) #prevent character data from being read as factor data

ls()
rm(list = ls())

options(max.print=1000000)

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
library(ggplot2)


#### Read in 2018 data missing DR outlier ####

# read data
FM_MET_18NDRout <- read.csv("C:/Users/Treberg's Lab PC/Desktop/Chapter 2 Analysis/FM_MET_18NDRout.csv")
View(FM_MET_18NDRout)

# CODE NOT USED CARRIED OVER FROM LAST WHOLE PCA
# upload data again
# X2018_MET_FM <- read.csv("C:/Users/user/Desktop/Regression/2018_MET_FM.csv")
# View(X2018_MET_FM)   
# raw_data<- X2018_MET_FM
# load in object with met_groups
met_groups <- readxl::read_excel("Class_2018.xlsx") %>% 
  print()
# met_groups <- met_groups[-1,]
# Read in the same raw data in Excel format like Matt the whole mettabolites
raw_data <- readxl::read_excel("Metabolites_NMR_MS_ALL_2018.xlsx") %>% 
  mutate(Sex = na_if(Sex, "NA"), Spine_Age = na_if(Spine_Age, "NA")) %>% 
  mutate(Sex = as.factor(Sex), Site = as.factor(Site), Mortality = as.factor(Mortality), ) # did not work for me
names(raw_data)

raw_data2<- FM_MET_18NDRout

# as data.frame
DM_data.frame2 <- data.frame(raw_data2)
names(DM_data.frame2)
dim(DM_data.frame2)
# 53 172

# dataframe with just met data columns
DM_data_pca2 <- DM_data.frame2[,c(10:172)] %>%
  print()

dim(DM_data_pca2) 
#53 163

# metabolite groups does not work as is code from last PCA the complete one
# load in object with met_groups
met_groups <- readxl::read_excel("Metabolite_groups_18.xlsx") %>% 
  print()

summary(DM_data_pca2)
# Serine        Putrescine    Trans_Hydroxyproline   Asparagine     
# Min.   : 23.5   Min.   : 3.62   Min.   :  0.242      Min.   :  0.000  
# 1st Qu.: 71.4   1st Qu.:14.60   1st Qu.:  1.240      1st Qu.:  0.000  
# Median : 93.2   Median :19.20   Median :  3.780      Median :  0.935  
# Mean   :105.2   Mean   :20.95   Mean   : 16.890      Mean   : 13.884  
# 3rd Qu.:133.0   3rd Qu.:23.70   3rd Qu.: 26.200      3rd Qu.: 23.400  
# Max.   :268.0   Max.   :65.20   Max.   :118.000      Max.   :130.000  
# Glutamine     Alpha_Aminoadipic_Acid Methionine_Sulfoxide
# Min.   : 0.00   Min.   : 0.542         Min.   : 1.030      
# 1st Qu.: 8.30   1st Qu.: 1.500         1st Qu.: 6.050      
# Median :15.70   Median : 2.900         Median : 7.990      
# Mean   :19.56   Mean   : 4.375         Mean   : 8.732      
# 3rd Qu.:25.70   3rd Qu.: 6.140         3rd Qu.:11.300      
# Max.   :84.00   Max.   :24.200         Max.   :20.900      
# Acetyl_Ornithine   Citrulline     Asymmetric_Dimethylarginine
# Min.   :0.305    Min.   :  2.01   Min.   :0.0000             
# 1st Qu.:1.440    1st Qu.:  3.43   1st Qu.:0.4390             
# Median :1.770    Median : 14.70   Median :0.5720             
# Mean   :2.130    Mean   : 46.54   Mean   :0.6021             
# 3rd Qu.:2.660    3rd Qu.: 42.30   3rd Qu.:0.7650             
# Max.   :5.960    Max.   :294.00   Max.   :1.9900             
# Total_Dimethylarginine   Tryptophan      Kynurenine       Ornithine    
# Min.   :0.263          Min.   : 5.98   Min.   :0.0000   Min.   : 0.00  
# 1st Qu.:0.850          1st Qu.:12.50   1st Qu.:0.0000   1st Qu.: 1.41  
# Median :1.090          Median :17.40   Median :0.0000   Median : 4.10  
# Mean   :1.147          Mean   :18.83   Mean   :0.2814   Mean   : 5.96  
# 3rd Qu.:1.400          3rd Qu.:21.80   3rd Qu.:0.3040   3rd Qu.: 7.04  
# Max.   :3.290          Max.   :45.20   Max.   :2.1900   Max.   :29.50  
# Lysine         Spermidine       Spermine       Sarcosine    
# Min.   :  0.82   Min.   : 3.23   Min.   :26.50   Min.   : 3.46  
# 1st Qu.:  4.56   1st Qu.: 9.49   1st Qu.:36.20   1st Qu.: 8.76  
# Median : 20.70   Median :13.70   Median :45.70   Median :15.10  
# Mean   : 27.86   Mean   :14.99   Mean   :46.48   Mean   :20.59  
# 3rd Qu.: 41.60   3rd Qu.:17.50   3rd Qu.:53.50   3rd Qu.:26.90  
# Max.   :110.00   Max.   :61.00   Max.   :82.70   Max.   :90.80  
# Methylhistidine Beta_Hydroxybutyric_Acid Alpha_Ketoglutaric_Acid
# Min.   : 2.38   Min.   : 0.736           Min.   :20.80          
# 1st Qu.: 6.94   1st Qu.: 4.570           1st Qu.:55.00          
# Median :13.30   Median : 6.680           Median :60.70          
# Mean   :19.95   Mean   : 8.094           Mean   :58.75          
# 3rd Qu.:31.70   3rd Qu.: 9.250           3rd Qu.:64.00          
# Max.   :85.80   Max.   :38.900           Max.   :71.20          
# Butyric_acid  Propionic_Acid   Fumaric_Acid   Isobutyric_Acid
# Min.   :1.28   Min.   :25.40   Min.   :21.20   Min.   :0.890  
# 1st Qu.:1.89   1st Qu.:35.20   1st Qu.:40.60   1st Qu.:1.720  
# Median :2.43   Median :42.70   Median :53.10   Median :2.520  
# Mean   :2.48   Mean   :44.95   Mean   :52.17   Mean   :2.867  
# 3rd Qu.:3.00   3rd Qu.:52.90   3rd Qu.:62.40   3rd Qu.:3.520  
# Max.   :4.58   Max.   :85.80   Max.   :99.40   Max.   :6.720  
# Hippuric_Acid     Methylmalonic_Acid   LYSOC14_0       LYSOC16_1     
# Min.   :0.01110   Min.   :0.02050    Min.   :1.852   Min.   : 4.323  
# 1st Qu.:0.03100   1st Qu.:0.04130    1st Qu.:3.029   1st Qu.: 7.990  
# Median :0.03670   Median :0.06410    Median :3.522   Median :10.371  
# Mean   :0.03744   Mean   :0.07654    Mean   :3.634   Mean   :10.769  
# 3rd Qu.:0.04390   3rd Qu.:0.08200    3rd Qu.:4.293   3rd Qu.:12.023  
# Max.   :0.07400   Max.   :0.73400    Max.   :6.116   Max.   :26.612  
# LYSOC16_0       LYSOC17_0       LYSOC18_2        LYSOC18_1     
# Min.   :18.97   Min.   :1.195   Min.   :0.3734   Min.   : 3.478  
# 1st Qu.:31.57   1st Qu.:2.095   1st Qu.:0.9065   1st Qu.: 7.836  
# Median :41.28   Median :2.579   Median :1.1541   Median : 9.397  
# Mean   :42.08   Mean   :2.760   Mean   :1.1893   Mean   : 9.706  
# 3rd Qu.:47.64   3rd Qu.:3.207   3rd Qu.:1.4132   3rd Qu.:11.978  
# Max.   :95.67   Max.   :6.248   Max.   :2.0876   Max.   :16.755  
# LYSOC18_0       LYSOC20_4        LYSOC20_3        LYSOC24_0    
# Min.   :1.389   Min.   : 0.985   Min.   : 6.953   Min.   :1.343  
# 1st Qu.:2.938   1st Qu.: 3.324   1st Qu.: 7.902   1st Qu.:2.518  
# Median :4.075   Median : 5.063   Median : 8.961   Median :2.876  
# Mean   :3.974   Mean   : 5.392   Mean   : 9.017   Mean   :3.048  
# 3rd Qu.:4.877   3rd Qu.: 6.244   3rd Qu.: 9.944   3rd Qu.:3.469  
# Max.   :7.102   Max.   :14.805   Max.   :12.636   Max.   :5.966  
# LYSOC26_1        LYSOC26_0        LYSOC28_1        LYSOC28_0     
# Min.   :0.7214   Min.   :0.8358   Min.   :0.9302   Min.   :0.5939  
# 1st Qu.:1.1996   1st Qu.:1.3596   1st Qu.:1.4988   1st Qu.:1.0115  
# Median :1.4203   Median :1.5663   Median :1.6560   Median :1.1970  
# Mean   :1.4411   Mean   :1.6417   Mean   :1.6750   Mean   :1.1949  
# 3rd Qu.:1.6799   3rd Qu.:1.8344   3rd Qu.:1.8642   3rd Qu.:1.3156  
# Max.   :2.2097   Max.   :2.8986   Max.   :2.6944   Max.   :1.9853  
# X14_1SMOH        X16_1SM          X16_0SM         X16_1SMOH     
# Min.   :1.762   Min.   : 3.842   Min.   : 2.759   Min.   : 2.322  
# 1st Qu.:3.263   1st Qu.: 6.019   1st Qu.: 5.197   1st Qu.: 4.324  
# Median :4.125   Median : 7.864   Median : 6.448   Median : 6.586  
# Mean   :4.191   Mean   : 7.989   Mean   : 6.993   Mean   : 6.560  
# 3rd Qu.:4.826   3rd Qu.: 9.738   3rd Qu.: 8.357   3rd Qu.: 8.137  
# Max.   :7.163   Max.   :14.747   Max.   :14.174   Max.   :14.324  
# X18_1SM         PC32_2AA        X18_0SM          X20_2SM      
# Min.   :1.446   Min.   :15.19   Min.   : 6.595   Min.   : 3.103  
# 1st Qu.:2.372   1st Qu.:23.10   1st Qu.: 9.742   1st Qu.: 6.911  
# Median :2.775   Median :30.06   Median :12.050   Median : 8.950  
# Mean   :2.995   Mean   :33.84   Mean   :13.148   Mean   : 9.277  
# 3rd Qu.:3.532   3rd Qu.:41.26   3rd Qu.:14.801   3rd Qu.:11.015  
# Max.   :5.857   Max.   :70.59   Max.   :29.556   Max.   :21.746  
# PC36_0AE         PC36_6AA         PC36_0AA        X22_2SMOH     
# Min.   : 4.782   Min.   : 33.39   Min.   : 18.22   Min.   : 3.798  
# 1st Qu.: 9.873   1st Qu.: 61.92   1st Qu.: 46.25   1st Qu.: 6.086  
# Median :13.353   Median : 78.68   Median : 80.24   Median : 8.284  
# Mean   :13.680   Mean   : 81.98   Mean   : 87.39   Mean   : 8.653  
# 3rd Qu.:17.052   3rd Qu.: 97.76   3rd Qu.:126.37   3rd Qu.:11.124  
# Max.   :25.719   Max.   :187.88   Max.   :209.69   Max.   :15.168  
# X22_1SMOH        PC38_6AA         PC38_0AA         PC40_6AE     
# Min.   :2.071   Min.   : 297.1   Min.   : 18.31   Min.   : 28.17  
# 1st Qu.:3.357   1st Qu.: 577.5   1st Qu.: 40.12   1st Qu.: 56.90  
# Median :4.213   Median : 697.8   Median : 54.80   Median : 68.45  
# Mean   :4.256   Mean   : 703.9   Mean   : 55.28   Mean   : 70.83  
# 3rd Qu.:5.151   3rd Qu.: 830.4   3rd Qu.: 70.40   3rd Qu.: 82.99  
# Max.   :6.856   Max.   :1337.0   Max.   :104.45   Max.   :143.91  
# X24_1SMOH         PC40_6AA         PC40_2AA        PC401AA     
# Min.   : 4.490   Min.   : 27.53   Min.   :1.489   Min.   :1.169  
# 1st Qu.: 9.049   1st Qu.: 78.78   1st Qu.:2.349   1st Qu.:2.220  
# Median :10.216   Median :113.36   Median :2.812   Median :2.764  
# Mean   :11.210   Mean   :120.07   Mean   :2.763   Mean   :2.679  
# 3rd Qu.:14.130   3rd Qu.:157.85   3rd Qu.:3.158   3rd Qu.:3.085  
# Max.   :21.248   Max.   :245.69   Max.   :4.721   Max.   :4.745  
# C2              C3_1               C3              C4_1        
# Min.   : 3.505   Min.   :0.01720   Min.   :0.4332   Min.   :0.03720  
# 1st Qu.:11.328   1st Qu.:0.03840   1st Qu.:0.9670   1st Qu.:0.05030  
# Median :14.913   Median :0.04490   Median :1.3656   Median :0.06340  
# Mean   :16.130   Mean   :0.04879   Mean   :1.5755   Mean   :0.06676  
# 3rd Qu.:19.889   3rd Qu.:0.06200   3rd Qu.:1.8675   3rd Qu.:0.07800  
# Max.   :37.159   Max.   :0.10650   Max.   :4.6204   Max.   :0.12330  
# C4              C3OH              C5_1               C5        
# Min.   :0.0617   Min.   :0.01140   Min.   :0.00860   Min.   :0.0303  
# 1st Qu.:0.1354   1st Qu.:0.02000   1st Qu.:0.02020   1st Qu.:0.0860  
# Median :0.1799   Median :0.02450   Median :0.02220   Median :0.1473  
# Mean   :0.2078   Mean   :0.02433   Mean   :0.02344   Mean   :0.2311  
# 3rd Qu.:0.2523   3rd Qu.:0.02800   3rd Qu.:0.02650   3rd Qu.:0.3018  
# Max.   :0.6779   Max.   :0.04350   Max.   :0.03890   Max.   :1.3979  
# C4OH              C6_1              C6              C5OH       
# Min.   :0.02170   Min.   :0.1171   Min.   :0.5102   Min.   :0.0662  
# 1st Qu.:0.03540   1st Qu.:0.1596   1st Qu.:0.6592   1st Qu.:0.1492  
# Median :0.04450   Median :0.1815   Median :0.7274   Median :0.1719  
# Mean   :0.05738   Mean   :0.1890   Mean   :0.7780   Mean   :0.1797  
# 3rd Qu.:0.06080   3rd Qu.:0.2149   3rd Qu.:0.8291   3rd Qu.:0.1995  
# Max.   :0.45600   Max.   :0.2802   Max.   :1.5829   Max.   :0.3586  
# C5_1DC            C5DC              C8             C5MDC        
# Min.   :0.0671   Min.   :0.0735   Min.   :0.1043   Min.   :0.02280  
# 1st Qu.:0.0859   1st Qu.:0.1652   1st Qu.:0.1783   1st Qu.:0.02730  
# Median :0.0972   Median :0.1926   Median :0.2074   Median :0.02930  
# Mean   :0.1041   Mean   :0.2051   Mean   :0.2294   Mean   :0.02998  
# 3rd Qu.:0.1173   3rd Qu.:0.2406   3rd Qu.:0.2565   3rd Qu.:0.03340  
# Max.   :0.2765   Max.   :0.5097   Max.   :0.5370   Max.   :0.03910  
# C9               C7DC             C10_2             C10_1       
# Min.   :0.01490   Min.   :0.00610   Min.   :0.03560   Min.   :0.1439  
# 1st Qu.:0.01940   1st Qu.:0.01340   1st Qu.:0.06520   1st Qu.:0.1745  
# Median :0.02400   Median :0.01570   Median :0.07600   Median :0.1921  
# Mean   :0.02461   Mean   :0.03919   Mean   :0.08816   Mean   :0.1903  
# 3rd Qu.:0.02850   3rd Qu.:0.05220   3rd Qu.:0.10080   3rd Qu.:0.2031  
# Max.   :0.03990   Max.   :0.17810   Max.   :0.21070   Max.   :0.2503  
# C10             C12_1              C12              C14_2        
# Min.   :0.1260   Min.   :0.02490   Min.   :0.01780   Min.   :0.01000  
# 1st Qu.:0.1318   1st Qu.:0.04940   1st Qu.:0.02690   1st Qu.:0.01430  
# Median :0.1356   Median :0.05390   Median :0.03060   Median :0.01570  
# Mean   :0.1379   Mean   :0.05514   Mean   :0.03405   Mean   :0.01677  
# 3rd Qu.:0.1419   3rd Qu.:0.06200   3rd Qu.:0.03780   3rd Qu.:0.01810  
# Max.   :0.1707   Max.   :0.09420   Max.   :0.08220   Max.   :0.02670  
# C14_1              C14              C12DC            C14_2OH       
# Min.   :0.03590   Min.   :0.01090   Min.   :0.00860   Min.   :0.00730  
# 1st Qu.:0.06060   1st Qu.:0.02710   1st Qu.:0.01400   1st Qu.:0.01390  
# Median :0.08340   Median :0.03860   Median :0.01750   Median :0.02000  
# Mean   :0.09107   Mean   :0.05097   Mean   :0.01883   Mean   :0.02356  
# 3rd Qu.:0.10100   3rd Qu.:0.06180   3rd Qu.:0.02170   3rd Qu.:0.02920  
# Max.   :0.25080   Max.   :0.17930   Max.   :0.03990   Max.   :0.06600  
# C14_1OH            C16_2             C16_1             C16        
# Min.   :0.00960   Min.   :0.00820   Min.   :0.0241   Min.   :0.0201  
# 1st Qu.:0.01960   1st Qu.:0.01070   1st Qu.:0.0720   1st Qu.:0.0741  
# Median :0.02690   Median :0.01380   Median :0.1125   Median :0.0982  
# Mean   :0.03354   Mean   :0.01512   Mean   :0.1436   Mean   :0.1216  
# 3rd Qu.:0.04330   3rd Qu.:0.01760   3rd Qu.:0.1723   3rd Qu.:0.1688  
# Max.   :0.11250   Max.   :0.03420   Max.   :0.4753   Max.   :0.3432  
# C16_2OH           C16_1OH            C16OH             C18_2        
# Min.   :0.00840   Min.   :0.01230   Min.   :0.00460   Min.   :0.00790  
# 1st Qu.:0.01450   1st Qu.:0.02400   1st Qu.:0.00790   1st Qu.:0.01610  
# Median :0.02020   Median :0.03540   Median :0.01030   Median :0.02070  
# Mean   :0.02368   Mean   :0.03891   Mean   :0.01117   Mean   :0.02643  
# 3rd Qu.:0.02790   3rd Qu.:0.04700   3rd Qu.:0.01330   3rd Qu.:0.03010  
# Max.   :0.06790   Max.   :0.10360   Max.   :0.02490   Max.   :0.09290  
# C18_1             C18            C18_1OH        Arginine_Average
# Min.   :0.0251   Min.   :0.0092   Min.   :0.00770   Min.   : 15.01  
# 1st Qu.:0.0852   1st Qu.:0.0272   1st Qu.:0.01320   1st Qu.: 34.76  
# Median :0.1362   Median :0.0353   Median :0.01790   Median : 64.06  
# Mean   :0.1747   Mean   :0.0425   Mean   :0.01955   Mean   : 65.23  
# 3rd Qu.:0.2366   3rd Qu.:0.0539   3rd Qu.:0.02410   3rd Qu.: 86.91  
# Max.   :0.5893   Max.   :0.1107   Max.   :0.04540   Max.   :181.06  
# Betaine_Average    C0_Average    Choline_Average  Citrate_Average 
# Min.   : 76.99   Min.   :13.79   Min.   : 24.64   Min.   : 51.34  
# 1st Qu.:145.75   1st Qu.:22.33   1st Qu.: 36.25   1st Qu.:139.44  
# Median :192.00   Median :31.13   Median : 44.02   Median :180.69  
# Mean   :201.12   Mean   :32.69   Mean   : 51.53   Mean   :228.48  
# 3rd Qu.:260.00   3rd Qu.:38.52   3rd Qu.: 56.33   3rd Qu.:289.88  
# Max.   :501.12   Max.   :70.55   Max.   :203.31   Max.   :823.38  
# Creatine_Average  Creatinine_Average Glucose_Average Glutamate_Average
# Min.   :  46.09   Min.   : 0.5705    Min.   : 1934   Min.   : 26.80   
# 1st Qu.: 249.38   1st Qu.: 5.8825    1st Qu.: 3867   1st Qu.: 69.44   
# Median : 355.12   Median : 8.8575    Median : 5790   Median : 82.67   
# Mean   : 470.96   Mean   :10.1706    Mean   : 6042   Mean   : 92.32   
# 3rd Qu.: 667.00   3rd Qu.:13.8250    3rd Qu.: 7627   3rd Qu.:108.38   
# Max.   :1351.12   Max.   :24.0750    Max.   :13459   Max.   :247.38   
# Glycine_Average Histidine_Average Isoleucine_Average Lactate_Average
# Min.   :127.8   Min.   :10.27     Min.   : 25.06     Min.   :1345   
# 1st Qu.:330.7   1st Qu.:35.70     1st Qu.: 55.48     1st Qu.:4242   
# Median :446.6   Median :45.95     Median : 77.95     Median :4841   
# Mean   :460.1   Mean   :49.53     Mean   : 96.28     Mean   :4750   
# 3rd Qu.:543.9   3rd Qu.:66.05     3rd Qu.:115.38     3rd Qu.:5539   
# Max.   :988.7   Max.   :92.21     Max.   :403.06     Max.   :6482   
# Leucine         Methionine     Phenylalanine_Average    Proline      
# Min.   : 18.18   Min.   : 3.795   Min.   : 16.73        Min.   : 10.19  
# 1st Qu.: 71.50   1st Qu.:18.900   1st Qu.: 46.27        1st Qu.: 30.40  
# Median :108.44   Median :27.113   Median : 60.54        Median : 50.21  
# Mean   :128.66   Mean   :29.146   Mean   : 66.39        Mean   : 69.08  
# 3rd Qu.:159.75   3rd Qu.:31.550   3rd Qu.: 82.33        3rd Qu.:101.69  
# Max.   :554.75   Max.   :76.463   Max.   :172.88        Max.   :249.31  
# Pyruvate_Average Succinate_Average Threonine_Average Tyrosine_Average
# Min.   : 47.36   Min.   : 3.587    Min.   : 23.27    Min.   : 10.13  
# 1st Qu.:122.81   1st Qu.: 7.710    1st Qu.: 69.16    1st Qu.: 35.76  
# Median :223.88   Median :10.188    Median :109.25    Median : 47.54  
# Mean   :196.77   Mean   :11.252    Mean   :117.95    Mean   : 58.75  
# 3rd Qu.:245.56   3rd Qu.:13.137    3rd Qu.:158.62    3rd Qu.: 60.36  
# Max.   :307.19   Max.   :24.900    Max.   :270.56    Max.   :381.19  
# Valine_Average   X1_Methylhistidine X2_Hydroxybutyrate
# Min.   : 49.99   Min.   : 1.000     Min.   : 2.25     
# 1st Qu.: 91.47   1st Qu.: 6.250     1st Qu.: 4.63     
# Median :137.88   Median : 7.000     Median : 6.88     
# Mean   :172.51   Mean   : 6.894     Mean   :10.88     
# 3rd Qu.:207.31   3rd Qu.: 7.500     3rd Qu.:13.38     
# Max.   :627.25   Max.   :14.500     Max.   :52.13     
# X2_Hydroxyisovaleric_Acid X3_Hydroxybutyrate      ADP        
# Min.   : 0.130            Min.   : 1.500     Min.   :  1.00  
# 1st Qu.: 0.880            1st Qu.: 4.630     1st Qu.: 12.00  
# Median : 1.130            Median : 6.380     Median : 42.00  
# Mean   : 1.689            Mean   : 8.056     Mean   : 56.02  
# 3rd Qu.: 1.750            3rd Qu.: 8.880     3rd Qu.: 75.00  
# Max.   :11.500            Max.   :40.880     Max.   :225.00  
# AMP             ATP            Acetamide        Acetate      
# Min.   : 23.0   Min.   :  0.800   Min.   : 1.80   Min.   : 31.00  
# 1st Qu.:150.0   1st Qu.:  2.000   1st Qu.: 7.30   1st Qu.: 63.00  
# Median :184.0   Median :  2.100   Median :11.00   Median : 80.00  
# Mean   :193.1   Mean   :  9.721   Mean   :12.09   Mean   : 81.58  
# 3rd Qu.:243.0   3rd Qu.:  3.100   3rd Qu.:17.30   3rd Qu.: 91.00  
# Max.   :382.0   Max.   :204.600   Max.   :25.10   Max.   :162.00  
# Acetoacetate      Acetone        Adenosine        Alanine     
# Min.   : 1.50   Min.   :0.880   Min.   : 10.0   Min.   : 51.0  
# 1st Qu.: 3.13   1st Qu.:1.250   1st Qu.: 66.0   1st Qu.:180.0  
# Median : 5.63   Median :1.750   Median :132.0   Median :244.0  
# Mean   :11.54   Mean   :1.875   Mean   :148.7   Mean   :271.4  
# 3rd Qu.: 8.25   3rd Qu.:2.000   3rd Qu.:195.0   3rd Qu.:343.0  
# Max.   :86.25   Max.   :7.130   Max.   :451.0   Max.   :857.0  
# Aspartate      Creatine_Phosphate    Cytidine    Dimethyl_Sulfone
# Min.   : 0.000   Min.   : 1.000     Min.   : 3.1   Min.   : 3.80   
# 1st Qu.: 6.875   1st Qu.: 3.100     1st Qu.: 8.9   1st Qu.: 8.80   
# Median : 9.875   Median : 3.100     Median :11.9   Median :11.30   
# Mean   :13.070   Mean   : 5.209     Mean   :12.3   Mean   :12.45   
# 3rd Qu.:15.500   3rd Qu.: 4.500     3rd Qu.:15.4   3rd Qu.:15.90   
# Max.   :44.250   Max.   :28.400     Max.   :23.5   Max.   :26.40   
# Ethanol        Ethanolamine      Formate         Glycerol     
# Min.   : 1.000   Min.   : 24.0   Min.   :14.60   Min.   :  66.0  
# 1st Qu.: 3.130   1st Qu.: 60.0   1st Qu.:20.40   1st Qu.: 125.0  
# Median : 3.750   Median : 91.0   Median :23.60   Median : 143.0  
# Mean   : 4.584   Mean   :105.9   Mean   :26.23   Mean   : 174.1  
# 3rd Qu.: 4.250   3rd Qu.:138.0   3rd Qu.:31.50   3rd Qu.: 177.0  
# Max.   :24.750   Max.   :258.0   Max.   :53.90   Max.   :1318.0  
# Guanidoacetate   Guanosine      Hypoxanthine        IMP       
# Min.   : 69    Min.   : 12.0   Min.   :11.10   Min.   : 21.9  
# 1st Qu.:123    1st Qu.: 62.0   1st Qu.:18.40   1st Qu.: 60.0  
# Median :167    Median :106.0   Median :27.00   Median : 84.4  
# Mean   :166    Mean   :127.5   Mean   :30.72   Mean   :118.2  
# 3rd Qu.:202    3rd Qu.:184.0   3rd Qu.:41.10   3rd Qu.:152.4  
# Max.   :280    Max.   :505.0   Max.   :73.30   Max.   :420.1  
# Inosine       Isopropanol        Malate          Malonate    
# Min.   : 42.0   Min.   :0.500   Min.   :  3.10   Min.   :14.00  
# 1st Qu.: 99.0   1st Qu.:1.000   1st Qu.:  8.80   1st Qu.:30.80  
# Median :140.0   Median :1.130   Median :  9.40   Median :39.50  
# Mean   :193.2   Mean   :1.161   Mean   : 16.87   Mean   :42.64  
# 3rd Qu.:234.0   3rd Qu.:1.250   3rd Qu.: 16.50   3rd Qu.:54.50  
# Max.   :717.0   Max.   :1.750   Max.   :102.90   Max.   :92.40  
# Mannose         Methanol     N_N_Dimethylglycine  Nicotinurate  
# Min.   : 3.30   Min.   :104.0   Min.   : 7.1        Min.   : 5.30  
# 1st Qu.:16.90   1st Qu.:126.0   1st Qu.:14.6        1st Qu.:17.00  
# Median :25.90   Median :142.0   Median :23.0        Median :27.90  
# Mean   :26.61   Mean   :153.7   Mean   :23.0        Mean   :27.53  
# 3rd Qu.:37.30   3rd Qu.:199.0   3rd Qu.:29.1        3rd Qu.:33.90  
# Max.   :61.30   Max.   :225.0   Max.   :48.1        Max.   :52.10  
# O_Acetylcarnitine O_Phosphocholine Propylene_Glycol    Taurine     
# Min.   : 3.1      Min.   :15.50    Min.   : 0.000   Min.   : 2787  
# 1st Qu.:11.3      1st Qu.:37.50    1st Qu.: 4.300   1st Qu.: 5352  
# Median :14.4      Median :45.60    Median : 6.100   Median : 6265  
# Mean   :15.9      Mean   :45.12    Mean   : 8.302   Mean   : 6601  
# 3rd Qu.:20.3      3rd Qu.:50.90    3rd Qu.:10.300   3rd Qu.: 7980  
# Max.   :33.0      Max.   :72.30    Max.   :50.400   Max.   :11645  
# Uridine        Xanthine      sn_Glycero_3_Phosphocholine
# Min.   : 4.1   Min.   :0.0000   Min.   : 16.00             
# 1st Qu.:13.1   1st Qu.:0.0000   1st Qu.: 39.00             
# Median :17.3   Median :0.7500   Median : 56.00             
# Mean   :18.1   Mean   :0.7524   Mean   : 73.23             
# 3rd Qu.:22.3   3rd Qu.:1.3750   3rd Qu.: 91.00             
# Max.   :36.8   Max.   :2.8750   Max.   :355.00             
# Beta.Alanine   
# Min.   :  31.0  
# 1st Qu.: 132.0  
# Median : 210.0  
# Mean   : 241.9  
# 3rd Qu.: 279.0  
# Max.   :1065.0  

#### Correlations plot ####

DM_cor2 <- cor(DM_data_pca2)
print(DM_data_pca2)
summary(DM_cor2)
# Serine            Putrescine        Trans_Hydroxyproline
# Min.   :-0.427345   Min.   :-0.471191   Min.   :-0.54141    
# 1st Qu.:-0.006915   1st Qu.:-0.006191   1st Qu.:-0.21043    
# Median : 0.176327   Median : 0.198896   Median : 0.06828    
# Mean   : 0.190669   Mean   : 0.182159   Mean   : 0.07688    
# 3rd Qu.: 0.377409   3rd Qu.: 0.389179   3rd Qu.: 0.33505    
# Max.   : 1.000000   Max.   : 1.000000   Max.   : 1.00000    
# Asparagine         Glutamine        Alpha_Aminoadipic_Acid
# Min.   :-0.50113   Min.   :-0.30755   Min.   :-0.4320       
# 1st Qu.:-0.19884   1st Qu.: 0.03814   1st Qu.:-0.0504       
# Median : 0.05857   Median : 0.15847   Median : 0.1247       
# Mean   : 0.08135   Mean   : 0.19454   Mean   : 0.1196       
# 3rd Qu.: 0.34891   3rd Qu.: 0.34225   3rd Qu.: 0.2797       
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000       
# Methionine_Sulfoxide Acetyl_Ornithine      Citrulline      
# Min.   :-0.3272      Min.   :-0.353166   Min.   :-0.50073  
# 1st Qu.: 0.0328      1st Qu.: 0.009821   1st Qu.:-0.21756  
# Median : 0.1440      Median : 0.138613   Median :-0.09785  
# Mean   : 0.1889      Mean   : 0.149658   Mean   :-0.06338  
# 3rd Qu.: 0.3213      3rd Qu.: 0.283209   3rd Qu.: 0.06695  
# Max.   : 1.0000      Max.   : 1.000000   Max.   : 1.00000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine   Tryptophan      
# Min.   :-0.49120            Min.   :-0.520292      Min.   :-0.54192  
# 1st Qu.:-0.04228            1st Qu.:-0.009831      1st Qu.:-0.05842  
# Median : 0.18606            Median : 0.202052      Median : 0.15089  
# Mean   : 0.20040            Mean   : 0.212490      Mean   : 0.18222  
# 3rd Qu.: 0.42575            3rd Qu.: 0.440275      3rd Qu.: 0.42831  
# Max.   : 1.00000            Max.   : 1.000000      Max.   : 1.00000  
# Kynurenine          Ornithine           Lysine        
# Min.   :-0.331759   Min.   :-0.4713   Min.   :-0.53802  
# 1st Qu.:-0.008328   1st Qu.:-0.0847   1st Qu.:-0.10712  
# Median : 0.119102   Median : 0.1076   Median : 0.07009  
# Mean   : 0.149080   Mean   : 0.1198   Mean   : 0.07456  
# 3rd Qu.: 0.297580   3rd Qu.: 0.2898   3rd Qu.: 0.24254  
# Max.   : 1.000000   Max.   : 1.0000   Max.   : 1.00000  
# Spermidine          Spermine         Sarcosine       
# Min.   :-0.53370   Min.   :-0.5432   Min.   :-0.30340  
# 1st Qu.:-0.09854   1st Qu.:-0.0699   1st Qu.: 0.00353  
# Median : 0.12775   Median : 0.1252   Median : 0.14122  
# Mean   : 0.13223   Mean   : 0.1245   Mean   : 0.17779  
# 3rd Qu.: 0.35268   3rd Qu.: 0.3233   3rd Qu.: 0.35007  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# Methylhistidine     Beta_Hydroxybutyric_Acid Alpha_Ketoglutaric_Acid
# Min.   :-0.547749   Min.   :-0.22975         Min.   :-0.2738        
# 1st Qu.:-0.144161   1st Qu.: 0.06705         1st Qu.: 0.1137        
# Median : 0.006971   Median : 0.18687         Median : 0.2393        
# Mean   :-0.002739   Mean   : 0.20357         Mean   : 0.2458        
# 3rd Qu.: 0.137285   3rd Qu.: 0.32071         3rd Qu.: 0.3756        
# Max.   : 1.000000   Max.   : 1.00000         Max.   : 1.0000        
# Butyric_acid      Propionic_Acid      Fumaric_Acid     
# Min.   :-0.48077   Min.   :-0.40098   Min.   :-0.38157  
# 1st Qu.:-0.13275   1st Qu.:-0.07779   1st Qu.: 0.07636  
# Median : 0.11761   Median : 0.02869   Median : 0.21892  
# Mean   : 0.09173   Mean   : 0.06217   Mean   : 0.21731  
# 3rd Qu.: 0.27595   3rd Qu.: 0.15692   3rd Qu.: 0.34098  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Isobutyric_Acid    Hippuric_Acid      Methylmalonic_Acid
# Min.   :-0.56184   Min.   :-0.39070   Min.   :-0.28720  
# 1st Qu.:-0.23372   1st Qu.:-0.11108   1st Qu.:-0.03047  
# Median : 0.02849   Median : 0.06090   Median : 0.04737  
# Mean   : 0.05613   Mean   : 0.06699   Mean   : 0.07028  
# 3rd Qu.: 0.32799   3rd Qu.: 0.24376   3rd Qu.: 0.16582  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC14_0          LYSOC16_1          LYSOC16_0       
# Min.   :-0.29907   Min.   :-0.32048   Min.   :-0.36971  
# 1st Qu.:-0.02588   1st Qu.:-0.04629   1st Qu.:-0.11573  
# Median : 0.11120   Median : 0.06141   Median : 0.03058  
# Mean   : 0.15368   Mean   : 0.11147   Mean   : 0.06689  
# 3rd Qu.: 0.28426   3rd Qu.: 0.24867   3rd Qu.: 0.23532  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC17_0          LYSOC18_2         LYSOC18_1         LYSOC18_0      
# Min.   :-0.38253   Min.   :-0.4224   Min.   :-0.3857   Min.   :-0.3205  
# 1st Qu.:-0.09604   1st Qu.: 0.0700   1st Qu.: 0.1147   1st Qu.: 0.1029  
# Median : 0.05008   Median : 0.2137   Median : 0.2520   Median : 0.2762  
# Mean   : 0.11937   Mean   : 0.2241   Mean   : 0.2473   Mean   : 0.2576  
# 3rd Qu.: 0.31775   3rd Qu.: 0.3180   3rd Qu.: 0.3619   3rd Qu.: 0.3981  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.0000  
# LYSOC20_4          LYSOC20_3           LYSOC24_0       
# Min.   :-0.35462   Min.   :-0.264122   Min.   :-0.34340  
# 1st Qu.:-0.08121   1st Qu.:-0.002536   1st Qu.:-0.06042  
# Median : 0.15807   Median : 0.155280   Median : 0.06823  
# Mean   : 0.17811   Mean   : 0.151826   Mean   : 0.13325  
# 3rd Qu.: 0.43256   3rd Qu.: 0.279178   3rd Qu.: 0.29573  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000  
# LYSOC26_1           LYSOC26_0          LYSOC28_1       
# Min.   :-0.353388   Min.   :-0.41049   Min.   :-0.25112  
# 1st Qu.: 0.003027   1st Qu.:-0.07910   1st Qu.: 0.04425  
# Median : 0.189886   Median : 0.08019   Median : 0.19480  
# Mean   : 0.196424   Mean   : 0.11450   Mean   : 0.20968  
# 3rd Qu.: 0.350701   3rd Qu.: 0.30276   3rd Qu.: 0.31482  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC28_0          X14_1SMOH           X16_1SM       
# Min.   :-0.35311   Min.   :-0.32102   Min.   :-0.5288  
# 1st Qu.:-0.02057   1st Qu.:-0.06495   1st Qu.:-0.1060  
# Median : 0.10662   Median : 0.11902   Median : 0.1382  
# Mean   : 0.14727   Mean   : 0.14604   Mean   : 0.1571  
# 3rd Qu.: 0.29435   3rd Qu.: 0.30786   3rd Qu.: 0.4053  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# X16_0SM           X16_1SMOH           X18_1SM       
# Min.   :-0.55132   Min.   :-0.54775   Min.   :-0.3099  
# 1st Qu.:-0.08681   1st Qu.:-0.19693   1st Qu.: 0.1013  
# Median : 0.14826   Median : 0.02536   Median : 0.2267  
# Mean   : 0.15832   Mean   : 0.06706   Mean   : 0.2535  
# 3rd Qu.: 0.40014   3rd Qu.: 0.30129   3rd Qu.: 0.4202  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# PC32_2AA           X18_0SM             X20_2SM        
# Min.   :-0.41063   Min.   :-0.390731   Min.   :-0.40444  
# 1st Qu.:-0.02721   1st Qu.: 0.007305   1st Qu.: 0.05269  
# Median : 0.18568   Median : 0.195269   Median : 0.23278  
# Mean   : 0.17985   Mean   : 0.209125   Mean   : 0.24051  
# 3rd Qu.: 0.36381   3rd Qu.: 0.406854   3rd Qu.: 0.41043  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000  
# PC36_0AE          PC36_6AA           PC36_0AA          X22_2SMOH      
# Min.   :-0.5434   Min.   :-0.34240   Min.   :-0.58839   Min.   :-0.5073  
# 1st Qu.:-0.1058   1st Qu.: 0.03023   1st Qu.:-0.17305   1st Qu.:-0.1168  
# Median : 0.1433   Median : 0.15144   Median : 0.09999   Median : 0.1141  
# Mean   : 0.1518   Mean   : 0.20782   Mean   : 0.10102   Mean   : 0.1452  
# 3rd Qu.: 0.3903   3rd Qu.: 0.35224   3rd Qu.: 0.36832   3rd Qu.: 0.3924  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# X22_1SMOH            PC38_6AA          PC38_0AA       
# Min.   :-0.322033   Min.   :-0.2757   Min.   :-0.41278  
# 1st Qu.:-0.004332   1st Qu.: 0.0488   1st Qu.:-0.08797  
# Median : 0.191621   Median : 0.1761   Median : 0.15743  
# Mean   : 0.223408   Mean   : 0.2170   Mean   : 0.18553  
# 3rd Qu.: 0.426157   3rd Qu.: 0.3566   3rd Qu.: 0.42568  
# Max.   : 1.000000   Max.   : 1.0000   Max.   : 1.00000  
# PC40_6AE          X24_1SMOH           PC40_6AA       
# Min.   :-0.26818   Min.   :-0.34711   Min.   :-0.44896  
# 1st Qu.: 0.04655   1st Qu.: 0.02437   1st Qu.: 0.08098  
# Median : 0.15938   Median : 0.16167   Median : 0.23873  
# Mean   : 0.21071   Mean   : 0.20031   Mean   : 0.24387  
# 3rd Qu.: 0.35507   3rd Qu.: 0.38921   3rd Qu.: 0.39515  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# PC40_2AA           PC401AA               C2          
# Min.   :-0.29866   Min.   :-0.43267   Min.   :-0.26067  
# 1st Qu.: 0.02819   1st Qu.: 0.07754   1st Qu.:-0.02255  
# Median : 0.12844   Median : 0.20913   Median : 0.12689  
# Mean   : 0.20085   Mean   : 0.23621   Mean   : 0.13551  
# 3rd Qu.: 0.31645   3rd Qu.: 0.36227   3rd Qu.: 0.29732  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C3_1                C3                C4_1        
# Min.   :-0.44531   Min.   :-0.57886   Min.   :-0.4073  
# 1st Qu.:-0.03308   1st Qu.:-0.09599   1st Qu.:-0.0774  
# Median : 0.10438   Median : 0.16863   Median : 0.1001  
# Mean   : 0.12093   Mean   : 0.16364   Mean   : 0.1247  
# 3rd Qu.: 0.26795   3rd Qu.: 0.39904   3rd Qu.: 0.3420  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# C4                C3OH               C5_1         
# Min.   :-0.54366   Min.   :-0.29486   Min.   :-0.47800  
# 1st Qu.:-0.04976   1st Qu.: 0.01386   1st Qu.:-0.08834  
# Median : 0.15207   Median : 0.10522   Median : 0.15413  
# Mean   : 0.18814   Mean   : 0.13471   Mean   : 0.12540  
# 3rd Qu.: 0.40914   3rd Qu.: 0.24020   3rd Qu.: 0.32224  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C5               C4OH               C6_1         
# Min.   :-0.5507   Min.   :-0.35788   Min.   :-0.51830  
# 1st Qu.:-0.1121   1st Qu.:-0.01839   1st Qu.:-0.12209  
# Median : 0.1015   Median : 0.07258   Median : 0.03895  
# Mean   : 0.1299   Mean   : 0.07358   Mean   : 0.02575  
# 3rd Qu.: 0.3717   3rd Qu.: 0.14146   3rd Qu.: 0.16797  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# C6                 C5OH              C5_1DC         
# Min.   :-0.268314   Min.   :-0.23316   Min.   :-0.383318  
# 1st Qu.: 0.004645   1st Qu.: 0.09246   1st Qu.: 0.005827  
# Median : 0.109406   Median : 0.19519   Median : 0.148490  
# Mean   : 0.159331   Mean   : 0.20119   Mean   : 0.155129  
# 3rd Qu.: 0.261937   3rd Qu.: 0.30142   3rd Qu.: 0.277764  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.000000  
# C5DC                C8                C5MDC         
# Min.   :-0.34816   Min.   :-0.314840   Min.   :-0.23622  
# 1st Qu.: 0.02926   1st Qu.: 0.006066   1st Qu.:-0.03537  
# Median : 0.13420   Median : 0.192218   Median : 0.11119  
# Mean   : 0.13675   Mean   : 0.199214   Mean   : 0.12238  
# 3rd Qu.: 0.22557   3rd Qu.: 0.365690   3rd Qu.: 0.25876  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000  
# C9                 C7DC              C10_2         
# Min.   :-0.340932   Min.   :-0.30912   Min.   :-0.43104  
# 1st Qu.: 0.002329   1st Qu.:-0.05273   1st Qu.:-0.10799  
# Median : 0.203603   Median : 0.02438   Median : 0.08272  
# Mean   : 0.216074   Mean   : 0.03705   Mean   : 0.10987  
# 3rd Qu.: 0.386955   3rd Qu.: 0.13022   3rd Qu.: 0.25529  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# C10_1               C10               C12_1         
# Min.   :-0.23800   Min.   :-0.34946   Min.   :-0.17470  
# 1st Qu.: 0.05396   1st Qu.:-0.05958   1st Qu.: 0.05869  
# Median : 0.14585   Median : 0.15726   Median : 0.17603  
# Mean   : 0.13034   Mean   : 0.17692   Mean   : 0.20149  
# 3rd Qu.: 0.22970   3rd Qu.: 0.34023   3rd Qu.: 0.31430  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C12               C14_2              C14_1         
# Min.   :-0.35744   Min.   :-0.37625   Min.   :-0.34065  
# 1st Qu.:-0.08137   1st Qu.:-0.01558   1st Qu.: 0.02676  
# Median : 0.10477   Median : 0.13842   Median : 0.22499  
# Mean   : 0.15981   Mean   : 0.18330   Mean   : 0.25589  
# 3rd Qu.: 0.34407   3rd Qu.: 0.37542   3rd Qu.: 0.42549  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C14               C12DC             C14_2OH           C14_1OH       
# Min.   :-0.42223   Min.   :-0.50064   Min.   :-0.4941   Min.   :-0.4567  
# 1st Qu.:-0.06956   1st Qu.:-0.09502   1st Qu.:-0.1054   1st Qu.:-0.1092  
# Median : 0.15296   Median : 0.16604   Median : 0.1666   Median : 0.1384  
# Mean   : 0.18552   Mean   : 0.16601   Mean   : 0.1744   Mean   : 0.1712  
# 3rd Qu.: 0.37657   3rd Qu.: 0.38594   3rd Qu.: 0.4054   3rd Qu.: 0.3486  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.0000  
# C16_2             C16_1              C16              C16_2OH       
# Min.   :-0.4468   Min.   :-0.4757   Min.   :-0.53108   Min.   :-0.4608  
# 1st Qu.:-0.1282   1st Qu.:-0.1189   1st Qu.:-0.15367   1st Qu.:-0.1310  
# Median : 0.1443   Median : 0.1307   Median : 0.09288   Median : 0.1161  
# Mean   : 0.1624   Mean   : 0.1627   Mean   : 0.12618   Mean   : 0.1532  
# 3rd Qu.: 0.3872   3rd Qu.: 0.3885   3rd Qu.: 0.33452   3rd Qu.: 0.3654  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# C16_1OH            C16OH               C18_2         
# Min.   :-0.4791   Min.   :-0.334291   Min.   :-0.43400  
# 1st Qu.:-0.1473   1st Qu.: 0.005014   1st Qu.:-0.08193  
# Median : 0.1347   Median : 0.169720   Median : 0.15422  
# Mean   : 0.1521   Mean   : 0.193170   Mean   : 0.19634  
# 3rd Qu.: 0.3516   3rd Qu.: 0.337505   3rd Qu.: 0.41424  
# Max.   : 1.0000   Max.   : 1.000000   Max.   : 1.00000  
# C18_1              C18              C18_1OH        
# Min.   :-0.5210   Min.   :-0.54141   Min.   :-0.46087  
# 1st Qu.:-0.1410   1st Qu.:-0.16790   1st Qu.:-0.09091  
# Median : 0.1304   Median : 0.09915   Median : 0.12904  
# Mean   : 0.1645   Mean   : 0.11884   Mean   : 0.15626  
# 3rd Qu.: 0.3878   3rd Qu.: 0.30637   3rd Qu.: 0.35218  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# Arginine_Average   Betaine_Average      C0_Average      
# Min.   :-0.48654   Min.   :-0.44161   Min.   :-0.25981  
# 1st Qu.:-0.01797   1st Qu.:-0.01653   1st Qu.: 0.02059  
# Median : 0.13109   Median : 0.09042   Median : 0.18664  
# Mean   : 0.15667   Mean   : 0.09472   Mean   : 0.19254  
# 3rd Qu.: 0.34045   3rd Qu.: 0.23755   3rd Qu.: 0.36319  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Choline_Average    Citrate_Average    Creatine_Average  
# Min.   :-0.40466   Min.   :-0.46193   Min.   :-0.37400  
# 1st Qu.:-0.05129   1st Qu.:-0.04362   1st Qu.:-0.05971  
# Median : 0.14625   Median : 0.13348   Median : 0.09258  
# Mean   : 0.15234   Mean   : 0.16874   Mean   : 0.11990  
# 3rd Qu.: 0.36886   3rd Qu.: 0.39484   3rd Qu.: 0.30589  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Creatinine_Average Glucose_Average    Glutamate_Average 
# Min.   :-0.4910    Min.   :-0.32288   Min.   :-0.24866  
# 1st Qu.:-0.1056    1st Qu.:-0.02224   1st Qu.: 0.05522  
# Median : 0.1149    Median : 0.12577   Median : 0.20167  
# Mean   : 0.1198    Mean   : 0.15658   Mean   : 0.21443  
# 3rd Qu.: 0.3323    3rd Qu.: 0.34683   3rd Qu.: 0.35124  
# Max.   : 1.0000    Max.   : 1.00000   Max.   : 1.00000  
# Glycine_Average    Histidine_Average  Isoleucine_Average
# Min.   :-0.32458   Min.   :-0.25966   Min.   :-0.5147   
# 1st Qu.: 0.01982   1st Qu.: 0.02273   1st Qu.:-0.0660   
# Median : 0.14856   Median : 0.22927   Median : 0.1108   
# Mean   : 0.18608   Mean   : 0.21194   Mean   : 0.1569   
# 3rd Qu.: 0.35022   3rd Qu.: 0.36030   3rd Qu.: 0.3766   
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000   
# Lactate_Average       Leucine           Methionine      
# Min.   :-0.41222   Min.   :-0.36055   Min.   :-0.39057  
# 1st Qu.: 0.03378   1st Qu.:-0.05186   1st Qu.: 0.05157  
# Median : 0.27966   Median : 0.10466   Median : 0.18018  
# Mean   : 0.23667   Mean   : 0.13150   Mean   : 0.22142  
# 3rd Qu.: 0.41588   3rd Qu.: 0.28639   3rd Qu.: 0.36677  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Phenylalanine_Average    Proline        Pyruvate_Average  
# Min.   :-0.371619     Min.   :-0.5101   Min.   :-0.40777  
# 1st Qu.: 0.005945     1st Qu.:-0.1242   1st Qu.:-0.03408  
# Median : 0.158427     Median : 0.1404   Median : 0.15212  
# Mean   : 0.205686     Mean   : 0.1453   Mean   : 0.12997  
# 3rd Qu.: 0.391671     3rd Qu.: 0.4479   3rd Qu.: 0.29753  
# Max.   : 1.000000     Max.   : 1.0000   Max.   : 1.00000  
# Succinate_Average  Threonine_Average  Tyrosine_Average  
# Min.   :-0.42918   Min.   :-0.50073   Min.   :-0.27403  
# 1st Qu.: 0.02183   1st Qu.: 0.09624   1st Qu.:-0.01554  
# Median : 0.18775   Median : 0.19024   Median : 0.11930  
# Mean   : 0.18336   Mean   : 0.20859   Mean   : 0.14538  
# 3rd Qu.: 0.35526   3rd Qu.: 0.34319   3rd Qu.: 0.26986  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Valine_Average     X1_Methylhistidine X2_Hydroxybutyrate
# Min.   :-0.53218   Min.   :-0.43198   Min.   :-0.36937  
# 1st Qu.:-0.09994   1st Qu.:-0.15548   1st Qu.:-0.16024  
# Median : 0.12027   Median :-0.03157   Median : 0.03393  
# Mean   : 0.16094   Mean   :-0.04337   Mean   : 0.05937  
# 3rd Qu.: 0.41425   3rd Qu.: 0.05164   3rd Qu.: 0.25219  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# X2_Hydroxyisovaleric_Acid X3_Hydroxybutyrate      ADP          
# Min.   :-0.257749         Min.   :-0.25808   Min.   :-0.56184  
# 1st Qu.:-0.092657         1st Qu.: 0.05523   1st Qu.:-0.28571  
# Median : 0.005605         Median : 0.17539   Median :-0.12027  
# Mean   : 0.006542         Mean   : 0.19478   Mean   :-0.10289  
# 3rd Qu.: 0.093500         3rd Qu.: 0.28887   3rd Qu.: 0.06949  
# Max.   : 1.000000         Max.   : 1.00000   Max.   : 1.00000  
# AMP                ATP             Acetamide       
# Min.   :-0.30912   Min.   :-0.32741   Min.   :-0.58839  
# 1st Qu.:-0.09507   1st Qu.:-0.18757   1st Qu.:-0.28531  
# Median : 0.00614   Median :-0.11128   Median :-0.07275  
# Mean   : 0.01279   Mean   :-0.08895   Mean   :-0.07453  
# 3rd Qu.: 0.11016   3rd Qu.:-0.02549   3rd Qu.: 0.14904  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Acetate          Acetoacetate         Acetone        
# Min.   :-0.18819   Min.   :-0.34240   Min.   :-0.27981  
# 1st Qu.:-0.01057   1st Qu.:-0.14956   1st Qu.:-0.01582  
# Median : 0.08392   Median :-0.06248   Median : 0.07809  
# Mean   : 0.11691   Mean   :-0.04443   Mean   : 0.13092  
# 3rd Qu.: 0.21218   3rd Qu.: 0.02522   3rd Qu.: 0.25814  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Adenosine           Alanine           Aspartate       
# Min.   :-0.41886   Min.   :-0.18120   Min.   :-0.43769  
# 1st Qu.:-0.07966   1st Qu.: 0.02043   1st Qu.:-0.13379  
# Median : 0.08704   Median : 0.11270   Median : 0.09107  
# Mean   : 0.09626   Mean   : 0.11988   Mean   : 0.10754  
# 3rd Qu.: 0.25554   3rd Qu.: 0.21022   3rd Qu.: 0.32025  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Creatine_Phosphate    Cytidine        Dimethyl_Sulfone  
# Min.   :-0.25902   Min.   :-0.38797   Min.   :-0.38359  
# 1st Qu.:-0.06514   1st Qu.:-0.11719   1st Qu.:-0.09115  
# Median : 0.05907   Median : 0.02626   Median : 0.09664  
# Mean   : 0.05287   Mean   : 0.03510   Mean   : 0.10715  
# 3rd Qu.: 0.15973   3rd Qu.: 0.15247   3rd Qu.: 0.28169  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Ethanol           Ethanolamine          Formate         
# Min.   :-0.388482   Min.   :-0.543238   Min.   :-0.321789  
# 1st Qu.:-0.168447   1st Qu.:-0.193439   1st Qu.:-0.132049  
# Median :-0.095347   Median :-0.002125   Median :-0.012035  
# Mean   :-0.064595   Mean   : 0.012499   Mean   :-0.006031  
# 3rd Qu.: 0.005542   3rd Qu.: 0.232765   3rd Qu.: 0.088470  
# Max.   : 1.000000   Max.   : 1.000000   Max.   : 1.000000  
# Glycerol         Guanidoacetate       Guanosine       
# Min.   :-0.177312   Min.   :-0.39258   Min.   :-0.40589  
# 1st Qu.:-0.038788   1st Qu.:-0.10587   1st Qu.:-0.15023  
# Median : 0.005588   Median : 0.05994   Median : 0.04171  
# Mean   : 0.027709   Mean   : 0.05043   Mean   : 0.04646  
# 3rd Qu.: 0.084159   3rd Qu.: 0.17640   3rd Qu.: 0.22469  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# Hypoxanthine           IMP              Inosine        
# Min.   :-0.34934   Min.   :-0.44236   Min.   :-0.41886  
# 1st Qu.:-0.12126   1st Qu.:-0.09520   1st Qu.:-0.18674  
# Median : 0.04735   Median : 0.06825   Median :-0.01795  
# Mean   : 0.07448   Mean   : 0.09971   Mean   : 0.02436  
# 3rd Qu.: 0.27818   3rd Qu.: 0.29471   3rd Qu.: 0.26576  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Isopropanol           Malate             Malonate        
# Min.   :-0.47605   Min.   :-0.192388   Min.   :-0.230575  
# 1st Qu.:-0.20069   1st Qu.:-0.008674   1st Qu.: 0.001764  
# Median :-0.04599   Median : 0.113895   Median : 0.071195  
# Mean   :-0.05142   Mean   : 0.138409   Mean   : 0.093913  
# 3rd Qu.: 0.06297   3rd Qu.: 0.236282   3rd Qu.: 0.161516  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.000000  
# Mannose            Methanol         N_N_Dimethylglycine
# Min.   :-0.49484   Min.   :-0.431384   Min.   :-0.39068   
# 1st Qu.:-0.02688   1st Qu.:-0.122063   1st Qu.:-0.08709   
# Median : 0.11361   Median : 0.017199   Median : 0.04877   
# Mean   : 0.10891   Mean   :-0.002519   Mean   : 0.07009   
# 3rd Qu.: 0.25422   3rd Qu.: 0.115819   3rd Qu.: 0.22074   
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000   
# Nicotinurate      O_Acetylcarnitine  O_Phosphocholine  
# Min.   :-0.27604   Min.   :-0.25992   Min.   :-0.26380  
# 1st Qu.: 0.01133   1st Qu.:-0.03171   1st Qu.:-0.04873  
# Median : 0.13862   Median : 0.11147   Median : 0.13983  
# Mean   : 0.15668   Mean   : 0.11716   Mean   : 0.13581  
# 3rd Qu.: 0.30068   3rd Qu.: 0.26854   3rd Qu.: 0.29864  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Propylene_Glycol      Taurine            Uridine        
# Min.   :-0.37654   Min.   :-0.23568   Min.   :-0.33752  
# 1st Qu.:-0.09834   1st Qu.:-0.04864   1st Qu.: 0.01171  
# Median : 0.03215   Median : 0.03683   Median : 0.17205  
# Mean   : 0.05134   Mean   : 0.07029   Mean   : 0.19062  
# 3rd Qu.: 0.21277   3rd Qu.: 0.14606   3rd Qu.: 0.34921  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Xanthine        sn_Glycero_3_Phosphocholine  Beta.Alanine     
# Min.   :-0.44161   Min.   :-0.36580            Min.   :-0.27547  
# 1st Qu.:-0.27922   1st Qu.:-0.16068            1st Qu.:-0.06905  
# Median :-0.13882   Median : 0.07206            Median : 0.02619  
# Mean   :-0.10638   Mean   : 0.03476            Mean   : 0.02245  
# 3rd Qu.: 0.07176   3rd Qu.: 0.18331            3rd Qu.: 0.09210  
# Max.   : 1.00000   Max.   : 1.00000            Max.   : 1.00000 


corrplot(DM_cor2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
plot_correlation182no43 <- corrplot(DM_cor2, type = "upper", order = "hclust", 
                               tl.col = "black", tl.srt = 60, 
                               tl.cex = 1, cl.cex = 1)
plot_correlation182no43
dev.list()
# Open the graphics device
png("plot_correlation182no43.png")
print("plot_correlation182no43")
dev.off()


#### PCA ####

# PCA with FactoMineR package

DM_pca_FM2 <- PCA(DM_data_pca2, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(DM_pca_FM2)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  29.8164044      18.29227266                    18.29227
# Dim.2  29.5020600      18.09942331                    36.39170
# Dim.3  15.1600431       9.30063994                    45.69234
# Dim.4  11.3646899       6.97220241                    52.66454
# Dim.5   7.9784321       4.89474364                    57.55928
# Dim.6   6.1259746       3.75826664                    61.31755
# Dim.7   4.3151605       2.64733770                    63.96489
# Dim.8   4.1536098       2.54822686                    66.51311
# Dim.9   3.8429673       2.35764867                    68.87076
# Dim.10  3.4816962       2.13600995                    71.00677
# Dim.11  3.2692331       2.00566445                    73.01244
# Dim.12  3.2042365       1.96578928                    74.97823
# Dim.13  2.8382056       1.74123045                    76.71946
# Dim.14  2.4386246       1.49608868                    78.21554
# Dim.15  2.2627212       1.38817250                    79.60372
# Dim.16  2.1278946       1.30545682                    80.90917
# Dim.17  2.0383915       1.25054696                    82.15972
# Dim.18  1.9668989       1.20668645                    83.36641
# Dim.19  1.8459992       1.13251484                    84.49892
# Dim.20  1.7340861       1.06385651                    85.56278
# Dim.21  1.5888706       0.97476727                    86.53755
# Dim.22  1.4608634       0.89623523                    87.43378
# Dim.23  1.3797107       0.84644829                    88.28023
# Dim.24  1.3338018       0.81828328                    89.09851
# Dim.25  1.3016362       0.79854984                    89.89706
# Dim.26  1.2128062       0.74405290                    90.64112
# Dim.27  1.1694980       0.71748341                    91.35860
# Dim.28  1.1011653       0.67556151                    92.03416
# Dim.29  1.0949539       0.67175089                    92.70591
# Dim.30  0.9679144       0.59381249                    93.29972
# Dim.31  0.9530056       0.58466602                    93.88439
# Dim.32  0.8816075       0.54086347                    94.42525
# Dim.33  0.8206931       0.50349272                    94.92875
# Dim.34  0.8078434       0.49560944                    95.42436
# Dim.35  0.7031409       0.43137477                    95.85573
# Dim.36  0.6750277       0.41412742                    96.26986
# Dim.37  0.6476531       0.39733316                    96.66719
# Dim.38  0.5953202       0.36522709                    97.03242
# Dim.39  0.5369166       0.32939669                    97.36181
# Dim.40  0.5201657       0.31912006                    97.68093
# Dim.41  0.4862691       0.29832458                    97.97926
# Dim.42  0.4529122       0.27786027                    98.25712
# Dim.43  0.4295696       0.26353964                    98.52066
# Dim.44  0.4107212       0.25197619                    98.77264
# Dim.45  0.3455951       0.21202152                    98.98466
# Dim.46  0.3410703       0.20924560                    99.19390
# Dim.47  0.2825768       0.17336003                    99.36726
# Dim.48  0.2692495       0.16518372                    99.53245
# Dim.49  0.2205493       0.13530633                    99.66775
# Dim.50  0.2145332       0.13161547                    99.79937
# Dim.51  0.1783822       0.10943692                    99.90880
# Dim.52  0.1486480       0.09119506                   100.00000              

fviz_screeplot(DM_pca_FM2, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
screeplot18_all2no43 <- fviz_screeplot(DM_pca_FM2, addlabels = TRUE, 
                                  ggtheme = theme_classic(),
                                  main = "",
                                  font.x = c(14, "bold"), font.y = c(14, "bold"),
                                  font.tickslab = 12,
                                  barfill = "#99d8c9", barcolor = "#66c2a4",
                                  font.submain = 16,
                                  ncp = Inf
)
screeplot18_all2no43

# Open the graphics device
png("screeplot18_all2no43")
dev.off()

fviz_screeplot(DM_pca_FM2, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 10,
               ncp = 16
)

screeplot18_partial2no43<- fviz_screeplot(DM_pca_FM2, addlabels = TRUE, 
                                     ggtheme = theme_classic(),
                                     main = "",
                                     font.x = c(14, "bold"), font.y = c(14, "bold"),
                                     font.tickslab = 12,
                                     barfill = "#99d8c9", barcolor = "#66c2a4",
                                     font.submain = 10,
                                     ncp = 16
)
screeplot18_partial2no43


# Open the graphics device
png("screeplot18_partial2no43")
dev.off()

# extract the results for variables
var2 <- get_pca_var(DM_pca_FM2)
var2
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
head(var2$coord)  # coordinates of variables
head(var2$contrib) # contributions of variables to the PCs
head(var2$cos2)  # Cos2: quality on the factore map
var2$coord
# Dim.1       Dim.2         Dim.3        Dim.4
# Serine                       0.556411004  0.33127822  0.5519716401 -0.374320521
# Putrescine                  -0.144991976  0.72049273 -0.0898051525 -0.061662493
# Trans_Hydroxyproline        -0.541296331  0.61885136  0.0870491472  0.143322782
# Asparagine                  -0.500331867  0.60522794  0.1092515433  0.119042335
# Glutamine                    0.218714607  0.55099656  0.4565055524 -0.269475355
# Alpha_Aminoadipic_Acid      -0.085393861  0.49186416  0.1290863341 -0.262518084
# Methionine_Sulfoxide         0.272471977  0.48501826  0.5339057703 -0.203027570
# Acetyl_Ornithine             0.001099744  0.52188630  0.0457176916 -0.025872362
# Citrulline                   0.024112348 -0.28620788 -0.4259257153  0.103231954
# Asymmetric_Dimethylarginine -0.135045640  0.78295528  0.2780184042 -0.033840694
# Total_Dimethylarginine      -0.024355679  0.77825839  0.4418556089 -0.073269173
# Tryptophan                  -0.198687151  0.77643135  0.3346937144 -0.007757709
# Kynurenine                   0.123104940  0.47297272  0.3390283408 -0.224818171
# Ornithine                    0.487520281  0.15772727  0.5148162608 -0.516362420
# Lysine                       0.446078764  0.03559760  0.4861389815 -0.517912751
# Spermidine                  -0.309391794  0.66985463 -0.0704923637  0.002176697
# Spermine                    -0.235212495  0.62009601 -0.0990818703 -0.052109631
# Sarcosine                    0.495504232  0.30055045  0.4489645918 -0.214393346
# Methylhistidine              0.282650240 -0.15273784  0.2731860387 -0.123138702
# Beta_Hydroxybutyric_Acid     0.345410793  0.36955877  0.1874744864  0.633309573
# Alpha_Ketoglutaric_Acid      0.303816770  0.63637176  0.2919029216 -0.112618254
# Butyric_acid                -0.245724075  0.46745323 -0.2214902604  0.217852284
# Propionic_Acid               0.178762154  0.10754013 -0.3836069650 -0.059688386
# Fumaric_Acid                 0.374140515  0.49283153  0.3435498724 -0.038623449
# Isobutyric_Acid             -0.569220604  0.55678531  0.0196478148  0.341640149
# Hippuric_Acid               -0.275601657  0.36153738 -0.1736921638  0.329191891
# Methylmalonic_Acid          -0.043595443  0.29527328 -0.0423611984 -0.022140522
# LYSOC14_0                    0.318211960  0.30774008 -0.5924816497 -0.211141495
# LYSOC16_1                    0.264370676  0.19370858 -0.5235227537 -0.156335739
# LYSOC16_0                    0.316525201 -0.01043467 -0.5240893116 -0.103115609
# LYSOC17_0                    0.414200104  0.11866416 -0.7027428810 -0.141973973
# LYSOC18_2                    0.429772573  0.51769343 -0.3121138123 -0.218788046
# LYSOC18_1                    0.463163703  0.55938470 -0.1975623354 -0.172168635
# LYSOC18_0                    0.498163497  0.57671255  0.0313687497 -0.244414763
# LYSOC20_4                   -0.197768679  0.74792811 -0.1928117234 -0.046700554
# LYSOC20_3                    0.296059505  0.35665848 -0.3338755811 -0.269303798
# LYSOC24_0                    0.556153077  0.05883561 -0.4641499439  0.124179762
# LYSOC26_1                    0.445228516  0.37560974 -0.5188920974  0.033951911
# LYSOC26_0                    0.487796947  0.05863877 -0.5409817423 -0.023619366
# LYSOC28_1                    0.303102617  0.49956409 -0.6064495915  0.036496601
# LYSOC28_0                    0.403944045  0.22738099 -0.5095813858  0.062301266
# X14_1SMOH                    0.006681616  0.47600218 -0.6918345801  0.030203969
# X16_1SM                     -0.314340193  0.72979871 -0.2708552559  0.284891044
# X16_0SM                     -0.324971166  0.73483768 -0.3220133770  0.234857010
# X16_1SMOH                   -0.441470679  0.49929843 -0.4325766440  0.163302339
# X18_1SM                      0.314456527  0.70214605 -0.1450885464 -0.342344564
# PC32_2AA                     0.648276173  0.24227247 -0.0935028568 -0.537180188
# X18_0SM                      0.682347886  0.31041166 -0.0143795424 -0.462036444
# X20_2SM                      0.546819407  0.52522478 -0.1842676448 -0.490396486
# PC36_0AE                    -0.290560105  0.70621284 -0.4129246425  0.046234578
# PC36_6AA                     0.263490214  0.57588081 -0.5480305083 -0.394570494
# PC36_0AA                    -0.455126677  0.63842528 -0.3579151328  0.107540740
# X22_2SMOH                   -0.277892157  0.66363461 -0.5675754129  0.085374724
# X22_1SMOH                    0.041894107  0.74837781 -0.5374614917 -0.057732744
# PC38_6AA                     0.335549977  0.54777019 -0.5383601700 -0.303640332
# PC38_0AA                    -0.125133798  0.72841072 -0.5627877729 -0.039085151
# PC40_6AE                     0.309462960  0.55041992 -0.5511347478 -0.301963785
# X24_1SMOH                    0.510126725  0.37497289 -0.5212797442 -0.147695994
# PC40_6AA                     0.497805879  0.55681038  0.0069446930 -0.302403739
# PC40_2AA                     0.283089099  0.51193509 -0.6697613140 -0.141401632
# PC401AA                      0.527515191  0.47561254 -0.4333932523 -0.195576658
# C2                           0.440385966  0.18535155  0.0261173804  0.103594458
# C3_1                        -0.127644221  0.46519812 -0.1127168962  0.273505976
# C3                          -0.266518409  0.76082387  0.0475737842  0.081195022
# C4_1                        -0.233079066  0.59960082 -0.3881371626 -0.039409296
# C4                          -0.200336851  0.74423970  0.0982635570  0.314646962
# C3OH                         0.079837467  0.40706869 -0.3204248345 -0.119400439
# C5_1                        -0.263520909  0.60233962 -0.0113764666  0.164166978
# C5                          -0.359492532  0.67279529  0.0837092708  0.158367142
# C4OH                        -0.023025694  0.23288409 -0.0538776846  0.080880555
# C6_1                        -0.140197760  0.13059387 -0.5287380178  0.346733145
# C6                           0.487310921  0.15458209 -0.1864071920  0.499568825
# C5OH                         0.260422228  0.49884054 -0.0375249183 -0.160132946
# C5_1DC                      -0.078197405  0.54181174 -0.0449745659  0.252604932
# C5DC                         0.064688671  0.41047987 -0.0353045453  0.183924558
# C8                           0.689842731  0.19042452  0.1167748661  0.197272132
# C5MDC                        0.495696693  0.05426124  0.1715214683  0.159884200
# C9                           0.764234970  0.17391756  0.0820783075  0.227375436
# C7DC                         0.179342563  0.01393773 -0.1223298637 -0.102136383
# C10_2                        0.697365207 -0.10857066  0.0227623257  0.196038502
# C10_1                        0.091195473  0.32296842 -0.1679154724  0.347138106
# C10                          0.758515325  0.04876658  0.1209442323  0.330031157
# C12_1                        0.452150441  0.31409584 -0.0593467870  0.427631980
# C12                          0.804449134 -0.05120477  0.0253871440  0.396000725
# C14_2                        0.719371308  0.09726115 -0.2339362237  0.214920190
# C14_1                        0.845399719  0.26080725 -0.0164412729  0.278298523
# C14                          0.893352404 -0.02046020  0.0328316625  0.356944862
# C12DC                        0.878241199 -0.03858655  0.0857566276  0.054517696
# C14_2OH                      0.929207890 -0.05871467  0.0713060691  0.185778306
# C14_1OH                      0.881563314 -0.07065857 -0.0399374048  0.373375110
# C16_2                        0.883757811 -0.05790067  0.1591499173  0.197144284
# C16_1                        0.902816775 -0.09084007  0.0798457840  0.264567634
# C16                          0.819976899 -0.16383447 -0.1880770109  0.336119531
# C16_2OH                      0.864080493 -0.10486139 -0.1629626022  0.327717297
# C16_1OH                      0.856829673 -0.09015711 -0.2681529431  0.224762728
# C16OH                        0.732537831  0.11543823  0.0657228121  0.322069669
# C18_2                        0.926181919  0.01885683  0.0671940201  0.247629088
# C18_1                        0.922542450 -0.09679035 -0.0518771005  0.273683527
# C18                          0.821678624 -0.16396102 -0.2414725136  0.177878840
# C18_1OH                      0.839036280 -0.08935645 -0.0778811821  0.393607938
# Arginine_Average             0.518213800  0.24932778  0.4646189915 -0.514875529
# Betaine_Average              0.354402424  0.07056341 -0.0508668615 -0.131551339
# C0_Average                   0.552913182  0.29280524  0.0672156436  0.119301263
# Choline_Average              0.599230044  0.15859178 -0.0307227713 -0.385451289
# Citrate_Average             -0.138147101  0.65900688  0.1683827908  0.278813237
# Creatine_Average            -0.195768114  0.50102939 -0.1261397511  0.298271843
# Creatinine_Average          -0.253751425  0.56583012 -0.0438274661  0.265548822
# Glucose_Average              0.134638872  0.50865318 -0.0238973344 -0.520702916
# Glutamate_Average            0.245124483  0.58281765  0.4238143315 -0.131520410
# Glycine_Average              0.323686458  0.46979241  0.5597711634 -0.193476251
# Histidine_Average            0.101717334  0.67380817  0.1848180286 -0.308749019
# Isoleucine_Average          -0.223338010  0.71548317  0.3300862229 -0.019879805
# Lactate_Average              0.032802934  0.78163342  0.1143995431  0.087012418
# Leucine                     -0.022161637  0.50590052  0.4930079580 -0.256168240
# Methionine                   0.244424303  0.60535734  0.5977661643 -0.228743227
# Phenylalanine_Average        0.177823582  0.64199426  0.5463019212 -0.300009028
# Proline                     -0.381049345  0.75987186  0.2284927467 -0.016473514
# Pyruvate_Average             0.543610252  0.11869647  0.3140880455 -0.347983505
# Succinate_Average            0.554783428  0.23511437  0.2785517080 -0.028531056
# Threonine_Average            0.235291178  0.55843311  0.3293505925 -0.159559693
# Tyrosine_Average             0.059086195  0.47534298  0.4678035610 -0.091951348
# Valine_Average              -0.291691903  0.76799699  0.3311619630 -0.002785009
# X1_Methylhistidine           0.093355017 -0.17649319  0.1510433690 -0.111308522
# X2_Hydroxybutyrate          -0.339202391  0.38082062  0.2449402202  0.361105098
# X2_Hydroxyisovaleric_Acid   -0.144084821  0.02714724  0.0396039197  0.250458354
# X3_Hydroxybutyrate           0.393762637  0.30615575  0.2018551844  0.653266196
# ADP                          0.234460878 -0.51817646 -0.0152066281 -0.207215430
# AMP                          0.175183819 -0.12473546  0.0932918436  0.133437227
# ATP                         -0.082198446 -0.30427796  0.0316826630  0.062749890
# Acetamide                    0.316079884 -0.54465729  0.0108671620 -0.096773325
# Acetate                      0.194509063  0.12750832  0.3095701112  0.417962445
# Acetoacetate                -0.093592197 -0.10406528  0.0450088496  0.340685435
# Acetone                      0.310196688  0.10867409  0.1391390542  0.667136149
# Adenosine                    0.432093408 -0.03387076  0.2991315172  0.212439705
# Alanine                      0.191151809  0.22233009  0.2810060385 -0.023358043
# Aspartate                   -0.308924085  0.47748649  0.2456231455  0.410274837
# Creatine_Phosphate           0.264647358 -0.04032803  0.0948609387 -0.033857597
# Cytidine                    -0.194794775  0.16070911  0.2949731342  0.200824322
# Dimethyl_Sulfone            -0.251695137  0.41664173  0.1488658546  0.484548784
# Ethanol                     -0.229312774 -0.14096697 -0.0474058121 -0.012525449
# Ethanolamine                 0.450159133 -0.38058267  0.1226645393  0.037595406
# Formate                     -0.209691716  0.04937571 -0.0687286999 -0.099265768
# Glycerol                     0.066266915  0.01306469  0.0201142516 -0.100370789
# Guanidoacetate               0.292765950 -0.13271594  0.2127260993  0.099375797
# Guanosine                    0.381256993 -0.15593623  0.2696150285  0.013794143
# Hypoxanthine                 0.484209733 -0.15355514 -0.0153563864  0.107928389
# IMP                         -0.234090179  0.43058753  0.2043669975  0.312427867
# Inosine                     -0.443662018  0.31370156  0.0466765787  0.279308982
# Isopropanol                  0.108081224 -0.30843735 -0.0242851477 -0.335784117
# Malate                       0.354553852  0.13860302  0.2794859837  0.497584276
# Malonate                     0.010079653  0.18260448  0.2030990879  0.372977817
# Mannose                      0.078217904  0.30886873  0.5048262527 -0.035034990
# Methanol                    -0.253138067  0.14503496 -0.0493342343 -0.157186634
# N_N_Dimethylglycine          0.464193571 -0.14415183  0.1718999425  0.176685520
# Nicotinurate                 0.086632115  0.39142952  0.4504828990  0.220232344
# O_Acetylcarnitine            0.354170099  0.16244405  0.0002316738  0.199356466
# O_Phosphocholine            -0.117221817  0.44757823  0.3285210701  0.210413395
# Propylene_Glycol            -0.185846062  0.22067125  0.1477574411  0.666062049
# Taurine                      0.040610271  0.05830368  0.2434041316  0.310607564
# Uridine                      0.040510685  0.54597838  0.4891299912  0.217135677
# Xanthine                    -0.485136717 -0.10800746  0.1050896402  0.213224634
# sn_Glycero_3_Phosphocholine  0.385355953 -0.17231654  0.0768988295 -0.257806837
# Beta.Alanine                 0.103199191 -0.04595353 -0.0356203842 -0.089278403
# Dim.5
# Serine                      -3.529456e-02
# Putrescine                  -2.385017e-01
# Trans_Hydroxyproline         1.326439e-01
# Asparagine                   1.323980e-01
# Glutamine                    2.941962e-02
# Alpha_Aminoadipic_Acid       6.514001e-02
# Methionine_Sulfoxide        -4.159108e-02
# Acetyl_Ornithine             2.343011e-01
# Citrulline                   6.671114e-01
# Asymmetric_Dimethylarginine -3.173886e-02
# Total_Dimethylarginine      -6.964762e-02
# Tryptophan                   2.505068e-01
# Kynurenine                   1.749992e-02
# Ornithine                   -1.479488e-01
# Lysine                      -2.104642e-01
# Spermidine                  -9.165500e-02
# Spermine                    -9.589155e-02
# Sarcosine                    1.153743e-01
# Methylhistidine             -1.534272e-01
# Beta_Hydroxybutyric_Acid     6.735647e-03
# Alpha_Ketoglutaric_Acid      2.110477e-01
# Butyric_acid                -3.066998e-01
# Propionic_Acid              -5.252086e-01
# Fumaric_Acid                 3.791472e-01
# Isobutyric_Acid              2.697174e-02
# Hippuric_Acid               -1.104884e-01
# Methylmalonic_Acid          -1.014934e-01
# LYSOC14_0                    3.385359e-01
# LYSOC16_1                    5.717302e-01
# LYSOC16_0                    5.213679e-01
# LYSOC17_0                    3.677419e-01
# LYSOC18_2                    1.366428e-01
# LYSOC18_1                    1.431692e-01
# LYSOC18_0                   -2.619043e-02
# LYSOC20_4                    1.530170e-01
# LYSOC20_3                   -2.320705e-01
# LYSOC24_0                   -2.048535e-01
# LYSOC26_1                   -2.862876e-01
# LYSOC26_0                   -1.637251e-01
# LYSOC28_1                   -1.860465e-01
# LYSOC28_0                   -3.320157e-01
# X14_1SMOH                   -1.624651e-02
# X16_1SM                      6.408117e-02
# X16_0SM                      6.127146e-02
# X16_1SMOH                    2.534326e-01
# X18_1SM                      1.473281e-02
# PC32_2AA                    -2.032513e-02
# X18_0SM                     -1.991494e-02
# X20_2SM                     -6.108211e-02
# PC36_0AE                     3.145857e-01
# PC36_6AA                     1.265628e-01
# PC36_0AA                     2.843460e-01
# X22_2SMOH                    1.094908e-01
# X22_1SMOH                    1.051926e-01
# PC38_6AA                     6.912096e-02
# PC38_0AA                     1.679212e-01
# PC40_6AE                    -6.026407e-03
# X24_1SMOH                    8.279924e-02
# PC40_6AA                    -1.995593e-01
# PC40_2AA                    -5.436816e-02
# PC401AA                     -1.123340e-01
# C2                          -5.042893e-01
# C3_1                        -1.229614e-01
# C3                           2.968345e-05
# C4_1                        -1.549624e-01
# C4                          -1.683224e-02
# C3OH                         7.130383e-02
# C5_1                        -3.123551e-01
# C5                          -1.332592e-02
# C4OH                        -7.816307e-02
# C6_1                         9.940360e-02
# C6                          -1.800676e-01
# C5OH                        -1.957481e-01
# C5_1DC                      -3.907124e-02
# C5DC                        -4.499681e-01
# C8                          -1.501367e-01
# C5MDC                       -3.177414e-02
# C9                          -2.353630e-02
# C7DC                         2.225795e-01
# C10_2                       -1.921829e-01
# C10_1                       -1.490539e-01
# C10                         -1.248971e-01
# C12_1                       -2.015227e-01
# C12                         -5.489169e-02
# C14_2                        1.902838e-01
# C14_1                       -4.716641e-02
# C14                         -3.501372e-02
# C12DC                       -1.116118e-03
# C14_2OH                     -1.765610e-02
# C14_1OH                      3.702819e-02
# C16_2                       -3.807664e-02
# C16_1                        2.592935e-02
# C16                         -1.318897e-02
# C16_2OH                      1.575490e-02
# C16_1OH                      7.450141e-02
# C16OH                       -1.061153e-01
# C18_2                       -4.714423e-02
# C18_1                       -1.171503e-02
# C18                          5.517637e-02
# C18_1OH                     -2.972659e-02
# Arginine_Average            -2.039402e-01
# Betaine_Average              5.509305e-01
# C0_Average                  -3.072262e-01
# Choline_Average              1.723609e-01
# Citrate_Average              1.748467e-01
# Creatine_Average            -1.214028e-01
# Creatinine_Average          -2.421873e-01
# Glucose_Average             -7.592552e-02
# Glutamate_Average            1.892447e-01
# Glycine_Average              3.220310e-03
# Histidine_Average           -1.225612e-01
# Isoleucine_Average           1.229399e-01
# Lactate_Average              9.377505e-03
# Leucine                      1.278876e-02
# Methionine                   2.381702e-02
# Phenylalanine_Average        9.176470e-02
# Proline                      1.410960e-01
# Pyruvate_Average            -2.432578e-01
# Succinate_Average            4.713492e-01
# Threonine_Average           -3.306160e-01
# Tyrosine_Average             2.059915e-01
# Valine_Average               1.746854e-01
# X1_Methylhistidine          -2.037913e-01
# X2_Hydroxybutyrate           2.712998e-01
# X2_Hydroxyisovaleric_Acid    1.986327e-01
# X3_Hydroxybutyrate          -4.515451e-02
# ADP                          3.424900e-01
# AMP                          3.870172e-04
# ATP                          2.203552e-01
# Acetamide                    1.690715e-01
# Acetate                      2.097006e-01
# Acetoacetate                -3.735849e-01
# Acetone                      1.936204e-01
# Adenosine                    3.789276e-01
# Alanine                      1.254397e-01
# Aspartate                    4.740577e-02
# Creatine_Phosphate           1.079441e-01
# Cytidine                    -2.796157e-02
# Dimethyl_Sulfone             1.744966e-01
# Ethanol                     -3.141073e-01
# Ethanolamine                 4.704277e-01
# Formate                     -3.254605e-01
# Glycerol                     5.938529e-02
# Guanidoacetate               2.228813e-01
# Guanosine                    5.771245e-01
# Hypoxanthine                 4.100356e-01
# IMP                          2.379970e-02
# Inosine                      1.293358e-02
# Isopropanol                  2.342803e-03
# Malate                       9.255765e-03
# Malonate                     1.480197e-01
# Mannose                     -3.658219e-01
# Methanol                    -5.196097e-01
# N_N_Dimethylglycine         -8.120879e-02
# Nicotinurate                 3.780405e-01
# O_Acetylcarnitine           -5.591922e-01
# O_Phosphocholine             2.223154e-01
# Propylene_Glycol            -1.290254e-01
# Taurine                      2.075958e-01
# Uridine                      1.220825e-01
# Xanthine                    -1.518620e-01
# sn_Glycero_3_Phosphocholine  2.580126e-01
# Beta.Alanine                -1.077883e-01



var2$contrib
# Dim.1        Dim.2        Dim.3        Dim.4
# Serine                      1.038332e+00 0.3719918571 2.009709e+00 1.232905e+00
# Putrescine                  7.050707e-02 1.7595712755 5.319883e-02 3.345681e-02
# Trans_Hydroxyproline        9.826863e-01 1.2981364836 4.998372e-02 1.807477e-01
# Asparagine                  8.395780e-01 1.2416111207 7.873262e-02 1.246939e-01
# Glutamine                   1.604354e-01 1.0290712051 1.374649e+00 6.389701e-01
# Alpha_Aminoadipic_Acid      2.445671e-02 0.8200456351 1.099158e-01 6.064023e-01
# Methionine_Sulfoxide        2.489937e-01 0.7973772349 1.880307e+00 3.627041e-01
# Acetyl_Ornithine            4.056276e-06 0.9232077808 1.378695e-02 5.889990e-03
# Citrulline                  1.949951e-03 0.2776584151 1.196650e+00 9.377147e-02
# Asymmetric_Dimethylarginine 6.116541e-02 2.0778852943 5.098550e-01 1.007676e-02
# Total_Dimethylarginine      1.989506e-03 2.0530299549 1.287835e+00 4.723729e-02
# Tryptophan                  1.323989e-01 2.0434018587 7.389153e-01 5.295529e-04
# Kynurenine                  5.082714e-02 0.7582629680 7.581787e-01 4.447390e-01
# Ornithine                   7.971317e-01 0.0843259476 1.748252e+00 2.346128e+00
# Lysine                      6.673718e-01 0.0042952562 1.558908e+00 2.360237e+00
# Spermidine                  3.210423e-01 1.5209284751 3.277810e-02 4.169061e-05
# Spermine                    1.855519e-01 1.3033634451 6.475718e-02 2.389342e-02
# Sarcosine                   8.234542e-01 0.3061839518 1.329608e+00 4.044502e-01
# Methylhistidine             2.679436e-01 0.0790753178 4.922850e-01 1.334233e-01
# Beta_Hydroxybutyric_Acid    4.001442e-01 0.4629293261 2.318376e-01 3.529186e+00
# Alpha_Ketoglutaric_Acid     3.095767e-01 1.3726804649 5.620519e-01 1.115989e-01
# Butyric_acid                2.025071e-01 0.7406686854 3.236002e-01 4.176059e-01
# Propionic_Acid              1.071756e-01 0.0392002407 9.706721e-01 3.134888e-02
# Fumaric_Acid                4.694769e-01 0.8232744449 7.785368e-01 1.312637e-02
# Isobutyric_Acid             1.086691e+00 1.0508075867 2.546409e-03 1.027023e+00
# Hippuric_Acid               2.547466e-01 0.4430513457 1.990032e-01 9.535438e-01
# Methylmalonic_Acid          6.374218e-03 0.2955261787 1.183685e-02 4.313384e-03
# LYSOC14_0                   3.396079e-01 0.3210079541 2.315524e+00 3.922741e-01
# LYSOC16_1                   2.344074e-01 0.1271877791 1.807885e+00 2.150597e-01
# LYSOC16_0                   3.360171e-01 0.0003690667 1.811800e+00 9.356022e-02
# LYSOC17_0                   5.753937e-01 0.0477294890 3.257560e+00 1.773617e-01
# LYSOC18_2                   6.194726e-01 0.9084331442 6.425775e-01 4.212012e-01
# LYSOC18_1                   7.194718e-01 1.0606420032 2.574589e-01 2.608258e-01
# LYSOC18_0                   8.323166e-01 1.1273699722 6.490737e-03 5.256507e-01
# LYSOC20_4                   1.311776e-01 1.8961267751 2.452260e-01 1.919051e-02
# LYSOC20_3                   2.939698e-01 0.4311742091 7.353073e-01 6.381567e-01
# LYSOC24_0                   1.037369e+00 0.0117335167 1.421072e+00 1.356888e-01
# LYSOC26_1                   6.648301e-01 0.4782129602 1.776044e+00 1.014310e-02
# LYSOC26_0                   7.980367e-01 0.0116551347 1.930478e+00 4.908840e-03
# LYSOC28_1                   3.081230e-01 0.8459215298 2.425990e+00 1.172053e-02
# LYSOC28_0                   5.472517e-01 0.1752491650 1.712879e+00 3.415357e-02
# X14_1SMOH                   1.497296e-04 0.7680076378 3.157215e+00 8.027317e-03
# X16_1SM                     3.313939e-01 1.8053185363 4.839206e-01 7.141674e-01
# X16_0SM                     3.541884e-01 1.8303345921 6.839863e-01 4.853438e-01
# X16_1SMOH                   6.536548e-01 0.8450220785 1.234314e+00 2.346536e-01
# X18_1SM                     3.316393e-01 1.6711005051 1.388564e-01 1.031263e+00
# PC32_2AA                    1.409499e+00 0.1989554247 5.766992e-02 2.539115e+00
# X18_0SM                     1.561552e+00 0.3266056619 1.363922e-03 1.878429e+00
# X20_2SM                     1.002842e+00 0.9350569629 2.239741e-01 2.116104e+00
# PC36_0AE                    2.831501e-01 1.6905144027 1.124712e+00 1.880945e-02
# PC36_6AA                    2.328486e-01 1.1241204923 1.981112e+00 1.369909e+00
# PC36_0AA                    6.947192e-01 1.3815538169 8.450058e-01 1.017627e-01
# X22_2SMOH                   2.589985e-01 1.4928140301 2.124940e+00 6.413588e-02
# X22_1SMOH                   5.886411e-03 1.8984076055 1.905436e+00 2.932829e-02
# PC38_6AA                    3.776236e-01 1.0170550096 1.911813e+00 8.112624e-01
# PC38_0AA                    5.251628e-02 1.7984580607 2.089243e+00 1.344207e-02
# PC40_6AE                    3.211900e-01 1.0269184053 2.003619e+00 8.023283e-01
# X24_1SMOH                   8.727722e-01 0.4765927032 1.792426e+00 1.919463e-01
# PC40_6AA                    8.311220e-01 1.0509022146 3.181308e-04 8.046680e-01
# PC40_2AA                    2.687763e-01 0.8883364053 2.958964e+00 1.759346e-01
# PC401AA                     9.332858e-01 0.7667508328 1.238979e+00 3.365708e-01
# C2                          6.504466e-01 0.1164501580 4.499443e-03 9.443119e-02
# C3_1                        5.464457e-02 0.7335396058 8.380648e-02 6.582275e-01
# C3                          2.382315e-01 1.9620763893 1.492915e-02 5.800978e-02
# C4_1                        1.822012e-01 1.2186306281 9.937337e-01 1.366595e-02
# C4                          1.346066e-01 1.8774713612 6.369195e-02 8.711431e-01
# C3OH                        2.137756e-02 0.5616723642 6.772545e-01 1.254453e-01
# C5_1                        2.329029e-01 1.2297887545 8.537178e-04 2.371450e-01
# C5                          4.334355e-01 1.5343115314 4.622178e-02 2.206849e-01
# C4OH                        1.778157e-03 0.1838346117 1.914774e-02 5.756131e-02
# C6_1                        6.592147e-02 0.0578087025 1.844084e+00 1.057872e+00
# C6                          7.964473e-01 0.0809964555 2.292054e-01 2.196004e+00
# C5OH                        2.274578e-01 0.8434728983 9.288361e-03 2.256336e-01
# C5_1DC                      2.050829e-02 0.9950490436 1.334239e-02 5.614694e-01
# C5DC                        1.403464e-02 0.5711252747 8.221685e-03 2.976609e-01
# C8                          1.596044e+00 0.1229117531 8.994941e-02 3.424316e-01
# C5MDC                       8.240940e-01 0.0099799205 1.940602e-01 2.249332e-01
# C9                          1.958838e+00 0.1025261233 4.443819e-02 4.549142e-01
# C7DC                        1.078727e-01 0.0006584638 9.871077e-02 9.179169e-02
# C10_2                       1.631043e+00 0.0399551325 3.417691e-03 3.381623e-01
# C10_1                       2.789275e-02 0.3535637918 1.859863e-01 1.060344e+00
# C10                         1.929627e+00 0.0080610610 9.648724e-02 9.584121e-01
# C12_1                       6.856629e-01 0.3344044428 2.323240e-02 1.609099e+00
# C12                         2.170411e+00 0.0088872739 4.251354e-03 1.379858e+00
# C14_2                       1.735605e+00 0.0320646476 3.609895e-01 4.064404e-01
# C14_1                       2.397005e+00 0.2305616074 1.783078e-03 6.814974e-01
# C14                         2.676642e+00 0.0014189510 7.110257e-03 1.121101e+00
# C12DC                       2.586857e+00 0.0050468414 4.851041e-02 2.615275e-02
# C14_2OH                     2.895813e+00 0.0116853277 3.353919e-02 3.036913e-01
# C14_1OH                     2.606464e+00 0.0169229979 1.052105e-02 1.226685e+00
# C16_2                       2.619457e+00 0.0113635705 1.670754e-01 3.419879e-01
# C16_1                       2.733657e+00 0.0279706516 4.205364e-02 6.159080e-01
# C16                         2.255007e+00 0.0909825740 2.333302e-01 9.940996e-01
# C16_2OH                     2.504108e+00 0.0372716697 1.751763e-01 9.450203e-01
# C16_1OH                     2.462259e+00 0.0275516501 4.743126e-01 4.445197e-01
# C16OH                       1.799720e+00 0.0451696768 2.849258e-02 9.127294e-01
# C18_2                       2.876983e+00 0.0012052720 2.978248e-02 5.395674e-01
# C18_1                       2.854417e+00 0.0317549744 1.775215e-02 6.590824e-01
# C18                         2.264377e+00 0.0911231825 3.846228e-01 2.784139e-01
# C18_1OH                     2.361056e+00 0.0270644680 4.000964e-02 1.363233e+00
# Arginine_Average            9.006637e-01 0.2107118692 1.423946e+00 2.332636e+00
# Betaine_Average             4.212482e-01 0.0168774482 1.706748e-02 1.522765e-01
# C0_Average                  1.025318e+00 0.2906065149 2.980165e-02 1.252370e-01
# Choline_Average             1.204292e+00 0.0852528663 6.226161e-03 1.307319e+00
# Citrate_Average             6.400712e-02 1.4720669443 1.870230e-01 6.840206e-01
# Creatine_Average            1.285371e-01 0.8508912659 1.049551e-01 7.828290e-01
# Creatinine_Average          2.159542e-01 1.0852249868 1.267046e-02 6.204848e-01
# Glucose_Average             6.079749e-02 0.8769830215 3.767025e-03 2.385736e+00
# Glutamate_Average           2.015200e-01 1.1513650822 1.184816e+00 1.522049e-01
# Glycine_Average             3.513936e-01 0.7480999932 2.066905e+00 3.293804e-01
# Histidine_Average           3.470042e-02 1.5389347611 2.253140e-01 8.387907e-01
# Isoleucine_Average          1.672900e-01 1.7351878581 7.187111e-01 3.477496e-03
# Lactate_Average             3.608861e-03 2.0708750509 8.632730e-02 6.662004e-02
# Leucine                     1.647208e-03 0.8675168454 1.603273e+00 5.774215e-01
# Methionine                  2.003704e-01 1.2421420945 2.357014e+00 4.604038e-01
# Phenylalanine_Average       1.060531e-01 1.3970435413 1.968634e+00 7.919742e-01
# Proline                     4.869756e-01 1.9571692589 3.443851e-01 2.387893e-03
# Pyruvate_Average            9.911058e-01 0.0477554869 6.507323e-01 1.065515e+00
# Succinate_Average           1.032266e+00 0.1873725626 5.118129e-01 7.162722e-03
# Threonine_Average           1.856761e-01 1.0570364700 7.155112e-01 2.240210e-01
# Tyrosine_Average            1.170892e-02 0.7658819512 1.443533e+00 7.439755e-02
# Valine_Average              2.853603e-01 1.9992480999 7.234033e-01 6.824889e-05
# X1_Methylhistidine          2.922941e-02 0.1055853305 1.504884e-01 1.090183e-01
# X2_Hydroxybutyrate          3.858891e-01 0.4915736097 3.957490e-01 1.147386e+00
# X2_Hydroxyisovaleric_Acid   6.962756e-02 0.0024980383 1.034608e-02 5.519674e-01
# X3_Hydroxybutyrate          5.200124e-01 0.3177111816 2.687691e-01 3.755111e+00
# ADP                         1.843680e-01 0.9101291291 1.525336e-03 3.778214e-01
# AMP                         1.029278e-01 0.0527384718 5.740992e-02 1.566738e-01
# ATP                         2.266063e-02 0.3138257995 6.621295e-03 3.464722e-02
# Acetamide                   3.350722e-01 1.0055282864 7.789899e-04 8.240503e-02
# Acetate                     1.268891e-01 0.0551092749 6.321463e-01 1.537152e+00
# Acetoacetate                2.937812e-02 0.0367078885 1.336274e-02 1.021291e+00
# Acetone                     3.227149e-01 0.0400312972 1.277020e-01 3.916259e+00
# Adenosine                   6.261812e-01 0.0038886375 5.902336e-01 3.971127e-01
# Alanine                     1.225467e-01 0.1675498938 5.208718e-01 4.800819e-03
# Aspartate                   3.200724e-01 0.7728048560 3.979588e-01 1.481127e+00
# Creatine_Phosphate          2.348983e-01 0.0055126660 5.935734e-02 1.008683e-02
# Cytidine                    1.272622e-01 0.0875444567 5.739374e-01 3.548747e-01
# Dimethyl_Sulfone            2.124684e-01 0.5884007029 1.461806e-01 2.065939e+00
# Ethanol                     1.763605e-01 0.0673569422 1.482391e-02 1.380476e-03
# Ethanolamine                6.796368e-01 0.4909595132 9.925163e-02 1.243689e-02
# Formate                     1.474712e-01 0.0082636967 3.115845e-02 8.670446e-02
# Glycerol                    1.472781e-02 0.0005785567 2.668746e-03 8.864558e-02
# Guanidoacetate              2.874656e-01 0.0597026814 2.984978e-01 8.689677e-02
# Guanosine                   4.875064e-01 0.0824217275 4.794991e-01 1.674295e-03
# Hypoxanthine                7.863425e-01 0.0799238423 1.555527e-03 1.024976e-01
# IMP                         1.837854e-01 0.6284497387 2.754997e-01 8.588987e-01
# Inosine                     6.601600e-01 0.3335654172 1.437135e-02 6.864552e-01
# Isopropanol                 3.917827e-02 0.3224642669 3.890282e-03 9.921166e-01
# Malate                      4.216083e-01 0.0651168026 5.152519e-01 2.178591e+00
# Malonate                    3.407500e-04 0.1130239627 2.720918e-01 1.224076e+00
# Mannose                     2.051904e-02 0.3233668899 1.681061e+00 1.080056e-02
# Methanol                    2.149115e-01 0.0713005795 1.605448e-02 2.174071e-01
# N_N_Dimethylglycine         7.226749e-01 0.0704349134 1.949176e-01 2.746909e-01
# Nicotinurate                2.517112e-02 0.5193436197 1.338617e+00 4.267805e-01
# O_Acetylcarnitine           4.206961e-01 0.0894448397 3.540408e-07 3.497060e-01
# O_Phosphocholine            4.608522e-02 0.6790246881 7.119115e-01 3.895733e-01
# Propylene_Glycol            1.158381e-01 0.1650589882 1.440119e-01 3.903658e+00
# Taurine                     5.531164e-03 0.0115223119 3.908008e-01 8.489194e-01
# Uridine                     5.504069e-03 1.0104121364 1.578150e+00 4.148631e-01
# Xanthine                    7.893562e-01 0.0395416857 7.284829e-02 4.000527e-01
# sn_Glycero_3_Phosphocholine 4.980453e-01 0.1006471757 3.900668e-02 5.848322e-01
# Beta.Alanine                3.571884e-02 0.0071578973 8.369447e-03 7.013507e-02
# Dim.5
# Serine                      1.561342e-02
# Putrescine                  7.129601e-01
# Trans_Hydroxyproline        2.205245e-01
# Asparagine                  2.197076e-01
# Glutamine                   1.084817e-02
# Alpha_Aminoadipic_Acid      5.318364e-02
# Methionine_Sulfoxide        2.168118e-02
# Acetyl_Ornithine            6.880678e-01
# Citrulline                  5.578009e+00
# Asymmetric_Dimethylarginine 1.262598e-02
# Total_Dimethylarginine      6.079881e-02
# Tryptophan                  7.865409e-01
# Kynurenine                  3.838441e-03
# Ornithine                   2.743503e-01
# Lysine                      5.551865e-01
# Spermidine                  1.052919e-01
# Spermine                    1.152506e-01
# Sarcosine                   1.668402e-01
# Methylhistidine             2.950441e-01
# Beta_Hydroxybutyric_Acid    5.686449e-04
# Alpha_Ketoglutaric_Acid     5.582691e-01
# Butyric_acid                1.178988e+00
# Propionic_Acid              3.457371e+00
# Fumaric_Acid                1.801765e+00
# Isobutyric_Acid             9.118018e-03
# Hippuric_Acid               1.530087e-01
# Methylmalonic_Acid          1.291095e-01
# LYSOC14_0                   1.436454e+00
# LYSOC16_1                   4.096988e+00
# LYSOC16_0                   3.406992e+00
# LYSOC17_0                   1.694996e+00
# LYSOC18_2                   2.340215e-01
# LYSOC18_1                   2.569104e-01
# LYSOC18_0                   8.597409e-03
# LYSOC20_4                   2.934688e-01
# LYSOC20_3                   6.750290e-01
# LYSOC24_0                   5.259800e-01
# LYSOC26_1                   1.027277e+00
# LYSOC26_0                   3.359795e-01
# LYSOC28_1                   4.338359e-01
# LYSOC28_0                   1.381655e+00
# X14_1SMOH                   3.308283e-03
# X16_1SM                     5.146872e-02
# X16_0SM                     4.705425e-02
# X16_1SMOH                   8.050216e-01
# X18_1SM                     2.720531e-03
# PC32_2AA                    5.177847e-03
# X18_0SM                     4.970963e-03
# X20_2SM                     4.676388e-02
# PC36_0AE                    1.240396e+00
# PC36_6AA                    2.007679e-01
# PC36_0AA                    1.013390e+00
# X22_2SMOH                   1.502579e-01
# X22_1SMOH                   1.386924e-01
# PC38_6AA                    5.988279e-02
# PC38_0AA                    3.534220e-01
# PC40_6AE                    4.551970e-04
# X24_1SMOH                   8.592809e-02
# PC40_6AA                    4.991445e-01
# PC40_2AA                    3.704859e-02
# PC401AA                     1.581631e-01
# C2                          3.187439e+00
# C3_1                        1.895047e-01
# C3                          1.104361e-08
# C4_1                        3.009783e-01
# C4                          3.551128e-03
# C3OH                        6.372475e-02
# C5_1                        1.222868e+00
# C5                          2.225754e-03
# C4OH                        7.657477e-02
# C6_1                        1.238473e-01
# C6                          4.064000e-01
# C5OH                        4.802613e-01
# C5_1DC                      1.913361e-02
# C5DC                        2.537733e+00
# C8                          2.825247e-01
# C5MDC                       1.265406e-02
# C9                          6.943188e-03
# C7DC                        6.209445e-01
# C10_2                       4.629265e-01
# C10_1                       2.784640e-01
# C10                         1.955183e-01
# C12_1                       5.090148e-01
# C12                         3.776553e-02
# C14_2                       4.538227e-01
# C14_1                       2.788355e-02
# C14                         1.536593e-02
# C12DC                       1.561359e-05
# C14_2OH                     3.907257e-03
# C14_1OH                     1.718491e-02
# C16_2                       1.817187e-02
# C16_1                       8.426858e-03
# C16                         2.180241e-03
# C16_2OH                     3.111098e-03
# C16_1OH                     6.956831e-02
# C16OH                       1.411363e-01
# C18_2                       2.785733e-02
# C18_1                       1.720162e-03
# C18                         3.815828e-02
# C18_1OH                     1.107574e-02
# Arginine_Average            5.213005e-01
# Betaine_Average             3.804312e+00
# C0_Average                  1.183039e+00
# Choline_Average             3.723572e-01
# Citrate_Average             3.831749e-01
# Creatine_Average            1.847310e-01
# Creatinine_Average          7.351654e-01
# Glucose_Average             7.225334e-02
# Glutamate_Average           4.488795e-01
# Glycine_Average             1.299804e-04
# Histidine_Average           1.882733e-01
# Isoleucine_Average          1.894384e-01
# Lactate_Average             1.102192e-03
# Leucine                     2.049932e-03
# Methionine                  7.109800e-03
# Phenylalanine_Average       1.055440e-01
# Proline                     2.495239e-01
# Pyruvate_Average            7.416791e-01
# Succinate_Average           2.784634e+00
# Threonine_Average           1.370030e+00
# Tyrosine_Average            5.318400e-01
# Valine_Average              3.824685e-01
# X1_Methylhistidine          5.205398e-01
# X2_Hydroxybutyrate          9.225317e-01
# X2_Hydroxyisovaleric_Acid   4.945201e-01
# X3_Hydroxybutyrate          2.555552e-02
# ADP                         1.470206e+00
# AMP                         1.877341e-06
# ATP                         6.085958e-01
# Acetamide                   3.582807e-01
# Acetate                     5.511651e-01
# Acetoacetate                1.749287e+00
# Acetone                     4.698774e-01
# Adenosine                   1.799679e+00
# Alanine                     1.972205e-01
# Aspartate                   2.816728e-02
# Creatine_Phosphate          1.460428e-01
# Cytidine                    9.799535e-03
# Dimethyl_Sulfone            3.816422e-01
# Ethanol                     1.236626e+00
# Ethanolamine                2.773756e+00
# Formate                     1.327636e+00
# Glycerol                    4.420182e-02
# Guanidoacetate              6.226295e-01
# Guanosine                   4.174664e+00
# Hypoxanthine                2.107296e+00
# IMP                         7.099463e-03
# Inosine                     2.096620e-03
# Isopropanol                 6.879457e-05
# Malate                      1.073760e-03
# Malonate                    2.746132e-01
# Mannose                     1.677343e+00
# Methanol                    3.384051e+00
# N_N_Dimethylglycine         8.265869e-02
# Nicotinurate                1.791262e+00
# O_Acetylcarnitine           3.919265e+00
# O_Phosphocholine            6.194719e-01
# Propylene_Glycol            2.086568e-01
# Taurine                     5.401562e-01
# Uridine                     1.868054e-01
# Xanthine                    2.890551e-01
# sn_Glycero_3_Phosphocholine 8.343808e-01
# Beta.Alanine                1.456215e-01
var2$cos2
# Dim.1        Dim.2        Dim.3        Dim.4
# Serine                      1.038332e+00 0.3719918571 2.009709e+00 1.232905e+00
# Putrescine                  7.050707e-02 1.7595712755 5.319883e-02 3.345681e-02
# Trans_Hydroxyproline        9.826863e-01 1.2981364836 4.998372e-02 1.807477e-01
# Asparagine                  8.395780e-01 1.2416111207 7.873262e-02 1.246939e-01
# Glutamine                   1.604354e-01 1.0290712051 1.374649e+00 6.389701e-01
# Alpha_Aminoadipic_Acid      2.445671e-02 0.8200456351 1.099158e-01 6.064023e-01
# Methionine_Sulfoxide        2.489937e-01 0.7973772349 1.880307e+00 3.627041e-01
# Acetyl_Ornithine            4.056276e-06 0.9232077808 1.378695e-02 5.889990e-03
# Citrulline                  1.949951e-03 0.2776584151 1.196650e+00 9.377147e-02
# Asymmetric_Dimethylarginine 6.116541e-02 2.0778852943 5.098550e-01 1.007676e-02
# Total_Dimethylarginine      1.989506e-03 2.0530299549 1.287835e+00 4.723729e-02
# Tryptophan                  1.323989e-01 2.0434018587 7.389153e-01 5.295529e-04
# Kynurenine                  5.082714e-02 0.7582629680 7.581787e-01 4.447390e-01
# Ornithine                   7.971317e-01 0.0843259476 1.748252e+00 2.346128e+00
# Lysine                      6.673718e-01 0.0042952562 1.558908e+00 2.360237e+00
# Spermidine                  3.210423e-01 1.5209284751 3.277810e-02 4.169061e-05
# Spermine                    1.855519e-01 1.3033634451 6.475718e-02 2.389342e-02
# Sarcosine                   8.234542e-01 0.3061839518 1.329608e+00 4.044502e-01
# Methylhistidine             2.679436e-01 0.0790753178 4.922850e-01 1.334233e-01
# Beta_Hydroxybutyric_Acid    4.001442e-01 0.4629293261 2.318376e-01 3.529186e+00
# Alpha_Ketoglutaric_Acid     3.095767e-01 1.3726804649 5.620519e-01 1.115989e-01
# Butyric_acid                2.025071e-01 0.7406686854 3.236002e-01 4.176059e-01
# Propionic_Acid              1.071756e-01 0.0392002407 9.706721e-01 3.134888e-02
# Fumaric_Acid                4.694769e-01 0.8232744449 7.785368e-01 1.312637e-02
# Isobutyric_Acid             1.086691e+00 1.0508075867 2.546409e-03 1.027023e+00
# Hippuric_Acid               2.547466e-01 0.4430513457 1.990032e-01 9.535438e-01
# Methylmalonic_Acid          6.374218e-03 0.2955261787 1.183685e-02 4.313384e-03
# LYSOC14_0                   3.396079e-01 0.3210079541 2.315524e+00 3.922741e-01
# LYSOC16_1                   2.344074e-01 0.1271877791 1.807885e+00 2.150597e-01
# LYSOC16_0                   3.360171e-01 0.0003690667 1.811800e+00 9.356022e-02
# LYSOC17_0                   5.753937e-01 0.0477294890 3.257560e+00 1.773617e-01
# LYSOC18_2                   6.194726e-01 0.9084331442 6.425775e-01 4.212012e-01
# LYSOC18_1                   7.194718e-01 1.0606420032 2.574589e-01 2.608258e-01
# LYSOC18_0                   8.323166e-01 1.1273699722 6.490737e-03 5.256507e-01
# LYSOC20_4                   1.311776e-01 1.8961267751 2.452260e-01 1.919051e-02
# LYSOC20_3                   2.939698e-01 0.4311742091 7.353073e-01 6.381567e-01
# LYSOC24_0                   1.037369e+00 0.0117335167 1.421072e+00 1.356888e-01
# LYSOC26_1                   6.648301e-01 0.4782129602 1.776044e+00 1.014310e-02
# LYSOC26_0                   7.980367e-01 0.0116551347 1.930478e+00 4.908840e-03
# LYSOC28_1                   3.081230e-01 0.8459215298 2.425990e+00 1.172053e-02
# LYSOC28_0                   5.472517e-01 0.1752491650 1.712879e+00 3.415357e-02
# X14_1SMOH                   1.497296e-04 0.7680076378 3.157215e+00 8.027317e-03
# X16_1SM                     3.313939e-01 1.8053185363 4.839206e-01 7.141674e-01
# X16_0SM                     3.541884e-01 1.8303345921 6.839863e-01 4.853438e-01
# X16_1SMOH                   6.536548e-01 0.8450220785 1.234314e+00 2.346536e-01
# X18_1SM                     3.316393e-01 1.6711005051 1.388564e-01 1.031263e+00
# PC32_2AA                    1.409499e+00 0.1989554247 5.766992e-02 2.539115e+00
# X18_0SM                     1.561552e+00 0.3266056619 1.363922e-03 1.878429e+00
# X20_2SM                     1.002842e+00 0.9350569629 2.239741e-01 2.116104e+00
# PC36_0AE                    2.831501e-01 1.6905144027 1.124712e+00 1.880945e-02
# PC36_6AA                    2.328486e-01 1.1241204923 1.981112e+00 1.369909e+00
# PC36_0AA                    6.947192e-01 1.3815538169 8.450058e-01 1.017627e-01
# X22_2SMOH                   2.589985e-01 1.4928140301 2.124940e+00 6.413588e-02
# X22_1SMOH                   5.886411e-03 1.8984076055 1.905436e+00 2.932829e-02
# PC38_6AA                    3.776236e-01 1.0170550096 1.911813e+00 8.112624e-01
# PC38_0AA                    5.251628e-02 1.7984580607 2.089243e+00 1.344207e-02
# PC40_6AE                    3.211900e-01 1.0269184053 2.003619e+00 8.023283e-01
# X24_1SMOH                   8.727722e-01 0.4765927032 1.792426e+00 1.919463e-01
# PC40_6AA                    8.311220e-01 1.0509022146 3.181308e-04 8.046680e-01
# PC40_2AA                    2.687763e-01 0.8883364053 2.958964e+00 1.759346e-01
# PC401AA                     9.332858e-01 0.7667508328 1.238979e+00 3.365708e-01
# C2                          6.504466e-01 0.1164501580 4.499443e-03 9.443119e-02
# C3_1                        5.464457e-02 0.7335396058 8.380648e-02 6.582275e-01
# C3                          2.382315e-01 1.9620763893 1.492915e-02 5.800978e-02
# C4_1                        1.822012e-01 1.2186306281 9.937337e-01 1.366595e-02
# C4                          1.346066e-01 1.8774713612 6.369195e-02 8.711431e-01
# C3OH                        2.137756e-02 0.5616723642 6.772545e-01 1.254453e-01
# C5_1                        2.329029e-01 1.2297887545 8.537178e-04 2.371450e-01
# C5                          4.334355e-01 1.5343115314 4.622178e-02 2.206849e-01
# C4OH                        1.778157e-03 0.1838346117 1.914774e-02 5.756131e-02
# C6_1                        6.592147e-02 0.0578087025 1.844084e+00 1.057872e+00
# C6                          7.964473e-01 0.0809964555 2.292054e-01 2.196004e+00
# C5OH                        2.274578e-01 0.8434728983 9.288361e-03 2.256336e-01
# C5_1DC                      2.050829e-02 0.9950490436 1.334239e-02 5.614694e-01
# C5DC                        1.403464e-02 0.5711252747 8.221685e-03 2.976609e-01
# C8                          1.596044e+00 0.1229117531 8.994941e-02 3.424316e-01
# C5MDC                       8.240940e-01 0.0099799205 1.940602e-01 2.249332e-01
# C9                          1.958838e+00 0.1025261233 4.443819e-02 4.549142e-01
# C7DC                        1.078727e-01 0.0006584638 9.871077e-02 9.179169e-02
# C10_2                       1.631043e+00 0.0399551325 3.417691e-03 3.381623e-01
# C10_1                       2.789275e-02 0.3535637918 1.859863e-01 1.060344e+00
# C10                         1.929627e+00 0.0080610610 9.648724e-02 9.584121e-01
# C12_1                       6.856629e-01 0.3344044428 2.323240e-02 1.609099e+00
# C12                         2.170411e+00 0.0088872739 4.251354e-03 1.379858e+00
# C14_2                       1.735605e+00 0.0320646476 3.609895e-01 4.064404e-01
# C14_1                       2.397005e+00 0.2305616074 1.783078e-03 6.814974e-01
# C14                         2.676642e+00 0.0014189510 7.110257e-03 1.121101e+00
# C12DC                       2.586857e+00 0.0050468414 4.851041e-02 2.615275e-02
# C14_2OH                     2.895813e+00 0.0116853277 3.353919e-02 3.036913e-01
# C14_1OH                     2.606464e+00 0.0169229979 1.052105e-02 1.226685e+00
# C16_2                       2.619457e+00 0.0113635705 1.670754e-01 3.419879e-01
# C16_1                       2.733657e+00 0.0279706516 4.205364e-02 6.159080e-01
# C16                         2.255007e+00 0.0909825740 2.333302e-01 9.940996e-01
# C16_2OH                     2.504108e+00 0.0372716697 1.751763e-01 9.450203e-01
# C16_1OH                     2.462259e+00 0.0275516501 4.743126e-01 4.445197e-01
# C16OH                       1.799720e+00 0.0451696768 2.849258e-02 9.127294e-01
# C18_2                       2.876983e+00 0.0012052720 2.978248e-02 5.395674e-01
# C18_1                       2.854417e+00 0.0317549744 1.775215e-02 6.590824e-01
# C18                         2.264377e+00 0.0911231825 3.846228e-01 2.784139e-01
# C18_1OH                     2.361056e+00 0.0270644680 4.000964e-02 1.363233e+00
# Arginine_Average            9.006637e-01 0.2107118692 1.423946e+00 2.332636e+00
# Betaine_Average             4.212482e-01 0.0168774482 1.706748e-02 1.522765e-01
# C0_Average                  1.025318e+00 0.2906065149 2.980165e-02 1.252370e-01
# Choline_Average             1.204292e+00 0.0852528663 6.226161e-03 1.307319e+00
# Citrate_Average             6.400712e-02 1.4720669443 1.870230e-01 6.840206e-01
# Creatine_Average            1.285371e-01 0.8508912659 1.049551e-01 7.828290e-01
# Creatinine_Average          2.159542e-01 1.0852249868 1.267046e-02 6.204848e-01
# Glucose_Average             6.079749e-02 0.8769830215 3.767025e-03 2.385736e+00
# Glutamate_Average           2.015200e-01 1.1513650822 1.184816e+00 1.522049e-01
# Glycine_Average             3.513936e-01 0.7480999932 2.066905e+00 3.293804e-01
# Histidine_Average           3.470042e-02 1.5389347611 2.253140e-01 8.387907e-01
# Isoleucine_Average          1.672900e-01 1.7351878581 7.187111e-01 3.477496e-03
# Lactate_Average             3.608861e-03 2.0708750509 8.632730e-02 6.662004e-02
# Leucine                     1.647208e-03 0.8675168454 1.603273e+00 5.774215e-01
# Methionine                  2.003704e-01 1.2421420945 2.357014e+00 4.604038e-01
# Phenylalanine_Average       1.060531e-01 1.3970435413 1.968634e+00 7.919742e-01
# Proline                     4.869756e-01 1.9571692589 3.443851e-01 2.387893e-03
# Pyruvate_Average            9.911058e-01 0.0477554869 6.507323e-01 1.065515e+00
# Succinate_Average           1.032266e+00 0.1873725626 5.118129e-01 7.162722e-03
# Threonine_Average           1.856761e-01 1.0570364700 7.155112e-01 2.240210e-01
# Tyrosine_Average            1.170892e-02 0.7658819512 1.443533e+00 7.439755e-02
# Valine_Average              2.853603e-01 1.9992480999 7.234033e-01 6.824889e-05
# X1_Methylhistidine          2.922941e-02 0.1055853305 1.504884e-01 1.090183e-01
# X2_Hydroxybutyrate          3.858891e-01 0.4915736097 3.957490e-01 1.147386e+00
# X2_Hydroxyisovaleric_Acid   6.962756e-02 0.0024980383 1.034608e-02 5.519674e-01
# X3_Hydroxybutyrate          5.200124e-01 0.3177111816 2.687691e-01 3.755111e+00
# ADP                         1.843680e-01 0.9101291291 1.525336e-03 3.778214e-01
# AMP                         1.029278e-01 0.0527384718 5.740992e-02 1.566738e-01
# ATP                         2.266063e-02 0.3138257995 6.621295e-03 3.464722e-02
# Acetamide                   3.350722e-01 1.0055282864 7.789899e-04 8.240503e-02
# Acetate                     1.268891e-01 0.0551092749 6.321463e-01 1.537152e+00
# Acetoacetate                2.937812e-02 0.0367078885 1.336274e-02 1.021291e+00
# Acetone                     3.227149e-01 0.0400312972 1.277020e-01 3.916259e+00
# Adenosine                   6.261812e-01 0.0038886375 5.902336e-01 3.971127e-01
# Alanine                     1.225467e-01 0.1675498938 5.208718e-01 4.800819e-03
# Aspartate                   3.200724e-01 0.7728048560 3.979588e-01 1.481127e+00
# Creatine_Phosphate          2.348983e-01 0.0055126660 5.935734e-02 1.008683e-02
# Cytidine                    1.272622e-01 0.0875444567 5.739374e-01 3.548747e-01
# Dimethyl_Sulfone            2.124684e-01 0.5884007029 1.461806e-01 2.065939e+00
# Ethanol                     1.763605e-01 0.0673569422 1.482391e-02 1.380476e-03
# Ethanolamine                6.796368e-01 0.4909595132 9.925163e-02 1.243689e-02
# Formate                     1.474712e-01 0.0082636967 3.115845e-02 8.670446e-02
# Glycerol                    1.472781e-02 0.0005785567 2.668746e-03 8.864558e-02
# Guanidoacetate              2.874656e-01 0.0597026814 2.984978e-01 8.689677e-02
# Guanosine                   4.875064e-01 0.0824217275 4.794991e-01 1.674295e-03
# Hypoxanthine                7.863425e-01 0.0799238423 1.555527e-03 1.024976e-01
# IMP                         1.837854e-01 0.6284497387 2.754997e-01 8.588987e-01
# Inosine                     6.601600e-01 0.3335654172 1.437135e-02 6.864552e-01
# Isopropanol                 3.917827e-02 0.3224642669 3.890282e-03 9.921166e-01
# Malate                      4.216083e-01 0.0651168026 5.152519e-01 2.178591e+00
# Malonate                    3.407500e-04 0.1130239627 2.720918e-01 1.224076e+00
# Mannose                     2.051904e-02 0.3233668899 1.681061e+00 1.080056e-02
# Methanol                    2.149115e-01 0.0713005795 1.605448e-02 2.174071e-01
# N_N_Dimethylglycine         7.226749e-01 0.0704349134 1.949176e-01 2.746909e-01
# Nicotinurate                2.517112e-02 0.5193436197 1.338617e+00 4.267805e-01
# O_Acetylcarnitine           4.206961e-01 0.0894448397 3.540408e-07 3.497060e-01
# O_Phosphocholine            4.608522e-02 0.6790246881 7.119115e-01 3.895733e-01
# Propylene_Glycol            1.158381e-01 0.1650589882 1.440119e-01 3.903658e+00
# Taurine                     5.531164e-03 0.0115223119 3.908008e-01 8.489194e-01
# Uridine                     5.504069e-03 1.0104121364 1.578150e+00 4.148631e-01
# Xanthine                    7.893562e-01 0.0395416857 7.284829e-02 4.000527e-01
# sn_Glycero_3_Phosphocholine 4.980453e-01 0.1006471757 3.900668e-02 5.848322e-01
# Beta.Alanine                3.571884e-02 0.0071578973 8.369447e-03 7.013507e-02
# Dim.5
# Serine                      1.561342e-02
# Putrescine                  7.129601e-01
# Trans_Hydroxyproline        2.205245e-01
# Asparagine                  2.197076e-01
# Glutamine                   1.084817e-02
# Alpha_Aminoadipic_Acid      5.318364e-02
# Methionine_Sulfoxide        2.168118e-02
# Acetyl_Ornithine            6.880678e-01
# Citrulline                  5.578009e+00
# Asymmetric_Dimethylarginine 1.262598e-02
# Total_Dimethylarginine      6.079881e-02
# Tryptophan                  7.865409e-01
# Kynurenine                  3.838441e-03
# Ornithine                   2.743503e-01
# Lysine                      5.551865e-01
# Spermidine                  1.052919e-01
# Spermine                    1.152506e-01
# Sarcosine                   1.668402e-01
# Methylhistidine             2.950441e-01
# Beta_Hydroxybutyric_Acid    5.686449e-04
# Alpha_Ketoglutaric_Acid     5.582691e-01
# Butyric_acid                1.178988e+00
# Propionic_Acid              3.457371e+00
# Fumaric_Acid                1.801765e+00
# Isobutyric_Acid             9.118018e-03
# Hippuric_Acid               1.530087e-01
# Methylmalonic_Acid          1.291095e-01
# LYSOC14_0                   1.436454e+00
# LYSOC16_1                   4.096988e+00
# LYSOC16_0                   3.406992e+00
# LYSOC17_0                   1.694996e+00
# LYSOC18_2                   2.340215e-01
# LYSOC18_1                   2.569104e-01
# LYSOC18_0                   8.597409e-03
# LYSOC20_4                   2.934688e-01
# LYSOC20_3                   6.750290e-01
# LYSOC24_0                   5.259800e-01
# LYSOC26_1                   1.027277e+00
# LYSOC26_0                   3.359795e-01
# LYSOC28_1                   4.338359e-01
# LYSOC28_0                   1.381655e+00
# X14_1SMOH                   3.308283e-03
# X16_1SM                     5.146872e-02
# X16_0SM                     4.705425e-02
# X16_1SMOH                   8.050216e-01
# X18_1SM                     2.720531e-03
# PC32_2AA                    5.177847e-03
# X18_0SM                     4.970963e-03
# X20_2SM                     4.676388e-02
# PC36_0AE                    1.240396e+00
# PC36_6AA                    2.007679e-01
# PC36_0AA                    1.013390e+00
# X22_2SMOH                   1.502579e-01
# X22_1SMOH                   1.386924e-01
# PC38_6AA                    5.988279e-02
# PC38_0AA                    3.534220e-01
# PC40_6AE                    4.551970e-04
# X24_1SMOH                   8.592809e-02
# PC40_6AA                    4.991445e-01
# PC40_2AA                    3.704859e-02
# PC401AA                     1.581631e-01
# C2                          3.187439e+00
# C3_1                        1.895047e-01
# C3                          1.104361e-08
# C4_1                        3.009783e-01
# C4                          3.551128e-03
# C3OH                        6.372475e-02
# C5_1                        1.222868e+00
# C5                          2.225754e-03
# C4OH                        7.657477e-02
# C6_1                        1.238473e-01
# C6                          4.064000e-01
# C5OH                        4.802613e-01
# C5_1DC                      1.913361e-02
# C5DC                        2.537733e+00
# C8                          2.825247e-01
# C5MDC                       1.265406e-02
# C9                          6.943188e-03
# C7DC                        6.209445e-01
# C10_2                       4.629265e-01
# C10_1                       2.784640e-01
# C10                         1.955183e-01
# C12_1                       5.090148e-01
# C12                         3.776553e-02
# C14_2                       4.538227e-01
# C14_1                       2.788355e-02
# C14                         1.536593e-02
# C12DC                       1.561359e-05
# C14_2OH                     3.907257e-03
# C14_1OH                     1.718491e-02
# C16_2                       1.817187e-02
# C16_1                       8.426858e-03
# C16                         2.180241e-03
# C16_2OH                     3.111098e-03
# C16_1OH                     6.956831e-02
# C16OH                       1.411363e-01
# C18_2                       2.785733e-02
# C18_1                       1.720162e-03
# C18                         3.815828e-02
# C18_1OH                     1.107574e-02
# Arginine_Average            5.213005e-01
# Betaine_Average             3.804312e+00
# C0_Average                  1.183039e+00
# Choline_Average             3.723572e-01
# Citrate_Average             3.831749e-01
# Creatine_Average            1.847310e-01
# Creatinine_Average          7.351654e-01
# Glucose_Average             7.225334e-02
# Glutamate_Average           4.488795e-01
# Glycine_Average             1.299804e-04
# Histidine_Average           1.882733e-01
# Isoleucine_Average          1.894384e-01
# Lactate_Average             1.102192e-03
# Leucine                     2.049932e-03
# Methionine                  7.109800e-03
# Phenylalanine_Average       1.055440e-01
# Proline                     2.495239e-01
# Pyruvate_Average            7.416791e-01
# Succinate_Average           2.784634e+00
# Threonine_Average           1.370030e+00
# Tyrosine_Average            5.318400e-01
# Valine_Average              3.824685e-01
# X1_Methylhistidine          5.205398e-01
# X2_Hydroxybutyrate          9.225317e-01
# X2_Hydroxyisovaleric_Acid   4.945201e-01
# X3_Hydroxybutyrate          2.555552e-02
# ADP                         1.470206e+00
# AMP                         1.877341e-06
# ATP                         6.085958e-01
# Acetamide                   3.582807e-01
# Acetate                     5.511651e-01
# Acetoacetate                1.749287e+00
# Acetone                     4.698774e-01
# Adenosine                   1.799679e+00
# Alanine                     1.972205e-01
# Aspartate                   2.816728e-02
# Creatine_Phosphate          1.460428e-01
# Cytidine                    9.799535e-03
# Dimethyl_Sulfone            3.816422e-01
# Ethanol                     1.236626e+00
# Ethanolamine                2.773756e+00
# Formate                     1.327636e+00
# Glycerol                    4.420182e-02
# Guanidoacetate              6.226295e-01
# Guanosine                   4.174664e+00
# Hypoxanthine                2.107296e+00
# IMP                         7.099463e-03
# Inosine                     2.096620e-03
# Isopropanol                 6.879457e-05
# Malate                      1.073760e-03
# Malonate                    2.746132e-01
# Mannose                     1.677343e+00
# Methanol                    3.384051e+00
# N_N_Dimethylglycine         8.265869e-02
# Nicotinurate                1.791262e+00
# O_Acetylcarnitine           3.919265e+00
# O_Phosphocholine            6.194719e-01
# Propylene_Glycol            2.086568e-01
# Taurine                     5.401562e-01
# Uridine                     1.868054e-01
# Xanthine                    2.890551e-01
# sn_Glycero_3_Phosphocholine 8.343808e-01
# Beta.Alanine                1.456215e-01




#### graph of variables default plot ####
fviz_pca_var(DM_pca_FM2, col.var = "black")
pca18_biplot_black2no43<- fviz_pca_var(DM_pca_FM2, col.var = "black")
pca18_biplot_black2no43
# Open the graphics device
png("pca18_biplot_black2")
dev.off()
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(DM_pca_FM2, col.var = "black", axes = c(2,3))
pca18_biplot_black2_3no43<- fviz_pca_var(DM_pca_FM2, col.var = "black", axes= c(2,3))
pca18_biplot_black2_3no43
# Open the graphics device
png("pca18_biplot_black2_3no43")
dev.off()

fviz_pca_var(DM_pca_FM2, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4no43<- fviz_pca_var(DM_pca_FM2, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4no43
png("pca18_biplot_black3_4no43")
dev.off()

fviz_pca_var(DM_pca_FM2, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5no43<- fviz_pca_var(DM_pca_FM2, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5no43
png("pca18_biplot_black4_5no43")
dev.off()


#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var2$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(DM_pca_FM2, choice = "var", axes = 1:2)
bargraph18_cos2no43<- fviz_cos2(DM_pca_FM2, choice = "var", axes = 1:2)
bargraph18_cos2no43
png("bargraph18_cos2no43")
dev.off()
# Color by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM2, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2no43<-fviz_pca_var(DM_pca_FM2, col.var = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                   repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2no43
png("biplot18_colour_cos2no43")
dev.off()
# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM2, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2no43<- fviz_pca_var(DM_pca_FM2, alpha.var = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                   repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2no43
png("biplot18_black_cos2no43")
dev.off()

#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.



# visualize contributions as corrplot

corrplot(var2$contrib, is.corr=FALSE,
         tl.col = "black", tl.cex = 1, 
         cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
         mar = c(1, 1, 1, 1)
) 
corrplot18_contrib2no43 <- corrplot(var2$contrib, is.corr=FALSE,
                               tl.col = "black", tl.cex = 1, 
                               cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
                               mar = c(1, 1, 1, 1)
)
corrplot18_contrib2no43
pdf("mcorrplot18_contrib2no43", height = 7, width =5)
dev.off()
# Barplot of contributions to a component
# The red dashed line on the graph above indicates the expected average contribution. 
# If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/32 = 3.125%.
# For a given component, a variable with a contribution larger than this cutoff 
# could be considered as important in contributing to the component.

#### contributions of variables ####
# to PC1
contrib18_PC1no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = 1, top = 80,
                              fill = "#6baed6", color = "#2171b5",
                              title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 19, label = expression(bold("PC1")), size = 5)
contrib18_PC1no43
png("contrib18_PC1no43")
dev.off()

# contributions of variables to PC2
contrib18_PC2no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = 2, top = 80,
                              fill = "#6baed6", color = "#2171b5",
                              title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 8.5, x = 19, label = expression(bold("PC2")), size = 5)
contrib18_PC2no43
png("contrib18_PC2no43")
dev.off()
# contributions of variables to PC3
contrib18_PC3no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = 3, top = 40,
                              fill = "#6baed6", color = "#2171b5",
                              title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 10, x = 19, label = expression(bold("PC3")), size = 5)
contrib18_PC3no43
png("contrib18_PC3no43")
dev.off()

# contributions of variables to PC4
contrib18_PC4no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = 4, top = 40,
                              fill = "#6baed6", color = "#2171b5",
                              title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 15, x = 19, label = expression(bold("PC4")), size = 5)
contrib18_PC4no43
png("contrib18_PC4no43")
dev.off()
# contributions of variables to PC5
contrib18_PC5no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = 5, top = 40,
                              fill = "#6baed6", color = "#2171b5",
                              title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 25, x = 19, label = expression(bold("PC5")), size = 5)
contrib18_PC5no43
png("contrib18_PC5no43")
dev.off()
# contributions of variables to PC1 and PC2
contrib18_PC1_PC2no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = c(1,2), top = 80,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 5, x = 18, label = expression(bold("PC1 & PC2")), size = 5)
contrib18_PC1_PC2no43
png("contrib18_PC1_PC2no43")
dev.off()
# contributions of variables to PC2 and PC3
contrib18_PC2_PC3no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = c(2,3), top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 6.5, x = 18, label = expression(bold("PC2 & PC3")), size = 5)
contrib18_PC2_PC3no43
png("contrib18_PC2_PC3no43")
dev.off()
# contributions of variables to PC3 and PC4
contrib18_PC3_PC4no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = c(3,4), top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 18, label = expression(bold("PC3 & PC4")), size = 5)
contrib18_PC3_PC4no43
png("contrib18_PC3_PC4no43")
dev.off()
# contributions of variables to PC4 and PC5
contrib18_PC4_PC5no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = c(4,5), top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 13, x = 18, label = expression(bold("PC4 & PC5")), size = 5)
contrib18_PC4_PC5no43
png("contrib18_PC4_PC5no43")
dev.off()
# contributions of variables to PC1:PC3
contrib18_PC1_to_PC3no43 <- fviz_contrib(DM_pca_FM2, choice = "var", axes = c(1:3), top = 40,
                                     fill = "#6baed6", color = "#2171b5",
                                     title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 4.5, x = 18, label = expression(bold("PC1 to PC3")), size = 5)
contrib18_PC1_to_PC3no43
png("contrib18_PC1_PC3no43")
dev.off()
#### colour variable colors using their contributions ####

# PC1 and PC2
varplot18_contrib_PC1_PC2no43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                          repel = TRUE,  # avoid text overlapping
                                          title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2no43
pdf("varplot18_contrib_PC1_PC2no43")
print()
png("varplot18_contrib_PC1_PC2no43")
print(varplot18_contrib_PC1_PC2no43)
dev.off()

# PC1 and PC2 select var
varplot18_contrib_PC1_PC2_selectno43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE,  # avoid text overlapping
                                                 title = "",
                                                 select.var = list(contrib = 60) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2_selectno43
pdf("varplot18_contrib_PC1_PC2_selectno43")
dev.off()
png("varplot18_contrib_PC1_PC2_selectno43")
print(varplot18_contrib_PC1_PC2_selectno43)
dev.off()

# PC2 and PC3
varplot18_contrib_PC2_PC3no43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                          axes = c(2,3),
                                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                          repel = TRUE,  # avoid text overlapping
                                          title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3no43
pdf("varplot18_contrib_PC2_PC3no43")
dev.off()
png("varplot18_contrib_PC2_PC3no43")
print(varplot18_contrib_PC2_PC3no43)
dev.off()

# PC2 and PC3 select var
varplot18_contrib_PC2_PC3_selectno43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                                 axes = c(2,3),
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE,  # avoid text overlapping
                                                 title = "",
                                                 select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3_selectno43
pdf("varplot18_contrib_PC2_PC3_selectno43")
dev.off()
png("varplot18_contrib_PC2_PC3_selectno43")
print(varplot18_contrib_PC2_PC3_selectno43)
dev.off()
# PC4 and PC3
varplot18_contrib_PC4_PC3no43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                          axes = c(3,4),
                                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                          repel = TRUE,  # avoid text overlapping
                                          title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3no43
pdf("varplot18_contrib_PC4_PC3no43")
dev.off()
png("varplot18_contrib_PC4_PC3no43")
print(varplot18_contrib_PC4_PC3no43)
dev.off()

# PC4 and PC3 select var
varplot18_contrib_PC4_PC3_selectno43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                                 axes = c(3,4),
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE,  # avoid text overlapping
                                                 title = "",
                                                 select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3_selectno43
pdf("varplot18_contrib_PC4_PC3_selectno43")
dev.off()
png("varplot18_contrib_PC4_PC3_selectno43")
print(varplot18_contrib_PC4_PC3_selectno43)
dev.off()


# PC4 and PC5
varplot18_contrib_PC4_PC5no43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                          axes = c(4,5),
                                          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                          repel = TRUE,  # avoid text overlapping
                                          title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5no43
pdf("varplot18_contrib_PC4_PC5no43")
dev.off()
png("varplot18_contrib_PC4_PC5no43")
print(varplot18_contrib_PC4_PC5no43)
dev.off()

# PC4 and PC5 select var
varplot18_contrib_PC4_PC5_selectno43 <- fviz_pca_var(DM_pca_FM2, col.var = "contrib",  # colour by contributions to the PC
                                                 axes = c(4,5),
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE,  # avoid text overlapping
                                                 title = "",
                                                 select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5_selectno43
pdf("varplot18_contrib_PC4_PC5_selectno43")
dev.off()
png("varplot18_contrib_PC4_PC5_selectno43")
print(varplot18_contrib_PC4_PC5_selectno43)
dev.off()

# install the colour brewer
library("RColorBrewer") # did not use this colours

#### cos2 selected variance trying plot ####

fviz_pca_var(DM_pca_FM2, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             select.var = list(cos2 = 50)
)
cos2_18_selectno43 <- fviz_pca_var(DM_pca_FM2, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               select.var = list(cos2 = 50)
)
pdf("cos2_18_select")
dev.off()
png("cos2_18_select")
print(cos2_18_selectno43)
dev.off()



#### colour variables by groups ####
#did not load the sheet of the metabolites
# PC1-2
#options(ggrepel.max.overlaps = Inf)
varplot18_met_PC1_PC2_select <- fviz_pca_var(DM_pca_FM, col.var = met_groups$Category, 
                                             # palette = c("dodgerblue2", "#E31A1C", # red
                                             #             "green4",
                                             #             "#6A3D9A", # purple
                                             #             "#FF7F00", # orange
                                             #             "black", "gold1",
                                             #             "skyblue2", "#FB9A99", # lt pink
                                             #             "palegreen2",
                                             #             "#CAB2D6", # lt purple
                                             #             "#FDBF6F", # lt orange
                                             #             "gray70", "khaki2",
                                             #             "maroon", "orchid1", "deeppink1", "blue1", "steelblue4"),
                                             select.var = list(contrib = 40),
                                             legend.title = "Metabolite function",
                                             title = "",
                                             repel = TRUE) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13, face = "bold"))
varplot18_met_PC1_PC2_select
pdf("varplot18_met_PC1_PC2_select")
dev.off()
png("varplot18_met_PC1_PC2_select")
dev.off()

# PC2-3 select
#options(ggrepel.max.overlaps = Inf)
varplot18_met_PC3_PC2_select <- fviz_pca_var(DM_pca_FM, col.var = met_groups$Category,
                                             axes = c(2,3),
                                             # palette = c("dodgerblue2", "#E31A1C", # red
                                             #             "green4",
                                             #             "#6A3D9A", # purple
                                             #             "#FF7F00", # orange
                                             #             "black", "gold1",
                                             #             "skyblue2", "#FB9A99", # lt pink
                                             #             "palegreen2",
                                             #             "#CAB2D6", # lt purple
                                             #             "#FDBF6F", # lt orange
                                             #             "gray70", "khaki2",
                                             #             "maroon", "orchid1", "deeppink1", "blue1", "steelblue4"),
                                             select.var = list(contrib = 50),
                                             legend.title = "Metabolite function",
                                             title = "",
                                             repel = TRUE) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13, face = "bold"))

varplot18_met_PC3_PC2_select
pdf("varplot18_met_PC3_PC2_select")
dev.off()
png("varplot18_met_PC3_PC2_select")
dev.off()

# cPC3-4 select
#options(ggrepel.max.overlaps = Inf)
varplot18_met_PC3_PC4_select <- fviz_pca_var(DM_pca_FM, col.var = met_groups$Category,
                                             axes = c(3,4),
                                             # palette = c("dodgerblue2", "#E31A1C", # red
                                             #             "green4",
                                             #             "#6A3D9A", # purple
                                             #             "#FF7F00", # orange
                                             #             "black", "gold1",
                                             #             "skyblue2", "#FB9A99", # lt pink
                                             #             "palegreen2",
                                             #             "#CAB2D6", # lt purple
                                             #             "#FDBF6F", # lt orange
                                             #             "gray70", "khaki2",
                                             #             "maroon", "orchid1", "deeppink1", "blue1", "steelblue4"),
                                             select.var = list(contrib = 80),
                                             legend.title = "Metabolite function",
                                             title = "",
                                             repel = TRUE) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13, face = "bold"))
varplot18_met_PC3_PC4_select
pdf("varplot18_met_PC3_PC4_select")
dev.off()
png("varplot18_met_PC3_PC4_select")
dev.off()

# PC4-5
#options(ggrepel.max.overlaps = Inf)
varplot18_met_PC4_PC5_select <- fviz_pca_var(DM_pca_FM, col.var = met_groups$Category,
                                             axes = c(4,5),
                                             # palette = c("dodgerblue2", "#E31A1C", # red
                                             #             "green4",
                                             #             "#6A3D9A", # purple
                                             #             "#FF7F00", # orange
                                             #             "black", "gold1",
                                             #             "skyblue2", "#FB9A99", # lt pink
                                             #             "palegreen2",
                                             #             "#CAB2D6", # lt purple
                                             #             "#FDBF6F", # lt orange
                                             #             "gray70", "khaki2",
                                             #             "maroon", "orchid1", "deeppink1", "blue1", "steelblue4"),
                                             select.var = list(contrib = 80),
                                             legend.title = "Metabolite function",
                                             title = "",
                                             repel = TRUE) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13, face = "bold"))
varplot18_met_PC4_PC5_select
pdf("varplot18_met_PC4_PC5_select")
dev.off()
png("varplot18_met_PC4_PC5_select")
dev.off()
#### Dimension description ####

# shows the variables that are signficantly correlated to the dimensions (pos or neg)
# (correlation cooeficient is signficantly different from 0)
# by default, the first 3 dimentions are characterized, use "axes =" to set dims

dim.desc2 <- dimdesc(DM_pca_FM2, axes = c(1:5))

dim.desc2$Dim.1
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# C14_2OH                       0.9292079 1.067592e-23
# C18_2                         0.9261819 2.991665e-23
# C18_1                         0.9225424 9.762914e-23
# C16_1                         0.9028168 2.493373e-20
# C14                           0.8933524 2.371962e-19
# C16_2                         0.8837578 1.893824e-18
# C14_1OH                       0.8815633 2.968810e-18
# C12DC                         0.8782412 5.766782e-18
# C16_2OH                       0.8640805 7.984859e-17
# C16_1OH                       0.8568297 2.742529e-16
# C14_1                         0.8453997 1.682359e-15
# C18_1OH                       0.8390363 4.339899e-15
# C18                           0.8216786 4.733481e-14
# C16                           0.8199769 5.900187e-14
# C12                           0.8044491 3.981774e-13
# C9                            0.7642350 2.770472e-11
# C10                           0.7585153 4.732762e-11
# C16OH                         0.7325378 4.528735e-10
# C14_2                         0.7193713 1.291295e-09
# C10_2                         0.6973652 6.561312e-09
# C8                            0.6898427 1.106834e-08
# X18_0SM                       0.6823479 1.835197e-08
# PC32_2AA                      0.6482762 1.534632e-07
# Choline_Average               0.5992300 2.130051e-06
# Serine                        0.5564110 1.517378e-05
# LYSOC24_0                     0.5561531 1.534181e-05
# Succinate_Average             0.5547834 1.626325e-05
# C0_Average                    0.5529132 1.760422e-05
# X20_2SM                       0.5468194 2.271483e-05
# Pyruvate_Average              0.5436103 2.592637e-05
# PC401AA                       0.5275152 4.933089e-05
# Arginine_Average              0.5182138 7.050605e-05
# X24_1SMOH                     0.5101267 9.538203e-05
# LYSOC18_0                     0.4981635 1.471166e-04
# PC40_6AA                      0.4978059 1.489982e-04
# C5MDC                         0.4956967 1.605494e-04
# Sarcosine                     0.4955042 1.616431e-04
# LYSOC26_0                     0.4877969 2.114500e-04
# Ornithine                     0.4875203 2.134734e-04
# C6                            0.4873109 2.150163e-04
# Hypoxanthine                  0.4842097 2.390919e-04
# N_N_Dimethylglycine           0.4641936 4.631993e-04
# LYSOC18_1                     0.4631637 4.787134e-04
# C12_1                         0.4521504 6.765640e-04
# Ethanolamine                  0.4501591 7.193491e-04
# Lysine                        0.4460788 8.147195e-04
# LYSOC26_1                     0.4452285 8.359697e-04
# C2                            0.4403860 9.667911e-04
# Adenosine                     0.4320934 1.234040e-03
# LYSOC18_2                     0.4297726 1.319837e-03
# LYSOC17_0                     0.4142001 2.047236e-03
# LYSOC28_0                     0.4039440 2.703353e-03
# X3_Hydroxybutyrate            0.3937626 3.532641e-03
# sn_Glycero_3_Phosphocholine   0.3853560 4.378937e-03
# Guanosine                     0.3812570 4.852736e-03
# Fumaric_Acid                  0.3741405 5.782804e-03
# Malate                        0.3545539 9.191041e-03
# Betaine_Average               0.3544024 9.223033e-03
# O_Acetylcarnitine             0.3541701 9.272302e-03
# Beta_Hydroxybutyric_Acid      0.3454108 1.130386e-02
# PC38_6AA                      0.3355500 1.403896e-02
# Glycine_Average               0.3236865 1.806312e-02
# LYSOC14_0                     0.3182120 2.022777e-02
# LYSOC16_0                     0.3165252 2.093743e-02
# Acetamide                     0.3160799 2.112826e-02
# X18_1SM                       0.3144565 2.183641e-02
# Acetone                       0.3101967 2.379054e-02
# PC40_6AE                      0.3094630 2.414156e-02
# Alpha_Ketoglutaric_Acid       0.3038168 2.699116e-02
# LYSOC28_1                     0.3031026 2.737088e-02
# LYSOC20_3                     0.2960595 3.136102e-02
# Guanidoacetate                0.2927659 3.338671e-02
# PC40_2AA                      0.2830891 3.997552e-02
# Methylhistidine               0.2826502 4.029804e-02
# Methionine_Sulfoxide          0.2724720 4.839905e-02
# Hippuric_Acid                -0.2756017 4.577792e-02
# X22_2SMOH                    -0.2778922 4.393409e-02
# PC36_0AE                     -0.2905601 3.480323e-02
# Valine_Average               -0.2916919 3.407033e-02
# Aspartate                    -0.3089241 2.440212e-02
# Spermidine                   -0.3093918 2.417583e-02
# X16_1SM                      -0.3143402 2.188792e-02
# X16_0SM                      -0.3249712 1.758469e-02
# X2_Hydroxybutyrate           -0.3392024 1.296621e-02
# C5                           -0.3594925 8.199009e-03
# Proline                      -0.3810493 4.877892e-03
# X16_1SMOH                    -0.4414707 9.359864e-04
# Inosine                      -0.4436620 8.764310e-04
# PC36_0AA                     -0.4551267 6.168855e-04
# Xanthine                     -0.4851367 2.316503e-04
# Asparagine                   -0.5003319 1.361652e-04
# Trans_Hydroxyproline         -0.5412963 2.849637e-05
# Isobutyric_Acid              -0.5692206 8.678984e-06
dim.desc2$Dim.2
# 
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Asymmetric_Dimethylarginine   0.7829553 4.298760e-12
# Lactate_Average               0.7816334 4.932451e-12
# Total_Dimethylarginine        0.7782584 6.977321e-12
# Tryptophan                    0.7764314 8.397302e-12
# Valine_Average                0.7679970 1.931983e-11
# C3                            0.7608239 3.819668e-11
# Proline                       0.7598719 4.173877e-11
# X22_1SMOH                     0.7483778 1.180267e-10
# LYSOC20_4                     0.7479281 1.227866e-10
# C4                            0.7442397 1.692936e-10
# X16_0SM                       0.7348377 3.747628e-10
# X16_1SM                       0.7297987 5.659963e-10
# PC38_0AA                      0.7284107 6.330513e-10
# Putrescine                    0.7204927 1.183782e-09
# Isoleucine_Average            0.7154832 1.739936e-09
# PC36_0AE                      0.7062128 3.474425e-09
# X18_1SM                       0.7021460 4.666991e-09
# Histidine_Average             0.6738082 3.207601e-08
# C5                            0.6727953 3.423070e-08
# Spermidine                    0.6698546 4.128318e-08
# X22_2SMOH                     0.6636346 6.093847e-08
# Citrate_Average               0.6590069 8.094138e-08
# Phenylalanine_Average         0.6419943 2.206238e-07
# PC36_0AA                      0.6384253 2.701677e-07
# Alpha_Ketoglutaric_Acid       0.6363718 3.032122e-07
# Spermine                      0.6200960 7.349552e-07
# Trans_Hydroxyproline          0.6188514 7.848199e-07
# Methionine                    0.6053573 1.570716e-06
# Asparagine                    0.6052279 1.580957e-06
# C5_1                          0.6023396 1.826404e-06
# C4_1                          0.5996008 2.091524e-06
# Glutamate_Average             0.5828177 4.672712e-06
# LYSOC18_0                     0.5767126 6.191499e-06
# PC36_6AA                      0.5758808 6.430695e-06
# Creatinine_Average            0.5658301 1.008494e-05
# LYSOC18_1                     0.5593847 1.335562e-05
# Threonine_Average             0.5584331 1.391425e-05
# PC40_6AA                      0.5568104 1.491696e-05
# Isobutyric_Acid               0.5567853 1.493296e-05
# Glutamine                     0.5509966 1.908383e-05
# PC40_6AE                      0.5504199 1.955096e-05
# PC38_6AA                      0.5477702 2.183638e-05
# Uridine                       0.5459784 2.351894e-05
# C5_1DC                        0.5418117 2.790439e-05
# X20_2SM                       0.5252248 5.391807e-05
# Acetyl_Ornithine              0.5218863 6.130904e-05
# LYSOC18_2                     0.5176934 7.190700e-05
# PC40_2AA                      0.5119351 8.920847e-05
# Glucose_Average               0.5086532 1.006992e-04
# Leucine                       0.5059005 1.113649e-04
# Creatine_Average              0.5010294 1.328039e-04
# LYSOC28_1                     0.4995641 1.399547e-04
# X16_1SMOH                     0.4992984 1.412882e-04
# C5OH                          0.4988405 1.436138e-04
# Fumaric_Acid                  0.4928315 1.775508e-04
# Alpha_Aminoadipic_Acid        0.4918642 1.836517e-04
# Methionine_Sulfoxide          0.4850183 2.325894e-04
# Aspartate                     0.4774865 2.999167e-04
# X14_1SMOH                     0.4760022 3.151091e-04
# PC401AA                       0.4756125 3.192113e-04
# Tyrosine_Average              0.4753430 3.220775e-04
# Kynurenine                    0.4729727 3.483051e-04
# Glycine_Average               0.4697924 3.865394e-04
# Butyric_acid                  0.4674532 4.170471e-04
# C3_1                          0.4651981 4.485063e-04
# O_Phosphocholine              0.4475782 7.784249e-04
# IMP                           0.4305875 1.289120e-03
# Dimethyl_Sulfone              0.4166417 1.913672e-03
# C5DC                          0.4104799 2.266717e-03
# C3OH                          0.4070687 2.486086e-03
# Nicotinurate                  0.3914295 3.751666e-03
# X2_Hydroxybutyrate            0.3808206 4.905734e-03
# LYSOC26_1                     0.3756097 5.578943e-03
# X24_1SMOH                     0.3749729 5.666523e-03
# Beta_Hydroxybutyric_Acid      0.3695588 6.460971e-03
# Hippuric_Acid                 0.3615374 7.816309e-03
# LYSOC20_3                     0.3566585 8.756252e-03
# Serine                        0.3312782 1.538919e-02
# C10_1                         0.3229684 1.833531e-02
# C12_1                         0.3140958 2.199644e-02
# Inosine                       0.3137016 2.217251e-02
# X18_0SM                       0.3104117 2.368852e-02
# Mannose                       0.3088687 2.442902e-02
# LYSOC14_0                     0.3077401 2.498293e-02
# X3_Hydroxybutyrate            0.3061558 2.577826e-02
# Sarcosine                     0.3005505 2.876463e-02
# Methylmalonic_Acid            0.2952733 3.183504e-02
# C0_Average                    0.2928052 3.336192e-02
# Citrulline                   -0.2862079 3.774424e-02
# ATP                          -0.3042780 2.674830e-02
# Isopropanol                  -0.3084374 2.463950e-02
# Ethanolamine                 -0.3805827 4.934846e-03
# ADP                          -0.5181765 7.060575e-05
# Acetamide                    -0.5446573 2.483523e-05
dim.desc2$Dim.3
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Methionine                    0.5977662 2.288712e-06
# Glycine_Average               0.5597712 1.313473e-05
# Serine                        0.5519716 1.831730e-05
# Phenylalanine_Average         0.5463019 2.320654e-05
# Methionine_Sulfoxide          0.5339058 3.836117e-05
# Ornithine                     0.5148163 8.012468e-05
# Mannose                       0.5048263 1.158005e-04
# Leucine                       0.4930080 1.764583e-04
# Uridine                       0.4891300 2.019432e-04
# Lysine                        0.4861390 2.238421e-04
# Tyrosine_Average              0.4678036 4.123434e-04
# Arginine_Average              0.4646190 4.569246e-04
# Glutamine                     0.4565056 5.908825e-04
# Nicotinurate                  0.4504829 7.122308e-04
# Sarcosine                     0.4489646 7.461700e-04
# Total_Dimethylarginine        0.4418556 9.252687e-04
# Glutamate_Average             0.4238143 1.565053e-03
# Fumaric_Acid                  0.3435499 1.178169e-02
# Kynurenine                    0.3390283 1.301568e-02
# Tryptophan                    0.3346937 1.430117e-02
# Valine_Average                0.3311620 1.542744e-02
# Isoleucine_Average            0.3300862 1.578519e-02
# Threonine_Average             0.3293506 1.603389e-02
# O_Phosphocholine              0.3285211 1.631833e-02
# Pyruvate_Average              0.3140880 2.199991e-02
# Acetate                       0.3095701 2.409003e-02
# Adenosine                     0.2991315 2.956482e-02
# Cytidine                      0.2949731 3.201757e-02
# Alpha_Ketoglutaric_Acid       0.2919029 3.393511e-02
# Alanine                       0.2810060 4.152543e-02
# Malate                        0.2794860 4.268731e-02
# Succinate_Average             0.2785517 4.341458e-02
# Asymmetric_Dimethylarginine   0.2780184 4.383425e-02
# Methylhistidine               0.2731860 4.779050e-02
# X16_1SM                      -0.2708553 4.980023e-02
# LYSOC18_2                    -0.3121138 2.289361e-02
# C3OH                         -0.3204248 1.932771e-02
# X16_0SM                      -0.3220134 1.870273e-02
# LYSOC20_3                    -0.3338756 1.455560e-02
# PC36_0AA                     -0.3579151 8.505239e-03
# Propionic_Acid               -0.3836070 4.575883e-03
# C4_1                         -0.3881372 4.081087e-03
# PC36_0AE                     -0.4129246 2.120254e-03
# Citrulline                   -0.4259257 1.473859e-03
# X16_1SMOH                    -0.4325766 1.216817e-03
# PC401AA                      -0.4333933 1.188201e-03
# LYSOC24_0                    -0.4641499 4.638472e-04
# LYSOC28_0                    -0.5095814 9.731911e-05
# LYSOC26_1                    -0.5188921 6.871752e-05
# X24_1SMOH                    -0.5212797 6.274778e-05
# LYSOC16_1                    -0.5235228 5.757718e-05
# LYSOC16_0                    -0.5240893 5.633458e-05
# C6_1                         -0.5287380 4.703145e-05
# X22_1SMOH                    -0.5374615 3.327792e-05
# PC38_6AA                     -0.5383602 3.209538e-05
# LYSOC26_0                    -0.5409817 2.886336e-05
# PC36_6AA                     -0.5480305 2.160139e-05
# PC40_6AE                     -0.5511347 1.897342e-05
# PC38_0AA                     -0.5627878 1.152320e-05
# X22_2SMOH                    -0.5675754 9.336812e-06
# LYSOC14_0                    -0.5924816 2.957682e-06
# LYSOC28_1                    -0.6064496 1.486702e-06
# PC40_2AA                     -0.6697613 4.152791e-08
# X14_1SMOH                    -0.6918346 9.651978e-09
# LYSOC17_0                    -0.7027429 4.470574e-09

dim.desc2$Dim.4
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Acetone                    0.6671361 4.899787e-08
# Propylene_Glycol           0.6660620 5.240397e-08
# X3_Hydroxybutyrate         0.6532662 1.143334e-07
# Beta_Hydroxybutyric_Acid   0.6333096 3.595789e-07
# C6                         0.4995688 1.399310e-04
# Malate                     0.4975843 1.501751e-04
# Dimethyl_Sulfone           0.4845488 2.363452e-04
# C12_1                      0.4276320 1.403664e-03
# Acetate                    0.4179624 1.844709e-03
# Aspartate                  0.4102748 2.279399e-03
# C12                        0.3960007 3.333222e-03
# C18_1OH                    0.3936079 3.546806e-03
# C14_1OH                    0.3733751 5.891566e-03
# Malonate                   0.3729778 5.948723e-03
# X2_Hydroxybutyrate         0.3611051 7.895891e-03
# C14                        0.3569449 8.698490e-03
# C10_1                      0.3471381 1.087536e-02
# C6_1                       0.3467331 1.097455e-02
# Isobutyric_Acid            0.3416401 1.229003e-02
# Acetoacetate               0.3406854 1.255114e-02
# C16                        0.3361195 1.386683e-02
# C10                        0.3300312 1.580369e-02
# Hippuric_Acid              0.3291919 1.608798e-02
# C16_2OH                    0.3277173 1.659804e-02
# C16OH                      0.3220697 1.868091e-02
# C4                         0.3146470 2.175231e-02
# IMP                        0.3124279 2.274943e-02
# Taurine                    0.3106076 2.359586e-02
# Creatine_Average           0.2982718 3.005860e-02
# X16_1SM                    0.2848910 3.867348e-02
# Inosine                    0.2793090 4.282432e-02
# Citrate_Average            0.2788132 4.320998e-02
# C14_1                      0.2782985 4.361341e-02
# C18_1                      0.2736835 4.737023e-02
# C3_1                       0.2735060 4.751987e-02
# Phenylalanine_Average     -0.3000090 2.906779e-02
# PC40_6AE                  -0.3019638 2.798566e-02
# PC40_6AA                  -0.3024037 2.774680e-02
# PC38_6AA                  -0.3036403 2.708456e-02
# Histidine_Average         -0.3087490 2.448728e-02
# Isopropanol               -0.3357841 1.396798e-02
# X18_1SM                   -0.3423446 1.210037e-02
# Pyruvate_Average          -0.3479835 1.067080e-02
# Serine                    -0.3743205 5.757482e-03
# Choline_Average           -0.3854513 4.368419e-03
# PC36_6AA                  -0.3945705 3.459477e-03
# X18_0SM                   -0.4620364 4.962323e-04
# X20_2SM                   -0.4903965 1.932735e-04
# Arginine_Average          -0.5148755 7.994706e-05
# Ornithine                 -0.5163624 7.560737e-05
# Lysine                    -0.5179128 7.131346e-05
# Glucose_Average           -0.5207029 6.414468e-05
# PC32_2AA                  -0.5371802 3.365624e-05
dim.desc2$Dim.5
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Citrulline           0.6671114 4.907386e-08
# Guanosine            0.5771245 6.076100e-06
# LYSOC16_1            0.5717302 7.757795e-06
# Betaine_Average      0.5509305 1.913679e-05
# LYSOC16_0            0.5213679 6.253671e-05
# Succinate_Average    0.4713492 3.673716e-04
# Ethanolamine         0.4704277 3.786111e-04
# Hypoxanthine         0.4100356 2.294275e-03
# Fumaric_Acid         0.3791472 5.113705e-03
# Adenosine            0.3789276 5.141564e-03
# Nicotinurate         0.3780405 5.255459e-03
# LYSOC17_0            0.3677419 6.748534e-03
# ADP                  0.3424900 1.206153e-02
# LYSOC14_0            0.3385359 1.315655e-02
# PC36_0AE             0.3145857 2.177932e-02
# PC36_0AA             0.2843460 3.906360e-02
# X2_Hydroxybutyrate   0.2712998 4.941173e-02
# LYSOC26_1           -0.2862876 3.768862e-02
# Butyric_acid        -0.3066998 2.550281e-02
# C0_Average          -0.3072262 2.523859e-02
# C5_1                -0.3123551 2.278278e-02
# Ethanol             -0.3141073 2.199136e-02
# Formate             -0.3254605 1.740531e-02
# Threonine_Average   -0.3306160 1.560814e-02
# LYSOC28_0           -0.3320157 1.514848e-02
# Mannose             -0.3658219 7.064463e-03
# Acetoacetate        -0.3735849 5.861581e-03
# C5DC                -0.4499681 7.235794e-04
# C2                  -0.5042893 1.180778e-04
# Methanol            -0.5196097 6.687082e-05
# Propionic_Acid      -0.5252086 5.395192e-05
# O_Acetylcarnitine   -0.5591922 1.346691e-05


#### extract the results for individuals ####
ind2 <- get_pca_ind(DM_pca_FM2)
ind2
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"
ind2$coord
# Dim.1        Dim.2        Dim.3       Dim.4       Dim.5
# 1  15.27677740   6.74811184   2.86444294 -4.28013141  1.36612307
# 2   2.26302098  -5.28362818  -1.96331988  1.04511617  3.39559858
# 3   6.49440445  -2.87490835  -3.27593713  4.21986268  7.34897915
# 4   4.32885441  -3.36534775   1.65269020 -2.22586811 -2.64322449
# 5   5.26602648   4.73995021   8.21308832 -6.60964148 -0.08259330
# 6   4.36808628  -0.32688872   7.26676117 -2.42135124 -0.75893389
# 7   2.83064775  -1.85484452   1.24844281 -2.35969413  1.25418817
# 8  13.51237476   3.66776822   2.80150295 -1.59238124 -0.36794939
# 9   9.75340198  -0.29563744   5.41831673  1.07615518  1.89462766
# 10 14.89972405   0.46939713   1.88275317 12.01642337 -0.46980809
# 11 -1.73244159  -7.21182841   6.20553765 -1.42044524  1.30144766
# 12  0.02890034  -6.20435543  -3.94140253 -1.01379096  3.82481899
# 13  3.92521432  -0.98503474  -7.33550350 -1.05137824  5.31523480
# 14 -3.10486116  -8.27766507   0.78664760  0.24853558  2.46256500
# 15 -2.34145062 -10.35547912   0.88583910  1.81982871  2.54396818
# 16  4.74086995   0.61745693   0.74204692 -3.58769396 -0.26102206
# 17 -3.11495726  -9.76244289  -2.02719789 -0.34219442  4.29019884
# 18  1.19822982  -5.64358563  -5.51532546 -0.74261901  1.64181795
# 19 -7.75276951 -13.17611476   0.98255347 -0.25445903  0.55804365
# 20  5.62410853   0.97614406 -13.31272205 -2.26227006  1.34877646
# 21 -1.90742199  -6.25244379   4.72143307 -6.04056281 -1.91500485
# 22  5.96647831  -0.33192398   2.65596150 -3.87554185  0.58965146
# 23 -2.86742408  -0.20253517   0.48058878 -3.05616258  0.32907876
# 24  1.54052835  -3.30834298   0.91261908 -2.61564401  0.17153446
# 25 -3.52008114   2.04527174  -1.98221416 -2.60387521  0.50905774
# 26 -3.72251319  -2.38813334   1.01225393  0.68965352 -3.06553809
# 27  1.23432057   6.46374821  -8.72963882 -3.97462930 -4.09182571
# 28 -0.64322109   0.04820938  -6.58604779 -1.83550209 -4.36473347
# 29 -5.53235017  -6.61349712   1.98821988  0.21171878 -5.13663350
# 30 -1.83575151  -0.69114999   0.07280662 -0.06593171 -0.47890966
# 31 -3.86082888   0.56201845  -0.92543856 -1.89225303 -4.87491142
# 32 -1.70599339   2.47431902  -2.05950153 -1.06450837 -2.66471538
# 33  3.65058661   1.99662741  -1.59316603  2.42414161 -3.69892670
# 34 -1.10588888  -0.19218675   0.69787663 -3.06806001 -1.04838423
# 35  2.74889696  -3.58695503   0.93937680  9.65805705 -7.67718049
# 36  5.61257503   0.68463122  -2.98408217  1.01539338 -2.45154467
# 37 -3.12912874  -3.75446928   1.64597872  2.34629153 -2.37081118
# 38  0.20317923  -2.31243820  -0.23670227  4.60947605 -2.08435659
# 39  1.98258928  -2.75748551  -0.68483837  2.35003405 -3.84036131
# 40 -6.04860122   6.13955294   3.57609051 -0.41444758  0.91953130
# 41 -3.45365310   6.57134386  -1.77062335 -1.18810881 -1.93222134
# 42 -3.43609568  11.35830615  -1.40477960  0.25931868  0.60846658
# 43 -4.02748703  11.16384608  -2.58136816  2.36420591  1.30351371
# 44 -3.11824345  -0.35420062   0.79704050  2.75556759  1.89423113
# 45 -4.22067049   7.10527621   1.87608793  7.42287249  1.83592620
# 46 -5.25588275  -1.95133641  -1.31924866 -0.04531593 -0.96114521
# 47 -5.49048439   2.58815863   0.92231748 -0.36089093 -0.73862012
# 48 -5.88537742   1.14933394   2.43280167  1.37071376  2.03514416
# 49 -6.27391577   7.39056369  -0.08965333  3.58517150  3.37269255
# 50 -5.42803255   0.19484943  -2.53280761  0.02753816 -0.12494819
# 51 -3.39802341   7.27420772  -3.16448233 -1.29012959  0.01238801
# 52 -7.71052755   8.76589517   6.86299324  1.85951586  3.97763381
# 53 -5.82571783   9.11987154   3.47093179  0.17989073  1.999065                                                                               
ind2$cos2
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# 1  6.180385e-01 1.205915e-01 2.172863e-02 4.851391e-02 4.942330e-03
# 2  5.124640e-02 2.793513e-01 3.857166e-02 1.092988e-02 1.153769e-01
# 3  2.288936e-01 4.485420e-02 5.824065e-02 9.663884e-02 2.930955e-01
# 4  1.971617e-01 1.191616e-01 2.873819e-02 5.212851e-02 7.350969e-02
# 5  1.283142e-01 1.039576e-01 3.121201e-01 2.021457e-01 3.156446e-05
# 6  1.539356e-01 8.620979e-04 4.260283e-01 4.730121e-02 4.646914e-03
# 7  1.090754e-01 4.683498e-02 2.121742e-02 7.579956e-02 2.141316e-02
# 8  6.856881e-01 5.052043e-02 2.947443e-02 9.522639e-03 5.084397e-04
# 9  5.259247e-01 4.832036e-04 1.623081e-01 6.402670e-03 1.984536e-02
# 10 4.983495e-01 4.946050e-04 7.957275e-03 3.241366e-01 4.954714e-04
# 11 1.933719e-02 3.350943e-01 2.481047e-01 1.299946e-02 1.091264e-02
# 12 5.799046e-06 2.672662e-01 1.078580e-01 7.135880e-03 1.015716e-01
# 13 9.702646e-02 6.110356e-03 3.388623e-01 6.961157e-03 1.779133e-01
# 14 7.690567e-02 5.466253e-01 4.936674e-03 4.927778e-04 4.837817e-02
# 15 2.756604e-02 5.391936e-01 3.945609e-03 1.665196e-02 3.254076e-02
# 16 2.258583e-01 3.831187e-03 5.533281e-03 1.293453e-01 6.846591e-04
# 17 4.890872e-02 4.803957e-01 2.071450e-02 5.902391e-04 9.277628e-02
# 18 1.241916e-02 2.755003e-01 2.631201e-01 4.770281e-03 2.331642e-02
# 19 2.146045e-01 6.198687e-01 3.446969e-03 2.311856e-04 1.111889e-03
# 20 1.101264e-01 3.317507e-03 6.170463e-01 1.781857e-02 6.333794e-03
# 21 2.328599e-02 2.502078e-01 1.426753e-01 2.335372e-01 2.347151e-02
# 22 2.672367e-01 8.270603e-04 5.295457e-02 1.127521e-01 2.610061e-03
# 23 6.194545e-02 3.090483e-04 1.740094e-03 7.036852e-02 8.158780e-04
# 24 2.871792e-02 1.324447e-01 1.007842e-02 8.278866e-02 3.560544e-04
# 25 6.751247e-02 2.279191e-02 2.140818e-02 3.694188e-02 1.411929e-03
# 26 1.092897e-01 4.498046e-02 8.081382e-03 3.751186e-03 7.411739e-02
# 27 7.424047e-03 2.035887e-01 3.713448e-01 7.698002e-02 8.158663e-02
# 28 3.849766e-03 2.162604e-05 4.036115e-01 3.134899e-02 1.772675e-01
# 29 2.185220e-01 3.122757e-01 2.822310e-02 3.200329e-04 1.883792e-01
# 30 4.192307e-02 5.942509e-03 6.594273e-05 5.407718e-05 2.853202e-03
# 31 1.482228e-01 3.140906e-03 8.516262e-03 3.560510e-02 2.363128e-01
# 32 3.857460e-02 8.114432e-02 5.621743e-02 1.501913e-02 9.411273e-02
# 33 9.162498e-02 2.740834e-02 1.745061e-02 4.040216e-02 9.406760e-02
# 34 1.453204e-02 4.388848e-04 5.787100e-03 1.118488e-01 1.306004e-02
# 35 2.723129e-02 4.636638e-02 3.180028e-03 3.361482e-01 2.124001e-01
# 36 2.712289e-01 4.035757e-03 7.667136e-02 8.877287e-03 5.174772e-02
# 37 8.906980e-02 1.282274e-01 2.464518e-02 5.007807e-02 5.113021e-02
# 38 3.947729e-04 5.113626e-02 5.357887e-04 2.031851e-01 4.154634e-02
# 39 5.066214e-02 9.800415e-02 6.044970e-03 7.118137e-02 1.900913e-01
# 40 2.429139e-01 2.502741e-01 8.491013e-02 1.140464e-03 5.614040e-03
# 41 7.047579e-02 2.551473e-01 1.852401e-02 8.340556e-03 2.205953e-02
# 42 5.556107e-02 6.071098e-01 9.286594e-03 3.164519e-04 1.742262e-03
# 43 6.130482e-02 4.710359e-01 2.518409e-02 2.112501e-02 6.421811e-03
# 44 1.094180e-01 1.411781e-03 7.148740e-03 8.544583e-02 4.037706e-02
# 45 9.893064e-02 2.803689e-01 1.954674e-02 3.059933e-01 1.871882e-02
# 46 2.697122e-01 3.717693e-02 1.699271e-02 2.004985e-05 9.019600e-03
# 47 3.146461e-01 6.991701e-02 8.878957e-03 1.359418e-03 5.694339e-03
# 48 3.186002e-01 1.215038e-02 5.443909e-02 1.728188e-02 3.809674e-02
# 49 2.037387e-01 2.827166e-01 4.160343e-05 6.652978e-02 5.887755e-02
# 50 2.538430e-01 3.270985e-04 5.526940e-02 6.533562e-06 1.345056e-04
# 51 8.438296e-02 3.866993e-01 7.318252e-02 1.216378e-02 1.121515e-06
# 52 1.961620e-01 2.535357e-01 1.554081e-01 1.140897e-02 5.220309e-02
# 53 1.350678e-01 3.310018e-01 4.794519e-02 1.287865e-04 1.590402e-02
ind2$contrib
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# 1  1.476836e+01 2.912301e+00 1.021184e+00 3.041450e+00 4.413531e-01
# 2  3.240754e-01 1.785403e+00 4.797399e-01 1.813409e-01 2.726712e+00
# 3  2.668994e+00 5.285910e-01 1.335657e+00 2.956400e+00 1.277205e+01
# 4  1.185809e+00 7.243220e-01 3.399434e-01 8.225558e-01 1.652246e+00
# 5  1.754830e+00 1.436876e+00 8.395309e+00 7.253078e+00 1.613230e-03
# 6  1.207400e+00 6.833948e-03 6.572120e+00 9.733793e-01 1.362117e-01
# 7  5.070380e-01 2.200325e-01 1.939818e-01 9.244384e-01 3.719906e-01
# 8  1.155400e+01 8.603508e-01 9.768005e-01 4.209792e-01 3.201716e-02
# 9  6.019787e+00 5.589728e-03 3.653865e+00 1.922721e-01 8.488957e-01
# 10 1.404835e+01 1.409135e-02 4.411746e-01 2.397270e+01 5.219724e-02
# 11 1.899267e-01 3.326309e+00 4.792725e+00 3.349774e-01 4.005530e-01
# 12 5.285362e-05 2.461870e+00 1.933416e+00 1.706331e-01 3.459620e+00
# 13 9.749798e-01 6.205473e-02 6.697050e+00 1.835205e-01 6.681154e+00
# 14 6.100329e-01 4.382153e+00 7.701656e-02 1.025519e-02 1.434108e+00
# 15 3.469276e-01 6.858232e+00 9.766375e-02 5.498289e-01 1.530488e+00
# 16 1.422279e+00 2.438289e-02 6.853090e-02 2.136965e+00 1.611243e-02
# 17 6.140066e-01 6.095212e+00 5.114651e-01 1.944072e-02 4.352727e+00
# 18 9.085506e-02 2.036958e+00 3.785873e+00 9.155850e-02 6.374653e-01
# 19 3.803493e+00 1.110316e+01 1.201534e-01 1.074985e-02 7.364494e-02
# 20 2.001595e+00 6.093960e-02 2.205756e+01 8.496801e-01 4.302160e-01
# 21 2.302303e-01 2.500181e+00 2.774415e+00 6.057890e+00 8.672540e-01
# 22 2.252708e+00 7.046104e-03 8.779449e-01 2.493626e+00 8.222376e-02
# 23 5.202987e-01 2.623449e-03 2.874557e-02 1.550670e+00 2.560981e-02
# 24 1.501787e-01 6.999916e-01 1.036580e-01 1.135857e+00 6.958399e-03
# 25 7.841050e-01 2.675308e-01 4.890180e-01 1.125659e+00 6.128309e-02
# 26 8.768824e-01 3.647447e-01 1.275271e-01 7.896388e-02 2.222389e+00
# 27 9.641060e-02 2.672026e+00 9.484539e+00 2.622768e+00 3.959504e+00
# 28 2.618119e-02 1.486397e-04 5.398504e+00 5.593406e-01 4.505283e+00
# 29 1.936815e+00 2.797268e+00 4.919857e-01 7.441926e-03 6.239700e+00
# 30 2.132537e-01 3.055036e-02 6.597287e-04 7.216976e-04 5.423926e-02
# 31 9.432568e-01 2.020100e-02 1.065906e-01 5.944632e-01 5.620048e+00
# 32 1.841720e-01 3.915463e-01 5.278955e-01 1.881329e-01 1.679222e+00
# 33 8.433234e-01 2.549564e-01 3.158972e-01 9.756241e-01 3.235624e+00
# 34 7.739125e-02 2.362211e-03 6.061512e-02 1.562767e+00 2.599247e-01
# 35 4.781738e-01 8.228556e-01 1.098255e-01 1.548624e+01 1.393831e+01
# 36 1.993394e+00 2.997679e-02 1.108269e+00 1.711730e-01 1.421301e+00
# 37 6.196062e-01 9.015066e-01 3.371880e-01 9.139669e-01 1.329231e+00
# 38 2.612327e-03 3.419889e-01 6.973142e-03 3.527521e+00 1.027426e+00
# 39 2.487335e-01 4.862933e-01 5.837137e-02 9.168850e-01 3.487793e+00
# 40 2.315148e+00 2.410712e+00 1.591625e+00 2.851713e-02 1.999584e-01
# 41 7.547903e-01 2.761723e+00 3.901899e-01 2.343576e-01 8.829179e-01
# 42 7.471355e-01 8.250855e+00 2.456066e-01 1.116437e-02 8.755481e-02
# 43 1.026448e+00 7.970755e+00 8.293228e-01 9.279768e-01 4.018258e-01
# 44 6.153028e-01 8.023621e-03 7.906504e-02 1.260633e+00 8.485404e-01
# 45 1.127280e+00 3.228745e+00 4.380564e-01 9.147671e+00 7.971077e-01
# 46 1.748076e+00 2.435208e-01 2.166093e-01 3.409325e-04 2.184664e-01
# 47 1.907613e+00 4.284041e-01 1.058728e-01 2.162313e-02 1.290176e-01
# 48 2.191884e+00 8.448201e-02 7.366091e-01 3.119321e-01 9.794831e-01
# 49 2.490842e+00 3.493228e+00 1.000361e-03 2.133961e+00 2.690048e+00
# 50 1.864463e+00 2.428120e-03 7.984139e-01 1.259030e-04 3.692041e-03
# 51 7.306705e-01 3.384100e+00 1.246319e+00 2.763333e-01 3.629190e-05
# 52 3.762158e+00 4.914330e+00 5.862068e+00 5.740719e-01 3.741590e+00
# 53 2.147671e+00 5.319236e+00 1.499394e+00 5.372595e-03 9.450625e-01



fviz_pca_ind(DM_pca_FM2, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE  # avoid text overlapping, slow if many points
)
pca18_indno43 <- 
  fviz_pca_ind(DM_pca_FM2, col.ind = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE  # avoid text overlapping, slow if many points
  )
pca18_indno43
png("pca18_indno43")
print(pca18_indno43)
dev.off()
# change point size according to cos2
fviz_pca_ind(DM_pca_FM2, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_pointno43<- fviz_pca_ind(DM_pca_FM2, pointsize = "cos2", 
                               pointshape = 21, fill = "#E7B800",
                               repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_pointno43
png("pca18_ind_pointno43")
print(pca18_ind_pointno43)
dev.off()
# to change boht pointsize and colour according to cos2
fviz_pca_ind(DM_pca_FM2, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_sizeno43<- fviz_pca_ind(DM_pca_FM2, col.ind = "cos2", pointsize = "cos2",
                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                    repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_sizeno43
png("pca18_ind_point_sizeno43")
print(pca18_ind_point_sizeno43)
dev.off()
# To visualize the contribution of individuals to the first two principal components
fviz_contrib(DM_pca_FM2, choice = "ind", axes = 1:2)
contr18_indno43<- fviz_contrib(DM_pca_FM2, choice = "ind", axes = 1:2)
contr18_indno43
png("contr18_indno43")
print(contr18_indno43)
dev.off()
#### Color by a custom continuous variable ####
# Spine age

fviz_pca_ind(DM_pca_FM2, col.ind = DM_data.frame2$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# trying to do the comparisson of age with the PCA loadins

# Get the scores of the first principal component
pc1_scores18no43 <- DM_pca_FM2$ind$coord[, "Dim.1"]
# second dim
pc1_scores18PC2no43 <- DM_pca_FM2$ind$coord[, "Dim.2"]

# Calculate the correlation
correlation18no43 <- cor(DM_data.frame2$Spine.Age, pc1_scores18no43)

correlation18no43
#[1] 0.5076787
# cortests
correlation18testno43 <- cor.test(DM_data.frame2$Spine.Age, pc1_scores18no43)

correlation18testno43
# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$Spine.Age and pc1_scores18no43
# t = 4.2082, df = 51, p-value = 0.0001044
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2751403 0.6840988
# sample estimates:
#       cor 
# 0.5076787 


correlation18PC2no43 <- cor(DM_data.frame2$Spine.Age, pc1_scores18PC2no43)

correlation18PC2no43
#[1] 0.04068033
correlation18PC2testno43 <- cor.test(DM_data.frame2$Spine.Age, pc1_scores18PC2no43)

correlation18PC2testno43
# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$Spine.Age and pc1_scores18PC2no43
# t = 0.29076, df = 51, p-value = 0.7724
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2321663  0.3075920
# sample estimates:
#        cor 
# 0.04068033 

# weight now
# Calculate the correlation
correlation18wno43 <- cor(DM_data.frame2$W, pc1_scores18no43)

correlation18wno43
#
correlation18wtestno43 <- cor.test(DM_data.frame2$W, pc1_scores18no43)

correlation18wtestno43
# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$W and pc1_scores18no43
# t = 4.6943, df = 50, p-value = 2.117e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3299992 0.7176668
# sample estimates:
#     cor 
# 0.55309 

# Calculate the correlation
correlation18wPC2no43 <- cor(DM_data.frame2$W, pc1_scores18PC2no43)

correlation18wPC2no43
#[1] 0.2850698
correlation18wPC2testno43 <- cor.test(DM_data.frame2$W, pc1_scores18PC2no43)

correlation18wPC2testno43

# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$W and pc1_scores18PC2no43
# t = 1.0141, df = 50, p-value = 0.3154
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1362211  0.3993846
# sample estimates:
#       cor 
# # 0.1419565 


# fl now
# Calculate the correlation
cor(DM_data.frame2$FL, pc1_scores18no43)
correlation18flno43 <- cor(DM_data.frame2$FL, pc1_scores18no43)

correlation18flno43
#[1] 0.4588664
cor.test(DM_data.frame2$FL, pc1_scores18no43)

correlation18fltestno43 <- cor.test(DM_data.frame2$FL, pc1_scores18no43)

correlation18fltestno43

# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$FL and pc1_scores18no43
# t = 3.6882, df = 51, p-value = 0.0005486
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2152726 0.6487025
# sample estimates:
#       cor 
# 0.4588664 

correlation18flPC2no43 <- cor(DM_data.frame2$FL, pc1_scores18PC2no43)

correlation18flPC2no43
#[1] 0.215295

correlation18flPC2testno43 <- cor.test(DM_data.frame2$FL, pc1_scores18PC2no43)

correlation18flPC2testno43
# Pearson's product-moment correlation
# 
# data:  DM_data.frame2$FL and pc1_scores18PC2no43
# t = 0.81731, df = 51, p-value = 0.4176
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1615554  0.3725478
# sample estimates:
#       cor 
# 0.1137036 




fviz_pca_ind(DM_pca_FM2, col.ind = DM_data.frame2$Spine_Age, axes = c(2,3),
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# W

fviz_pca_ind(DM_pca_FM2, col.ind = DM_data.frame2$W,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Weight") 
#### Anova sex and site ####

# Perform the Two-Way ANOVA

# Get the scores of the first principal component
pc1_scores18 <- DM_pca_FM$ind$coord[, "Dim.1"]

# Perform the ANOVA
modelssno43 <- aov(pc1_scores18no43 ~ DM_data.frame2$Sex * DM_data.frame2$Site, data = DM_data.frame2)

# Print the summary of the model
summary(modelssno43)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame2$Sex                      2   83.7   41.86   2.755   0.0746 .  
# DM_data.frame2$Site                     3  762.9  254.30  16.737 2.13e-07 ***
#   DM_data.frame2$Sex:DM_data.frame2$Site  3   65.1   21.70   1.428   0.2473    
# Residuals                              44  668.5   15.19                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# only site model anova
# Perform the ANOVA
modelsno43site <- aov(pc1_scores18no43 ~ DM_data.frame2$Site, data = DM_data.frame2)

# Print the summary of the model
summary(modelsno43site)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame2$Site  3  735.9  245.32   14.24 8.47e-07 ***
#   Residuals           49  844.3   17.23                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsno43site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18no43 ~ DM_data.frame2$Site, data = DM_data.frame2)
# 
# $`DM_data.frame2$Site`
# diff        lwr         upr     p adj
# Matheson-Dauphin River   4.327906  0.2255255  8.43028558 0.0350230
# Red River-Dauphin River 10.020725  5.9183447 14.12310485 0.0000002
# Sandy Bar-Dauphin River  5.350329  0.6337705 10.06688667 0.0204091
# Red River-Matheson       5.692819  1.6617904  9.72384815 0.0025134
# Sandy Bar-Matheson       1.022423 -3.6322081  5.67705429 0.9363894
# Sandy Bar-Red River     -4.670396 -9.3250274 -0.01576499 0.0489370
library(report)
report(models)

# pc 2
# Perform the ANOVA
modelsPC2no43site <- aov(pc1_scores18PC2no43 ~ DM_data.frame2$Site, data = DM_data.frame2)

# Print the summary of the model
summary(modelsPC2no43site)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame2$Site  3  641.1  213.68   11.35 9.08e-06 ***
#   Residuals           49  922.6   18.83                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
report(modelsPC2no43site)
TukeyHSD(modelsPC2no43site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2no43 ~ DM_data.frame2$Site, data = DM_data.frame2)
# 
# $`DM_data.frame2$Site`
# diff        lwr       upr     p adj
# Matheson-Dauphin River  -6.000171 -10.288397 -1.711944 0.0027898
# Red River-Dauphin River -7.559431 -11.847657 -3.271204 0.0001283
# Sandy Bar-Dauphin River -9.585826 -14.516054 -4.655598 0.0000250
# Red River-Matheson      -1.559260  -5.772903  2.654383 0.7590940
# Sandy Bar-Matheson      -3.585655  -8.451151  1.279841 0.2172122
# Sandy Bar-Red River     -2.026395  -6.891891  2.839101 0.6866068


# only sex model anova
# Perform the ANOVA
modelsexno43 <- aov(pc1_scores18no43 ~ DM_data.frame2$Sex, data = DM_data.frame2)

# Print the summary of the model
summary(modelsexno43)
# Df Sum Sq Mean Sq F value Pr(>F)
# DM_data.frame2$Sex  2   83.7   41.86   1.399  0.256
# Residuals          50 1496.5   29.93               

TukeyHSD(modelsexno43)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18no43 ~ DM_data.frame2$Sex, data = DM_data.frame2)
# 
# $`DM_data.frame2$Sex`
# diff       lwr      upr     p adj
# M-F  -0.1629032 -4.600024 4.274218 0.9956742
# UN-F -3.7495508 -9.234289 1.735187 0.2340810
# UN-M -3.5866476 -9.871420 2.698125 0.3597624

# pc2
# Perform the ANOVA
modelsexPC2no43 <- aov(pc1_scores18PC2no43 ~ DM_data.frame2$Sex, data = DM_data.frame2)

# Print the summary ofPC2 the model
summary(modelsexPC2no43)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame2$Sex  2  501.8  250.89   11.81 6.28e-05 ***
#   Residuals          50 1061.8   21.24                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsexPC2no43)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2no43 ~ DM_data.frame2$Sex, data = DM_data.frame2)
# 
# $`DM_data.frame2$Sex`
# diff        lwr       upr     p adj
# M-F  -6.573955 -10.311482 -2.836428 0.0002721
# UN-F  2.700772  -1.919196  7.320740 0.3425092
# UN-M  9.274727   3.980864 14.568590 0.0002873

#### Color by groups Site ####
pca18_PC1_PC2no43 <- fviz_pca_ind(DM_pca_FM2,
                              geom.ind = "point", # show points only (but not "text"),
                              pointsize = 2.75,
                              pointshape = 21,
                              fill.ind = DM_data.frame2$Site, # color by groups
                              col.ind = DM_data.frame2$Site,
                              palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                              legend.title = "Site",
                              mean.point = FALSE,  # removes group mean point
                              title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2no43
png("pca18_PC1_PC2no43")
print(pca18_PC1_PC2no43)
dev.off()

pca18_PC1_PC2_meanno43 <- fviz_pca_ind(DM_pca_FM2,
                                   geom.ind = "point", # show points only (but not "text"),
                                   pointsize = 2.75,
                                   pointshape = 21,
                                   fill.ind = DM_data.frame2$Site, # color by groups
                                   col.ind = DM_data.frame2$Site,
                                   palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                   legend.title = "Site",
                                   title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_meanno43
png("pca18_PC1_PC2_meanno43")
print(pca18_PC1_PC2_meanno43)
dev.off()

pca18_PC2_PC3no43 <- fviz_pca_ind(DM_pca_FM2,
                              axes = c(2,3),
                              geom.ind = "point", # show points only (but not "text"),
                              pointsize = 2.75,
                              pointshape = 21,
                              fill.ind = DM_data.frame2$Site, # color by groups
                              col.ind = DM_data.frame2$Site,
                              palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                              legend.title = "Site",
                              mean.point = FALSE,  # removes group mean point
                              title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3no43
png("pca18_PC2_PC3no43")
print(pca18_PC2_PC3no43)
dev.off()


pca18_PC2_PC3_meanno43 <- fviz_pca_ind(DM_pca_FM2,
                                   axes = c(2,3),
                                   geom.ind = "point", # show points only (but not "text"),
                                   pointshape = 21,
                                   pointsize = 2.75,
                                   fill.ind = DM_data.frame2$Site, # color by groups
                                   col.ind = DM_data.frame2$Site,
                                   palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                   legend.title = "Site",
                                   title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3_meanno43
png("pca18_PC2_PC3_meanno43")
print(pca18_PC2_PC3_meanno43)
dev.off()

pca18_PC1_PC3_meanno43 <- fviz_pca_ind(DM_pca_FM2,
                                   axes = c(1,3),
                                   geom.ind = "point", # show points only (but not "text"),
                                   pointshape = 21,
                                   pointsize = 2.75,
                                   fill.ind = DM_data.frame2$Site, # color by groups
                                   col.ind = DM_data.frame2$Site,
                                   palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                   legend.title = "Site",
                                   title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC3_meanno43
png("pca18_PC1_PC3_meanno43")
print(pca18_PC1_PC3_meanno43)
dev.off()


pca18_PC3_PC4_meanno43 <- fviz_pca_ind(DM_pca_FM2,
                                   axes = c(3,4),
                                   geom.ind = "point", # show points only (but not "text"),
                                   pointshape = 21,
                                   pointsize = 2.75,
                                   fill.ind = DM_data.frame2$Site, # color by groups
                                   col.ind = DM_data.frame2$Site,
                                   palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                   legend.title = "Site",
                                   title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_meanno43
png("pca18_PC3_PC4_meanno43")
print(pca18_PC3_PC4_meanno43)
dev.off()


pca18_PC4_PC5_meanno43 <- fviz_pca_ind(DM_pca_FM2,
                                   axes = c(4,5),
                                   geom.ind = "point", # show points only (but not "text"),
                                   pointshape = 21,
                                   pointsize = 2.75,
                                   fill.ind = DM_data.frame2$Site, # color by groups
                                   col.ind = DM_data.frame2$Site,
                                   palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                   legend.title = "Site",
                                   title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_meanno43




#### biplot of individuals and variables ####
fviz_pca_biplot(DM_pca_FM2, repel = TRUE,
                col.var = "#2E9FDF",  # variables colour
                col.ind = "#696969"  # individuals colour
)


#### colour individuals by site ####
fviz_pca_biplot(DM_pca_FM2,
                geom.ind = "point",
                label = "var",
                col.ind = DM_data.frame2$Site,  # colour by groups
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                axes = c(1,2),
                legend.title = "Site",
                col.var = "black",
                repel = TRUE,
                mean.point = FALSE  # removes group mean point
)

# colour individuals by Site
fviz_pca_biplot(DM_pca_FM2, axes = c(1,2),
                geom.ind = "point",
                col.ind = DM_data.frame2$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 50),  # top contributing var
                legend.title = "Site"
)
pca18_ind_by_siteno43 <- fviz_pca_biplot(DM_pca_FM2, axes = c(1,2),
                                     geom.ind = "point",
                                     col.ind = DM_data.frame2$Site,  # colour by groups
                                     mean.point = FALSE,  # removes group mean point
                                     pointsize = 2.5,
                                     palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                     label = "var",
                                     alpha.var ="contrib",col.var = "black",
                                     repel = TRUE,
                                     select.var = list(contrib = 50),  # top contributing var
                                     legend.title = "Site"
)
pca18_ind_by_siteno43
png("pca18_ind_by_siteno43")

print(pca18_ind_by_siteno43)
dev.off()
# colour individuals by Site
fviz_pca_biplot(DM_pca_FM2, axes = c(2,3),
                geom.ind = "point",
                col.ind = DM_data.frame2$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 10),  # top contributing var
                legend.title = "Site"
)
# colour individuals by Site
fviz_pca_biplot(DM_pca_FM2, axes = c(3,4),
                geom.ind = "point",
                col.ind = DM_data.frame2$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 10),  # top contributing var
                legend.title = "Site"
)
# colour individuals by Site
fviz_pca_biplot(DM_pca_FM2, axes = c(4,5),
                geom.ind = "point",
                col.ind = DM_data.frame2$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 10),  # top contributing var
                legend.title = "Site"
)
#### by sex ####
# Color by groups
pca18_PC1_PC2_sexno43 <- fviz_pca_ind(DM_pca_FM2,
                                  geom.ind = "point", # show points only (but not "text"),
                                  pointsize = 2.75,
                                  pointshape = 21,
                                  fill.ind = DM_data.frame2$Sex, # color by groups
                                  col.ind = DM_data.frame2$Sex,
                                  palette = c("#99d8c9", "#08519c","#E7B800"),
                                  legend.title = "Sex",
                                  mean.point = FALSE,  # removes group mean point
                                  title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_sexno43

pca18_PC1_PC2_mean_sexno43 <- fviz_pca_ind(DM_pca_FM2,
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointsize = 2.75,
                                       pointshape = 21,
                                       fill.ind = DM_data.frame2$Sex, # color by groups
                                       col.ind = DM_data.frame2$Sex,
                                       palette = c("#99d8c9", "#08519c","#E7B800"),
                                       legend.title = "Sex",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_mean_sexno43
png("pca18_PC1_PC2_mean_sexno43")
print(pca18_PC1_PC2_mean_sexno43)
dev.off()

pca18_PC3_PC2_mean_sexno43 <- fviz_pca_ind(DM_pca_FM2, axes = c(2,3),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointsize = 2.75,
                                       pointshape = 21,
                                       fill.ind = DM_data.frame2$Sex, # color by groups
                                       col.ind = DM_data.frame2$Sex,
                                       palette = c("#99d8c9", "#08519c","#E7B800"),
                                       legend.title = "Sex",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC2_mean_sexno43
png("pca18_PC3_PC2_mean_sexno43")
print(pca18_PC3_PC2_mean_sexno43)
dev.off()

pca18_PC3_PC4_mean_sexno43 <- fviz_pca_ind(DM_pca_FM2, axes = c(3,4),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointsize = 2.75,
                                       pointshape = 21,
                                       fill.ind = DM_data.frame2$Sex, # color by groups
                                       col.ind = DM_data.frame2$Sex,
                                       palette = c("#99d8c9", "#08519c","#E7B800"),
                                       legend.title = "Sex",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_mean_sexno43
png("pca18_PC3_PC4_mean_sexno43")
print(pca18_PC3_PC4_mean_sexno43)
dev.off()

pca18_PC4_PC5_mean_sexno43 <- fviz_pca_ind(DM_pca_FM2, axes = c(4,5),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointsize = 2.75,
                                       pointshape = 21,
                                       fill.ind = DM_data.frame2$Sex, # color by groups
                                       col.ind = DM_data.frame2$Sex,
                                       palette = c("#99d8c9", "#08519c","#E7B800"),
                                       legend.title = "Sex",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_mean_sexno43

png("pca18_PC4_PC5_mean_sexno43")
dev.off()

install.packages("harrypotter")
library("harrypotter")
hp_palette <- hp(128, option = "Ravenclaw")  # Choose your house option: "Gryffindor," "Slytherin," "Ravenclaw," or "Hufflepuff"

#### Some graphs ####
library(tidyverse)
theme_set(theme_light())
# All missing values can be filtered out by filtering the `sex` variable
# dat <- palmerpenguins::penguins %>% filter(!is.na(sex))
# by sex
point_plot_FL_W <- raw_data2 %>% 
  ggplot(aes(FL, W, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21) +
  scale_fill_hp(discrete = TRUE, option = "Ravenclaw")
point_plot_FL_W

ggplot(raw_data2, aes(factor(Sex), fill=factor(vs))) +
  geom_bar() +
  scale_fill_hp(discrete = TRUE, option = "Ravenclaw")

ggplot(aes(FL, W, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21) +
  scale_fill_hp(discrete = TRUE, option = "Ravenclaw")

ggplot(data = raw_data2, aes(FL, W, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21) +
  scale_fill_hp(discrete = TRUE, option = "Ravenclaw")



point_plot_W_FL <- raw_data %>% 
  ggplot(aes(W, FL, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_W_FL

point_plot_A_W <- raw_data %>% 
  ggplot(aes(Spine_Age, W, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_W

point_plot_A_FL <- raw_data %>% 
  ggplot(aes(Spine_Age, FL, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_FL

# by site
point_plot_FL_W_site <- raw_data %>% 
  ggplot(aes(FL, W, fill = Site)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_FL_W_site

point_plot_W_FL_site <- raw_data %>% 
  ggplot(aes(W, FL, fill = Site)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_W_FL_site

point_plot_A_W_site <- raw_data %>% 
  ggplot(aes(Spine_Age, W, fill = Site)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_W_site

point_plot_A_FL_site <- raw_data %>% 
  ggplot(aes(Spine_Age, FL, fill = Site)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_FL_site


# by mortality
point_plot_FL_W_M <- raw_data %>% 
  ggplot(aes(FL, W, fill = Mortality)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_FL_W_M

point_plot_W_FL_M <- raw_data %>% 
  ggplot(aes(W, FL, fill = Mortality)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_W_FL_M

point_plot_A_W_M <- raw_data %>% 
  ggplot(aes(Spine_Age, W, fill = Mortality)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_W_M

point_plot_A_FL_M <- raw_data %>% 
  ggplot(aes(Spine_Age, FL, fill = Mortality)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_A_FL_M

# plot_plot is obviously a fun name
boxplot_plot_W_Sex <- raw_data %>% 
  ggplot(aes(x = W, fill = Sex)) +
  geom_boxplot()
boxplot_plot_W_Sex

# plot_plot is obviously a fun name
boxplot_plot_FL_Sex <- raw_data %>% 
  ggplot(aes(x = FL, fill = Sex)) +
  geom_boxplot()
boxplot_plot_FL_Sex

# plot_plot is obviously a fun name
boxplot_plot_A_Sex <- raw_data %>% 
  ggplot(aes(x = Spine_Age, fill = Sex)) +
  geom_boxplot()
boxplot_plot_A_Sex

# plot_plot is obviously a fun name
boxplot_plot_W_Site <- raw_data %>% 
  ggplot(aes(x = W, fill = Site)) +
  geom_boxplot()
boxplot_plot_W_Site

# plot_plot is obviously a fun name
boxplot_plot_FL_Site <- raw_data %>% 
  ggplot(aes(x = FL, fill = Site)) +
  geom_boxplot()
boxplot_plot_FL_Site

# plot_plot is obviously a fun name
boxplot_plot_A_Site <- raw_data %>% 
  ggplot(aes(x = Spine_Age, fill = Site)) +
  geom_boxplot()
boxplot_plot_A_Site

# plot_plot is obviously a fun name
boxplot_plot_W_M <- raw_data %>% 
  ggplot(aes(x = W, fill = Mortality)) +
  geom_boxplot()
boxplot_plot_W_M

# plot_plot is obviously a fun name
boxplot_plot_FL_M <- raw_data %>% 
  ggplot(aes(x = FL, fill = Mortality)) +
  geom_boxplot()
boxplot_plot_FL_M

# plot_plot is obviously a fun name
boxplot_plot_A_M <- raw_data %>% 
  ggplot(aes(x = Spine_Age, fill = Mortality)) +
  geom_boxplot()
boxplot_plot_A_M

library(patchwork)
# by site
p_A_Site <- (point_plot_A_W_site + point_plot_A_FL_site) / boxplot_plot_A_Site
p_A_Site 

# by sex
p_A_Sex <- (point_plot_A_W + point_plot_A_FL) / boxplot_plot_A_Sex
p_A_Sex 

# by Mortality
p_A_M <- (point_plot_A_W_M + point_plot_A_FL_M) / boxplot_plot_A_M
p_A_M

#### random forest try ####
# get data
RF_metabolites18 <- read.csv("C:/Users/user/Desktop/Regression/RF_metabolites18.csv", header=FALSE)
View(RF_metabolites18)
# Assume 'metabolites' is your dataframe
standardized_data <- scale(metabolites)

# Convert back to data frame
standardized_data <- as.data.frame(standardized_data)


# Load necessary libraries
library(randomForest)
set.seed(123)
# Assume 'metabolites' is your dataframe and 'Age' is the column with age data
# Remove 'Age' column from the predictors dataframe
predictors <- RF_metabolites18[, !names(RF_metabolites18) %in% "Age"]

# Run Random Forest
rf18 <- randomForest(predictors, RF_metabolites18$Age, ntree = 500, importance=TRUE,
                     proximity=TRUE)

# Print summary of the model
print(rf18)
# Call:
#   randomForest(x = predictors, y = RF_metabolites18$Age, ntree = 500) 
# Type of random forest: unsupervised
# Number of trees: 500
# No. of variables tried at each split: 12

# Assume 'rf' is your random forest model
randomForest::varImpPlot(rf18)

## Look at variable importance:
round(importance(rf18), 2)
# 1     2 MeanDecreaseAccuracy MeanDecreaseGini
# V1   -1.62 -0.68                -1.88             0.17
# V2   -0.38 -0.96                -1.01             0.27
# V3    2.22 -2.19                -0.38             0.40
# V4    0.60  1.22                 1.11             0.34
# V5    1.36 -0.60                 0.60             0.34
# V6   -2.00  0.69                -1.03             0.29
# V7    0.83  1.86                 1.59             0.31
# V8   -0.53 -0.63                -0.73             0.25
# V9    0.59 -2.03                -1.05             0.30
# V10  -1.61 -1.77                -2.31             0.35
# V11   0.91  0.19                 0.74             0.35
# V12  -0.85 -0.23                -0.67             0.35
# V13   0.28 -0.76                -0.74             0.28
# V14   0.72  0.36                 0.65             0.16
# V15  -0.02 -0.96                -0.65             0.27
# V16  -0.08  0.75                 0.50             0.39
# V17   0.69 -0.81                -0.32             0.30
# V18  -0.76 -1.50                -1.31             0.30
# V19  -0.58 -0.09                -0.33             0.41
# V20  -0.15 -0.53                -0.50             0.33
# V21   0.70 -0.05                 0.70             0.30
# V22  -1.05 -0.36                -0.63             0.36
# V23   1.32 -1.02                 0.14             0.28
# V24  -0.74 -1.82                -1.94             0.28
# V25  -0.29  0.08                -0.18             0.26
# V26   2.17 -0.29                 1.46             0.35
# V27   0.80 -1.07                -0.36             0.40
# V28  -2.01  0.70                -0.88             0.26
# V29   0.16 -0.08                -0.19             0.36
# V30  -0.11 -0.78                -0.57             0.27
# V31  -1.07  0.24                -0.81             0.26
# V32   0.84 -0.77                -0.01             0.41
# V33   2.54 -0.63                 1.41             0.35
# V34   1.36  2.00                 2.67             0.23
# V35   0.11  0.54                 0.46             0.39
# V36   0.20 -0.69                -0.61             0.32
# V37   1.49  0.27                 1.08             0.38
# V38  -0.06 -2.31                -1.74             0.40
# V39   1.91  0.62                 1.79             0.48
# V40   2.26 -0.72                 1.58             0.32
# V41   1.42 -0.85                 0.18             0.39
# V42  -0.16 -2.77                -2.03             0.30
# V43   0.60 -0.38                 0.31             0.34
# V44   0.48  1.64                 1.22             0.32
# V45   0.42 -0.26                 0.30             0.43
# V46  -0.28 -1.05                -0.84             0.31
# V47  -0.13  0.77                 0.59             0.40
# V48   1.24  0.36                 0.99             0.42
# V49   0.73  0.46                 0.82             0.36
# V50   0.57 -0.23                 0.32             0.31
# V51  -0.17 -0.58                -0.44             0.31
# V52  -0.09 -0.12                 0.01             0.41
# V53   0.94 -0.62                 0.15             0.26
# V54  -0.21 -1.41                -1.05             0.39
# V55   1.22 -0.57                 0.86             0.36
# V56   1.38 -0.52                 0.57             0.32
# V57   1.15  0.43                 0.82             0.34
# V58  -0.24 -2.87                -2.32             0.29
# V59   0.12 -0.62                -0.25             0.36
# V60  -0.42 -0.61                -0.80             0.27
# V61   1.36 -1.58                -0.20             0.40
# V62   1.71 -0.56                 0.74             0.24
# V63   0.68  0.72                 0.70             0.30
# V64  -0.08  0.69                 0.32             0.29
# V65   2.90  0.35                 2.79             0.31
# V66  -0.34 -0.69                -0.77             0.31
# V67   1.24 -1.35                -0.09             0.28
# V68  -1.20  1.55                 0.22             0.34
# V69   2.39 -0.01                 1.44             0.47
# V70   1.82 -0.39                 1.18             0.42
# V71   2.05 -1.19                 0.12             0.30
# V72  -0.27 -0.50                -0.61             0.53
# V73   1.65 -0.32                 1.17             0.34
# V74  -1.77 -1.14                -2.13             0.37
# V75   0.97 -0.56                 0.41             0.27
# V76   0.59 -0.26                 0.16             0.30
# V77  -0.54 -1.46                -1.33             0.46
# V78  -0.10 -0.85                -0.86             0.35
# V79  -1.27 -0.77                -1.44             0.31
# V80  -0.17 -0.78                -0.55             0.22
# V81   0.67 -0.62                 0.15             0.27
# V82   0.60  1.31                 1.16             0.31
# V83   2.42 -0.22                 1.74             0.39
# V84  -1.99 -0.91                -1.73             0.28
# V85   2.18 -0.23                 1.49             0.34
# V86   0.84 -0.36                 0.53             0.40
# V87   3.04 -0.77                 2.18             0.39
# V88   3.04 -2.24                 1.56             0.45
# V89   3.38 -1.85                 1.84             0.47
# V90   2.05  0.37                 1.62             0.42
# V91   2.99  0.95                 2.97             0.41
# V92   3.60  0.39                 2.87             0.45
# V93   2.82 -1.02                 2.07             0.29
# V94   2.62  1.11                 2.32             0.50
# V95   3.45  0.13                 2.52             0.49
# V96   3.40 -0.25                 2.88             0.46
# V97   3.32  0.68                 3.12             0.39
# V98   2.84 -0.48                 2.06             0.46
# V99   3.16  0.31                 2.72             0.49
# V100  2.79 -0.75                 2.26             0.42
# V101  2.81 -0.76                 2.06             0.47
# V102 -0.91 -0.33                -0.82             0.42
# V103  0.05  0.18                 0.15             0.33
# V104  2.78 -0.41                 1.72             0.36
# V105 -0.13  0.45                 0.03             0.36
# V106  0.80 -1.46                -0.19             0.38
# V107 -1.52  0.56                -0.65             0.27
# V108 -0.40 -1.69                -1.42             0.35
# V109  1.27 -1.58                -0.22             0.29
# V110  0.41  0.26                 0.61             0.36
# V111  1.88 -0.94                 0.74             0.31
# V112  1.65  0.84                 1.49             0.41
# V113  0.37 -0.35                -0.08             0.28
# V114  1.50  0.61                 1.63             0.32
# V115 -0.58 -1.13                -0.96             0.41
# V116 -0.80  0.58                -0.53             0.31
# V117  0.60 -0.67                -0.30             0.33
# V118 -1.21 -0.25                -1.06             0.39
# V119  1.24 -1.09                 0.44             0.27
# V120 -2.26 -0.31                -2.07             0.36
# V121 -0.72  0.17                -0.35             0.27
# V122  1.09 -0.90                 0.19             0.45
# V123 -1.17 -1.58                -1.81             0.25
# V124 -0.87 -1.47                -1.49             0.31
# V125  0.30  0.16                 0.36             0.37
# V126  1.41 -0.59                 0.59             0.30
# V127 -0.25  0.27                 0.17             0.29
# V128 -2.04 -0.10                -1.81             0.32
# V129 -0.02 -0.33                -0.08             0.22
# V130  0.06 -1.55                -0.93             0.18
# V131 -2.07 -0.78                -1.82             0.33
# V132 -0.84 -0.86                -1.10             0.35
# V133  1.27  0.68                 1.21             0.30
# V134  0.86 -0.72                 0.14             0.22
# V135 -0.78 -1.61                -1.66             0.25
# V136 -1.05 -0.28                -0.70             0.38
# V137  0.75 -0.73                 0.02             0.26
# V138 -0.47  0.19                -0.26             0.20
# V139  0.94 -1.70                -0.59             0.34
# V140  0.06  0.43                 0.35             0.33
# V141 -0.33  2.17                 1.37             0.26
# V142 -0.44 -2.08                -2.10             0.24
# V143  0.96  0.27                 1.20             0.31
# V144 -0.46 -1.45                -1.13             0.26
# V145  0.59  0.74                 0.92             0.28
# V146 -1.58  0.65                -0.81             0.37
# V147 -1.09 -0.83                -1.17             0.41
# V148 -0.51 -1.90                -1.51             0.27
# V149 -1.30 -0.39                -0.96             0.25
# V150  1.99  0.66                 2.05             0.18
# V151  0.23 -0.09                -0.01             0.20
# V152 -0.40 -1.34                -1.21             0.38
# V153 -0.21 -1.47                -0.93             0.28
# V154  0.01 -0.35                -0.07             0.29
# V155  0.76  0.80                 0.90             0.39
# V156  1.60 -1.93                -0.28             0.38
# V157  0.23 -0.60                 0.03             0.23
# V158 -0.44 -1.93                -1.67             0.35
# V159  0.29  0.75                 0.69             0.36
# V160 -0.13 -2.24                -1.51             0.30
# V161 -1.39 -0.31                -1.16             0.23
# V162  0.65 -0.30                 0.17             0.23
# V163 -1.13  0.07                -0.74             0.26
# V164  1.87 -1.44                 0.61             0.23

# Assuming 'rf' is your random forest model
plot(rf18, type = "l", main = deparse(substitute(rf)))

# Assume 'rf' is your random forest model and 'metabolites' is your data
# This will create a partial dependence plot for the variable 'Serine'
randomForest::partialPlot(rf18, predictors, "C6")

MDSplot(rf18, RF_metabolites18$Age)

# Now the proximity matrix is stored in the 'proximity' element of the model
prox_matrix18 <- rf18$proximity

# Convert proximity matrix to a distance matrix
dist_matrix18 <- as.dist(1 - prox_matrix18)

# Create a heatmap
heatmap(as.matrix(dist_matrix18))

# Load necessary library
library(rpart.plot)

# Assume 'rf' is your random forest model
# This will get the first tree from the random forest
tree18 <- randomForest::getTree(rf18, k = 1, labelVar=TRUE)

# Print the tree
print(tree18)


# Assume 'tree' is your decision tree model
rpart.plot::rpart.plot(tree, type = 4, extra = 101)


library(rfviz)
rf_viz(rf18, input = TRUE, imp = TRUE, cmd = TRUE, hl_color = "orange")

capabilities("tcltk")


# Assume 'rf' is your random forest model and 'metabolites' is your data
# This will create a partial dependence plot for the variable 'Serine'
randomForest::partialPlot(rf18, RF_metabolites18, "Serine")

getTree(rf18, k=12, labelVar=FALSE)
