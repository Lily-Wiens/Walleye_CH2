#### NO O Chapter 2 rerun of PCA ####
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

# on my laptop
getwd()
setwd("C:/Users/user/Desktop/Chapter 2 Analysis")


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
No_zeroth18 <- read.csv("C:/Users/Treberg's Lab PC/Desktop/Chapter 2 Analysis/No_zeroth18.csv")
   View(No_zeroth18)

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

raw_data3<- No_zeroth18

# as data.frame
DM_data.frame3 <- data.frame(raw_data3)
names(DM_data.frame3)
dim(DM_data.frame3)
#  54 169

# dataframe with just met data columns
DM_data_pca3 <- DM_data.frame3[,c(10:169)] %>%
  print()

dim(DM_data_pca3) 
#54 160

# metabolite groups does not work as is code from last PCA the complete one
# load in object with met_groups
met_groups <- readxl::read_excel("Metabolite_groups_18.xlsx") %>% 
  print()

summary(DM_data_pca3)
# Serine         Putrescine    Trans_Hydroxyproline   Glutamine    
# Min.   : 23.50   Min.   : 3.62   Min.   :  0.242      Min.   : 0.00  
# 1st Qu.: 71.67   1st Qu.:14.62   1st Qu.:  1.248      1st Qu.: 8.46  
# Median : 93.65   Median :19.50   Median :  4.195      Median :15.85  
# Mean   :105.26   Mean   :21.18   Mean   : 17.229      Mean   :19.80  
# 3rd Qu.:128.00   3rd Qu.:24.45   3rd Qu.: 26.800      3rd Qu.:27.43  
# Max.   :268.00   Max.   :65.20   Max.   :118.000      Max.   :84.00  
# Alpha_Aminoadipic_Acid Methionine_Sulfoxide Acetyl_Ornithine   Citrulline     
# Min.   : 0.542         Min.   : 1.030       Min.   :0.305    Min.   :  2.010  
# 1st Qu.: 1.508         1st Qu.: 5.930       1st Qu.:1.458    1st Qu.:  3.513  
# Median : 3.000         Median : 7.850       Median :1.820    Median : 14.900  
# Mean   : 4.456         Mean   : 8.646       Mean   :2.154    Mean   : 46.011  
# 3rd Qu.: 6.478         3rd Qu.:11.175       3rd Qu.:2.743    3rd Qu.: 39.500  
# Max.   :24.200         Max.   :20.900       Max.   :5.960    Max.   :294.000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine   Tryptophan      Ornithine     
# Min.   :0.0000              Min.   :0.2630         Min.   : 5.98   Min.   : 0.000  
# 1st Qu.:0.4420              1st Qu.:0.8555         1st Qu.:12.65   1st Qu.: 1.450  
# Median :0.5735              Median :1.0950         Median :17.55   Median : 4.065  
# Mean   :0.6091              Mean   :1.1694         Mean   :19.27   Mean   : 5.916  
# 3rd Qu.:0.7943              3rd Qu.:1.4000         3rd Qu.:22.70   3rd Qu.: 7.015  
# Max.   :1.9900              Max.   :3.2900         Max.   :45.20   Max.   :29.500  
# Lysine          Spermidine        Spermine       Sarcosine     
# Min.   :  0.820   Min.   : 3.230   Min.   :26.50   Min.   : 3.460  
# 1st Qu.:  4.755   1st Qu.: 9.643   1st Qu.:36.23   1st Qu.: 9.195  
# Median : 19.200   Median :13.750   Median :45.70   Median :15.500  
# Mean   : 27.468   Mean   :15.433   Mean   :47.13   Mean   :20.583  
# 3rd Qu.: 41.025   3rd Qu.:17.950   3rd Qu.:54.62   3rd Qu.:26.575  
# Max.   :110.000   Max.   :61.000   Max.   :82.70   Max.   :90.800  
# Methylhistidine  Beta_Hydroxybutyric_Acid Alpha_Ketoglutaric_Acid  Butyric_acid  
# Min.   : 2.380   Min.   : 0.736           Min.   :20.80           Min.   :1.280  
# 1st Qu.: 6.955   1st Qu.: 4.575           1st Qu.:55.27           1st Qu.:1.900  
# Median :11.950   Median : 6.695           Median :60.75           Median :2.445  
# Mean   :19.753   Mean   : 8.883           Mean   :58.89           Mean   :2.524  
# 3rd Qu.:29.800   3rd Qu.: 9.595           3rd Qu.:64.22           3rd Qu.:3.053  
# Max.   :85.800   Max.   :50.700           Max.   :71.20           Max.   :4.900  
# Propionic_Acid   Fumaric_Acid   Isobutyric_Acid Hippuric_Acid    
# Min.   :25.40   Min.   :21.20   Min.   :0.890   Min.   :0.01110  
# 1st Qu.:34.67   1st Qu.:40.75   1st Qu.:1.750   1st Qu.:0.03062  
# Median :42.65   Median :53.10   Median :2.530   Median :0.03650  
# Mean   :44.71   Mean   :52.42   Mean   :2.935   Mean   :0.03721  
# 3rd Qu.:52.88   3rd Qu.:62.77   3rd Qu.:3.715   3rd Qu.:0.04353  
# Max.   :85.80   Max.   :99.40   Max.   :6.720   Max.   :0.07400  
# Methylmalonic_Acid   LYSOC14_0       LYSOC16_1        LYSOC16_0       LYSOC17_0    
# Min.   :0.02050    Min.   :1.852   Min.   : 4.323   Min.   :18.97   Min.   :1.195  
# 1st Qu.:0.04188    1st Qu.:3.036   1st Qu.: 8.038   1st Qu.:31.57   1st Qu.:2.106  
# Median :0.06485    Median :3.539   Median :10.398   Median :41.44   Median :2.608  
# Mean   :0.07671    Mean   :3.646   Mean   :10.843   Mean   :42.15   Mean   :2.771  
# 3rd Qu.:0.08380    3rd Qu.:4.294   3rd Qu.:12.441   3rd Qu.:47.61   3rd Qu.:3.279  
# Max.   :0.73400    Max.   :6.116   Max.   :26.612   Max.   :95.67   Max.   :6.248  
# LYSOC18_2        LYSOC18_1        LYSOC18_0       LYSOC20_4        LYSOC20_3     
# Min.   :0.3734   Min.   : 3.478   Min.   :1.389   Min.   : 0.985   Min.   : 6.953  
# 1st Qu.:0.9103   1st Qu.: 7.837   1st Qu.:2.954   1st Qu.: 3.347   1st Qu.: 7.989  
# Median :1.1622   Median : 9.475   Median :4.135   Median : 5.082   Median : 8.931  
# Mean   :1.1990   Mean   : 9.742   Mean   :3.979   Mean   : 5.546   Mean   : 9.007  
# 3rd Qu.:1.4203   3rd Qu.:11.953   3rd Qu.:4.874   3rd Qu.: 6.721   3rd Qu.: 9.942  
# Max.   :2.0876   Max.   :16.755   Max.   :7.102   Max.   :14.805   Max.   :12.636  
# LYSOC24_0       LYSOC26_1        LYSOC26_0        LYSOC28_1        LYSOC28_0     
# Min.   :1.343   Min.   :0.7214   Min.   :0.8358   Min.   :0.9302   Min.   :0.5939  
# 1st Qu.:2.532   1st Qu.:1.2015   1st Qu.:1.3633   1st Qu.:1.5009   1st Qu.:1.0157  
# Median :2.888   Median :1.4402   Median :1.5753   Median :1.6607   Median :1.2035  
# Mean   :3.049   Mean   :1.4443   Mean   :1.6407   Mean   :1.6879   Mean   :1.1967  
# 3rd Qu.:3.462   3rd Qu.:1.6796   3rd Qu.:1.8186   3rd Qu.:1.8913   3rd Qu.:1.3107  
# Max.   :5.966   Max.   :2.2097   Max.   :2.8986   Max.   :2.6944   Max.   :1.9853  
# X14_1SMOH        X16_1SM          X16_0SM         X16_1SMOH         X18_1SM     
# Min.   :1.762   Min.   : 3.842   Min.   : 2.759   Min.   : 2.322   Min.   :1.446  
# 1st Qu.:3.288   1st Qu.: 6.063   1st Qu.: 5.203   1st Qu.: 4.333   1st Qu.:2.401  
# Median :4.163   Median : 7.876   Median : 6.506   Median : 6.615   Median :2.797  
# Mean   :4.271   Mean   : 8.229   Mean   : 7.182   Mean   : 6.677   Mean   :3.026  
# 3rd Qu.:4.865   3rd Qu.: 9.940   3rd Qu.: 8.398   3rd Qu.: 8.251   3rd Qu.:3.582  
# Max.   :8.506   Max.   :20.919   Max.   :17.205   Max.   :14.324   Max.   :5.857  
# PC32_2AA        X18_0SM          X20_2SM          PC36_0AE         PC36_6AA     
# Min.   :15.19   Min.   : 6.595   Min.   : 3.103   Min.   : 4.782   Min.   : 33.39  
# 1st Qu.:23.56   1st Qu.: 9.783   1st Qu.: 6.921   1st Qu.: 9.883   1st Qu.: 62.69  
# Median :30.48   Median :12.134   Median : 9.034   Median :13.371   Median : 79.88  
# Mean   :34.00   Mean   :13.229   Mean   : 9.354   Mean   :13.925   Mean   : 82.93  
# 3rd Qu.:41.68   3rd Qu.:15.823   3rd Qu.:11.081   3rd Qu.:17.455   3rd Qu.:100.71  
# Max.   :70.59   Max.   :29.556   Max.   :21.746   Max.   :26.898   Max.   :187.88  
# PC36_0AA        X22_2SMOH        X22_1SMOH        PC38_6AA         PC38_0AA     
# Min.   : 18.22   Min.   : 3.798   Min.   :2.071   Min.   : 297.1   Min.   : 18.31  
# 1st Qu.: 46.35   1st Qu.: 6.093   1st Qu.:3.363   1st Qu.: 578.1   1st Qu.: 40.32  
# Median : 81.88   Median : 8.335   Median :4.269   Median : 700.6   Median : 55.16  
# Mean   : 88.53   Mean   : 8.991   Mean   :4.361   Mean   : 713.2   Mean   : 56.84  
# 3rd Qu.:127.87   3rd Qu.:11.318   3rd Qu.:5.151   3rd Qu.: 840.6   3rd Qu.: 71.04  
# Max.   :209.69   Max.   :26.864   Max.   :9.969   Max.   :1337.0   Max.   :139.60  
# PC40_6AE        X24_1SMOH         PC40_6AA         PC40_2AA        PC401AA     
# Min.   : 28.17   Min.   : 4.490   Min.   : 27.53   Min.   :1.489   Min.   :1.169  
# 1st Qu.: 56.93   1st Qu.: 9.051   1st Qu.: 79.26   1st Qu.:2.374   1st Qu.:2.264  
# Median : 69.31   Median :10.227   Median :114.12   Median :2.816   Median :2.794  
# Mean   : 72.31   Mean   :11.331   Mean   :121.08   Mean   :2.795   Mean   :2.710  
# 3rd Qu.: 83.23   3rd Qu.:14.213   3rd Qu.:159.76   3rd Qu.:3.159   3rd Qu.:3.097  
# Max.   :150.51   Max.   :21.248   Max.   :245.69   Max.   :4.721   Max.   :4.745  
# C2              C3_1               C3              C4_1        
# Min.   : 3.505   Min.   :0.01720   Min.   :0.4332   Min.   :0.03720  
# 1st Qu.:11.411   1st Qu.:0.03857   1st Qu.:0.9799   1st Qu.:0.05042  
# Median :15.130   Median :0.04515   Median :1.3726   Median :0.06380  
# Mean   :16.520   Mean   :0.04956   Mean   :1.6131   Mean   :0.06774  
# 3rd Qu.:20.302   3rd Qu.:0.06237   3rd Qu.:1.9258   3rd Qu.:0.07852  
# Max.   :37.195   Max.   :0.10650   Max.   :4.6204   Max.   :0.12330  
# C4              C3OH              C5_1               C5         
# Min.   :0.0617   Min.   :0.01140   Min.   :0.00860   Min.   :0.03030  
# 1st Qu.:0.1359   1st Qu.:0.02010   1st Qu.:0.02022   1st Qu.:0.08697  
# Median :0.1817   Median :0.02455   Median :0.02240   Median :0.14810  
# Mean   :0.2239   Mean   :0.02446   Mean   :0.02362   Mean   :0.23941  
# 3rd Qu.:0.2545   3rd Qu.:0.02823   3rd Qu.:0.02792   3rd Qu.:0.30443  
# Max.   :1.0766   Max.   :0.04350   Max.   :0.03890   Max.   :1.39790  
# C4OH              C6_1              C6              C5OH       
# Min.   :0.02170   Min.   :0.1171   Min.   :0.5102   Min.   :0.0662  
# 1st Qu.:0.03560   1st Qu.:0.1603   1st Qu.:0.6616   1st Qu.:0.1496  
# Median :0.04455   Median :0.1818   Median :0.7370   Median :0.1722  
# Mean   :0.06464   Mean   :0.1914   Mean   :0.7994   Mean   :0.1807  
# 3rd Qu.:0.06110   3rd Qu.:0.2150   3rd Qu.:0.8328   3rd Qu.:0.2005  
# Max.   :0.45600   Max.   :0.3187   Max.   :1.9323   Max.   :0.3586  
# C5_1DC             C5DC              C8             C5MDC        
# Min.   :0.06710   Min.   :0.0735   Min.   :0.1043   Min.   :0.02280  
# 1st Qu.:0.08605   1st Qu.:0.1664   1st Qu.:0.1803   1st Qu.:0.02730  
# Median :0.09720   Median :0.1930   Median :0.2081   Median :0.02935  
# Mean   :0.10748   Mean   :0.2098   Mean   :0.2405   Mean   :0.03024  
# 3rd Qu.:0.11790   3rd Qu.:0.2427   3rd Qu.:0.2644   3rd Qu.:0.03347  
# Max.   :0.28850   Max.   :0.5097   Max.   :0.8334   Max.   :0.04400  
# C9               C7DC             C10_2             C10_1       
# Min.   :0.01490   Min.   :0.00610   Min.   :0.03560   Min.   :0.1439  
# 1st Qu.:0.01945   1st Qu.:0.01342   1st Qu.:0.06535   1st Qu.:0.1760  
# Median :0.02405   Median :0.01585   Median :0.07680   Median :0.1926  
# Mean   :0.02546   Mean   :0.03934   Mean   :0.08996   Mean   :0.1909  
# 3rd Qu.:0.02902   3rd Qu.:0.05097   3rd Qu.:0.10410   3rd Qu.:0.2032  
# Max.   :0.07090   Max.   :0.17810   Max.   :0.21070   Max.   :0.2503  
# C10             C12_1              C12              C14_2        
# Min.   :0.1260   Min.   :0.02490   Min.   :0.01780   Min.   :0.01000  
# 1st Qu.:0.1318   1st Qu.:0.04942   1st Qu.:0.02712   1st Qu.:0.01432  
# Median :0.1358   Median :0.05445   Median :0.03095   Median :0.01600  
# Mean   :0.1389   Mean   :0.05732   Mean   :0.03594   Mean   :0.01750  
# 3rd Qu.:0.1424   3rd Qu.:0.06237   3rd Qu.:0.03870   3rd Qu.:0.01847  
# Max.   :0.1930   Max.   :0.17290   Max.   :0.13590   Max.   :0.05620  
# C14_1              C14              C12DC            C14_2OH       
# Min.   :0.03590   Min.   :0.01090   Min.   :0.00860   Min.   :0.00730  
# 1st Qu.:0.06103   1st Qu.:0.02712   1st Qu.:0.01405   1st Qu.:0.01413  
# Median :0.08355   Median :0.03960   Median :0.01760   Median :0.02035  
# Mean   :0.09990   Mean   :0.06046   Mean   :0.01944   Mean   :0.02557  
# 3rd Qu.:0.10287   3rd Qu.:0.06323   3rd Qu.:0.02170   3rd Qu.:0.02950  
# Max.   :0.56820   Max.   :0.56350   Max.   :0.05170   Max.   :0.13220  
# C14_1OH            C16_2             C16_1              C16         
# Min.   :0.00960   Min.   :0.00820   Min.   :0.02410   Min.   :0.02010  
# 1st Qu.:0.01960   1st Qu.:0.01075   1st Qu.:0.07233   1st Qu.:0.07425  
# Median :0.02735   Median :0.01405   Median :0.11485   Median :0.10055  
# Mean   :0.03942   Mean   :0.01590   Mean   :0.16115   Mean   :0.13835  
# 3rd Qu.:0.04330   3rd Qu.:0.01812   3rd Qu.:0.17433   3rd Qu.:0.16910  
# Max.   :0.35090   Max.   :0.05710   Max.   :1.09100   Max.   :1.02600  
# C16_2OH           C16_1OH            C16OH              C18_2        
# Min.   :0.00840   Min.   :0.01230   Min.   :0.004600   Min.   :0.00790  
# 1st Qu.:0.01455   1st Qu.:0.02425   1st Qu.:0.007925   1st Qu.:0.01617  
# Median :0.02045   Median :0.03635   Median :0.010600   Median :0.02075  
# Mean   :0.02648   Mean   :0.04385   Mean   :0.012157   Mean   :0.02930  
# 3rd Qu.:0.02933   3rd Qu.:0.04888   3rd Qu.:0.013375   3rd Qu.:0.03243  
# Max.   :0.17480   Max.   :0.30590   Max.   :0.064300   Max.   :0.18110  
# C18_1              C18             C18_1OH        Arginine_Average
# Min.   :0.02510   Min.   :0.00920   Min.   :0.00770   Min.   : 15.01  
# 1st Qu.:0.08568   1st Qu.:0.02725   1st Qu.:0.01323   1st Qu.: 34.93  
# Median :0.13795   Median :0.03600   Median :0.01820   Median : 61.05  
# Mean   :0.19660   Mean   :0.04684   Mean   :0.02116   Mean   : 64.76  
# 3rd Qu.:0.24088   3rd Qu.:0.05435   3rd Qu.:0.02448   3rd Qu.: 85.88  
# Max.   :1.35790   Max.   :0.27680   Max.   :0.10650   Max.   :181.06  
# Betaine_Average    C0_Average    Choline_Average  Citrate_Average 
# Min.   : 76.99   Min.   :13.79   Min.   : 24.64   Min.   : 51.34  
# 1st Qu.:145.78   1st Qu.:22.38   1st Qu.: 36.25   1st Qu.:139.84  
# Median :189.81   Median :31.41   Median : 44.69   Median :181.28  
# Mean   :200.16   Mean   :33.22   Mean   : 51.43   Mean   :239.50  
# 3rd Qu.:259.12   3rd Qu.:39.31   3rd Qu.: 56.25   3rd Qu.:292.59  
# Max.   :501.12   Max.   :70.55   Max.   :203.31   Max.   :823.62  
# Creatine_Average  Creatinine_Average Glucose_Average   Glutamate_Average
# Min.   :  46.09   Min.   : 0.5705    Min.   :  737.4   Min.   : 26.80   
# 1st Qu.: 250.39   1st Qu.: 5.8881    1st Qu.: 3725.3   1st Qu.: 69.51   
# Median : 381.09   Median : 8.8575    Median : 5749.0   Median : 83.09   
# Mean   : 485.99   Mean   :10.3945    Mean   : 5943.9   Mean   : 92.49   
# 3rd Qu.: 668.27   3rd Qu.:14.0687    3rd Qu.: 7596.0   3rd Qu.:107.41   
# Max.   :1351.12   Max.   :24.0750    Max.   :13459.3   Max.   :247.38   
# Glycine_Average Histidine_Average Isoleucine_Average Lactate_Average
# Min.   :127.8   Min.   :10.27     Min.   : 25.06     Min.   :1345   
# 1st Qu.:332.8   1st Qu.:36.12     1st Qu.: 55.56     1st Qu.:4245   
# Median :450.9   Median :46.17     Median : 80.43     Median :4864   
# Mean   :464.1   Mean   :49.81     Mean   : 98.25     Mean   :4773   
# 3rd Qu.:544.5   3rd Qu.:65.72     3rd Qu.:115.84     3rd Qu.:5597   
# Max.   :988.7   Max.   :92.21     Max.   :403.06     Max.   :6482   
# Leucine         Methionine     Phenylalanine_Average    Proline      
# Min.   : 18.18   Min.   : 3.795   Min.   : 16.73        Min.   : 10.19  
# 1st Qu.: 71.55   1st Qu.:19.231   1st Qu.: 46.50        1st Qu.: 30.41  
# Median :109.88   Median :26.962   Median : 61.59        Median : 50.23  
# Mean   :129.27   Mean   :29.057   Mean   : 66.55        Mean   : 69.92  
# 3rd Qu.:161.44   3rd Qu.:31.406   3rd Qu.: 81.93        3rd Qu.:102.06  
# Max.   :554.75   Max.   :76.463   Max.   :172.88        Max.   :249.31  
# Pyruvate_Average Succinate_Average Threonine_Average Tyrosine_Average
# Min.   : 47.36   Min.   : 3.587    Min.   : 23.27    Min.   : 10.13  
# 1st Qu.:122.91   1st Qu.: 7.726    1st Qu.: 69.42    1st Qu.: 36.11  
# Median :224.06   Median :10.181    Median :109.50    Median : 47.70  
# Mean   :197.47   Mean   :11.230    Mean   :118.85    Mean   : 59.24  
# 3rd Qu.:245.33   3rd Qu.:13.109    3rd Qu.:164.91    3rd Qu.: 61.94  
# Max.   :307.19   Max.   :24.900    Max.   :270.56    Max.   :381.19  
# Valine_Average   X1_Methylhistidine X2_Hydroxybutyrate X2_Hydroxyisovaleric_Acid
# Min.   : 49.99   Min.   : 1.00      Min.   : 2.250     Min.   : 0.130           
# 1st Qu.: 91.88   1st Qu.: 6.25      1st Qu.: 4.817     1st Qu.: 0.880           
# Median :142.88   Median : 7.00      Median : 7.005     Median : 1.130           
# Mean   :176.41   Mean   : 6.90      Mean   :11.132     Mean   : 1.760           
# 3rd Qu.:220.81   3rd Qu.: 7.47      3rd Qu.:14.595     3rd Qu.: 1.938           
# Max.   :627.25   Max.   :14.50      Max.   :52.130     Max.   :11.500           
# X3_Hydroxybutyrate      ADP              AMP             ATP         
# Min.   : 1.500     Min.   :  1.00   Min.   : 23.0   Min.   :  0.800  
# 1st Qu.: 4.660     1st Qu.: 12.00   1st Qu.:143.2   1st Qu.:  2.000  
# Median : 6.380     Median : 39.50   Median :179.0   Median :  2.100  
# Mean   : 8.886     Mean   : 55.00   Mean   :192.1   Mean   :  9.578  
# 3rd Qu.: 9.630     3rd Qu.: 74.25   3rd Qu.:241.5   3rd Qu.:  3.050  
# Max.   :52.880     Max.   :225.00   Max.   :382.0   Max.   :204.600  
# Acetamide         Acetate        Acetoacetate      Acetone        Adenosine    
# Min.   : 1.800   Min.   : 31.00   Min.   : 1.50   Min.   :0.880   Min.   : 10.0  
# 1st Qu.: 7.075   1st Qu.: 63.25   1st Qu.: 3.16   1st Qu.:1.250   1st Qu.: 66.0  
# Median :11.000   Median : 80.00   Median : 5.88   Median :1.750   Median :132.0  
# Mean   :11.911   Mean   : 83.09   Mean   :11.46   Mean   :1.949   Mean   :149.2  
# 3rd Qu.:17.100   3rd Qu.: 92.50   3rd Qu.: 8.25   3rd Qu.:2.098   3rd Qu.:192.0  
# Max.   :25.100   Max.   :163.00   Max.   :86.25   Max.   :7.130   Max.   :451.0  
# Alanine        Aspartate      Creatine_Phosphate    Cytidine    
# Min.   : 51.0   Min.   : 0.000   Min.   : 1.00      Min.   : 3.10  
# 1st Qu.:180.2   1st Qu.: 7.062   1st Qu.: 3.10      1st Qu.: 8.90  
# Median :244.5   Median :10.168   Median : 3.10      Median :11.90  
# Mean   :272.1   Mean   :13.053   Mean   : 5.17      Mean   :12.21  
# 3rd Qu.:340.5   3rd Qu.:15.281   3rd Qu.: 4.50      3rd Qu.:15.32  
# Max.   :857.0   Max.   :44.250   Max.   :28.40      Max.   :23.50  
# Dimethyl_Sulfone    Ethanol        Ethanolamine       Formate         Glycerol     
# Min.   : 3.800   Min.   : 1.000   Min.   : 24.00   Min.   :14.60   Min.   :  66.0  
# 1st Qu.: 8.875   1st Qu.: 3.033   1st Qu.: 59.25   1st Qu.:20.45   1st Qu.: 125.0  
# Median :11.650   Median : 3.750   Median : 90.00   Median :23.95   Median : 145.0  
# Mean   :12.513   Mean   : 4.555   Mean   :104.91   Mean   :26.20   Mean   : 174.3  
# 3rd Qu.:15.875   3rd Qu.: 4.250   3rd Qu.:137.75   3rd Qu.:31.35   3rd Qu.: 177.8  
# Max.   :26.400   Max.   :24.750   Max.   :258.00   Max.   :53.90   Max.   :1318.0  
# Guanidoacetate    Guanosine      Hypoxanthine        IMP            Inosine     
# Min.   : 69.0   Min.   : 12.0   Min.   :11.10   Min.   : 21.90   Min.   : 42.0  
# 1st Qu.:121.5   1st Qu.: 59.0   1st Qu.:18.45   1st Qu.: 60.00   1st Qu.: 99.0  
# Median :167.0   Median :106.0   Median :26.65   Median : 86.85   Median :140.0  
# Mean   :164.9   Mean   :126.2   Mean   :30.64   Mean   :119.34   Mean   :192.2  
# 3rd Qu.:201.8   3rd Qu.:182.0   3rd Qu.:40.70   3rd Qu.:155.10   3rd Qu.:232.0  
# Max.   :280.0   Max.   :505.0   Max.   :73.30   Max.   :420.10   Max.   :717.0  
# Isopropanol       Malate          Malonate        Mannose         Methanol    
# Min.   :0.50   Min.   :  3.10   Min.   :14.00   Min.   : 3.30   Min.   :104.0  
# 1st Qu.:1.00   1st Qu.:  8.80   1st Qu.:30.95   1st Qu.:16.15   1st Qu.:126.0  
# Median :1.13   Median :  9.40   Median :39.80   Median :25.85   Median :141.5  
# Mean   :1.16   Mean   : 16.73   Mean   :42.90   Mean   :26.33   Mean   :153.4  
# 3rd Qu.:1.25   3rd Qu.: 16.40   3rd Qu.:54.88   3rd Qu.:37.30   3rd Qu.:191.2  
# Max.   :1.75   Max.   :102.90   Max.   :92.40   Max.   :61.30   Max.   :225.0  
# N_N_Dimethylglycine  Nicotinurate   O_Acetylcarnitine O_Phosphocholine
# Min.   : 7.10       Min.   : 5.30   Min.   : 3.10     Min.   :15.50   
# 1st Qu.:14.78       1st Qu.:17.23   1st Qu.:11.32     1st Qu.:36.67   
# Median :22.25       Median :28.10   Median :14.50     Median :45.05   
# Mean   :22.95       Mean   :27.79   Mean   :16.21     Mean   :44.76   
# 3rd Qu.:29.02       3rd Qu.:33.98   3rd Qu.:20.30     3rd Qu.:50.75   
# Max.   :48.10       Max.   :52.10   Max.   :33.00     Max.   :72.30   
# Propylene_Glycol    Taurine         Uridine      sn_Glycero_3_Phosphocholine
# Min.   : 0.000   Min.   : 2787   Min.   : 4.10   Min.   : 16.00             
# 1st Qu.: 4.300   1st Qu.: 5383   1st Qu.:13.10   1st Qu.: 38.25             
# Median : 6.100   Median : 6200   Median :17.35   Median : 54.50             
# Mean   : 8.444   Mean   : 6593   Mean   :18.19   Mean   : 72.33             
# 3rd Qu.:10.450   3rd Qu.: 7956   3rd Qu.:22.45   3rd Qu.: 90.25             
# Max.   :50.400   Max.   :11645   Max.   :36.80   Max.   :355.00             
# Beta.Alanine   
# Min.   :  31.0  
# 1st Qu.: 132.5  
# Median : 209.5  
# Mean   : 240.7  
# 3rd Qu.: 277.2  
# Max.   :1065.0  

#### Correlations plot ####

DM_cor3 <- cor(DM_data_pca3)
print(DM_data_pca3)
summary(DM_cor3)
# Serine            Putrescine       Trans_Hydroxyproline   Glutamine       
# Min.   :-0.386034   Min.   :-0.48154   Min.   :-0.51380     Min.   :-0.26786  
# 1st Qu.:-0.003761   1st Qu.: 0.03823   1st Qu.:-0.13546     1st Qu.: 0.06591  
# Median : 0.141126   Median : 0.21402   Median : 0.05812     Median : 0.17812  
# Mean   : 0.177338   Mean   : 0.21317   Mean   : 0.10675     Mean   : 0.20632  
# 3rd Qu.: 0.321278   3rd Qu.: 0.40812   3rd Qu.: 0.33665     3rd Qu.: 0.33369  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000     Max.   : 1.00000  
# Alpha_Aminoadipic_Acid Methionine_Sulfoxide Acetyl_Ornithine     Citrulline      
# Min.   :-0.424566      Min.   :-0.33609     Min.   :-0.37548   Min.   :-0.50286  
# 1st Qu.: 0.008064      1st Qu.:-0.03216     1st Qu.: 0.07943   1st Qu.:-0.22475  
# Median : 0.131629      Median : 0.07444     Median : 0.16611   Median :-0.10054  
# Mean   : 0.151091      Mean   : 0.14018     Mean   : 0.17906   Mean   :-0.07492  
# 3rd Qu.: 0.284638      3rd Qu.: 0.28109     3rd Qu.: 0.30448   3rd Qu.: 0.02880  
# Max.   : 1.000000      Max.   : 1.00000     Max.   : 1.00000   Max.   : 1.00000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine   Tryptophan      
# Min.   :-0.50143            Min.   :-0.5480        Min.   :-0.57112  
# 1st Qu.: 0.02838            1st Qu.: 0.1154        1st Qu.: 0.07541  
# Median : 0.19505            Median : 0.2584        Median : 0.23771  
# Mean   : 0.22509            Mean   : 0.2602        Mean   : 0.24512  
# 3rd Qu.: 0.43318            3rd Qu.: 0.4538        3rd Qu.: 0.46565  
# Max.   : 1.00000            Max.   : 1.0000        Max.   : 1.00000  
# Ornithine            Lysine           Spermidine          Spermine      
# Min.   :-0.44796   Min.   :-0.54111   Min.   :-0.54135   Min.   :-0.5488  
# 1st Qu.:-0.08847   1st Qu.:-0.11068   1st Qu.: 0.02119   1st Qu.: 0.0420  
# Median : 0.07680   Median : 0.03084   Median : 0.19740   Median : 0.1954  
# Mean   : 0.09973   Mean   : 0.04767   Mean   : 0.20036   Mean   : 0.1936  
# 3rd Qu.: 0.22080   3rd Qu.: 0.18268   3rd Qu.: 0.39172   3rd Qu.: 0.3622  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# Sarcosine         Methylhistidine    Beta_Hydroxybutyric_Acid
# Min.   :-0.277369   Min.   :-0.54561   Min.   :-0.33828        
# 1st Qu.: 0.003232   1st Qu.:-0.15077   1st Qu.: 0.08556        
# Median : 0.122199   Median :-0.00401   Median : 0.25856        
# Mean   : 0.160036   Mean   :-0.02288   Mean   : 0.28491        
# 3rd Qu.: 0.284849   3rd Qu.: 0.07610   3rd Qu.: 0.43003        
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000        
# Alpha_Ketoglutaric_Acid  Butyric_acid       Propionic_Acid      Fumaric_Acid    
# Min.   :-0.2785         Min.   :-0.487771   Min.   :-0.38071   Min.   :-0.3852  
# 1st Qu.: 0.1457         1st Qu.: 0.002168   1st Qu.:-0.09943   1st Qu.: 0.0938  
# Median : 0.2493         Median : 0.207143   Median :-0.02205   Median : 0.2232  
# Mean   : 0.2554         Mean   : 0.173959   Mean   : 0.02367   Mean   : 0.2229  
# 3rd Qu.: 0.3757         3rd Qu.: 0.340402   3rd Qu.: 0.10965   3rd Qu.: 0.3389  
# Max.   : 1.0000         Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.0000  
# Isobutyric_Acid    Hippuric_Acid      Methylmalonic_Acid   LYSOC14_0       
# Min.   :-0.57244   Min.   :-0.39532   Min.   :-0.28320   Min.   :-0.30027  
# 1st Qu.:-0.07405   1st Qu.:-0.13846   1st Qu.:-0.02196   1st Qu.:-0.01383  
# Median : 0.09591   Median : 0.02876   Median : 0.04460   Median : 0.13169  
# Mean   : 0.13365   Mean   : 0.03687   Mean   : 0.06850   Mean   : 0.16329  
# 3rd Qu.: 0.37943   3rd Qu.: 0.20473   3rd Qu.: 0.15209   3rd Qu.: 0.28512  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC16_1          LYSOC16_0          LYSOC17_0          LYSOC18_2       
# Min.   :-0.33400   Min.   :-0.36973   Min.   :-0.38989   Min.   :-0.25866  
# 1st Qu.:-0.02307   1st Qu.:-0.11327   1st Qu.:-0.07298   1st Qu.: 0.09511  
# Median : 0.08994   Median : 0.04145   Median : 0.08226   Median : 0.24959  
# Mean   : 0.13128   Mean   : 0.06413   Mean   : 0.12496   Mean   : 0.24801  
# 3rd Qu.: 0.26185   3rd Qu.: 0.19099   3rd Qu.: 0.27686   3rd Qu.: 0.33997  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC18_1         LYSOC18_0         LYSOC20_4          LYSOC20_3        
# Min.   :-0.2270   Min.   :-0.3212   Min.   :-0.40107   Min.   :-0.252566  
# 1st Qu.: 0.1345   1st Qu.: 0.1091   1st Qu.: 0.08196   1st Qu.:-0.008618  
# Median : 0.2447   Median : 0.2314   Median : 0.23651   Median : 0.102721  
# Mean   : 0.2482   Mean   : 0.2427   Mean   : 0.24955   Mean   : 0.126623  
# 3rd Qu.: 0.3494   3rd Qu.: 0.3661   3rd Qu.: 0.45815   3rd Qu.: 0.248172  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.000000  
# LYSOC24_0          LYSOC26_1         LYSOC26_0          LYSOC28_1       
# Min.   :-0.27781   Min.   :-0.2951   Min.   :-0.29045   Min.   :-0.25108  
# 1st Qu.:-0.05198   1st Qu.: 0.0112   1st Qu.:-0.06997   1st Qu.: 0.06226  
# Median : 0.07443   Median : 0.2010   Median : 0.07598   Median : 0.25098  
# Mean   : 0.11268   Mean   : 0.1926   Mean   : 0.09674   Mean   : 0.24468  
# 3rd Qu.: 0.25157   3rd Qu.: 0.2885   3rd Qu.: 0.20756   3rd Qu.: 0.36751  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC28_0           X14_1SMOH            X16_1SM            X16_0SM        
# Min.   :-0.266982   Min.   :-0.328043   Min.   :-0.55117   Min.   :-0.57470  
# 1st Qu.:-0.008469   1st Qu.:-0.004088   1st Qu.:-0.01017   1st Qu.: 0.01725  
# Median : 0.119016   Median : 0.231346   Median : 0.27262   Median : 0.25615  
# Mean   : 0.141458   Mean   : 0.226345   Mean   : 0.25124   Mean   : 0.24610  
# 3rd Qu.: 0.242403   3rd Qu.: 0.431488   3rd Qu.: 0.48283   3rd Qu.: 0.46059  
# Max.   : 1.000000   Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# X16_1SMOH          X18_1SM           PC32_2AA            X18_0SM        
# Min.   :-0.5456   Min.   :-0.2305   Min.   :-0.358308   Min.   :-0.31392  
# 1st Qu.:-0.0586   1st Qu.: 0.1532   1st Qu.: 0.005022   1st Qu.: 0.05831  
# Median : 0.1138   Median : 0.2862   Median : 0.189435   Median : 0.20806  
# Mean   : 0.1306   Mean   : 0.2859   Mean   : 0.185059   Mean   : 0.21630  
# 3rd Qu.: 0.3336   3rd Qu.: 0.4419   3rd Qu.: 0.344561   3rd Qu.: 0.39595  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.000000   Max.   : 1.00000  
# X20_2SM            PC36_0AE           PC36_6AA           PC36_0AA      
# Min.   :-0.30991   Min.   :-0.57279   Min.   :-0.33964   Min.   :-0.6029  
# 1st Qu.: 0.09098   1st Qu.: 0.01689   1st Qu.: 0.07184   1st Qu.:-0.1008  
# Median : 0.25262   Median : 0.20595   Median : 0.23069   Median : 0.1119  
# Mean   : 0.25672   Mean   : 0.22100   Mean   : 0.24516   Mean   : 0.1410  
# 3rd Qu.: 0.39549   3rd Qu.: 0.43347   3rd Qu.: 0.38861   3rd Qu.: 0.3890  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# X22_2SMOH          X22_1SMOH          PC38_6AA           PC38_0AA      
# Min.   :-0.44440   Min.   :-0.3815   Min.   :-0.27067   Min.   :-0.4575  
# 1st Qu.: 0.05437   1st Qu.: 0.1121   1st Qu.: 0.07341   1st Qu.: 0.0881  
# Median : 0.27046   Median : 0.3354   Median : 0.25997   Median : 0.2852  
# Mean   : 0.25646   Mean   : 0.3023   Mean   : 0.26354   Mean   : 0.2673  
# 3rd Qu.: 0.48318   3rd Qu.: 0.4783   3rd Qu.: 0.41848   3rd Qu.: 0.4740  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# PC40_6AE          X24_1SMOH           PC40_6AA          PC40_2AA       
# Min.   :-0.25545   Min.   :-0.28122   Min.   :-0.4514   Min.   :-0.25066  
# 1st Qu.: 0.08587   1st Qu.: 0.05087   1st Qu.: 0.1063   1st Qu.: 0.04538  
# Median : 0.29083   Median : 0.20001   Median : 0.2520   Median : 0.26633  
# Mean   : 0.27507   Mean   : 0.22735   Mean   : 0.2554   Mean   : 0.25412  
# 3rd Qu.: 0.43815   3rd Qu.: 0.39504   3rd Qu.: 0.3930   3rd Qu.: 0.40248  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# PC401AA               C2                C3_1                C3           
# Min.   :-0.27623   Min.   :-0.23804   Min.   :-0.46355   Min.   :-0.603335  
# 1st Qu.: 0.09355   1st Qu.: 0.02102   1st Qu.: 0.04199   1st Qu.: 0.005741  
# Median : 0.25232   Median : 0.18411   Median : 0.19821   Median : 0.220022  
# Mean   : 0.27693   Mean   : 0.19014   Mean   : 0.18439   Mean   : 0.217495  
# 3rd Qu.: 0.42552   3rd Qu.: 0.34393   3rd Qu.: 0.32530   3rd Qu.: 0.406689  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000  
# C4_1                C4                C3OH               C5_1          
# Min.   :-0.42795   Min.   :-0.52913   Min.   :-0.31143   Min.   :-0.492808  
# 1st Qu.: 0.02189   1st Qu.: 0.05597   1st Qu.: 0.05157   1st Qu.: 0.006111  
# Median : 0.18970   Median : 0.31250   Median : 0.15254   Median : 0.174408  
# Mean   : 0.19384   Mean   : 0.28838   Mean   : 0.16256   Mean   : 0.168813  
# 3rd Qu.: 0.37547   3rd Qu.: 0.54636   3rd Qu.: 0.27324   3rd Qu.: 0.342370  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000  
# C5                C4OH               C6_1                C6          
# Min.   :-0.57424   Min.   :-0.28991   Min.   :-0.51031   Min.   :-0.24049  
# 1st Qu.:-0.01002   1st Qu.: 0.03888   1st Qu.:-0.08518   1st Qu.: 0.04368  
# Median : 0.15179   Median : 0.20216   Median : 0.14463   Median : 0.16886  
# Mean   : 0.17941   Mean   : 0.20596   Mean   : 0.11457   Mean   : 0.24574  
# 3rd Qu.: 0.40104   3rd Qu.: 0.37839   3rd Qu.: 0.31913   3rd Qu.: 0.41725  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C5OH             C5_1DC              C5DC                C8         
# Min.   :-0.2440   Min.   :-0.38970   Min.   :-0.38651   Min.   :-0.2348  
# 1st Qu.: 0.1251   1st Qu.: 0.06737   1st Qu.: 0.08054   1st Qu.: 0.1038  
# Median : 0.2162   Median : 0.25003   Median : 0.22949   Median : 0.2280  
# Mean   : 0.2149   Mean   : 0.25293   Mean   : 0.21343   Mean   : 0.2827  
# 3rd Qu.: 0.3076   3rd Qu.: 0.45835   3rd Qu.: 0.38793   3rd Qu.: 0.3772  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# C5MDC                C9               C7DC              C10_2         
# Min.   :-0.23834   Min.   :-0.2322   Min.   :-0.30993   Min.   :-0.36773  
# 1st Qu.: 0.03753   1st Qu.: 0.1086   1st Qu.:-0.04636   1st Qu.:-0.03503  
# Median : 0.14345   Median : 0.2455   Median : 0.03070   Median : 0.11804  
# Mean   : 0.18338   Mean   : 0.2921   Mean   : 0.03770   Mean   : 0.15587  
# 3rd Qu.: 0.27143   3rd Qu.: 0.4087   3rd Qu.: 0.12209   3rd Qu.: 0.28848  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# C10_1               C10               C12_1               C12          
# Min.   :-0.25501   Min.   :-0.21481   Min.   :-0.22622   Min.   :-0.25986  
# 1st Qu.: 0.06977   1st Qu.: 0.06431   1st Qu.: 0.06811   1st Qu.: 0.04019  
# Median : 0.19409   Median : 0.19593   Median : 0.24234   Median : 0.17244  
# Mean   : 0.16168   Mean   : 0.25888   Mean   : 0.28857   Mean   : 0.25845  
# 3rd Qu.: 0.27081   3rd Qu.: 0.34707   3rd Qu.: 0.46327   3rd Qu.: 0.35941  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C14_2              C14_1              C14               C12DC         
# Min.   :-0.26850   Min.   :-0.1729   Min.   :-0.26319   Min.   :-0.36353  
# 1st Qu.: 0.06929   1st Qu.: 0.1100   1st Qu.: 0.05376   1st Qu.: 0.02244  
# Median : 0.20938   Median : 0.2475   Median : 0.18088   Median : 0.19301  
# Mean   : 0.27622   Mean   : 0.3117   Mean   : 0.27235   Mean   : 0.23697  
# 3rd Qu.: 0.43403   3rd Qu.: 0.4569   3rd Qu.: 0.42013   3rd Qu.: 0.40339  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# C14_2OH            C14_1OH             C16_2              C16_1         
# Min.   :-0.30363   Min.   :-0.27917   Min.   :-0.29934   Min.   :-0.25688  
# 1st Qu.: 0.03525   1st Qu.: 0.04596   1st Qu.: 0.02952   1st Qu.: 0.04088  
# Median : 0.19650   Median : 0.17191   Median : 0.19020   Median : 0.18127  
# Mean   : 0.26569   Mean   : 0.26619   Mean   : 0.25126   Mean   : 0.26007  
# 3rd Qu.: 0.37031   3rd Qu.: 0.40714   3rd Qu.: 0.35628   3rd Qu.: 0.35913  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C16               C16_2OH            C16_1OH             C16OH         
# Min.   :-0.289144   Min.   :-0.27602   Min.   :-0.28035   Min.   :-0.24956  
# 1st Qu.: 0.009578   1st Qu.: 0.01963   1st Qu.: 0.03058   1st Qu.: 0.08589  
# Median : 0.160962   Median : 0.17617   Median : 0.17418   Median : 0.20072  
# Mean   : 0.245306   Mean   : 0.25772   Mean   : 0.25950   Mean   : 0.27993  
# 3rd Qu.: 0.400566   3rd Qu.: 0.41341   3rd Qu.: 0.42202   3rd Qu.: 0.43338  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C18_2              C18_1               C18               C18_1OH       
# Min.   :-0.24006   Min.   :-0.25820   Min.   :-0.302980   Min.   :-0.2664  
# 1st Qu.: 0.05139   1st Qu.: 0.03312   1st Qu.: 0.004634   1st Qu.: 0.0279  
# Median : 0.20687   Median : 0.18803   Median : 0.158510   Median : 0.1676  
# Mean   : 0.27992   Mean   : 0.26190   Mean   : 0.238156   Mean   : 0.2591  
# 3rd Qu.: 0.37539   3rd Qu.: 0.38343   3rd Qu.: 0.404416   3rd Qu.: 0.4008  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.0000  
# Arginine_Average   Betaine_Average      C0_Average       Choline_Average   
# Min.   :-0.47774   Min.   :-0.38868   Min.   :-0.23009   Min.   :-0.38147  
# 1st Qu.:-0.03336   1st Qu.:-0.03401   1st Qu.: 0.04219   1st Qu.:-0.04164  
# Median : 0.09012   Median : 0.05010   Median : 0.22594   Median : 0.13387  
# Mean   : 0.12530   Mean   : 0.06714   Mean   : 0.22736   Mean   : 0.13115  
# 3rd Qu.: 0.26569   3rd Qu.: 0.17037   3rd Qu.: 0.36100   3rd Qu.: 0.29763  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Citrate_Average   Creatine_Average   Creatinine_Average Glucose_Average   
# Min.   :-0.4857   Min.   :-0.40498   Min.   :-0.52157   Min.   :-0.36054  
# 1st Qu.: 0.0350   1st Qu.: 0.03295   1st Qu.: 0.01738   1st Qu.:-0.09705  
# Median : 0.2976   Median : 0.20146   Median : 0.18521   Median : 0.08677  
# Mean   : 0.2516   Mean   : 0.18847   Mean   : 0.17793   Mean   : 0.08917  
# 3rd Qu.: 0.4368   3rd Qu.: 0.36124   3rd Qu.: 0.37810   3rd Qu.: 0.25585  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Glutamate_Average  Glycine_Average   Histidine_Average  Isoleucine_Average
# Min.   :-0.24993   Min.   :-0.3285   Min.   :-0.25286   Min.   :-0.53638  
# 1st Qu.: 0.06057   1st Qu.: 0.0679   1st Qu.: 0.06007   1st Qu.: 0.01943  
# Median : 0.18395   Median : 0.1947   Median : 0.22930   Median : 0.16180  
# Mean   : 0.20823   Mean   : 0.2105   Mean   : 0.22729   Mean   : 0.19893  
# 3rd Qu.: 0.34191   3rd Qu.: 0.3380   3rd Qu.: 0.36360   3rd Qu.: 0.40935  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# Lactate_Average      Leucine           Methionine       Phenylalanine_Average
# Min.   :-0.4297   Min.   :-0.36218   Min.   :-0.38762   Min.   :-0.33167     
# 1st Qu.: 0.0939   1st Qu.:-0.01687   1st Qu.: 0.02905   1st Qu.: 0.02377     
# Median : 0.2821   Median : 0.10272   Median : 0.15233   Median : 0.15192     
# Mean   : 0.2557   Mean   : 0.14296   Mean   : 0.19859   Mean   : 0.20510     
# 3rd Qu.: 0.4164   3rd Qu.: 0.28586   3rd Qu.: 0.34507   3rd Qu.: 0.38203     
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000     
# Proline         Pyruvate_Average  Succinate_Average  Threonine_Average
# Min.   :-0.51920   Min.   :-0.3953   Min.   :-0.42636   Min.   :-0.5029  
# 1st Qu.:-0.07782   1st Qu.:-0.0108   1st Qu.: 0.01886   1st Qu.: 0.1129  
# Median : 0.16737   Median : 0.1565   Median : 0.15086   Median : 0.2082  
# Mean   : 0.17211   Mean   : 0.1338   Mean   : 0.15726   Mean   : 0.2183  
# 3rd Qu.: 0.42472   3rd Qu.: 0.2744   3rd Qu.: 0.30208   3rd Qu.: 0.3365  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# Tyrosine_Average   Valine_Average      X1_Methylhistidine X2_Hydroxybutyrate
# Min.   :-0.27447   Min.   :-0.556598   Min.   :-0.42457   Min.   :-0.35034  
# 1st Qu.: 0.01771   1st Qu.: 0.007498   1st Qu.:-0.14331   1st Qu.:-0.09362  
# Median : 0.12332   Median : 0.174492   Median :-0.02288   Median : 0.05999  
# Mean   : 0.15438   Mean   : 0.212389   Mean   :-0.03730   Mean   : 0.10027  
# 3rd Qu.: 0.26530   3rd Qu.: 0.454069   3rd Qu.: 0.05115   3rd Qu.: 0.27863  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# X2_Hydroxyisovaleric_Acid X3_Hydroxybutyrate      ADP                 AMP          
# Min.   :-0.24252          Min.   :-0.3587    Min.   :-0.572444   Min.   :-0.30993  
# 1st Qu.:-0.06316          1st Qu.: 0.0729    1st Qu.:-0.294440   1st Qu.:-0.10373  
# Median : 0.08973          Median : 0.2461    Median :-0.137938   Median :-0.01512  
# Mean   : 0.06665          Mean   : 0.2783    Mean   :-0.129910   Mean   :-0.01370  
# 3rd Qu.: 0.16626          3rd Qu.: 0.4264    3rd Qu.: 0.004719   3rd Qu.: 0.05376  
# Max.   : 1.00000          Max.   : 1.0000    Max.   : 1.000000   Max.   : 1.00000  
# ATP             Acetamide           Acetate          Acetoacetate     
# Min.   :-0.32877   Min.   :-0.60333   Min.   :-0.27396   Min.   :-0.33964  
# 1st Qu.:-0.18565   1st Qu.:-0.28267   1st Qu.: 0.04897   1st Qu.:-0.15159  
# Median :-0.11324   Median :-0.10693   Median : 0.17529   Median :-0.06400  
# Mean   :-0.09323   Mean   :-0.12018   Mean   : 0.19178   Mean   :-0.05388  
# 3rd Qu.:-0.03924   3rd Qu.: 0.02133   3rd Qu.: 0.30384   3rd Qu.: 0.00957  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Acetone           Adenosine           Alanine           Aspartate       
# Min.   :-0.36054   Min.   :-0.41978   Min.   :-0.18010   Min.   :-0.43109  
# 1st Qu.: 0.04478   1st Qu.:-0.05937   1st Qu.: 0.02725   1st Qu.:-0.10826  
# Median : 0.15414   Median : 0.09349   Median : 0.10438   Median : 0.06753  
# Mean   : 0.20467   Mean   : 0.08842   Mean   : 0.12078   Mean   : 0.10041  
# 3rd Qu.: 0.31899   3rd Qu.: 0.22204   3rd Qu.: 0.20379   3rd Qu.: 0.29799  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Creatine_Phosphate    Cytidine          Dimethyl_Sulfone      Ethanol         
# Min.   :-0.25800   Min.   :-0.3814700   Min.   :-0.38909   Min.   :-0.385297  
# 1st Qu.:-0.07229   1st Qu.:-0.1491086   1st Qu.:-0.05188   1st Qu.:-0.158086  
# Median : 0.03471   Median : 0.0006488   Median : 0.09112   Median :-0.106751  
# Mean   : 0.03514   Mean   : 0.0096985   Mean   : 0.11861   Mean   :-0.069940  
# 3rd Qu.: 0.11133   3rd Qu.: 0.1185031   3rd Qu.: 0.27929   3rd Qu.:-0.003025  
# Max.   : 1.00000   Max.   : 1.0000000   Max.   : 1.00000   Max.   : 1.000000  
# Ethanolamine         Formate             Glycerol         Guanidoacetate     
# Min.   :-0.54881   Min.   :-0.320546   Min.   :-0.177262   Min.   :-0.417646  
# 1st Qu.:-0.20337   1st Qu.:-0.125611   1st Qu.:-0.035192   1st Qu.:-0.132523  
# Median :-0.02285   Median :-0.021935   Median : 0.007469   Median : 0.001254  
# Mean   :-0.02507   Mean   :-0.004978   Mean   : 0.028765   Mean   : 0.007974  
# 3rd Qu.: 0.12969   3rd Qu.: 0.085000   3rd Qu.: 0.078516   3rd Qu.: 0.112288  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.000000   Max.   : 1.000000  
# Guanosine         Hypoxanthine           IMP              Inosine        
# Min.   :-0.41547   Min.   :-0.35127   Min.   :-0.44915   Min.   :-0.41978  
# 1st Qu.:-0.15554   1st Qu.:-0.11536   1st Qu.:-0.02809   1st Qu.:-0.18006  
# Median : 0.01527   Median : 0.04649   Median : 0.07865   Median :-0.02712  
# Mean   : 0.01267   Mean   : 0.05474   Mean   : 0.11916   Mean   : 0.01935  
# 3rd Qu.: 0.15869   3rd Qu.: 0.20506   3rd Qu.: 0.29506   3rd Qu.: 0.22539  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Isopropanol           Malate            Malonate           Mannose        
# Min.   :-0.46007   Min.   :-0.18871   Min.   :-0.23515   Min.   :-0.47985  
# 1st Qu.:-0.18718   1st Qu.:-0.02421   1st Qu.: 0.02455   1st Qu.:-0.08915  
# Median :-0.04526   Median : 0.10082   Median : 0.11273   Median : 0.05504  
# Mean   :-0.04878   Mean   : 0.10158   Mean   : 0.11326   Mean   : 0.06702  
# 3rd Qu.: 0.04659   3rd Qu.: 0.18343   3rd Qu.: 0.18880   3rd Qu.: 0.23173  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Methanol         N_N_Dimethylglycine  Nicotinurate      O_Acetylcarnitine  
# Min.   :-0.427228   Min.   :-0.39154    Min.   :-0.28084   Min.   :-0.252116  
# 1st Qu.:-0.117090   1st Qu.:-0.09306    1st Qu.: 0.06457   1st Qu.: 0.003107  
# Median : 0.002353   Median : 0.04730    Median : 0.18408   Median : 0.175496  
# Mean   :-0.003713   Mean   : 0.04301    Mean   : 0.18342   Mean   : 0.165500  
# 3rd Qu.: 0.097549   3rd Qu.: 0.18542    3rd Qu.: 0.30967   3rd Qu.: 0.306387  
# Max.   : 1.000000   Max.   : 1.00000    Max.   : 1.00000   Max.   : 1.000000  
# O_Phosphocholine   Propylene_Glycol      Taurine            Uridine        
# Min.   :-0.30298   Min.   :-0.37081   Min.   :-0.23627   Min.   :-0.34058  
# 1st Qu.:-0.10253   1st Qu.:-0.09003   1st Qu.:-0.05506   1st Qu.: 0.04448  
# Median : 0.07349   Median : 0.07519   Median : 0.02173   Median : 0.18056  
# Mean   : 0.08485   Mean   : 0.07334   Mean   : 0.05715   Mean   : 0.20100  
# 3rd Qu.: 0.24876   3rd Qu.: 0.22664   3rd Qu.: 0.10801   3rd Qu.: 0.35857  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# sn_Glycero_3_Phosphocholine  Beta.Alanine       
# Min.   :-0.380609           Min.   :-0.2724350  
# 1st Qu.:-0.178012           1st Qu.:-0.0709378  
# Median : 0.019231           Median : 0.0003351  
# Mean   : 0.004626           Mean   : 0.0123529  
# 3rd Qu.: 0.119182           3rd Qu.: 0.0796982  
# Max.   : 1.000000           Max.   : 1.0000000  

corrplot(DM_cor3, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
plot_correlation18no0 <- corrplot(DM_cor3, type = "upper", order = "hclust", 
                                    tl.col = "black", tl.srt = 60, 
                                    tl.cex = 1, cl.cex = 1)
plot_correlation18no0
dev.list()
# Open the graphics device
png("plot_correlation18no0")
print(plot_correlation18no0)
dev.off()


#### PCA ####

# PCA with FactoMineR package

DM_pca_FM3 <- PCA(DM_data_pca3, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(DM_pca_FM3)
# eigenvalue variance.percent cumulative.variance.percent
# Dim.1  37.1958078      23.24737987                    23.24738
# Dim.2  21.6851874      13.55324213                    36.80062
# Dim.3  16.2343513      10.14646955                    46.94709
# Dim.4  12.8399680       8.02497998                    54.97207
# Dim.5   7.6012923       4.75080767                    59.72288
# Dim.6   6.2372815       3.89830097                    63.62118
# Dim.7   4.4301577       2.76884859                    66.39003
# Dim.8   3.6537949       2.28362181                    68.67365
# Dim.9   3.4298974       2.14368589                    70.81734
# Dim.10  3.3273185       2.07957405                    72.89691
# Dim.11  3.0320851       1.89505318                    74.79196
# Dim.12  2.9757161       1.85982255                    76.65179
# Dim.13  2.4701479       1.54384246                    78.19563
# Dim.14  2.2498993       1.40618705                    79.60182
# Dim.15  2.0594337       1.28714604                    80.88896
# Dim.16  1.9252265       1.20326655                    82.09223
# Dim.17  1.9037942       1.18987135                    83.28210
# Dim.18  1.7987116       1.12419474                    84.40629
# Dim.19  1.7449686       1.09060540                    85.49690
# Dim.20  1.6455068       1.02844173                    86.52534
# Dim.21  1.4527269       0.90795432                    87.43330
# Dim.22  1.3838245       0.86489031                    88.29819
# Dim.23  1.3184385       0.82402409                    89.12221
# Dim.24  1.2140344       0.75877148                    89.88098
# Dim.25  1.2013572       0.75084825                    90.63183
# Dim.26  1.1576508       0.72353176                    91.35536
# Dim.27  1.0912954       0.68205961                    92.03742
# Dim.28  0.9920799       0.62004995                    92.65747
# Dim.29  0.9799172       0.61244827                    93.26992
# Dim.30  0.9345609       0.58410058                    93.85402
# Dim.31  0.8406819       0.52542618                    94.37945
# Dim.32  0.7553150       0.47207187                    94.85152
# Dim.33  0.7074422       0.44215137                    95.29367
# Dim.34  0.6982144       0.43638402                    95.73005
# Dim.35  0.6545505       0.40909405                    96.13915
# Dim.36  0.6037797       0.37736233                    96.51651
# Dim.37  0.5688654       0.35554086                    96.87205
# Dim.38  0.5138191       0.32113694                    97.19319
# Dim.39  0.4844305       0.30276907                    97.49596
# Dim.40  0.4519307       0.28245666                    97.77841
# Dim.41  0.4443792       0.27773698                    98.05615
# Dim.42  0.4274337       0.26714605                    98.32330
# Dim.43  0.3682660       0.23016625                    98.55346
# Dim.44  0.3538328       0.22114550                    98.77461
# Dim.45  0.3468749       0.21679684                    98.99141
# Dim.46  0.2853876       0.17836727                    99.16977
# Dim.47  0.2648099       0.16550620                    99.33528
# Dim.48  0.2395349       0.14970930                    99.48499
# Dim.49  0.2104227       0.13151418                    99.61650
# Dim.50  0.1858375       0.11614845                    99.73265
# Dim.51  0.1559617       0.09747609                    99.83013
# Dim.52  0.1420940       0.08880876                    99.91894
# Dim.53  0.1297033       0.08106458                   100.00000             

fviz_screeplot(DM_pca_FM3, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
screeplot18_allno0 <- fviz_screeplot(DM_pca_FM3, addlabels = TRUE, 
                                       ggtheme = theme_classic(),
                                       main = "",
                                       font.x = c(14, "bold"), font.y = c(14, "bold"),
                                       font.tickslab = 12,
                                       barfill = "#99d8c9", barcolor = "#66c2a4",
                                       font.submain = 16,
                                       ncp = Inf
)
screeplot18_allno0

# Open the graphics device
png("screeplot18_allno0")
print(screeplot18_allno0)
dev.off()

fviz_screeplot(DM_pca_FM3, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 10,
               ncp = 16
)

screeplot18_partialno0<- fviz_screeplot(DM_pca_FM3, addlabels = TRUE, 
                                          ggtheme = theme_classic(),
                                          main = "",
                                          font.x = c(14, "bold"), font.y = c(14, "bold"),
                                          font.tickslab = 12,
                                          barfill = "#99d8c9", barcolor = "#66c2a4",
                                          font.submain = 10,
                                          ncp = 16
)
screeplot18_partialno0


# Open the graphics device
png("screeplot18_partialno0")
print(screeplot18_partialno0)
dev.off()

# extract the results for variables
var30 <- get_pca_var(DM_pca_FM3)
var30
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
head(var30$coord)  # coordinates of variables
head(var30$contrib) # contributions of variables to the PCs
head(var30$cos2)  # Cos2: quality on the factore map
var30$coord
# Dim.1        Dim.2         Dim.3        Dim.4
# Serine                       0.331101645  0.201066568  0.8264140530 -0.135711731
# Putrescine                   0.525783164 -0.503749796  0.0966336006  0.188401639
# Trans_Hydroxyproline         0.293793941 -0.769157021 -0.0642484951 -0.104485685
# Glutamine                    0.436141231 -0.168315062  0.6161848573 -0.134244903
# Alpha_Aminoadipic_Acid       0.352440227 -0.327952546  0.2755399786  0.081273008
# Methionine_Sulfoxide         0.237586734 -0.149467636  0.7260839236 -0.156560440
# Acetyl_Ornithine             0.436391107 -0.272606788  0.1420870381  0.037307366
# Citrulline                  -0.183133649  0.213393607 -0.3424371873  0.244910164
# Asymmetric_Dimethylarginine  0.536118158 -0.567121207  0.3185357602 -0.110241674
# Total_Dimethylarginine       0.635152165 -0.446686196  0.3989426415 -0.241657845
# Tryptophan                   0.619057655 -0.535940893  0.2068025711 -0.222132290
# Ornithine                    0.155354913  0.249828786  0.7909366618 -0.074746088
# Lysine                       0.028818456  0.282906653  0.7522272787 -0.057572734
# Spermidine                   0.543629187 -0.541283246 -0.0812066979  0.058637861
# Spermine                     0.551929177 -0.443153143 -0.0827519324  0.086419094
# Sarcosine                    0.301613798  0.168603738  0.6643976299 -0.142766670
# Methylhistidine             -0.086505955  0.273114179  0.3086875965 -0.134210184
# Beta_Hydroxybutyric_Acid     0.751203178  0.139306623 -0.2274228210 -0.425316693
# Alpha_Ketoglutaric_Acid      0.550963676 -0.145621795  0.5375405635 -0.044685806
# Butyric_acid                 0.497263905 -0.344748186 -0.2884364288  0.048906489
# Propionic_Acid               0.047289698  0.067986955 -0.0067111738  0.437496510
# Fumaric_Acid                 0.483548535 -0.016001048  0.5027697617 -0.153012015
# Isobutyric_Acid              0.414266065 -0.677999990 -0.3168506574 -0.208270740
# Hippuric_Acid                0.088756356 -0.469591893 -0.1644599295  0.066964016
# Methylmalonic_Acid           0.183768318 -0.212920473  0.0385569724  0.091207675
# LYSOC14_0                    0.383885890  0.141089466 -0.0186889555  0.658916807
# LYSOC16_1                    0.320287066  0.168332361 -0.0936682499  0.518963099
# LYSOC16_0                    0.142835368  0.304288801 -0.1269311570  0.495717344
# LYSOC17_0                    0.302533292  0.324412987 -0.1447173613  0.677022300
# LYSOC18_2                    0.594973833  0.091766806  0.2189796075  0.454496721
# LYSOC18_1                    0.571424928  0.061823730  0.3229122658  0.385187233
# LYSOC18_0                    0.528226292  0.043802852  0.5312111183  0.267837389
# LYSOC20_4                    0.644741553 -0.464937288 -0.0394763369  0.194328087
# LYSOC20_3                    0.287599042  0.028083207  0.1948782934  0.521988734
# LYSOC24_0                    0.250237515  0.418988706 -0.0472767624  0.418319892
# LYSOC26_1                    0.464449283  0.167885002 -0.0040531972  0.531897478
# LYSOC26_0                    0.205965498  0.382263702 -0.0262683896  0.555880892
# LYSOC28_1                    0.617191021  0.030043906 -0.1442004883  0.536613238
# LYSOC28_0                    0.335806975  0.216938824 -0.0588770585  0.490576848
# X14_1SMOH                    0.636307224 -0.109480567 -0.4094967067  0.431441644
# X16_1SM                      0.714430141 -0.452064856 -0.3617842946  0.012976782
# X16_0SM                      0.692837293 -0.486912631 -0.3571097659  0.085385259
# X16_1SMOH                    0.407329309 -0.532051662 -0.4487215104  0.219140171
# X18_1SM                      0.681437821 -0.100301495  0.3500436239  0.382956105
# PC32_2AA                     0.396974416  0.403514392  0.5045757173  0.423672907
# X18_0SM                      0.472648674  0.380383113  0.5317205660  0.329192779
# X20_2SM                      0.590052415  0.172050792  0.4683580651  0.501680176
# PC36_0AE                     0.605358890 -0.509154820 -0.2586900127  0.292314559
# PC36_6AA                     0.604978159 -0.029124241  0.0956460723  0.692774389
# PC36_0AA                     0.405044156 -0.670802775 -0.2852756464  0.255439144
# X22_2SMOH                    0.729276881 -0.303108826 -0.4444364074  0.227870937
# X22_1SMOH                    0.810312120 -0.202334044 -0.2119580586  0.375124535
# PC38_6AA                     0.665167162  0.046488564  0.0335451280  0.603372213
# PC38_0AA                     0.733167287 -0.331554516 -0.2768631778  0.398112921
# PC40_6AE                     0.712292027  0.051621438 -0.0308538052  0.550552197
# X24_1SMOH                    0.577445030  0.250456313  0.0172756694  0.534881277
# PC40_6AA                     0.587039795  0.084640331  0.4905133851  0.274243152
# PC40_2AA                     0.652646899  0.039842871 -0.1264888540  0.623112043
# PC401AA                      0.667909262  0.218886077  0.1244461289  0.505896796
# C2                           0.502200777  0.273480675  0.0251837351 -0.089500686
# C3_1                         0.493878444 -0.289416504 -0.2057331665 -0.047922949
# C3                           0.583194377 -0.581914378 -0.0004277435 -0.037656116
# C4_1                         0.533830472 -0.406949420 -0.2071835719  0.324927194
# C4                           0.789259255 -0.295201465 -0.2404023960 -0.259497525
# C3OH                         0.405000268 -0.129446923 -0.0210149273  0.350570233
# C5_1                         0.462021742 -0.517473006 -0.0862741309 -0.037564677
# C5                           0.481631645 -0.624250737 -0.0641075736 -0.115930345
# C4OH                         0.570954622  0.024799724 -0.2588376185 -0.126582175
# C6_1                         0.382691626 -0.063755723 -0.6060676748  0.126036795
# C6                           0.668500495  0.335647085 -0.3331979832 -0.169102877
# C5OH                         0.476905223 -0.063893220  0.2961040575  0.207604153
# C5_1DC                       0.690769374 -0.190358251 -0.2525361205 -0.146795355
# C5DC                         0.581027159 -0.090451142 -0.1461019502 -0.091647951
# C8                           0.737809493  0.424148480 -0.0802903716 -0.242181712
# C5MDC                        0.454816535  0.384469678  0.0571064128 -0.240834565
# C9                           0.748752183  0.473499834 -0.0736863148 -0.228873550
# C7DC                         0.090137943  0.140809058  0.0317298518  0.138079217
# C10_2                        0.397126120  0.611200148 -0.0073637218 -0.155870384
# C10_1                        0.414281575 -0.084273068 -0.1873379824 -0.011954645
# C10                          0.668650032  0.531873323 -0.0844363027 -0.291536661
# C12_1                        0.777967609  0.249894298 -0.2834345992 -0.230765695
# C12                          0.696414606  0.547975642 -0.2389638013 -0.285283197
# C14_2                        0.748366092  0.454041342 -0.2829041453 -0.114342586
# C14_1                        0.827561541  0.429051781 -0.1863266576 -0.208253165
# C14                          0.747218085  0.478278971 -0.2772450065 -0.278899939
# C12DC                        0.589078074  0.681488605  0.0765824797 -0.153000372
# C14_2OH                      0.699304808  0.628571977 -0.1125643059 -0.237511635
# C14_1OH                      0.734904672  0.487682003 -0.3103550870 -0.263496087
# C16_2                        0.651093193  0.639576060 -0.0383554469 -0.280565688
# C16_1                        0.691433595  0.611248206 -0.1634631760 -0.274573616
# C16                          0.689593414  0.547807188 -0.3450497447 -0.201501094
# C16_2OH                      0.711302068  0.558568635 -0.3110942203 -0.197792385
# C16_1OH                      0.724647471  0.516155504 -0.3349398033 -0.148600916
# C16OH                        0.759409919  0.417461942 -0.2376230737 -0.276589245
# C18_2                        0.738643961  0.581261567 -0.1311552730 -0.245692425
# C18_1                        0.705807743  0.610957917 -0.2216947894 -0.215904881
# C18                          0.668002365  0.594375477 -0.3027985425 -0.121709833
# C18_1OH                      0.709058799  0.546547101 -0.2951340778 -0.247835514
# Arginine_Average             0.198886377  0.214009597  0.8320144178  0.001671720
# Betaine_Average              0.088723408  0.231543418  0.2535314667  0.177505657
# C0_Average                   0.565516798  0.279520772  0.1575340591 -0.082611889
# Choline_Average              0.253662730  0.377150645  0.4707220549  0.304904147
# Citrate_Average              0.678258204 -0.359105446 -0.0858739313 -0.265745425
# Creatine_Average             0.517637243 -0.351821195 -0.2564683152 -0.066972692
# Creatinine_Average           0.499170215 -0.463683652 -0.1906367863 -0.093460536
# Glucose_Average              0.137133831 -0.229771088  0.5485655085  0.438299608
# Glutamate_Average            0.432404895 -0.186859836  0.5986721304 -0.145829169
# Glycine_Average              0.470789200 -0.035858538  0.6068090944 -0.282031417
# Histidine_Average            0.493044889 -0.301145018  0.4818665624  0.116975516
# Isoleucine_Average           0.502701268 -0.567041781  0.2425642545 -0.197538414
# Lactate_Average              0.606538733 -0.428159744  0.2512353012 -0.010908469
# Leucine                      0.305573752 -0.336182941  0.5433414584 -0.196185877
# Methionine                   0.383044031 -0.224539948  0.7644746193 -0.217118729
# Phenylalanine_Average        0.437964692 -0.267300887  0.7144757041 -0.161224167
# Proline                      0.423652792 -0.733721196  0.1941898842 -0.095789439
# Pyruvate_Average             0.253886972  0.350852207  0.5816443506 -0.008304074
# Succinate_Average            0.281286133  0.260992663  0.5222000816 -0.086849547
# Threonine_Average            0.475107590 -0.157983780  0.4828347436 -0.083336956
# Tyrosine_Average             0.334275545 -0.263972028  0.4651865924 -0.262082086
# Valine_Average               0.534667814 -0.637932352  0.2199002843 -0.209645568
# X1_Methylhistidine          -0.081910679  0.188531050  0.0817430609 -0.114053904
# X2_Hydroxybutyrate           0.265964202 -0.483474692 -0.1092552445 -0.375773267
# X2_Hydroxyisovaleric_Acid    0.167889357 -0.076927682 -0.2416872446 -0.226773571
# X3_Hydroxybutyrate           0.736581298  0.193472526 -0.2350782650 -0.445329414
# ADP                         -0.371521693  0.482940145  0.0771149893  0.070186452
# AMP                         -0.089762111  0.176937871  0.0766729198 -0.107286059
# ATP                         -0.243048408  0.106961579 -0.1470813738 -0.123154552
# Acetamide                   -0.407185231  0.533150958  0.1249451526  0.044502577
# Acetate                      0.426569611  0.124599355 -0.0479151887 -0.449567527
# Acetoacetate                -0.085435669 -0.034814389 -0.2220597989 -0.233405506
# Acetone                      0.508337380  0.209199399 -0.2379879529 -0.438822712
# Adenosine                    0.133447405  0.323663452  0.2433355073 -0.309166294
# Alanine                      0.204468165  0.002321407  0.3387533103 -0.142997563
# Aspartate                    0.203271269 -0.574986116  0.0184068425 -0.308791051
# Creatine_Phosphate           0.001337456  0.216000417  0.2188130108 -0.005355832
# Cytidine                    -0.064355757 -0.305519274  0.1591586089 -0.243841230
# Dimethyl_Sulfone             0.268701148 -0.464067233 -0.1386198855 -0.319770291
# Ethanol                     -0.202542212 -0.096299382 -0.1389383819  0.015721023
# Ethanolamine                -0.201303018  0.547094275  0.1910044786 -0.092884799
# Formate                     -0.058382183 -0.182670845 -0.0315717715  0.112364370
# Glycerol                     0.026690137  0.059226905  0.0981745199  0.041770854
# Guanidoacetate              -0.116341605  0.262725790  0.2567787322 -0.148122508
# Guanosine                   -0.071860591  0.351732441  0.3070842885 -0.187891512
# Hypoxanthine                 0.028810052  0.459316878  0.1552846276  0.018526873
# IMP                          0.270586350 -0.461987117  0.0182361080 -0.266866971
# Inosine                      0.041219841 -0.561774761 -0.1482870796 -0.137058903
# Isopropanol                 -0.200201985  0.290797993  0.1028486350  0.147176870
# Malate                       0.182657904  0.115227381  0.1410455178 -0.375750351
# Malonate                     0.210984824 -0.105081773 -0.0087733051 -0.319951571
# Mannose                      0.077341208 -0.198410533  0.5064850999 -0.263366822
# Methanol                    -0.032839222 -0.274531180  0.0264927610  0.148927608
# N_N_Dimethylglycine          0.030494813  0.408121787  0.1748090394 -0.180895003
# Nicotinurate                 0.369874310 -0.172289087  0.2759914608 -0.404657158
# O_Acetylcarnitine            0.442512373  0.215039987 -0.0457641490 -0.110028517
# O_Phosphocholine             0.094925009 -0.428868971  0.3197227504 -0.208311960
# Propylene_Glycol             0.225969821 -0.288492952 -0.3044038750 -0.456514570
# Taurine                      0.028049233 -0.047652014  0.0880794106 -0.300814596
# Uridine                      0.403645747 -0.330633117  0.3767410486 -0.386804740
# sn_Glycero_3_Phosphocholine -0.097253938  0.391068470  0.3380706520  0.111628238
# Beta.Alanine                -0.040474472  0.106354515  0.1112527783  0.107177906
# Dim.5
# Serine                      -0.075688528
# Putrescine                  -0.214584521
# Trans_Hydroxyproline         0.125828791
# Glutamine                   -0.010491092
# Alpha_Aminoadipic_Acid       0.047965373
# Methionine_Sulfoxide        -0.063952409
# Acetyl_Ornithine             0.218038533
# Citrulline                   0.680053989
# Asymmetric_Dimethylarginine -0.030514734
# Total_Dimethylarginine      -0.080002320
# Tryptophan                   0.205584766
# Ornithine                   -0.210125322
# Lysine                      -0.264095407
# Spermidine                  -0.083287744
# Spermine                    -0.105425029
# Sarcosine                    0.088081558
# Methylhistidine             -0.167899047
# Beta_Hydroxybutyric_Acid     0.024670299
# Alpha_Ketoglutaric_Acid      0.217098391
# Butyric_acid                -0.250017059
# Propionic_Acid              -0.474997431
# Fumaric_Acid                 0.384112000
# Isobutyric_Acid              0.039008550
# Hippuric_Acid               -0.034189789
# Methylmalonic_Acid          -0.113469515
# LYSOC14_0                    0.365671653
# LYSOC16_1                    0.576167739
# LYSOC16_0                    0.540529953
# LYSOC17_0                    0.388785515
# LYSOC18_2                    0.158547925
# LYSOC18_1                    0.172215293
# LYSOC18_0                   -0.005691351
# LYSOC20_4                    0.159461739
# LYSOC20_3                   -0.234818622
# LYSOC24_0                   -0.133489062
# LYSOC26_1                   -0.218417764
# LYSOC26_0                   -0.100813219
# LYSOC28_1                   -0.132540035
# LYSOC28_0                   -0.266396007
# X14_1SMOH                   -0.006885507
# X16_1SM                      0.066876269
# X16_0SM                      0.067716644
# X16_1SMOH                    0.249550375
# X18_1SM                      0.004985942
# PC32_2AA                    -0.051580735
# X18_0SM                     -0.046178072
# X20_2SM                     -0.077357159
# PC36_0AE                     0.305120463
# PC36_6AA                     0.120139379
# PC36_0AA                     0.295053507
# X22_2SMOH                    0.080802448
# X22_1SMOH                    0.092843463
# PC38_6AA                     0.069331289
# PC38_0AA                     0.154062817
# PC40_6AE                    -0.013752657
# X24_1SMOH                    0.092725391
# PC40_6AA                    -0.186892407
# PC40_2AA                    -0.022751936
# PC401AA                     -0.081667372
# C2                          -0.453355985
# C3_1                        -0.082935109
# C3                          -0.009461014
# C4_1                        -0.130689535
# C4                          -0.015584410
# C3OH                         0.065165257
# C5_1                        -0.280108771
# C5                          -0.004669109
# C4OH                        -0.055927426
# C6_1                         0.123992658
# C6                          -0.105448014
# C5OH                        -0.193181673
# C5_1DC                      -0.020011010
# C5DC                        -0.369462254
# C8                          -0.109185129
# C5MDC                       -0.014510882
# C9                          -0.015513063
# C7DC                         0.207072866
# C10_2                       -0.167668895
# C10_1                       -0.068716813
# C10                         -0.079459882
# C12_1                       -0.104052323
# C12                         -0.034732471
# C14_2                        0.089216008
# C14_1                       -0.035162115
# C14                         -0.034641026
# C12DC                       -0.011764100
# C14_2OH                     -0.027645674
# C14_1OH                     -0.004037208
# C16_2                       -0.039640910
# C16_1                        0.005360863
# C16                         -0.020779716
# C16_2OH                     -0.004972031
# C16_1OH                      0.011185637
# C16OH                       -0.068493446
# C18_2                       -0.039854020
# C18_1                       -0.017240351
# C18                          0.011260436
# C18_1OH                     -0.021391521
# Arginine_Average            -0.247786486
# Betaine_Average              0.546287207
# C0_Average                  -0.282966142
# Choline_Average              0.147515633
# Citrate_Average              0.159294849
# Creatine_Average            -0.087693389
# Creatinine_Average          -0.207492794
# Glucose_Average             -0.090478948
# Glutamate_Average            0.184857416
# Glycine_Average             -0.031859474
# Histidine_Average           -0.145460836
# Isoleucine_Average           0.093732972
# Lactate_Average              0.041133988
# Leucine                     -0.025516683
# Methionine                  -0.006728918
# Phenylalanine_Average        0.051138300
# Proline                      0.126267399
# Pyruvate_Average            -0.267530258
# Succinate_Average            0.493887569
# Threonine_Average           -0.346197292
# Tyrosine_Average             0.185878258
# Valine_Average               0.147774442
# X1_Methylhistidine          -0.247521741
# X2_Hydroxybutyrate           0.283113653
# X2_Hydroxyisovaleric_Acid    0.203307563
# X3_Hydroxybutyrate          -0.012239810
# ADP                          0.312235163
# AMP                          0.041448426
# ATP                          0.211980832
# Acetamide                    0.182369136
# Acetate                      0.221411120
# Acetoacetate                -0.335613132
# Acetone                      0.217618962
# Adenosine                    0.396220709
# Alanine                      0.148151740
# Aspartate                    0.099884935
# Creatine_Phosphate           0.123918527
# Cytidine                     0.032449277
# Dimethyl_Sulfone             0.225406804
# Ethanol                     -0.315830970
# Ethanolamine                 0.490125858
# Formate                     -0.321793668
# Glycerol                     0.043635158
# Guanidoacetate               0.264604795
# Guanosine                    0.567920799
# Hypoxanthine                 0.456059581
# IMP                          0.071613229
# Inosine                      0.055399884
# Isopropanol                 -0.026920457
# Malate                       0.069669579
# Malonate                     0.197854231
# Mannose                     -0.332875552
# Methanol                    -0.514646263
# N_N_Dimethylglycine         -0.049815246
# Nicotinurate                 0.398445074
# O_Acetylcarnitine           -0.499073386
# O_Phosphocholine             0.279364162
# Propylene_Glycol            -0.056167110
# Taurine                      0.263712047
# Uridine                      0.170813663
# sn_Glycero_3_Phosphocholine  0.238833434
# Beta.Alanine                -0.078363498
 



var30$contrib
# Dim.1        Dim.2        Dim.3        Dim.4
# Serine                      2.947329e-01 1.864303e-01 4.206883e+00 1.434402e-01
# Putrescine                  7.432234e-01 1.170217e+00 5.752033e-02 2.764429e-01
# Trans_Hydroxyproline        2.320554e-01 2.728141e+00 2.542676e-02 8.502559e-02
# Glutamine                   5.113995e-01 1.306420e-01 2.338768e+00 1.403562e-01
# Alpha_Aminoadipic_Acid      3.339465e-01 4.959739e-01 4.676644e-01 5.144329e-02
# Methionine_Sulfoxide        1.517576e-01 1.030223e-01 3.247422e+00 1.908974e-01
# Acetyl_Ornithine            5.119857e-01 3.426969e-01 1.243581e-01 1.083990e-02
# Citrulline                  9.016590e-02 2.099905e-01 7.223154e-01 4.671428e-01
# Asymmetric_Dimethylarginine 7.727287e-01 1.483162e+00 6.250021e-01 9.465153e-02
# Total_Dimethylarginine      1.084580e+00 9.201145e-01 9.803609e-01 4.548182e-01
# Tryptophan                  1.030311e+00 1.324557e+00 2.634371e-01 3.842903e-01
# Ornithine                   6.488675e-02 2.878205e-01 3.853439e+00 4.351240e-02
# Lysine                      2.232788e-03 3.690822e-01 3.485485e+00 2.581486e-02
# Spermidine                  7.945323e-01 1.351095e+00 4.062083e-02 2.677887e-02
# Spermine                    8.189789e-01 9.056168e-01 4.218143e-02 5.816416e-02
# Sarcosine                   2.445729e-01 1.310905e-01 2.719075e+00 1.587412e-01
# Methylhistidine             2.011861e-02 3.439738e-01 5.869531e-01 1.402836e-01
# Beta_Hydroxybutyric_Acid    1.517123e+00 8.949120e-02 3.185907e-01 1.408838e+00
# Alpha_Ketoglutaric_Acid     8.161161e-01 9.778890e-02 1.779867e+00 1.555161e-02
# Butyric_acid                6.647830e-01 5.480760e-01 5.124663e-01 1.862812e-02
# Propionic_Acid              6.012279e-03 2.131513e-02 2.774355e-04 1.490683e+00
# Fumaric_Acid                6.286170e-01 1.180684e-03 1.557053e+00 1.823422e-01
# Isobutyric_Acid             4.613863e-01 2.119806e+00 6.184068e-01 3.378256e-01
# Hippuric_Acid               2.117897e-02 1.016899e+00 1.666039e-01 3.492360e-02
# Methylmalonic_Acid          9.079194e-02 2.090603e-01 9.157373e-03 6.478864e-02
# LYSOC14_0                   3.961962e-01 9.179647e-02 2.151469e-03 3.381405e+00
# LYSOC16_1                   2.757940e-01 1.306688e-01 5.404430e-02 2.097534e+00
# LYSOC16_0                   5.485011e-02 4.269812e-01 9.924338e-02 1.913834e+00
# LYSOC17_0                   2.460664e-01 4.853257e-01 1.290049e-01 3.569785e+00
# LYSOC18_2                   9.517037e-01 3.883364e-02 2.953741e-01 1.608783e+00
# LYSOC18_1                   8.778582e-01 1.762573e-02 6.422944e-01 1.155526e+00
# LYSOC18_0                   7.501464e-01 8.847928e-03 1.738199e+00 5.586997e-01
# LYSOC20_4                   1.117577e+00 9.968403e-01 9.599282e-03 2.941083e-01
# LYSOC20_3                   2.223724e-01 3.636891e-03 2.339333e-01 2.122063e+00
# LYSOC24_0                   1.683491e-01 8.095459e-01 1.376767e-02 1.362866e+00
# LYSOC26_1                   5.799394e-01 1.299752e-01 1.011953e-04 2.203393e+00
# LYSOC26_0                   1.140499e-01 6.738496e-01 4.250421e-03 2.406576e+00
# LYSOC28_1                   1.024107e+00 4.162456e-03 1.280851e-01 2.242636e+00
# LYSOC28_0                   3.031694e-01 2.170258e-01 2.135292e-02 1.874348e+00
# X14_1SMOH                   1.088528e+00 5.527273e-02 1.032918e+00 1.449707e+00
# X16_1SM                     1.372226e+00 9.424066e-01 8.062403e-01 1.311505e-03
# X16_0SM                     1.290531e+00 1.093299e+00 7.855404e-01 5.678085e-02
# X16_1SMOH                   4.460642e-01 1.305402e+00 1.240277e+00 3.740073e-01
# X18_1SM                     1.248414e+00 4.639291e-02 7.547609e-01 1.142179e+00
# PC32_2AA                    4.236732e-01 7.508529e-01 1.568259e+00 1.397969e+00
# X18_0SM                     6.005966e-01 6.672357e-01 1.741534e+00 8.439888e-01
# X20_2SM                     9.360244e-01 1.365055e-01 1.351204e+00 1.960153e+00
# PC36_0AE                    9.852169e-01 1.195464e+00 4.122156e-01 6.654830e-01
# PC36_6AA                    9.839780e-01 3.911524e-03 5.635070e-02 3.737831e+00
# PC36_0AA                    4.410733e-01 2.075040e+00 5.012963e-01 5.081723e-01
# X22_2SMOH                   1.429851e+00 4.236761e-01 1.216702e+00 4.044026e-01
# X22_1SMOH                   1.765268e+00 1.887882e-01 2.767355e-01 1.095941e+00
# PC38_6AA                    1.189509e+00 9.966188e-03 6.931448e-03 2.835350e+00
# PC38_0AA                    1.445147e+00 5.069285e-01 4.721668e-01 1.234379e+00
# PC40_6AE                    1.364024e+00 1.228845e-02 5.863846e-03 2.360658e+00
# X24_1SMOH                   8.964525e-01 2.892683e-01 1.838378e-03 2.228183e+00
# PC40_6AA                    9.264908e-01 3.303631e-02 1.482063e+00 5.857437e-01
# PC40_2AA                    1.145150e+00 7.320455e-03 9.855294e-02 3.023906e+00
# PC401AA                     1.199336e+00 2.209394e-01 9.539549e-02 1.993241e+00
# C2                          6.780485e-01 3.448975e-01 3.906658e-03 6.238624e-02
# C3_1                        6.557619e-01 3.862633e-01 2.607196e-01 1.788641e-02
# C3                          9.143925e-01 1.561547e+00 1.127021e-06 1.104351e-02
# C4_1                        7.661481e-01 7.636910e-01 2.644087e-01 8.222581e-01
# C4                          1.674732e+00 4.018591e-01 3.559940e-01 5.244481e-01
# C3OH                        4.409777e-01 7.727167e-02 2.720325e-03 9.571635e-01
# C5_1                        5.738929e-01 1.234844e+00 4.584862e-02 1.098994e-02
# C5                          6.236430e-01 1.797028e+00 2.531534e-02 1.046719e-01
# C4OH                        8.764138e-01 2.836159e-03 4.126861e-01 1.247904e-01
# C6_1                        3.937349e-01 1.874456e-02 2.262598e+00 1.237174e-01
# C6                          1.201460e+00 5.195204e-01 6.838641e-01 2.227091e-01
# C5OH                        6.114630e-01 1.882549e-02 5.400746e-01 3.356666e-01
# C5_1DC                      1.282839e+00 1.671015e-01 3.928367e-01 1.678266e-01
# C5DC                        9.076092e-01 3.772810e-02 1.314853e-01 6.541564e-02
# C8                          1.463506e+00 8.296075e-01 3.970928e-02 4.567923e-01
# C5MDC                       5.561328e-01 6.816493e-01 2.008791e-02 4.517246e-01
# C9                          1.507239e+00 1.033895e+00 3.344558e-02 4.079691e-01
# C7DC                        2.184345e-02 9.143195e-02 6.201563e-03 1.484885e-01
# C10_2                       4.239971e-01 1.722676e+00 3.340103e-04 1.892184e-01
# C10_1                       4.614209e-01 3.275024e-02 2.161806e-01 1.113037e-03
# C10                         1.201998e+00 1.304527e+00 4.391607e-02 6.619458e-01
# C12_1                       1.627155e+00 2.879715e-01 4.948468e-01 4.147425e-01
# C12                         1.303892e+00 1.384712e+00 3.517461e-01 6.338528e-01
# C14_2                       1.505685e+00 9.506652e-01 4.929963e-01 1.018245e-01
# C14_1                       1.841224e+00 8.488994e-01 2.138529e-01 3.377686e-01
# C14                         1.501069e+00 1.054871e+00 4.734701e-01 6.058051e-01
# C12DC                       9.329357e-01 2.141677e+00 3.612634e-02 1.823144e-01
# C14_2OH                     1.314737e+00 1.821994e+00 7.804884e-02 4.393452e-01
# C14_1OH                     1.452005e+00 1.096757e+00 5.933115e-01 5.407349e-01
# C16_2                       1.139705e+00 1.886345e+00 9.061898e-03 6.130631e-01
# C16_1                       1.285307e+00 1.722947e+00 1.645906e-01 5.871562e-01
# C16                         1.278475e+00 1.383860e+00 7.333790e-01 3.162211e-01
# C16_2OH                     1.360236e+00 1.438765e+00 5.961409e-01 3.046879e-01
# C16_1OH                     1.411756e+00 1.228564e+00 6.910327e-01 1.719804e-01
# C16OH                       1.550453e+00 8.036568e-01 3.478102e-01 5.958084e-01
# C18_2                       1.466818e+00 1.558045e+00 1.059587e-01 4.701318e-01
# C18_1                       1.339303e+00 1.721311e+00 3.027443e-01 3.630454e-01
# C18                         1.199671e+00 1.629141e+00 5.647713e-01 1.153685e-01
# C18_1OH                     1.351669e+00 1.377501e+00 5.365421e-01 4.783691e-01
# Arginine_Average            1.063448e-01 2.112046e-01 4.264094e+00 2.176521e-05
# Betaine_Average             2.116325e-02 2.472303e-01 3.959395e-01 2.453920e-01
# C0_Average                  8.597992e-01 3.603006e-01 1.528671e-01 5.315219e-02
# Choline_Average             1.729893e-01 6.559436e-01 1.364879e+00 7.240403e-01
# Citrate_Average             1.236790e+00 5.946765e-01 4.542425e-02 5.500063e-01
# Creatine_Average            7.203724e-01 5.707959e-01 4.051655e-01 3.493265e-02
# Creatinine_Average          6.698897e-01 9.914719e-01 2.238610e-01 6.802877e-02
# Glucose_Average             5.055862e-02 2.434600e-01 1.853626e+00 1.496161e+00
# Glutamate_Average           5.026749e-01 1.610159e-01 2.207716e+00 1.656246e-01
# Glycine_Average             5.958802e-01 5.929553e-03 2.268137e+00 6.194853e-01
# Histidine_Average           6.535502e-01 4.182040e-01 1.430272e+00 1.065678e-01
# Isoleucine_Average          6.794007e-01 1.482747e+00 3.624254e-01 3.039059e-01
# Lactate_Average             9.890610e-01 8.453732e-01 3.888001e-01 9.267523e-04
# Leucine                     2.510372e-01 5.211805e-01 1.818489e+00 2.997585e-01
# Methionine                  3.944604e-01 2.325006e-01 3.599906e+00 3.671391e-01
# Phenylalanine_Average       5.156847e-01 3.294865e-01 3.144416e+00 2.024400e-01
# Proline                     4.825320e-01 2.482555e+00 2.322834e-01 7.146137e-02
# Pyruvate_Average            1.732953e-01 5.676560e-01 2.083915e+00 5.370547e-04
# Succinate_Average           2.127172e-01 3.141184e-01 1.679728e+00 5.874504e-02
# Threonine_Average           6.068620e-01 1.150964e-01 1.436025e+00 5.408930e-02
# Tyrosine_Average            3.004106e-01 3.213310e-01 1.332967e+00 5.349470e-01
# Valine_Average              7.685535e-01 1.876662e+00 2.978631e-01 3.423004e-01
# X1_Methylhistidine          1.803794e-02 1.639089e-01 4.115919e-02 1.013109e-01
# X2_Hydroxybutyrate          1.901745e-01 1.077914e+00 7.352747e-02 1.099734e+00
# X2_Hydroxyisovaleric_Acid   7.577960e-02 2.728991e-02 3.598094e-01 4.005170e-01
# X3_Hydroxybutyrate          1.458638e+00 1.726138e-01 3.404004e-01 1.544539e+00
# ADP                         3.710858e-01 1.075532e+00 3.663048e-02 3.836566e-02
# AMP                         2.166168e-02 1.443705e-01 3.621171e-02 8.964429e-02
# ATP                         1.588150e-01 5.275850e-02 1.332540e-01 1.181237e-01
# Acetamide                   4.457487e-01 1.310802e+00 9.616209e-02 1.542433e-02
# Acetate                     4.891993e-01 7.159265e-02 1.414202e-02 1.574077e+00
# Acetoacetate                1.962386e-02 5.589261e-03 3.037421e-01 4.242856e-01
# Acetone                     6.947205e-01 2.018170e-01 3.488791e-01 1.499734e+00
# Adenosine                   4.787693e-02 4.830857e-01 3.647338e-01 7.444240e-01
# Alanine                     1.123977e-01 2.485074e-05 7.068580e-01 1.592551e-01
# Aspartate                   1.110857e-01 1.524585e+00 2.087006e-03 7.426180e-01
# Creatine_Phosphate          4.809110e-06 2.151523e-01 2.949248e-01 2.234035e-04
# Cytidine                    1.113476e-02 4.304414e-01 1.560362e-01 4.630739e-01
# Dimethyl_Sulfone            1.941087e-01 9.931129e-01 1.183630e-01 7.963652e-01
# Ethanol                     1.102902e-01 4.276454e-02 1.189076e-01 1.924854e-03
# Ethanolamine                1.089448e-01 1.380261e+00 2.247254e-01 6.719320e-02
# Formate                     9.163611e-03 1.538776e-01 6.139924e-03 9.833165e-02
# Glycerol                    1.915171e-03 1.617614e-02 5.936940e-02 1.358885e-02
# Guanidoacetate              3.638950e-02 3.183041e-01 4.061469e-01 1.708749e-01
# Guanosine                   1.388314e-02 5.705079e-01 5.808717e-01 2.749479e-01
# Hypoxanthine                2.231486e-03 9.728853e-01 1.485327e-01 2.673255e-03
# IMP                         1.968420e-01 9.842299e-01 2.048469e-03 5.546586e-01
# Inosine                     4.567921e-03 1.455329e+00 1.354477e-01 1.463021e-01
# Isopropanol                 1.077563e-01 3.899596e-01 6.515716e-02 1.687000e-01
# Malate                      8.969804e-02 6.122774e-02 1.225416e-01 1.099600e+00
# Malonate                    1.196764e-01 5.092038e-02 4.741235e-04 7.972684e-01
# Mannose                     1.608155e-02 1.815375e-01 1.580150e+00 5.402045e-01
# Methanol                    2.899290e-03 3.475523e-01 4.323341e-03 1.727374e-01
# N_N_Dimethylglycine         2.500103e-03 7.680976e-01 1.882317e-01 2.548527e-01
# Nicotinurate                3.678022e-01 1.368839e-01 4.691982e-01 1.275295e+00
# O_Acetylcarnitine           5.264497e-01 2.132432e-01 1.290078e-02 9.428586e-02
# O_Phosphocholine            2.422520e-02 8.481762e-01 6.296688e-01 3.379594e-01
# Propylene_Glycol            1.372799e-01 3.838020e-01 5.707756e-01 1.623100e+00
# Taurine                     2.115183e-03 1.047127e-02 4.778745e-02 7.047480e-01
# Uridine                     4.380329e-01 5.041149e-01 8.742808e-01 1.165251e+00
# sn_Glycero_3_Phosphocholine 2.542848e-02 7.052489e-01 7.040119e-01 9.704747e-02
# Beta.Alanine                4.404214e-03 5.216133e-02 7.624069e-02 8.946365e-02
# Dim.5
# Serine                      0.0753655180
# Putrescine                  0.6057722140
# Trans_Hydroxyproline        0.2082920135
# Glutamine                   0.0014479514
# Alpha_Aminoadipic_Acid      0.0302669190
# Methionine_Sulfoxide        0.0538054645
# Acetyl_Ornithine            0.6254305235
# Citrulline                  6.0841421662
# Asymmetric_Dimethylarginine 0.0122498776
# Total_Dimethylarginine      0.0842010928
# Tryptophan                  0.5560251394
# Ornithine                   0.5808571652
# Lysine                      0.9175595597
# Spermidine                  0.0912588023
# Spermine                    0.1462177257
# Sarcosine                   0.1020663400
# Methylhistidine             0.3708591772
# Beta_Hydroxybutyric_Acid    0.0080068448
# Alpha_Ketoglutaric_Acid     0.6200486693
# Butyric_acid                0.8223408268
# Propionic_Acid              2.9682131905
# Fumaric_Acid                1.9410124369
# Isobutyric_Acid             0.0200185298
# Hippuric_Acid               0.0153781964
# Methylmalonic_Acid          0.1693834468
# LYSOC14_0                   1.7591187492
# LYSOC16_1                   4.3672740350
# LYSOC16_0                   3.8437231422
# LYSOC17_0                   1.9885326276
# LYSOC18_2                   0.3306996198
# LYSOC18_1                   0.3901719079
# LYSOC18_0                   0.0004261311
# LYSOC20_4                   0.3345226727
# LYSOC20_3                   0.7254001493
# LYSOC24_0                   0.2344250037
# LYSOC26_1                   0.6276080146
# LYSOC26_0                   0.1337049654
# LYSOC28_1                   0.2311036115
# LYSOC28_0                   0.9336153618
# X14_1SMOH                   0.0006237125
# X16_1SM                     0.0588378292
# X16_0SM                     0.0603258455
# X16_1SMOH                   0.8192737130
# X18_1SM                     0.0003270446
# PC32_2AA                    0.0350015775
# X18_0SM                     0.0280533127
# X20_2SM                     0.0787251672
# PC36_0AE                    1.2247719663
# PC36_6AA                    0.1898817964
# PC36_0AA                    1.1452864695
# X22_2SMOH                   0.0858937572
# X22_1SMOH                   0.1134005686
# PC38_6AA                    0.0632369797
# PC38_0AA                    0.3122541646
# PC40_6AE                    0.0024882028
# X24_1SMOH                   0.1131123220
# PC40_6AA                    0.4595109686
# PC40_2AA                    0.0068100339
# PC401AA                     0.0877424441
# C2                          2.7039040419
# C3_1                        0.0904876714
# C3                          0.0011775733
# C4_1                        0.2246954072
# C4                          0.0031951649
# C3OH                        0.0558656426
# C5_1                        1.0322050652
# C5                          0.0002868010
# C4OH                        0.0411492791
# C6_1                        0.2022574400
# C6                          0.1462814905
# C5OH                        0.4909580840
# C5_1DC                      0.0052680587
# C5DC                        1.7957782973
# C8                          0.1568337595
# C5MDC                       0.0027701303
# C9                          0.0031659766
# C7DC                        0.5641037092
# C10_2                       0.3698431439
# C10_1                       0.0621210215
# C10                         0.0830631507
# C12_1                       0.1424348067
# C12                         0.0158702561
# C14_2                       0.1047124067
# C14_1                       0.0162653176
# C14                         0.0157867986
# C12DC                       0.0018206646
# C14_2OH                     0.0100546491
# C14_1OH                     0.0002144247
# C16_2                       0.0206728240
# C16_1                       0.0003780785
# C16                         0.0056805685
# C16_2OH                     0.0003252223
# C16_1OH                     0.0016460159
# C16OH                       0.0617178238
# C18_2                       0.0208956962
# C18_1                       0.0039102522
# C18                         0.0016681034
# C18_1OH                     0.0060199917
# Arginine_Average            0.8077329533
# Betaine_Average             3.9260391731
# C0_Average                  1.0533713858
# Choline_Average             0.2862784537
# Citrate_Average             0.3338228289
# Creatine_Average            0.1011687249
# Creatinine_Average          0.5663939518
# Glucose_Average             0.1076980035
# Glutamate_Average           0.4495586167
# Glycine_Average             0.0133533355
# Histidine_Average           0.2783586542
# Isoleucine_Average          0.1155839002
# Lactate_Average             0.0222594380
# Leucine                     0.0085656632
# Methionine                  0.0005956663
# Phenylalanine_Average       0.0344036998
# Proline                     0.2097466535
# Pyruvate_Average            0.9415825119
# Succinate_Average           3.2089929176
# Threonine_Average           1.5767393308
# Tyrosine_Average            0.4545375390
# Valine_Average              0.2872838598
# X1_Methylhistidine          0.8060078509
# X2_Hydroxybutyrate          1.0544699231
# X2_Hydroxyisovaleric_Acid   0.5437755034
# X3_Hydroxybutyrate          0.0019708879
# ADP                         1.2825555614
# AMP                         0.0226010523
# ATP                         0.5911609710
# Acetamide                   0.4375374671
# Acetate                     0.6449282858
# Acetoacetate                1.4818029700
# Acetone                     0.6230258058
# Adenosine                   2.0653179050
# Alanine                     0.2887527175
# Aspartate                   0.1312540023
# Creatine_Phosphate          0.2020156667
# Cytidine                    0.0138523233
# Dimethyl_Sulfone            0.6684156526
# Ethanol                     1.3122663656
# Ethanolamine                3.1602962733
# Formate                     1.3622836873
# Glycerol                    0.0250487278
# Guanidoacetate              0.9211025573
# Guanosine                   4.2431473772
# Hypoxanthine                2.7362497581
# IMP                         0.0674681934
# Inosine                     0.0403766489
# Isopropanol                 0.0095340498
# Malate                      0.0638555923
# Malonate                    0.5149952850
# Mannose                     1.4577275710
# Methanol                    3.4844177366
# N_N_Dimethylglycine         0.0326465370
# Nicotinurate                2.0885722005
# O_Acetylcarnitine           3.2767355267
# O_Phosphocholine            1.0267245635
# Propylene_Glycol            0.0415027363
# Taurine                     0.9148976395
# Uridine                     0.3838466723
# sn_Glycero_3_Phosphocholine 0.7504172624
# Beta.Alanine                0.0807867618
var30$cos2
# Dim.1        Dim.2        Dim.3        Dim.4
# Serine                      1.096283e-01 4.042776e-02 6.829602e-01 1.841767e-02
# Putrescine                  2.764479e-01 2.537639e-01 9.338053e-03 3.549518e-02
# Trans_Hydroxyproline        8.631488e-02 5.916025e-01 4.127869e-03 1.091726e-02
# Glutamine                   1.902192e-01 2.832996e-02 3.796838e-01 1.802169e-02
# Alpha_Aminoadipic_Acid      1.242141e-01 1.075529e-01 7.592228e-02 6.605302e-03
# Methionine_Sulfoxide        5.644746e-02 2.234057e-02 5.271979e-01 2.451117e-02
# Acetyl_Ornithine            1.904372e-01 7.431446e-02 2.018873e-02 1.391840e-03
# Citrulline                  3.353793e-02 4.553683e-02 1.172632e-01 5.998099e-02
# Asymmetric_Dimethylarginine 2.874227e-01 3.216265e-01 1.014650e-01 1.215323e-02
# Total_Dimethylarginine      4.034183e-01 1.995286e-01 1.591552e-01 5.839851e-02
# Tryptophan                  3.832324e-01 2.872326e-01 4.276730e-02 4.934275e-02
# Ornithine                   2.413515e-02 6.241442e-02 6.255808e-01 5.586978e-03
# Lysine                      8.305034e-04 8.003617e-02 5.658459e-01 3.314620e-03
# Spermidine                  2.955327e-01 2.929876e-01 6.594528e-03 3.438399e-03
# Spermine                    3.046258e-01 1.963847e-01 6.847882e-03 7.468260e-03
# Sarcosine                   9.097088e-02 2.842722e-02 4.414242e-01 2.038232e-02
# Methylhistidine             7.483280e-03 7.459135e-02 9.528803e-02 1.801237e-02
# Beta_Hydroxybutyric_Acid    5.643062e-01 1.940634e-02 5.172114e-02 1.808943e-01
# Alpha_Ketoglutaric_Acid     3.035610e-01 2.120571e-02 2.889499e-01 1.996821e-03
# Butyric_acid                2.472714e-01 1.188513e-01 8.319557e-02 2.391845e-03
# Propionic_Acid              2.236316e-03 4.622226e-03 4.503985e-05 1.914032e-01
# Fumaric_Acid                2.338192e-01 2.560335e-04 2.527774e-01 2.341268e-02
# Isobutyric_Acid             1.716164e-01 4.596840e-01 1.003943e-01 4.337670e-02
# Hippuric_Acid               7.877691e-03 2.205165e-01 2.704707e-02 4.484179e-03
# Methylmalonic_Acid          3.377079e-02 4.533513e-02 1.486640e-03 8.318840e-03
# LYSOC14_0                   1.473684e-01 1.990624e-02 3.492771e-04 4.341714e-01
# LYSOC16_1                   1.025838e-01 2.833578e-02 8.773741e-03 2.693227e-01
# LYSOC16_0                   2.040194e-02 9.259167e-02 1.611152e-02 2.457357e-01
# LYSOC17_0                   9.152639e-02 1.052438e-01 2.094311e-02 4.583592e-01
# LYSOC18_2                   3.539939e-01 8.421147e-03 4.795207e-02 2.065673e-01
# LYSOC18_1                   3.265264e-01 3.822174e-03 1.042723e-01 1.483692e-01
# LYSOC18_0                   2.790230e-01 1.918690e-03 2.821853e-01 7.173687e-02
# LYSOC20_4                   4.156917e-01 2.161667e-01 1.558381e-03 3.776341e-02
# LYSOC20_3                   8.271321e-02 7.886665e-04 3.797755e-02 2.724722e-01
# LYSOC24_0                   6.261881e-02 1.755515e-01 2.235092e-03 1.749915e-01
# LYSOC26_1                   2.157131e-01 2.818537e-02 1.642841e-05 2.829149e-01
# LYSOC26_0                   4.242179e-02 1.461255e-01 6.900283e-04 3.090036e-01
# LYSOC28_1                   3.809248e-01 9.026363e-04 2.079378e-02 2.879538e-01
# LYSOC28_0                   1.127663e-01 4.706245e-02 3.466508e-03 2.406656e-01
# X14_1SMOH                   4.048869e-01 1.198599e-02 1.676876e-01 1.861419e-01
# X16_1SM                     5.104104e-01 2.043626e-01 1.308879e-01 1.683969e-04
# X16_0SM                     4.800235e-01 2.370839e-01 1.275274e-01 7.290643e-03
# X16_1SMOH                   1.659172e-01 2.830790e-01 2.013510e-01 4.802241e-02
# X18_1SM                     4.643575e-01 1.006039e-02 1.225305e-01 1.466554e-01
# PC32_2AA                    1.575887e-01 1.628239e-01 2.545967e-01 1.794987e-01
# X18_0SM                     2.233968e-01 1.446913e-01 2.827268e-01 1.083679e-01
# X20_2SM                     3.481619e-01 2.960147e-02 2.193593e-01 2.516830e-01
# PC36_0AE                    3.664594e-01 2.592386e-01 6.692052e-02 8.544780e-02
# PC36_6AA                    3.659986e-01 8.482214e-04 9.148171e-03 4.799364e-01
# PC36_0AA                    1.640608e-01 4.499764e-01 8.138219e-02 6.524916e-02
# X22_2SMOH                   5.318448e-01 9.187496e-02 1.975237e-01 5.192516e-02
# X22_1SMOH                   6.566057e-01 4.093907e-02 4.492622e-02 1.407184e-01
# PC38_6AA                    4.424474e-01 2.161187e-03 1.125276e-03 3.640580e-01
# PC38_0AA                    5.375343e-01 1.099284e-01 7.665322e-02 1.584939e-01
# PC40_6AE                    5.073599e-01 2.664773e-03 9.519573e-04 3.031077e-01
# X24_1SMOH                   3.334428e-01 6.272836e-02 2.984488e-04 2.860980e-01
# PC40_6AA                    3.446157e-01 7.163986e-03 2.406034e-01 7.520931e-02
# PC40_2AA                    4.259480e-01 1.587454e-03 1.599943e-02 3.882686e-01
# PC401AA                     4.461028e-01 4.791111e-02 1.548684e-02 2.559316e-01
# C2                          2.522056e-01 7.479168e-02 6.342205e-04 8.010373e-03
# C3_1                        2.439159e-01 8.376191e-02 4.232614e-02 2.296609e-03
# C3                          3.401157e-01 3.386243e-01 1.829645e-07 1.417983e-03
# C4_1                        2.849750e-01 1.656078e-01 4.292503e-02 1.055777e-01
# C4                          6.229302e-01 8.714390e-02 5.779331e-02 6.733897e-02
# C3OH                        1.640252e-01 1.675651e-02 4.416272e-04 1.228995e-01
# C5_1                        2.134641e-01 2.677783e-01 7.443226e-03 1.411105e-03
# C5                          2.319690e-01 3.896890e-01 4.109781e-03 1.343984e-02
# C4OH                        3.259892e-01 6.150263e-04 6.699691e-02 1.602305e-02
# C6_1                        1.464529e-01 4.064792e-03 3.673180e-01 1.588527e-02
# C6                          4.468929e-01 1.126590e-01 1.110209e-01 2.859578e-02
# C5OH                        2.274386e-01 4.082344e-03 8.767761e-02 4.309948e-02
# C5_1DC                      4.771623e-01 3.623626e-02 6.377449e-02 2.154888e-02
# C5DC                        3.375926e-01 8.181409e-03 2.134578e-02 8.399347e-03
# C8                          5.443628e-01 1.799019e-01 6.446544e-03 5.865198e-02
# C5MDC                       2.068581e-01 1.478169e-01 3.261142e-03 5.800129e-02
# C9                          5.606298e-01 2.242021e-01 5.429673e-03 5.238310e-02
# C7DC                        8.124849e-03 1.982719e-02 1.006783e-03 1.906587e-02
# C10_2                       1.577092e-01 3.735656e-01 5.422440e-05 2.429558e-02
# C10_1                       1.716292e-01 7.101950e-03 3.509552e-02 1.429135e-04
# C10                         4.470929e-01 2.828892e-01 7.129489e-03 8.499362e-02
# C12_1                       6.052336e-01 6.244716e-02 8.033517e-02 5.325281e-02
# C12                         4.849933e-01 3.002773e-01 5.710370e-02 8.138650e-02
# C14_2                       5.600518e-01 2.061535e-01 8.003476e-02 1.307423e-02
# C14_1                       6.848581e-01 1.840854e-01 3.471762e-02 4.336938e-02
# C14                         5.583349e-01 2.287508e-01 7.686479e-02 7.778518e-02
# C12DC                       3.470130e-01 4.644267e-01 5.864876e-03 2.340911e-02
# C14_2OH                     4.890272e-01 3.951027e-01 1.267072e-02 5.641178e-02
# C14_1OH                     5.400849e-01 2.378337e-01 9.632028e-02 6.943019e-02
# C16_2                       4.239223e-01 4.090575e-01 1.471140e-03 7.871711e-02
# C16_1                       4.780804e-01 3.736244e-01 2.672021e-02 7.539067e-02
# C16                         4.755391e-01 3.000927e-01 1.190593e-01 4.060269e-02
# C16_2OH                     5.059506e-01 3.119989e-01 9.677961e-02 3.912183e-02
# C16_1OH                     5.251140e-01 2.664165e-01 1.121847e-01 2.208223e-02
# C16OH                       5.767034e-01 1.742745e-01 5.646473e-02 7.650161e-02
# C18_2                       5.455949e-01 3.378650e-01 1.720171e-02 6.036477e-02
# C18_1                       4.981646e-01 3.732696e-01 4.914858e-02 4.661492e-02
# C18                         4.462272e-01 3.532822e-01 9.168696e-02 1.481328e-02
# C18_1OH                     5.027644e-01 2.987137e-01 8.710412e-02 6.142244e-02
# Arginine_Average            3.955579e-02 4.580011e-02 6.922480e-01 2.794646e-06
# Betaine_Average             7.871843e-03 5.361235e-02 6.427820e-02 3.150826e-02
# C0_Average                  3.198092e-01 7.813186e-02 2.481698e-02 6.824724e-03
# Choline_Average             6.434478e-02 1.422426e-01 2.215793e-01 9.296654e-02
# Citrate_Average             4.600342e-01 1.289567e-01 7.374332e-03 7.062063e-02
# Creatine_Average            2.679483e-01 1.237782e-01 6.577600e-02 4.485341e-03
# Creatinine_Average          2.491709e-01 2.150025e-01 3.634238e-02 8.734872e-03
# Glucose_Average             1.880569e-02 5.279475e-02 3.009241e-01 1.921065e-01
# Glutamate_Average           1.869740e-01 3.491660e-02 3.584083e-01 2.126615e-02
# Glycine_Average             2.216425e-01 1.285835e-03 3.682173e-01 7.954172e-02
# Histidine_Average           2.430933e-01 9.068832e-02 2.321954e-01 1.368327e-02
# Isoleucine_Average          2.527086e-01 3.215364e-01 5.883742e-02 3.902142e-02
# Lactate_Average             3.678892e-01 1.833208e-01 6.311918e-02 1.189947e-04
# Leucine                     9.337532e-02 1.130190e-01 2.952199e-01 3.848890e-02
# Methionine                  1.467227e-01 5.041819e-02 5.844214e-01 4.714054e-02
# Phenylalanine_Average       1.918131e-01 7.144976e-02 5.104755e-01 2.599323e-02
# Proline                     1.794817e-01 5.383468e-01 3.770971e-02 9.175617e-03
# Pyruvate_Average            6.445859e-02 1.230973e-01 3.383102e-01 6.895765e-05
# Succinate_Average           7.912189e-02 6.811717e-02 2.726929e-01 7.542844e-03
# Threonine_Average           2.257272e-01 2.495887e-02 2.331294e-01 6.945048e-03
# Tyrosine_Average            1.117401e-01 6.968123e-02 2.163986e-01 6.868702e-02
# Valine_Average              2.858697e-01 4.069577e-01 4.835614e-02 4.395126e-02
# X1_Methylhistidine          6.709359e-03 3.554396e-02 6.681928e-03 1.300829e-02
# X2_Hydroxybutyrate          7.073696e-02 2.337478e-01 1.193671e-02 1.412055e-01
# X2_Hydroxyisovaleric_Acid   2.818684e-02 5.917868e-03 5.841272e-02 5.142625e-02
# X3_Hydroxybutyrate          5.425520e-01 3.743162e-02 5.526179e-02 1.983183e-01
# ADP                         1.380284e-01 2.332312e-01 5.946722e-03 4.926138e-03
# AMP                         8.057236e-03 3.130701e-02 5.878737e-03 1.151030e-02
# ATP                         5.907253e-02 1.144078e-02 2.163293e-02 1.516704e-02
# Acetamide                   1.657998e-01 2.842499e-01 1.561129e-02 1.980479e-03
# Acetate                     1.819616e-01 1.552500e-02 2.295865e-03 2.021110e-01
# Acetoacetate                7.299253e-03 1.212042e-03 4.931055e-02 5.447813e-02
# Acetone                     2.584069e-01 4.376439e-02 5.663827e-02 1.925654e-01
# Adenosine                   1.780821e-02 1.047580e-01 5.921217e-02 9.558380e-02
# Alanine                     4.180723e-02 5.388929e-06 1.147538e-01 2.044830e-02
# Aspartate                   4.131921e-02 3.306090e-01 3.388119e-04 9.535191e-02
# Creatine_Phosphate          1.788787e-06 4.665618e-02 4.787913e-02 2.868493e-05
# Cytidine                    4.141663e-03 9.334203e-02 2.533146e-02 5.945855e-02
# Dimethyl_Sulfone            7.220031e-02 2.153584e-01 1.921547e-02 1.022530e-01
# Ethanol                     4.102335e-02 9.273571e-03 1.930387e-02 2.471506e-04
# Ethanolamine                4.052291e-02 2.993121e-01 3.648271e-02 8.627586e-03
# Formate                     3.408479e-03 3.336864e-02 9.967768e-04 1.262575e-02
# Glycerol                    7.123634e-04 3.507826e-03 9.638236e-03 1.744804e-03
# Guanidoacetate              1.353537e-02 6.902484e-02 6.593532e-02 2.194028e-02
# Guanosine                   5.163945e-03 1.237157e-01 9.430076e-02 3.530322e-02
# Hypoxanthine                8.300191e-04 2.109720e-01 2.411332e-02 3.432450e-04
# IMP                         7.321697e-02 2.134321e-01 3.325556e-04 7.121798e-02
# Inosine                     1.699075e-03 3.155909e-01 2.198906e-02 1.878514e-02
# Isopropanol                 4.008083e-02 8.456347e-02 1.057784e-02 2.166103e-02
# Malate                      3.336391e-02 1.327735e-02 1.989384e-02 1.411883e-01
# Malonate                    4.451460e-02 1.104218e-02 7.697088e-05 1.023690e-01
# Mannose                     5.981662e-03 3.936674e-02 2.565272e-01 6.936208e-02
# Methanol                    1.078414e-03 7.536737e-02 7.018664e-04 2.217943e-02
# N_N_Dimethylglycine         9.299336e-04 1.665634e-01 3.055820e-02 3.272300e-02
# Nicotinurate                1.368070e-01 2.968353e-02 7.617129e-02 1.637474e-01
# O_Acetylcarnitine           1.958172e-01 4.624220e-02 2.094357e-03 1.210627e-02
# O_Phosphocholine            9.010757e-03 1.839286e-01 1.022226e-01 4.339387e-02
# Propylene_Glycol            5.106236e-02 8.322818e-02 9.266172e-02 2.084056e-01
# Taurine                     7.867594e-04 2.270714e-03 7.757983e-03 9.048942e-02
# Uridine                     1.629299e-01 1.093183e-01 1.419338e-01 1.496179e-01
# sn_Glycero_3_Phosphocholine 9.458329e-03 1.529345e-01 1.142918e-01 1.246086e-02
# Beta.Alanine                1.638183e-03 1.131128e-02 1.237718e-02 1.148710e-02
# Dim.5
# Serine                      5.728753e-03
# Putrescine                  4.604652e-02
# Trans_Hydroxyproline        1.583288e-02
# Glutamine                   1.100630e-04
# Alpha_Aminoadipic_Acid      2.300677e-03
# Methionine_Sulfoxide        4.089911e-03
# Acetyl_Ornithine            4.754080e-02
# Citrulline                  4.624734e-01
# Asymmetric_Dimethylarginine 9.311490e-04
# Total_Dimethylarginine      6.400371e-03
# Tryptophan                  4.226510e-02
# Ornithine                   4.415265e-02
# Lysine                      6.974638e-02
# Spermidine                  6.936848e-03
# Spermine                    1.111444e-02
# Sarcosine                   7.758361e-03
# Methylhistidine             2.819009e-02
# Beta_Hydroxybutyric_Acid    6.086237e-04
# Alpha_Ketoglutaric_Acid     4.713171e-02
# Butyric_acid                6.250853e-02
# Propionic_Acid              2.256226e-01
# Fumaric_Acid                1.475420e-01
# Isobutyric_Acid             1.521667e-03
# Hippuric_Acid               1.168942e-03
# Methylmalonic_Acid          1.287533e-02
# LYSOC14_0                   1.337158e-01
# LYSOC16_1                   3.319693e-01
# LYSOC16_0                   2.921726e-01
# LYSOC17_0                   1.511542e-01
# LYSOC18_2                   2.513744e-02
# LYSOC18_1                   2.965811e-02
# LYSOC18_0                   3.239147e-05
# LYSOC20_4                   2.542805e-02
# LYSOC20_3                   5.513979e-02
# LYSOC24_0                   1.781933e-02
# LYSOC26_1                   4.770632e-02
# LYSOC26_0                   1.016331e-02
# LYSOC28_1                   1.756686e-02
# LYSOC28_0                   7.096683e-02
# X14_1SMOH                   4.741021e-05
# X16_1SM                     4.472435e-03
# X16_0SM                     4.585544e-03
# X16_1SMOH                   6.227539e-02
# X18_1SM                     2.485961e-05
# PC32_2AA                    2.660572e-03
# X18_0SM                     2.132414e-03
# X20_2SM                     5.984130e-03
# PC36_0AE                    9.309850e-02
# PC36_6AA                    1.443347e-02
# PC36_0AA                    8.705657e-02
# X22_2SMOH                   6.529036e-03
# X22_1SMOH                   8.619909e-03
# PC38_6AA                    4.806828e-03
# PC38_0AA                    2.373535e-02
# PC40_6AE                    1.891356e-04
# X24_1SMOH                   8.597998e-03
# PC40_6AA                    3.492877e-02
# PC40_2AA                    5.176506e-04
# PC401AA                     6.669560e-03
# C2                          2.055316e-01
# C3_1                        6.878232e-03
# C3                          8.951079e-05
# C4_1                        1.707975e-02
# C4                          2.428738e-04
# C3OH                        4.246511e-03
# C5_1                        7.846092e-02
# C5                          2.180058e-05
# C4OH                        3.127877e-03
# C6_1                        1.537418e-02
# C6                          1.111928e-02
# C5OH                        3.731916e-02
# C5_1DC                      4.004405e-04
# C5DC                        1.365024e-01
# C8                          1.192139e-02
# C5MDC                       2.105657e-04
# C9                          2.406551e-04
# C7DC                        4.287917e-02
# C10_2                       2.811286e-02
# C10_1                       4.722000e-03
# C10                         6.313873e-03
# C12_1                       1.082689e-02
# C12                         1.206345e-03
# C14_2                       7.959496e-03
# C14_1                       1.236374e-03
# C14                         1.200001e-03
# C12DC                       1.383940e-04
# C14_2OH                     7.642833e-04
# C14_1OH                     1.629905e-05
# C16_2                       1.571402e-03
# C16_1                       2.873885e-05
# C16                         4.317966e-04
# C16_2OH                     2.472110e-05
# C16_1OH                     1.251185e-04
# C16OH                       4.691352e-03
# C18_2                       1.588343e-03
# C18_1                       2.972297e-04
# C18                         1.267974e-04
# C18_1OH                     4.575972e-04
# Arginine_Average            6.139814e-02
# Betaine_Average             2.984297e-01
# C0_Average                  8.006984e-02
# Choline_Average             2.176086e-02
# Citrate_Average             2.537485e-02
# Creatine_Average            7.690130e-03
# Creatinine_Average          4.305326e-02
# Glucose_Average             8.186440e-03
# Glutamate_Average           3.417226e-02
# Glycine_Average             1.015026e-03
# Histidine_Average           2.115885e-02
# Isoleucine_Average          8.785870e-03
# Lactate_Average             1.692005e-03
# Leucine                     6.511011e-04
# Methionine                  4.527834e-05
# Phenylalanine_Average       2.615126e-03
# Proline                     1.594346e-02
# Pyruvate_Average            7.157244e-02
# Succinate_Average           2.439249e-01
# Threonine_Average           1.198526e-01
# Tyrosine_Average            3.455073e-02
# Valine_Average              2.183729e-02
# X1_Methylhistidine          6.126701e-02
# X2_Hydroxybutyrate          8.015334e-02
# X2_Hydroxyisovaleric_Acid   4.133397e-02
# X3_Hydroxybutyrate          1.498130e-04
# ADP                         9.749080e-02
# AMP                         1.717972e-03
# ATP                         4.493587e-02
# Acetamide                   3.325850e-02
# Acetate                     4.902288e-02
# Acetoacetate                1.126362e-01
# Acetone                     4.735801e-02
# Adenosine                   1.569909e-01
# Alanine                     2.194894e-02
# Aspartate                   9.977000e-03
# Creatine_Phosphate          1.535580e-02
# Cytidine                    1.052956e-03
# Dimethyl_Sulfone            5.080823e-02
# Ethanol                     9.974920e-02
# Ethanolamine                2.402234e-01
# Formate                     1.035512e-01
# Glycerol                    1.904027e-03
# Guanidoacetate              7.001570e-02
# Guanosine                   3.225340e-01
# Hypoxanthine                2.079903e-01
# IMP                         5.128455e-03
# Inosine                     3.069147e-03
# Isopropanol                 7.247110e-04
# Malate                      4.853850e-03
# Malonate                    3.914630e-02
# Mannose                     1.108061e-01
# Methanol                    2.648608e-01
# N_N_Dimethylglycine         2.481559e-03
# Nicotinurate                1.587585e-01
# O_Acetylcarnitine           2.490742e-01
# O_Phosphocholine            7.804433e-02
# Propylene_Glycol            3.154744e-03
# Taurine                     6.954404e-02
# Uridine                     2.917731e-02
# sn_Glycero_3_Phosphocholine 5.704141e-02
# Beta.Alanine                6.140838e-03


#### graph of variables default plot ####
fviz_pca_var(DM_pca_FM3, col.var = "black")
pca18_biplot_blackno0<- fviz_pca_var(DM_pca_FM3, col.var = "black")
pca18_biplot_blackno0
# Open the graphics device
png("pca18_biplot_blackno0")
print(pca18_biplot_blackno0)
dev.off()
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(DM_pca_FM3, col.var = "black", axes = c(2,3))
pca18_biplot_blackpc2_3no0<- fviz_pca_var(DM_pca_FM3, col.var = "black", axes= c(2,3))
pca18_biplot_blackpc2_3no0
# Open the graphics device
png("pca18_biplot_blackpc2_3no0")
print(pca18_biplot_blackpc2_3no0)
dev.off()

fviz_pca_var(DM_pca_FM3, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4no0<- fviz_pca_var(DM_pca_FM3, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4no0
png("pca18_biplot_black3_4no0")
print(pca18_biplot_black3_4no0)
dev.off()

fviz_pca_var(DM_pca_FM3, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5no0<- fviz_pca_var(DM_pca_FM3, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5no0
png("pca18_biplot_black4_5no0")
print("pca18_biplot_black4_5no0")
dev.off()


#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var30$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(DM_pca_FM3, choice = "var", axes = 1:2)
bargraph18_cos2no0<- fviz_cos2(DM_pca_FM3, choice = "var", axes = 1:2)
bargraph18_cos2no0
png("bargraph18_cos2no0")
print(bargraph18_cos2no0)
dev.off()
# Color by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM3, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2no0<-fviz_pca_var(DM_pca_FM3, col.var = "cos2",
                                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                       repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2no0
png("biplot18_colour_cos2no0")
print(biplot18_colour_cos2no0)
dev.off()
# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM3, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2no0<- fviz_pca_var(DM_pca_FM3, alpha.var = "cos2",
                                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                       repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2no0
png("biplot18_black_cos2no0")
print(biplot18_black_cos2no0)
dev.off()

#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.



# visualize contributions as corrplot

corrplot(var30$contrib, is.corr=FALSE,
         tl.col = "black", tl.cex = 1, 
         cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
         mar = c(1, 1, 1, 1)
) 
corrplot18_contrib2no0 <- corrplot(var30$contrib, is.corr=FALSE,
                                    tl.col = "black", tl.cex = 1, 
                                    cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
                                    mar = c(1, 1, 1, 1)
)
corrplot18_contrib2no0
pdf("mcorrplot18_contrib2no43", height = 7, width =5)
png("corrplot18_contrib2no0")
print(corrplot18_contrib2no0)
dev.off()
# Barplot of contributions to a component
# The red dashed line on the graph above indicates the expected average contribution. 
# If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/32 = 3.125%.
# For a given component, a variable with a contribution larger than this cutoff 
# could be considered as important in contributing to the component.

#### contributions of variables ####
# to PC1
contrib18_PC1no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = 1, top = 80,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 19, label = expression(bold("PC1")), size = 5)
contrib18_PC1no0
png("contrib18_PC1no43")
print(contrib18_PC1no0)
dev.off()

# contributions of variables to PC2
contrib18_PC2no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = 2, top = 80,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 8.5, x = 19, label = expression(bold("PC2")), size = 5)
contrib18_PC2no0
png("contrib18_PC2no0")
print(contrib18_PC2no0)
dev.off()
# contributions of variables to PC3
contrib18_PC3no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = 3, top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 10, x = 19, label = expression(bold("PC3")), size = 5)
contrib18_PC3no0
png("contrib18_PC3no0")
print(contrib18_PC3no0)
dev.off()

# contributions of variables to PC4
contrib18_PC4no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = 4, top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 15, x = 19, label = expression(bold("PC4")), size = 5)
contrib18_PC4no0
png("contrib18_PC4no0")
print(contrib18_PC4no0)
dev.off()
# contributions of variables to PC5
contrib18_PC5no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = 5, top = 40,
                                  fill = "#6baed6", color = "#2171b5",
                                  title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 25, x = 19, label = expression(bold("PC5")), size = 5)
contrib18_PC5no0
png("contrib18_PC5no0")
print(contrib18_PC5no0)
dev.off()
# contributions of variables to PC1 and PC2
contrib18_PC1_PC2no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = c(1,2), top = 80,
                                      fill = "#6baed6", color = "#2171b5",
                                      title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 5, x = 18, label = expression(bold("PC1 & PC2")), size = 5)
contrib18_PC1_PC2no0
png("contrib18_PC1_PC2no0")
print(contrib18_PC1_PC2no0)
dev.off()
# contributions of variables to PC2 and PC3
contrib18_PC2_PC3no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = c(2,3), top = 40,
                                      fill = "#6baed6", color = "#2171b5",
                                      title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 6.5, x = 18, label = expression(bold("PC2 & PC3")), size = 5)
contrib18_PC2_PC3no0
png("contrib18_PC2_PC3no0")
print(contrib18_PC2_PC3no0)
dev.off()
# contributions of variables to PC3 and PC4
contrib18_PC3_PC4no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = c(3,4), top = 40,
                                      fill = "#6baed6", color = "#2171b5",
                                      title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 18, label = expression(bold("PC3 & PC4")), size = 5)
contrib18_PC3_PC4no0
png("contrib18_PC3_PC4no0")
print(contrib18_PC3_PC4no0)
dev.off()
# contributions of variables to PC4 and PC5
contrib18_PC4_PC5no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = c(4,5), top = 40,
                                      fill = "#6baed6", color = "#2171b5",
                                      title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 13, x = 18, label = expression(bold("PC4 & PC5")), size = 5)
contrib18_PC4_PC5no0
png("contrib18_PC4_PC5no0")
print(contrib18_PC4_PC5no0)
dev.off()
# contributions of variables to PC1:PC3
contrib18_PC1_to_PC3no0 <- fviz_contrib(DM_pca_FM3, choice = "var", axes = c(1:3), top = 40,
                                         fill = "#6baed6", color = "#2171b5",
                                         title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 4.5, x = 18, label = expression(bold("PC1 to PC3")), size = 5)
contrib18_PC1_to_PC3no0
png("contrib18_PC1_PC3no0")
print(contrib18_PC1_to_PC3no0)
dev.off()
#### colour variable colors using their contributions ####

# PC1 and PC2
varplot18_contrib_PC1_PC2no0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                              repel = TRUE,  # avoid text overlapping
                                              title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2no0
# pdf("varplot18_contrib_PC1_PC2no43")
# print()
png("varplot18_contrib_PC1_PC2no0")
print(varplot18_contrib_PC1_PC2no0)
dev.off()

# PC1 and PC2 select var
varplot18_contrib_PC1_PC2_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                     repel = TRUE,  # avoid text overlapping
                                                     title = "",
                                                     select.var = list(contrib = 60) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2_selectno0
# pdf("varplot18_contrib_PC1_PC2_selectno43")
# dev.off()
png("varplot18_contrib_PC1_PC2_selectno0")
print(varplot18_contrib_PC1_PC2_selectno0)
dev.off()

# PC2 and PC3
varplot18_contrib_PC2_PC3no0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                              axes = c(2,3),
                                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                              repel = TRUE,  # avoid text overlapping
                                              title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3no0
# pdf("varplot18_contrib_PC2_PC3no0")
# dev.off()
png("varplot18_contrib_PC2_PC3no0")
print(varplot18_contrib_PC2_PC3no0)
dev.off()

# PC2 and PC3 select var
varplot18_contrib_PC2_PC3_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                                     axes = c(2,3),
                                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                     repel = TRUE,  # avoid text overlapping
                                                     title = "",
                                                     select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3_selectno0
# pdf("varplot18_contrib_PC2_PC3_selectno0")
# dev.off()
png("varplot18_contrib_PC2_PC3_selectno0")
print(varplot18_contrib_PC2_PC3_selectno0)
dev.off()
# PC4 and PC3
varplot18_contrib_PC4_PC3no0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                              axes = c(3,4),
                                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                              repel = TRUE,  # avoid text overlapping
                                              title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3no43
# pdf("varplot18_contrib_PC4_PC3no0")
# dev.off()
png("varplot18_contrib_PC4_PC3no0")
print(varplot18_contrib_PC4_PC3no0)
dev.off()

# PC4 and PC3 select var
varplot18_contrib_PC4_PC3_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                                     axes = c(3,4),
                                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                     repel = TRUE,  # avoid text overlapping
                                                     title = "",
                                                     select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3_selectno0
# pdf("varplot18_contrib_PC4_PC3_selectno0")
# dev.off()
png("varplot18_contrib_PC4_PC3_selectno0")
print(varplot18_contrib_PC4_PC3_selectno0)
dev.off()


# PC4 and PC5
varplot18_contrib_PC4_PC5no0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                              axes = c(4,5),
                                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                              repel = TRUE,  # avoid text overlapping
                                              title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5no0
# pdf("varplot18_contrib_PC4_PC5no0")
# dev.off()
png("varplot18_contrib_PC4_PC5no0")
print(varplot18_contrib_PC4_PC5no0)
dev.off()

# PC4 and PC5 select var
varplot18_contrib_PC4_PC5_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = "contrib",  # colour by contributions to the PC
                                                     axes = c(4,5),
                                                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                     repel = TRUE,  # avoid text overlapping
                                                     title = "",
                                                     select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5_selectno0
# pdf("varplot18_contrib_PC4_PC5_selectno0")
# dev.off()
png("varplot18_contrib_PC4_PC5_selectno0")
print(varplot18_contrib_PC4_PC5_selectno0)
dev.off()

# install the colour brewer
library("RColorBrewer") # did not use this colours

#### cos2 selected variance trying plot ####

fviz_pca_var(DM_pca_FM3, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             select.var = list(cos2 = 50)
)
cos2_18_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                   repel = TRUE, # Avoid text overlapping
                                   select.var = list(cos2 = 50)
)
# pdf("cos2_18_select")
# dev.off()
png("cos2_18_selectno0")
print(cos2_18_selectno0)
dev.off()



#### colour variables by groups ####
#did not load the sheet of the metabolites
# PC1-2
#options(ggrepel.max.overlaps = Inf)
varplot18_met_PC1_PC2_selectno0 <- fviz_pca_var(DM_pca_FM3, col.var = met_groups$Category, 
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

dim.desc30 <- dimdesc(DM_pca_FM3, axes = c(1:5))

dim.desc30$Dim.1
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# C14_1                         0.8275615 1.206652e-14
# X22_1SMOH                     0.8103121 1.147373e-13
# C4                            0.7892593 1.339224e-12
# C12_1                         0.7779676 4.473340e-12
# C16OH                         0.7594099 2.808121e-11
# Beta_Hydroxybutyric_Acid      0.7512032 6.008954e-11
# C9                            0.7487522 7.499112e-11
# C14_2                         0.7483661 7.763637e-11
# C14                           0.7472181 8.603430e-11
# C18_2                         0.7386440 1.821678e-10
# C8                            0.7378095 1.956616e-10
# X3_Hydroxybutyrate            0.7365813 2.172555e-10
# C14_1OH                       0.7349047 2.503989e-10
# PC38_0AA                      0.7331673 2.897630e-10
# X22_2SMOH                     0.7292769 4.001963e-10
# C16_1OH                       0.7246475 5.835075e-10
# X16_1SM                       0.7144301 1.306521e-09
# PC40_6AE                      0.7122920 1.539812e-09
# C16_2OH                       0.7113021 1.660684e-09
# C18_1OH                       0.7090588 1.968588e-09
# C18_1                         0.7058077 2.511833e-09
# C14_2OH                       0.6993048 4.050243e-09
# C12                           0.6964146 4.988236e-09
# X16_0SM                       0.6928373 6.433805e-09
# C16_1                         0.6914336 7.102373e-09
# C5_1DC                        0.6907694 7.441090e-09
# C16                           0.6895934 8.078427e-09
# X18_1SM                       0.6814378 1.413615e-08
# Citrate_Average               0.6782582 1.749792e-08
# C10                           0.6686500 3.282147e-08
# C6                            0.6685005 3.313833e-08
# C18                           0.6680024 3.421474e-08
# PC401AA                       0.6679093 3.441954e-08
# PC38_6AA                      0.6651672 4.099503e-08
# PC40_2AA                      0.6526469 8.905319e-08
# C16_2                         0.6510932 9.780878e-08
# LYSOC20_4                     0.6447416 1.427205e-07
# X14_1SMOH                     0.6363072 2.326154e-07
# Total_Dimethylarginine        0.6351522 2.484244e-07
# Tryptophan                    0.6190577 6.041626e-07
# LYSOC28_1                     0.6171910 6.676088e-07
# Lactate_Average               0.6065387 1.166059e-06
# PC36_0AE                      0.6053589 1.238811e-06
# PC36_6AA                      0.6049782 1.263177e-06
# LYSOC18_2                     0.5949738 2.088642e-06
# X20_2SM                       0.5900524 2.658417e-06
# C12DC                         0.5890781 2.787141e-06
# PC40_6AA                      0.5870398 3.075411e-06
# C3                            0.5831944 3.696220e-06
# C5DC                          0.5810272 4.095575e-06
# X24_1SMOH                     0.5774450 4.844675e-06
# LYSOC18_1                     0.5714249 6.396459e-06
# C4OH                          0.5709546 6.535321e-06
# C0_Average                    0.5655168 8.357517e-06
# Spermine                      0.5519292 1.516573e-05
# Alpha_Ketoglutaric_Acid       0.5509637 1.580607e-05
# Spermidine                    0.5436292 2.155074e-05
# Asymmetric_Dimethylarginine   0.5361182 2.938426e-05
# Valine_Average                0.5346678 3.117078e-05
# C4_1                          0.5338305 3.224727e-05
# LYSOC18_0                     0.5282263 4.038152e-05
# Putrescine                    0.5257832 4.448696e-05
# Creatine_Average              0.5176372 6.111219e-05
# Acetone                       0.5083374 8.696047e-05
# Isoleucine_Average            0.5027013 1.071591e-04
# C2                            0.5022008 1.091459e-04
# Creatinine_Average            0.4991702 1.219145e-04
# Butyric_acid                  0.4972639 1.306317e-04
# C3_1                          0.4938784 1.475302e-04
# Histidine_Average             0.4930449 1.519863e-04
# Fumaric_Acid                  0.4835485 2.121875e-04
# C5                            0.4816316 2.267055e-04
# C5OH                          0.4769052 2.664483e-04
# Threonine_Average             0.4751076 2.831576e-04
# X18_0SM                       0.4726487 3.075574e-04
# Glycine_Average               0.4707892 3.272609e-04
# LYSOC26_1                     0.4644493 4.033653e-04
# C5_1                          0.4620217 4.365206e-04
# C5MDC                         0.4548165 5.499803e-04
# O_Acetylcarnitine             0.4425124 8.066657e-04
# Phenylalanine_Average         0.4379647 9.260162e-04
# Acetyl_Ornithine              0.4363911 9.708715e-04
# Glutamine                     0.4361412 9.781712e-04
# Glutamate_Average             0.4324049 1.093364e-03
# Acetate                       0.4265696 1.297740e-03
# Proline                       0.4236528 1.412204e-03
# C10_1                         0.4142816 1.843587e-03
# Isobutyric_Acid               0.4142661 1.844389e-03
# X16_1SMOH                     0.4073293 2.235858e-03
# PC36_0AA                      0.4050442 2.380116e-03
# C3OH                          0.4050003 2.382965e-03
# Uridine                       0.4036457 2.472420e-03
# C10_2                         0.3971261 2.946050e-03
# PC32_2AA                      0.3969744 2.957967e-03
# LYSOC14_0                     0.3838859 4.161492e-03
# Methionine                    0.3830440 4.251908e-03
# C6_1                          0.3826916 4.290268e-03
# Nicotinurate                  0.3698743 5.909208e-03
# Alpha_Aminoadipic_Acid        0.3524402 8.955464e-03
# LYSOC28_0                     0.3358070 1.304718e-02
# Tyrosine_Average              0.3342755 1.349407e-02
# Serine                        0.3311016 1.446219e-02
# LYSOC16_1                     0.3202871 1.821990e-02
# Leucine                       0.3055738 2.463996e-02
# LYSOC17_0                     0.3025333 2.618027e-02
# Sarcosine                     0.3016138 2.666171e-02
# Trans_Hydroxyproline          0.2937939 3.106393e-02
# LYSOC20_3                     0.2875990 3.496716e-02
# Succinate_Average             0.2812861 3.935438e-02
# IMP                           0.2705863 4.782090e-02
# Dimethyl_Sulfone              0.2687011 4.945635e-02
# ADP                          -0.3715217 5.674952e-03
# Acetamide                    -0.4071852 2.244718e-03
dim.desc30$Dim.2
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# C12DC                         0.6814886 1.408777e-08
# C16_2                         0.6395761 1.928328e-07
# C14_2OH                       0.6285720 3.594581e-07
# C16_1                         0.6112482 9.135633e-07
# C10_2                         0.6112001 9.158594e-07
# C18_1                         0.6109579 9.275152e-07
# C18                           0.5943755 2.151263e-06
# C18_2                         0.5812616 4.050522e-06
# C16_2OH                       0.5585686 1.137184e-05
# C12                           0.5479756 1.794968e-05
# C16                           0.5478072 1.807819e-05
# Ethanolamine                  0.5470943 1.863150e-05
# C18_1OH                       0.5465471 1.906675e-05
# Acetamide                     0.5331510 3.314594e-05
# C10                           0.5318733 3.489843e-05
# C16_1OH                       0.5161555 6.468966e-05
# C14_1OH                       0.4876820 1.837227e-04
# ADP                           0.4829401 2.167008e-04
# C14                           0.4782790 2.542893e-04
# C9                            0.4734998 2.989032e-04
# Hypoxanthine                  0.4593169 4.763566e-04
# C14_2                         0.4540413 5.636534e-04
# C14_1                         0.4290518 1.206955e-03
# C8                            0.4241485 1.392136e-03
# LYSOC24_0                     0.4189887 1.614084e-03
# C16OH                         0.4174619 1.685559e-03
# N_N_Dimethylglycine           0.4081218 2.187677e-03
# PC32_2AA                      0.4035144 2.481251e-03
# sn_Glycero_3_Phosphocholine   0.3910685 3.456396e-03
# C5MDC                         0.3844697 4.099791e-03
# LYSOC26_0                     0.3822637 4.337256e-03
# X18_0SM                       0.3803831 4.549166e-03
# Choline_Average               0.3771506 4.934735e-03
# Guanosine                     0.3517324 9.103647e-03
# Pyruvate_Average              0.3508522 9.290894e-03
# C6                            0.3356471 1.309324e-02
# LYSOC17_0                     0.3244130 1.669848e-02
# Adenosine                     0.3236635 1.696653e-02
# LYSOC16_0                     0.3042888 2.528136e-02
# Isopropanol                   0.2907980 3.290369e-02
# Lysine                        0.2829067 3.818707e-02
# C0_Average                    0.2795208 4.065933e-02
# C2                            0.2734807 4.539612e-02
# Methylhistidine               0.2731142 4.569748e-02
# Acetyl_Ornithine             -0.2726068 4.611739e-02
# Methanol                     -0.2745312 4.454131e-02
# Propylene_Glycol             -0.2884930 3.438003e-02
# C3_1                         -0.2894165 3.378205e-02
# C4                           -0.2952015 3.022966e-02
# Histidine_Average            -0.3011450 2.691000e-02
# X22_2SMOH                    -0.3031088 2.588266e-02
# Cytidine                     -0.3055193 2.466687e-02
# Alpha_Aminoadipic_Acid       -0.3279525 1.548090e-02
# Uridine                      -0.3306331 1.461002e-02
# PC38_0AA                     -0.3315545 1.432052e-02
# Leucine                      -0.3361829 1.293943e-02
# Butyric_acid                 -0.3447482 1.068354e-02
# Creatine_Average             -0.3518212 9.084950e-03
# Citrate_Average              -0.3591054 7.659505e-03
# C4_1                         -0.4069494 2.259287e-03
# Lactate_Average              -0.4281597 1.238904e-03
# O_Phosphocholine             -0.4288690 1.213442e-03
# Spermine                     -0.4431531 7.910126e-04
# Total_Dimethylarginine       -0.4466862 7.095149e-04
# X16_1SM                      -0.4520649 5.999165e-04
# IMP                          -0.4619871 4.370109e-04
# Creatinine_Average           -0.4636837 4.135671e-04
# Dimethyl_Sulfone             -0.4640672 4.084271e-04
# LYSOC20_4                    -0.4649373 3.969824e-04
# Hippuric_Acid                -0.4695919 3.405475e-04
# X2_Hydroxybutyrate           -0.4834747 2.127307e-04
# X16_0SM                      -0.4869126 1.887413e-04
# Putrescine                   -0.5037498 1.031037e-04
# PC36_0AE                     -0.5091548 8.434001e-05
# C5_1                         -0.5174730 6.149956e-05
# X16_1SMOH                    -0.5320517 3.464878e-05
# Tryptophan                   -0.5359409 2.959743e-05
# Spermidine                   -0.5412832 2.376100e-05
# Inosine                      -0.5617748 9.873859e-06
# Isoleucine_Average           -0.5670418 7.803950e-06
# Asymmetric_Dimethylarginine  -0.5671212 7.776072e-06
# Aspartate                    -0.5749861 5.430543e-06
# C3                           -0.5819144 3.927469e-06
# C5                           -0.6242507 4.560404e-07
# Valine_Average               -0.6379324 2.119633e-07
# PC36_0AA                     -0.6708028 2.856436e-08
# Isobutyric_Acid              -0.6780000 1.780167e-08
# Proline                      -0.7337212 2.766176e-10
# Trans_Hydroxyproline         -0.7691570 1.092958e-11
dim.desc30$Dim.3
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Arginine_Average              0.8320144 6.477919e-15
# Serine                        0.8264141 1.412382e-14
# Ornithine                     0.7909367 1.112571e-12
# Methionine                    0.7644746 1.729498e-11
# Lysine                        0.7522273 5.473564e-11
# Methionine_Sulfoxide          0.7260839 5.194972e-10
# Phenylalanine_Average         0.7144757 1.301934e-09
# Sarcosine                     0.6643976 4.304253e-08
# Glutamine                     0.6161849 7.043398e-07
# Glycine_Average               0.6068091 1.149960e-06
# Glutamate_Average             0.5986721 1.737761e-06
# Pyruvate_Average              0.5816444 3.977941e-06
# Glucose_Average               0.5485655 1.750633e-05
# Leucine                       0.5433415 2.181122e-05
# Alpha_Ketoglutaric_Acid       0.5375406 2.772433e-05
# X18_0SM                       0.5317206 3.511357e-05
# LYSOC18_0                     0.5312111 3.583997e-05
# Succinate_Average             0.5222001 5.120612e-05
# Mannose                       0.5064851 9.317677e-05
# PC32_2AA                      0.5045757 1.000087e-04
# Fumaric_Acid                  0.5027698 1.068898e-04
# PC40_6AA                      0.4905134 1.662840e-04
# Threonine_Average             0.4828347 2.174915e-04
# Histidine_Average             0.4818666 2.248789e-04
# Choline_Average               0.4707221 3.279934e-04
# X20_2SM                       0.4683581 3.547500e-04
# Tyrosine_Average              0.4651866 3.937571e-04
# Total_Dimethylarginine        0.3989426 2.806618e-03
# Uridine                       0.3767410 4.985583e-03
# X18_1SM                       0.3500436 9.465830e-03
# Alanine                       0.3387533 1.222306e-02
# sn_Glycero_3_Phosphocholine   0.3380707 1.240991e-02
# LYSOC18_1                     0.3229123 1.723884e-02
# O_Phosphocholine              0.3197228 1.843688e-02
# Asymmetric_Dimethylarginine   0.3185358 1.890044e-02
# Methylhistidine               0.3086876 2.314218e-02
# Guanosine                     0.3070843 2.390350e-02
# C5OH                          0.2961041 2.970456e-02
# Nicotinurate                  0.2759915 4.337500e-02
# Alpha_Aminoadipic_Acid        0.2755400 4.373289e-02
# PC38_0AA                     -0.2768632 4.269080e-02
# C14                          -0.2772450 4.239390e-02
# C14_2                        -0.2829041 3.818886e-02
# C12_1                        -0.2834346 3.781300e-02
# PC36_0AA                     -0.2852756 3.653216e-02
# Butyric_acid                 -0.2884364 3.441691e-02
# C18_1OH                      -0.2951341 3.026917e-02
# C18                          -0.3027985 2.604276e-02
# Propylene_Glycol             -0.3044039 2.522335e-02
# C14_1OH                      -0.3103551 2.237216e-02
# C16_2OH                      -0.3110942 2.203782e-02
# Isobutyric_Acid              -0.3168507 1.957550e-02
# C6                           -0.3331980 1.381632e-02
# C16_1OH                      -0.3349398 1.329864e-02
# Citrulline                   -0.3424372 1.125595e-02
# C16                          -0.3450497 1.061073e-02
# X16_0SM                      -0.3571098 8.029231e-03
# X16_1SM                      -0.3617843 7.186548e-03
# X14_1SMOH                    -0.4094967 2.106274e-03
# X22_2SMOH                    -0.4444364 7.604842e-04
# X16_1SMOH                    -0.4487215 6.660808e-04
# C6_1                         -0.6060677 1.194614e-06
dim.desc30$Dim.4
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# PC36_6AA                   0.6927744 6.462449e-09
# LYSOC17_0                  0.6770223 1.899737e-08
# LYSOC14_0                  0.6589168 6.065885e-08
# PC40_2AA                   0.6231120 4.852638e-07
# PC38_6AA                   0.6033722 1.370960e-06
# LYSOC26_0                  0.5558809 1.278699e-05
# PC40_6AE                   0.5505522 1.608650e-05
# LYSOC28_1                  0.5366132 2.879638e-05
# X24_1SMOH                  0.5348813 3.090170e-05
# LYSOC26_1                  0.5318975 3.486451e-05
# LYSOC20_3                  0.5219887 5.163020e-05
# LYSOC16_1                  0.5189631 5.806612e-05
# PC401AA                    0.5058968 9.523474e-05
# X20_2SM                    0.5016802 1.112483e-04
# LYSOC16_0                  0.4957173 1.381186e-04
# LYSOC28_0                  0.4905768 1.659111e-04
# LYSOC18_2                  0.4544967 5.555845e-04
# Glucose_Average            0.4382996 9.167133e-04
# Propionic_Acid             0.4374965 9.391630e-04
# X14_1SMOH                  0.4314416 1.124970e-03
# PC32_2AA                   0.4236729 1.411384e-03
# LYSOC24_0                  0.4183199 1.645054e-03
# PC38_0AA                   0.3981129 2.869563e-03
# LYSOC18_1                  0.3851872 4.025057e-03
# X18_1SM                    0.3829561 4.261451e-03
# X22_1SMOH                  0.3751245 5.190778e-03
# C3OH                       0.3505702 9.351578e-03
# X18_0SM                    0.3291928 1.507260e-02
# C4_1                       0.3249272 1.651669e-02
# Choline_Average            0.3049041 2.497247e-02
# PC36_0AE                   0.2923146 3.196139e-02
# PC40_6AA                   0.2742432 4.477436e-02
# C16_1                     -0.2745736 4.450706e-02
# C16OH                     -0.2765892 4.290484e-02
# C14                       -0.2788999 4.112663e-02
# C16_2                     -0.2805657 3.988270e-02
# Glycine_Average           -0.2820314 3.881393e-02
# C12                       -0.2852832 3.652698e-02
# C10                       -0.2915367 3.244188e-02
# Taurine                   -0.3008146 2.708616e-02
# Aspartate                 -0.3087911 2.309376e-02
# Adenosine                 -0.3091663 2.291887e-02
# Dimethyl_Sulfone          -0.3197703 1.841852e-02
# Malonate                  -0.3199516 1.834864e-02
# Malate                    -0.3757504 5.110477e-03
# X2_Hydroxybutyrate        -0.3757733 5.107557e-03
# Uridine                   -0.3868047 3.860979e-03
# Nicotinurate              -0.4046572 2.405349e-03
# Beta_Hydroxybutyric_Acid  -0.4253167 1.345848e-03
# Acetone                   -0.4388227 9.023513e-04
# X3_Hydroxybutyrate        -0.4453294 7.398708e-04
# Acetate                   -0.4495675 6.487437e-04
# Propylene_Glycol          -0.4565146 5.210719e-04
dim.desc30$Dim.5
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Citrulline           0.6800540 1.551664e-08
# LYSOC16_1            0.5761677 5.141248e-06
# Guanosine            0.5679208 7.500522e-06
# Betaine_Average      0.5462872 1.927676e-05
# LYSOC16_0            0.5405300 2.451395e-05
# Succinate_Average    0.4938876 1.474820e-04
# Ethanolamine         0.4901259 1.685780e-04
# Hypoxanthine         0.4560596 5.286799e-04
# Nicotinurate         0.3984451 2.844216e-03
# Adenosine            0.3962207 3.017811e-03
# LYSOC17_0            0.3887855 3.668105e-03
# Fumaric_Acid         0.3841120 4.137498e-03
# LYSOC14_0            0.3656717 6.545527e-03
# ADP                  0.3122352 2.153002e-02
# PC36_0AE             0.3051205 2.486464e-02
# PC36_0AA             0.2950535 3.031647e-02
# X2_Hydroxybutyrate   0.2831137 3.804004e-02
# O_Phosphocholine     0.2793642 4.077680e-02
# C5_1                -0.2801088 4.022079e-02
# C0_Average          -0.2829661 3.814477e-02
# Ethanol             -0.3158310 1.999384e-02
# Formate             -0.3217937 1.765123e-02
# Mannose             -0.3328756 1.391402e-02
# Acetoacetate        -0.3356131 1.310303e-02
# Threonine_Average   -0.3461973 1.033757e-02
# C5DC                -0.3694623 5.969108e-03
# C2                  -0.4533560 5.759971e-04
# Propionic_Acid      -0.4749974 2.842119e-04
# O_Acetylcarnitine   -0.4990734 1.223441e-04
# Methanol            -0.5146463 6.853037e-05


#### extract the results for individuals ####
ind30 <- get_pca_ind(DM_pca_FM3)
ind30
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"
ind30$coord
# Dim.1       Dim.2       Dim.3       Dim.4       Dim.5
# 1    9.2116333   6.2910780  10.2291741  1.86340370  1.06983549
# 2   -3.4137569   4.6073486  -1.7986937  0.63956119  3.66790799
# 3    0.4707443   5.6712781  -2.6870369  0.17233609  7.73417966
# 4   -2.0643010   4.3637829   2.7023771 -0.07821887 -2.89775573
# 5    3.1214893   0.7388319  10.9433300 -1.82617275 -0.81878123
# 6   -0.4539689   2.3682227   7.2885519 -3.74385104 -0.97433749
# 7   -1.4646399   2.8631609   3.0068901  0.62476492  1.12160982
# 8    6.6019316   6.4413704   7.5913729  0.02488315 -0.41683878
# 9    2.2690953   5.8449272   5.9081401 -3.89021671  2.09311083
# 10   6.4747118   8.0558665  -0.6491482 -5.70009410  0.71962143
# 11  -7.2059263   2.1639135   2.3287920 -5.02501726  0.61592956
# 12  -5.4115148   4.1706010  -1.9381619  3.36121660  3.97669479
# 13   0.7022370   3.8130435  -2.3147214  6.45044757  5.31885082
# 14  -7.8897138   2.0281416  -2.6340079 -1.99550038  2.06553226
# 15  -8.8801126   3.2372919  -3.9741191 -3.19353591  2.29667114
# 16   0.9536137   2.9965751   4.9139089  2.23722751 -0.17950490
# 17  -8.5343895   3.6089459  -4.4759106  0.20020680  4.08648142
# 18  -3.6643447   4.0883762  -4.0036465  3.62573332  1.41019828
# 19 -13.0420261   1.5143057  -5.4439720 -3.06561719 -0.12725115
# 20   3.5712104   3.6928472  -5.2726503 11.08861111  1.05653638
# 21  -6.9781251   1.5652691   3.8421680 -1.41482261 -3.06216924
# 22   0.3831574   4.1637649   6.2428386  0.82468009  0.29298552
# 23  -2.1133693  -1.3621190   1.5842089  1.10709243 -0.16217922
# 24  -2.7660416   2.9321336   2.0595454  0.73154448  0.13899921
# 25  -0.5334030  -3.0717235   0.5842253  3.37091391  0.90084320
# 26  -3.7829875  -1.6308516  -0.8257395 -1.14619535 -2.78295750
# 27   4.8211304  -1.4966092  -0.1001565  9.90353920 -3.69326851
# 28  -0.5201724   0.1111678  -2.8122386  6.08393422 -4.29582015
# 29  -7.3520779  -0.6619852  -2.8174385 -2.98233512 -5.44144165
# 30  -1.8591115  -1.1564114  -0.2701211 -0.24746346 -0.56609628
# 31  -1.4725389  -2.3722050  -0.6636488  1.35434067 -5.07293010
# 32   0.5751293  -2.2561379  -0.1417630  2.58064672 -2.65811744
# 33   3.1967429   1.2284919  -1.1172700  0.40397723 -3.36355077
# 34  -1.5473269  -0.4924011   2.0952803  1.10817284 -1.28195355
# 35  -0.4712374   2.7524519  -5.7081833 -6.11785553 -6.43272781
# 36   2.6800828   3.3138851  -0.3214744  2.15622990 -2.03189002
# 37  -4.6440215  -0.8606392  -1.3271870 -2.57200815 -1.82236153
# 38  -1.4275183   0.5836021  -2.9248633 -2.43314862 -1.55273609
# 39  -1.5127389   2.2153265  -1.6835304 -0.69605722 -3.40460804
# 40   0.9038152  -7.7869837   2.0140100 -1.84573998  0.74678969
# 41   2.7713425  -5.8581903  -0.1775224  2.57246833 -2.29212216
# 42   6.8846211  -8.2039745   0.2687651  2.18839542  0.61063549
# 43  30.3243838   7.3287541 -11.6876200 -6.66919947 -0.96182799
# 44   6.6767112  -8.5563683  -0.8633369  2.13647824  2.02407906
# 45  -1.6392361  -2.0680614  -1.3865348 -1.85009918  2.47384749
# 46   4.2634773  -7.5785682  -3.1634094 -5.13179344  1.97149737
# 47  -3.8655320  -2.3970200  -2.2657507  0.71353216 -0.78087472
# 48  -0.7467098  -4.9310402   0.1991646 -0.31277333 -0.69359625
# 49  -2.0952031  -5.1693901  -0.7624226 -2.73104986  1.79510454
# 50   3.3528699  -8.0871838  -2.0026788 -1.33067051  3.50640483
# 51  -2.2580163  -3.0282433  -2.7614582  1.79427331 -0.03559218
# 52   3.6035662  -6.0077117   0.2160616  3.79293377  0.29825025
# 53   2.6798669 -10.8723153   3.3954950 -5.64058166  3.99115123
# 54   3.1164984  -8.8486225   3.5621168 -1.47152717  1.81954270
                                                                              
ind30$cos2
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# 1  0.2802853100 0.1307305814 3.456274e-01 1.146941e-02 3.780603e-03
# 2  0.1208471327 0.2201270021 3.354947e-02 4.241655e-03 1.395109e-01
# 3  0.0014842892 0.2154320462 4.836104e-02 1.989303e-04 4.006609e-01
# 4  0.0512444500 0.2289955550 8.781985e-02 7.357391e-05 1.009774e-01
# 5  0.0508889475 0.0028509592 6.254589e-01 1.741741e-02 3.501350e-03
# 6  0.0017832551 0.0485295578 4.596674e-01 1.212826e-01 8.214480e-03
# 7  0.0305813370 0.1168654882 1.288932e-01 5.564524e-03 1.793405e-02
# 8  0.2239446334 0.2131842728 2.961006e-01 3.181336e-06 8.927603e-04
# 9  0.0361544639 0.2398913980 2.451083e-01 1.062684e-01 3.076387e-02
# 10 0.1570368843 0.2431002015 1.578514e-03 1.217096e-01 1.939854e-03
# 11 0.3628981982 0.0327254001 3.790239e-02 1.764738e-01 2.651351e-03
# 12 0.2098080616 0.1246182818 2.691315e-02 8.094266e-02 1.132998e-01
# 13 0.0033731330 0.0994510965 3.664912e-02 2.846071e-01 1.935092e-01
# 14 0.5206685907 0.0344060992 5.803273e-02 3.330754e-02 3.568641e-02
# 15 0.4061834501 0.0539820040 8.135178e-02 5.253260e-02 2.716958e-02
# 16 0.0097187932 0.0959660734 2.580605e-01 5.349186e-02 3.443658e-04
# 17 0.3787255897 0.0677236706 1.041700e-01 2.084190e-04 8.683180e-02
# 18 0.1257749040 0.1565680232 1.501457e-01 1.231383e-01 1.862781e-02
# 19 0.6457622147 0.0087058188 1.125160e-01 3.567947e-02 6.147605e-05
# 20 0.0539244766 0.0576604130 1.175474e-01 5.198875e-01 4.719804e-03
# 21 0.3320916967 0.0167093115 1.006777e-01 1.365163e-02 6.394976e-02
# 22 0.0011692005 0.1380723971 3.103835e-01 5.416335e-03 6.836389e-04
# 23 0.0365856509 0.0151981280 2.055819e-02 1.003986e-02 2.154519e-04
# 24 0.0999410050 0.1123036345 5.540756e-02 6.990493e-03 2.523777e-04
# 25 0.0018596524 0.0616715773 2.230907e-03 7.427046e-02 5.304192e-03
# 26 0.1175887039 0.0218536836 5.602501e-03 1.079476e-02 6.363689e-02
# 27 0.1241034676 0.0119592229 5.356048e-05 5.236812e-01 7.282960e-02
# 28 0.0027218893 0.0001243181 7.955733e-02 3.723445e-01 1.856385e-01
# 29 0.4051534256 0.0032847009 5.949880e-02 6.666720e-02 2.219358e-01
# 30 0.0445851840 0.0172506092 9.412323e-04 7.899541e-04 4.133903e-03
# 31 0.0242114971 0.0628336827 4.917730e-03 2.048066e-02 2.873465e-01
# 32 0.0047614308 0.0732719452 2.892891e-04 9.586577e-02 1.017079e-01
# 33 0.0879879763 0.0129942877 1.074791e-02 1.405147e-03 9.741008e-02
# 34 0.0314868692 0.0031886249 5.773639e-02 1.615029e-02 2.161274e-02
# 35 0.0009966107 0.0340005155 1.462318e-01 1.679749e-01 1.857105e-01
# 36 0.0817934964 0.1250539085 1.176834e-03 5.294348e-02 4.701352e-02
# 37 0.2037773370 0.0069985753 1.664300e-02 6.250462e-02 3.137881e-02
# 38 0.0226114282 0.0037791869 9.492394e-02 6.569040e-02 2.675222e-02
# 39 0.0336554097 0.0721775980 4.168395e-02 7.125524e-03 1.704752e-01
# 40 0.0066793094 0.4958048956 3.316619e-02 2.785565e-02 4.560043e-03
# 41 0.0515053142 0.2301434260 2.113382e-04 4.437841e-02 3.523280e-02
# 42 0.2623531383 0.3725416061 3.998270e-04 2.650805e-02 2.063908e-03
# 43 0.7710460696 0.0450356733 1.145377e-01 3.729441e-02 7.756953e-04
# 44 0.1985759164 0.3261221836 3.320183e-03 2.033285e-02 1.824973e-02
# 45 0.0352223727 0.0560612033 2.519980e-02 4.486685e-02 8.021977e-02
# 46 0.1179767972 0.3727717793 6.495005e-02 1.709256e-01 2.522674e-02
# 47 0.1664230194 0.0639938402 5.717670e-02 5.670508e-03 6.791373e-03
# 48 0.0072506563 0.3161918916 5.158198e-04 1.272135e-03 6.255862e-03
# 49 0.0454146622 0.2764537224 6.013611e-03 7.716193e-02 3.333675e-02
# 50 0.0698908733 0.4066133984 2.493502e-02 1.100850e-02 7.643833e-02
# 51 0.0518027289 0.0931708397 7.747747e-02 3.270963e-02 1.287086e-05
# 52 0.1109344286 0.3083326486 3.988012e-04 1.229000e-01 7.599114e-04
# 53 0.0259955883 0.4278753610 4.173298e-02 1.151651e-01 5.765933e-02
# 54 0.0467794099 0.3771136503 6.111351e-02 1.042936e-02 1.594577e-02
ind30$contrib
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# 1   4.224599347  3.379817099 11.935824932 5.007906e-01 2.788385e-01
# 2   0.580199586  1.812780408  0.369050442 5.899382e-02 3.277598e+00
# 3   0.011032714  2.746660319  0.823604191 4.283460e-03 1.457294e+01
# 4   0.212157460  1.626182706  0.833034916 8.823998e-04 2.045702e+00
# 5   0.485105229  0.046615925 13.660616295 4.809789e-01 1.633257e-01
# 6   0.010260404  0.478947771  6.059737428 2.021526e+00 2.312798e-01
# 7   0.106800663  0.700058888  1.031350801 5.629574e-02 3.064802e-01
# 8   2.169972799  3.543232151  6.573731423 8.930035e-05 4.233066e-02
# 9   0.256340782  2.917435543  3.981739743 2.182679e+00 1.067341e+00
# 10  2.087147409  5.542013598  0.048068327 4.686044e+00 1.261615e-01
# 11  2.585186476  0.399873703  0.618632462 3.641811e+00 9.242318e-02
# 12  1.457974547  1.485389484  0.428500553 1.629428e+00 3.852683e+00
# 13  0.024551592  1.241613948  0.611179456 6.000978e+00 6.892145e+00
# 14  3.099093982  0.351268834  0.791417145 5.743094e-01 1.039400e+00
# 15  3.925989997  0.894966507  1.801577966 1.470909e+00 1.285039e+00
# 16  0.045274868  0.766819931  2.754390413 7.218760e-01 7.850032e-03
# 17  3.626245148  1.112253586  2.285252078 5.780960e-03 4.068346e+00
# 18  0.668505331  1.427396949  1.828448593 1.895981e+00 4.844837e-01
# 19  8.468419750  0.195825925  3.380672219 1.355433e+00 3.944957e-03
# 20  0.634955239  1.164570370  3.171240884 1.773359e+01 2.719492e-01
# 21  2.424318918  0.209228626  1.683929866 2.886997e-01 2.284426e+00
# 22  0.007309148  1.480524044  4.445652466 9.808742e-02 2.091275e-02
# 23  0.222363256  0.158443034  0.286283416 1.767706e-01 6.407804e-03
# 24  0.380916385  0.734193557  0.483853689 7.718338e-02 4.706992e-03
# 25  0.014165216  0.805762964  0.038934267 1.638844e+00 1.977048e-01
# 26  0.712495391  0.227128652  0.077778120 1.894783e-01 1.886826e+00
# 27  1.157204207  0.191275729  0.001144273 1.414566e+01 3.323078e+00
# 28  0.013471217  0.001055362  0.902143489 5.338403e+00 4.495841e+00
# 29  2.691116137  0.037423084  0.905482755 1.282789e+00 7.213512e+00
# 30  0.172077109  0.114200443  0.008323163 8.832110e-03 7.807274e-02
# 31  0.107955754  0.480559848  0.050239796 2.645441e-01 6.269553e+00
# 32  0.016468091  0.434684670  0.002292434 9.605045e-01 1.721342e+00
# 33  0.508777223  0.128880622  0.142392643 2.353727e-02 2.756226e+00
# 34  0.119200034  0.020705282  0.500789915 1.771158e-01 4.003717e-01
# 35  0.011055840  0.646967600  3.716782337 5.398098e+00 1.008113e+01
# 36  0.357609175  0.937816707  0.011788650 6.705519e-01 1.005818e+00
# 37  1.073743871  0.063253603  0.200925734 9.540848e-01 8.090733e-01
# 38  0.101455509  0.029085519  0.975848606 8.538460e-01 5.873732e-01
# 39  0.113930551  0.419100865  0.323305604 6.987667e-02 2.823924e+00
# 40  0.040669752  5.178233018  0.462694715 4.913414e-01 1.358675e-01
# 41  0.382377791  2.930690848  0.003594824 9.544262e-01 1.279954e+00
# 42  2.359784417  5.747667853  0.008239810 6.907071e-01 9.084121e-02
# 43 45.782153092  4.586731954 15.582009583 6.414897e+00 2.253791e-01
# 44  2.219409404  6.252043673  0.085022106 6.583233e-01 9.980993e-01
# 45  0.133781259  0.365232908  0.219297081 4.936650e-01 1.490955e+00
# 46  0.904982431  4.904755780  1.141516280 3.798225e+00 9.469155e-01
# 47  0.743927796  0.490666465  0.585592564 7.342930e-02 1.485530e-01
# 48  0.027759773  2.076443599  0.004524761 1.410918e-02 1.172013e-01
# 49  0.218556891  2.282031701  0.066307527 1.075726e+00 7.850518e-01
# 50  0.559687564  5.585186633  0.457502974 2.553779e-01 2.995318e+00
# 51  0.253843704  0.783113267  0.869857781 4.643222e-01 3.086227e-04
# 52  0.646512982  3.082202923  0.005325083 2.074879e+00 2.167107e-02
# 53  0.357551574 10.094554018  1.315156635 4.588704e+00 3.880746e+00
# 54  0.483555215  6.686431507  1.447396786 3.123050e-01 8.065723e-01



fviz_pca_ind(DM_pca_FM3, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE  # avoid text overlapping, slow if many points
)
pca18_indno0 <- 
  fviz_pca_ind(DM_pca_FM3, col.ind = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE  # avoid text overlapping, slow if many points
  )
pca18_indno0
png("pca18_indno0")
print(pca18_indno0)
dev.off()
# change point size according to cos2
fviz_pca_ind(DM_pca_FM3, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_pointno0<- fviz_pca_ind(DM_pca_FM3, pointsize = "cos2", 
                                   pointshape = 21, fill = "#E7B800",
                                   repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_pointno0
png("pca18_ind_pointno0")
print(pca18_ind_pointno0)
dev.off()
# to change boht pointsize and colour according to cos2
fviz_pca_ind(DM_pca_FM3, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_sizeno0<- fviz_pca_ind(DM_pca_FM3, col.ind = "cos2", pointsize = "cos2",
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_sizeno0
png("pca18_ind_point_sizeno0")
print(pca18_ind_point_sizeno0)
dev.off()
# To visualize the contribution of individuals to the first two principal components
fviz_contrib(DM_pca_FM3, choice = "ind", axes = 1:2)
contr18_indno0<- fviz_contrib(DM_pca_FM3, choice = "ind", axes = 1:2)
contr18_indno0
png("contr18_indno0")
print(contr18_indno0)
dev.off()
#### Color by a custom continuous variable ####
# Spine age

fviz_pca_ind(DM_pca_FM3, col.ind = DM_data.frame3$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# trying to do the comparisson of age with the PCA loadins

# Get the scores of the first principal component
pc1_scores18no0 <- DM_pca_FM3$ind$coord[, "Dim.1"]
# second dim
pc1_scores18PC2no0 <- DM_pca_FM3$ind$coord[, "Dim.2"]

# Calculate the correlation
correlation18no0 <- cor(DM_data.frame3$Spine.Age, pc1_scores18no0)

correlation18no0
#[1] 0.06079824
# cortests
correlation18testno0 <- cor.test(DM_data.frame3$Spine.Age, pc1_scores18no0)

correlation18testno0


# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$Spine.Age and pc1_scores18PC2no0
# t = 2.277, df = 52, p-value = 0.02693
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.0362669 0.5264220
# sample estimates:
#       cor 
# 0.3011035 


correlation18PC2no0 <- cor(DM_data.frame3$Spine.Age, pc1_scores18PC2no0)

correlation18PC2no0
#[1] 0.3011035
correlation18PC2testno0 <- cor.test(DM_data.frame3$Spine.Age, pc1_scores18PC2no0)

correlation18PC2testno0
# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$Spine.Age and pc1_scores18PC2no0
# t = 2.277, df = 52, p-value = 0.02693
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.0362669 0.5264220
# sample estimates:
#       cor 
# 0.3011035 

# weight now
# Calculate the correlation
correlation18wno0 <- cor(DM_data.frame3$W, pc1_scores18no0)

correlation18wno0
#0.1600918
correlation18wtestno0 <- cor.test(DM_data.frame3$W, pc1_scores18no0)

correlation18wtestno0

# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$W and pc1_scores18no0
# t = 1.1695, df = 52, p-value = 0.2475
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1124908  0.4102658
# sample estimates:
#       cor 
# 0.1600918 

# Calculate the correlation
correlation18wPC2no0 <- cor(DM_data.frame3$W, pc1_scores18PC2no0)

correlation18wPC2no0
#[1] 0.2807467
correlation18wPC2testno0 <- cor.test(DM_data.frame3$W, pc1_scores18PC2no0)

correlation18wPC2testno0

# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$W and pc1_scores18PC2no0
# t = 2.1093, df = 52, p-value = 0.03975
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.01404175 0.51015730
# sample estimates:
#       cor 
# 0.2807467 


# fl now
# Calculate the correlation
cor(DM_data.frame3$FL, pc1_scores18no0)
correlation18flno0 <- cor(DM_data.frame3$FL, pc1_scores18no0)

correlation18flno0
#[1] 0.1109496
cor.test(DM_data.frame2$FL, pc1_scores18no43)

correlation18fltestno0 <- cor.test(DM_data.frame3$FL, pc1_scores18no0)

correlation18fltestno0

# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$FL and pc1_scores18no0
# t = 0.80504, df = 52, p-value = 0.4245
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1616121  0.3677841
# sample estimates:
#       cor 
# 0.1109496 

correlation18flPC2no0 <- cor(DM_data.frame3$FL, pc1_scores18PC2no0)

correlation18flPC2no0
#[1] 0.2117646

correlation18flPC2testno0 <- cor.test(DM_data.frame3$FL, pc1_scores18PC2no0)

correlation18flPC2testno0
# Pearson's product-moment correlation
# 
# data:  DM_data.frame3$FL and pc1_scores18PC2no0
# t = 1.5625, df = 52, p-value = 0.1242
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.05936187  0.45379406
# sample estimates:
#       cor 
# 0.2117646 




fviz_pca_ind(DM_pca_FM3, col.ind = DM_data.frame3$Spine.Age, axes = c(2,3),
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# W

fviz_pca_ind(DM_pca_FM3, col.ind = DM_data.frame3$W,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Weight") 
#### Anova sex and site ####

# Perform the Two-Way ANOVA

# Get the scores of the first principal component
pc1_scores18no0 <- DM_pca_FM3$ind$coord[, "Dim.1"]

# Perform the ANOVA
modelssno0 <- aov(pc1_scores18no0 ~ DM_data.frame3$Sex * DM_data.frame3$Site, data = DM_data.frame3)

# Print the summary of the model
summary(modelssno0)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# DM_data.frame3$Sex                      2  417.6  208.78   6.636 0.00298 **
#   DM_data.frame3$Site                     3  158.9   52.96   1.683 0.18404   
# DM_data.frame3$Sex:DM_data.frame3$Site  3   16.4    5.47   0.174 0.91352   
# Residuals                              45 1415.7   31.46                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# only site model anova
# Perform the ANOVA
modelsno0site <- aov(pc1_scores18no0 ~ DM_data.frame3$Site, data = DM_data.frame3)

# Print the summary of the model
summary(modelsno0site)
# Df Sum Sq Mean Sq F value Pr(>F)  
# DM_data.frame3$Site  3  326.3  108.77   3.233 0.0299 *
#   Residuals           50 1682.3   33.65                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsno0site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18no0 ~ DM_data.frame3$Site, data = DM_data.frame3)
# 
# $`DM_data.frame3$Site`
# diff        lwr        upr     p adj
# Matheson-Dauphin River  -4.5215003 -10.150323  1.1073222 0.1562532
# Red River-Dauphin River -4.1269698  -9.755792  1.5018527 0.2214056
# Sandy Bar-Dauphin River -7.1748653 -13.674470 -0.6752609 0.0251680
# Red River-Matheson       0.3945305  -5.234292  6.0233530 0.9976762
# Sandy Bar-Matheson      -2.6533650  -9.152969  3.8462393 0.7002641
# Sandy Bar-Red River     -3.0478955  -9.547500  3.4517088 0.6008705
library(report)
report(models)

# pc 2
# Perform the ANOVA
modelsPC2no0site <- aov(pc1_scores18PC2no0 ~ DM_data.frame3$Site, data = DM_data.frame3)

# Print the summary of the model
summary(modelsPC2no0site)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame3$Site  3  771.5  257.16   32.18 9.91e-12 ***
#   Residuals           50  399.5    7.99                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
report(modelsPC2no43site)
TukeyHSD(modelsPC2no0site)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2no0 ~ DM_data.frame3$Site, data = DM_data.frame3)
# 
# $`DM_data.frame3$Site`
# diff        lwr       upr     p adj
# Matheson-Dauphin River   5.218059  2.4749741  7.961143 0.0000356
# Red River-Dauphin River  9.648252  6.9051673 12.391336 0.0000000
# Sandy Bar-Dauphin River  8.048783  4.8813421 11.216225 0.0000001
# Red River-Matheson       4.430193  1.6871086  7.173278 0.0004590
# Sandy Bar-Matheson       2.830725 -0.3367166  5.998166 0.0951807
# Sandy Bar-Red River     -1.599468 -4.7669097  1.567973 0.5410910


# only sex model anova
# Perform the ANOVA
modelsexno0 <- aov(pc1_scores18no0 ~ DM_data.frame3$Sex, data = DM_data.frame3)

# Print the summary of the model
summary(modelsexno0)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# DM_data.frame3$Sex  2  417.6   208.8   6.692 0.00262 **
#   Residuals          51 1591.0    31.2                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1             

TukeyHSD(modelsexno0)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18no0 ~ DM_data.frame3$Sex, data = DM_data.frame3)
# 
# $`DM_data.frame3$Sex`
# diff        lwr        upr     p adj
# M-F  -4.412739 -8.9399881  0.1145098 0.0575351
# UN-F  4.842922 -0.4552465 10.1410897 0.0796463
# UN-M  9.255661  3.1015554 15.4097661 0.0018659

# pc2
# Perform the ANOVA
modelsexPC2no0 <- aov(pc1_scores18PC2no0 ~ DM_data.frame3$Sex, data = DM_data.frame3)

# Print the summary ofPC2 the model
summary(modelsexPC2no0)
# Df Sum Sq Mean Sq F value Pr(>F)  
# DM_data.frame3$Sex  2  192.5   96.24   5.016 0.0103 *
#   Residuals          51  978.5   19.19                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsexPC2no0)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2no0 ~ DM_data.frame3$Sex, data = DM_data.frame3)
# 
# $`DM_data.frame3$Sex`
# diff         lwr       upr     p adj
# M-F   3.727438   0.1770034  7.277872 0.0376191
# UN-F -2.215397  -6.3704140  1.939621 0.4088780
# UN-M -5.942834 -10.7691089 -1.116560 0.0123038

#### Color by groups Site ####
pca18_PC1_PC2no0 <- fviz_pca_ind(DM_pca_FM3,
                                  geom.ind = "point", # show points only (but not "text"),
                                  pointsize = 2.75,
                                  pointshape = 21,
                                  fill.ind = DM_data.frame3$Site, # color by groups
                                  col.ind = DM_data.frame3$Site,
                                  palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                  legend.title = "Site",
                                  mean.point = FALSE,  # removes group mean point
                                  title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2no0
png("pca18_PC1_PC2no0")
print(pca18_PC1_PC2no0)
dev.off()

pca18_PC1_PC2_meanno0 <- fviz_pca_ind(DM_pca_FM3,
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointsize = 2.75,
                                       pointshape = 21,
                                       fill.ind = DM_data.frame3$Site, # color by groups
                                       col.ind = DM_data.frame3$Site,
                                       palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                       legend.title = "Site",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_meanno0
png("pca18_PC1_PC2_meanno0")
print(pca18_PC1_PC2_meanno0)
dev.off()

pca18_PC2_PC3no0 <- fviz_pca_ind(DM_pca_FM3,
                                  axes = c(2,3),
                                  geom.ind = "point", # show points only (but not "text"),
                                  pointsize = 2.75,
                                  pointshape = 21,
                                  fill.ind = DM_data.frame3$Site, # color by groups
                                  col.ind = DM_data.frame3$Site,
                                  palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                  legend.title = "Site",
                                  mean.point = FALSE,  # removes group mean point
                                  title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3no0
png("pca18_PC2_PC3no0")
print(pca18_PC2_PC3no0)
dev.off()


pca18_PC2_PC3_meanno0 <- fviz_pca_ind(DM_pca_FM3,
                                       axes = c(2,3),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointshape = 21,
                                       pointsize = 2.75,
                                       fill.ind = DM_data.frame3$Site, # color by groups
                                       col.ind = DM_data.frame3$Site,
                                       palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                       legend.title = "Site",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3_meanno0
png("pca18_PC2_PC3_meanno0")
print(pca18_PC2_PC3_meanno0)
dev.off()

pca18_PC1_PC3_meanno0 <- fviz_pca_ind(DM_pca_FM3,
                                       axes = c(1,3),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointshape = 21,
                                       pointsize = 2.75,
                                       fill.ind = DM_data.frame3$Site, # color by groups
                                       col.ind = DM_data.frame3$Site,
                                       palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                       legend.title = "Site",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC3_meanno0
png("pca18_PC1_PC3_meanno0")
print(pca18_PC1_PC3_meanno0)
dev.off()


pca18_PC3_PC4_meanno0 <- fviz_pca_ind(DM_pca_FM3,
                                       axes = c(3,4),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointshape = 21,
                                       pointsize = 2.75,
                                       fill.ind = DM_data.frame3$Site, # color by groups
                                       col.ind = DM_data.frame3$Site,
                                       palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                       legend.title = "Site",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_meanno0
png("pca18_PC3_PC4_meanno0")
print(pca18_PC3_PC4_meanno0)
dev.off()


pca18_PC4_PC5_meanno0 <- fviz_pca_ind(DM_pca_FM3,
                                       axes = c(4,5),
                                       geom.ind = "point", # show points only (but not "text"),
                                       pointshape = 21,
                                       pointsize = 2.75,
                                       fill.ind = DM_data.frame3$Site, # color by groups
                                       col.ind = DM_data.frame3$Site,
                                       palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                       legend.title = "Site",
                                       title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_meanno0




#### biplot of individuals and variables ####
fviz_pca_biplot(DM_pca_FM3, repel = TRUE,
                col.var = "#2E9FDF",  # variables colour
                col.ind = "#696969"  # individuals colour
)


#### colour individuals by site ####
fviz_pca_biplot(DM_pca_FM3,
                geom.ind = "point",
                label = "var",
                col.ind = DM_data.frame3$Site,  # colour by groups
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                axes = c(1,2),
                legend.title = "Site",
                col.var = "black",
                repel = TRUE,
                mean.point = FALSE  # removes group mean point
)

# colour individuals by Site
fviz_pca_biplot(DM_pca_FM3, axes = c(1,2),
                geom.ind = "point",
                col.ind = DM_data.frame3$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 50),  # top contributing var
                legend.title = "Site"
)
pca18_ind_by_siteno0 <- fviz_pca_biplot(DM_pca_FM3, axes = c(1,2),
                                         geom.ind = "point",
                                         col.ind = DM_data.frame3$Site,  # colour by groups
                                         mean.point = FALSE,  # removes group mean point
                                         pointsize = 2.5,
                                         palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                         label = "var",
                                         alpha.var ="contrib",col.var = "black",
                                         repel = TRUE,
                                         select.var = list(contrib = 50),  # top contributing var
                                         legend.title = "Site"
)
pca18_ind_by_siteno0
png("pca18_ind_by_siteno0")

print(pca18_ind_by_siteno0)
dev.off()
# colour individuals by Site
fviz_pca_biplot(DM_pca_FM3, axes = c(2,3),
                geom.ind = "point",
                col.ind = DM_data.frame3$Site,  # colour by groups
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
fviz_pca_biplot(DM_pca_FM3, axes = c(3,4),
                geom.ind = "point",
                col.ind = DM_data.frame3$Site,  # colour by groups
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
fviz_pca_biplot(DM_pca_FM3, axes = c(4,5),
                geom.ind = "point",
                col.ind = DM_data.frame3$Site,  # colour by groups
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
pca18_PC1_PC2_sexno0 <- fviz_pca_ind(DM_pca_FM3,
                                      geom.ind = "point", # show points only (but not "text"),
                                      pointsize = 2.75,
                                      pointshape = 21,
                                      fill.ind = DM_data.frame3$Sex, # color by groups
                                      col.ind = DM_data.frame3$Sex,
                                      palette = c("#99d8c9", "#08519c","#E7B800"),
                                      legend.title = "Sex",
                                      mean.point = FALSE,  # removes group mean point
                                      title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_sexno0

pca18_PC1_PC2_mean_sexno0 <- fviz_pca_ind(DM_pca_FM3,
                                           geom.ind = "point", # show points only (but not "text"),
                                           pointsize = 2.75,
                                           pointshape = 21,
                                           fill.ind = DM_data.frame3$Sex, # color by groups
                                           col.ind = DM_data.frame3$Sex,
                                           palette = c("#99d8c9", "#08519c","#E7B800"),
                                           legend.title = "Sex",
                                           title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_mean_sexno0
png("pca18_PC1_PC2_mean_sexno0")
print(pca18_PC1_PC2_mean_sexno0)
dev.off()

pca18_PC3_PC2_mean_sexno0 <- fviz_pca_ind(DM_pca_FM3, axes = c(2,3),
                                           geom.ind = "point", # show points only (but not "text"),
                                           pointsize = 2.75,
                                           pointshape = 21,
                                           fill.ind = DM_data.frame3$Sex, # color by groups
                                           col.ind = DM_data.frame3$Sex,
                                           palette = c("#99d8c9", "#08519c","#E7B800"),
                                           legend.title = "Sex",
                                           title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC2_mean_sexno0
png("pca18_PC3_PC2_mean_sexno0")
print(pca18_PC3_PC2_mean_sexno0)
dev.off()

pca18_PC3_PC4_mean_sexno0 <- fviz_pca_ind(DM_pca_FM3, axes = c(3,4),
                                           geom.ind = "point", # show points only (but not "text"),
                                           pointsize = 2.75,
                                           pointshape = 21,
                                           fill.ind = DM_data.frame3$Sex, # color by groups
                                           col.ind = DM_data.frame3$Sex,
                                           palette = c("#99d8c9", "#08519c","#E7B800"),
                                           legend.title = "Sex",
                                           title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_mean_sexno0
png("pca18_PC3_PC4_mean_sexno0")
print(pca18_PC3_PC4_mean_sexno0)
dev.off()

pca18_PC4_PC5_mean_sexno0 <- fviz_pca_ind(DM_pca_FM3, axes = c(4,5),
                                           geom.ind = "point", # show points only (but not "text"),
                                           pointsize = 2.75,
                                           pointshape = 21,
                                           fill.ind = DM_data.frame3$Sex, # color by groups
                                           col.ind = DM_data.frame3$Sex,
                                           palette = c("#99d8c9", "#08519c","#E7B800"),
                                           legend.title = "Sex",
                                           title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_mean_sexno0

png("pca18_PC4_PC5_mean_sexno0")
print(pca18_PC4_PC5_mean_sexno0)
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
