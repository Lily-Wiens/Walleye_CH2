#------------------------------------------#
#              Walleye Project             #
#------------------------------------------#

setwd("C:/Users/user/Dropbox/PC/Desktop/Regression")

#### 2018 Metabolite Only PCA ####
#### set up ####
getwd()
setwd("C:/Users/user/Dropbox/PC/Desktop/Regression")
# setting for the desktop
setwd("C:/Users/user/Desktop/Regression")

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


#### Read in 2018 data total ####

# read data
# W_2018_MET <- read.csv("C:/Users/user/Desktop/Wallaye new stat/2018_MET_FM.csv", row.names=1)
# View(W_2018_MET)
 x2018_MET_FM <- read.csv("C:/Users/user/Dropbox/PC/Desktop/Regression/2018_MET_FM.csv", quote="")
   View(x2018_MET_FM)
   
# upload data again
X2018_MET_FM <- read.csv("C:/Users/user/Desktop/Regression/2018_MET_FM.csv")
View(X2018_MET_FM)   
raw_data<- X2018_MET_FM
# load in object with met_groups
met_groups <- readxl::read_excel("Class_2018.xlsx") %>% 
  print()
# met_groups <- met_groups[-1,]
# Read in the same raw data in Excel format like Matt the whole mettabolites
raw_data <- readxl::read_excel("Metabolites_NMR_MS_ALL_2018.xlsx") %>% 
  mutate(Sex = na_if(Sex, "NA"), Spine_Age = na_if(Spine_Age, "NA")) %>% 
  mutate(Sex = as.factor(Sex), Site = as.factor(Site), Mortality = as.factor(Mortality), ) # did not work for me
names(raw_data)

# as data.frame
DM_data.frame <- data.frame(raw_data)
names(DM_data.frame)
dim(DM_data.frame)
# 54 172

# dataframe with just met data columns
DM_data_pca <- DM_data.frame[,c(10:172)] %>%
print()

dim(DM_data_pca) 
#54 163

# load in object with met_groups
met_groups <- readxl::read_excel("Metabolite_groups_18.xlsx") %>% 
  print()
# A tibble: 163 x 2
# Category                 Metabolite                 
# <chr>                    <chr>                      
#   1 Nonessential_Amino_Acids Serine                     
# 2 Biogenic_Amines            Putrescine                 
# 3 Biogenic_Amines          trans_Hydroxyproline       
# 4 Nonessential_Amino_Acids Asparagine                 
# 5 Nonessential_Amino_Acids Glutamine                  
# 6 Biogenic_Amines          alpha_Aminoadipic_acid     
# 7 Biogenic_Amines          Methionine_sulfoxide       
# 8 Biogenic_Amines          Acetyl_ornithine           
# 9 Biogenic_Amines          Citrulline                 
# 10 Biogenic_Amines          Asymmetric_dimethylarginine
# ... with 153 more rows
summary(DM_data_pca)
# Serine         Putrescine    Trans_Hydroxyproline
# Min.   : 23.50   Min.   : 3.62   Min.   :  0.242     
# 1st Qu.: 71.67   1st Qu.:14.62   1st Qu.:  1.248     
# Median : 93.65   Median :19.50   Median :  4.195     
# Mean   :105.26   Mean   :21.18   Mean   : 17.229     
# 3rd Qu.:128.00   3rd Qu.:24.45   3rd Qu.: 26.800     
# Max.   :268.00   Max.   :65.20   Max.   :118.000     
# Asparagine        Glutamine     Alpha_Aminoadipic_Acid
# Min.   :  0.000   Min.   : 0.00   Min.   : 0.542        
# 1st Qu.:  0.000   1st Qu.: 8.46   1st Qu.: 1.508        
# Median :  1.087   Median :15.85   Median : 3.000        
# Mean   : 14.208   Mean   :19.80   Mean   : 4.456        
# 3rd Qu.: 23.400   3rd Qu.:27.43   3rd Qu.: 6.478        
# Max.   :130.000   Max.   :84.00   Max.   :24.200        
# Methionine_Sulfoxide Acetyl_Ornithine   Citrulline     
# Min.   : 1.030       Min.   :0.305    Min.   :  2.010  
# 1st Qu.: 5.930       1st Qu.:1.458    1st Qu.:  3.513  
# Median : 7.850       Median :1.820    Median : 14.900  
# Mean   : 8.646       Mean   :2.154    Mean   : 46.011  
# 3rd Qu.:11.175       3rd Qu.:2.743    3rd Qu.: 39.500  
# Max.   :20.900       Max.   :5.960    Max.   :294.000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine
# Min.   :0.0000              Min.   :0.2630        
# 1st Qu.:0.4420              1st Qu.:0.8555        
# Median :0.5735              Median :1.0950        
# Mean   :0.6091              Mean   :1.1694        
# 3rd Qu.:0.7943              3rd Qu.:1.4000        
# Max.   :1.9900              Max.   :3.2900        
# Tryptophan      Kynurenine        Ornithine     
# Min.   : 5.98   Min.   :0.00000   Min.   : 0.000  
# 1st Qu.:12.65   1st Qu.:0.00000   1st Qu.: 1.450  
# Median :17.55   Median :0.00735   Median : 4.065  
# Mean   :19.27   Mean   :0.28816   Mean   : 5.916  
# 3rd Qu.:22.70   3rd Qu.:0.31375   3rd Qu.: 7.015  
# Max.   :45.20   Max.   :2.19000   Max.   :29.500  
# Lysine          Spermidine        Spermine    
# Min.   :  0.820   Min.   : 3.230   Min.   :26.50  
# 1st Qu.:  4.755   1st Qu.: 9.643   1st Qu.:36.23  
# Median : 19.200   Median :13.750   Median :45.70  
# Mean   : 27.468   Mean   :15.433   Mean   :47.13  
# 3rd Qu.: 41.025   3rd Qu.:17.950   3rd Qu.:54.62  
# Max.   :110.000   Max.   :61.000   Max.   :82.70  
# Sarcosine      Methylhistidine  Beta_Hydroxybutyric_Acid
# Min.   : 3.460   Min.   : 2.380   Min.   : 0.736          
# 1st Qu.: 9.195   1st Qu.: 6.955   1st Qu.: 4.575          
# Median :15.500   Median :11.950   Median : 6.695          
# Mean   :20.583   Mean   :19.753   Mean   : 8.883          
# 3rd Qu.:26.575   3rd Qu.:29.800   3rd Qu.: 9.595          
# Max.   :90.800   Max.   :85.800   Max.   :50.700          
# Alpha_Ketoglutaric_Acid  Butyric_acid   Propionic_Acid 
# Min.   :20.80           Min.   :1.280   Min.   :25.40  
# 1st Qu.:55.27           1st Qu.:1.900   1st Qu.:34.67  
# Median :60.75           Median :2.445   Median :42.65  
# Mean   :58.89           Mean   :2.524   Mean   :44.71  
# 3rd Qu.:64.22           3rd Qu.:3.053   3rd Qu.:52.88  
# Max.   :71.20           Max.   :4.900   Max.   :85.80  
# Fumaric_Acid   Isobutyric_Acid Hippuric_Acid    
# Min.   :21.20   Min.   :0.890   Min.   :0.01110  
# 1st Qu.:40.75   1st Qu.:1.750   1st Qu.:0.03062  
# Median :53.10   Median :2.530   Median :0.03650  
# Mean   :52.42   Mean   :2.935   Mean   :0.03721  
# 3rd Qu.:62.77   3rd Qu.:3.715   3rd Qu.:0.04353  
# Max.   :99.40   Max.   :6.720   Max.   :0.07400  
# Methylmalonic_Acid   LYSOC14_0       LYSOC16_1     
# Min.   :0.02050    Min.   :1.852   Min.   : 4.323  
# 1st Qu.:0.04188    1st Qu.:3.036   1st Qu.: 8.038  
# Median :0.06485    Median :3.539   Median :10.398  
# Mean   :0.07671    Mean   :3.646   Mean   :10.843  
# 3rd Qu.:0.08380    3rd Qu.:4.294   3rd Qu.:12.441  
# Max.   :0.73400    Max.   :6.116   Max.   :26.612  
# LYSOC16_0       LYSOC17_0       LYSOC18_2     
# Min.   :18.97   Min.   :1.195   Min.   :0.3734  
# 1st Qu.:31.57   1st Qu.:2.106   1st Qu.:0.9103  
# Median :41.44   Median :2.608   Median :1.1622  
# Mean   :42.15   Mean   :2.771   Mean   :1.1990  
# 3rd Qu.:47.61   3rd Qu.:3.279   3rd Qu.:1.4203  
# Max.   :95.67   Max.   :6.248   Max.   :2.0876  
# LYSOC18_1        LYSOC18_0       LYSOC20_4     
# Min.   : 3.478   Min.   :1.389   Min.   : 0.985  
# 1st Qu.: 7.837   1st Qu.:2.954   1st Qu.: 3.347  
# Median : 9.475   Median :4.135   Median : 5.082  
# Mean   : 9.742   Mean   :3.979   Mean   : 5.546  
# 3rd Qu.:11.953   3rd Qu.:4.874   3rd Qu.: 6.721  
# Max.   :16.755   Max.   :7.102   Max.   :14.805  
# LYSOC20_3        LYSOC24_0       LYSOC26_1     
# Min.   : 6.953   Min.   :1.343   Min.   :0.7214  
# 1st Qu.: 7.989   1st Qu.:2.532   1st Qu.:1.2015  
# Median : 8.931   Median :2.888   Median :1.4402  
# Mean   : 9.007   Mean   :3.049   Mean   :1.4443  
# 3rd Qu.: 9.942   3rd Qu.:3.462   3rd Qu.:1.6796  
# Max.   :12.636   Max.   :5.966   Max.   :2.2097  
# LYSOC26_0        LYSOC28_1        LYSOC28_0     
# Min.   :0.8358   Min.   :0.9302   Min.   :0.5939  
# 1st Qu.:1.3633   1st Qu.:1.5009   1st Qu.:1.0157  
# Median :1.5753   Median :1.6607   Median :1.2035  
# Mean   :1.6407   Mean   :1.6879   Mean   :1.1967  
# 3rd Qu.:1.8186   3rd Qu.:1.8913   3rd Qu.:1.3107  
# Max.   :2.8986   Max.   :2.6944   Max.   :1.9853  
# X14_1SMOH        X16_1SM          X16_0SM      
# Min.   :1.762   Min.   : 3.842   Min.   : 2.759  
# 1st Qu.:3.288   1st Qu.: 6.063   1st Qu.: 5.203  
# Median :4.163   Median : 7.876   Median : 6.506  
# Mean   :4.271   Mean   : 8.229   Mean   : 7.182  
# 3rd Qu.:4.865   3rd Qu.: 9.940   3rd Qu.: 8.398  
# Max.   :8.506   Max.   :20.919   Max.   :17.205  
# X16_1SMOH         X18_1SM         PC32_2AA    
# Min.   : 2.322   Min.   :1.446   Min.   :15.19  
# 1st Qu.: 4.333   1st Qu.:2.401   1st Qu.:23.56  
# Median : 6.615   Median :2.797   Median :30.48  
# Mean   : 6.677   Mean   :3.026   Mean   :34.00  
# 3rd Qu.: 8.251   3rd Qu.:3.582   3rd Qu.:41.68  
# Max.   :14.324   Max.   :5.857   Max.   :70.59  
# X18_0SM          X20_2SM          PC36_0AE     
# Min.   : 6.595   Min.   : 3.103   Min.   : 4.782  
# 1st Qu.: 9.783   1st Qu.: 6.921   1st Qu.: 9.883  
# Median :12.134   Median : 9.034   Median :13.371  
# Mean   :13.229   Mean   : 9.354   Mean   :13.925  
# 3rd Qu.:15.823   3rd Qu.:11.081   3rd Qu.:17.455  
# Max.   :29.556   Max.   :21.746   Max.   :26.898  
# PC36_6AA         PC36_0AA        X22_2SMOH     
# Min.   : 33.39   Min.   : 18.22   Min.   : 3.798  
# 1st Qu.: 62.69   1st Qu.: 46.35   1st Qu.: 6.093  
# Median : 79.88   Median : 81.88   Median : 8.335  
# Mean   : 82.93   Mean   : 88.53   Mean   : 8.991  
# 3rd Qu.:100.71   3rd Qu.:127.87   3rd Qu.:11.318  
# Max.   :187.88   Max.   :209.69   Max.   :26.864  
# X22_1SMOH        PC38_6AA         PC38_0AA     
# Min.   :2.071   Min.   : 297.1   Min.   : 18.31  
# 1st Qu.:3.363   1st Qu.: 578.1   1st Qu.: 40.32  
# Median :4.269   Median : 700.6   Median : 55.16  
# Mean   :4.361   Mean   : 713.2   Mean   : 56.84  
# 3rd Qu.:5.151   3rd Qu.: 840.6   3rd Qu.: 71.04  
# Max.   :9.969   Max.   :1337.0   Max.   :139.60  
# PC40_6AE        X24_1SMOH         PC40_6AA     
# Min.   : 28.17   Min.   : 4.490   Min.   : 27.53  
# 1st Qu.: 56.93   1st Qu.: 9.051   1st Qu.: 79.26  
# Median : 69.31   Median :10.227   Median :114.12  
# Mean   : 72.31   Mean   :11.331   Mean   :121.08  
# 3rd Qu.: 83.23   3rd Qu.:14.213   3rd Qu.:159.76  
# Max.   :150.51   Max.   :21.248   Max.   :245.69  
# PC40_2AA        PC401AA            C2        
# Min.   :1.489   Min.   :1.169   Min.   : 3.505  
# 1st Qu.:2.374   1st Qu.:2.264   1st Qu.:11.411  
# Median :2.816   Median :2.794   Median :15.130  
# Mean   :2.795   Mean   :2.710   Mean   :16.520  
# 3rd Qu.:3.159   3rd Qu.:3.097   3rd Qu.:20.302  
# Max.   :4.721   Max.   :4.745   Max.   :37.195  
# C3_1               C3              C4_1        
# Min.   :0.01720   Min.   :0.4332   Min.   :0.03720  
# 1st Qu.:0.03857   1st Qu.:0.9799   1st Qu.:0.05042  
# Median :0.04515   Median :1.3726   Median :0.06380  
# Mean   :0.04956   Mean   :1.6131   Mean   :0.06774  
# 3rd Qu.:0.06237   3rd Qu.:1.9258   3rd Qu.:0.07852  
# Max.   :0.10650   Max.   :4.6204   Max.   :0.12330  
# C4              C3OH              C5_1        
# Min.   :0.0617   Min.   :0.01140   Min.   :0.00860  
# 1st Qu.:0.1359   1st Qu.:0.02010   1st Qu.:0.02022  
# Median :0.1817   Median :0.02455   Median :0.02240  
# Mean   :0.2239   Mean   :0.02446   Mean   :0.02362  
# 3rd Qu.:0.2545   3rd Qu.:0.02823   3rd Qu.:0.02792  
# Max.   :1.0766   Max.   :0.04350   Max.   :0.03890  
# C5               C4OH              C6_1       
# Min.   :0.03030   Min.   :0.02170   Min.   :0.1171  
# 1st Qu.:0.08697   1st Qu.:0.03560   1st Qu.:0.1603  
# Median :0.14810   Median :0.04455   Median :0.1818  
# Mean   :0.23941   Mean   :0.06464   Mean   :0.1914  
# 3rd Qu.:0.30443   3rd Qu.:0.06110   3rd Qu.:0.2150  
# Max.   :1.39790   Max.   :0.45600   Max.   :0.3187  
# C6              C5OH            C5_1DC       
# Min.   :0.5102   Min.   :0.0662   Min.   :0.06710  
# 1st Qu.:0.6616   1st Qu.:0.1496   1st Qu.:0.08605  
# Median :0.7370   Median :0.1722   Median :0.09720  
# Mean   :0.7994   Mean   :0.1807   Mean   :0.10748  
# 3rd Qu.:0.8328   3rd Qu.:0.2005   3rd Qu.:0.11790  
# Max.   :1.9323   Max.   :0.3586   Max.   :0.28850  
# C5DC              C8             C5MDC        
# Min.   :0.0735   Min.   :0.1043   Min.   :0.02280  
# 1st Qu.:0.1664   1st Qu.:0.1803   1st Qu.:0.02730  
# Median :0.1930   Median :0.2081   Median :0.02935  
# Mean   :0.2098   Mean   :0.2405   Mean   :0.03024  
# 3rd Qu.:0.2427   3rd Qu.:0.2644   3rd Qu.:0.03347  
# Max.   :0.5097   Max.   :0.8334   Max.   :0.04400  
# C9               C7DC             C10_2        
# Min.   :0.01490   Min.   :0.00610   Min.   :0.03560  
# 1st Qu.:0.01945   1st Qu.:0.01342   1st Qu.:0.06535  
# Median :0.02405   Median :0.01585   Median :0.07680  
# Mean   :0.02546   Mean   :0.03934   Mean   :0.08996  
# 3rd Qu.:0.02902   3rd Qu.:0.05097   3rd Qu.:0.10410  
# Max.   :0.07090   Max.   :0.17810   Max.   :0.21070  
# C10_1             C10             C12_1        
# Min.   :0.1439   Min.   :0.1260   Min.   :0.02490  
# 1st Qu.:0.1760   1st Qu.:0.1318   1st Qu.:0.04942  
# Median :0.1926   Median :0.1358   Median :0.05445  
# Mean   :0.1909   Mean   :0.1389   Mean   :0.05732  
# 3rd Qu.:0.2032   3rd Qu.:0.1424   3rd Qu.:0.06237  
# Max.   :0.2503   Max.   :0.1930   Max.   :0.17290  
# C12              C14_2             C14_1        
# Min.   :0.01780   Min.   :0.01000   Min.   :0.03590  
# 1st Qu.:0.02712   1st Qu.:0.01432   1st Qu.:0.06103  
# Median :0.03095   Median :0.01600   Median :0.08355  
# Mean   :0.03594   Mean   :0.01750   Mean   :0.09990  
# 3rd Qu.:0.03870   3rd Qu.:0.01847   3rd Qu.:0.10287  
# Max.   :0.13590   Max.   :0.05620   Max.   :0.56820  
# C14              C12DC            C14_2OH       
# Min.   :0.01090   Min.   :0.00860   Min.   :0.00730  
# 1st Qu.:0.02712   1st Qu.:0.01405   1st Qu.:0.01413  
# Median :0.03960   Median :0.01760   Median :0.02035  
# Mean   :0.06046   Mean   :0.01944   Mean   :0.02557  
# 3rd Qu.:0.06323   3rd Qu.:0.02170   3rd Qu.:0.02950  
# Max.   :0.56350   Max.   :0.05170   Max.   :0.13220  
# C14_1OH            C16_2             C16_1        
# Min.   :0.00960   Min.   :0.00820   Min.   :0.02410  
# 1st Qu.:0.01960   1st Qu.:0.01075   1st Qu.:0.07233  
# Median :0.02735   Median :0.01405   Median :0.11485  
# Mean   :0.03942   Mean   :0.01590   Mean   :0.16115  
# 3rd Qu.:0.04330   3rd Qu.:0.01812   3rd Qu.:0.17433  
# Max.   :0.35090   Max.   :0.05710   Max.   :1.09100  
# C16             C16_2OH           C16_1OH       
# Min.   :0.02010   Min.   :0.00840   Min.   :0.01230  
# 1st Qu.:0.07425   1st Qu.:0.01455   1st Qu.:0.02425  
# Median :0.10055   Median :0.02045   Median :0.03635  
# Mean   :0.13835   Mean   :0.02648   Mean   :0.04385  
# 3rd Qu.:0.16910   3rd Qu.:0.02933   3rd Qu.:0.04888  
# Max.   :1.02600   Max.   :0.17480   Max.   :0.30590  
# C16OH              C18_2             C18_1        
# Min.   :0.004600   Min.   :0.00790   Min.   :0.02510  
# 1st Qu.:0.007925   1st Qu.:0.01617   1st Qu.:0.08568  
# Median :0.010600   Median :0.02075   Median :0.13795  
# Mean   :0.012157   Mean   :0.02930   Mean   :0.19660  
# 3rd Qu.:0.013375   3rd Qu.:0.03243   3rd Qu.:0.24088  
# Max.   :0.064300   Max.   :0.18110   Max.   :1.35790  
# C18             C18_1OH        Arginine_Average
# Min.   :0.00920   Min.   :0.00770   Min.   : 15.01  
# 1st Qu.:0.02725   1st Qu.:0.01323   1st Qu.: 34.93  
# Median :0.03600   Median :0.01820   Median : 61.05  
# Mean   :0.04684   Mean   :0.02116   Mean   : 64.76  
# 3rd Qu.:0.05435   3rd Qu.:0.02448   3rd Qu.: 85.88  
# Max.   :0.27680   Max.   :0.10650   Max.   :181.06  
# Betaine_Average    C0_Average    Choline_Average 
# Min.   : 76.99   Min.   :13.79   Min.   : 24.64  
# 1st Qu.:145.78   1st Qu.:22.38   1st Qu.: 36.25  
# Median :189.81   Median :31.41   Median : 44.69  
# Mean   :200.16   Mean   :33.22   Mean   : 51.43  
# 3rd Qu.:259.12   3rd Qu.:39.31   3rd Qu.: 56.25  
# Max.   :501.12   Max.   :70.55   Max.   :203.31  
# Citrate_Average  Creatine_Average  Creatinine_Average
# Min.   : 51.34   Min.   :  46.09   Min.   : 0.5705   
# 1st Qu.:139.84   1st Qu.: 250.39   1st Qu.: 5.8881   
# Median :181.28   Median : 381.09   Median : 8.8575   
# Mean   :239.50   Mean   : 485.99   Mean   :10.3945   
# 3rd Qu.:292.59   3rd Qu.: 668.27   3rd Qu.:14.0687   
# Max.   :823.62   Max.   :1351.12   Max.   :24.0750   
# Glucose_Average   Glutamate_Average Glycine_Average
# Min.   :  737.4   Min.   : 26.80    Min.   :127.8  
# 1st Qu.: 3725.3   1st Qu.: 69.51    1st Qu.:332.8  
# Median : 5749.0   Median : 83.09    Median :450.9  
# Mean   : 5943.9   Mean   : 92.49    Mean   :464.1  
# 3rd Qu.: 7596.0   3rd Qu.:107.41    3rd Qu.:544.5  
# Max.   :13459.3   Max.   :247.38    Max.   :988.7  
# Histidine_Average Isoleucine_Average Lactate_Average
# Min.   :10.27     Min.   : 25.06     Min.   :1345   
# 1st Qu.:36.12     1st Qu.: 55.56     1st Qu.:4245   
# Median :46.17     Median : 80.43     Median :4864   
# Mean   :49.81     Mean   : 98.25     Mean   :4773   
# 3rd Qu.:65.72     3rd Qu.:115.84     3rd Qu.:5597   
# Max.   :92.21     Max.   :403.06     Max.   :6482   
# Leucine         Methionine     Phenylalanine_Average
# Min.   : 18.18   Min.   : 3.795   Min.   : 16.73       
# 1st Qu.: 71.55   1st Qu.:19.231   1st Qu.: 46.50       
# Median :109.88   Median :26.962   Median : 61.59       
# Mean   :129.27   Mean   :29.057   Mean   : 66.55       
# 3rd Qu.:161.44   3rd Qu.:31.406   3rd Qu.: 81.93       
# Max.   :554.75   Max.   :76.463   Max.   :172.88       
# Proline       Pyruvate_Average Succinate_Average
# Min.   : 10.19   Min.   : 47.36   Min.   : 3.587   
# 1st Qu.: 30.41   1st Qu.:122.91   1st Qu.: 7.726   
# Median : 50.23   Median :224.06   Median :10.181   
# Mean   : 69.92   Mean   :197.47   Mean   :11.230   
# 3rd Qu.:102.06   3rd Qu.:245.33   3rd Qu.:13.109   
# Max.   :249.31   Max.   :307.19   Max.   :24.900   
# Threonine_Average Tyrosine_Average Valine_Average  
# Min.   : 23.27    Min.   : 10.13   Min.   : 49.99  
# 1st Qu.: 69.42    1st Qu.: 36.11   1st Qu.: 91.88  
# Median :109.50    Median : 47.70   Median :142.88  
# Mean   :118.85    Mean   : 59.24   Mean   :176.41  
# 3rd Qu.:164.91    3rd Qu.: 61.94   3rd Qu.:220.81  
# Max.   :270.56    Max.   :381.19   Max.   :627.25  
# X1_Methylhistidine X2_Hydroxybutyrate
# Min.   : 1.00      Min.   : 2.250    
# 1st Qu.: 6.25      1st Qu.: 4.817    
# Median : 7.00      Median : 7.005    
# Mean   : 6.90      Mean   :11.132    
# 3rd Qu.: 7.47      3rd Qu.:14.595    
# Max.   :14.50      Max.   :52.130    
# X2_Hydroxyisovaleric_Acid X3_Hydroxybutyrate
# Min.   : 0.130            Min.   : 1.500    
# 1st Qu.: 0.880            1st Qu.: 4.660    
# Median : 1.130            Median : 6.380    
# Mean   : 1.760            Mean   : 8.886    
# 3rd Qu.: 1.938            3rd Qu.: 9.630    
# Max.   :11.500            Max.   :52.880    
# ADP              AMP             ATP         
# Min.   :  1.00   Min.   : 23.0   Min.   :  0.800  
# 1st Qu.: 12.00   1st Qu.:143.2   1st Qu.:  2.000  
# Median : 39.50   Median :179.0   Median :  2.100  
# Mean   : 55.00   Mean   :192.1   Mean   :  9.578  
# 3rd Qu.: 74.25   3rd Qu.:241.5   3rd Qu.:  3.050  
# Max.   :225.00   Max.   :382.0   Max.   :204.600  
# Acetamide         Acetate        Acetoacetate  
# Min.   : 1.800   Min.   : 31.00   Min.   : 1.50  
# 1st Qu.: 7.075   1st Qu.: 63.25   1st Qu.: 3.16  
# Median :11.000   Median : 80.00   Median : 5.88  
# Mean   :11.911   Mean   : 83.09   Mean   :11.46  
# 3rd Qu.:17.100   3rd Qu.: 92.50   3rd Qu.: 8.25  
# Max.   :25.100   Max.   :163.00   Max.   :86.25  
# Acetone        Adenosine        Alanine     
# Min.   :0.880   Min.   : 10.0   Min.   : 51.0  
# 1st Qu.:1.250   1st Qu.: 66.0   1st Qu.:180.2  
# Median :1.750   Median :132.0   Median :244.5  
# Mean   :1.949   Mean   :149.2   Mean   :272.1  
# 3rd Qu.:2.098   3rd Qu.:192.0   3rd Qu.:340.5  
# Max.   :7.130   Max.   :451.0   Max.   :857.0  
# Aspartate      Creatine_Phosphate    Cytidine    
# Min.   : 0.000   Min.   : 1.00      Min.   : 3.10  
# 1st Qu.: 7.062   1st Qu.: 3.10      1st Qu.: 8.90  
# Median :10.168   Median : 3.10      Median :11.90  
# Mean   :13.053   Mean   : 5.17      Mean   :12.21  
# 3rd Qu.:15.281   3rd Qu.: 4.50      3rd Qu.:15.32  
# Max.   :44.250   Max.   :28.40      Max.   :23.50  
# Dimethyl_Sulfone    Ethanol        Ethanolamine   
# Min.   : 3.800   Min.   : 1.000   Min.   : 24.00  
# 1st Qu.: 8.875   1st Qu.: 3.033   1st Qu.: 59.25  
# Median :11.650   Median : 3.750   Median : 90.00  
# Mean   :12.513   Mean   : 4.555   Mean   :104.91  
# 3rd Qu.:15.875   3rd Qu.: 4.250   3rd Qu.:137.75  
# Max.   :26.400   Max.   :24.750   Max.   :258.00  
# Formate         Glycerol      Guanidoacetate 
# Min.   :14.60   Min.   :  66.0   Min.   : 69.0  
# 1st Qu.:20.45   1st Qu.: 125.0   1st Qu.:121.5  
# Median :23.95   Median : 145.0   Median :167.0  
# Mean   :26.20   Mean   : 174.3   Mean   :164.9  
# 3rd Qu.:31.35   3rd Qu.: 177.8   3rd Qu.:201.8  
# Max.   :53.90   Max.   :1318.0   Max.   :280.0  
# Guanosine      Hypoxanthine        IMP        
# Min.   : 12.0   Min.   :11.10   Min.   : 21.90  
# 1st Qu.: 59.0   1st Qu.:18.45   1st Qu.: 60.00  
# Median :106.0   Median :26.65   Median : 86.85  
# Mean   :126.2   Mean   :30.64   Mean   :119.34  
# 3rd Qu.:182.0   3rd Qu.:40.70   3rd Qu.:155.10  
# Max.   :505.0   Max.   :73.30   Max.   :420.10  
# Inosine       Isopropanol       Malate      
# Min.   : 42.0   Min.   :0.50   Min.   :  3.10  
# 1st Qu.: 99.0   1st Qu.:1.00   1st Qu.:  8.80  
# Median :140.0   Median :1.13   Median :  9.40  
# Mean   :192.2   Mean   :1.16   Mean   : 16.73  
# 3rd Qu.:232.0   3rd Qu.:1.25   3rd Qu.: 16.40  
# Max.   :717.0   Max.   :1.75   Max.   :102.90  
# Malonate        Mannose         Methanol    
# Min.   :14.00   Min.   : 3.30   Min.   :104.0  
# 1st Qu.:30.95   1st Qu.:16.15   1st Qu.:126.0  
# Median :39.80   Median :25.85   Median :141.5  
# Mean   :42.90   Mean   :26.33   Mean   :153.4  
# 3rd Qu.:54.88   3rd Qu.:37.30   3rd Qu.:191.2  
# Max.   :92.40   Max.   :61.30   Max.   :225.0  
# N_N_Dimethylglycine  Nicotinurate   O_Acetylcarnitine
# Min.   : 7.10       Min.   : 5.30   Min.   : 3.10    
# 1st Qu.:14.78       1st Qu.:17.23   1st Qu.:11.32    
# Median :22.25       Median :28.10   Median :14.50    
# Mean   :22.95       Mean   :27.79   Mean   :16.21    
# 3rd Qu.:29.02       3rd Qu.:33.98   3rd Qu.:20.30    
# Max.   :48.10       Max.   :52.10   Max.   :33.00    
# O_Phosphocholine Propylene_Glycol    Taurine     
# Min.   :15.50    Min.   : 0.000   Min.   : 2787  
# 1st Qu.:36.67    1st Qu.: 4.300   1st Qu.: 5383  
# Median :45.05    Median : 6.100   Median : 6200  
# Mean   :44.76    Mean   : 8.444   Mean   : 6593  
# 3rd Qu.:50.75    3rd Qu.:10.450   3rd Qu.: 7956  
# Max.   :72.30    Max.   :50.400   Max.   :11645  
# Uridine         Xanthine     
# Min.   : 4.10   Min.   :0.0000  
# 1st Qu.:13.10   1st Qu.:0.0000  
# Median :17.35   Median :0.6875  
# Mean   :18.19   Mean   :0.7384  
# 3rd Qu.:22.45   3rd Qu.:1.3438  
# Max.   :36.80   Max.   :2.8750  
# sn_Glycero_3_Phosphocholine  Beta.Alanine   
# Min.   : 16.00              Min.   :  31.0  
# 1st Qu.: 38.25              1st Qu.: 132.5  
# Median : 54.50              Median : 209.5  
# Mean   : 72.33              Mean   : 240.7  
# 3rd Qu.: 90.25              3rd Qu.: 277.2  
# Max.   :355.00              Max.   :1065.0 

#### Correlations plot ####

DM_cor <- cor(DM_data_pca)

summary(DM_cor)
# Serine         Putrescine    Trans_Hydroxyproline
# Min.   : 23.50   Min.   : 3.62   Min.   :  0.242     
# 1st Qu.: 71.67   1st Qu.:14.62   1st Qu.:  1.248     
# Median : 93.65   Median :19.50   Median :  4.195     
# Mean   :105.26   Mean   :21.18   Mean   : 17.229     
# 3rd Qu.:128.00   3rd Qu.:24.45   3rd Qu.: 26.800     
# Max.   :268.00   Max.   :65.20   Max.   :118.000     
# Asparagine        Glutamine     Alpha_Aminoadipic_Acid
# Min.   :  0.000   Min.   : 0.00   Min.   : 0.542        
# 1st Qu.:  0.000   1st Qu.: 8.46   1st Qu.: 1.508        
# Median :  1.087   Median :15.85   Median : 3.000        
# Mean   : 14.208   Mean   :19.80   Mean   : 4.456        
# 3rd Qu.: 23.400   3rd Qu.:27.43   3rd Qu.: 6.478        
# Max.   :130.000   Max.   :84.00   Max.   :24.200        
# Methionine_Sulfoxide Acetyl_Ornithine   Citrulline     
# Min.   : 1.030       Min.   :0.305    Min.   :  2.010  
# 1st Qu.: 5.930       1st Qu.:1.458    1st Qu.:  3.513  
# Median : 7.850       Median :1.820    Median : 14.900  
# Mean   : 8.646       Mean   :2.154    Mean   : 46.011  
# 3rd Qu.:11.175       3rd Qu.:2.743    3rd Qu.: 39.500  
# Max.   :20.900       Max.   :5.960    Max.   :294.000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine
# Min.   :0.0000              Min.   :0.2630        
# 1st Qu.:0.4420              1st Qu.:0.8555        
# Median :0.5735              Median :1.0950        
# Mean   :0.6091              Mean   :1.1694        
# 3rd Qu.:0.7943              3rd Qu.:1.4000        
# Max.   :1.9900              Max.   :3.2900        
# Tryptophan      Kynurenine        Ornithine     
# Min.   : 5.98   Min.   :0.00000   Min.   : 0.000  
# 1st Qu.:12.65   1st Qu.:0.00000   1st Qu.: 1.450  
# Median :17.55   Median :0.00735   Median : 4.065  
# Mean   :19.27   Mean   :0.28816   Mean   : 5.916  
# 3rd Qu.:22.70   3rd Qu.:0.31375   3rd Qu.: 7.015  
# Max.   :45.20   Max.   :2.19000   Max.   :29.500  
# Lysine          Spermidine        Spermine    
# Min.   :  0.820   Min.   : 3.230   Min.   :26.50  
# 1st Qu.:  4.755   1st Qu.: 9.643   1st Qu.:36.23  
# Median : 19.200   Median :13.750   Median :45.70  
# Mean   : 27.468   Mean   :15.433   Mean   :47.13  
# 3rd Qu.: 41.025   3rd Qu.:17.950   3rd Qu.:54.62  
# Max.   :110.000   Max.   :61.000   Max.   :82.70  
# Sarcosine      Methylhistidine  Beta_Hydroxybutyric_Acid
# Min.   : 3.460   Min.   : 2.380   Min.   : 0.736          
# 1st Qu.: 9.195   1st Qu.: 6.955   1st Qu.: 4.575          
# Median :15.500   Median :11.950   Median : 6.695          
# Mean   :20.583   Mean   :19.753   Mean   : 8.883          
# 3rd Qu.:26.575   3rd Qu.:29.800   3rd Qu.: 9.595          
# Max.   :90.800   Max.   :85.800   Max.   :50.700          
# Alpha_Ketoglutaric_Acid  Butyric_acid   Propionic_Acid 
# Min.   :20.80           Min.   :1.280   Min.   :25.40  
# 1st Qu.:55.27           1st Qu.:1.900   1st Qu.:34.67  
# Median :60.75           Median :2.445   Median :42.65  
# Mean   :58.89           Mean   :2.524   Mean   :44.71  
# 3rd Qu.:64.22           3rd Qu.:3.053   3rd Qu.:52.88  
# Max.   :71.20           Max.   :4.900   Max.   :85.80  
# Fumaric_Acid   Isobutyric_Acid Hippuric_Acid    
# Min.   :21.20   Min.   :0.890   Min.   :0.01110  
# 1st Qu.:40.75   1st Qu.:1.750   1st Qu.:0.03062  
# Median :53.10   Median :2.530   Median :0.03650  
# Mean   :52.42   Mean   :2.935   Mean   :0.03721  
# 3rd Qu.:62.77   3rd Qu.:3.715   3rd Qu.:0.04353  
# Max.   :99.40   Max.   :6.720   Max.   :0.07400  
# Methylmalonic_Acid   LYSOC14_0       LYSOC16_1     
# Min.   :0.02050    Min.   :1.852   Min.   : 4.323  
# 1st Qu.:0.04188    1st Qu.:3.036   1st Qu.: 8.038  
# Median :0.06485    Median :3.539   Median :10.398  
# Mean   :0.07671    Mean   :3.646   Mean   :10.843  
# 3rd Qu.:0.08380    3rd Qu.:4.294   3rd Qu.:12.441  
# Max.   :0.73400    Max.   :6.116   Max.   :26.612  
# LYSOC16_0       LYSOC17_0       LYSOC18_2     
# Min.   :18.97   Min.   :1.195   Min.   :0.3734  
# 1st Qu.:31.57   1st Qu.:2.106   1st Qu.:0.9103  
# Median :41.44   Median :2.608   Median :1.1622  
# Mean   :42.15   Mean   :2.771   Mean   :1.1990  
# 3rd Qu.:47.61   3rd Qu.:3.279   3rd Qu.:1.4203  
# Max.   :95.67   Max.   :6.248   Max.   :2.0876  
# LYSOC18_1        LYSOC18_0       LYSOC20_4     
# Min.   : 3.478   Min.   :1.389   Min.   : 0.985  
# 1st Qu.: 7.837   1st Qu.:2.954   1st Qu.: 3.347  
# Median : 9.475   Median :4.135   Median : 5.082  
# Mean   : 9.742   Mean   :3.979   Mean   : 5.546  
# 3rd Qu.:11.953   3rd Qu.:4.874   3rd Qu.: 6.721  
# Max.   :16.755   Max.   :7.102   Max.   :14.805  
# LYSOC20_3        LYSOC24_0       LYSOC26_1     
# Min.   : 6.953   Min.   :1.343   Min.   :0.7214  
# 1st Qu.: 7.989   1st Qu.:2.532   1st Qu.:1.2015  
# Median : 8.931   Median :2.888   Median :1.4402  
# Mean   : 9.007   Mean   :3.049   Mean   :1.4443  
# 3rd Qu.: 9.942   3rd Qu.:3.462   3rd Qu.:1.6796  
# Max.   :12.636   Max.   :5.966   Max.   :2.2097  
# LYSOC26_0        LYSOC28_1        LYSOC28_0     
# Min.   :0.8358   Min.   :0.9302   Min.   :0.5939  
# 1st Qu.:1.3633   1st Qu.:1.5009   1st Qu.:1.0157  
# Median :1.5753   Median :1.6607   Median :1.2035  
# Mean   :1.6407   Mean   :1.6879   Mean   :1.1967  
# 3rd Qu.:1.8186   3rd Qu.:1.8913   3rd Qu.:1.3107  
# Max.   :2.8986   Max.   :2.6944   Max.   :1.9853  
# X14_1SMOH        X16_1SM          X16_0SM      
# Min.   :1.762   Min.   : 3.842   Min.   : 2.759  
# 1st Qu.:3.288   1st Qu.: 6.063   1st Qu.: 5.203  
# Median :4.163   Median : 7.876   Median : 6.506  
# Mean   :4.271   Mean   : 8.229   Mean   : 7.182  
# 3rd Qu.:4.865   3rd Qu.: 9.940   3rd Qu.: 8.398  
# Max.   :8.506   Max.   :20.919   Max.   :17.205  
# X16_1SMOH         X18_1SM         PC32_2AA    
# Min.   : 2.322   Min.   :1.446   Min.   :15.19  
# 1st Qu.: 4.333   1st Qu.:2.401   1st Qu.:23.56  
# Median : 6.615   Median :2.797   Median :30.48  
# Mean   : 6.677   Mean   :3.026   Mean   :34.00  
# 3rd Qu.: 8.251   3rd Qu.:3.582   3rd Qu.:41.68  
# Max.   :14.324   Max.   :5.857   Max.   :70.59  
# X18_0SM          X20_2SM          PC36_0AE     
# Min.   : 6.595   Min.   : 3.103   Min.   : 4.782  
# 1st Qu.: 9.783   1st Qu.: 6.921   1st Qu.: 9.883  
# Median :12.134   Median : 9.034   Median :13.371  
# Mean   :13.229   Mean   : 9.354   Mean   :13.925  
# 3rd Qu.:15.823   3rd Qu.:11.081   3rd Qu.:17.455  
# Max.   :29.556   Max.   :21.746   Max.   :26.898  
# PC36_6AA         PC36_0AA        X22_2SMOH     
# Min.   : 33.39   Min.   : 18.22   Min.   : 3.798  
# 1st Qu.: 62.69   1st Qu.: 46.35   1st Qu.: 6.093  
# Median : 79.88   Median : 81.88   Median : 8.335  
# Mean   : 82.93   Mean   : 88.53   Mean   : 8.991  
# 3rd Qu.:100.71   3rd Qu.:127.87   3rd Qu.:11.318  
# Max.   :187.88   Max.   :209.69   Max.   :26.864  
# X22_1SMOH        PC38_6AA         PC38_0AA     
# Min.   :2.071   Min.   : 297.1   Min.   : 18.31  
# 1st Qu.:3.363   1st Qu.: 578.1   1st Qu.: 40.32  
# Median :4.269   Median : 700.6   Median : 55.16  
# Mean   :4.361   Mean   : 713.2   Mean   : 56.84  
# 3rd Qu.:5.151   3rd Qu.: 840.6   3rd Qu.: 71.04  
# Max.   :9.969   Max.   :1337.0   Max.   :139.60  
# PC40_6AE        X24_1SMOH         PC40_6AA     
# Min.   : 28.17   Min.   : 4.490   Min.   : 27.53  
# 1st Qu.: 56.93   1st Qu.: 9.051   1st Qu.: 79.26  
# Median : 69.31   Median :10.227   Median :114.12  
# Mean   : 72.31   Mean   :11.331   Mean   :121.08  
# 3rd Qu.: 83.23   3rd Qu.:14.213   3rd Qu.:159.76  
# Max.   :150.51   Max.   :21.248   Max.   :245.69  
# PC40_2AA        PC401AA            C2        
# Min.   :1.489   Min.   :1.169   Min.   : 3.505  
# 1st Qu.:2.374   1st Qu.:2.264   1st Qu.:11.411  
# Median :2.816   Median :2.794   Median :15.130  
# Mean   :2.795   Mean   :2.710   Mean   :16.520  
# 3rd Qu.:3.159   3rd Qu.:3.097   3rd Qu.:20.302  
# Max.   :4.721   Max.   :4.745   Max.   :37.195  
# C3_1               C3              C4_1        
# Min.   :0.01720   Min.   :0.4332   Min.   :0.03720  
# 1st Qu.:0.03857   1st Qu.:0.9799   1st Qu.:0.05042  
# Median :0.04515   Median :1.3726   Median :0.06380  
# Mean   :0.04956   Mean   :1.6131   Mean   :0.06774  
# 3rd Qu.:0.06237   3rd Qu.:1.9258   3rd Qu.:0.07852  
# Max.   :0.10650   Max.   :4.6204   Max.   :0.12330  
# C4              C3OH              C5_1        
# Min.   :0.0617   Min.   :0.01140   Min.   :0.00860  
# 1st Qu.:0.1359   1st Qu.:0.02010   1st Qu.:0.02022  
# Median :0.1817   Median :0.02455   Median :0.02240  
# Mean   :0.2239   Mean   :0.02446   Mean   :0.02362  
# 3rd Qu.:0.2545   3rd Qu.:0.02823   3rd Qu.:0.02792  
# Max.   :1.0766   Max.   :0.04350   Max.   :0.03890  
# C5               C4OH              C6_1       
# Min.   :0.03030   Min.   :0.02170   Min.   :0.1171  
# 1st Qu.:0.08697   1st Qu.:0.03560   1st Qu.:0.1603  
# Median :0.14810   Median :0.04455   Median :0.1818  
# Mean   :0.23941   Mean   :0.06464   Mean   :0.1914  
# 3rd Qu.:0.30443   3rd Qu.:0.06110   3rd Qu.:0.2150  
# Max.   :1.39790   Max.   :0.45600   Max.   :0.3187  
# C6              C5OH            C5_1DC       
# Min.   :0.5102   Min.   :0.0662   Min.   :0.06710  
# 1st Qu.:0.6616   1st Qu.:0.1496   1st Qu.:0.08605  
# Median :0.7370   Median :0.1722   Median :0.09720  
# Mean   :0.7994   Mean   :0.1807   Mean   :0.10748  
# 3rd Qu.:0.8328   3rd Qu.:0.2005   3rd Qu.:0.11790  
# Max.   :1.9323   Max.   :0.3586   Max.   :0.28850  
# C5DC              C8             C5MDC        
# Min.   :0.0735   Min.   :0.1043   Min.   :0.02280  
# 1st Qu.:0.1664   1st Qu.:0.1803   1st Qu.:0.02730  
# Median :0.1930   Median :0.2081   Median :0.02935  
# Mean   :0.2098   Mean   :0.2405   Mean   :0.03024  
# 3rd Qu.:0.2427   3rd Qu.:0.2644   3rd Qu.:0.03347  
# Max.   :0.5097   Max.   :0.8334   Max.   :0.04400  
# C9               C7DC             C10_2        
# Min.   :0.01490   Min.   :0.00610   Min.   :0.03560  
# 1st Qu.:0.01945   1st Qu.:0.01342   1st Qu.:0.06535  
# Median :0.02405   Median :0.01585   Median :0.07680  
# Mean   :0.02546   Mean   :0.03934   Mean   :0.08996  
# 3rd Qu.:0.02902   3rd Qu.:0.05097   3rd Qu.:0.10410  
# Max.   :0.07090   Max.   :0.17810   Max.   :0.21070  
# C10_1             C10             C12_1        
# Min.   :0.1439   Min.   :0.1260   Min.   :0.02490  
# 1st Qu.:0.1760   1st Qu.:0.1318   1st Qu.:0.04942  
# Median :0.1926   Median :0.1358   Median :0.05445  
# Mean   :0.1909   Mean   :0.1389   Mean   :0.05732  
# 3rd Qu.:0.2032   3rd Qu.:0.1424   3rd Qu.:0.06237  
# Max.   :0.2503   Max.   :0.1930   Max.   :0.17290  
# C12              C14_2             C14_1        
# Min.   :0.01780   Min.   :0.01000   Min.   :0.03590  
# 1st Qu.:0.02712   1st Qu.:0.01432   1st Qu.:0.06103  
# Median :0.03095   Median :0.01600   Median :0.08355  
# Mean   :0.03594   Mean   :0.01750   Mean   :0.09990  
# 3rd Qu.:0.03870   3rd Qu.:0.01847   3rd Qu.:0.10287  
# Max.   :0.13590   Max.   :0.05620   Max.   :0.56820  
# C14              C12DC            C14_2OH       
# Min.   :0.01090   Min.   :0.00860   Min.   :0.00730  
# 1st Qu.:0.02712   1st Qu.:0.01405   1st Qu.:0.01413  
# Median :0.03960   Median :0.01760   Median :0.02035  
# Mean   :0.06046   Mean   :0.01944   Mean   :0.02557  
# 3rd Qu.:0.06323   3rd Qu.:0.02170   3rd Qu.:0.02950  
# Max.   :0.56350   Max.   :0.05170   Max.   :0.13220  
# C14_1OH            C16_2             C16_1        
# Min.   :0.00960   Min.   :0.00820   Min.   :0.02410  
# 1st Qu.:0.01960   1st Qu.:0.01075   1st Qu.:0.07233  
# Median :0.02735   Median :0.01405   Median :0.11485  
# Mean   :0.03942   Mean   :0.01590   Mean   :0.16115  
# 3rd Qu.:0.04330   3rd Qu.:0.01812   3rd Qu.:0.17433  
# Max.   :0.35090   Max.   :0.05710   Max.   :1.09100  
# C16             C16_2OH           C16_1OH       
# Min.   :0.02010   Min.   :0.00840   Min.   :0.01230  
# 1st Qu.:0.07425   1st Qu.:0.01455   1st Qu.:0.02425  
# Median :0.10055   Median :0.02045   Median :0.03635  
# Mean   :0.13835   Mean   :0.02648   Mean   :0.04385  
# 3rd Qu.:0.16910   3rd Qu.:0.02933   3rd Qu.:0.04888  
# Max.   :1.02600   Max.   :0.17480   Max.   :0.30590  
# C16OH              C18_2             C18_1        
# Min.   :0.004600   Min.   :0.00790   Min.   :0.02510  
# 1st Qu.:0.007925   1st Qu.:0.01617   1st Qu.:0.08568  
# Median :0.010600   Median :0.02075   Median :0.13795  
# Mean   :0.012157   Mean   :0.02930   Mean   :0.19660  
# 3rd Qu.:0.013375   3rd Qu.:0.03243   3rd Qu.:0.24088  
# Max.   :0.064300   Max.   :0.18110   Max.   :1.35790  
# C18             C18_1OH        Arginine_Average
# Min.   :0.00920   Min.   :0.00770   Min.   : 15.01  
# 1st Qu.:0.02725   1st Qu.:0.01323   1st Qu.: 34.93  
# Median :0.03600   Median :0.01820   Median : 61.05  
# Mean   :0.04684   Mean   :0.02116   Mean   : 64.76  
# 3rd Qu.:0.05435   3rd Qu.:0.02448   3rd Qu.: 85.88  
# Max.   :0.27680   Max.   :0.10650   Max.   :181.06  
# Betaine_Average    C0_Average    Choline_Average 
# Min.   : 76.99   Min.   :13.79   Min.   : 24.64  
# 1st Qu.:145.78   1st Qu.:22.38   1st Qu.: 36.25  
# Median :189.81   Median :31.41   Median : 44.69  
# Mean   :200.16   Mean   :33.22   Mean   : 51.43  
# 3rd Qu.:259.12   3rd Qu.:39.31   3rd Qu.: 56.25  
# Max.   :501.12   Max.   :70.55   Max.   :203.31  
# Citrate_Average  Creatine_Average  Creatinine_Average
# Min.   : 51.34   Min.   :  46.09   Min.   : 0.5705   
# 1st Qu.:139.84   1st Qu.: 250.39   1st Qu.: 5.8881   
# Median :181.28   Median : 381.09   Median : 8.8575   
# Mean   :239.50   Mean   : 485.99   Mean   :10.3945   
# 3rd Qu.:292.59   3rd Qu.: 668.27   3rd Qu.:14.0687   
# Max.   :823.62   Max.   :1351.12   Max.   :24.0750   
# Glucose_Average   Glutamate_Average Glycine_Average
# Min.   :  737.4   Min.   : 26.80    Min.   :127.8  
# 1st Qu.: 3725.3   1st Qu.: 69.51    1st Qu.:332.8  
# Median : 5749.0   Median : 83.09    Median :450.9  
# Mean   : 5943.9   Mean   : 92.49    Mean   :464.1  
# 3rd Qu.: 7596.0   3rd Qu.:107.41    3rd Qu.:544.5  
# Max.   :13459.3   Max.   :247.38    Max.   :988.7  
# Histidine_Average Isoleucine_Average Lactate_Average
# Min.   :10.27     Min.   : 25.06     Min.   :1345   
# 1st Qu.:36.12     1st Qu.: 55.56     1st Qu.:4245   
# Median :46.17     Median : 80.43     Median :4864   
# Mean   :49.81     Mean   : 98.25     Mean   :4773   
# 3rd Qu.:65.72     3rd Qu.:115.84     3rd Qu.:5597   
# Max.   :92.21     Max.   :403.06     Max.   :6482   
# Leucine         Methionine     Phenylalanine_Average
# Min.   : 18.18   Min.   : 3.795   Min.   : 16.73       
# 1st Qu.: 71.55   1st Qu.:19.231   1st Qu.: 46.50       
# Median :109.88   Median :26.962   Median : 61.59       
# Mean   :129.27   Mean   :29.057   Mean   : 66.55       
# 3rd Qu.:161.44   3rd Qu.:31.406   3rd Qu.: 81.93       
# Max.   :554.75   Max.   :76.463   Max.   :172.88       
# Proline       Pyruvate_Average Succinate_Average
# Min.   : 10.19   Min.   : 47.36   Min.   : 3.587   
# 1st Qu.: 30.41   1st Qu.:122.91   1st Qu.: 7.726   
# Median : 50.23   Median :224.06   Median :10.181   
# Mean   : 69.92   Mean   :197.47   Mean   :11.230   
# 3rd Qu.:102.06   3rd Qu.:245.33   3rd Qu.:13.109   
# Max.   :249.31   Max.   :307.19   Max.   :24.900   
# Threonine_Average Tyrosine_Average Valine_Average  
# Min.   : 23.27    Min.   : 10.13   Min.   : 49.99  
# 1st Qu.: 69.42    1st Qu.: 36.11   1st Qu.: 91.88  
# Median :109.50    Median : 47.70   Median :142.88  
# Mean   :118.85    Mean   : 59.24   Mean   :176.41  
# 3rd Qu.:164.91    3rd Qu.: 61.94   3rd Qu.:220.81  
# Max.   :270.56    Max.   :381.19   Max.   :627.25  
# X1_Methylhistidine X2_Hydroxybutyrate
# Min.   : 1.00      Min.   : 2.250    
# 1st Qu.: 6.25      1st Qu.: 4.817    
# Median : 7.00      Median : 7.005    
# Mean   : 6.90      Mean   :11.132    
# 3rd Qu.: 7.47      3rd Qu.:14.595    
# Max.   :14.50      Max.   :52.130    
# X2_Hydroxyisovaleric_Acid X3_Hydroxybutyrate
# Min.   : 0.130            Min.   : 1.500    
# 1st Qu.: 0.880            1st Qu.: 4.660    
# Median : 1.130            Median : 6.380    
# Mean   : 1.760            Mean   : 8.886    
# 3rd Qu.: 1.938            3rd Qu.: 9.630    
# Max.   :11.500            Max.   :52.880    
# ADP              AMP             ATP         
# Min.   :  1.00   Min.   : 23.0   Min.   :  0.800  
# 1st Qu.: 12.00   1st Qu.:143.2   1st Qu.:  2.000  
# Median : 39.50   Median :179.0   Median :  2.100  
# Mean   : 55.00   Mean   :192.1   Mean   :  9.578  
# 3rd Qu.: 74.25   3rd Qu.:241.5   3rd Qu.:  3.050  
# Max.   :225.00   Max.   :382.0   Max.   :204.600  
# Acetamide         Acetate        Acetoacetate  
# Min.   : 1.800   Min.   : 31.00   Min.   : 1.50  
# 1st Qu.: 7.075   1st Qu.: 63.25   1st Qu.: 3.16  
# Median :11.000   Median : 80.00   Median : 5.88  
# Mean   :11.911   Mean   : 83.09   Mean   :11.46  
# 3rd Qu.:17.100   3rd Qu.: 92.50   3rd Qu.: 8.25  
# Max.   :25.100   Max.   :163.00   Max.   :86.25  
# Acetone        Adenosine        Alanine     
# Min.   :0.880   Min.   : 10.0   Min.   : 51.0  
# 1st Qu.:1.250   1st Qu.: 66.0   1st Qu.:180.2  
# Median :1.750   Median :132.0   Median :244.5  
# Mean   :1.949   Mean   :149.2   Mean   :272.1  
# 3rd Qu.:2.098   3rd Qu.:192.0   3rd Qu.:340.5  
# Max.   :7.130   Max.   :451.0   Max.   :857.0  
# Aspartate      Creatine_Phosphate    Cytidine    
# Min.   : 0.000   Min.   : 1.00      Min.   : 3.10  
# 1st Qu.: 7.062   1st Qu.: 3.10      1st Qu.: 8.90  
# Median :10.168   Median : 3.10      Median :11.90  
# Mean   :13.053   Mean   : 5.17      Mean   :12.21  
# 3rd Qu.:15.281   3rd Qu.: 4.50      3rd Qu.:15.32  
# Max.   :44.250   Max.   :28.40      Max.   :23.50  
# Dimethyl_Sulfone    Ethanol        Ethanolamine   
# Min.   : 3.800   Min.   : 1.000   Min.   : 24.00  
# 1st Qu.: 8.875   1st Qu.: 3.033   1st Qu.: 59.25  
# Median :11.650   Median : 3.750   Median : 90.00  
# Mean   :12.513   Mean   : 4.555   Mean   :104.91  
# 3rd Qu.:15.875   3rd Qu.: 4.250   3rd Qu.:137.75  
# Max.   :26.400   Max.   :24.750   Max.   :258.00  
# Formate         Glycerol      Guanidoacetate 
# Min.   :14.60   Min.   :  66.0   Min.   : 69.0  
# 1st Qu.:20.45   1st Qu.: 125.0   1st Qu.:121.5  
# Median :23.95   Median : 145.0   Median :167.0  
# Mean   :26.20   Mean   : 174.3   Mean   :164.9  
# 3rd Qu.:31.35   3rd Qu.: 177.8   3rd Qu.:201.8  
# Max.   :53.90   Max.   :1318.0   Max.   :280.0  
# Guanosine      Hypoxanthine        IMP        
# Min.   : 12.0   Min.   :11.10   Min.   : 21.90  
# 1st Qu.: 59.0   1st Qu.:18.45   1st Qu.: 60.00  
# Median :106.0   Median :26.65   Median : 86.85  
# Mean   :126.2   Mean   :30.64   Mean   :119.34  
# 3rd Qu.:182.0   3rd Qu.:40.70   3rd Qu.:155.10  
# Max.   :505.0   Max.   :73.30   Max.   :420.10  
# Inosine       Isopropanol       Malate      
# Min.   : 42.0   Min.   :0.50   Min.   :  3.10  
# 1st Qu.: 99.0   1st Qu.:1.00   1st Qu.:  8.80  
# Median :140.0   Median :1.13   Median :  9.40  
# Mean   :192.2   Mean   :1.16   Mean   : 16.73  
# 3rd Qu.:232.0   3rd Qu.:1.25   3rd Qu.: 16.40  
# Max.   :717.0   Max.   :1.75   Max.   :102.90  
# Malonate        Mannose         Methanol    
# Min.   :14.00   Min.   : 3.30   Min.   :104.0  
# 1st Qu.:30.95   1st Qu.:16.15   1st Qu.:126.0  
# Median :39.80   Median :25.85   Median :141.5  
# Mean   :42.90   Mean   :26.33   Mean   :153.4  
# 3rd Qu.:54.88   3rd Qu.:37.30   3rd Qu.:191.2  
# Max.   :92.40   Max.   :61.30   Max.   :225.0  
# N_N_Dimethylglycine  Nicotinurate   O_Acetylcarnitine
# Min.   : 7.10       Min.   : 5.30   Min.   : 3.10    
# 1st Qu.:14.78       1st Qu.:17.23   1st Qu.:11.32    
# Median :22.25       Median :28.10   Median :14.50    
# Mean   :22.95       Mean   :27.79   Mean   :16.21    
# 3rd Qu.:29.02       3rd Qu.:33.98   3rd Qu.:20.30    
# Max.   :48.10       Max.   :52.10   Max.   :33.00    
# O_Phosphocholine Propylene_Glycol    Taurine     
# Min.   :15.50    Min.   : 0.000   Min.   : 2787  
# 1st Qu.:36.67    1st Qu.: 4.300   1st Qu.: 5383  
# Median :45.05    Median : 6.100   Median : 6200  
# Mean   :44.76    Mean   : 8.444   Mean   : 6593  
# 3rd Qu.:50.75    3rd Qu.:10.450   3rd Qu.: 7956  
# Max.   :72.30    Max.   :50.400   Max.   :11645  
# Uridine         Xanthine     
# Min.   : 4.10   Min.   :0.0000  
# 1st Qu.:13.10   1st Qu.:0.0000  
# Median :17.35   Median :0.6875  
# Mean   :18.19   Mean   :0.7384  
# 3rd Qu.:22.45   3rd Qu.:1.3438  
# Max.   :36.80   Max.   :2.8750  
# sn_Glycero_3_Phosphocholine  Beta.Alanine   
# Min.   : 16.00              Min.   :  31.0  
# 1st Qu.: 38.25              1st Qu.: 132.5  
# Median : 54.50              Median : 209.5  
# Mean   : 72.33              Mean   : 240.7  
# 3rd Qu.: 90.25              3rd Qu.: 277.2  
# Max.   :355.00              Max.   :1065.0  
# > DM_cor <- cor(DM_data_pca)
# > summary(DM_cor)
# Serine            Putrescine       Trans_Hydroxyproline
# Min.   :-0.396603   Min.   :-0.48154   Min.   :-0.51380    
# 1st Qu.:-0.006163   1st Qu.: 0.03771   1st Qu.:-0.13300    
# Median : 0.139882   Median : 0.21944   Median : 0.07324    
# Mean   : 0.174635   Mean   : 0.21264   Mean   : 0.11260    
# 3rd Qu.: 0.327823   3rd Qu.: 0.40671   3rd Qu.: 0.33782    
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000    
# Asparagine         Glutamine       
# Min.   :-0.50812   Min.   :-0.27134  
# 1st Qu.:-0.12739   1st Qu.: 0.06585  
# Median : 0.05558   Median : 0.17900  
# Mean   : 0.11526   Mean   : 0.20724  
# 3rd Qu.: 0.34433   3rd Qu.: 0.35419  
# Max.   : 1.00000   Max.   : 1.00000  
# Alpha_Aminoadipic_Acid Methionine_Sulfoxide
# Min.   :-0.424566      Min.   :-0.33609    
# 1st Qu.: 0.007831      1st Qu.:-0.03356    
# Median : 0.132037      Median : 0.07540    
# Mean   : 0.151367      Mean   : 0.14186    
# 3rd Qu.: 0.286182      3rd Qu.: 0.28971    
# Max.   : 1.000000      Max.   : 1.00000    
# Acetyl_Ornithine     Citrulline      
# Min.   :-0.37548   Min.   :-0.50286  
# 1st Qu.: 0.07867   1st Qu.:-0.22221  
# Median : 0.16655   Median :-0.10136  
# Mean   : 0.17893   Mean   :-0.07612  
# 3rd Qu.: 0.30573   3rd Qu.: 0.02742  
# Max.   : 1.00000   Max.   : 1.00000  
# Asymmetric_Dimethylarginine Total_Dimethylarginine
# Min.   :-0.50143            Min.   :-0.5480       
# 1st Qu.: 0.03624            1st Qu.: 0.1152       
# Median : 0.19510            Median : 0.2597       
# Mean   : 0.22665            Mean   : 0.2609       
# 3rd Qu.: 0.43079            3rd Qu.: 0.4545       
# Max.   : 1.00000            Max.   : 1.0000       
# Tryptophan         Kynurenine         Ornithine       
# Min.   :-0.57112   Min.   :-0.34335   Min.   :-0.44796  
# 1st Qu.: 0.07466   1st Qu.: 0.03163   1st Qu.:-0.09449  
# Median : 0.24224   Median : 0.13685   Median : 0.07353  
# Mean   : 0.24776   Mean   : 0.16280   Mean   : 0.09769  
# 3rd Qu.: 0.46676   3rd Qu.: 0.28962   3rd Qu.: 0.22129  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Lysine           Spermidine          Spermine      
# Min.   :-0.54111   Min.   :-0.54135   Min.   :-0.5488  
# 1st Qu.:-0.11734   1st Qu.: 0.02882   1st Qu.: 0.0448  
# Median : 0.03079   Median : 0.19817   Median : 0.1958  
# Mean   : 0.04501   Mean   : 0.20091   Mean   : 0.1952  
# 3rd Qu.: 0.18494   3rd Qu.: 0.39399   3rd Qu.: 0.3607  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# Sarcosine         Methylhistidine   
# Min.   :-0.300303   Min.   :-0.54561  
# 1st Qu.: 0.002097   1st Qu.:-0.15632  
# Median : 0.120844   Median :-0.00448  
# Mean   : 0.157637   Mean   :-0.02578  
# 3rd Qu.: 0.286426   3rd Qu.: 0.07567  
# Max.   : 1.000000   Max.   : 1.00000  
# Beta_Hydroxybutyric_Acid Alpha_Ketoglutaric_Acid
# Min.   :-0.33828         Min.   :-0.2785        
# 1st Qu.: 0.08495         1st Qu.: 0.1442        
# Median : 0.25511         Median : 0.2494        
# Mean   : 0.28070         Mean   : 0.2531        
# 3rd Qu.: 0.42832         3rd Qu.: 0.3742        
# Max.   : 1.00000         Max.   : 1.0000        
# Butyric_acid        Propionic_Acid      Fumaric_Acid     
# Min.   :-0.4877707   Min.   :-0.38071   Min.   :-0.38521  
# 1st Qu.:-0.0008598   1st Qu.:-0.10042   1st Qu.: 0.09324  
# Median : 0.2061905   Median :-0.02476   Median : 0.22269  
# Mean   : 0.1739917   Mean   : 0.02189   Mean   : 0.22043  
# 3rd Qu.: 0.3416129   3rd Qu.: 0.10467   3rd Qu.: 0.32468  
# Max.   : 1.0000000   Max.   : 1.00000   Max.   : 1.00000  
# Isobutyric_Acid   Hippuric_Acid      Methylmalonic_Acid
# Min.   :-0.5724   Min.   :-0.39532   Min.   :-0.28320  
# 1st Qu.:-0.0717   1st Qu.:-0.13829   1st Qu.:-0.02202  
# Median : 0.1018   Median : 0.02960   Median : 0.04834  
# Mean   : 0.1373   Mean   : 0.03868   Mean   : 0.06986  
# 3rd Qu.: 0.3796   3rd Qu.: 0.20836   3rd Qu.: 0.15688  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC14_0          LYSOC16_1          LYSOC16_0       
# Min.   :-0.30027   Min.   :-0.33400   Min.   :-0.36973  
# 1st Qu.:-0.01772   1st Qu.:-0.02663   1st Qu.:-0.11419  
# Median : 0.12570   Median : 0.08898   Median : 0.03441  
# Mean   : 0.15814   Mean   : 0.12690   Mean   : 0.05978  
# 3rd Qu.: 0.27669   3rd Qu.: 0.25812   3rd Qu.: 0.18663  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# LYSOC17_0          LYSOC18_2          LYSOC18_1      
# Min.   :-0.38989   Min.   :-0.43561   Min.   :-0.3924  
# 1st Qu.:-0.07968   1st Qu.: 0.09311   1st Qu.: 0.1290  
# Median : 0.07458   Median : 0.24604   Median : 0.2411  
# Mean   : 0.11871   Mean   : 0.24142   Mean   : 0.2426  
# 3rd Qu.: 0.27451   3rd Qu.: 0.33819   3rd Qu.: 0.3386  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# LYSOC18_0         LYSOC20_4          LYSOC20_3        
# Min.   :-0.3212   Min.   :-0.40107   Min.   :-0.254436  
# 1st Qu.: 0.1021   1st Qu.: 0.07662   1st Qu.:-0.009316  
# Median : 0.2304   Median : 0.23204   Median : 0.101552  
# Mean   : 0.2387   Mean   : 0.24812   Mean   : 0.123719  
# 3rd Qu.: 0.3667   3rd Qu.: 0.45883   3rd Qu.: 0.247010  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.000000  
# LYSOC24_0          LYSOC26_1           LYSOC26_0       
# Min.   :-0.34076   Min.   :-0.358866   Min.   :-0.40484  
# 1st Qu.:-0.05630   1st Qu.: 0.005739   1st Qu.:-0.07679  
# Median : 0.06371   Median : 0.193299   Median : 0.07291  
# Mean   : 0.10648   Mean   : 0.187086   Mean   : 0.09015  
# 3rd Qu.: 0.25089   3rd Qu.: 0.285774   3rd Qu.: 0.20623  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000  
# LYSOC28_1          LYSOC28_0          X14_1SMOH        
# Min.   :-0.25108   Min.   :-0.35578   Min.   :-0.328043  
# 1st Qu.: 0.06136   1st Qu.:-0.02143   1st Qu.:-0.005511  
# Median : 0.24580   Median : 0.11366   Median : 0.228263  
# Mean   : 0.23994   Mean   : 0.13597   Mean   : 0.223680  
# 3rd Qu.: 0.36731   3rd Qu.: 0.24074   3rd Qu.: 0.429325  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000  
# X16_1SM             X16_0SM           X16_1SMOH       
# Min.   :-0.551174   Min.   :-0.57470   Min.   :-0.54561  
# 1st Qu.:-0.005581   1st Qu.: 0.02794   1st Qu.:-0.05481  
# Median : 0.269637   Median : 0.25498   Median : 0.11559  
# Mean   : 0.252077   Mean   : 0.24725   Mean   : 0.13389  
# 3rd Qu.: 0.483511   3rd Qu.: 0.46148   3rd Qu.: 0.33476  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# X18_1SM           PC32_2AA             X18_0SM        
# Min.   :-0.3292   Min.   :-0.4175603   Min.   :-0.40054  
# 1st Qu.: 0.1520   1st Qu.: 0.0002709   1st Qu.: 0.03343  
# Median : 0.2864   Median : 0.1870406   Median : 0.20752  
# Mean   : 0.2829   Mean   : 0.1793133   Mean   : 0.21060  
# 3rd Qu.: 0.4361   3rd Qu.: 0.3386638   3rd Qu.: 0.39371  
# Max.   : 1.0000   Max.   : 1.0000000   Max.   : 1.00000  
# X20_2SM            PC36_0AE           PC36_6AA       
# Min.   :-0.41624   Min.   :-0.57279   Min.   :-0.33964  
# 1st Qu.: 0.08593   1st Qu.: 0.01679   1st Qu.: 0.07053  
# Median : 0.25066   Median : 0.20460   Median : 0.22783  
# Mean   : 0.25126   Mean   : 0.22146   Mean   : 0.24082  
# 3rd Qu.: 0.39089   3rd Qu.: 0.43377   3rd Qu.: 0.37998  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# PC36_0AA         X22_2SMOH          X22_1SMOH      
# Min.   :-0.6029   Min.   :-0.44440   Min.   :-0.3815  
# 1st Qu.:-0.1004   1st Qu.: 0.04966   1st Qu.: 0.1102  
# Median : 0.1190   Median : 0.26324   Median : 0.3304  
# Mean   : 0.1441   Mean   : 0.25531   Mean   : 0.2991  
# 3rd Qu.: 0.3898   3rd Qu.: 0.48274   3rd Qu.: 0.4765  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# PC38_6AA          PC38_0AA          PC40_6AE       
# Min.   :-0.2707   Min.   :-0.4575   Min.   :-0.26576  
# 1st Qu.: 0.0714   1st Qu.: 0.0871   1st Qu.: 0.08505  
# Median : 0.2547   Median : 0.2808   Median : 0.28352  
# Mean   : 0.2588   Mean   : 0.2653   Mean   : 0.26979  
# 3rd Qu.: 0.4148   3rd Qu.: 0.4690   3rd Qu.: 0.43664  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000  
# X24_1SMOH          PC40_6AA          PC40_2AA       
# Min.   :-0.3648   Min.   :-0.4514   Min.   :-0.32259  
# 1st Qu.: 0.0427   1st Qu.: 0.1007   1st Qu.: 0.03588  
# Median : 0.1946   Median : 0.2520   Median : 0.26493  
# Mean   : 0.2210   Mean   : 0.2508   Mean   : 0.24827  
# 3rd Qu.: 0.3910   3rd Qu.: 0.3842   3rd Qu.: 0.39821  
# Max.   : 1.0000   Max.   : 1.0000   Max.   : 1.00000  
# PC401AA               C2                C3_1         
# Min.   :-0.44826   Min.   :-0.23804   Min.   :-0.46355  
# 1st Qu.: 0.09145   1st Qu.: 0.01603   1st Qu.: 0.03684  
# Median : 0.24831   Median : 0.18373   Median : 0.19551  
# Mean   : 0.26961   Mean   : 0.18581   Mean   : 0.18212  
# 3rd Qu.: 0.42345   3rd Qu.: 0.34341   3rd Qu.: 0.32175  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C3                 C4_1                C4         
# Min.   :-0.603335   Min.   :-0.42795   Min.   :-0.5291  
# 1st Qu.: 0.003915   1st Qu.: 0.02058   1st Qu.: 0.0529  
# Median : 0.221548   Median : 0.19002   Median : 0.3150  
# Mean   : 0.220385   Mean   : 0.19351   Mean   : 0.2880  
# 3rd Qu.: 0.415687   3rd Qu.: 0.37834   3rd Qu.: 0.5459  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.0000  
# C3OH               C5_1                 C5           
# Min.   :-0.31143   Min.   :-0.492808   Min.   :-0.574241  
# 1st Qu.: 0.05143   1st Qu.: 0.008959   1st Qu.:-0.008348  
# Median : 0.15540   Median : 0.177397   Median : 0.155314  
# Mean   : 0.16063   Mean   : 0.170137   Mean   : 0.183100  
# 3rd Qu.: 0.26769   3rd Qu.: 0.343327   3rd Qu.: 0.405651  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.000000  
# C4OH               C6_1               C6          
# Min.   :-0.28991   Min.   :-0.5103   Min.   :-0.28909  
# 1st Qu.: 0.03888   1st Qu.:-0.0858   1st Qu.: 0.04023  
# Median : 0.19461   Median : 0.1397   Median : 0.16572  
# Mean   : 0.20302   Mean   : 0.1121   Mean   : 0.23978  
# 3rd Qu.: 0.36986   3rd Qu.: 0.3180   3rd Qu.: 0.41337  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# C5OH             C5_1DC              C5DC         
# Min.   :-0.2440   Min.   :-0.38970   Min.   :-0.38651  
# 1st Qu.: 0.1190   1st Qu.: 0.06537   1st Qu.: 0.07876  
# Median : 0.2152   Median : 0.25015   Median : 0.22489  
# Mean   : 0.2129   Mean   : 0.25201   Mean   : 0.21153  
# 3rd Qu.: 0.3056   3rd Qu.: 0.45627   3rd Qu.: 0.38309  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.00000  
# C8              C5MDC                C9         
# Min.   :-0.3118   Min.   :-0.23834   Min.   :-0.2787  
# 1st Qu.: 0.1024   1st Qu.: 0.03666   1st Qu.: 0.1012  
# Median : 0.2274   Median : 0.14344   Median : 0.2393  
# Mean   : 0.2768   Mean   : 0.17971   Mean   : 0.2855  
# 3rd Qu.: 0.3743   3rd Qu.: 0.26975   3rd Qu.: 0.3966  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# C7DC              C10_2             C10_1         
# Min.   :-0.30993   Min.   :-0.3677   Min.   :-0.25501  
# 1st Qu.:-0.04885   1st Qu.:-0.0359   1st Qu.: 0.06346  
# Median : 0.02918   Median : 0.1146   Median : 0.18970  
# Mean   : 0.03501   Mean   : 0.1491   Mean   : 0.15943  
# 3rd Qu.: 0.12035   3rd Qu.: 0.2857   3rd Qu.: 0.26975  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# C10               C12_1               C12          
# Min.   :-0.26610   Min.   :-0.22622   Min.   :-0.27893  
# 1st Qu.: 0.06103   1st Qu.: 0.06025   1st Qu.: 0.03148  
# Median : 0.18985   Median : 0.23565   Median : 0.16754  
# Mean   : 0.25195   Mean   : 0.28333   Mean   : 0.25179  
# 3rd Qu.: 0.33743   3rd Qu.: 0.46280   3rd Qu.: 0.34858  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C14_2              C14_1               C14          
# Min.   :-0.32128   Min.   :-0.29209   Min.   :-0.26319  
# 1st Qu.: 0.06089   1st Qu.: 0.09752   1st Qu.: 0.04636  
# Median : 0.20707   Median : 0.24215   Median : 0.17296  
# Mean   : 0.26948   Mean   : 0.30492   Mean   : 0.26607  
# 3rd Qu.: 0.42555   3rd Qu.: 0.45303   3rd Qu.: 0.41588  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C12DC             C14_2OH            C14_1OH        
# Min.   :-0.36832   Min.   :-0.32309   Min.   :-0.27917  
# 1st Qu.: 0.01588   1st Qu.: 0.03506   1st Qu.: 0.03648  
# Median : 0.19110   Median : 0.19597   Median : 0.16984  
# Mean   : 0.22947   Mean   : 0.25817   Mean   : 0.25990  
# 3rd Qu.: 0.40082   3rd Qu.: 0.35980   3rd Qu.: 0.40504  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C16_2              C16_1               C16           
# Min.   :-0.39651   Min.   :-0.32578   Min.   :-0.289144  
# 1st Qu.: 0.01767   1st Qu.: 0.02732   1st Qu.: 0.009065  
# Median : 0.18793   Median : 0.17531   Median : 0.154081  
# Mean   : 0.24393   Mean   : 0.25285   Mean   : 0.238623  
# 3rd Qu.: 0.35059   3rd Qu.: 0.34998   3rd Qu.: 0.396124  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000  
# C16_2OH            C16_1OH            C16OH         
# Min.   :-0.29537   Min.   :-0.2960   Min.   :-0.24956  
# 1st Qu.: 0.01661   1st Qu.: 0.0255   1st Qu.: 0.08498  
# Median : 0.16680   Median : 0.1717   Median : 0.19721  
# Mean   : 0.25071   Mean   : 0.2526   Mean   : 0.27391  
# 3rd Qu.: 0.40863   3rd Qu.: 0.4204   3rd Qu.: 0.43026  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# C18_2              C18_1               C18           
# Min.   :-0.36627   Min.   :-0.33909   Min.   :-0.302980  
# 1st Qu.: 0.05101   1st Qu.: 0.02472   1st Qu.:-0.003973  
# Median : 0.19542   Median : 0.18008   Median : 0.154035  
# Mean   : 0.27239   Mean   : 0.25440   Mean   : 0.230967  
# 3rd Qu.: 0.37469   3rd Qu.: 0.37925   3rd Qu.: 0.397986  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.000000  
# C18_1OH         Arginine_Average   Betaine_Average   
# Min.   :-0.26637   Min.   :-0.47774   Min.   :-0.42578  
# 1st Qu.: 0.02605   1st Qu.:-0.04264   1st Qu.:-0.04205  
# Median : 0.16331   Median : 0.08971   Median : 0.05008  
# Mean   : 0.25214   Mean   : 0.12282   Mean   : 0.06362  
# 3rd Qu.: 0.39494   3rd Qu.: 0.26828   3rd Qu.: 0.16729  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# C0_Average       Choline_Average    Citrate_Average   
# Min.   :-0.28460   Min.   :-0.39795   Min.   :-0.48572  
# 1st Qu.: 0.04185   1st Qu.:-0.05167   1st Qu.: 0.03184  
# Median : 0.21271   Median : 0.13122   Median : 0.29915  
# Mean   : 0.22246   Mean   : 0.12648   Mean   : 0.25141  
# 3rd Qu.: 0.35778   3rd Qu.: 0.29684   3rd Qu.: 0.43104  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Creatine_Average   Creatinine_Average Glucose_Average   
# Min.   :-0.40498   Min.   :-0.52157   Min.   :-0.36054  
# 1st Qu.: 0.02827   1st Qu.: 0.02099   1st Qu.:-0.09891  
# Median : 0.19950   Median : 0.18605   Median : 0.08709  
# Mean   : 0.18728   Mean   : 0.17837   Mean   : 0.08969  
# 3rd Qu.: 0.36500   3rd Qu.: 0.37629   3rd Qu.: 0.25713  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Glutamate_Average Glycine_Average    Histidine_Average
# Min.   :-0.2499   Min.   :-0.32847   Min.   :-0.2529  
# 1st Qu.: 0.0601   1st Qu.: 0.06695   1st Qu.: 0.0588  
# Median : 0.1844   Median : 0.19551   Median : 0.2294  
# Mean   : 0.2081   Mean   : 0.20918   Mean   : 0.2275  
# 3rd Qu.: 0.3460   3rd Qu.: 0.33828   3rd Qu.: 0.3649  
# Max.   : 1.0000   Max.   : 1.00000   Max.   : 1.0000  
# Isoleucine_Average Lactate_Average      Leucine        
# Min.   :-0.53638   Min.   :-0.4297   Min.   :-0.36218  
# 1st Qu.: 0.01795   1st Qu.: 0.0926   1st Qu.:-0.01434  
# Median : 0.16275   Median : 0.2846   Median : 0.10364  
# Mean   : 0.20191   Mean   : 0.2550   Mean   : 0.14390  
# 3rd Qu.: 0.41363   3rd Qu.: 0.4159   3rd Qu.: 0.28586  
# Max.   : 1.00000   Max.   : 1.0000   Max.   : 1.00000  
# Methionine       Phenylalanine_Average    Proline        
# Min.   :-0.38762   Min.   :-0.33167      Min.   :-0.51920  
# 1st Qu.: 0.02846   1st Qu.: 0.02357      1st Qu.:-0.07607  
# Median : 0.15289   Median : 0.15198      Median : 0.17222  
# Mean   : 0.20047   Mean   : 0.20706      Mean   : 0.17685  
# 3rd Qu.: 0.34624   3rd Qu.: 0.38704      3rd Qu.: 0.42984  
# Max.   : 1.00000   Max.   : 1.00000      Max.   : 1.00000  
# Pyruvate_Average   Succinate_Average  Threonine_Average
# Min.   :-0.39532   Min.   :-0.42636   Min.   :-0.5029  
# 1st Qu.:-0.02185   1st Qu.: 0.01573   1st Qu.: 0.1128  
# Median : 0.15524   Median : 0.15058   Median : 0.2083  
# Mean   : 0.12969   Mean   : 0.15277   Mean   : 0.2185  
# 3rd Qu.: 0.27369   3rd Qu.: 0.29864   3rd Qu.: 0.3394  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.0000  
# Tyrosine_Average   Valine_Average      X1_Methylhistidine
# Min.   :-0.27447   Min.   :-0.556598   Min.   :-0.42457  
# 1st Qu.: 0.01664   1st Qu.: 0.006041   1st Qu.:-0.13847  
# Median : 0.12447   Median : 0.176300   Median :-0.02398  
# Mean   : 0.15499   Mean   : 0.215124   Mean   :-0.03652  
# 3rd Qu.: 0.27043   3rd Qu.: 0.454843   3rd Qu.: 0.05178  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.00000  
# X2_Hydroxybutyrate X2_Hydroxyisovaleric_Acid
# Min.   :-0.35034   Min.   :-0.24252         
# 1st Qu.:-0.08738   1st Qu.:-0.05663         
# Median : 0.06050   Median : 0.08965         
# Mean   : 0.10261   Mean   : 0.06645         
# 3rd Qu.: 0.27864   3rd Qu.: 0.16514         
# Max.   : 1.00000   Max.   : 1.00000         
# X3_Hydroxybutyrate      ADP                  AMP          
# Min.   :-0.3587    Min.   :-0.5724437   Min.   :-0.30993  
# 1st Qu.: 0.0729    1st Qu.:-0.2989716   1st Qu.:-0.10441  
# Median : 0.2319    Median :-0.1396655   Median :-0.01404  
# Mean   : 0.2738    Mean   :-0.1322011   Mean   :-0.01396  
# 3rd Qu.: 0.4201    3rd Qu.:-0.0001898   3rd Qu.: 0.05305  
# Max.   : 1.0000    Max.   : 1.0000000   Max.   : 1.00000  
# ATP             Acetamide           Acetate        
# Min.   :-0.32877   Min.   :-0.60333   Min.   :-0.27396  
# 1st Qu.:-0.18336   1st Qu.:-0.29024   1st Qu.: 0.04896  
# Median :-0.10636   Median :-0.10808   Median : 0.17020  
# Mean   :-0.09172   Mean   :-0.12360   Mean   : 0.18946  
# 3rd Qu.:-0.03759   3rd Qu.: 0.01764   3rd Qu.: 0.29999  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Acetoacetate         Acetone           Adenosine       
# Min.   :-0.33964   Min.   :-0.36054   Min.   :-0.41978  
# 1st Qu.:-0.15105   1st Qu.: 0.03583   1st Qu.:-0.07066  
# Median :-0.06195   Median : 0.15129   Median : 0.08642  
# Mean   :-0.05225   Mean   : 0.20039   Mean   : 0.08441  
# 3rd Qu.: 0.01512   3rd Qu.: 0.31792   3rd Qu.: 0.21987  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Alanine           Aspartate        Creatine_Phosphate
# Min.   :-0.18010   Min.   :-0.43109   Min.   :-0.25800  
# 1st Qu.: 0.02577   1st Qu.:-0.10618   1st Qu.:-0.07546  
# Median : 0.10725   Median : 0.08306   Median : 0.02929  
# Mean   : 0.12015   Mean   : 0.10476   Mean   : 0.03239  
# 3rd Qu.: 0.20342   3rd Qu.: 0.29994   3rd Qu.: 0.11055  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Cytidine          Dimethyl_Sulfone      Ethanol         
# Min.   :-0.3814700   Min.   :-0.38909   Min.   :-0.385297  
# 1st Qu.:-0.1479331   1st Qu.:-0.04869   1st Qu.:-0.157202  
# Median : 0.0007744   Median : 0.09766   Median :-0.099691  
# Mean   : 0.0119877   Mean   : 0.12266   Mean   :-0.067424  
# 3rd Qu.: 0.1265201   3rd Qu.: 0.28027   3rd Qu.:-0.002905  
# Max.   : 1.0000000   Max.   : 1.00000   Max.   : 1.000000  
# Ethanolamine         Formate             Glycerol        
# Min.   :-0.54881   Min.   :-0.320546   Min.   :-0.177262  
# 1st Qu.:-0.21065   1st Qu.:-0.126181   1st Qu.:-0.034609  
# Median :-0.02502   Median :-0.018828   Median : 0.007125  
# Mean   :-0.02791   Mean   :-0.002905   Mean   : 0.028673  
# 3rd Qu.: 0.12761   3rd Qu.: 0.085817   3rd Qu.: 0.079076  
# Max.   : 1.00000   Max.   : 1.000000   Max.   : 1.000000  
# Guanidoacetate        Guanosine         Hypoxanthine     
# Min.   :-0.417646   Min.   :-0.41547   Min.   :-0.35127  
# 1st Qu.:-0.130011   1st Qu.:-0.15988   1st Qu.:-0.12348  
# Median :-0.003061   Median : 0.01526   Median : 0.03978  
# Mean   : 0.005801   Mean   : 0.01044   Mean   : 0.04833  
# 3rd Qu.: 0.107987   3rd Qu.: 0.15778   3rd Qu.: 0.20439  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# IMP              Inosine          Isopropanol      
# Min.   :-0.44915   Min.   :-0.41978   Min.   :-0.46007  
# 1st Qu.:-0.02848   1st Qu.:-0.17964   1st Qu.:-0.18776  
# Median : 0.08218   Median :-0.02175   Median :-0.04564  
# Mean   : 0.12032   Mean   : 0.02330   Mean   :-0.04955  
# 3rd Qu.: 0.29709   3rd Qu.: 0.23197   3rd Qu.: 0.04698  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Malate            Malonate           Mannose        
# Min.   :-0.18871   Min.   :-0.23515   Min.   :-0.47985  
# 1st Qu.:-0.02428   1st Qu.: 0.02416   1st Qu.:-0.08698  
# Median : 0.09753   Median : 0.11236   Median : 0.06373  
# Mean   : 0.09888   Mean   : 0.11376   Mean   : 0.06879  
# 3rd Qu.: 0.18090   3rd Qu.: 0.18941   3rd Qu.: 0.23204  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# Methanol         N_N_Dimethylglycine  Nicotinurate     
# Min.   :-0.427228   Min.   :-0.39154    Min.   :-0.28084  
# 1st Qu.:-0.117802   1st Qu.:-0.09314    1st Qu.: 0.06443  
# Median : 0.004117   Median : 0.04570    Median : 0.18492  
# Mean   :-0.002587   Mean   : 0.04066    Mean   : 0.18336  
# 3rd Qu.: 0.101671   3rd Qu.: 0.18359    3rd Qu.: 0.31093  
# Max.   : 1.000000   Max.   : 1.00000    Max.   : 1.00000  
# O_Acetylcarnitine   O_Phosphocholine   Propylene_Glycol  
# Min.   :-0.252116   Min.   :-0.30298   Min.   :-0.37081  
# 1st Qu.: 0.003015   1st Qu.:-0.09475   1st Qu.:-0.08792  
# Median : 0.169691   Median : 0.07910   Median : 0.07588  
# Mean   : 0.161285   Mean   : 0.08692   Mean   : 0.07489  
# 3rd Qu.: 0.306028   3rd Qu.: 0.24944   3rd Qu.: 0.22728  
# Max.   : 1.000000   Max.   : 1.00000   Max.   : 1.00000  
# Taurine            Uridine            Xanthine       
# Min.   :-0.23627   Min.   :-0.34058   Min.   :-0.44826  
# 1st Qu.:-0.05613   1st Qu.: 0.05309   1st Qu.:-0.27990  
# Median : 0.02144   Median : 0.18206   Median :-0.15767  
# Mean   : 0.05691   Mean   : 0.20194   Mean   :-0.11837  
# 3rd Qu.: 0.10979   3rd Qu.: 0.35532   3rd Qu.: 0.01406  
# Max.   : 1.00000   Max.   : 1.00000   Max.   : 1.00000  
# sn_Glycero_3_Phosphocholine  Beta.Alanine       
# Min.   :-0.380609           Min.   :-0.2724350  
# 1st Qu.:-0.189013           1st Qu.:-0.0775420  
# Median : 0.012453           Median :-0.0002643  
# Mean   : 0.001548           Mean   : 0.0105928  
# 3rd Qu.: 0.117228           3rd Qu.: 0.0806783  
# Max.   : 1.000000           Max.   : 1.0000000  

corrplot(DM_cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 60, 
         tl.cex = 1, cl.cex = 1)
plot_correlation18 <- corrplot(DM_cor, type = "upper", order = "hclust", 
                          tl.col = "black", tl.srt = 60, 
                          tl.cex = 1, cl.cex = 1)
plot_correlation18
# Open the graphics device
png("plot_correlation18.png")
dev.off()


#### PCA ####

# PCA with FactoMineR package

DM_pca_FM <- PCA(DM_data_pca, graph = FALSE)


#### Eigenvalues / Variances ####

# extract eigenvalues/varianes
get_eig(DM_pca_FM)
# eigenvalue variance.percent
# Dim.1  37.5128906      23.01404328
# Dim.2  22.3902305      13.73633773
# Dim.3  16.4728889      10.10606680
# Dim.4  12.9089060       7.91957426
# Dim.5   7.6355911       4.68441172
# Dim.6   6.3521377       3.89701697
# Dim.7   4.5350610       2.78224602
# Dim.8   3.8173049       2.34190486
# Dim.9   3.4950015       2.14417268
# Dim.10  3.4061504       2.08966284
# Dim.11  3.0736151       1.88565344
# Dim.12  3.0121307       1.84793296
# Dim.13  2.5562552       1.56825471
# Dim.14  2.3065133       1.41503882
# Dim.15  2.0707038       1.27037042
# Dim.16  1.9801633       1.21482410
# Dim.17  1.9335574       1.18623151
# Dim.18  1.8248909       1.11956494
# Dim.19  1.7714681       1.08679025
# Dim.20  1.6796392       1.03045352
# Dim.21  1.4706978       0.90226856
# Dim.22  1.4106331       0.86541906
# Dim.23  1.3538345       0.83057334
# Dim.24  1.2414541       0.76162829
# Dim.25  1.2219668       0.74967286
# Dim.26  1.1586548       0.71083116
# Dim.27  1.0962716       0.67255926
# Dim.28  1.0421076       0.63932980
# Dim.29  0.9861697       0.60501209
# Dim.30  0.9576700       0.58752762
# Dim.31  0.8792031       0.53938837
# Dim.32  0.8314265       0.51007761
# Dim.33  0.7490166       0.45951937
# Dim.34  0.7426690       0.45562513
# Dim.35  0.7015282       0.43038537
# Dim.36  0.6179649       0.37911960
# Dim.37  0.5854380       0.35916441
# Dim.38  0.5361584       0.32893153
# Dim.39  0.5010210       0.30737484
# Dim.40  0.4730543       0.29021736
# Dim.41  0.4480202       0.27485900
# Dim.42  0.4303495       0.26401810
# Dim.43  0.3829127       0.23491577
# Dim.44  0.3695006       0.22668751
# Dim.45  0.3568867       0.21894889
# Dim.46  0.3120355       0.19143282
# Dim.47  0.2865869       0.17582019
# Dim.48  0.2451416       0.15039360
# Dim.49  0.2307491       0.14156388
# Dim.50  0.1932637       0.11856668
# Dim.51  0.1714599       0.10519011
# Dim.52  0.1503885       0.09226286
# Dim.53  0.1306659       0.08016314
# cumulative.variance.percent
# Dim.1                     23.01404
# Dim.2                     36.75038
# Dim.3                     46.85645
# Dim.4                     54.77602
# Dim.5                     59.46043
# Dim.6                     63.35745
# Dim.7                     66.13970
# Dim.8                     68.48160
# Dim.9                     70.62577
# Dim.10                    72.71544
# Dim.11                    74.60109
# Dim.12                    76.44902
# Dim.13                    78.01728
# Dim.14                    79.43232
# Dim.15                    80.70269
# Dim.16                    81.91751
# Dim.17                    83.10374
# Dim.18                    84.22331
# Dim.19                    85.31010
# Dim.20                    86.34055
# Dim.21                    87.24282
# Dim.22                    88.10824
# Dim.23                    88.93881
# Dim.24                    89.70044
# Dim.25                    90.45011
# Dim.26                    91.16095
# Dim.27                    91.83350
# Dim.28                    92.47283
# Dim.29                    93.07785
# Dim.30                    93.66537
# Dim.31                    94.20476
# Dim.32                    94.71484
# Dim.33                    95.17436
# Dim.34                    95.62998
# Dim.35                    96.06037
# Dim.36                    96.43949
# Dim.37                    96.79865
# Dim.38                    97.12759
# Dim.39                    97.43496
# Dim.40                    97.72518
# Dim.41                    98.00004
# Dim.42                    98.26405
# Dim.43                    98.49897
# Dim.44                    98.72566
# Dim.45                    98.94461
# Dim.46                    99.13604
# Dim.47                    99.31186
# Dim.48                    99.46225
# Dim.49                    99.60382
# Dim.50                    99.72238
# Dim.51                    99.82757
# Dim.52                    99.91984
# Dim.53                   100.00000                  

fviz_screeplot(DM_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 16,
               ncp = Inf
)
screeplot18_all <- fviz_screeplot(DM_pca_FM, addlabels = TRUE, 
                                  ggtheme = theme_classic(),
                                  main = "",
                                  font.x = c(14, "bold"), font.y = c(14, "bold"),
                                  font.tickslab = 12,
                                  barfill = "#99d8c9", barcolor = "#66c2a4",
                                  font.submain = 16,
                                  ncp = Inf
)
screeplot18_all

# Open the graphics device
png("screeplot18_all")
dev.off()

fviz_screeplot(DM_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 10,
               ncp = 16
)

screeplot18_partial<- fviz_screeplot(DM_pca_FM, addlabels = TRUE, 
               ggtheme = theme_classic(),
               main = "",
               font.x = c(14, "bold"), font.y = c(14, "bold"),
               font.tickslab = 12,
               barfill = "#99d8c9", barcolor = "#66c2a4",
               font.submain = 10,
               ncp = 16
)
screeplot18_partial


# Open the graphics device
png("screeplot18_partial")
dev.off()

# extract the results for variables
var <- get_pca_var(DM_pca_FM)
var
# Principal Component Analysis Results for variables
#   Name       Description                                    
# 1 "$coord"   "Coordinates for the variables"                
# 2 "$cor"     "Correlations between variables and dimensions"
# 3 "$cos2"    "Cos2 for the variables"                       
# 4 "$contrib" "contributions of the variables"
head(var$coord)  # coordinates of variables
head(var$contrib) # contributions of variables to the PCs
head(var$cos2)  # Cos2: quality on the factore map
var$coord
# Dim.1         Dim.2        Dim.3
# Serine                       0.3388464093  2.038905e-01  0.824648578
# Putrescine                   0.5302696188 -4.867777e-01  0.090331769
# Trans_Hydroxyproline         0.3059433757 -7.766448e-01 -0.064470213
# Asparagine                   0.3038074802 -7.418040e-01 -0.025614733
# Glutamine                    0.4484044329 -1.718089e-01  0.616048692
# Alpha_Aminoadipic_Acid       0.3586028289 -3.194118e-01  0.272574403
# Methionine_Sulfoxide         0.2505179712 -1.536242e-01  0.728665703
# Acetyl_Ornithine             0.4453616497 -2.666102e-01  0.147391232
# Citrulline                  -0.1863719675  2.112130e-01 -0.335280461
# Asymmetric_Dimethylarginine  0.5437495169 -5.571181e-01  0.312053591
# Total_Dimethylarginine       0.6430978503 -4.392563e-01  0.391183820
# Tryptophan                   0.6308594242 -5.345513e-01  0.208028417
# Kynurenine                   0.3807518377 -1.890290e-01  0.442630248
# Ornithine                    0.1606198384  2.539382e-01  0.791607410
# Lysine                       0.0311485661  2.871536e-01  0.747960076
# Spermidine                   0.5476740505 -5.284207e-01 -0.085688110
# Spermine                     0.5543752384 -4.343337e-01 -0.087263214
# Sarcosine                    0.3071636897  1.683414e-01  0.659074987
# Methylhistidine             -0.0868544020  2.759798e-01  0.307504930
# Beta_Hydroxybutyric_Acid     0.7470116317  1.406963e-01 -0.237891251
# Alpha_Ketoglutaric_Acid      0.5566267081 -1.371849e-01  0.527120130
# Butyric_acid                 0.4957848677 -3.393116e-01 -0.299948815
# Propionic_Acid               0.0444735953  7.733669e-02 -0.009665286
# Fumaric_Acid                 0.4870625475 -9.659349e-03  0.493918257
# Isobutyric_Acid              0.4195623138 -6.779555e-01 -0.320503231
# Hippuric_Acid                0.0895676894 -4.651227e-01 -0.171968448
# Methylmalonic_Acid           0.1913479735 -2.087259e-01  0.049623722
# LYSOC14_0                    0.3814819058  1.589151e-01 -0.021284814
# LYSOC16_1                    0.3186804838  1.800642e-01 -0.093056193
# LYSOC16_0                    0.1385309916  3.142170e-01 -0.125292346
# LYSOC17_0                    0.2976529997  3.398759e-01 -0.145726096
# LYSOC18_2                    0.5945654628  1.128200e-01  0.210015321
# LYSOC18_1                    0.5727274058  8.071487e-02  0.315049500
# LYSOC18_0                    0.5311310650  5.994161e-02  0.523806460
# LYSOC20_4                    0.6476245787 -4.498801e-01 -0.049741892
# LYSOC20_3                    0.2894596547  4.290819e-02  0.195481346
# LYSOC24_0                    0.2446404778  4.292578e-01 -0.052179017
# LYSOC26_1                    0.4622372784  1.834427e-01 -0.010242433
# LYSOC26_0                    0.2010356582  3.956970e-01 -0.031045679
# LYSOC28_1                    0.6145798674  4.544715e-02 -0.152298150
# LYSOC28_0                    0.3326020824  2.289852e-01 -0.065632199
# X14_1SMOH                    0.6331417540 -9.701015e-02 -0.413795331
# X16_1SM                      0.7169647395 -4.480244e-01 -0.367418617
# X16_0SM                      0.6958427452 -4.813834e-01 -0.361224679
# X16_1SMOH                    0.4110879892 -5.284521e-01 -0.444593056
# X18_1SM                      0.6867497936 -8.883087e-02  0.342036204
# PC32_2AA                     0.3971957044  4.165011e-01  0.500908990
# X18_0SM                      0.4730653210  3.921227e-01  0.525803673
# X20_2SM                      0.5921260610  1.897404e-01  0.461854488
# PC36_0AE                     0.6091197098 -4.983897e-01 -0.261971389
# PC36_6AA                     0.6063220264 -9.804504e-03  0.091887948
# PC36_0AA                     0.4114987048 -6.637201e-01 -0.283604237
# X22_2SMOH                    0.7288230470 -2.944640e-01 -0.449813354
# X22_1SMOH                    0.8106063823 -1.875888e-01 -0.219155232
# PC38_6AA                     0.6636876851  6.535763e-02  0.027454235
# PC38_0AA                     0.7342564351 -3.167873e-01 -0.282470870
# PC40_6AE                     0.7098100615  7.015405e-02 -0.039285266
# X24_1SMOH                    0.5740478623  2.672553e-01  0.011066246
# PC40_6AA                     0.5884018537  1.014887e-01  0.480341707
# PC40_2AA                     0.6496039482  5.880178e-02 -0.136261557
# PC401AA                      0.6653513756  2.376208e-01  0.113466235
# C2                           0.4991459432  2.809680e-01  0.021562374
# C3_1                         0.4944874433 -2.776340e-01 -0.211485919
# C3                           0.5944766069 -5.784021e-01  0.004049056
# C4_1                         0.5373410343 -3.947879e-01 -0.208670988
# C4                           0.7916022675 -2.919006e-01 -0.244733071
# C3OH                         0.4074218557 -1.191333e-01 -0.021489683
# C5_1                         0.4645470086 -5.082754e-01 -0.092734528
# C5                           0.4897598110 -6.207907e-01 -0.062248263
# C4OH                         0.5660167619  3.009566e-02 -0.268775723
# C6_1                         0.3759414722 -5.813131e-02 -0.613038119
# C6                           0.6618472401  3.400685e-01 -0.340300570
# C5OH                         0.4783318633 -5.885792e-02  0.283107790
# C5_1DC                       0.6906948939 -1.855131e-01 -0.256520439
# C5DC                         0.5788484571 -8.183641e-02 -0.151990145
# C8                           0.7329648128  4.279972e-01 -0.085923217
# C5MDC                        0.4496815724  3.833098e-01  0.049138835
# C9                           0.7419645702  4.782717e-01 -0.082678441
# C7DC                         0.0899586554  1.471902e-01  0.034311336
# C10_2                        0.3897866249  6.146910e-01 -0.012686157
# C10_1                        0.4101403219 -7.870300e-02 -0.199357273
# C10                          0.6603984453  5.348487e-01 -0.095432808
# C12_1                        0.7704083320  2.544655e-01 -0.295785705
# C12                          0.6877571749  5.486462e-01 -0.248081134
# C14_2                        0.7408389055  4.588685e-01 -0.290253173
# C14_1                        0.8201879961  4.338776e-01 -0.197332377
# C14                          0.7389157704  4.797589e-01 -0.286380747
# C12DC                        0.5817953428  6.855649e-01  0.069423788
# C14_2OH                      0.6903732747  6.310981e-01 -0.122692319
# C14_1OH                      0.7261232077  4.891200e-01 -0.319186669
# C16_2                        0.6441052735  6.409279e-01 -0.045838107
# C16_1                        0.6828226355  6.126920e-01 -0.172321243
# C16                          0.6793946429  5.500985e-01 -0.353798576
# C16_2OH                      0.7018928302  5.608283e-01 -0.319491432
# C16_1OH                      0.7155825723  5.196280e-01 -0.342779681
# C16OH                        0.7514440092  4.204216e-01 -0.247837541
# C18_2                        0.7307245274  5.847884e-01 -0.140355324
# C18_1                        0.6966550519  6.140602e-01 -0.230088595
# C18                          0.6578760230  5.979459e-01 -0.310826375
# C18_1OH                      0.6993482341  5.485659e-01 -0.304820424
# Arginine_Average             0.2033255607  2.213717e-01  0.829041283
# Betaine_Average              0.0920296894  2.358516e-01  0.257754706
# C0_Average                   0.5632399475  2.835352e-01  0.147090028
# Choline_Average              0.2551280535  3.861413e-01  0.470887767
# Citrate_Average              0.6809853346 -3.525476e-01 -0.092181234
# Creatine_Average             0.5179974112 -3.434798e-01 -0.264345371
# Creatinine_Average           0.5018302496 -4.538630e-01 -0.193884547
# Glucose_Average              0.1484779978 -2.176406e-01  0.557668244
# Glutamate_Average            0.4424847020 -1.854787e-01  0.595076543
# Glycine_Average              0.4791644095 -3.363134e-02  0.604802040
# Histidine_Average            0.5034163087 -2.931094e-01  0.482277506
# Isoleucine_Average           0.5141469854 -5.658942e-01  0.242031471
# Lactate_Average              0.6135054180 -4.157303e-01  0.245389195
# Leucine                      0.3122922382 -3.310391e-01  0.535662259
# Methionine                   0.3947595880 -2.262349e-01  0.762575552
# Phenylalanine_Average        0.4516106571 -2.667886e-01  0.717181252
# Proline                      0.4364906077 -7.350098e-01  0.192696283
# Pyruvate_Average             0.2530927234  3.589106e-01  0.575390543
# Succinate_Average            0.2830452435  2.658768e-01  0.515704188
# Threonine_Average            0.4837490663 -1.554742e-01  0.482015937
# Tyrosine_Average             0.3414637183 -2.617159e-01  0.460227740
# Valine_Average               0.5459043242 -6.353091e-01  0.216771076
# X1_Methylhistidine          -0.0804240608  1.808143e-01  0.087524124
# X2_Hydroxybutyrate           0.2721620834 -4.868825e-01 -0.111726058
# X2_Hydroxyisovaleric_Acid    0.1656597728 -7.918181e-02 -0.245765484
# X3_Hydroxybutyrate           0.7313471814  1.944505e-01 -0.246128089
# ADP                         -0.3750485302  4.784334e-01  0.081909699
# AMP                         -0.0909221413  1.731015e-01  0.077279890
# ATP                         -0.2452139380  9.792143e-02 -0.143644262
# Acetamide                   -0.4140893109  5.294524e-01  0.121589507
# Acetate                      0.4231940011  1.212919e-01 -0.057653642
# Acetoacetate                -0.0876412598 -4.073374e-02 -0.222245302
# Acetone                      0.5026999222  2.085743e-01 -0.248204191
# Adenosine                    0.1316730409  3.197266e-01  0.234288506
# Alanine                      0.2074709870 -3.326186e-05  0.332375503
# Aspartate                    0.2101296176 -5.824379e-01  0.011693200
# Creatine_Phosphate           0.0007150375  2.178862e-01  0.216040381
# Cytidine                    -0.0636554674 -3.089267e-01  0.147826481
# Dimethyl_Sulfone             0.2744668778 -4.714872e-01 -0.138804101
# Ethanol                     -0.2045567861 -1.013160e-01 -0.139556450
# Ethanolamine                -0.2043794503  5.389669e-01  0.192163517
# Formate                     -0.0614294023 -1.841225e-01 -0.040709584
# Glycerol                     0.0254736829  5.675480e-02  0.093624476
# Guanidoacetate              -0.1171218232  2.572557e-01  0.251728061
# Guanosine                   -0.0706625348  3.454734e-01  0.308589852
# Hypoxanthine                 0.0236272310  4.640300e-01  0.147809752
# IMP                          0.2714798668 -4.605100e-01  0.002859496
# Inosine                      0.0452340400 -5.596950e-01 -0.147473638
# Isopropanol                 -0.2051232904  2.876024e-01  0.099508126
# Malate                       0.1811247476  1.142811e-01  0.131985627
# Malonate                     0.2116734450 -1.132763e-01 -0.017668753
# Mannose                      0.0829232049 -1.991459e-01  0.501689589
# Methanol                    -0.0324237789 -2.704743e-01  0.020256734
# N_N_Dimethylglycine          0.0266771463  4.050526e-01  0.172218357
# Nicotinurate                 0.3745107550 -1.737080e-01  0.270295113
# O_Acetylcarnitine            0.4381890299  2.234834e-01 -0.051099468
# O_Phosphocholine             0.0997167489 -4.290136e-01  0.308871515
# Propylene_Glycol             0.2252959034 -2.927906e-01 -0.309557275
# Taurine                      0.0251608463 -5.043755e-02  0.076932559
# Uridine                      0.4080026957 -3.298225e-01  0.364696018
# Xanthine                    -0.3014628836 -3.612490e-01 -0.228064247
# sn_Glycero_3_Phosphocholine -0.0978326915  3.902570e-01  0.337639963
# Beta.Alanine                -0.0457184323  1.099851e-01  0.100532207
# Dim.4        Dim.5
# Serine                       0.141707923 -0.072110979
# Putrescine                  -0.198054645 -0.224959148
# Trans_Hydroxyproline         0.089337327  0.135352890
# Asparagine                   0.092267751  0.132101856
# Glutamine                    0.134433428 -0.003097557
# Alpha_Aminoadipic_Acid      -0.086629498  0.042200412
# Methionine_Sulfoxide         0.157217479 -0.057435746
# Acetyl_Ornithine            -0.045211770  0.223969212
# Citrulline                  -0.244060610  0.683875877
# Asymmetric_Dimethylarginine  0.101089166 -0.038802072
# Total_Dimethylarginine       0.234538449 -0.083828897
# Tryptophan                   0.212365810  0.208080249
# Kynurenine                   0.095075464 -0.027800938
# Ornithine                    0.083792306 -0.212176013
# Lysine                       0.066094784 -0.266329897
# Spermidine                  -0.069431599 -0.092154926
# Spermine                    -0.090917196 -0.116913340
# Sarcosine                    0.147903220  0.090962189
# Methylhistidine              0.136112641 -0.163331064
# Beta_Hydroxybutyric_Acid     0.426006495  0.025823913
# Alpha_Ketoglutaric_Acid      0.042457961  0.212717211
# Butyric_acid                -0.055005250 -0.252883062
# Propionic_Acid              -0.433745761 -0.478725219
# Fumaric_Acid                 0.153260394  0.378776506
# Isobutyric_Acid              0.193064148  0.041052492
# Hippuric_Acid               -0.077296822 -0.038149710
# Methylmalonic_Acid          -0.095775615 -0.110142317
# LYSOC14_0                   -0.656895485  0.357997010
# LYSOC16_1                   -0.517767756  0.574548239
# LYSOC16_0                   -0.490964891  0.535219140
# LYSOC17_0                   -0.672787241  0.386110277
# LYSOC18_2                   -0.455712288  0.154886100
# LYSOC18_1                   -0.385638861  0.167200321
# LYSOC18_0                   -0.264601064 -0.015623480
# LYSOC20_4                   -0.205203551  0.152310083
# LYSOC20_3                   -0.519527410 -0.239809893
# LYSOC24_0                   -0.411591171 -0.132004395
# LYSOC26_1                   -0.529776154 -0.217783615
# LYSOC26_0                   -0.551105765 -0.099382254
# LYSOC28_1                   -0.535446572 -0.136003634
# LYSOC28_0                   -0.488220487 -0.262940397
# X14_1SMOH                   -0.432192067 -0.012662074
# X16_1SM                     -0.021477841  0.066672303
# X16_0SM                     -0.093915481  0.065186131
# X16_1SMOH                   -0.227550946  0.244401800
# X18_1SM                     -0.382566131  0.004875196
# PC32_2AA                    -0.412734485 -0.055073043
# X18_0SM                     -0.318298606 -0.049998967
# X20_2SM                     -0.495816293 -0.082241963
# PC36_0AE                    -0.302562069  0.300529817
# PC36_6AA                    -0.692037913  0.114620242
# PC36_0AA                    -0.268587243  0.291080325
# X22_2SMOH                   -0.233421371  0.078239565
# X22_1SMOH                   -0.379222207  0.089800261
# PC38_6AA                    -0.599903495  0.059852547
# PC38_0AA                    -0.405144778  0.150140184
# PC40_6AE                    -0.548066342 -0.020486290
# X24_1SMOH                   -0.530313212  0.090371428
# PC40_6AA                    -0.270115914 -0.195519281
# PC40_2AA                    -0.623425249 -0.025590754
# PC401AA                     -0.502653829 -0.083951452
# C2                           0.094916219 -0.451202412
# C3_1                         0.036666434 -0.080751686
# C3                           0.026322320 -0.004643628
# C4_1                        -0.333267096 -0.131699120
# C4                           0.253635942 -0.015504451
# C3OH                        -0.353174360  0.064103556
# C5_1                         0.028368931 -0.285839130
# C5                           0.105542820 -0.008976080
# C4OH                         0.127071855 -0.061655117
# C6_1                        -0.130399561  0.123338929
# C6                           0.172868645 -0.100151908
# C5OH                        -0.204402477 -0.195453064
# C5_1DC                       0.144071529 -0.023647914
# C5DC                         0.091211042 -0.373736535
# C8                           0.250582232 -0.106670806
# C5MDC                        0.251398454 -0.018307488
# C9                           0.238465756 -0.016962283
# C7DC                        -0.138255005  0.209594290
# C10_2                        0.166567227 -0.163395950
# C10_1                        0.009902967 -0.071983467
# C10                          0.301523852 -0.078595127
# C12_1                        0.235819483 -0.106257174
# C12                          0.295380827 -0.031959770
# C14_2                        0.122038071  0.090632623
# C14_1                        0.216388059 -0.035172202
# C14                          0.287876706 -0.033172396
# C12DC                        0.167013424 -0.012317287
# C14_2OH                      0.249775876 -0.026578439
# C14_1OH                      0.272669204 -0.002878986
# C16_2                        0.292209587 -0.035379015
# C16_1                        0.285898651  0.007754861
# C16                          0.211811847 -0.020240412
# C16_2OH                      0.207864919 -0.002479764
# C16_1OH                      0.157803008  0.012934858
# C16OH                        0.284862386 -0.069996579
# C18_2                        0.256114669 -0.037660276
# C18_1                        0.226909595 -0.015327566
# C18                          0.132818400  0.012129758
# C18_1OH                      0.257730290 -0.019980881
# Arginine_Average             0.007316440 -0.253789375
# Betaine_Average             -0.176840384  0.551875543
# C0_Average                   0.088271351 -0.276892978
# Choline_Average             -0.296141858  0.147653604
# Citrate_Average              0.257479242  0.155117775
# Creatine_Average             0.056278307 -0.085473268
# Creatinine_Average           0.082223412 -0.208769066
# Glucose_Average             -0.440289279 -0.093368012
# Glutamate_Average            0.143291397  0.188026266
# Glycine_Average              0.282033406 -0.026430568
# Histidine_Average           -0.119715964 -0.146141203
# Isoleucine_Average           0.186861934  0.097914581
# Lactate_Average              0.001580538  0.035936584
# Leucine                      0.191611865 -0.029674196
# Methionine                   0.218601497 -0.008793314
# Phenylalanine_Average        0.160255610  0.052469349
# Proline                      0.083141004  0.128199944
# Pyruvate_Average             0.018641003 -0.274371990
# Succinate_Average            0.089967722  0.493544536
# Threonine_Average            0.084886210 -0.345785483
# Tyrosine_Average             0.257237091  0.185851119
# Valine_Average               0.196472378  0.150898183
# X1_Methylhistidine           0.119522714 -0.237365996
# X2_Hydroxybutyrate           0.362198392  0.287684261
# X2_Hydroxyisovaleric_Acid    0.224286271  0.199090903
# X3_Hydroxybutyrate           0.447200341 -0.011056826
# ADP                         -0.062182864  0.315443933
# AMP                          0.112136273  0.038035575
# ATP                          0.126890118  0.210165785
# Acetamide                   -0.034751621  0.179893292
# Acetate                      0.452412191  0.214621328
# Acetoacetate                 0.232685397 -0.331163265
# Acetone                      0.439487219  0.216403778
# Adenosine                    0.312260013  0.399910824
# Alanine                      0.145225243  0.144223716
# Aspartate                    0.298889303  0.097278636
# Creatine_Phosphate           0.008578474  0.121030586
# Cytidine                     0.240005851  0.020171264
# Dimethyl_Sulfone             0.312551082  0.221828841
# Ethanol                     -0.012997664 -0.322268361
# Ethanolamine                 0.104257529  0.488790780
# Formate                     -0.108869700 -0.336004872
# Glycerol                    -0.036939686  0.039220113
# Guanidoacetate               0.152485887  0.264361136
# Guanosine                    0.193204568  0.571802101
# Hypoxanthine                -0.014623642  0.455229317
# IMP                          0.257093633  0.064339972
# Inosine                      0.126135495  0.045802679
# Isopropanol                 -0.135300771 -0.038911645
# Malate                       0.375498159  0.068322560
# Malonate                     0.318454733  0.197388684
# Mannose                      0.263251698 -0.338204743
# Methanol                    -0.152435151 -0.518511513
# N_N_Dimethylglycine          0.191564666 -0.055427840
# Nicotinurate                 0.401602157  0.391571193
# O_Acetylcarnitine            0.113295718 -0.497703641
# O_Phosphocholine             0.200934604  0.268582337
# Propylene_Glycol             0.448414351 -0.053146594
# Taurine                      0.301468719  0.249201943
# Uridine                      0.383150458  0.158159150
# Xanthine                     0.226596361 -0.142819951
# sn_Glycero_3_Phosphocholine -0.103479177  0.240171629
# Beta.Alanine                -0.102312176 -0.089427228



var$contrib
# Dim.1        Dim.2        Dim.3
# Serine                      3.060732e-01 1.856674e-01 4.128270e+00
# Putrescine                  7.495713e-01 1.058285e+00 4.953490e-02
# Trans_Hydroxyproline        2.495178e-01 2.693930e+00 2.523181e-02
# Asparagine                  2.460460e-01 2.457649e+00 3.982996e-03
# Glutamine                   5.359932e-01 1.318357e-01 2.303882e+00
# Alpha_Aminoadipic_Acid      3.428048e-01 4.556626e-01 4.510247e-01
# Methionine_Sulfoxide        1.673005e-01 1.054048e-01 3.223197e+00
# Acetyl_Ornithine            5.287436e-01 3.174643e-01 1.318784e-01
# Citrulline                  9.259353e-02 1.992428e-01 6.824121e-01
# Asymmetric_Dimethylarginine 7.881652e-01 1.386232e+00 5.911376e-01
# Total_Dimethylarginine      1.102487e+00 8.617422e-01 9.289493e-01
# Tryptophan                  1.060925e+00 1.276204e+00 2.627094e-01
# Kynurenine                  3.864591e-01 1.595873e-01 1.189357e+00
# Ornithine                   6.877298e-02 2.880034e-01 3.804083e+00
# Lysine                      2.586399e-03 3.682732e-01 3.396152e+00
# Spermidine                  7.995835e-01 1.247099e+00 4.457295e-02
# Spermine                    8.192701e-01 8.425361e-01 4.622667e-02
# Sarcosine                   2.515123e-01 1.265678e-01 2.636938e+00
# Methylhistidine             2.010959e-02 3.401699e-01 5.740297e-01
# Beta_Hydroxybutyric_Acid    1.487559e+00 8.841111e-02 3.435478e-01
# Alpha_Ketoglutaric_Acid     8.259382e-01 8.405318e-02 1.686745e+00
# Butyric_acid                6.552485e-01 5.142079e-01 5.461658e-01
# Propionic_Acid              5.272589e-03 2.671238e-02 5.671000e-04
# Fumaric_Acid                6.323957e-01 4.167131e-04 1.480950e+00
# Isobutyric_Acid             4.692588e-01 2.052787e+00 6.235841e-01
# Hippuric_Acid               2.138564e-02 9.662212e-01 1.795262e-01
# Methylmalonic_Acid          9.760391e-02 1.945781e-01 1.494889e-02
# LYSOC14_0                   3.879425e-01 1.127903e-01 2.750236e-03
# LYSOC16_1                   2.707263e-01 1.448092e-01 5.256792e-02
# LYSOC16_0                   5.115798e-02 4.409615e-01 9.529702e-02
# LYSOC17_0                   2.361783e-01 5.159197e-01 1.289154e-01
# LYSOC18_2                   9.423643e-01 5.684782e-02 2.677517e-01
# LYSOC18_1                   8.744106e-01 2.909702e-02 6.025427e-01
# LYSOC18_0                   7.520087e-01 1.604716e-02 1.665605e+00
# LYSOC20_4                   1.118063e+00 9.039302e-01 1.502017e-02
# LYSOC20_3                   2.233549e-01 8.222841e-03 2.319748e-01
# LYSOC24_0                   1.595424e-01 8.229583e-01 1.652807e-02
# LYSOC26_1                   5.695730e-01 1.502942e-01 6.368490e-04
# LYSOC26_0                   1.077372e-01 6.993057e-01 5.851033e-03
# LYSOC28_1                   1.006876e+00 9.224752e-03 1.408055e-01
# LYSOC28_0                   2.948964e-01 2.341834e-01 2.614955e-02
# X14_1SMOH                   1.068615e+00 4.203159e-02 1.039445e+00
# X16_1SM                     1.370298e+00 8.964886e-01 8.195068e-01
# X16_0SM                     1.290749e+00 1.034960e+00 7.921092e-01
# X16_1SMOH                   4.504940e-01 1.247247e+00 1.199929e+00
# X18_1SM                     1.257235e+00 3.524271e-02 7.101897e-01
# PC32_2AA                    4.205606e-01 7.747717e-01 1.523168e+00
# X18_0SM                     5.965704e-01 6.867289e-01 1.678330e+00
# X20_2SM                     9.346474e-01 1.607908e-01 1.294913e+00
# PC36_0AE                    9.890649e-01 1.109378e+00 4.166179e-01
# PC36_6AA                    9.800002e-01 4.293315e-04 5.125631e-02
# PC36_0AA                    4.513947e-01 1.967485e+00 4.882651e-01
# X22_2SMOH                   1.416001e+00 3.872629e-01 1.228273e+00
# X22_1SMOH                   1.751618e+00 1.571648e-01 2.915640e-01
# PC38_6AA                    1.174213e+00 1.907805e-02 4.575609e-03
# PC38_0AA                    1.437193e+00 4.482051e-01 4.843704e-01
# PC40_6AE                    1.343086e+00 2.198098e-02 9.368922e-03
# X24_1SMOH                   8.784472e-01 3.190026e-01 7.434142e-04
# PC40_6AA                    9.229274e-01 4.600200e-02 1.400654e+00
# PC40_2AA                    1.124907e+00 1.544267e-02 1.127138e-01
# PC401AA                     1.180108e+00 2.521798e-01 7.815622e-02
# C2                          6.641628e-01 3.525778e-01 2.822431e-03
# C3_1                        6.518235e-01 3.442601e-01 2.715146e-01
# C3                          9.420827e-01 1.494174e+00 9.952631e-05
# C4_1                        7.696965e-01 6.960959e-01 2.643348e-01
# C4                          1.670450e+00 3.805497e-01 3.635930e-01
# C3OH                        4.424947e-01 6.338810e-02 2.803433e-03
# C5_1                        5.752794e-01 1.153824e+00 5.220513e-02
# C5                          6.394193e-01 1.721202e+00 2.352257e-02
# C4OH                        8.540397e-01 4.045285e-03 4.385411e-01
# C6_1                        3.767558e-01 1.509252e-02 2.281420e+00
# C6                          1.167710e+00 5.165046e-01 7.030004e-01
# C5OH                        6.099273e-01 1.547217e-02 4.865572e-01
# C5_1DC                      1.271721e+00 1.537060e-01 3.994608e-01
# C5DC                        8.932011e-01 2.991125e-02 1.402365e-01
# C8                          1.432141e+00 8.181319e-01 4.481788e-02
# C5MDC                       5.390507e-01 6.562076e-01 1.465818e-02
# C9                          1.467526e+00 1.021623e+00 4.149682e-02
# C7DC                        2.157274e-02 9.676077e-02 7.146699e-03
# C10_2                       4.050171e-01 1.687544e+00 9.769907e-04
# C10_1                       4.484194e-01 2.766457e-02 2.412650e-01
# C10                         1.162603e+00 1.277625e+00 5.528733e-02
# C12_1                       1.582200e+00 2.892006e-01 5.311101e-01
# C12                         1.260926e+00 1.344393e+00 3.736093e-01
# C14_2                       1.463076e+00 9.404115e-01 5.114276e-01
# C14_1                       1.793272e+00 8.407676e-01 2.363888e-01
# C14                         1.455490e+00 1.027987e+00 4.978722e-01
# C12DC                       9.023187e-01 2.099127e+00 2.925815e-02
# C14_2OH                     1.270537e+00 1.778833e+00 9.138291e-02
# C14_1OH                     1.405530e+00 1.068494e+00 6.184715e-01
# C16_2                       1.105944e+00 1.834678e+00 1.275509e-02
# C16_1                       1.242897e+00 1.676586e+00 1.802635e-01
# C16                         1.230449e+00 1.351519e+00 7.598754e-01
# C16_2OH                     1.313291e+00 1.404757e+00 6.196532e-01
# C16_1OH                     1.365020e+00 1.205942e+00 7.132805e-01
# C16OH                       1.505264e+00 7.894259e-01 3.728760e-01
# C18_2                       1.423400e+00 1.527351e+00 1.195881e-01
# C18_1                       1.293764e+00 1.684082e+00 3.213812e-01
# C18                         1.153739e+00 1.596854e+00 5.864972e-01
# C18_1OH                     1.303786e+00 1.343999e+00 5.640510e-01
# Arginine_Average            1.102055e-01 2.188697e-01 4.172367e+00
# Betaine_Average             2.257748e-02 2.484385e-01 4.033141e-01
# C0_Average                  8.456806e-01 3.590503e-01 1.313399e-01
# Choline_Average             1.735146e-01 6.659381e-01 1.346062e+00
# Citrate_Average             1.236218e+00 5.551073e-01 5.158403e-02
# Creatine_Average            7.152776e-01 5.269190e-01 4.242029e-01
# Creatinine_Average          6.713255e-01 9.200068e-01 2.282005e-01
# Glucose_Average             5.876837e-02 2.115540e-01 1.887913e+00
# Glutamate_Average           5.219345e-01 1.536489e-01 2.149690e+00
# Glycine_Average             6.120524e-01 5.051608e-03 2.220530e+00
# Histidine_Average           6.755757e-01 3.837082e-01 1.411966e+00
# Isoleucine_Average          7.046834e-01 1.430250e+00 3.556100e-01
# Lactate_Average             1.003359e+00 7.719068e-01 3.655452e-01
# Leucine                     2.599811e-01 4.894407e-01 1.741856e+00
# Methionine                  4.154176e-01 2.285918e-01 3.530173e+00
# Phenylalanine_Average       5.436856e-01 3.178893e-01 3.122397e+00
# Proline                     5.078895e-01 2.412835e+00 2.254119e-01
# Pyruvate_Average            1.707571e-01 5.753260e-01 2.009813e+00
# Succinate_Average           2.135655e-01 3.157203e-01 1.614476e+00
# Threonine_Average           6.238207e-01 1.079588e-01 1.410435e+00
# Tyrosine_Average            3.108197e-01 3.059156e-01 1.285807e+00
# Valine_Average              7.944243e-01 1.802651e+00 2.852548e-01
# X1_Methylhistidine          1.724215e-02 1.460182e-01 4.650351e-02
# X2_Hydroxybutyrate          1.974580e-01 1.058741e+00 7.577731e-02
# X2_Hydroxyisovaleric_Acid   7.315661e-02 2.800221e-02 3.666672e-01
# X3_Hydroxybutyrate          1.425826e+00 1.688728e-01 3.677499e-01
# ADP                         3.749682e-01 1.022314e+00 4.072873e-02
# AMP                         2.203732e-02 1.338268e-01 3.625461e-02
# ATP                         1.602912e-01 4.282496e-02 1.252584e-01
# Acetamide                   4.570961e-01 1.251974e+00 8.974751e-02
# Acetate                     4.774177e-01 6.570604e-02 2.017826e-02
# Acetoacetate                2.047560e-02 7.410542e-03 2.998440e-01
# Acetone                     6.736543e-01 1.942956e-01 3.739801e-01
# Adenosine                   4.621822e-02 4.565610e-01 3.332209e-01
# Alanine                     1.147451e-01 4.941222e-09 6.706381e-01
# Aspartate                   1.177048e-01 1.515098e+00 8.300360e-04
# Creatine_Phosphate          1.362941e-06 2.120317e-01 2.833349e-01
# Cytidine                    1.080167e-02 4.262380e-01 1.326584e-01
# Dimethyl_Sulfone            2.008165e-01 9.928445e-01 1.169593e-01
# Ethanol                     1.115443e-01 4.584557e-02 1.182306e-01
# Ethanolamine                1.113509e-01 1.297375e+00 2.241672e-01
# Formate                     1.005940e-02 1.514102e-01 1.006059e-02
# Glycerol                    1.729828e-03 1.438622e-02 5.321193e-02
# Guanidoacetate              3.656749e-02 2.955775e-01 3.846746e-01
# Guanosine                   1.331061e-02 5.330535e-01 5.780874e-01
# Hypoxanthine                1.488145e-03 9.616868e-01 1.326284e-01
# IMP                         1.964693e-01 9.471519e-01 4.963744e-05
# Inosine                     5.454441e-03 1.399085e+00 1.320259e-01
# Isopropanol                 1.121629e-01 3.694252e-01 6.011008e-02
# Malate                      8.745307e-02 5.832975e-02 1.057508e-01
# Malonate                    1.194407e-01 5.730858e-02 1.895143e-03
# Mannose                     1.833039e-02 1.771267e-01 1.527919e+00
# Methanol                    2.802507e-03 3.267333e-01 2.490973e-03
# N_N_Dimethylglycine         1.897135e-03 7.327643e-01 1.800483e-01
# Nicotinurate                3.738936e-01 1.347662e-01 4.435133e-01
# O_Acetylcarnitine           5.118497e-01 2.230652e-01 1.585123e-02
# O_Phosphocholine            2.650670e-02 8.220221e-01 5.791432e-01
# Propylene_Glycol            1.353088e-01 3.828739e-01 5.817177e-01
# Taurine                     1.687602e-03 1.136186e-02 3.592945e-02
# Uridine                     4.437573e-01 4.858496e-01 8.074066e-01
# Xanthine                    2.422630e-01 5.828473e-01 3.157509e-01
# sn_Glycero_3_Phosphocholine 2.551452e-02 6.802097e-01 6.920507e-01
# Beta.Alanine                5.571885e-03 5.402675e-02 6.135369e-02
# Dim.4        Dim.5
# Serine                      1.555603e-01 6.810204e-02
# Putrescine                  3.038650e-01 6.627728e-01
# Trans_Hydroxyproline        6.182676e-02 2.399343e-01
# Asparagine                  6.594934e-02 2.285468e-01
# Glutamine                   1.399991e-01 1.256597e-04
# Alpha_Aminoadipic_Acid      5.813560e-02 2.332334e-02
# Methionine_Sulfoxide        1.914751e-01 4.320379e-02
# Acetyl_Ornithine            1.583484e-02 6.569525e-01
# Citrulline                  4.614301e-01 6.125082e+00
# Asymmetric_Dimethylarginine 7.916255e-02 1.971820e-02
# Total_Dimethylarginine      4.261266e-01 9.203327e-02
# Tryptophan                  3.493653e-01 5.670470e-01
# Kynurenine                  7.002409e-02 1.012223e-02
# Ornithine                   5.438997e-02 5.895897e-01
# Lysine                      3.384114e-02 9.289604e-01
# Spermidine                  3.734435e-02 1.112230e-01
# Spermine                    6.403282e-02 1.790134e-01
# Sarcosine                   1.694595e-01 1.083625e-01
# Methylhistidine             1.435184e-01 3.493775e-01
# Beta_Hydroxybutyric_Acid    1.405863e+00 8.733764e-03
# Alpha_Ketoglutaric_Acid     1.396461e-02 5.926013e-01
# Butyric_acid                2.343791e-02 8.375231e-01
# Propionic_Acid              1.457408e+00 3.001442e+00
# Fumaric_Acid                1.819577e-01 1.878985e+00
# Isobutyric_Acid             2.887446e-01 2.207173e-02
# Hippuric_Acid               4.628431e-02 1.906074e-02
# Methylmalonic_Acid          7.105922e-02 1.588787e-01
# LYSOC14_0                   3.342744e+00 1.678480e+00
# LYSOC16_1                   2.076733e+00 4.323250e+00
# LYSOC16_0                   1.867289e+00 3.751635e+00
# LYSOC17_0                   3.506437e+00 1.952451e+00
# LYSOC18_2                   1.608763e+00 3.141827e-01
# LYSOC18_1                   1.152052e+00 3.661268e-01
# LYSOC18_0                   5.423676e-01 3.196781e-03
# LYSOC20_4                   3.261973e-01 3.038188e-01
# LYSOC20_3                   2.090872e+00 7.531674e-01
# LYSOC24_0                   1.312329e+00 2.282097e-01
# LYSOC26_1                   2.174179e+00 6.211661e-01
# LYSOC26_0                   2.352775e+00 1.293526e-01
# LYSOC28_1                   2.220971e+00 2.422470e-01
# LYSOC28_0                   1.846471e+00 9.054656e-01
# X14_1SMOH                   1.446985e+00 2.099747e-03
# X16_1SM                     3.573484e-03 5.821679e-02
# X16_0SM                     6.832583e-02 5.565033e-02
# X16_1SMOH                   4.011140e-01 7.822870e-01
# X18_1SM                     1.133766e+00 3.112731e-04
# PC32_2AA                    1.319630e+00 3.972240e-02
# X18_0SM                     7.848380e-01 3.274005e-02
# X20_2SM                     1.904374e+00 8.858175e-02
# PC36_0AE                    7.091523e-01 1.182858e+00
# PC36_6AA                    3.709969e+00 1.720600e-01
# PC36_0AA                    5.588321e-01 1.109642e+00
# X22_2SMOH                   4.220771e-01 8.016969e-02
# X22_1SMOH                   1.114033e+00 1.056118e-01
# PC38_6AA                    2.787875e+00 4.691618e-02
# PC38_0AA                    1.271543e+00 2.952237e-01
# PC40_6AE                    2.326895e+00 5.496471e-03
# X24_1SMOH                   2.178590e+00 1.069596e-01
# PC40_6AA                    5.652114e-01 5.006526e-01
# PC40_2AA                    3.010782e+00 8.576765e-03
# PC401AA                     1.957260e+00 9.230256e-02
# C2                          6.978971e-02 2.666246e+00
# C3_1                        1.041473e-02 8.540052e-02
# C3                          5.367337e-03 2.824049e-04
# C4_1                        8.603902e-01 2.271554e-01
# C4                          4.983474e-01 3.148256e-03
# C3OH                        9.662486e-01 5.381726e-02
# C5_1                        6.234426e-03 1.070042e+00
# C5                          8.629149e-02 1.055190e-03
# C4OH                        1.250862e-01 4.978467e-02
# C6_1                        1.317234e-01 1.992314e-01
# C6                          2.314957e-01 1.313638e-01
# C5OH                        3.236554e-01 5.003136e-01
# C5_1DC                      1.607929e-01 7.323910e-03
# C5DC                        6.444740e-02 1.829315e+00
# C8                          4.864196e-01 1.490213e-01
# C5MDC                       4.895936e-01 4.389498e-03
# C9                          4.405169e-01 3.768131e-03
# C7DC                        1.480718e-01 5.753290e-01
# C10_2                       2.149264e-01 3.496551e-01
# C10_1                       7.596983e-04 6.786141e-02
# C10                         7.042939e-01 8.090001e-02
# C12_1                       4.307943e-01 1.478679e-01
# C12                         6.758887e-01 1.337718e-02
# C14_2                       1.153722e-01 1.075787e-01
# C14_1                       3.627247e-01 1.620155e-02
# C14                         6.419831e-01 1.441156e-02
# C12DC                       2.160794e-01 1.986952e-03
# C14_2OH                     4.832942e-01 9.251588e-03
# C14_1OH                     5.759473e-01 1.085517e-04
# C16_2                       6.614537e-01 1.639264e-02
# C16_1                       6.331911e-01 7.875993e-04
# C16                         3.475450e-01 5.365325e-03
# C16_2OH                     3.347133e-01 8.053378e-05
# C16_1OH                     1.929039e-01 2.191193e-03
# C16OH                       6.286093e-01 6.416689e-02
# C18_2                       5.081354e-01 1.857481e-02
# C18_1                       3.988561e-01 3.076832e-03
# C18                         1.366555e-01 1.926911e-03
# C18_1OH                     5.145665e-01 5.228614e-03
# Arginine_Average            4.146773e-04 8.435371e-01
# Betaine_Average             2.422554e-01 3.988776e+00
# C0_Average                  6.036012e-02 1.004110e+00
# Choline_Average             6.793759e-01 2.855259e-01
# Citrate_Average             5.135645e-01 3.151233e-01
# Creatine_Average            2.453537e-02 9.567929e-02
# Creatinine_Average          5.237229e-02 5.708075e-01
# Glucose_Average             1.501712e+00 1.141704e-01
# Glutamate_Average           1.590563e-01 4.630143e-01
# Glycine_Average             6.161858e-01 9.148931e-03
# Histidine_Average           1.110234e-01 2.797066e-01
# Isoleucine_Average          2.704906e-01 1.255602e-01
# Lactate_Average             1.935175e-05 1.691340e-02
# Leucine                     2.844169e-01 1.153228e-02
# Methionine                  3.701833e-01 1.012657e-03
# Phenylalanine_Average       1.989468e-01 3.605527e-02
# Proline                     5.354773e-02 2.152450e-01
# Pyruvate_Average            2.691839e-03 9.859091e-01
# Succinate_Average           6.270238e-02 3.190142e+00
# Threonine_Average           5.581936e-02 1.565925e+00
# Tyrosine_Average            5.125990e-01 4.523636e-01
# Valine_Average              2.990292e-01 2.982122e-01
# X1_Methylhistidine          1.106653e-01 7.378946e-01
# X2_Hydroxybutyrate          1.016257e+00 1.083901e+00
# X2_Hydroxyisovaleric_Acid   3.896870e-01 5.191109e-01
# X3_Hydroxybutyrate          1.549226e+00 1.601099e-03
# ADP                         2.995380e-02 1.303172e+00
# AMP                         9.740983e-02 1.894686e-02
# ATP                         1.247286e-01 5.784707e-01
# Acetamide                   9.355364e-03 4.238257e-01
# Acetate                     1.585547e+00 6.032580e-01
# Acetoacetate                4.194197e-01 1.436288e+00
# Acetone                     1.496246e+00 6.133198e-01
# Adenosine                   7.553414e-01 2.094516e+00
# Alanine                     1.633785e-01 2.724148e-01
# Aspartate                   6.920402e-01 1.239345e-01
# Creatine_Phosphate          5.700732e-04 1.918437e-01
# Cytidine                    4.462253e-01 5.328728e-03
# Dimethyl_Sulfone            7.567503e-01 6.444561e-01
# Ethanol                     1.308703e-03 1.360168e+00
# Ethanolamine                8.420258e-02 3.128984e+00
# Formate                     9.181732e-02 1.478592e+00
# Glycerol                    1.057053e-02 2.014536e-02
# Guanidoacetate              1.801233e-01 9.152770e-01
# Guanosine                   2.891647e-01 4.282021e+00
# Hypoxanthine                1.656615e-03 2.714050e+00
# IMP                         5.120274e-01 5.421495e-02
# Inosine                     1.232495e-01 2.747509e-02
# Isopropanol                 1.418114e-01 1.982972e-02
# Malate                      1.092260e+00 6.113439e-02
# Malonate                    7.856081e-01 5.102721e-01
# Mannose                     5.368500e-01 1.498017e+00
# Methanol                    1.800034e-01 3.521066e+00
# N_N_Dimethylglycine         2.842768e-01 4.023586e-02
# Nicotinurate                1.249403e+00 2.008070e+00
# O_Acetylcarnitine           9.943460e-02 3.244135e+00
# O_Phosphocholine            3.127664e-01 9.447399e-01
# Propylene_Glycol            1.557649e+00 3.699203e-02
# Taurine                     7.040363e-01 8.133176e-01
# Uridine                     1.137232e+00 3.276016e-01
# Xanthine                    3.977557e-01 2.671376e-01
# sn_Glycero_3_Phosphocholine 8.295002e-02 7.554413e-01
# Beta.Alanine                8.108961e-02 1.047362e-01

var$cos2
# Dim.1        Dim.2        Dim.3
# Serine                      1.148169e-01 4.157136e-02 6.800453e-01
# Putrescine                  2.811859e-01 2.369525e-01 8.159828e-03
# Trans_Hydroxyproline        9.360135e-02 6.031771e-01 4.156408e-03
# Asparagine                  9.229899e-02 5.502732e-01 6.561146e-04
# Glutamine                   2.010665e-01 2.951831e-02 3.795160e-01
# Alpha_Aminoadipic_Acid      1.285960e-01 1.020239e-01 7.429680e-02
# Methionine_Sulfoxide        6.275925e-02 2.360039e-02 5.309537e-01
# Acetyl_Ornithine            1.983470e-01 7.108099e-02 2.172418e-02
# Citrulline                  3.473451e-02 4.461093e-02 1.124130e-01
# Asymmetric_Dimethylarginine 2.956635e-01 3.103806e-01 9.737744e-02
# Total_Dimethylarginine      4.135748e-01 1.929461e-01 1.530248e-01
# Tryptophan                  3.979836e-01 2.857451e-01 4.327582e-02
# Kynurenine                  1.449720e-01 3.573197e-02 1.959215e-01
# Ornithine                   2.579873e-02 6.448463e-02 6.266423e-01
# Lysine                      9.702332e-04 8.245721e-02 5.594443e-01
# Spermidine                  2.999469e-01 2.792284e-01 7.342452e-03
# Spermine                    3.073319e-01 1.886458e-01 7.614869e-03
# Sarcosine                   9.434953e-02 2.833882e-02 4.343798e-01
# Methylhistidine             7.543687e-03 7.616482e-02 9.455928e-02
# Beta_Hydroxybutyric_Acid    5.580264e-01 1.979545e-02 5.659225e-02
# Alpha_Ketoglutaric_Acid     3.098333e-01 1.881970e-02 2.778556e-01
# Butyric_acid                2.458026e-01 1.151323e-01 8.996929e-02
# Propionic_Acid              1.977901e-03 5.980963e-03 9.341775e-05
# Fumaric_Acid                2.372299e-01 9.330303e-05 2.439552e-01
# Isobutyric_Acid             1.760325e-01 4.596237e-01 1.027223e-01
# Hippuric_Acid               8.022371e-03 2.163392e-01 2.957315e-02
# Methylmalonic_Acid          3.661405e-02 4.356649e-02 2.462514e-03
# LYSOC14_0                   1.455284e-01 2.525400e-02 4.530433e-04
# LYSOC16_1                   1.015573e-01 3.242311e-02 8.659455e-03
# LYSOC16_0                   1.919084e-02 9.873231e-02 1.569817e-02
# LYSOC17_0                   8.859731e-02 1.155156e-01 2.123609e-02
# LYSOC18_2                   3.535081e-01 1.272836e-02 4.410643e-02
# LYSOC18_1                   3.280167e-01 6.514890e-03 9.925619e-02
# LYSOC18_0                   2.821002e-01 3.592996e-03 2.743732e-01
# LYSOC20_4                   4.194176e-01 2.023921e-01 2.474256e-03
# LYSOC20_3                   8.378689e-02 1.841113e-03 3.821296e-02
# LYSOC24_0                   5.984896e-02 1.842623e-01 2.722650e-03
# LYSOC26_1                   2.136633e-01 3.365121e-02 1.049074e-04
# LYSOC26_0                   4.041534e-02 1.565761e-01 9.638342e-04
# LYSOC28_1                   3.777084e-01 2.065443e-03 2.319473e-02
# LYSOC28_0                   1.106241e-01 5.243420e-02 4.307585e-03
# X14_1SMOH                   4.008685e-01 9.410970e-03 1.712266e-01
# X16_1SM                     5.140384e-01 2.007259e-01 1.349964e-01
# X16_0SM                     4.841971e-01 2.317300e-01 1.304833e-01
# X16_1SMOH                   1.689933e-01 2.792616e-01 1.976630e-01
# X18_1SM                     4.716253e-01 7.890924e-03 1.169888e-01
# PC32_2AA                    1.577644e-01 1.734732e-01 2.509098e-01
# X18_0SM                     2.237908e-01 1.537602e-01 2.764695e-01
# X20_2SM                     3.506133e-01 3.600143e-02 2.133096e-01
# PC36_0AE                    3.710268e-01 2.483923e-01 6.862901e-02
# PC36_6AA                    3.676264e-01 9.612830e-05 8.443395e-03
# PC36_0AA                    1.693312e-01 4.405244e-01 8.043136e-02
# X22_2SMOH                   5.311830e-01 8.670906e-02 2.023321e-01
# X22_1SMOH                   6.570827e-01 3.518957e-02 4.802902e-02
# PC38_6AA                    4.404813e-01 4.271620e-03 7.537350e-04
# PC38_0AA                    5.391325e-01 1.003542e-01 7.978979e-02
# PC40_6AE                    5.038303e-01 4.921591e-03 1.543332e-03
# X24_1SMOH                   3.295309e-01 7.142542e-02 1.224618e-04
# PC40_6AA                    3.462167e-01 1.029995e-02 2.307282e-01
# PC40_2AA                    4.219853e-01 3.457650e-03 1.856721e-02
# PC401AA                     4.426925e-01 5.646364e-02 1.287459e-02
# C2                          2.491467e-01 7.894299e-02 4.649360e-04
# C3_1                        2.445178e-01 7.708064e-02 4.472629e-02
# C3                          3.534024e-01 3.345490e-01 1.639486e-05
# C4_1                        2.887354e-01 1.558575e-01 4.354358e-02
# C4                          6.266341e-01 8.520595e-02 5.989428e-02
# C3OH                        1.659926e-01 1.419274e-02 4.618065e-04
# C5_1                        2.158039e-01 2.583438e-01 8.599693e-03
# C5                          2.398647e-01 3.853811e-01 3.874846e-03
# C4OH                        3.203750e-01 9.057487e-04 7.224039e-02
# C6_1                        1.413320e-01 3.379249e-03 3.758157e-01
# C6                          4.380418e-01 1.156466e-01 1.158045e-01
# C5OH                        2.288014e-01 3.464254e-03 8.015002e-02
# C5_1DC                      4.770594e-01 3.441512e-02 6.580274e-02
# C5DC                        3.350655e-01 6.697198e-03 2.310100e-02
# C8                          5.372374e-01 1.831816e-01 7.382799e-03
# C5MDC                       2.022135e-01 1.469264e-01 2.414625e-03
# C9                          5.505114e-01 2.287438e-01 6.835725e-03
# C7DC                        8.092560e-03 2.166496e-02 1.177268e-03
# C10_2                       1.519336e-01 3.778451e-01 1.609386e-04
# C10_1                       1.682151e-01 6.194162e-03 3.974332e-02
# C10                         4.361261e-01 2.860631e-01 9.107421e-03
# C12_1                       5.935290e-01 6.475269e-02 8.748918e-02
# C12                         4.730099e-01 3.010127e-01 6.154425e-02
# C14_2                       5.488423e-01 2.105603e-01 8.424690e-02
# C14_1                       6.727083e-01 1.882498e-01 3.894007e-02
# C14                         5.459965e-01 2.301686e-01 8.201393e-02
# C12DC                       3.384858e-01 4.699993e-01 4.819662e-03
# C14_2OH                     4.766153e-01 3.982848e-01 1.505341e-02
# C14_1OH                     5.272549e-01 2.392383e-01 1.018801e-01
# C16_2                       4.148716e-01 4.107886e-01 2.101132e-03
# C16_1                       4.662468e-01 3.753915e-01 2.969461e-02
# C16                         4.615771e-01 3.026083e-01 1.251734e-01
# C16_2OH                     4.926535e-01 3.145283e-01 1.020748e-01
# C16_1OH                     5.120584e-01 2.700132e-01 1.174979e-01
# C16OH                       5.646681e-01 1.767543e-01 6.142345e-02
# C18_2                       5.339583e-01 3.419775e-01 1.969962e-02
# C18_1                       4.853283e-01 3.770699e-01 5.294076e-02
# C18                         4.328009e-01 3.575393e-01 9.661304e-02
# C18_1OH                     4.890880e-01 3.009246e-01 9.291549e-02
# Arginine_Average            4.134128e-02 4.900543e-02 6.873094e-01
# Betaine_Average             8.469464e-03 5.562596e-02 6.643749e-02
# C0_Average                  3.172392e-01 8.039220e-02 2.163548e-02
# Choline_Average             6.509032e-02 1.491051e-01 2.217353e-01
# Citrate_Average             4.637410e-01 1.242898e-01 8.497380e-03
# Creatine_Average            2.683213e-01 1.179784e-01 6.987848e-02
# Creatinine_Average          2.518336e-01 2.059917e-01 3.759122e-02
# Glucose_Average             2.204572e-02 4.736742e-02 3.109939e-01
# Glutamate_Average           1.957927e-01 3.440235e-02 3.541161e-01
# Glycine_Average             2.295985e-01 1.131067e-03 3.657855e-01
# Histidine_Average           2.534280e-01 8.591314e-02 2.325916e-01
# Isoleucine_Average          2.643471e-01 3.202363e-01 5.857923e-02
# Lactate_Average             3.763889e-01 1.728317e-01 6.021586e-02
# Leucine                     9.752644e-02 1.095869e-01 2.869341e-01
# Methionine                  1.558351e-01 5.118223e-02 5.815215e-01
# Phenylalanine_Average       2.039522e-01 7.117615e-02 5.143489e-01
# Proline                     1.905241e-01 5.402394e-01 3.713186e-02
# Pyruvate_Average            6.405593e-02 1.288168e-01 3.310743e-01
# Succinate_Average           8.011461e-02 7.069049e-02 2.659508e-01
# Threonine_Average           2.340132e-01 2.417224e-02 2.323394e-01
# Tyrosine_Average            1.165975e-01 6.849520e-02 2.118096e-01
# Valine_Average              2.980115e-01 4.036176e-01 4.698970e-02
# X1_Methylhistidine          6.468030e-03 3.269380e-02 7.660472e-03
# X2_Hydroxybutyrate          7.407220e-02 2.370545e-01 1.248271e-02
# X2_Hydroxyisovaleric_Acid   2.744316e-02 6.269759e-03 6.040067e-02
# X3_Hydroxybutyrate          5.348687e-01 3.781100e-02 6.057904e-02
# ADP                         1.406614e-01 2.288985e-01 6.709199e-03
# AMP                         8.266836e-03 2.996414e-02 5.972181e-03
# ATP                         6.012988e-02 9.588607e-03 2.063367e-02
# Acetamide                   1.714700e-01 2.803199e-01 1.478401e-02
# Acetate                     1.790932e-01 1.471173e-02 3.323942e-03
# Acetoacetate                7.680990e-03 1.659237e-03 4.939297e-02
# Acetone                     2.527072e-01 4.350324e-02 6.160532e-02
# Adenosine                   1.733779e-02 1.022251e-01 5.489110e-02
# Alanine                     4.304421e-02 1.106351e-09 1.104735e-01
# Aspartate                   4.415446e-02 3.392339e-01 1.367309e-04
# Creatine_Phosphate          5.112786e-07 4.747439e-02 4.667345e-02
# Cytidine                    4.052019e-03 9.543568e-02 2.185267e-02
# Dimethyl_Sulfone            7.533207e-02 2.223002e-01 1.926658e-02
# Ethanol                     4.184348e-02 1.026493e-02 1.947600e-02
# Ethanolamine                4.177096e-02 2.904854e-01 3.692682e-02
# Formate                     3.773571e-03 3.390110e-02 1.657270e-03
# Glycerol                    6.489085e-04 3.221107e-03 8.765543e-03
# Guanidoacetate              1.371752e-02 6.618048e-02 6.336702e-02
# Guanosine                   4.993194e-03 1.193519e-01 9.522770e-02
# Hypoxanthine                5.582460e-04 2.153239e-01 2.184772e-02
# IMP                         7.370132e-02 2.120695e-01 8.176720e-06
# Inosine                     2.046118e-03 3.132584e-01 2.174847e-02
# Isopropanol                 4.207556e-02 8.271516e-02 9.901867e-03
# Malate                      3.280617e-02 1.306017e-02 1.742021e-02
# Malonate                    4.480565e-02 1.283152e-02 3.121848e-04
# Mannose                     6.876258e-03 3.965909e-02 2.516924e-01
# Methanol                    1.051301e-03 7.315634e-02 4.103353e-04
# N_N_Dimethylglycine         7.116701e-04 1.640676e-01 2.965916e-02
# Nicotinurate                1.402583e-01 3.017445e-02 7.305945e-02
# O_Acetylcarnitine           1.920096e-01 4.994482e-02 2.611156e-03
# O_Phosphocholine            9.943430e-03 1.840526e-01 9.540161e-02
# Propylene_Glycol            5.075824e-02 8.572636e-02 9.582571e-02
# Taurine                     6.330682e-04 2.543946e-03 5.918619e-03
# Uridine                     1.664662e-01 1.087828e-01 1.330032e-01
# Xanthine                    9.087987e-02 1.305008e-01 5.201330e-02
# sn_Glycero_3_Phosphocholine 9.571236e-03 1.523005e-01 1.140007e-01
# Beta.Alanine                2.090175e-03 1.209671e-02 1.010672e-02
# Dim.4        Dim.5
# Serine                      2.008114e-02 5.199993e-03
# Putrescine                  3.922564e-02 5.060662e-02
# Trans_Hydroxyproline        7.981158e-03 1.832040e-02
# Asparagine                  8.513338e-03 1.745090e-02
# Glutamine                   1.807235e-02 9.594858e-06
# Alpha_Aminoadipic_Acid      7.504670e-03 1.780875e-03
# Methionine_Sulfoxide        2.471734e-02 3.298865e-03
# Acetyl_Ornithine            2.044104e-03 5.016221e-02
# Citrulline                  5.956558e-02 4.676862e-01
# Asymmetric_Dimethylarginine 1.021902e-02 1.505601e-03
# Total_Dimethylarginine      5.500828e-02 7.027284e-03
# Tryptophan                  4.509924e-02 4.329739e-02
# Kynurenine                  9.039344e-03 7.728922e-04
# Ornithine                   7.021151e-03 4.501866e-02
# Lysine                      4.368520e-03 7.093161e-02
# Spermidine                  4.820747e-03 8.492530e-03
# Spermine                    8.265937e-03 1.366873e-02
# Sarcosine                   2.187536e-02 8.274120e-03
# Methylhistidine             1.852665e-02 2.667704e-02
# Beta_Hydroxybutyric_Acid    1.814815e-01 6.668745e-04
# Alpha_Ketoglutaric_Acid     1.802678e-03 4.524861e-02
# Butyric_acid                3.025577e-03 6.394984e-02
# Propionic_Acid              1.881354e-01 2.291778e-01
# Fumaric_Acid                2.348875e-02 1.434716e-01
# Isobutyric_Acid             3.727377e-02 1.685307e-03
# Hippuric_Acid               5.974799e-03 1.455400e-03
# Methylmalonic_Acid          9.172968e-03 1.213133e-02
# LYSOC14_0                   4.315117e-01 1.281619e-01
# LYSOC16_1                   2.680834e-01 3.301057e-01
# LYSOC16_0                   2.410465e-01 2.864595e-01
# LYSOC17_0                   4.526427e-01 1.490811e-01
# LYSOC18_2                   2.076737e-01 2.398970e-02
# LYSOC18_1                   1.487173e-01 2.795595e-02
# LYSOC18_0                   7.001372e-02 2.440931e-04
# LYSOC20_4                   4.210850e-02 2.319836e-02
# LYSOC20_3                   2.699087e-01 5.750878e-02
# LYSOC24_0                   1.694073e-01 1.742516e-02
# LYSOC26_1                   2.806628e-01 4.742970e-02
# LYSOC26_0                   3.037176e-01 9.876832e-03
# LYSOC28_1                   2.867030e-01 1.849699e-02
# LYSOC28_0                   2.383592e-01 6.913765e-02
# X14_1SMOH                   1.867900e-01 1.603281e-04
# X16_1SM                     4.612976e-04 4.445196e-03
# X16_0SM                     8.820118e-03 4.249232e-03
# X16_1SMOH                   5.177943e-02 5.973224e-02
# X18_1SM                     1.463568e-01 2.376754e-05
# PC32_2AA                    1.703498e-01 3.033040e-03
# X18_0SM                     1.013140e-01 2.499897e-03
# X20_2SM                     2.458338e-01 6.763740e-03
# PC36_0AE                    9.154381e-02 9.031817e-02
# PC36_6AA                    4.789165e-01 1.313780e-02
# PC36_0AA                    7.213911e-02 8.472776e-02
# X22_2SMOH                   5.448554e-02 6.121429e-03
# X22_1SMOH                   1.438095e-01 8.064087e-03
# PC38_6AA                    3.598842e-01 3.582327e-03
# PC38_0AA                    1.641423e-01 2.254207e-02
# PC40_6AE                    3.003767e-01 4.196881e-04
# X24_1SMOH                   2.812321e-01 8.166995e-03
# PC40_6AA                    7.296261e-02 3.822779e-02
# PC40_2AA                    3.886590e-01 6.548867e-04
# PC401AA                     2.526609e-01 7.047846e-03
# C2                          9.009089e-03 2.035836e-01
# C3_1                        1.344427e-03 6.520835e-03
# C3                          6.928645e-04 2.156328e-05
# C4_1                        1.110670e-01 1.734466e-02
# C4                          6.433119e-02 2.403880e-04
# C3OH                        1.247321e-01 4.109266e-03
# C5_1                        8.047962e-04 8.170401e-02
# C5                          1.113929e-02 8.057001e-05
# C4OH                        1.614726e-02 3.801353e-03
# C6_1                        1.700405e-02 1.521249e-02
# C6                          2.988357e-02 1.003040e-02
# C5OH                        4.178037e-02 3.820190e-02
# C5_1DC                      2.075661e-02 5.592239e-04
# C5DC                        8.319454e-03 1.396790e-01
# C8                          6.279146e-02 1.137866e-02
# C5MDC                       6.320118e-02 3.351641e-04
# C9                          5.686592e-02 2.877191e-04
# C7DC                        1.911445e-02 4.392977e-02
# C10_2                       2.774464e-02 2.669824e-02
# C10_1                       9.806875e-05 5.181620e-03
# C10                         9.091663e-02 6.177194e-03
# C12_1                       5.561083e-02 1.129059e-02
# C12                         8.724983e-02 1.021427e-03
# C14_2                       1.489329e-02 8.214272e-03
# C14_1                       4.682379e-02 1.237084e-03
# C14                         8.287300e-02 1.100408e-03
# C12DC                       2.789348e-02 1.517156e-04
# C14_2OH                     6.238799e-02 7.064134e-04
# C14_1OH                     7.434849e-02 8.288561e-06
# C16_2                       8.538644e-02 1.251675e-03
# C16_1                       8.173804e-02 6.013786e-05
# C16                         4.486426e-02 4.096743e-04
# C16_2OH                     4.320782e-02 6.149230e-06
# C16_1OH                     2.490179e-02 1.673106e-04
# C16OH                       8.114658e-02 4.899521e-03
# C18_2                       6.559472e-02 1.418296e-03
# C18_1                       5.148796e-02 2.349343e-04
# C18                         1.764073e-02 1.471310e-04
# C18_1OH                     6.642490e-02 3.992356e-04
# Arginine_Average            5.353030e-05 6.440905e-02
# Betaine_Average             3.127252e-02 3.045666e-01
# C0_Average                  7.791831e-03 7.666972e-02
# Choline_Average             8.770000e-02 2.180159e-02
# Citrate_Average             6.629556e-02 2.406152e-02
# Creatine_Average            3.167248e-03 7.305680e-03
# Creatinine_Average          6.760689e-03 4.358452e-02
# Glucose_Average             1.938546e-01 8.717586e-03
# Glutamate_Average           2.053242e-02 3.535388e-02
# Glycine_Average             7.954284e-02 6.985749e-04
# Histidine_Average           1.433191e-02 2.135725e-02
# Isoleucine_Average          3.491738e-02 9.587265e-03
# Lactate_Average             2.498099e-06 1.291438e-03
# Leucine                     3.671511e-02 8.805579e-04
# Methionine                  4.778661e-02 7.732237e-05
# Phenylalanine_Average       2.568186e-02 2.753033e-03
# Proline                     6.912427e-03 1.643523e-02
# Pyruvate_Average            3.474870e-04 7.527999e-02
# Succinate_Average           8.094191e-03 2.435862e-01
# Threonine_Average           7.205669e-03 1.195676e-01
# Tyrosine_Average            6.617092e-02 3.454064e-02
# Valine_Average              3.860140e-02 2.277026e-02
# X1_Methylhistidine          1.428568e-02 5.634262e-02
# X2_Hydroxybutyrate          1.311877e-01 8.276223e-02
# X2_Hydroxyisovaleric_Acid   5.030433e-02 3.963719e-02
# X3_Hydroxybutyrate          1.999881e-01 1.222534e-04
# ADP                         3.866709e-03 9.950487e-02
# AMP                         1.257454e-02 1.446705e-03
# ATP                         1.610110e-02 4.416966e-02
# Acetamide                   1.207675e-03 3.236160e-02
# Acetate                     2.046768e-01 4.606231e-02
# Acetoacetate                5.414249e-02 1.096691e-01
# Acetone                     1.931490e-01 4.683060e-02
# Adenosine                   9.750632e-02 1.599287e-01
# Alanine                     2.109037e-02 2.080048e-02
# Aspartate                   8.933482e-02 9.463133e-03
# Creatine_Phosphate          7.359022e-05 1.464840e-02
# Cytidine                    5.760281e-02 4.068799e-04
# Dimethyl_Sulfone            9.768818e-02 4.920803e-02
# Ethanol                     1.689393e-04 1.038569e-01
# Ethanolamine                1.086963e-02 2.389164e-01
# Formate                     1.185261e-02 1.128993e-01
# Glycerol                    1.364540e-03 1.538217e-03
# Guanidoacetate              2.325195e-02 6.988681e-02
# Guanosine                   3.732801e-02 3.269576e-01
# Hypoxanthine                2.138509e-04 2.072337e-01
# IMP                         6.609714e-02 4.139632e-03
# Inosine                     1.591016e-02 2.097885e-03
# Isopropanol                 1.830630e-02 1.514116e-03
# Malate                      1.409989e-01 4.667972e-03
# Malonate                    1.014134e-01 3.896229e-02
# Mannose                     6.930146e-02 1.143824e-01
# Methanol                    2.323648e-02 2.688542e-01
# N_N_Dimethylglycine         3.669702e-02 3.072245e-03
# Nicotinurate                1.612843e-01 1.533280e-01
# O_Acetylcarnitine           1.283592e-02 2.477089e-01
# O_Phosphocholine            4.037472e-02 7.213647e-02
# Propylene_Glycol            2.010754e-01 2.824560e-03
# Taurine                     9.088339e-02 6.210161e-02
# Uridine                     1.468043e-01 2.501432e-02
# Xanthine                    5.134591e-02 2.039754e-02
# sn_Glycero_3_Phosphocholine 1.070794e-02 5.768241e-02
# Beta.Alanine                1.046778e-02 7.997229e-03



#### graph of variables default plot ####
fviz_pca_var(DM_pca_FM, col.var = "black")
pca18_biplot_black<- fviz_pca_var(DM_pca_FM, col.var = "black")
pca18_biplot_black
# Open the graphics device
png("pca18_biplot_black")
dev.off()
# plot aka variable correlation plot
# positively correlated: variables grouped together 
# negatively correlated: variables positioned on opposite sides of the origin (opposed quadrants)
# distance between variables and the origin measures the quality of the variables on the factor map. 
# variables that are away from the origin are well represented on the factor map
fviz_pca_var(DM_pca_FM, col.var = "black", axes = c(2,3))
pca18_biplot_black2_3<- fviz_pca_var(DM_pca_FM, col.var = "black", axes= c(2,3))
pca18_biplot_black2_3
# Open the graphics device
png("pca18_biplot_black2_3")
dev.off()

fviz_pca_var(DM_pca_FM, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4<- fviz_pca_var(DM_pca_FM, col.var = "black", axes = c(3,4))
pca18_biplot_black3_4
png("pca18_biplot_black3_4")
dev.off()

fviz_pca_var(DM_pca_FM, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5<- fviz_pca_var(DM_pca_FM, col.var = "black", axes = c(4,5))
pca18_biplot_black4_5
png("pca18_biplot_black4_5")
dev.off()


#### Quality of representation ####
#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.
# visualize cos2 as bar graph
fviz_cos2(DM_pca_FM, choice = "var", axes = 1:2)
bargraph18_cos2<- fviz_cos2(DM_pca_FM, choice = "var", axes = 1:2)
bargraph18_cos2
png("bargraph18_cos2")
dev.off()
# Color by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2<-fviz_pca_var(DM_pca_FM, col.var = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                   repel = TRUE # Avoid text overlapping
)
biplot18_colour_cos2
png("biplot18_colour_cos2")
dev.off()
# Adjust transparancy by cos2 values: quality on the factor map
fviz_pca_var(DM_pca_FM, alpha.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2<- fviz_pca_var(DM_pca_FM, alpha.var = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                                   repel = TRUE # Avoid text overlapping
)
biplot18_black_cos2
png("biplot18_black_cos2")
dev.off()

#### Contributions of variables to PCs ####

# contributions of variables in a given principal component are expressed in percentage.
# Variables that are correlated with PC1 and PC2 are the most important in explaining the variability in the data set.
# Variables that do not correlate with any PC or correlate with the last dimensions are variables with 
# low contribution and might be removed to simplify the overall analysis.

head(var$contrib) # contributions of variables to the PCs
# The larger the value of the contribution, the more the variable contributes to the component.

# visualize contributions as corrplot

corrplot(var$contrib, is.corr=FALSE,
         tl.col = "black", tl.cex = 1, 
         cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
         mar = c(1, 1, 1, 1)
) 
corrplot18_contrib <- corrplot(var$contrib, is.corr=FALSE,
                               tl.col = "black", tl.cex = 1, 
                               cl.cex = 1, cl.align.text = "l", cl.ratio = 0.75,
                               mar = c(1, 1, 1, 1)
)
corrplot18_contrib
pdf("mcorrplot18_contrib", height = 7, width =5)
dev.off()
# Barplot of contributions to a component
# The red dashed line on the graph above indicates the expected average contribution. 
# If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/32 = 3.125%.
# For a given component, a variable with a contribution larger than this cutoff 
# could be considered as important in contributing to the component.

#### contributions of variables ####
# to PC1
contrib18_PC1 <- fviz_contrib(DM_pca_FM, choice = "var", axes = 1, top = 80,
                            fill = "#6baed6", color = "#2171b5",
                            title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 19, label = expression(bold("PC1")), size = 5)
contrib18_PC1
png("contrib18_PC1")
dev.off()

# contributions of variables to PC2
contrib18_PC2 <- fviz_contrib(DM_pca_FM, choice = "var", axes = 2, top = 80,
                            fill = "#6baed6", color = "#2171b5",
                            title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 8.5, x = 19, label = expression(bold("PC2")), size = 5)
contrib18_PC2
png("contrib18_PC2")
dev.off()
# contributions of variables to PC3
contrib18_PC3 <- fviz_contrib(DM_pca_FM, choice = "var", axes = 3, top = 40,
                            fill = "#6baed6", color = "#2171b5",
                            title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 10, x = 19, label = expression(bold("PC3")), size = 5)
contrib18_PC3
png("contrib18_PC3")
dev.off()

# contributions of variables to PC4
contrib18_PC4 <- fviz_contrib(DM_pca_FM, choice = "var", axes = 4, top = 40,
                            fill = "#6baed6", color = "#2171b5",
                            title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 15, x = 19, label = expression(bold("PC4")), size = 5)
contrib18_PC4
png("contrib18_PC4")
dev.off()
# contributions of variables to PC5
contrib18_PC5 <- fviz_contrib(DM_pca_FM, choice = "var", axes = 5, top = 40,
                            fill = "#6baed6", color = "#2171b5",
                            title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 25, x = 19, label = expression(bold("PC5")), size = 5)
contrib18_PC5
png("contrib18_PC5")
dev.off()
# contributions of variables to PC1 and PC2
contrib18_PC1_PC2 <- fviz_contrib(DM_pca_FM, choice = "var", axes = c(1,2), top = 80,
                                fill = "#6baed6", color = "#2171b5",
                                title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 5, x = 18, label = expression(bold("PC1 & PC2")), size = 5)
contrib18_PC1_PC2
png("contrib18_PC1_PC2")
dev.off()
# contributions of variables to PC2 and PC3
contrib18_PC2_PC3 <- fviz_contrib(DM_pca_FM, choice = "var", axes = c(2,3), top = 40,
                                fill = "#6baed6", color = "#2171b5",
                                title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 6.5, x = 18, label = expression(bold("PC2 & PC3")), size = 5)
contrib18_PC2_PC3
png("contrib18_PC2_PC3")
dev.off()
# contributions of variables to PC3 and PC4
contrib18_PC3_PC4 <- fviz_contrib(DM_pca_FM, choice = "var", axes = c(3,4), top = 40,
                                fill = "#6baed6", color = "#2171b5",
                                title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 7.5, x = 18, label = expression(bold("PC3 & PC4")), size = 5)
contrib18_PC3_PC4
png("contrib18_PC3_PC4")
dev.off()
# contributions of variables to PC4 and PC5
contrib18_PC4_PC5 <- fviz_contrib(DM_pca_FM, choice = "var", axes = c(4,5), top = 40,
                                fill = "#6baed6", color = "#2171b5",
                                title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 13, x = 18, label = expression(bold("PC4 & PC5")), size = 5)
contrib18_PC4_PC5
png("contrib18_PC4_PC5")
dev.off()
# contributions of variables to PC1:PC3
contrib18_PC1_to_PC3 <- fviz_contrib(DM_pca_FM, choice = "var", axes = c(1:3), top = 40,
                                   fill = "#6baed6", color = "#2171b5",
                                   title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold")) +
  theme(aspect.ratio = 1/2) +
  xlab("Metabolite") +
  annotate("text", y = 4.5, x = 18, label = expression(bold("PC1 to PC3")), size = 5)
contrib18_PC1_to_PC3
png("contrib18_PC1_PC3")
dev.off()
#### colour variable colors using their contributions ####

# PC1 and PC2
varplot18_contrib_PC1_PC2 <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE,  # avoid text overlapping
                                        title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2
pdf("varplot18_contrib_PC1_PC2")
dev.off()
png("varplot18_contrib_PC1_PC2")
dev.off()

# PC1 and PC2 select var
varplot18_contrib_PC1_PC2_select <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,  # avoid text overlapping
                                               title = "",
                                               select.var = list(contrib = 60) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))

varplot18_contrib_PC1_PC2_select
pdf("varplot18_contrib_PC1_PC2_select")
dev.off()
png("varplot18_contrib_PC1_PC2_select")
dev.off()

# PC2 and PC3
varplot18_contrib_PC2_PC3 <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                        axes = c(2,3),
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE,  # avoid text overlapping
                                        title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3
pdf("varplot18_contrib_PC2_PC3")
dev.off()
png("varplot18_contrib_PC2_PC3")
dev.off()

# PC2 and PC3 select var
varplot18_contrib_PC2_PC3_select <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                               axes = c(2,3),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,  # avoid text overlapping
                                               title = "",
                                               select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC2_PC3_select
pdf("varplot18_contrib_PC2_PC3_select")
dev.off()
png("varplot18_contrib_PC2_PC3_select")
dev.off()
# PC4 and PC3
varplot18_contrib_PC4_PC3 <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                        axes = c(3,4),
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE,  # avoid text overlapping
                                        title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3
pdf("varplot18_contrib_PC4_PC3")
dev.off()
png("varplot18_contrib_PC4_PC3")
dev.off()

# PC4 and PC3 select var
varplot18_contrib_PC4_PC3_select <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                               axes = c(3,4),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,  # avoid text overlapping
                                               title = "",
                                               select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC3_select
pdf("varplot18_contrib_PC4_PC3_select")
dev.off()
png("varplot18_contrib_PC4_PC3_select")
dev.off()


# PC4 and PC5
varplot18_contrib_PC4_PC5 <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                        axes = c(4,5),
                                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                        repel = TRUE,  # avoid text overlapping
                                        title = ""
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5
pdf("varplot18_contrib_PC4_PC5")
dev.off()
png("varplot18_contrib_PC4_PC5")
dev.off()

# PC4 and PC5 select var
varplot18_contrib_PC4_PC5_select <- fviz_pca_var(DM_pca_FM, col.var = "contrib",  # colour by contributions to the PC
                                               axes = c(4,5),
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE,  # avoid text overlapping
                                               title = "",
                                               select.var = list(contrib = 30) # keeps only top var
) +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14, face = "bold"))
varplot18_contrib_PC4_PC5_select
pdf("varplot18_contrib_PC4_PC5_select")
dev.off()
png("varplot18_contrib_PC4_PC5_select")
dev.off()

# install the colour brewer
library("RColorBrewer")

#### cos2 selected variance trying plot ####

fviz_pca_var(DM_pca_FM, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE, # Avoid text overlapping
             select.var = list(cos2 = 50)
)
cos2_18_select <- fviz_pca_var(DM_pca_FM, col.var = "cos2",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE, # Avoid text overlapping
                               select.var = list(cos2 = 50)
)
pdf("cos2_18_select")
dev.off()
png("cos2_18_select")
dev.off()



#### colour variables by groups ####
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

dim.desc <- dimdesc(DM_pca_FM, axes = c(1:5))

dim.desc$Dim.1
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# C14_1                         0.8201880 3.253969e-14
# X22_1SMOH                     0.8106064 1.106268e-13
# C4                            0.7916023 1.033169e-12
# C12_1                         0.7704083 9.650477e-12
# C16OH                         0.7514440 5.878758e-11
# Beta_Hydroxybutyric_Acid      0.7470116 8.763312e-11
# C9                            0.7419646 1.367126e-10
# C14_2                         0.7408389 1.507577e-10
# C14                           0.7389158 1.779663e-10
# PC38_0AA                      0.7342564 2.644533e-10
# C8                            0.7329648 2.947141e-10
# X3_Hydroxybutyrate            0.7313472 3.372480e-10
# C18_2                         0.7307245 3.551191e-10
# X22_2SMOH                     0.7288230 4.154085e-10
# C14_1OH                       0.7261232 5.178437e-10
# X16_1SM                       0.7169647 1.073236e-09
# C16_1OH                       0.7155826 1.195063e-09
# PC40_6AE                      0.7098101 1.859918e-09
# C16_2OH                       0.7018928 3.354026e-09
# C18_1OH                       0.6993482 4.037511e-09
# C18_1                         0.6966551 4.902992e-09
# X16_0SM                       0.6958427 5.196629e-09
# C5_1DC                        0.6906949 7.480007e-09
# C14_2OH                       0.6903733 7.650270e-09
# C12                           0.6877572 9.177488e-09
# X18_1SM                       0.6867498 9.838853e-09
# C16_1                         0.6828226 1.287113e-08
# Citrate_Average               0.6809853 1.457428e-08
# C16                           0.6793946 1.621817e-08
# PC401AA                       0.6653514 4.051867e-08
# PC38_6AA                      0.6636877 4.501642e-08
# C6                            0.6618472 5.053722e-08
# C10                           0.6603984 5.532445e-08
# C18                           0.6578760 6.469108e-08
# PC40_2AA                      0.6496039 1.069544e-07
# LYSOC20_4                     0.6476246 1.203563e-07
# C16_2                         0.6441053 1.481554e-07
# Total_Dimethylarginine        0.6430979 1.571592e-07
# X14_1SMOH                     0.6331418 2.783659e-07
# Tryptophan                    0.6308594 3.164428e-07
# LYSOC28_1                     0.6145799 7.668636e-07
# Lactate_Average               0.6135054 8.115777e-07
# PC36_0AE                      0.6091197 1.020589e-06
# PC36_6AA                      0.6063220 1.179116e-06
# LYSOC18_2                     0.5945655 2.131193e-06
# C3                            0.5944766 2.140558e-06
# X20_2SM                       0.5921261 2.402668e-06
# PC40_6AA                      0.5884019 2.879864e-06
# C12DC                         0.5817953 3.949644e-06
# C5DC                          0.5788485 4.537189e-06
# X24_1SMOH                     0.5740479 5.670946e-06
# LYSOC18_1                     0.5727274 6.026054e-06
# C4OH                          0.5660168 8.172141e-06
# C0_Average                    0.5632399 9.252092e-06
# Alpha_Ketoglutaric_Acid       0.5566267 1.237879e-05
# Spermine                      0.5543752 1.364932e-05
# Spermidine                    0.5476741 1.818035e-05
# Valine_Average                0.5459043 1.959004e-05
# Asymmetric_Dimethylarginine   0.5437495 2.144267e-05
# C4_1                          0.5373410 2.795184e-05
# LYSOC18_0                     0.5311311 3.595537e-05
# Putrescine                    0.5302696 3.721896e-05
# Creatine_Average              0.5179974 6.027052e-05
# Isoleucine_Average            0.5141470 6.984631e-05
# Histidine_Average             0.5034163 1.043780e-04
# Acetone                       0.5026999 1.071644e-04
# Creatinine_Average            0.5018302 1.106385e-04
# C2                            0.4991459 1.220220e-04
# Butyric_acid                  0.4957849 1.377837e-04
# C3_1                          0.4944874 1.443503e-04
# C5                            0.4897598 1.707714e-04
# Fumaric_Acid                  0.4870625 1.877537e-04
# Threonine_Average             0.4837491 2.107188e-04
# Glycine_Average               0.4791644 2.467222e-04
# C5OH                          0.4783319 2.538314e-04
# X18_0SM                       0.4730653 3.032931e-04
# C5_1                          0.4645470 4.020797e-04
# LYSOC26_1                     0.4622373 4.334800e-04
# Phenylalanine_Average         0.4516107 6.085420e-04
# C5MDC                         0.4496816 6.464381e-04
# Glutamine                     0.4484044 6.726855e-04
# Acetyl_Ornithine              0.4453616 7.391363e-04
# Glutamate_Average             0.4424847 8.073479e-04
# O_Acetylcarnitine             0.4381890 9.197755e-04
# Proline                       0.4364906 9.679785e-04
# Acetate                       0.4231940 1.431008e-03
# Isobutyric_Acid               0.4195623 1.587939e-03
# PC36_0AA                      0.4114987 1.992563e-03
# X16_1SMOH                     0.4110880 2.015434e-03
# C10_1                         0.4101403 2.069103e-03
# Uridine                       0.4080027 2.194858e-03
# C3OH                          0.4074219 2.230183e-03
# PC32_2AA                      0.3971957 2.940598e-03
# Methionine                    0.3947596 3.136884e-03
# C10_2                         0.3897866 3.573898e-03
# LYSOC14_0                     0.3814819 4.424272e-03
# Kynurenine                    0.3807518 4.506914e-03
# C6_1                          0.3759415 5.086171e-03
# Nicotinurate                  0.3745108 5.270608e-03
# Alpha_Aminoadipic_Acid        0.3586028 7.751198e-03
# Tyrosine_Average              0.3414637 1.150487e-02
# Serine                        0.3388464 1.219777e-02
# LYSOC28_0                     0.3326021 1.399735e-02
# LYSOC16_1                     0.3186805 1.884340e-02
# Leucine                       0.3122922 2.150488e-02
# Sarcosine                     0.3071637 2.386531e-02
# Trans_Hydroxyproline          0.3059434 2.445801e-02
# Asparagine                    0.3038075 2.552520e-02
# LYSOC17_0                     0.2976530 2.882120e-02
# LYSOC20_3                     0.2894597 3.375432e-02
# Succinate_Average             0.2830452 3.808858e-02
# Dimethyl_Sulfone              0.2744669 4.459325e-02
# X2_Hydroxybutyrate            0.2721621 4.648802e-02
# IMP                           0.2714799 4.706131e-02
# Xanthine                     -0.3014629 2.674143e-02
# ADP                          -0.3750485 5.200606e-03
# Acetamide                    -0.4140893 1.853550e-03
dim.desc$Dim.2

# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# C12DC                         0.6855649 1.067421e-08
# C16_2                         0.6409279 1.783250e-07
# C14_2OH                       0.6310981 3.122442e-07
# C10_2                         0.6146910 7.623715e-07
# C18_1                         0.6140602 7.881942e-07
# C16_1                         0.6126920 8.470346e-07
# C18                           0.5979459 1.801994e-06
# C18_2                         0.5847884 3.425945e-06
# C16_2OH                       0.5608283 1.029587e-05
# C16                           0.5500985 1.640104e-05
# C12                           0.5486462 1.744645e-05
# C18_1OH                       0.5485659 1.750604e-05
# Ethanolamine                  0.5389669 2.614700e-05
# C10                           0.5348487 3.094265e-05
# Acetamide                     0.5294524 3.845536e-05
# C16_1OH                       0.5196280 5.659173e-05
# C14_1OH                       0.4891200 1.746680e-04
# C14                           0.4797589 2.417576e-04
# ADP                           0.4784334 2.529547e-04
# C9                            0.4782717 2.543520e-04
# Hypoxanthine                  0.4640300 4.089229e-04
# C14_2                         0.4588685 4.832692e-04
# C14_1                         0.4338776 1.046580e-03
# LYSOC24_0                     0.4292578 1.199682e-03
# C8                            0.4279972 1.244806e-03
# C16OH                         0.4204216 1.549482e-03
# PC32_2AA                      0.4165011 1.731975e-03
# N_N_Dimethylglycine           0.4050526 2.379568e-03
# LYSOC26_0                     0.3956970 3.060019e-03
# X18_0SM                       0.3921227 3.362332e-03
# sn_Glycero_3_Phosphocholine   0.3902570 3.530379e-03
# Choline_Average               0.3861413 3.927553e-03
# C5MDC                         0.3833098 4.223182e-03
# Pyruvate_Average              0.3589106 7.694940e-03
# Guanosine                     0.3454734 1.050916e-02
# C6                            0.3400685 1.186989e-02
# LYSOC17_0                     0.3398759 1.192106e-02
# Adenosine                     0.3197266 1.843541e-02
# LYSOC16_0                     0.3142170 2.067147e-02
# Isopropanol                   0.2876024 3.496491e-02
# Lysine                        0.2871536 3.526279e-02
# C0_Average                    0.2835352 3.774208e-02
# C2                            0.2809680 3.958700e-02
# Methylhistidine               0.2759798 4.338425e-02
# Methanol                     -0.2704743 4.791686e-02
# C3_1                         -0.2776340 4.209319e-02
# C4                           -0.2919006 3.221635e-02
# Propylene_Glycol             -0.2927906 3.167025e-02
# Histidine_Average            -0.2931094 3.147654e-02
# X22_2SMOH                    -0.2944640 3.066439e-02
# Cytidine                     -0.3089267 2.303043e-02
# PC38_0AA                     -0.3167873 1.960130e-02
# Alpha_Aminoadipic_Acid       -0.3194118 1.855737e-02
# Uridine                      -0.3298225 1.486886e-02
# Leucine                      -0.3310391 1.448184e-02
# Butyric_acid                 -0.3393116 1.207206e-02
# Creatine_Average             -0.3434798 1.099451e-02
# Citrate_Average              -0.3525476 8.933169e-03
# Xanthine                     -0.3612490 7.278965e-03
# C4_1                         -0.3947879 3.134539e-03
# Lactate_Average              -0.4157303 1.770031e-03
# O_Phosphocholine             -0.4290136 1.208309e-03
# Spermine                     -0.4343337 1.032460e-03
# Total_Dimethylarginine       -0.4392563 8.906020e-04
# X16_1SM                      -0.4480244 6.806792e-04
# LYSOC20_4                    -0.4498801 6.424430e-04
# Creatinine_Average           -0.4538630 5.668416e-04
# IMP                          -0.4605100 4.583953e-04
# Hippuric_Acid                -0.4651227 3.945810e-04
# Dimethyl_Sulfone             -0.4714872 3.197342e-04
# X16_0SM                      -0.4813834 2.286501e-04
# Putrescine                   -0.4867777 1.896341e-04
# X2_Hydroxybutyrate           -0.4868825 1.889407e-04
# PC36_0AE                     -0.4983897 1.254170e-04
# C5_1                         -0.5082754 8.716232e-05
# Spermidine                   -0.5284207 4.007029e-05
# X16_1SMOH                    -0.5284521 4.002031e-05
# Tryptophan                   -0.5345513 3.131853e-05
# Asymmetric_Dimethylarginine  -0.5571181 1.211647e-05
# Inosine                      -0.5596950 1.082313e-05
# Isoleucine_Average           -0.5658942 8.217221e-06
# C3                           -0.5784021 4.632957e-06
# Aspartate                    -0.5824379 3.831302e-06
# C5                           -0.6207907 5.503466e-07
# Valine_Average               -0.6353091 2.462194e-07
# PC36_0AA                     -0.6637201 4.492447e-08
# Isobutyric_Acid              -0.6779555 1.785448e-08
# Proline                      -0.7350098 2.481879e-10
# Asparagine                   -0.7418040 1.386368e-10
# Trans_Hydroxyproline         -0.7766448 5.128591e-12
dim.desc$Dim.3
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Arginine_Average              0.8290413 9.832643e-15
# Serine                        0.8246486 1.795458e-14
# Ornithine                     0.7916074 1.032578e-12
# Methionine                    0.7625756 2.077098e-11
# Lysine                        0.7479601 8.051344e-11
# Methionine_Sulfoxide          0.7286657 4.208092e-10
# Phenylalanine_Average         0.7171813 1.055250e-09
# Sarcosine                     0.6590750 6.006711e-08
# Glutamine                     0.6160487 7.094533e-07
# Glycine_Average               0.6048020 1.274599e-06
# Glutamate_Average             0.5950765 2.078065e-06
# Pyruvate_Average              0.5753905 5.329867e-06
# Glucose_Average               0.5576682 1.182890e-05
# Leucine                       0.5356623 2.993539e-05
# Alpha_Ketoglutaric_Acid       0.5271201 4.219505e-05
# X18_0SM                       0.5258037 4.445095e-05
# LYSOC18_0                     0.5238065 4.808606e-05
# Succinate_Average             0.5157042 6.581692e-05
# Mannose                       0.5016896 1.112100e-04
# PC32_2AA                      0.5009090 1.144309e-04
# Fumaric_Acid                  0.4939183 1.473203e-04
# Histidine_Average             0.4822775 2.217159e-04
# Threonine_Average             0.4820159 2.237244e-04
# PC40_6AA                      0.4803417 2.369783e-04
# Choline_Average               0.4708878 3.261884e-04
# X20_2SM                       0.4618545 4.388934e-04
# Tyrosine_Average              0.4602277 4.625887e-04
# Kynurenine                    0.4426302 8.037654e-04
# Total_Dimethylarginine        0.3911838 3.445992e-03
# Uridine                       0.3646960 6.701527e-03
# X18_1SM                       0.3420362 1.135791e-02
# sn_Glycero_3_Phosphocholine   0.3376400 1.252905e-02
# Alanine                       0.3323755 1.406671e-02
# LYSOC18_1                     0.3150495 2.031955e-02
# Asymmetric_Dimethylarginine   0.3120536 2.161017e-02
# O_Phosphocholine              0.3088715 2.305616e-02
# Guanosine                     0.3085899 2.318800e-02
# Methylhistidine               0.3075049 2.370175e-02
# C5OH                          0.2831078 3.804420e-02
# Alpha_Aminoadipic_Acid        0.2725744 4.614430e-02
# Nicotinurate                  0.2702951 4.807063e-02
# C4OH                         -0.2687757 4.939080e-02
# PC38_0AA                     -0.2824709 3.849813e-02
# PC36_0AA                     -0.2836042 3.769345e-02
# C14                          -0.2863807 3.578071e-02
# C14_2                        -0.2902532 3.324779e-02
# C12_1                        -0.2957857 2.988889e-02
# Butyric_acid                 -0.2999488 2.755236e-02
# C18_1OH                      -0.3048204 2.501431e-02
# Propylene_Glycol             -0.3095573 2.273783e-02
# C18                          -0.3108264 2.215849e-02
# C14_1OH                      -0.3191867 1.864503e-02
# C16_2OH                      -0.3194914 1.852645e-02
# Isobutyric_Acid              -0.3205032 1.813737e-02
# Citrulline                   -0.3352805 1.319937e-02
# C6                           -0.3403006 1.180849e-02
# C16_1OH                      -0.3427797 1.116948e-02
# C16                          -0.3537986 8.676928e-03
# X16_0SM                      -0.3612247 7.283191e-03
# X16_1SM                      -0.3674186 6.274102e-03
# X14_1SMOH                    -0.4137953 1.868878e-03
# X16_1SMOH                    -0.4445931 7.568313e-04
# X22_2SMOH                    -0.4498134 6.437830e-04
# C6_1                         -0.6130381 8.317750e-07
dim.desc$Dim.4
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Acetate                    0.4524122 5.933951e-04
# Propylene_Glycol           0.4484144 6.724780e-04
# X3_Hydroxybutyrate         0.4472003 6.983084e-04
# Acetone                    0.4394872 8.843990e-04
# Beta_Hydroxybutyric_Acid   0.4260065 1.319168e-03
# Nicotinurate               0.4016022 2.613028e-03
# Uridine                    0.3831505 4.240383e-03
# Malate                     0.3754982 5.142705e-03
# X2_Hydroxybutyrate         0.3621984 7.115757e-03
# Malonate                   0.3184547 1.893244e-02
# Dimethyl_Sulfone           0.3125511 2.139117e-02
# Adenosine                  0.3122600 2.151907e-02
# C10                        0.3015239 2.670920e-02
# Taurine                    0.3014687 2.673834e-02
# Aspartate                  0.2988893 2.813199e-02
# C12                        0.2953808 3.012470e-02
# C16_2                      0.2922096 3.202588e-02
# C14                        0.2878767 3.478391e-02
# C16_1                      0.2858987 3.610694e-02
# C16OH                      0.2848624 3.681650e-02
# Glycine_Average            0.2820334 3.881250e-02
# C14_1OH                    0.2726692 4.606557e-02
# PC36_0AA                  -0.2685872 4.955660e-02
# PC40_6AA                  -0.2701159 4.822482e-02
# Choline_Average           -0.2961419 2.968274e-02
# PC36_0AE                  -0.3025621 2.616533e-02
# X18_0SM                   -0.3182986 1.899424e-02
# C4_1                      -0.3332671 1.379546e-02
# C3OH                      -0.3531744 8.803979e-03
# X22_1SMOH                 -0.3792222 4.684477e-03
# X18_1SM                   -0.3825661 4.304001e-03
# LYSOC18_1                 -0.3856389 3.978637e-03
# PC38_0AA                  -0.4051448 2.373594e-03
# LYSOC24_0                 -0.4115912 1.987445e-03
# PC32_2AA                  -0.4127345 1.925132e-03
# X14_1SMOH                 -0.4321921 1.100278e-03
# Propionic_Acid            -0.4337458 1.050694e-03
# Glucose_Average           -0.4402893 8.631590e-04
# LYSOC18_2                 -0.4557123 5.345546e-04
# LYSOC28_0                 -0.4882205 1.802829e-04
# LYSOC16_0                 -0.4909649 1.636473e-04
# X20_2SM                   -0.4958163 1.376281e-04
# PC401AA                   -0.5026538 1.073460e-04
# LYSOC16_1                 -0.5177678 6.080596e-05
# LYSOC20_3                 -0.5195274 5.681250e-05
# LYSOC26_1                 -0.5297762 3.796110e-05
# X24_1SMOH                 -0.5303132 3.715405e-05
# LYSOC28_1                 -0.5354466 3.019943e-05
# PC40_6AE                  -0.5480663 1.788084e-05
# LYSOC26_0                 -0.5511058 1.571029e-05
# PC38_6AA                  -0.5999035 1.633703e-06
# PC40_2AA                  -0.6234252 4.770557e-07
# LYSOC14_0                 -0.6568955 6.871916e-08
# LYSOC17_0                 -0.6727872 2.510535e-08
# PC36_6AA                  -0.6920379 6.806886e-09
dim.desc$Dim.5
# Link between the variable and the continuous variables (R-square)
# =================================================================================
#   correlation      p.value
# Citrulline           0.6838759 1.198126e-08
# LYSOC16_1            0.5745482 5.541534e-06
# Guanosine            0.5718021 6.287079e-06
# Betaine_Average      0.5518755 1.520067e-05
# LYSOC16_0            0.5352191 3.048018e-05
# Succinate_Average    0.4935445 1.493007e-04
# Ethanolamine         0.4887908 1.767041e-04
# Hypoxanthine         0.4552293 5.428225e-04
# Adenosine            0.3999108 2.734719e-03
# Nicotinurate         0.3915712 3.411254e-03
# LYSOC17_0            0.3861103 3.930688e-03
# Fumaric_Acid         0.3787765 4.737358e-03
# LYSOC14_0            0.3579970 7.862985e-03
# ADP                  0.3154439 2.015460e-02
# PC36_0AE             0.3005298 2.723877e-02
# PC36_0AA             0.2910803 3.272654e-02
# X2_Hydroxybutyrate   0.2876843 3.491083e-02
# O_Phosphocholine     0.2685823 4.956093e-02
# Pyruvate_Average    -0.2743720 4.466999e-02
# C0_Average          -0.2768930 4.266756e-02
# C5_1                -0.2858391 3.614738e-02
# Ethanol             -0.3222684 1.747521e-02
# Acetoacetate        -0.3311633 1.444284e-02
# Formate             -0.3360049 1.299037e-02
# Mannose             -0.3382047 1.237302e-02
# Threonine_Average   -0.3457855 1.043489e-02
# C5DC                -0.3737365 5.372839e-03
# C2                  -0.4512024 6.163901e-04
# Propionic_Acid      -0.4787252 2.504496e-04
# O_Acetylcarnitine   -0.4977036 1.285717e-04
# Methanol            -0.5185115 5.908758e-05


#### extract the results for individuals ####
ind <- get_pca_ind(DM_pca_FM)
ind
# Principal Component Analysis Results for individuals
#   Name       Description                       
# 1 "$coord"   "Coordinates for the individuals" 
# 2 "$cos2"    "Cos2 for the individuals"        
# 3 "$contrib" "contributions of the individuals"
ind$coord
# Dim.1       Dim.2        Dim.3       Dim.4       Dim.5
# 1    9.3597234   6.4343669  10.33650930  -1.6682317  1.05047424
# 2   -3.4658494   4.6736754  -1.73950973  -0.6469672  3.75342829
# 3    0.3886967   5.7581805  -2.68510745  -0.1579891  7.78138355
# 4   -1.9825596   4.3871951   2.89773708   0.1533919 -2.83573422
# 5    3.4341172   0.7950058  11.31872417   1.9131770 -0.84007669
# 6   -0.3675978   2.4299482   7.35405499   3.7660707 -0.95655026
# 7   -1.5194714   2.8890206   2.97035833  -0.5398773  1.05107010
# 8    6.5811372   6.6195810   7.45541012   0.1007786 -0.45126122
# 9    2.2601840   5.9124937   5.83123138   3.9648167  2.09401452
# 10   6.3631696   8.1218182  -0.81475719   5.7814811  0.75502499
# 11  -7.1914411   2.1559964   2.42349034   4.9598306  0.74373086
# 12  -5.4538344   4.2438933  -1.85634399  -3.3698526  4.03443984
# 13   0.6410306   4.0383973  -2.27631091  -6.4723060  5.38036876
# 14  -7.9006635   2.0284875  -2.46960430   1.9222890  2.19093246
# 15  -9.0217435   2.9877396  -3.94262702   3.3060306  2.28424782
# 16   0.9341614   3.2112990   4.84781856  -2.2416815 -0.15797757
# 17  -8.6263440   3.4910515  -4.36682762  -0.1473403  4.10268107
# 18  -3.7702748   4.1101454  -3.94763584  -3.5637018  1.44864091
# 19 -13.1429694   1.2636036  -5.29373746   3.1085571 -0.07873734
# 20   3.4772071   4.0079671  -5.23562966 -11.0746944  1.07537484
# 21  -7.0182356   1.4998758   3.86819356   1.4684151 -3.00365564
# 22   0.4101114   4.2908908   6.24529207  -0.7699629  0.30600221
# 23  -2.1728917  -1.5249763   1.44685963  -1.0205683 -0.24514552
# 24  -2.9030764   2.8715951   1.92835309  -0.5849704  0.01721363
# 25  -0.5800668  -3.0076724   0.45222617  -3.3970954  0.73593709
# 26  -3.8829291  -1.7012221  -0.90870014   1.1993902 -2.96055607
# 27   4.8227846  -1.1479427  -0.16375152  -9.9668929 -3.74781214
# 28  -0.5314745   0.2824234  -2.77784573  -6.1403862 -4.22295093
# 29  -7.3955172  -0.7659861  -2.73872659   2.9606386 -5.37523504
# 30  -1.8630852  -1.1870798  -0.27905477   0.2352490 -0.60254605
# 31  -1.4915949  -2.4854224  -0.68633094  -1.3175437 -5.06432088
# 32   0.5798846  -2.2400730  -0.21449338  -2.6155793 -2.67317318
# 33   3.1375164   1.4835035  -1.18194176  -0.4816026 -3.32579412
# 34  -1.5234373  -0.2968694   2.13229981  -1.1985514 -1.26669526
# 35  -0.6389409   2.6706736  -5.77559917   6.1661655 -6.35470493
# 36   2.5194141   3.3319382  -0.51101592  -2.0168390 -2.13638077
# 37  -4.7334816  -0.8704905  -1.38322381   2.5406587 -1.87003438
# 38  -1.6135292   0.5002786  -3.09672189   2.5190498 -1.69168447
# 39  -1.5917653   2.2546387  -1.72033853   0.7052074 -3.37297565
# 40   1.2089256  -8.0832386   2.33213663   1.8543022  0.67164964
# 41   3.0937311  -5.9245412   0.18768808  -2.6598220 -2.20902096
# 42   7.1266901  -8.1828277   0.46210815  -2.3022244  0.60132365
# 43  30.1734923   7.4762600 -11.98097831   6.7874882 -0.93938311
# 44   6.7430370  -8.4633699  -0.96512845  -2.2806689  1.87013591
# 45  -1.7929783  -2.2074550  -1.51549250   1.9526299  2.24063031
# 46   4.3943177  -7.7884692  -3.14556651   5.0067465  1.99650547
# 47  -3.9792604  -2.4020185  -2.34682374  -0.7162311 -0.90528890
# 48  -0.8134112  -5.0864664   0.02902364   0.3278962 -0.79158196
# 49  -1.8787366  -5.3870015  -0.60614395   2.5866944  1.94717765
# 50   3.4348180  -8.1408716  -2.07638745   1.1572383  3.47922622
# 51  -2.1720366  -3.1708537  -2.68837096  -1.8657552  0.07760855
# 52   3.6513870  -5.6987147   0.15577467  -4.0083844  0.27039710
# 53   2.7904038 -10.9939704   3.21114663   5.4209273  3.99681786
# 54   3.4932572  -9.4644106   3.50429079   1.3605992  2.12283974                                                                               
ind$cos2
# Dim.1        Dim.2        Dim.3        Dim.4        Dim.5
# 1  0.282559317 0.1335350397 3.446128e-01 8.976279e-03 3.559217e-03
# 2  0.122538887 0.2228289720 3.086799e-02 4.269913e-03 1.437179e-01
# 3  0.001001290 0.2197398884 4.778171e-02 1.654218e-04 4.012841e-01
# 4  0.045745103 0.2240092729 9.772599e-02 2.738399e-04 9.358864e-02
# 5  0.057005660 0.0030551205 6.192738e-01 1.769286e-02 3.411345e-03
# 6  0.001152803 0.0503736179 4.613844e-01 1.210002e-01 7.805931e-03
# 7  0.032682020 0.1181476518 1.248940e-01 4.125855e-03 1.563822e-02
# 8  0.221091004 0.2236815635 2.837346e-01 5.184485e-05 1.039502e-03
# 9  0.035553512 0.2432967445 2.366549e-01 1.094060e-01 3.051788e-02
# 10 0.150926780 0.2458814894 2.474433e-03 1.245941e-01 2.124917e-03
# 11 0.357460403 0.0321286246 4.059557e-02 1.700316e-01 3.823206e-03
# 12 0.210719770 0.1275940148 2.441284e-02 8.044947e-02 1.153102e-01
# 13 0.002780458 0.1103511337 3.506080e-02 2.834502e-01 1.958767e-01
# 14 0.515248289 0.0339651927 5.034359e-02 3.050189e-02 3.962300e-02
# 15 0.414166675 0.0454234035 7.909802e-02 5.561705e-02 2.655098e-02
# 16 0.009170195 0.1083667941 2.469605e-01 5.280591e-02 2.622561e-04
# 17 0.385568477 0.0631481283 9.880526e-02 1.124841e-04 8.721345e-02
# 18 0.132289048 0.1572143704 1.450281e-01 1.181900e-01 1.952987e-02
# 19 0.651713719 0.0060240915 1.057291e-01 3.645756e-02 2.339006e-05
# 20 0.050794127 0.0674839781 1.151570e-01 5.152475e-01 4.858163e-03
# 21 0.333776163 0.0152443898 1.013948e-01 1.461158e-02 6.113638e-02
# 22 0.001325900 0.1451449314 3.074769e-01 4.673546e-03 7.381694e-04
# 23 0.037867335 0.0186515433 1.678964e-02 8.353581e-03 4.819883e-04
# 24 0.106779046 0.1044757529 4.711316e-02 4.335478e-03 3.754171e-06
# 25 0.002188433 0.0588353552 1.330113e-03 7.505728e-02 3.522562e-03
# 26 0.121494351 0.0233216456 6.653932e-03 1.159199e-02 7.062911e-02
# 27 0.123533808 0.0069989123 1.424167e-04 5.276059e-01 7.460118e-02
# 28 0.002806790 0.0007925858 7.667632e-02 3.746592e-01 1.772055e-01
# 29 0.408592214 0.0043832352 5.603383e-02 6.548228e-02 2.158479e-01
# 30 0.044738679 0.0181625803 1.003683e-03 7.133019e-04 4.679487e-03
# 31 0.024441273 0.0678613439 5.174748e-03 1.907006e-02 2.817503e-01
# 32 0.004794266 0.0715423876 6.559441e-04 9.753822e-02 1.018810e-01
# 33 0.083610879 0.0186925229 1.186541e-02 1.970010e-03 9.394669e-02
# 34 0.030020792 0.0011399972 5.881249e-02 1.858174e-02 2.075474e-02
# 35 0.001821039 0.0318155546 1.487962e-01 1.696009e-01 1.801311e-01
# 36 0.070879975 0.1239704928 2.916039e-03 4.542207e-02 5.096614e-02
# 37 0.209473893 0.0070842934 1.788766e-02 6.034776e-02 3.269397e-02
# 38 0.027713057 0.0026641253 1.020787e-01 6.754672e-02 3.046278e-02
# 39 0.036979484 0.0741919299 4.319471e-02 7.258327e-03 1.660463e-01
# 40 0.010805637 0.4830833227 4.021233e-02 2.542215e-02 3.335315e-03
# 41 0.060587735 0.2221923082 2.229939e-04 4.478419e-02 3.089010e-02
# 42 0.272572035 0.3593456278 1.146022e-03 2.844467e-02 1.940538e-03
# 43 0.762152837 0.0467907174 1.201643e-01 3.856640e-02 7.387141e-04
# 44 0.201466605 0.3173793209 4.127264e-03 2.304711e-02 1.549667e-02
# 45 0.038213429 0.0579228187 2.730069e-02 4.532167e-02 5.967692e-02
# 46 0.122738889 0.3855702396 6.289224e-02 1.593348e-01 2.533613e-02
# 47 0.171946540 0.0626529334 5.980668e-02 5.570511e-03 8.899449e-03
# 48 0.008075868 0.3157921129 1.028188e-05 1.312326e-03 7.648226e-03
# 49 0.035061070 0.2882621196 3.649597e-03 6.646357e-02 3.766210e-02
# 50 0.072538779 0.4074783988 2.650816e-02 8.233952e-03 7.442659e-02
# 51 0.047258688 0.1007163044 7.239787e-02 3.487038e-02 6.033461e-05
# 52 0.112717880 0.2745559771 2.051496e-04 1.358363e-01 6.181319e-04
# 53 0.027839345 0.4321495077 3.686763e-02 1.050684e-01 5.711545e-02
# 54 0.052481246 0.3852389144 5.281330e-02 7.961656e-03 1.938106e-02
ind$contrib
#Dim.1       Dim.2        Dim.3        Dim.4        Dim.5
# 1   4.324657742 3.424201524 1.201114e+01  0.399235849 2.676298e-01
# 2   0.592986869 1.806611489 3.401654e-01  0.060045615 3.416802e+00
# 3   0.007458430 2.742320592 8.105127e-01  0.003580728 1.468511e+01
# 4   0.194034169 1.591921185 9.439618e-01  0.003375371 1.950269e+00
# 5   0.582178202 0.052274305 1.440229e+01  0.525081952 1.711597e-01
# 6   0.006670703 0.488361812 6.079813e+00  2.034668808 2.219112e-01
# 7   0.113975311 0.690317607 9.918686e-01  0.041812577 2.679335e-01
# 8   2.138097941 3.624170987 6.248554e+00  0.001456979 4.938779e-02
# 9   0.252181543 2.891272707 3.822584e+00  2.255085527 1.063464e+00
# 10  1.998815504 5.455746745 7.462647e-02  4.795070801 1.382567e-01
# 11  2.553039714 0.384453421 6.602662e-01  3.528990446 1.341514e-01
# 12  1.468350079 1.489623780 3.873944e-01  1.629065780 3.947572e+00
# 13  0.020285387 1.348856545 5.825050e-01  6.009452227 7.020817e+00
# 14  3.081433804 0.340323809 6.856322e-01  0.530095570 1.164184e+00
# 15  4.017969759 0.738300479 1.747462e+00  1.567943954 1.265465e+00
# 16  0.043079388 0.852921704 2.641978e+00  0.720882724 6.052773e-03
# 17  3.673493352 1.007999198 2.143723e+00  0.003114295 4.082245e+00
# 18  0.701732719 1.397211145 1.751905e+00  1.821878878 5.089617e-01
# 19  8.527322838 0.132059415 3.150368e+00  1.386227437 1.503575e-03
# 20  0.596879733 1.328605273 3.081586e+00 17.594636757 2.804680e-01
# 21  2.431541026 0.186062261 1.682104e+00  0.309324637 2.188083e+00
# 22  0.008302893 1.522799056 4.384721e+00  0.085046484 2.270977e-02
# 23  0.233078325 0.192341891 2.353365e-01  0.149417321 1.457510e-02
# 24  0.416047231 0.682015233 4.180321e-01  0.049089053 7.186353e-05
# 25  0.016610461 0.748184544 2.299047e-02  1.655511872 1.313545e-01
# 26  0.744294208 0.239369995 9.282772e-02  0.206365832 2.125740e+00
# 27  1.148210310 0.108990370 3.014444e-03 14.250706013 3.406585e+00
# 28  0.013944102 0.006597038 8.674665e-01  5.408890273 4.325095e+00
# 29  2.699994121 0.048527674 8.432063e-01  1.257440965 7.007426e+00
# 30  0.171352766 0.116548712 8.754178e-03  0.007939122 8.805298e-02
# 31  0.109831650 0.510914359 5.295459e-02  0.249027252 6.220224e+00
# 32  0.016600027 0.415022875 5.172068e-03  0.981414725 1.733077e+00
# 33  0.485956865 0.182022392 1.570466e-01  0.033273184 2.682590e+00
# 34  0.114571041 0.007289176 5.111319e-01  0.206077279 3.891418e-01
# 35  0.020153345 0.589914821 3.749994e+00  5.454402155 9.793870e+00
# 36  0.313346475 0.918209007 2.935657e-02  0.583524706 1.106932e+00
# 37  1.106081437 0.062672310 2.150906e-01  0.925996772 8.481306e-01
# 38  0.128522563 0.020700058 1.078055e+00  0.910312101 6.940684e-01
# 39  0.125078825 0.420437637 3.327088e-01  0.071342869 2.759243e+00
# 40  0.072148091 5.404039030 6.114268e-01  0.493262210 1.094080e-01
# 41  0.472488070 2.903067399 3.960134e-03  1.014896959 1.183486e+00
# 42  2.507272015 5.538019592 2.400622e-02  0.760347465 8.769607e-02
# 43 44.944532370 4.622920929 1.613693e+01  6.608988227 2.140174e-01
# 44  2.244588811 5.924262842 1.047145e-01  0.746176028 8.482227e-01
# 45  0.158699572 0.403024449 2.581928e-01  0.546961416 1.217598e+00
# 46  0.953253934 5.017089958 1.112331e+00  3.596068967 9.667286e-01
# 47  0.781682580 0.477199926 6.191522e-01  0.073590734 1.987641e-01
# 48  0.032662243 2.139833726 9.469785e-05  0.015423732 1.519691e-01
# 49  0.174243865 2.400169276 4.130361e-02  0.959858162 9.195487e-01
# 50  0.582415832 5.481374630 4.846779e-01  0.192115483 2.935817e+00
# 51  0.232894907 0.831572417 8.124841e-01  0.499374198 1.460773e-03
# 52  0.658175094 2.685972129 2.727911e-03  2.304918277 1.773240e-02
# 53  0.384379152 9.996703213 1.159196e+00  4.215644350 3.874291e+00
# 54  0.602402603 7.408577353 1.380501e+00  0.265568906 1.092944e+00



fviz_pca_ind(DM_pca_FM, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE  # avoid text overlapping, slow if many points
)
pca18_ind <- 
  fviz_pca_ind(DM_pca_FM, col.ind = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE  # avoid text overlapping, slow if many points
  )
pca18_ind
png("pca18_ind")
dev.off()
 # change point size according to cos2
fviz_pca_ind(DM_pca_FM, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point<- fviz_pca_ind(DM_pca_FM, pointsize = "cos2", 
                               pointshape = 21, fill = "#E7B800",
                               repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point
png("pca18_ind_point")
dev.off()
# to change boht pointsize and colour according to cos2
fviz_pca_ind(DM_pca_FM, col.ind = "cos2", pointsize = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_size<- fviz_pca_ind(DM_pca_FM, col.ind = "cos2", pointsize = "cos2",
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE # Avoid text overlapping (slow if many points)
)
pca18_ind_point_size
png("pca18_ind_point_size")
dev.off()
# To visualize the contribution of individuals to the first two principal components
fviz_contrib(DM_pca_FM, choice = "ind", axes = 1:2)
contr18_ind<- fviz_contrib(DM_pca_FM, choice = "ind", axes = 1:2)
contr18_ind
png("contr18_ind")
dev.off()
#### Color by a custom continuous variable ####
# Spine age

fviz_pca_ind(DM_pca_FM, col.ind = DM_data.frame$Spine.Age,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# trying to do the comparisson of age with the PCA loadins

# Get the scores of the first principal component
pc1_scores18 <- DM_pca_FM$ind$coord[, "Dim.1"]
# second dim
pc1_scores18PC2 <- DM_pca_FM$ind$coord[, "Dim.2"]

# Calculate the correlation
correlation18 <- cor(DM_data.frame$Spine.Age, pc1_scores18)

correlation18
# [1] 0.06323247
# cortests
correlation18test <- cor.test(DM_data.frame$Spine.Age, pc1_scores18)

correlation18test
# Pearson's product-moment correlation
# 
# data:  DM_data.frame$Spine.Age and pc1_scores18
# t = 0.45689, df = 52, p-value = 0.6497
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2080506  0.3254822
# sample estimates:
#        cor 
# 0.06323247 


correlation18PC2 <- cor(DM_data.frame$Spine.Age, pc1_scores18PC2)

correlation18PC2
#[1] 0.3049553
correlation18PC2test <- cor.test(DM_data.frame$Spine.Age, pc1_scores18PC2)

correlation18PC2test
# Pearson's product-moment correlation
# 
# data:  DM_data.frame$Spine.Age and pc1_scores18PC2
# t = 2.3091, df = 52, p-value = 0.02495
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.04050198 0.52948111
# sample estimates:
#       cor 
# 0.3049553 

# weight now
# Calculate the correlation
correlation18w <- cor(DM_data.frame$W, pc1_scores18)

correlation18w
#[1] 0.1613928
# Calculate the correlation
correlation18wPC2 <- cor(DM_data.frame$W, pc1_scores18PC2)

correlation18wPC2
#[1] 0.2850698
correlation18wPC2test <- cor.test(DM_data.frame$W, pc1_scores18PC2)

correlation18wPC2test

# Pearson's product-moment correlation
# 
# data:  DM_data.frame$W and pc1_scores18PC2
# t = 2.1447, df = 52, p-value = 0.03667
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.01873963 0.51362511
# sample estimates:
#       cor 
# 0.2850698 


# fl now
# Calculate the correlation
correlation18fl <- cor(DM_data.frame$FL, pc1_scores18)

correlation18fl
#[1] 0.1116018

correlation18flPC2 <- cor(DM_data.frame$FL, pc1_scores18PC2)

correlation18flPC2
#[1] 0.215295

correlation18flPC2test <- cor.test(DM_data.frame$FL, pc1_scores18PC2)

correlation18flPC2test
# Pearson's product-moment correlation
# 
# data:  DM_data.frame$FL and pc1_scores18PC2
# t = 1.5898, df = 52, p-value = 0.1179
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.05567505  0.45672644
# sample estimates:
#      cor 
# 0.215295 




fviz_pca_ind(DM_pca_FM, col.ind = DM_data.frame$Spine_Age, axes = c(2,3),
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Age")  

# W

fviz_pca_ind(DM_pca_FM, col.ind = DM_data.frame$W,
             geom.ind = "point", # show points only (nbut not "text"),
             pointsize = 2.5,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             legend.title = "Weight") 
#### Anova sex and site ####

# Perform the Two-Way ANOVA
model1 <- aov(PC_scores ~ Sex * Other_Variable, data = my_data)

# Print the summary of the model
summary(model)

# Get the scores of the first principal component
pc1_scores18 <- DM_pca_FM$ind$coord[, "Dim.1"]

# Perform the ANOVA
modelss <- aov(pc1_scores18 ~ DM_data.frame$Sex * DM_data.frame$Site, data = DM_data.frame)

# Print the summary of the model
summary(modelss)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# DM_data.frame$Sex                     2  425.9  212.93   6.786 0.00266 **
#   DM_data.frame$Site                    3  171.9   57.30   1.826 0.15597   
# DM_data.frame$Sex:DM_data.frame$Site  3   16.0    5.33   0.170 0.91607   
# Residuals                            45 1411.9   31.38                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# only site model anova
# Perform the ANOVA
models <- aov(pc1_scores18 ~ DM_data.frame$Site, data = DM_data.frame)

# Print the summary of the model
summary(models)
# Df Sum Sq Mean Sq F value Pr(>F)  
# DM_data.frame$Site  3  343.5  114.50   3.403 0.0246 *
#   Residuals          50 1682.2   33.64                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(models)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18 ~ DM_data.frame$Site, data = DM_data.frame)
# 
# $`DM_data.frame$Site`
# diff        lwr        upr     p adj
# Matheson-Dauphin River  -4.6839906 -10.312729  0.9447477 0.1341048
# Red River-Dauphin River -4.2232493  -9.851988  1.4054891 0.2039547
# Sandy Bar-Dauphin River -7.3440549 -13.843562 -0.8445477 0.0210171
# Red River-Matheson       0.4607414  -5.167997  6.0894797 0.9963171
# Sandy Bar-Matheson      -2.6600643  -9.159571  3.8394429 0.6986087
# Sandy Bar-Red River     -3.1208056  -9.620313  3.3787015 0.5822174
library(report)
report(models)

# pc 2
# Perform the ANOVA
modelsPC2 <- aov(pc1_scores18PC2 ~ DM_data.frame$Site, data = DM_data.frame)

# Print the summary of the model
summary(modelsPC2)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# DM_data.frame$Site  3  794.2   264.7   31.91 1.14e-11 ***
#   Residuals          50  414.9     8.3                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
report(modelsPC2)
TukeyHSD(modelsPC2)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2 ~ DM_data.frame$Site, data = DM_data.frame)
# 
# $`DM_data.frame$Site`
# diff        lwr       upr     p adj
# Matheson-Dauphin River   5.355910  2.5606638  8.151156 0.0000314
# Red River-Dauphin River  9.799583  7.0043373 12.594829 0.0000000
# Sandy Bar-Dauphin River  8.148025  4.9203526 11.375697 0.0000001
# Red River-Matheson       4.443673  1.6484275  7.238919 0.0005704
# Sandy Bar-Matheson       2.792115 -0.4355572  6.019787 0.1119640
# Sandy Bar-Red River     -1.651559 -4.8792306  1.576113 0.5300245


# only sex model anova
# Perform the ANOVA
modelsex <- aov(pc1_scores18 ~ DM_data.frame$Sex, data = DM_data.frame)

# Print the summary of the model
summary(modelsex)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# DM_data.frame$Sex  2  425.9  212.93   6.788 0.00243 **
#   Residuals         51 1599.8   31.37                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsex)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18 ~ DM_data.frame$Sex, data = DM_data.frame)
# 
# $`DM_data.frame$Sex`
# diff        lwr         upr     p adj
# M-F  -4.507119 -9.0469123  0.03267414 0.0520520
# UN-F  4.828409 -0.4844396 10.14125731 0.0818612
# UN-M  9.335528  3.1643706 15.50668521 0.0017497

# pc2
# Perform the ANOVA
modelsexPC2 <- aov(pc1_scores18PC2 ~ DM_data.frame$Sex, data = DM_data.frame)

# Print the summary ofPC2 the model
summary(modelsexPC2)
# Df Sum Sq Mean Sq F value Pr(>F)  
# DM_data.frame$Sex  2  193.6   96.80   4.862 0.0117 *
#   Residuals         51 1015.5   19.91                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelsexPC2)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = pc1_scores18PC2 ~ DM_data.frame$Sex, data = DM_data.frame)
# 
# $`DM_data.frame$Sex`
# diff        lwr       upr     p adj
# M-F   3.728918   0.112074  7.345761 0.0419463
# UN-F -2.238554  -6.471289  1.994181 0.4147000
# UN-M -5.967471 -10.884019 -1.050923 0.0137846

#### Color by groups Site ####
pca18_PC1_PC2 <- fviz_pca_ind(DM_pca_FM,
                            geom.ind = "point", # show points only (but not "text"),
                            pointsize = 2.75,
                            pointshape = 21,
                            fill.ind = DM_data.frame$Site, # color by groups
                            col.ind = DM_data.frame$Site,
                            palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                            legend.title = "Site",
                            mean.point = FALSE,  # removes group mean point
                            title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2
png("pca18_PC1_PC2")

pca18_PC1_PC2_mean <- fviz_pca_ind(DM_pca_FM,
                                 geom.ind = "point", # show points only (but not "text"),
                                 pointsize = 2.75,
                                 pointshape = 21,
                                 fill.ind = DM_data.frame$Site, # color by groups
                                 col.ind = DM_data.frame$Site,
                                 palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                 legend.title = "Site",
                                 title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_mean
png("pca18_PC1_PC2_mean")

pca18_PC2_PC3 <- fviz_pca_ind(DM_pca_FM,
                            axes = c(2,3),
                            geom.ind = "point", # show points only (but not "text"),
                            pointsize = 2.75,
                            pointshape = 21,
                            fill.ind = DM_data.frame$Site, # color by groups
                            col.ind = DM_data.frame$Site,
                            palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                            legend.title = "Site",
                            mean.point = FALSE,  # removes group mean point
                            title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3
png("pca18_PC2_PC3")

pca18_PC2_PC3_mean <- fviz_pca_ind(DM_pca_FM,
                                 axes = c(2,3),
                                 geom.ind = "point", # show points only (but not "text"),
                                 pointshape = 21,
                                 pointsize = 2.75,
                                 fill.ind = DM_data.frame$Site, # color by groups
                                 col.ind = DM_data.frame$Site,
                                 palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                 legend.title = "Site",
                                 title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC2_PC3_mean


pca18_PC1_PC3_mean <- fviz_pca_ind(DM_pca_FM,
                                 axes = c(1,3),
                                 geom.ind = "point", # show points only (but not "text"),
                                 pointshape = 21,
                                 pointsize = 2.75,
                                 fill.ind = DM_data.frame$Site, # color by groups
                                 col.ind = DM_data.frame$Site,
                                 palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                 legend.title = "Site",
                                 title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC3_mean


pca18_PC3_PC4_mean <- fviz_pca_ind(DM_pca_FM,
                                 axes = c(3,4),
                                 geom.ind = "point", # show points only (but not "text"),
                                 pointshape = 21,
                                 pointsize = 2.75,
                                 fill.ind = DM_data.frame$Site, # color by groups
                                 col.ind = DM_data.frame$Site,
                                 palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                 legend.title = "Site",
                                 title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_mean


pca18_PC4_PC5_mean <- fviz_pca_ind(DM_pca_FM,
                                 axes = c(4,5),
                                 geom.ind = "point", # show points only (but not "text"),
                                 pointshape = 21,
                                 pointsize = 2.75,
                                 fill.ind = DM_data.frame$Site, # color by groups
                                 col.ind = DM_data.frame$Site,
                                 palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                 legend.title = "Site",
                                 title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_mean




#### biplot of individuals and variables ####
fviz_pca_biplot(DM_pca_FM, repel = TRUE,
                col.var = "#2E9FDF",  # variables colour
                col.ind = "#696969"  # individuals colour
)


#### colour individuals by site ####
fviz_pca_biplot(DM_pca_FM,
                geom.ind = "point",
                label = "var",
                col.ind = DM_data.frame$Site,  # colour by groups
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                axes = c(1,2),
                legend.title = "Site",
                col.var = "black",
                repel = TRUE,
                mean.point = FALSE  # removes group mean point
)

# colour individuals by Site
fviz_pca_biplot(DM_pca_FM, axes = c(1,2),
                geom.ind = "point",
                col.ind = DM_data.frame$Site,  # colour by groups
                mean.point = FALSE,  # removes group mean point
                pointsize = 2.5,
                palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                label = "var",
                alpha.var ="contrib",col.var = "black",
                repel = TRUE,
                select.var = list(contrib = 50),  # top contributing var
                legend.title = "Site"
)
pca18_ind_by_site <- fviz_pca_biplot(DM_pca_FM, axes = c(1,2),
                                     geom.ind = "point",
                                     col.ind = DM_data.frame$Site,  # colour by groups
                                     mean.point = FALSE,  # removes group mean point
                                     pointsize = 2.5,
                                     palette = c("#99d8c9", "#08519c","#E7B800", "#FC4E07"),
                                     label = "var",
                                     alpha.var ="contrib",col.var = "black",
                                     repel = TRUE,
                                     select.var = list(contrib = 50),  # top contributing var
                                     legend.title = "Site"
)
pca18_ind_by_site
png("pca18_ind_by_site")

# colour individuals by Site
fviz_pca_biplot(DM_pca_FM, axes = c(2,3),
                geom.ind = "point",
                col.ind = DM_data.frame$Site,  # colour by groups
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
fviz_pca_biplot(DM_pca_FM, axes = c(3,4),
                geom.ind = "point",
                col.ind = DM_data.frame$Site,  # colour by groups
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
fviz_pca_biplot(DM_pca_FM, axes = c(4,5),
                geom.ind = "point",
                col.ind = DM_data.frame$Site,  # colour by groups
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
pca18_PC1_PC2_sex <- fviz_pca_ind(DM_pca_FM,
                                geom.ind = "point", # show points only (but not "text"),
                                pointsize = 2.75,
                                pointshape = 21,
                                fill.ind = DM_data.frame$Sex, # color by groups
                                col.ind = DM_data.frame$Sex,
                                palette = c("#99d8c9", "#08519c","#E7B800"),
                                legend.title = "Sex",
                                mean.point = FALSE,  # removes group mean point
                                title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_sex

pca18_PC1_PC2_mean_sex <- fviz_pca_ind(DM_pca_FM,
                                     geom.ind = "point", # show points only (but not "text"),
                                     pointsize = 2.75,
                                     pointshape = 21,
                                     fill.ind = DM_data.frame$Sex, # color by groups
                                     col.ind = DM_data.frame$Sex,
                                     palette = c("#99d8c9", "#08519c","#E7B800"),
                                     legend.title = "Sex",
                                     title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC1_PC2_mean_sex
png("pca18_PC1_PC2_mean_sex")

pca18_PC3_PC2_mean_sex <- fviz_pca_ind(DM_pca_FM, axes = c(2,3),
                                     geom.ind = "point", # show points only (but not "text"),
                                     pointsize = 2.75,
                                     pointshape = 21,
                                     fill.ind = DM_data.frame$Sex, # color by groups
                                     col.ind = DM_data.frame$Sex,
                                     palette = c("#99d8c9", "#08519c","#E7B800"),
                                     legend.title = "Sex",
                                     title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC2_mean_sex
png("pca18_PC3_PC2_mean_sex")

pca18_PC3_PC4_mean_sex <- fviz_pca_ind(DM_pca_FM, axes = c(3,4),
                                     geom.ind = "point", # show points only (but not "text"),
                                     pointsize = 2.75,
                                     pointshape = 21,
                                     fill.ind = DM_data.frame$Sex, # color by groups
                                     col.ind = DM_data.frame$Sex,
                                     palette = c("#99d8c9", "#08519c","#E7B800"),
                                     legend.title = "Sex",
                                     title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC3_PC4_mean_sex
png("pca18_PC3_PC4_mean_sex")

pca18_PC4_PC5_mean_sex <- fviz_pca_ind(DM_pca_FM, axes = c(4,5),
                                     geom.ind = "point", # show points only (but not "text"),
                                     pointsize = 2.75,
                                     pointshape = 21,
                                     fill.ind = DM_data.frame$Sex, # color by groups
                                     col.ind = DM_data.frame$Sex,
                                     palette = c("#99d8c9", "#08519c","#E7B800"),
                                     legend.title = "Sex",
                                     title = "") +
  theme(aspect.ratio = 1/1) +
  theme(axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 14, face = "bold"))
pca18_PC4_PC5_mean_sex

png("pca18_PC4_PC5_mean_sex")



#### Some graphs ####
library(tidyverse)
theme_set(theme_light())
# All missing values can be filtered out by filtering the `sex` variable
# dat <- palmerpenguins::penguins %>% filter(!is.na(sex))
# by sex
point_plot_FL_W <- raw_data %>% 
  ggplot(aes(FL, W, fill = Sex)) +
  geom_jitter(size = 3, alpha = 0.5, shape = 21)
point_plot_FL_W

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
