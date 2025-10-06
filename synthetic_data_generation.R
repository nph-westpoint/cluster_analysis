library(tidyverse)
library(janitor)

cgm = read_csv('cgm_download_mg.csv') %>% clean_names()
clinical = read_csv('clinical.csv') %>% clean_names()
real_world = read_csv('real_world_worry_waves_dataset_with_extras.csv') %>% 
  clean_names() %>% 
  mutate(cluster=as.factor(cluster_numbers))

names(real_world)
intersect(names(cgm),names(clinical))

cgm_clinical <-
clinical %>% 
  left_join(cgm, by = c('user_id','diagnosis'))

names(cgm_clinical)
# Use the means of the Krakoff paper for  the total 180-minute glucose AUC model
# 180Glucose_AUC = -.32M(log)-.25AIR(log)+.04(Age)+.23(Female)-.09(Body Fat%)
# Weighted means from Krakoff across NGR and IGR categories.
m_avg <- 2.7*.683+2.1*.317
m_sd <- (3.6-2.2)/1.349*.683+(2.6-1.9)/1.349*.317
m_avg_log <- log(m_avg)

air_sd <- (308.5-125.4)/1.349*.683+(247-93.9)/1.349*.317
air_avg <- 200.5*.683+163.5*.317

air_avg_log <- log(air_avg)

body_fat_avg <- 30.1*.683+35.4*.317
body_fat_sd <- 9.3*.683+7.9*.317

age_avg <- 27*.683+29.4*.317
age_sd <- 6.7*.683+6.9*.317
female_avg <- .35*.683+.54*.317
female_sd <- sqrt(.35*(1-.35)/201)*.683+sqrt(.54*(1-.54)/145)*.317
auc_avg <- 17248*.683+19784*.317
auc_sd <- 1863.*.683+2166*.317




cgm_clinical %>% 
  lm(bgi_paper_hbgi~ age + auc_100,data=.) %>% 
  summary()

m_avg <- 2.7*.683+2.1*.317
m_sd <- (3.6-2.2)/1.349*.683+(2.6-1.9)/1.349*.317
m_log_sd <- sqrt(log(1 + m_sd^2 / m_avg^2))
m_avg_log <- log(m_avg)

air_sd <- (308.5-125.4)/1.349*.683+(247-93.9)/1.349*.317
air_avg <- 200.5*.683+163.5*.317
air_log_sd <- sqrt(log(1 + air_sd^2 / air_avg^2))

air_avg_log <- log(air_avg)

body_fat_avg <- 30.1*.683+35.4*.317
body_fat_sd <- 9.3*.683+7.9*.317

age_avg <- 27*.683+29.4*.317
age_sd <- 6.7*.683+6.9*.317
female_avg <- .35*.683+.54*.317
female_sd <- sqrt(.35*(1-.35))*.683+sqrt(.54*(1-.54))*.317
auc_avg <- 17248*.683+19784*.317
auc_sd <- 1863.*.683+2166*.317

tibble(feature = c('log(m)','log(air)','body_fat','age','sex','auc'),
       r_partial = c(-.38,-.35,-.01,.09,.19,NA),
       mean = c(log(m_avg),log(air_avg),body_fat_avg,age_avg,female_avg,auc_avg), 
       sd  = c(m_log_sd,air_log_sd,body_fat_sd,age_sd,female_sd,auc_sd)) %>% 
  mutate(beta= r_partial*(auc_sd/sd))


auc_intercept <- auc_avg-(-2202*m_avg_log-1110*air_avg_log+26.1 *age_avg+769*female_avg-2.21*body_fat_avg)

# Create model for CGM values. mage, coef_variation, bgi_paper_lbgi, bgi_paper_hbgi
mage_model <- cgm_clinical %>% 
  mutate(sex = if_else(rank(height, ties.method = "first") <= 32, "Female", "Male")) %>% 
  lm(mage~sex + age, data =.)
cv_model <- cgm_clinical %>% 
  mutate(sex = if_else(rank(height, ties.method = "first") <= 32, "Female", "Male")) %>% 
  lm(coef_variation~sex + age, data =.)
lbgi_model <- cgm_clinical %>% 
  mutate(sex = if_else(rank(height, ties.method = "first") <= 32, "Female", "Male")) %>% 
  lm(bgi_paper_lbgi~sex + age, data =.)
hbgi_model <- cgm_clinical %>% 
  mutate(sex = if_else(rank(height, ties.method = "first") <= 32, "Female", "Male")) %>% 
  lm(bgi_paper_hbgi~sex + age, data =.)
summary(hbgi_model)

# Use actual age and sex, then normal for AUC predictions in real_world
set.seed(456)
real_world_synthetic <-
real_world %>% 
  select(sex_wave1,age_wave1, cluster, clout,analytic,tone,authentic) %>% 
  filter(sex_wave1 %in% c("Female","Male")) %>% 
  rowwise() %>% 
  # mutate(auc_pred = auc_intercept -2202*log(abs(rnorm(1,m_avg,m_sd)))-1110*log(abs(rnorm(1,air_avg,air_sd)))+26.1*age_wave1+769*(sex_wave1=="Female")-2.21*abs(rnorm(1,body_fat_avg,body_fat_sd))) %>%
  mutate(auc_pred = auc_intercept -2202*log(m_avg)-1110*log(air_avg)+26.1*age_wave1+769*(sex_wave1=="Female")-2.21*body_fat_avg) %>%
  ungroup() %>% 
  # mutate(mage = predict(mage_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))) %>%
  # mutate(cv = predict(cv_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))) %>%
  # mutate(lbgi = predict(lbgi_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))) %>%
  # mutate(hbgi = predict(hbgi_model,newdata = rename(.,sex=sex_wave1,age = age_wave1)))
  mutate(mage = predict(mage_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))+rnorm(1150,mean=0,sd=summary(mage_model)$sigma/2)) %>%
  mutate(cv = predict(cv_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))+rnorm(1150,mean=0,sd=summary(cv_model)$sigma/2)) %>%
  mutate(lbgi = predict(lbgi_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))+rnorm(1150,mean=0,sd=summary(lbgi_model)$sigma/2)) %>%
  mutate(hbgi = predict(hbgi_model,newdata = rename(.,sex=sex_wave1,age = age_wave1))+rnorm(1150,mean=0,sd=summary(hbgi_model)$sigma/2))

write_csv(real_world_synthetic,'real_world_synthetic.csv')
names(real_world_synthetic)


library(randomForest)
library(iml)
# Load the dataset
data <- real_world_synthetic
# data <- read.csv("real_world_synthetic.csv")

# Select relevant columns
data_rf <- data %>% select(-age_wave1,-sex_wave1)
  
library(caret)
dummies <- dummyVars(~ cluster, data = data_rf)
dummy_df <- predict(dummies, newdata = data_rf)

# Remove rows with missing values
data_rf_dummy <- data_rf %>% 
  select(-cluster) %>% 
  bind_cols(dummy_df)
data_rf_dummy <- na.omit(data_rf_dummy)

str(data_rf$cluster)
set.seed(123)
rf_model <- randomForest(auc_pred ~ ., data = data_rf_dummy, importance = TRUE, mtry = 2)
print(rf_model)
rf_model$importance
varImpPlot(rf_model) 


# Create a Predictor object
X <- data_rf_dummy %>% select(-auc_pred)  # predictors
y <- data_rf_dummy$auc_pred
predictor <- Predictor$new(rf_model, data = X, y = y)

# Compute SHAP values
shap <- Shapley$new(predictor, x.interest = X[1, ])

# Plot SHAP values for the first observation
plot(shap)
