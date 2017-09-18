#### Ecoinformatics with pest/pesticide data from california citrus farms #####

# Objective: clean and aggregate raw data to use a GLMM to test if high pesticide use worsens pest infestations. we hypothesize that it might because pests can evolve resistance.  

#### Libraries ####
library(tidyverse)
library(data.table)
library(date)
library(Hmisc)
library(lme4)
library(r2glmm)

load("~/citrus_demo.RData")

#### 1. Loading Raw Data #### 
# Raw data: 114816 rows, 57 variables, including weekly sampling records of pests and pesticides from 300 citrus groves over a ten year period.  

pest <- read.csv("~/Documents/Desktop/R input files/citrus/Data_Dump_2.18.15.csv")[, 1:57]
pest <- as_tibble(pest)
# names(pest) 

# Additional data to determine if high pesticide use worsens pest problems

pesticides_full_ready <- read.csv("~/Documents/Desktop/R input files/citrus/pesticides_full_ready.csv")

pests_full <- read.csv("~/Documents/Desktop/R input files/citrus/pest_full_final.csv")

#### 2. Data Filtering #### 

# Cull the data only to include one especially damaging pest (mites) from the largest grower 

mites_raw <-
  pest %>% 
    select(grower_name, crop_year, replicate_num, field_num_ucd, sample_date, insect_type, density_index, infestation_sample, infestation_finding, infestation_percent) %>%
    filter(grower_name == "Sunwest" & insect_type == "CRM (Citrus Red Mite)") %>% 
    distinct() 

# Select pesticide spray data and filter result to grower of interest.  

pesticides <- 
  pest %>% 
  select(grower_name, crop_year, replicate_num, field_num_ucd, active_ingredient, Spray_Start_Date) %>%
  filter(grower_name == "Sunwest") %>% 
  distinct() 

#### 3. Data Cleaning #### 

# Use replace to go through and replace all of the count and other values with the appropriate infestation percents 

levels(mites_raw$density_index)
mites_raw$density_index <- recode_factor(mites_raw$density_index, "NONE" = '0', "+/-" = '0.026', "+/+"= '0.06', "++/++" = '0.104', "+++/+++" = '.1429', "sw Low" = '0.108', "sw Medium-Low" = '0.17', "sw Medium" = '0.224', "sw Medium-High" = '0.3925', "sw High" = '0.507')

#### 4. Handeling Missing Values #### 

# Find and fill missing infestation percent values 
mites<-
  mites_raw %>%
  mutate(infestation_percent = ifelse(is.na(infestation_percent), infestation_finding/infestation_sample, infestation_percent)) %>%
  # Move density_index values into infestation_percent 
  mutate(infestation_percent = ifelse(is.na(infestation_percent), as.numeric(as.character(density_index)), infestation_percent)) %>%   # clean pest data 
  select(-c(7:9))

names(mites)

#### 5. Associate pest and pesticide data to calculate impact of sprays ####

pesticides<-as.data.table(pesticides) 
pesticides$Spray_Start_Date<-as.Date(pesticides$Spray_Start_Date,"%m/%d/%y", tz = "PST") 
pesticides$Spray_Start_Date<-as.POSIXct(pesticides$Spray_Start_Date) 
# head(pesticides$Spray_Start_Date)

mites<-as.data.table(mites) 
mites$sample_date<-as.Date(mites$sample_date, "%m/%d/%y", tz = "PST")
mites$sample_date<-as.POSIXct(mites$sample_date) 

pesticides[, join_time:=Spray_Start_Date]
mites[, join_time:=sample_date]

setkey(pesticides, replicate_num, join_time)
setkey(mites, replicate_num, join_time)

# Set rolling window to one month in legnth 
one_month <- 60*60*24*30

## Backward rolling join ####

# Run a backward rolling join to associate pest data with preceding pesticide sprays (closest spray before a pest record). 

output_before<-pesticides[mites, roll = -one_month]

sprays_before <-
  output_before %>%
  drop_na(infestation_percent, Spray_Start_Date) %>%
  select(-starts_with("i."), -join_time) %>% 
  dplyr::rename(infestation_percent_before = infestation_percent, sample_date_before = sample_date)


#### Forward rolling join ####

# Run a forward rolling join to associate pest data with proceding sprays (closest spray after a pest record)

output_after<-pesticides[mites, roll = one_month]

sprays_after <-
  output_after %>%
  drop_na(infestation_percent, Spray_Start_Date) %>%
  select(-starts_with("i."), -join_time) %>% 
  dplyr::rename(infestation_percent_after = infestation_percent, sample_date_after = sample_date)

# Inner join before and after sprays 
before_after_sprays<-inner_join(sprays_before, sprays_after)

# Take the difference in pest population before and after a spray to determine efficacy 

before_after_clean <-
  before_after_sprays %>%
  filter(sample_date_before != sample_date_after) %>%
  mutate(efficacy = infestation_percent_after-infestation_percent_before) 
  

#### 6. Plot variance in pest infestation level by grove ####

ggplot(before_after_clean, aes(x=reorder(field_num_ucd, efficacy), y=efficacy)) + theme(axis.text.x=element_blank()) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="red", width=0.2) + stat_summary(fun.y=mean, geom="point", color="red") + geom_hline(aes(yintercept=0), color="green", linetype="dashed") + xlab("Grove") + ylab("Impact of spray on pest")

#### 7. Example analyses with generlaized linear mixed model ####

# Summarize data by grove and crop year
pest_summary <- 
  pests_full  %>% 
  group_by(grower_type, field_num_ucd,crop_year,species) %>% 
  select(grower_type, field_num_ucd,crop_year,species,infestation_percent_comb) %>% 
  summarise(
    mean_infes = mean(infestation_percent_comb, na.rm = TRUE),
    count = n()) 

pesticide_summary <- 
  pesticides_full_ready %>% 
  group_by(grower_type, field_num_ucd,crop_year) %>% 
  select(grower_type, field_num_ucd,crop_year,value) %>% 
  summarise(
    ingredient_count = n(),
    pesticide_impact = sum(value)) %>%
  group_by(grower_type,field_num_ucd)
 
pests_pesticides_sum <- inner_join(pest_summary, pesticide_summary)

# Filter data to focus on one particularly damaging pest 

target_pest<-filter(pests_pesticides_sum, species == "CRM (Citrus Red Mite)")

#### Visualize correlation between pesticides use and pest ####

ggplot(target_pest, aes(pesticide_impact, mean_infes)) + geom_point() + ylab("mean pest infestation proportion")  + xlab("pesticide impact score") + scale_y_continuous(limits = c(0, 0.5)) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + geom_smooth(method = 'loess')

#### Fit a GLMM using a binomial probability dist ####  
# For the response variable (proportions ranging from 0 to 1) and a random effect for grove (different growers manage their groves in diff. ways)

fit_me <- glmer(mean_infes~pesticide_impact + (1|field_num_ucd), binomial(link = logit), data=target_pest)
summary(fit_me)

# 1| specifies a random intercept for each grove

# Results are marginally significant. there is a tendency for heavy spraying to make pests worse, but it is not particularly strong. 

# save.image("~/citrus_demo.RData")