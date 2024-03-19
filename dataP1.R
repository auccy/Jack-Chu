library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caTools)
setwd("C:/Users/11373/Desktop/code")
d0 <- read.csv("C:/Users/11373/Desktop/code/mmse_medical.csv") 
d1 <- d0 %>%
  group_by(AIBL.Id) %>%
  mutate(Age_diff = Age - first(Age))
BLdataset <- d1
d2.0 <- d1[d1$Neuropsych.Simple.Classification != "AD",]
d2 <-  d1 %>%
  group_by(AIBL.Id ) %>% 
  subset(d1$Neuropsych.Simple.Classification == "AD") %>%
  distinct( AIBL.Id, .keep_all = TRUE)




d3 <- rbind(d1[d1$Neuropsych.Simple.Classification != "AD",], d2) %>% 
group_by(AIBL.Id)

d <- d3[order(d3$AIBL.Id),]


d5 <- d0 %>%
  group_by(AIBL.Id) %>%
  mutate(Age_diff = AGE-last(AGE))

# Group by participant_id and filter only the last record for each participant
last_records <- d0 %>%
  group_by(AIBL.Id) %>%
  slice_tail(n = 1)

mean(last_records$Age,na.rm = TRUE)
sd(last_records$Age,na.rm = TRUE)
max(last_records$Age_diff)
min(last_records$Age_diff)


state_counts <- table(last_records$DX)
print(state_counts)








