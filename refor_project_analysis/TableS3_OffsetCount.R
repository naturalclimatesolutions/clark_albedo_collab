
# install libs 
library(dplyr)
library(raster)

# setwd
setwd('...')

# read tabs
files <- list.files(pattern = '.csv')

albOffset <- read.csv('finalcount_albedo.csv') 
combOpp <- read.csv('finalcount_CompbinedOpportunityNCI.csv') 
combined_nc <- read.csv('finalcount_NetClimateImpact.csv')

# Count values for combined opp 
summary_table_verif <- combOpp %>%
  summarise(
    negative_lt0 = sum(Value < 0),
    positive_gt0 = sum(Value > 0),
    neutral_eq0 = sum(Value == 0), 
    na = sum(Value == "-nan(ind)")
  )
# View(summary_table_verif)

# Count values for Net Climate Imp
summary_table_verif <- albOffset %>%
  summarise(
    negative_lt0 = sum(Value < 0),
    positive_gt0 = sum(Value > 0),
    neutral_eq0 = sum(Value == 0), 
    na = sum(Value == "-nan(ind)")
  )
# View(summary_table_verif)

# Count values for Net Climate Imp
summary_table_verif <- combined_nc %>%
  summarise(
    negative_lt0 = sum(Value < 0),
    positive_gt0 = sum(Value > 0),
    neutral_eq0 = sum(Value == 0), 
    na = sum(Value == "-nan(ind)")
  )
# View(summary_table_verif)

###########################
###### Count Bin Vals #####
###########################

# Table S3: Proportion of maximum carbon storage offset by albedo in project-pixels
# Summary table with bins
breaks <-  c(-10000, 
             0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 
             10001)

## Overall Pixel Count
combined <- albOffset

combined$Value <- as.numeric(combined$Value)

summary_table <- combined %>%
  mutate(Range = cut(Value, breaks = breaks, labels = FALSE, right = FALSE),
         RangeLabel = paste(breaks[Range], breaks[Range + 1], sep = " to ")) %>%
  group_by(RangeLabel) %>%
  summarise(Count = n()) %>%
  ungroup()

print(summary_table)

## Combined Opportunity
combined <- combOpp 

combined[combined == "-nan(ind)"] <- NA
combined$Value <- as.numeric(combined$Value)

summary_table <- combined %>%
  mutate(Range = cut(Value, breaks = breaks, labels = FALSE, right = FALSE),
         RangeLabel = paste(breaks[Range], breaks[Range + 1], sep = " to ")) %>%
  group_by(RangeLabel) %>%
  summarise(Count = n()) %>%
  ungroup()

print(summary_table)

## Net Climate Impact
combined <- combined_nc

combined[combined == "-nan(ind)"] <- NA
combined$Value <- as.numeric(combined$Value)

summary_table <- combined %>%
  mutate(Range = cut(Value, breaks = breaks, labels = FALSE, right = FALSE),
         RangeLabel = paste(breaks[Range], breaks[Range + 1], sep = " to ")) %>%
  group_by(RangeLabel) %>%
  summarise(Count = n()) %>%
  ungroup()

print(summary_table)



