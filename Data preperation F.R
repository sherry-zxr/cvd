
# Import data

library(readr)
library(tidyverse)
data <- read.csv("/.../framingham_raw.csv", 
                   header = TRUE, sep = "\t", dec = ",", na.strings = " ",
                 stringsAsFactors = TRUE)

# Overview of the data
summary(data)

# Transform virables (count less than 4) as factor

col_names <- sapply(data, function(col)
                    length(unique(col)) < 4)
data[,col_names] <- lapply(data[col_names], factor)

# store the original data
data_original <- data

# make sure variable consistency
which(data$SEX == 'Male')
data$SEX[c(2, 18,  221,  991,  992,  993,  996, 1000)] <- 1

which(data$SEX == 'Man')
data$SEX[c(12,215)] <- 1

which(data$SEX == 'Female')
data$SEX[c(6281,  6284,  6288,  6290,  6736,  6742,  6748,  8798,  8801,
           8809,  8810,  8811,  8814,
           8815, 11427, 11544, 11545, 11548, 11549, 11618, 11623, 11626, 11627)] <- 2

data$SEX <- as.numeric(data$SEX)

which(data$SEX > 2)
data$SEX[c(1,     3,     4,  1002,  6744,  6750,  6752,  8803,
           8804,  8805,  8806, 11419, 11423)] <- NA
summary(data$SEX)

# find possible mistakes in data
summary(data)

# Possible mistake: max age is 600, max cigpday is 80, max bmi is 241.90, max HDLC is 189.

# find the rows with mistake and replace them to NA
which(data$AGE>120)
data$AGE[c(10,15)] <- NA
hist(data$AGE)

which(data$CIGPDAY>70)
data$CIGPDAY[c(89,637,769,1687,2121,3284,3506,4422)] <- NA
hist(data$CIGPDAY)

which(data$BMI>80)
data$BMI[11627] <- NA
hist(data$BMI)

which(data$HDLC>160)
data$HDLC[7714] <- NA
hist(data$HDLC)

# keep relevant columns
keeps <- c("RANDID","SEX","TIME","AGE","PREVAP","PREVCHD","PREVMI","PREVSTRK","PREVHYP",
           "ANYCHD","STROKE","CVD","DEATH","ANGINA","HOSPMI","MI_FCHD")
data = data[keeps]

# create standard columns for model input
data$state <- ifelse(data$DEATH == 1, 'Death',
                ifelse(data$HOSPMI == 1 |
                         data$ANGINA == 1 |
                         data$MI_FCHD == 1, 'Acute',
                     ifelse(data$ANYCHD == 1 |
                              data$STROKE == 1 |
                              data$CVD == 1 |
                               data$PREVAP == 1 |
                               data$PREVCHD == 1 |
                               data$PREVMI == 1 |
                               data$PREVSTRK == 1, 'Stable',
                             ifelse(data$PREVHYP == 1,'Hyper','NA'))))

data$state_id <- ifelse(data$state == "Death", 4,
                        ifelse(data$state == "Acute", 2,
                               ifelse(data$state == "Stable", 3, 1)))

data$strategy_id <- 1
data$strategy_name <- "SOC"
data$female <- ifelse(data$SEX == 1, 0, 1)
data$age <- data$AGE
data$patient_id <- data$RANDID
data$time <- data$TIME / 365

# keep standard columns only
keeps <- c("state","strategy_name","female","age","patient_id","time","strategy_id","state_id")
data = data[keeps]

# remove rows with missing values
data <- na.omit(data)
data <- data[data$state != "NA", ]

# standardize variable type
data <- data %>% 
  mutate_at(c('female', 'age', 'patient_id', 'time', 'strategy_id', 'state_id'), as.numeric)

view(data)
summary(data)

# export prepared SOC data
write.csv(data, "/.../datasoc_f.csv")

# SOC+box data integration
databox <- read_csv("Desktop/Thesis code/model data/databox.csv", 
                    +     col_types = cols(...1 = col_skip()))
dataall_f = rbind(data,databox)
summary(dataall_f)

# export prepared data for model
write.csv(dataall_f, "/.../dataall_f.csv")

