library(readr)
library(tidyverse)

##### death data
elan_pat0 <- read.csv("/.../PAT.csv", 
                     header = TRUE, sep = "\t", dec = ",", na.strings = " ",
                     stringsAsFactors = TRUE)

elan_pat <- read.csv("/.../PAT.csv", 
                 header = TRUE, sep = "\t", dec = ",", na.strings = " ",
                 stringsAsFactors = TRUE)

summary(elan_pat)

# replace gender "O" with NA
elan_pat$dGeslacht <- as.factor(elan_pat$dGeslacht)
elan_pat["dGeslacht"][elan_pat["dGeslacht"] == "O"] <- NA

# remove rows with missing values
elan_pat <- na.omit(elan_pat)

# keep relevant columns
keeps <- c("PATNR","iGeboortejaar","iOverlijdensjaar","dGeslacht")
elan_pat = elan_pat[keeps]

# create standard columns for model input
elan_pat$patient_id <- elan_pat$PATNR
elan_pat$age <- elan_pat$iOverlijdensjaar - elan_pat$iGeboortejaar
elan_pat$female <- ifelse(elan_pat$dGeslacht == "M", 0, 1)
elan_pat$state <- "Death"
elan_pat$state_id <- 4
elan_pat$strategy_name <- "SOC"
elan_pat$strategy_id <- 1
elan_pat$time <- elan_pat$iOverlijdensjaar - 2017

# keep standard columns only
keeps <- c("state","strategy_name","female","age","patient_id","time","strategy_id","state_id")
elan_death = elan_pat[keeps]
summary(elan_death)

##### acute event data
elan_acute <- read.csv("/Users/r/Desktop/Thesis code/elan/COR.csv", 
                     header = TRUE, sep = "\t", dec = ",", na.strings = " ",
                     stringsAsFactors = TRUE)

summary(elan_acute)

# keep relevant columns
keeps <- c("PATNR","dDatum","dCorrespondentieICPC")
elan_acute = elan_acute[keeps]

# remove rows with missing values
elan_acute <- na.omit(elan_acute)

elan_acute <- subset(elan_acute, dCorrespondentieICPC == c("K74","K75","K76","K76.01",
                    "K76.02","K90","K90.01","K90.02","K90.03","K89","K77","K77.01",
                    "K77.02"))

# fill in patient info
elan_acute <- merge(elan_acute,elan_pat0, by = "PATNR")


# create standard columns for model input
elan_acute$patient_id <- elan_acute$PATNR

time0 <- "2017-01-01 00:00:00.0000000"
elan_acute$time <- difftime(elan_acute$dDatum, time0) / 365
elan_acute$time <- as.numeric(elan_acute$time)

elan_acute$year <- as.Date(elan_acute$dDatum, format="%Y")
elan_acute$year <- format(elan_acute$year, "%Y")
elan_acute$year <- as.numeric(elan_acute$year)
elan_acute$age <- elan_acute$year - elan_acute$iGeboortejaar


elan_acute$female <- ifelse(elan_acute$dGeslacht == "M", 0, 1)
elan_acute$state <- "Acute"
elan_acute$state_id <- 2
elan_acute$strategy_name <- "SOC"
elan_acute$strategy_id <- 1

# keep standard columns only
keeps <- c("state","strategy_name","female","age","patient_id","time","strategy_id","state_id")
elan_acute = elan_acute[keeps]
summary(elan_acute)

##### Hyper and stable data
elan_eps <- read.csv("/Users/r/Desktop/Thesis code/elan/EPS.csv", 
                     header = TRUE, sep = "\t", dec = ",", na.strings = " ",
                     stringsAsFactors = TRUE)

summary(elan_eps)

# keep relevant columns
keeps <- c("PATNR","dBegindatum","dICPC")
elan_eps = elan_eps[keeps]

elan_eps <- subset(elan_eps, dICPC == c("K74","K75","K76","K76.01",
                                        "K76.02","K90","K90.01","K90.02","K90.03","K89","K77","K77.01",
                                        "K77.02","K86","K87","K85"))

# remove rows with missing values
elan_eps <- na.omit(elan_eps)

#fill in patient info
elan_eps <- merge(elan_eps,elan_pat0, by = "PATNR")

# create standard columns for model input
elan_eps$patient_id <- elan_eps$PATNR

time0 <- "2017-01-01 00:00:00.0000000"
elan_eps$dBegindatum <- as.Date(elan_eps$dBegindatum)
elan_eps$time <- difftime(elan_eps$dBegindatum, time0) / 365
elan_eps$time <- as.numeric(elan_eps$time)

elan_eps$year <- as.Date(elan_eps$dBegindatum, format="%Y")
elan_eps$year <- format(elan_eps$year, "%Y")
elan_eps$year <- as.numeric(elan_eps$year)
elan_eps$age <- elan_eps$year - elan_eps$iGeboortejaar

elan_eps$female <- ifelse(elan_eps$dGeslacht == "M", 0, 1)
elan_eps$state <- ifelse(elan_eps$dICPC == c("K86","K87","K85"), 'Hyper','Stable')
elan_eps$state_id <- ifelse(elan_eps$state == "Hyper", 1, 3)
elan_eps$strategy_name <- "SOC"
elan_eps$strategy_id <- 1

# keep standard columns only
keeps <- c("state","strategy_name","female","age","patient_id","time","strategy_id","state_id")
elan_eps = elan_eps[keeps]

elan_eps <- elan_eps[elan_eps$age>20, ]
elan_eps <- elan_eps[elan_eps$time>0, ]
elan_eps <- na.omit(elan_eps)
summary(elan_eps)

# acute patients become stable after 1 year
elan_acute_after <- elan_acute
elan_acute_after$time <- elan_acute_after$time + 1
elan_acute_after$age <- elan_acute_after$age + 1
elan_acute_after$state_id <- 3
elan_acute_after$state <- "Stable"

##### SOC data integration
datasoc_e <- rbind(elan_death,elan_acute)
datasoc_e <- rbind(datasoc_e,elan_eps)
datasoc_e <- rbind(datasoc_e,elan_acute_after)
summary(datasoc_e)

# export prepared SOC data
write.csv(datasoc_e, "/.../datasoc_e.csv")

##### SOC+box data integration
databox <- read_csv("Desktop/Thesis code/model data/databox.csv", 
                    col_types = cols(...1 = col_skip()))
dataall_e = rbind(datasoc_e,databox)
dataall_e<- dataall_e[order(dataall_e$patient_id, dataall_e$time),]
dataall_e <- distinct(dataall_e)

summary(dataall_e)

# export prepared data for model
write.csv(dataall_e, "/.../dataall_e.csv")
