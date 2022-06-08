#####---- Temperature model for the Phd project of Emma Chollet #####
# 14.01.2022, JW
library(nlme)
library(lme4)
library(lmerTest) 
library(MuMIn)
library(tidyverse)
library(modelr)
library(jtools)
library(lattice)
library(reshape2)
# directory and file definitions
dir.data                <- "../data/"
dir.output              <- "../output/"
file.stations           <- "temperature_stations.dat"
temp.data.2005.2015    <- "../data/temp_data_files_2005_2015/"
temp.data.2016.2019    <- "../data/temp_data_files_2016_2019/"
dir.plots               <- "../output/temperature_model"

# load station info
data.stations           <- read.table(paste(dir.data,file.stations,sep=""),header=T,sep="\t", quote = "" )
# Dataframes to safe outputs in
data.stations.yearly    <- data.stations
data.stations.yearly    <- data.stations.yearly[c("ID", "Station", "station_elevation.m.", "catchment_area.km2.",
                                               "catchment_elevation.m.", "glaciation...", "coordinates")]

colnames(data.stations.yearly)  <- c("ID", "Station", "station_elevation", "catchment_area",
                                    "catchment_elevation", "glaciation", "coordinates")
rownames(data.stations.yearly)  <- data.stations.yearly$ID.

data.stations.yearly.output     <- data.frame() # empty dataframe to fill with yearly data


filenames.2005.2015     <- list.files(paste(temp.data.2005.2015), pattern="*.dat", full.names=TRUE)
filenames.2016.2019     <- list.files(paste(temp.data.2016.2019), pattern="*.csv", full.names=TRUE)

# read the files of all the stations
##### 2005 to 2015 data ####

for (file in filenames.2005.2015)
{
  # file <- filenames.2005.2015[62]
  dataset <- read.table(file, header=F, sep="\t")
  colnames(dataset)  <- c("stationID", "time", "temperature")
  
  # split the colum with time
  split.time <- strsplit(as.character(dataset$time), "-")                                   # split the time in 2
  temp <- data.frame(do.call(rbind, split.time))
  split.time <- strsplit(as.character(temp$X1), " ")                                        # split off the hour
  temp <- data.frame(do.call(rbind, split.time))
  split.time <- strsplit(as.character(temp$X1), "\\.")                                      # split the date in year, month and day
  dataset <- data.frame(dataset$stationID, do.call(rbind, split.time),temp$X1, temp$X2, dataset$temperature)
  colnames(dataset)  <- c("stationID", "year", "month", "day", "date", "hour", "temperature")
  
  # select summer
  dataset$month <- as.numeric(dataset$month)
  dataset$day <- as.numeric(dataset$day)
  summer <- which((dataset$month == 6 & dataset$day >= 21) | (dataset$month == 7) | (dataset$month == 8) | (dataset$month == 9 & dataset$day <21))
  dataset.summer <- dataset[summer,]  
  
  # mean yearly temperature
  yearly.mean <- mean(dataset$temperature)
  
  # mean of daily max temperature in summer
  daily.max <- aggregate(as.numeric(as.character(dataset.summer$temperature)),by=list(date=dataset.summer$date), FUN=max)
  mean.max.sum.temp <- mean(daily.max$x)
  
  # max morning temperature (= min daily temperature) in summer
  daily.min <- aggregate(as.numeric(as.character(dataset.summer$temperature)),by=list(date=dataset.summer$date), FUN=min)
  split.time <- strsplit(as.character(daily.min$date), "\\.")
  daily.min <- data.frame(do.call(rbind, split.time), daily.min$date, daily.min$x)
  colnames(daily.min)  <- c("year", "month", "day", "date", "temperature")
  max.morning.sum.temp.yearly <- aggregate(as.numeric(as.character(daily.min$temperature)), by=list(year=daily.min$year), FUN=max)
  max.morning.sum.temp<-  mean(max.morning.sum.temp.yearly$x) #REMOVE LATER
  
  ID <- as.character(unique(dataset$stationID))
  
  #this part is to produce the output dataframe by combining the station name, the year, and the corresponding temperature
  for (j in 1 : dim(max.morning.sum.temp.yearly)[[1]]){
    #j = 1
    max.morning.summer.temp <- max.morning.sum.temp.yearly[j,"x"]
    year  <- max.morning.sum.temp.yearly[j,"year"]
    tmp <- data.stations.yearly[data.stations.yearly$ID == ID,]
    tmp$year <- year
    tmp$max.morning.summer.temp <- max.morning.summer.temp
    data.stations.yearly.output <- rbind(data.stations.yearly.output, tmp)
  }
}
##### 2016 - 2019 data ####
for (file in filenames.2016.2019)
{
  # file <- filenames.2016.2019[1]
  
  # Extract time information from the data
  dataset.new <- read.csv(file, header = T, sep = ";", skip = 8)
  dataset.new$year <- substr(dataset.new$Zeitstempel, 1, 4)
  dataset.new$month <- substr(dataset.new$Zeitstempel, 6, 7)
  dataset.new$day <- substr(dataset.new$Zeitstempel, 9, 10)
  dataset.new$date <- substr(dataset.new$Zeitstempel, 1, 10)
  dataset.new$hour <- substr(dataset.new$Zeitstempel, 12, 13)
  
  # Select necessary columns
  dataset.new <- subset(dataset.new, select = c(Stationsnummer, year, 
                                                  month, day, date,
                                                  hour, Wert))
  
  # Change colnames to fit with previous projects
  colnames(dataset.new) <- c("stationID", "year", "month", "day", "date",
                              "hour", "temperature")
  
  # We remove the first row because the data we have from previous 
  # projects starts with the temperatour between the hours 00 - 01
  dataset.new <- dataset.new[-1,]
  
  # Convert Datatype
  # str(dataset.2017)
  dataset.new$month <- as.numeric(dataset.new$month)
  dataset.new$day <- as.numeric(dataset.new$day)
  
  # select summer
  summer <- which((dataset.new$month == 6 & dataset.new$day >= 21) | (dataset.new$month == 7) | (dataset.new$month == 8) | (dataset.new$month == 9 & dataset.new$day <21))
  dataset.new.summer <- dataset.new[summer,]
  
  # mean yearly temperature
  yearly.mean <- mean(dataset.new$temperature)
  
  # mean of daily max temperature in summer
  daily.max <- aggregate(as.numeric(as.character(dataset.new.summer$temperature)),by=list(date=dataset.new.summer$date), FUN=max)
  mean.max.sum.temp <- mean(daily.max$x)
  
  # max morning temperature (= min daily temperature) in summer
  daily.min <- aggregate(as.numeric(as.character(dataset.new.summer$temperature)),by=list(date=dataset.new.summer$date), FUN=min)
  split.time <- strsplit(as.character(daily.min$date), "-")
  daily.min <- data.frame(do.call(rbind, split.time), daily.min$date, daily.min$x)
  colnames(daily.min)  <- c("year", "month", "day", "date", "temperature")
  max.morning.sum.temp.yearly <- aggregate(as.numeric(as.character(daily.min$temperature)), by=list(year=daily.min$year), FUN=max)
  max.morning.sum.temp <-  mean(max.morning.sum.temp.yearly$x)
  
  ID <- as.character(unique(dataset.new$stationID))
  
  #this part is to produce the output dataframe by combining the station name, the year, and the corresponding temperature
  for (j in 1 : dim(max.morning.sum.temp.yearly)[[1]]){
    # j = 2
    max.morning.summer.temp <- max.morning.sum.temp.yearly[j,"x"]
    year  <- max.morning.sum.temp.yearly[j,"year"]
    tmp <- data.stations.yearly[data.stations.yearly$ID == ID,]
    tmp$year <- year
    tmp$max.morning.summer.temp <- max.morning.summer.temp
    data.stations.yearly.output <- rbind(data.stations.yearly.output, tmp)
  }

}

# read and prepare data ####
# only retain data for which we have temperatures
temp.rind <- which(!is.na(data.stations.yearly.output))
data.stations.yearly.output <- data.stations.yearly.output[temp.rind,]

# replace NA in glaciation by 0
data.stations$glaciation[is.na(data.stations$glaciation)] <- 0
data.stations.yearly.output$glaciation[is.na(data.stations.yearly.output$glaciation)] <- 0

# check for NAs, and remove those data
any(is.na(data.stations.yearly.output$station_elevation))
any(is.na(data.stations.yearly.output$catchment_area))
any(is.na(data.stations.yearly.output$catchment_elevation))   #2 stations with missing catchment_elevation
catchelev.rind <- which(!is.na(data.stations.yearly.output$catchment_elevation))
data.stations.yearly.output <- data.stations.yearly.output[catchelev.rind,]

# remove outlier with high glaciation and very low temperature
ind <- which(data.stations$mean.max.s.temp>5)
data.stations <- data.stations[ind,]

ind <- which(data.stations.yearly.output$"max.morning.summer.temp">5)
data.stations.yearly.output <- data.stations.yearly.output[ind,]
data.stations.yearly.output$year <- as.factor(data.stations.yearly.output$year)
data.stations.yearly.output$ID <- as.factor(data.stations.yearly.output$ID)

data.stations.yearly.output.2006.2015 <- data.stations.yearly.output
data.stations.yearly.output.2006.2015$year <- as.numeric(as.character(data.stations.yearly.output$year))
data.stations.yearly.output.2006.2015 <- subset(data.stations.yearly.output.2006.2015, year <= 2015)
#### Mixed effects model ####
lme1 <- lme(max.morning.summer.temp ~ log10(catchment_area) + catchment_elevation, random = ~1|year, data = data.stations.yearly.output)

lm1 <- lm(max.morning.summer.temp ~ log10(catchment_area) + catchment_elevation, data=data.stations.yearly.output)
lm2 <- lm(max.morning.summer.temp ~ log10(catchment_area) + catchment_elevation, data=data.stations.yearly.output.2006.2015)

#(mean.max.s.temp ~ log10(catchment_area) + catchment_elevation, data=data.stations

data.stations.yearly.output$lme1 <- predict(lme1)
data.stations.yearly.output$lm1 <- predict(lm1)
saveRDS(lme1, file = paste0(dir.output,"lme_temperature_model.rds"))
r.squaredGLMM(lme1)
r.squaredGLMM(lm1)

fixEffect <- fixef(lme1)
randEffect <- ranef(lme1)

data.stations.yearly.output.table <- as.data.frame(data.stations.yearly.output %>% 
                                                     group_by(year) %>% 
                                                     summarise(temp = mean(max.morning.summer.temp)))

max.morning.summer.temp.melt <- data.stations.yearly.output[c("ID", "Station", "year", "max.morning.summer.temp", "lme1", "lm1")]

max.morning.summer.temp.melt <- melt(max.morning.summer.temp.melt, id.vars = c("ID", "Station", "year"))
colnames(max.morning.summer.temp.melt)[5] <- "Temperature.value"
colnames(max.morning.summer.temp.melt)[4] <- "Temperature.method"

ggplot(data=max.morning.summer.temp.melt) + 
  geom_boxplot( aes(x=factor(year), y=Temperature.value, fill=factor(Temperature.method)), position=position_dodge(1)) +
  theme_minimal() + labs(x = "Year", y = "Temperature") + labs(fill='Temperature method') 


g1 <- ggplot(data.stations.yearly.output, aes(x=year, y=max.morning.summer.temp)) + 
  geom_boxplot()
g1  
g2 <- ggplot(data.stations.yearly.output, aes(x=year, y=lme1)) + 
  geom_boxplot()
g2  

g3 <- ggplot(data.stations.yearly.output, aes(x=year, y=lm1)) + 
  geom_boxplot()
g3 

plot(lme1)
plot(lme1, type = c("p", "smooth"))
plot(lme1, sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"))
qqnorm(resid(lme1))
qqline(resid(lme1))


plot(lm1)


g4 <- ggplot(data.stations.yearly.output, aes(x=max.morning.summer.temp, y=lme1, color = year)) + 
  geom_point() + qqline()
g4 

g4 <- ggplot(data.stations.yearly.output, aes(x=lme1, y=lm1)) + 
  geom_point()
g4 
