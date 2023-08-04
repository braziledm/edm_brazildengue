# Installing the rEDM package
install.packages("rEDM")
library(rEDM)

#### Load complete original time series ####
original_series <- read.csv("CCM/Data/CompleteData_SP.csv")
original_series <- na.omit(original_series)

#Normalise data
norm_data <- original_series[,-c(1,2,6,7)]
head(norm_data)
for (i in 1:ncol(norm_data)){
  norm_data[,i]<-(norm_data[,i]-mean(norm_data[,i]))/sd(norm_data[,i])
}
head(norm_data)

# Get time series for each variable of interest
tot_rain <- norm_data[,1]
max_temp <- norm_data[,2]
min_temp <- norm_data[,3]

# Create Surrogate Series
tot_rain_surrogate <- make_surrogate_data(ts = tot_rain, 
                                          method = 'seasonal', 
                                          num_surr = 100, 
                                          T_period = 52, 
                                          alpha = 0.05)
max_temp_surrogate <- make_surrogate_data(max_temp, 'seasonal', 100, 52, 0.05)
min_temp_surrogate <- make_surrogate_data(min_temp, 'seasonal', 100, 52, 0.05)

# Export values to CSV
DataFile = "CCM/Data/.csv"
write.csv(tot_rain_surrogate,paste('Rain_Surrogates',basename(DataFile),sep=''),row.names = F)
write.csv(max_temp_surrogate,paste('MaxTemp_Surrogates',basename(DataFile),sep=''),row.names = F)
write.csv(min_temp_surrogate,paste('MinTemp_Surrogates',basename(DataFile),sep=''),row.names = F)

