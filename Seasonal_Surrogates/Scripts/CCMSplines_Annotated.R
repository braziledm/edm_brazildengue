# This function obtains the convergent-cross mapping between seasonal
#               drivers (target variables) and one affected variable, 
#               e.g. dengue cases.
# Adapted from: https://github.com/mathbio/malariaCCM/blob/master/CCM/CCMSplines.R

# please see: vignette("rEDM-tutorial")

# input
# OriginalFile: File with data to be analysed
# SeasonalityFile: File with seasonality data
# DateCol: Column with date
# LibraryCol: Variable column to be cross-mapped from
# TimeForPred: Time for prediction
# MaxE: Maximum embedding to be used in optimal search
# NSurr: Number of surrogate series
# Normalize: Normalize variables (T or F)
# alpha: significance for statistical test
# output
# Boxplot_tp=[TimeForPred].eps: plot with ccm analysis for the original data 
#                               and surrogates (boxplot). Red asterisks
# Boxplot_tp=[TimeForPred].csv: Spreadsheet with ccm values for [NSurr] surrogate series and 
#                               for the original series (first line of spreadsheet).

# Installing the rEDM package
install.packages("rEDM")
library(rEDM)

CCMSplines<-function(OriginalFile= "CCM/Data/CompleteData_SP.csv",
                     SeasonalityFile="CCM/Data/Seasonality_ClimateVariables.csv",
                     DateCol=1,
                     LibraryCol=2,
                     TimeForPred=0,
                     MaxE=10,
                     NSurr=500,
                     Normalize=T,
                     alpha=0.05){
  # Read in all data
  Original<-read.csv(OriginalFile)
  Seasonality<-read.csv(SeasonalityFile)
  AllVariablesNames<-names(Original)
  DriversNames<-AllVariablesNames[-c(DateCol,LibraryCol)]
  Drivers<-seq(1:ncol(Original))[-c(DateCol,LibraryCol)]
  Target<-Original[,LibraryCol]
  
  # Normalise the data
  if (Normalize==T){
    NormData<-Original[,-DateCol]
    for (i in 1:ncol(NormData)){
      NormData[,i]<-(NormData[,i]-mean(NormData[,i]))/sd(NormData[,i])
    }
    Original[,-DateCol]<-NormData
  }
  
  # For storing residuals, rho values for surrogates, and significance values
  Residuals<-as.data.frame(matrix(NA,NROW(Original),length(Drivers)))
  RhoSurr<-matrix(NA,(NSurr+1),length(Drivers))
  signf<-matrix(NA,1,length(Drivers))
  
  # Creating 500 surrogate time series
  count<-1
  for (j in Drivers){
    Residuals[,count]<-(Original[,j]-Seasonality[,j]) #remove seasonality
    block_temp<-as.data.frame(matrix(NA,NROW(Original),2))
    E_start<-(matrix(NA,1,NSurr))
    lib_ccm<-c(1,NROW(Original)-1) # library size
    for (i in 1:NSurr){
      block_temp[,1]<-Target #Cases
      block_temp[,2]<-sample(Residuals[,1])+Seasonality[,3] #series with reshuffled residuals
      out.temp <- do.call(rbind,lapply(1:MaxE, function(E_i){
        ccm(block=block_temp,E=E_i,lib=lib_ccm,pred=lib_ccm,
            lib_sizes = NROW(block_temp),exclusion_radius=0,
            random_libs = FALSE,num_sample=1,tp = 0,
            lib_column = 1,target_column = 2)
      }))
      E_start[i] <- out.temp$E[which.max(out.temp$rho[2:MaxE])+1]
      df.out.ccm <- ccm(block=block_temp,E=E_star[i],lib=lib_ccm,
                        pred = pred_ccm,lib_sizes = NROW(block_temp),
                        exclusion_radius=0,random_libs = FALSE,
                        num_sample=1,tp = TimeForPred)
      RhoSurr[(i+1),count]<-df.out.ccm$rho
    } # end loop surrogates
    
    ## Original
    block_temp[,2]<-Original[,j]
    # optimal E for original
    out.temp <- do.call(rbind,lapply(1:MaxE, function(E_i){
      # pred_ccm <- make_pred_nozero(Target,E_i)
      ccm(block=block_temp,E=E_i,lib=lib_ccm,pred=pred_ccm,
          lib_sizes = NROW(block_temp),exclusion_radius=0,
          random_libs = FALSE,num_sample=1,tp = TimeForPred,
          lib_column = 1,target_column = 2)
    }))
    E_startOriginal<- out.temp$E[which.max(out.temp$rho[2:MaxE])+1]
    # pred_ccm <- make_pred_nozero(Target,E_starOriginal)
    df.out.ccm <- ccm(block=block_temp,E=E_startOriginal,
                      lib=lib_ccm,pred = pred_ccm,
                      lib_sizes = NROW(block_temp),exclusion_radius=0,
                      random_libs = FALSE,num_sample=1,tp = TimeForPred)
    RhoSurr[1,count]<-df.out.ccm$rho
    count<-count+1
  } #end loop drivers
  
  # check significance
  for (k in 1:length(Drivers)){
    signf[k]<-(sum(RhoSurr[1,k]>RhoSurr[,k])/NSurr)>(1-alpha)
  }
  
  RhoSurr<-as.data.frame(RhoSurr,row.names = F)
  colnames(RhoSurr)<-DriversNames
  write.csv(RhoSurr,paste('Boxplot_tp=',TimeForPred,'.csv',sep=''),row.names = F)
  
  boxplot(RhoSurr)
  title(main =paste('tp= ',TimeForPred,sep=''),ylab='rho')
  points(seq(1,length(Drivers)),RhoSurr[1,],pch=16,cex=1.25,lwd=1.25)
  points(seq(1,length(Drivers))[signf],RhoSurr[1,signf],pch=8,col='red',cex=2,lwd=2)
  
} # end function



# Run function
CCMSplines(OriginalFile="CCM/Data/CompleteData_SP.csv",
                     SeasonalityFile="CCM/Data/Seasonality_ClimateVariables.csv",
                     DateCol=1,
                     LibraryCol=2,
                     TimeForPred=0,
                     MaxE=10,
                     NSurr=500,
                     Normalize=T,
                     alpha=0.05)

