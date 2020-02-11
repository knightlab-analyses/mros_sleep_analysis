library(dplyr)
library(pspearman)
taxaMWWTest <- function(df, varName){
  df$Current_Log_Ratio <- lapply(df$Current_Log_Ratio,as.character)
  df$Current_Log_Ratio <- lapply(df$Current_Log_Ratio,as.numeric)
  df$Current_Log_Ratio <- unlist(df$Current_Log_Ratio)
  df <- na.omit(df)
  variable <- select(df, varName)
  variable.extract <- as.numeric(gsub(".*?([0-9]+).*", "\\1", variable[,1]))
  yesInd <- which(variable.extract == 1)
  if(length(yesInd) <= 1 || length(yesInd) == nrow(df)) return(99)
  yes <- df[yesInd, ]
  no <- df[-yesInd, ]
  p <- wilcox.test(yes$Current_Log_Ratio, no$Current_Log_Ratio)
  return(p$p.value)
}

#exclude the ones that have sample sizes less than 25
taxaSpearman <- function(df, varName){
  df$Current_Log_Ratio <- lapply(df$Current_Log_Ratio,as.character)
  df$Current_Log_Ratio <- lapply(df$Current_Log_Ratio,as.numeric)
  df$Current_Log_Ratio <- unlist(df$Current_Log_Ratio)
  df <- na.omit(df)
  if(nrow(df) <= 25) return(99)
  #value <- select(df, varName)
  p <- spearman.test(df$Current_Log_Ratio, df[,varName])
  return(p$p.value)
}

#NEED to pass processed data(change to numeric values)
get.pvalue<- function(df, folder.name, file.name, var.test){
  taxa.path <- paste("~/Dropbox/lab/sleep/Results/taxa/", folder.name, ".csv", sep = "")
  taxa <- read.csv(taxa.path)
  taxa.name <- as.character(taxa$Taxon)
  test.len <- length(var.test)
  results <- setNames(data.frame(matrix(ncol = test.len, nrow = 0)), var.test)
  file.path <- paste("~/Dropbox/lab/sleep/Results/taxa/", folder.name, sep = "")
  setwd(file.path)
  for(i in 1:10){
    dat = read.table(file.name[i], header=TRUE, sep='\t')
    name <- intersect(colnames(df), colnames(dat))
    df.delete <- df[ , -which(names(df) %in% c(name))]
    d <- merge(dat, df.delete, by.x = "Sample.ID", by.y = "X.SampleID")
    d <- d[, c("Current_Log_Ratio", var.test)]
    p.val <- numeric(length(var.test))
    for (j in 1:length(var.test)){
      if(var.test[j] %in% c("AMFVT_C1", "AMAMPT_C1", "PQBADSLP", "SLEEPHRS")){
        p.val[j] <- taxaMWWTest(d, var.test[j])
      } else {
        p.val[j] <- taxaSpearman(d, var.test[j])
      }
    }
    results[nrow(results)+1,] <- p.val
  }
  results$Taxaname <- taxa.name
  return(results)
}
