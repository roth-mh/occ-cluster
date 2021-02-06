# CREATE A STATIC DATASET
library(gtools)

f_in <- "../../../ICB General/data generation/NORM_weta_data_baseline_merged.csv"
WETA_2017_more <- read.delim(f_in, header=TRUE, sep = ",")

df <- WETA_2017_more
df$observation_date <- as.character(df$observation_date)
df$formatted_date <- mdy(df$observation_date)
df$as_date <- as_date(df$formatted_date)

# df <- subset(df, select = c(to_keep))

# df$occupied <- rbinom(length(df$occupied_prob), 1, df$occupied_prob)
# df$species_observed_syn <- rbinom(length(df$species_observed_syn), 1, df$species_observed_syn)

df_2017 <- subset(df, df$as_date >= as_date("2017-01-01"))
df_2017 <- subset(df_2017, df_2017$as_date < as_date("2018-01-01"))

colnames(df_2017a)

df_2017a <- subset(df,select=seq(from=25, to=372))
df_2017a <- subset(df_2017a, select = -c(latitude_300, longitude_300, latitude_150, longitude_150, latitude_75, longitude_75,
                                         checklist_id_150, checklist_id_300, checklist_id_75, system.index_300, system.index_150))

isCorrelated <- TRUE
while(isCorrelated){
  siteCovs <- list()
  while(length(siteCovs) < 5){
    x <- sample(seq(from=1, to=ncol(df_2017a)), 1)
    if(!(colnames(df_2017a)[x] %in% siteCovs)){
      siteCovs <- append(siteCovs, colnames(df_2017a)[x])
    }
  }
  cor_mat <- cor(df_2017a[,unlist(siteCovs)])
  if(length(cor_mat[abs(cor_mat) > .5]) == 5){
    isCorrelated <- FALSE
  } 
}


to_keep <- c(paste(siteCovs, sep = " "), paste(covObj$det_cov, sep=" "), 
               "observation_date", "formatted_date", "latitude", "longitude", "locality_id", 
               "observer_id", "checklist_id", "as_date")  



df <- subset(df, select = c(to_keep))
df$occupied <- -1 
df$species_observed_syn <- -1
df$occupied_prob <- -1

df_2017 <- subset(df, df$as_date >= as_date("2017-01-01"))
df_2017 <- subset(df_2017, df_2017$as_date < as_date("2018-01-01"))

write.csv(df, file = "../../../ICB General/data generation/UPDATED_COVS_df.csv")
write.csv(df_2017, file = "../../../ICB General/data generation/2017_UPDATED_COVS_df.csv")


