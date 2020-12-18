######
# DBSC
# DOI: 10.1016/j.cageo.2011.12.017
######
setwd("/Users/MarkRoth/Documents/Oregon State/Year 1/Research/")
source("eBird/site def/helper/helpers.R")
source("eBird/site def/helper/DBSCHelper.R")


# load("eBird/site def/data from SLURM/clustSite1000.Rdata")
# 
# 
# i <- 0
# for(item in li_1000){
#   i <- i + 1
#   row <- as.data.frame(item)
#   if(i == 1){
#     out_df = row
#   } else {
#     out_df = rbind(out_df, row)
#   }
# }
# 
# write.csv(out_df, file = "eBird/site def/helper/DBSCresults.csv")


# TRUE_OCC_COEFF <- c(-.5, 0.9772027357609328, -0.028406966644699994, 0.946381131416163, 0.7298163568647422, 0.20654714056879497)
# TRUE_DET_COEFF <- c(-1, 0.3935988905626474, 0.6282290623116853, -0.21418691151584734, 1.0118138402429864, 0.8490091884815945)

# f_in_WETA <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_2017.csv"
# WETA_2017 <- read.delim(f_in_WETA, header=TRUE, sep = ",")
# # WETA_2017 <- groomDataframe(WETA_2017, covObj$det_cov, covObj$occ_cov, syn_spec = T)
# f_in_syn_spec_form <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
# covObj <- loadCovariates(f_in_syn_spec_form)
# WETA_2017_region <- subset(WETA_2017, WETA_2017$latitude <= 44.5)
# WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
# WETA_2017 <- WETA_2017_region

v_pts <- formatVert(WETA_2017, covObj)

DT = RTriangle::triangulate(RTriangle::pslg(v_pts[c("latitude", "longitude")]))

globalMeanDT <- calcGlobalMeanSD(DT, v_pts)
globalMean <- globalMeanDT$mean
globalSD <- globalMeanDT$sd
edgeWts <- as.data.frame(globalMeanDT$edgeWts)

DT_graph <- igraph::graph_from_data_frame(edgeWts, directed = F)
DT_subgraphs <- removeGlobalLongEdges(DT, globalMean, globalSD, edgeWts, DT_graph)
de <- decompose(DT_subgraphs$graph)

new_graphs <- removeLongLocalEdges(de)
de_new_graphs <- list()
g_num <- 1
for(n_g in new_graphs){
  decomp <- decompose(n_g)
  for(d in decomp){
    de_new_graphs[[g_num]] <- d
    g_num <- g_num + 1
  }
}

T1 <- calcT1(de_new_graphs)

# step 3; make clusters & assign noise to separate cluster
clust_li <- createClusters(de_new_graphs, v_pts, T1)
sorted_DI <- as.data.frame(clust_li[1])
max_clust <- max(sorted_DI$clust) + 1
for(i in 1:nrow(sorted_DI)){
  pt <- sorted_DI[i,]
  if(pt$clust == -1){
    sorted_DI[i,]$clust <- max_clust
    max_clust <- max_clust + 1
  }
}


# assign clusters to checklists with the same lat/long
WETA_2017$clust <- -1
for(j in 1:nrow(sorted_DI)){
  c_id <- as.character(sorted_DI[j,]$v_name)
  lat_long <- v_pts[v_pts$checklist_id == c_id,][c("latitude", "longitude")]
  WETA_2017[WETA_2017$latitude == lat_long$latitude & WETA_2017$longitude == lat_long$longitude,]$clust <- sorted_DI[j,]$clust
}

WETA_2017$site <- WETA_2017$clust
table(WETA_2017$site)
DBSC.MSE <- calcOccMSE(WETA_2017, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)

# TODO: 
#   1. assign clusters to checklists at the same locations (like SKATER)
#   2. test occupancy? (given kmSq shit ... what does recover these sites really mean?)



















