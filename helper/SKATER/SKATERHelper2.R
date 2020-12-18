# SKATERHelper2.0
# neighbors are nearby not Voronoi


# findNeighbors <- function(checklists_no_dups, rad){
#   # checklists_no_dups <- dplyr::distinct(checklists_df, latitude, longitude, .keep_all = TRUE)
#   ngbhrs <- list()
#   for(i in 1:nrow(checklists_no_dups)){
#     
#     list_i <- c()
#     
#     y_max <- checklists_no_dups[i, "latitude"] + rad
#     y_min <- checklists_no_dups[i, "latitude"] - rad
#     x_max <- checklists_no_dups[i, "longitude"] + rad
#     x_min <- checklists_no_dups[i, "longitude"] - rad
#     for(j in 1:nrow(checklists_no_dups)){
#       if(j != i){
#         if((checklists_no_dups[j, "latitude"] <= y_max) && (checklists_no_dups[j, "latitude"] >= y_min) && (checklists_no_dups[j, "longitude"] <= x_max) && (checklists_no_dups[j, "longitude"] >= x_min)){
#           list_i <- append(list_i, j)
#         }
#       }
#     }
#     ngbhrs <- append(ngbhrs, list(list_i))
#   }
#   return(ngbhrs)
# }

define_edges <- function(edges, v_covs){
  edge_weights <- list()
  for(i in 1:length(edges)){
    list_i <- c()
    edge1 <- i
    for(j in 1:length(edges[i]$nbs)){
      edge2 <- edges[[i]][j]
      
      diff <- dist(rbind(v_covs[edge1,], v_covs[edge2,]))[1]
      list_i <- append(list_i, diff)
    }
    edge_weights <- append(edge_weights, list(wts=list_i))
  }
  return(edge_weights)
}


adjListToMat <- function(adjList, edge_costs){
  for(q in 1:length(adjList)){
    row <- numeric(length(adjList))
    for(n in 1:length(adjList[[q]])){
      v2 <- adjList[[q]][n]
      e_wt <- edge_costs[[q]][n]
      row[v2] <- e_wt
    }
    # row <- as.data.frame(row)
    if(q==1){
      adj_mat = row
    } else {
      adj_mat = rbind(adj_mat, row)
    }
  }
  return(adj_mat)
}

calcSKATER.MSE2 <- function(checklists_df, nbs, covObj, num_sites, rad=.001, enforce_false_positives=FALSE){
  checklists_no_dups <- dplyr::distinct(checklists_df, latitude, longitude, .keep_all = TRUE)
  env_covs_df <- subset(checklists_no_dups, select = covObj$siteCovs)
  # ngbrs <- findNeighbors(checklists_no_dups, rad)
  edge_costs <- define_edges(nbs, env_covs_df)
  
  adjMat <- adjListToMat(nbs, edge_costs)
  adjMatInf <- replace(adjMat, adjMat == 0, Inf)
  
  G <- igraph::graph_from_adjacency_matrix(adjMatInf, weighted = T, mode = "undirected")
  mst_G <- igraph::mst(G)
  G_DF <- igraph::as_data_frame(mst_G)
  mst_m <- as.matrix(G_DF[,1:2])
  
  sites <- spdep::skater(mst_m, env_covs_df, ncuts = num_sites)
  
  checklists_no_dups$site <- sites$groups
  checklists_df$site <- numeric(nrow(checklists_df))
  for(row in 1:nrow(checklists_df)){
    long <- checklists_df[row, ]$longitude
    lat <- checklists_df[row, ]$latitude
    site_obj <- checklists_no_dups[checklists_no_dups$longitude == long & checklists_no_dups$latitude == lat,]
    checklists_df[row,]$site <- as.numeric(site_obj$site)
  }
  
  SKATER_MSE <- calcOccMSE(checklists_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = TRUE, enforce_false_positives = enforce_false_positives)
  return(SKATER_MSE)
}
