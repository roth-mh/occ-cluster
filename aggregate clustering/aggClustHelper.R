# AGG Clustering Helper Fcn

library(gtools)

##########
# checklist similarity
#  - returns 1 if two checklists are in the same site
#     in one of the clusterings but in different sites
#     in the other clustering
#  - o/w returns 0
##########
# clusterSimil <- function(site1, site2){
#   if(as.character(site1[[1]]) == as.character(site1[[2]]) & 
#      as.character(site2[[1]]) != as.character(site2[[2]])){
#     return(1)
#   }
#   if(as.character(site1[[1]]) != as.character(site1[[2]]) &
#      as.character(site2[[1]]) == as.character(site2[[2]])){
#     return(1)
#   }
#   return(0)
# }


checklistSimil <- function(site1, site2){
  if(as.character(site1[[1]]) == as.character(site1[[2]]) & 
     as.character(site2[[1]]) == as.character(site2[[2]])){
    return(0)
  }
  return(1)
}




##########
# cluster similarity
#  - counts the number of pairwise checklist disagreements
#     between two clusterings
#  - not scalable, O(n^2)



#  - calc the weights for 2 checklists given
#      a df with sites as columns
##############################
# clust.DF contains each object and 2 clusterings in
# the columns site1 and site2
##########
clusterSimil <- function(sites.DF, ignoreFirstCol){
  
  if(ignoreFirstCol){
    ignoreFirstCol <- 1
  } else {
    ignoreFirstCol <- 0
  }
  
  distMat <- matrix(nrow = nrow(sites.DF), ncol = nrow(sites.DF))
  
  chcklst_comb <- combinations(nrow(sites.DF), 2)
  # clust_comb <- combinations(ncol(sites.DF)-1, 2)
  for(i in 1:nrow(chcklst_comb)){
    # disp("on checklist row: ", as.character(i))
    chklsts <- chcklst_comb[i,]
    
    # sites.DF[chklsts,]
    tot_distance <- 0
    
    for(j in (1+ignoreFirstCol):ncol(sites.DF)){
      if(sites.DF[chklsts[[1]],j] != sites.DF[chklsts[[2]],j]){
        tot_distance <- tot_distance + 1  
      }
      # distance <- checklistSimil(site1_lst, site2_lst)
    }
    
    # calc fraction of disagreements
    tot_distance <- tot_distance/(ncol(sites.DF) - ignoreFirstCol)
    distMat[chklsts[1], chklsts[2]] <- tot_distance
  }
  diag(distMat) = 0
  distMat[lower.tri(distMat)] = t(distMat)[lower.tri(distMat)]
  return(distMat)
}



case1 <- function(sites){
  if(length(unique(sites)) > 1){
    return(FALSE)
  } else if(unique(sites) == -1){
    return(TRUE)
  }
  return(FALSE)
}

case2 <- function(sites){
  if(length(unique(sites)) > 1){
    return(FALSE)
  } else if(unique(sites) != -1){
    return(TRUE)
  }
  return(FALSE)
}

#########
# function to chain 
# different sites to be the same
#########
chainSites <- function(over, proj_cent, site_num){
  # dealing with unlabeled checklists
  unlab_chks_over <- over[over$site == -1,]$checklist_id
  if(length(unlab_chks_over) > 0){
    proj_cent[proj_cent$checklist_id %in% unlab_chks_over,]$site <- site_num
  }
  
  buffer_sites <- c(site_num, over[over$site != -1,]$site)
  
  proj_cent[proj_cent$site %in% unique(buffer_sites),]$site <- site_num
  return(proj_cent)
}




######
# BALLS CONSENSUS CLUSTERING ALGORITHM
######
# Step 1: sort the vertices in order of 
#   decreasing incident weight
# 
# Step 2: select the next unclustered vertex, v
# Step 3: find the set of vertices, B, that
#   are at most distance .5 from v
# Step 4: if d(u, B) <= alpha, B <- B U {v}
#   where d(u, B) is the avg distance of all
#   nodes in B to u
######


# ALPHA: threshold for mean distance from u; OR
#         the average fraction of disagreements
#         i.e., lower ==> tighter groups
#
#   paper proves 1/4 is a 2-approx, but 2/5 
#     works better in practice

# BETA: threshold for what is considered a neighbor
#       i.e., lower ==> looser groups
#
#     - this is not considered a parameter in the 
#       original paper but I don't see why it 
#       can't be       
ballsAggr <- function(dMat, v_list, ALPHA = .4, BETA = .5){
  v_wghts <- as.data.frame(matrix(nrow = nrow(dMat), ncol = 3))
  colnames(v_wghts) <- c("vertex", "total_weight", "site")
  
  v_wghts$vertex <- seq(1:nrow(dMat))
  v_wghts$pvs_vertex <- v_list
  v_wghts$total_weight <- rowSums(dMat)
  v_wghts$site <- -1
  v_sorted <- v_wghts[order(v_wghts$total_weight),]
  
  next_site <- 1
  for(v in 1:nrow(v_sorted)){
    if(v_sorted[v,]$site == -1){
      vert <- v_sorted[v,]$vertex
      # all 'neighbors' of ver
      B_with_v <- seq(1:nrow(v_sorted))[dMat[vert,] <= BETA]
      
      # removing the vertex itself
      B <- B_with_v[!(B_with_v %in% vert)]
      
      if(length(B) == 0){
        v_sorted[v,]$site <- next_site
      } else {
        d_v_B <- mean(dMat[vert,][B])
        if(d_v_B <= ALPHA){
          # disp("here")
          v_sorted[v_sorted$vertex %in% B_with_v,]$site <- next_site
        } else {
          v_sorted[v,]$site <- next_site
        }
      }
      next_site <- next_site + 1
    }
  }
  return(as.data.frame(v_sorted[c("pvs_vertex", "site")]))
  
}



###
#########
###############
#####################
#
# AGGLOMERATIVE ALGORITHM FUNCS
#
#####################
###############
#########
###
setClass("Cluster", slots = list(name="numeric", objects="numeric", closest="numeric", min_dist="numeric"))

agglCluster <- function(dMat){
  aggl_sites <- as.data.frame(matrix(data = c(seq(1:nrow(dMat)), seq(1:nrow(dMat))), nrow = nrow(dMat), ncol = 2))
  colnames(aggl_sites) <- c("vertex", "site")
  
  # initialize objects with each being in their own cluster
  CLUSTERS <- list()
  cluster_name <- 1
  for(v in aggl_sites$vertex){
    new_c <- new("Cluster", name=cluster_name, objects=c(v), closest=-1, min_dist=-1)
    CLUSTERS <- append(CLUSTERS, new_c)
    cluster_name <- cluster_name + 1
  }
  
  n_CLUSTERS <- list()
  cluster_name <- 1
  for(v in CLUSTERS){
    min_obj <- find_min_cluster(v, CLUSTERS, dMat) 
    pop_clst <- new("Cluster", name=cluster_name, objects=c(v@objects), closest=min_obj$clust, min_dist=min_obj$dist)
    n_CLUSTERS <- append(n_CLUSTERS, pop_clst)
    cluster_name <- cluster_name + 1
  }
  
  idx <- find_min_dist(n_CLUSTERS)
  while(n_CLUSTERS[[idx]]@min_dist < .5){
    new_CLUSTERS <- mergeClusters(n_CLUSTERS, dMat, n_CLUSTERS[[idx]]@name, n_CLUSTERS[[idx]]@closest, cluster_name)
    cluster_name <- cluster_name + 1
    stopifnot(length(new_CLUSTERS) < length(n_CLUSTERS))
    
    # recalculate the min distances ... untested
    new_CLUSTERS2 <- list()
    for(v in new_CLUSTERS){
      min_obj <- find_min_cluster(v, new_CLUSTERS, dMat) 
      v@closest <- min_obj$clust
      v@min_dist <- min_obj$dist
      new_CLUSTERS2 <- append(new_CLUSTERS2, v)
    }
    stopifnot(length(new_CLUSTERS) == length(new_CLUSTERS2))
    n_CLUSTERS <- new_CLUSTERS2
    idx <- find_min_dist(n_CLUSTERS)
    if(idx == -1){
      break
    }
  }
  
  return(n_CLUSTERS)
}






# Merges two groups in cluster_list w/ names specified
# creates new group with new_clust_name
#   - removes to-be-merged groups prior to merging
mergeClusters <- function(clusters, dMat, clust_name1, clust_name2, new_clust_name){
  
  clust_name1.idx <- find_clust_idx(clusters, clust_name1)
  clust_name2.idx <- find_clust_idx(clusters, clust_name2)
  new_clst <- new("Cluster", name=new_clust_name, objects=c(clusters[[clust_name1.idx]]@objects, clusters[[clust_name2.idx]]@objects), closest=-1, min_dist=-1)
  clusters <- removeClusters(clusters, clust_name1, clust_name2)
  min_obj <- find_min_cluster(new_clst, clusters, dMat)
  new_clst@closest <- min_obj$clust
  new_clst@min_dist <- min_obj$dist
  clusters <- append(clusters, new_clst)
  return(clusters)
}


# given a cluster name, finds the index in cluster list
find_clust_idx <- function(clusters, cls_name){
  for(i in 1:length(clusters)){
    if(clusters[[i]]@name == cls_name){
      return(i)
    }
  }
  disp("ERROR, clust name not found")
}

# removes given clusters from cluster list
removeClusters <- function(clusters, clust_name1, clust_name2){
  to_del <- c()
  for(i in 1:length(clusters)){
    if(clusters[[i]]@name %in% c(clust_name1, clust_name2)){
      to_del <- append(to_del, i)
    }
  }
  clusters <- clusters[-c(to_del)]
  return(clusters)
}

# finds the minimum distance in the current
# cluster list
find_min_dist <- function(clusters){
  m_dist <- Inf
  m_idx <- -1
  i <- 1
  for(clst in clusters){
    if(clst@min_dist < m_dist){
      m_idx <- i
      m_dist <- clst@min_dist
    }
    i <- i + 1
  }
  return(m_idx)
}


# finds the closest cluster from the
# given cluster
find_min_cluster <- function(c, clusters, dMat){
  min_dist <- Inf
  min_clust <- -1
  for(clst in clusters){
    if(c@name != clst@name){
      # mean(find all edges between clusters)
      mean_dist <- calcMeanDist(c, clst, dMat)
      if(mean_dist <= min_dist){
        min_dist <- mean_dist
        min_clust <- clst@name
      }
    }
  }
  return(list(clust=min_clust, dist=min_dist))
}


# assuming we have an undirected graph,
# therefore we only have to check the edges from
# cls1 to clst2

# calculates the average distance between two clusters
#   - i.e. the average of all edges between the two
#     clusters
calcMeanDist <- function(clst1, clst2, dMat){
  dist <- c()
  for(v in clst1@objects){
    dv <- dMat[v,clst2@objects]
    dist <- append(dist, dv)
  }
  return(mean(dist))
}


# # # # # # # # # # 
# this function prepares the WETA.sites df
# 
# first it groups checklists by territory size
# then uses the true sites to populate a syn species
# and loads those detections into the original dataset
# # # # # # # # # # 
prepWETA.sites <- function(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS){
  
  # Looking only at unique locations
  # bc if same location ==> same site (proof)
  WETA_uniq <- sqldf("SELECT * FROM WETA_2017 GROUP BY latitude, longitude")
  WETA_uniq$site <- -1
  
  ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/OR_2020_fittedImage.tif")
  prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ################
  # define the checklists as centroids of the buffer
  # reproject centroids into appropriate projection for meters
  # create buffer around each centroid of territory size
  ################
  centroids <- st_as_sf(WETA_uniq, coords = c("longitude", "latitude"), crs=prj)
  proj_cent <- st_transform(centroids, crs(ras.OR2))
  buffers <- st_buffer(proj_cent, dist = MDD)
  
  ########
  # for each buffer, find the checklists that overlap
  # if all are unlabeled, create new cluster and add them
  # if some are labeled (possibly some are unlabeled)
  # union them
  ########
  i <- nrow(WETA_uniq)
  for(b in 1:nrow(buffers)){
    buff <- buffers[b,]
    over <- proj_cent[st_intersection(buff, proj_cent),]
    
    ########
    # plot all overlapping point for a given buffer
    ########
    # ggplot() +
    #   geom_sf(data = buff, pch = 4, col = "gray40") +
    #   geom_sf(data = over, col = "red") +
    #   theme_minimal()
    ########
    
    if(length(over) > 1){
      sites <- proj_cent[proj_cent$checklist_id %in% over$checklist_id,]$site
      # over$site
      # if sites are both -1: site <- i
      if(case1(sites)){
        proj_cent[proj_cent$checklist_id %in% over$checklist_id,]$site <- i
        i <- i + 1
      }
      # all of the same site
      else if(case2(sites)){
        
      }
      # some mix of sites ==> union them
      else {
        proj_cent <- chainSites(over, proj_cent, i)
        i <- i + 1
      }
      
    } else {
      # what is this condition? single point
      proj_cent[proj_cent$checklist_id == over$checklist_id,]$site <- i
      i <- i + 1
    }
  }
  
  # rename sites
  j <- 1
  for(s in unique(proj_cent$site)){
    proj_cent[proj_cent$site == s,]$site <- j
    j <- j + 1
  }
  
  ############
  # HIST & TABLE
  ############
  # hist(proj_cent$site, breaks=seq(0,length(unique(proj_cent$site))), main = sprintf("MDD: %s, # Visits/Site, all OR", MDD) )
  # dt <- as.data.frame(table(table(proj_cent$site)))
  # colnames(dt) <- c("num visits", "freq")
  # grid.table(dt)
  ############
  
  # load ground truth data set
  f_name <- "../../../clusteredSites_2020-12-26_.csv"
  clusted_sites <- read.delim(f_name, sep=",")
  WETA_filtered <- WETA_2017
  WETA_filtered$site <- -1
  WETA_filtered$det_prob <- -1
  
  # link sites with their respective checklists
  for(row in 1:nrow(WETA_filtered)){
    # disp(row)
    long <- WETA_filtered[row, ]$longitude
    lat <- WETA_filtered[row, ]$latitude
    site_obj <- clusted_sites[clusted_sites$longitude == long & clusted_sites$latitude == lat,]
    WETA_filtered[row,]$site <- as.character(site_obj$site)
  }
  
  WETA_2017$vertex <- seq(1:nrow(WETA_2017))
  og_data <- WETA_2017
  
  truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
  og_data$species_observed_syn <- truth_df$species_observed_syn
  
  WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
  return(list(WETA_sites, og_data, proj_cent, truth_df))
}

# construct many sites and join them together
appendSites <- function(tests, WETA_sites, og_data){
  
  need_L_join <- FALSE
  
  df_to_join <- list()
  if(!is.na(tests$DBSC)){
    DBSC_df <- runDBSC(og_data, covObj)  
    df_to_join <- append(df_to_join, list(DBSC_df))
  }
  
  if(!is.na(tests$eBird)){
    eBird_df <- filter_repeat_visits(
      og_data,
      min_obs = 2,
      max_obs = 10,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id", "observer_id")
    )
    df_to_join <- append(df_to_join, list(eBird_df))
    need_L_join <- TRUE
  }
  
  if(!is.na(tests$eBird_simple)){
    eBird_simple_df <- filter_repeat_visits(
      og_data,
      min_obs = 1,
      max_obs = 1000000,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id")
    )
    df_to_join <- append(df_to_join, list(eBird_simple_df))
  }
  
  if(!is.na(tests$rounded)){
    for(i in tests$rounded){
      WETA_2017_i <- roundLatLong(og_data, i)
      eBird_rounded_df <- filter_repeat_visits(
        WETA_2017_i,
        min_obs = 1,
        max_obs = 1000000,
        annual_closure = TRUE,
        date_var = "formatted_date",
        site_vars = c("rounded_locality_id")
      )
      df_to_join <- append(df_to_join, list(WETA_2017_i))
    }
  }
  
  if(!is.na(tests$clustGeo)){
    for(i in tests$clustGeo){
      disp(i)
      clustGeo_df_i <- clustGeoSites(alpha = i[1], og_data, covObj, num_sites = i[2])
      df_to_join <- append(df_to_join, list(clustGeo_df_i))
    }
  }
  ######
  if(need_L_join){
    disp("implement nasty left join for eBird")
    # # max_site <- 1
    # for(row in 1:nrow(WETA_sites)){
    #   if(is.na(WETA_sites[row,]$site4)){
    #     WETA_sites[row,]$site4 <- max_site
    #     max_site <- max_site + 1
    #   }
    # }
    # WETA_sites <- left_join(WETA_sites, eBird_df[c("checklist_id", "site")], by="checklist_id", suffix=c("", "4"))
  }
  ######
  i <- 1
  for(df in df_to_join){
    # TODO: debug here; one doesn't have checklist_id
    WETA_sites <- inner_join(WETA_sites, df[c("checklist_id", "site")], by="checklist_id", suffix=c("", as.character(i)))  
    i <- i + 1
  }
  
  return(WETA_sites)
}


combineMethods <- function(proj_cent, WETA_sites, og_data){
  v_by_s_df <- data.frame(vertex = numeric(nrow(proj_cent)), site = numeric(nrow(proj_cent)), stringsAsFactors = FALSE)
  
  # run the balls aggregation technique
  # find the site that each checklist belongs to
  max_site <- 1
  row <- 1
  for(group in unique(proj_cent$site)){
    group_checklists <- proj_cent[proj_cent$site == group,]$checklist_id
    group_df <- WETA_sites[WETA_sites$checklist_id %in% group_checklists,]
    
    vertex_by_sites.DF <- subset(group_df, select = -c(checklist_id))
    
    if(length(group_checklists) > 1){
      dMat <- clusterSimil(subset(vertex_by_sites.DF, select = -c(vertex)), ignoreFirstCol = 0)
      clust_sites_df <- ballsAggr(dMat, vertex_by_sites.DF$vertex, ALPHA = .4, BETA = .5)
      for(s in unique(clust_sites_df$site)){
        clust_sites_df[clust_sites_df$site == s,]$site <- max_site
        max_site <- max_site + 1
      }
      
      # from vertices --> checklists
      for(i in 1:nrow(clust_sites_df)){
        v_by_s_df$vertex[row] <- clust_sites_df[i,]$pvs_vertex
        v_by_s_df$site[row] <- clust_sites_df[i,]$site
        row <- row + 1
      }
    } else {
      v_by_s_df$vertex[row] <- vertex_by_sites.DF$vertex
      v_by_s_df$site[row] <- max_site
      max_site <- max_site + 1
      row <- row + 1
    }
    
  }
  
  og_data$site <- -1
  # trying to join results of independent locations to checklists
  # with the same lat/long
  for(row in 1:nrow(v_by_s_df)){
    v <- v_by_s_df[row,]$vertex
    lat <- og_data[og_data$vertex == v, ]$latitude
    lon <- og_data[og_data$vertex == v, ]$longitude
    og_data[og_data$latitude == lat & og_data$longitude == lon,]$site <- v_by_s_df[row,]$site
  }
  return(og_data)
}


combineMethodsAgg <- function(proj_cent, WETA_sites, og_data){
  v_by_s_df <- data.frame(vertex = numeric(nrow(proj_cent)), site = numeric(nrow(proj_cent)), stringsAsFactors = FALSE)

  # run the agglomerative aggregation technique
  # find the site that each checklist belongs to
  max_site <- 1
  row <- 1
  for(group in unique(proj_cent$site)){
    # disp("group is: ", as.character(group))
    group_checklists <- proj_cent[proj_cent$site == group,]$checklist_id
    group_df <- WETA_sites[WETA_sites$checklist_id %in% group_checklists,]

    vertex_by_sites.DF <- subset(group_df, select = -c(checklist_id))
    if(length(group_checklists) > 1){
      
      dMat <- clusterSimil(subset(vertex_by_sites.DF, select = -c(vertex)), ignoreFirstCol = 0)
      v_map <- data.frame(pvs_vertex = vertex_by_sites.DF$vertex, aggl_v = seq(1:nrow(dMat)))
      agglom_sites <- agglCluster(dMat)
      
      for(cl in agglom_sites){
        
        verts <- v_map[v_map$aggl_v %in% cl@objects,]$pvs_vertex
        
        for(v in verts){
          v_by_s_df$vertex[row] <- v
          v_by_s_df$site[row] <- max_site
          row <- row + 1  
        }
      }
    } else {
      v_by_s_df$vertex[row] <- vertex_by_sites.DF$vertex
      v_by_s_df$site[row] <- max_site
      row <- row + 1  
    }
    max_site <- max_site + 1
  }
  return(v_by_s_df)
}

