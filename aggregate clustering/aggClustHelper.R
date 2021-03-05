# AGG Clustering Helper Fcn

library(gtools)
library(mclustcomp)
library(cluster)

##########
# checklist similarity
#  - returns 1 if two checklists are in the same site
#     in one of the clusterings but in different sites
#     in the other clustering
#  - o/w returns 0
##########
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

agglCluster <- function(dMat, run_modification){
  aggl_sites <- as.data.frame(matrix(data = c(seq(1:nrow(dMat)), seq(1:nrow(dMat))), nrow = nrow(dMat), ncol = 2))
  colnames(aggl_sites) <- c("vertex", "site")
  
  # initialize objects with each being in their own cluster
  CLUSTERS <- list()
  cluster_name <- 1
  
  for(v in aggl_sites$vertex){
    new_c <- new("Cluster", name=cluster_name, objects=c(v), closest=-1, min_dist=-1)
    # CLUSTERS <- append(CLUSTERS, new_c)
    CLUSTERS[[as.character(cluster_name)]] <- new_c
    cluster_name <- cluster_name + 1
  }
  
  n_CLUSTERS <- list()
  cluster_name <- 1
  for(v in CLUSTERS){
    min_obj <- find_min_cluster(v, CLUSTERS, dMat) 
    pop_clst <- new("Cluster", name=cluster_name, objects=c(v@objects), closest=min_obj$clust, min_dist=min_obj$dist)
    # n_CLUSTERS <- append(n_CLUSTERS, pop_clst)
    n_CLUSTERS[[as.character(cluster_name)]] <- pop_clst
    cluster_name <- cluster_name + 1
  }
  
  idx <- find_min_dist(n_CLUSTERS)
  while(n_CLUSTERS[[idx]]@min_dist <= .5){
    new_CLUSTERS <- mergeClusters(n_CLUSTERS, dMat, n_CLUSTERS[[idx]]@name, n_CLUSTERS[[idx]]@closest, cluster_name)
    cluster_name <- cluster_name + 1
    stopifnot(length(new_CLUSTERS) < length(n_CLUSTERS))
    
    
    # the code bwloe needs to be refactored, it is the
    # bottleneck for the rest of the algorithm it
    # takes about 1s to run and there are many runs (~100)
    # when # unique locations is large
    
    # this is a hacky-workaround that is not 100p
    # correct
    
    # recalculate the min distances
    new_CLUSTERS2 <- list()
    for(v in new_CLUSTERS){
      if(run_modification){
        if(is.null(new_CLUSTERS[[as.character(v@closest)]])){
        # if v's objects and the newly combined objects
        # have a non-0 intersection, do this... o/w
        # skip it
        min_obj <- find_min_cluster(v, new_CLUSTERS, dMat) 
        v@closest <- min_obj$clust
        v@min_dist <- min_obj$dist
        # new_CLUSTERS2 <- append(new_CLUSTERS2, v)
        # new_CLUSTERS2[[as.character(v@name)]] <- v
        }
      } else {
        min_obj <- find_min_cluster(v, new_CLUSTERS, dMat) 
        v@closest <- min_obj$clust
        v@min_dist <- min_obj$dist
      }
      new_CLUSTERS2[[as.character(v@name)]] <- v
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
  # clusters <- append(clusters, new_clst)
  clusters[[as.character(new_clust_name)]] <- new_clst
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
prepProj.centers <- function(WETA_2017, covObj, MDD, MIN_OBS, MAX_OBS){
  
  # Looking only at unique locations
  # bc if same location ==> same site (proof)
  # ALTHOUGH TRUE, NOT ALL ALGs MAKE THIS ASSUMPTION
  
  # WETA_uniq <- sqldf("SELECT * FROM WETA_2017 GROUP BY latitude, longitude")
  # WETA_2017$site <- -1
  
  ras.OR2 <- raster("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/occ-cluster/site construction/input data/land cover data/OR_2020_fittedImage.tif")
  prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ################
  # define the checklists as centroids of the buffer
  # reproject centroids into appropriate projection for meters
  # create buffer around each centroid of territory size
  ################
  centroids <- st_as_sf(WETA_2017, coords = c("longitude", "latitude"), crs=prj)
  proj_cent <- st_transform(centroids, crs(ras.OR2))
  buffers <- st_buffer(proj_cent, dist = MDD)
  
  ########
  # for each buffer, find the checklists that overlap
  # if all are unlabeled, create new cluster and add them
  # if some are labeled (possibly some are unlabeled)
  # union them
  ########
  i <- nrow(WETA_2017)
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
  
  return(list(proj_cent, WETA_filtered))
}

# construct many sites and join them together
appendSites <- function(tests, WETA_sites, og_data, covObj=NA, truth_df=NA){
  
  need_L_join <- FALSE
  
  df_to_join <- list()
  
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
      df_to_join <- append(df_to_join, list(eBird_rounded_df))
    }
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
  
  if(!is.na(tests$noisy_gt)){
    prob <- tests$noisy_gt
    nrows <- nrow(truth_df)
    for(row in 1:nrows){
      if(runif(1) <= prob){
        # neighbor <- as.integer(runif(1, -10, 10))
        # if((neighbor + row) > nrows){
        #   neighbor <- -1
        # }
        # if((neighbor + row) < 1){
        #   neighbor <- 1
        # }
        rand_site <- sample(seq(1:nrows), 1)
        truth_df[row,]$site <- truth_df[rand_site,]$site
      }  
    }
    df_to_join <- append(df_to_join, list(truth_df))
  }
  
  if(!is.na(tests$DBSC)){
    DBSC_df <- runDBSC(og_data, covObj)  
    df_to_join <- append(df_to_join, list(DBSC_df))
  }
  
  if(length(tests$clustGeo) > 0){
    for(i in tests$clustGeo){
      disp(i)
      clustGeo_df_i <- clustGeoSites(alpha = i[1], og_data, covObj, num_sites = i[2])
      df_to_join <- append(df_to_join, list(clustGeo_df_i))
    }
  }
  
  if(length(tests$agnes) > 0){
    for(i in 1:length(tests$agnes)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs))
      weta_clust <- agnes(x, method = "ward")
      agnes.clusters <- cutree(as.hclust(weta_clust), k = tests$agnes[[i]])
      cpy$site <- as.numeric(agnes.clusters)
      df_to_join <- append(df_to_join, list(cpy))
    }
  }
  
  if(length(tests$kmeans) > 0){
    for(i in 1:length(tests$kmeans)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs))
      kmean_res <- kmeans(x, tests$kmeans[[i]])
      cpy$site <- as.numeric(kmean_res$cluster)
      df_to_join <- append(df_to_join, list(cpy))  
    }
  }
  
  # need: overlap, interval
  if(length(tests$local) > 0){
    overlap_p <- tests$local$overlap
    interval <- tests$local$interval
    overlap <- as.integer(overlap_p * interval)
    stopifnot(overlap <= interval)
    for(i in seq(from=1, to=nrow(og_data), by=interval)){
      if(i - overlap > 0){
        start <- i - overlap
      } else {
        start <- 1
      }
    #   disp(as.character(start))
    # }
      df <- partialCorrectSites(truth_df, start, interval)
      df_to_join <- append(df_to_join, list(df))  
    }
  }
  
  ######
  # if(need_L_join){
  #   disp("implement nasty left join for eBird")
  #   # # max_site <- 1
  #   # for(row in 1:nrow(WETA_sites)){
  #   #   if(is.na(WETA_sites[row,]$site4)){
  #   #     WETA_sites[row,]$site4 <- max_site
  #   #     max_site <- max_site + 1
  #   #   }
  #   # }
  #   # WETA_sites <- left_join(WETA_sites, eBird_df[c("checklist_id", "site")], by="checklist_id", suffix=c("", "4"))
  # }
  ######
  
  i <- 1
  for(df in df_to_join){
    WETA_sites <- inner_join(WETA_sites, df[c("checklist_id", "site")], by="checklist_id", suffix=c("", as.character(i)))  
    i <- i + 1
  }
  
  return(list(WETA_sites, df_to_join))
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
  # og_data <- inner_join(og_data, WETA_sites[c("checklist_id", "vertex")], by="checklist_id")
  og_data <- inner_join(og_data, v_by_s_df[c("vertex", "site")], by="vertex")
  
  # og_data$site <- -1
  # joining results of independent locations to checklists
  # with the same lat/long
  # for(row in 1:nrow(v_by_s_df)){
  #   v <- v_by_s_df[row,]$vertex
  #   lat <- og_data[og_data$vertex == v, ]$latitude
  #   lon <- og_data[og_data$vertex == v, ]$longitude
  #   og_data[og_data$latitude == lat & og_data$longitude == lon,]$site <- v_by_s_df[row,]$site
  # }
  return(og_data)
}

combineMethodsAgg <- function(proj_cent, WETA_sites, og_data, run_mod=TRUE){
  v_by_s_df <- data.frame(vertex = numeric(nrow(proj_cent)), site = numeric(nrow(proj_cent)), stringsAsFactors = FALSE)
  
  # run the agglomerative aggregation technique
  # find the site that each checklist belongs to
  max_site <- 1
  row <- 1
  for(group in unique(proj_cent$site)){
    group_checklists <- proj_cent[proj_cent$site == group,]$checklist_id
    group_df <- WETA_sites[WETA_sites$checklist_id %in% group_checklists,]

    vertex_by_sites.DF <- subset(group_df, select = -c(checklist_id))
    if(length(group_checklists) > 1){
      # disp("length(group_checklists):", as.character(length(group_checklists)))
      
      dMat <- clusterSimil(subset(vertex_by_sites.DF, select = -c(vertex)), ignoreFirstCol = 0)
      v_map <- data.frame(pvs_vertex = vertex_by_sites.DF$vertex, aggl_v = seq(1:nrow(dMat)))
      agglom_sites <- agglCluster(dMat, run_mod)

      for(cl in agglom_sites){
        
        verts <- v_map[v_map$aggl_v %in% cl@objects,]$pvs_vertex
        
        for(v in verts){
          v_by_s_df$vertex[row] <- v
          v_by_s_df$site[row] <- max_site
          row <- row + 1  
        }
        max_site <- max_site + 1
      }
    } else {
      v_by_s_df$vertex[row] <- vertex_by_sites.DF$vertex
      v_by_s_df$site[row] <- max_site
      row <- row + 1
      max_site <- max_site + 1
    }
  }

  # og_data <- inner_join(og_data, WETA_sites[c("checklist_id", "vertex")], by="checklist_id")
  og_data <- inner_join(og_data, v_by_s_df[c("vertex", "site")], by="vertex")
    
  # og_data$site <- -1
  # for(row in 1:nrow(v_by_s_df)){
  #   v <- v_by_s_df[row,]$vertex
  #   lat <- og_data[og_data$vertex == v, ]$latitude
  #   lon <- og_data[og_data$vertex == v, ]$longitude
  #   og_data[og_data$latitude == lat & og_data$longitude == lon,]$site <- v_by_s_df[row,]$site
  # }
  # 
  return(og_data)
}

calcStats <- function(pred_sites, og_sites, mse){
  cmetrics = c("jaccard", "overlap", "mirkin", "f", "mmm", "nmi1", "nvi")
  ARI <- adjustedRandIndex(og_sites, pred_sites)
  
  m <- mclustcomp(og_sites, pred_sites, types=cmetrics)
  p <- ClusterPurity(as.factor(og_sites), as.factor(pred_sites))
  m_list <- as.list(as.numeric(as.matrix(m)[,2]))
  names(m_list) <- m$types
  return(list(purity=p, ARI=ARI, m_list))
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

runExp <- function(tests, covObj, WETA_sites, og_data, truth_df, proj_cent){
  
  # generate sites!
  sites_obj <- appendSites(tests, WETA_sites, og_data, covObj, truth_df)
  WETA_sites <- sites_obj[[1]]
  dfs_list <- sites_obj[[2]]

  ptm <- proc.time()
  comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
  end <- proc.time() - ptm
  disp("BALLS alg duration: ", as.character(end))

  ptm <- proc.time()
  comb_aggl_df_fast <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = TRUE)
  end <- proc.time() - ptm
  disp("AGGLOM-FAST alg duration: ", as.character(end))
  
  # ptm <- proc.time()
  # comb_aggl_df <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = FALSE)
  # end <- proc.time() - ptm
  # disp("AGGLOM alg duration: ", as.character(end))
  
  # initialize mse list
  mse_list <- list()
  
  # run this through an occupancy model
  consensus_mse <- calcOccMSE(comb_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  mse_list <- append(mse_list, list(consensus_mse))
  
  aggl_consensus_mse <- calcOccMSE(comb_aggl_df_fast, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  mse_list <- append(mse_list, list(aggl_consensus_mse))
  
  # aggl_consensus_mse <- calcOccMSE(comb_aggl_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  # mse_list <- append(mse_list, list(aggl_consensus_mse))
  
  base_mse <- calcOccMSE(truth_df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
  
  # calc mse for each individual method:
  for(df in dfs_list){
    mse.obj <- calcOccMSE(df, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)  
    mse_list <- append(mse_list, list(mse.obj))
  }
  
  return(list(base.mse=base_mse, mse.list=mse_list, pop.WETA.Sites=WETA_sites))
}


# calc similarity btwn each input cluster and the ground truth
makeSIMIL_TO_GT.DF <- function(res_obj, test_names){
  j <- 1
  for(df in res_obj$mse.list){
    o1 <- df$checklists[order(df$checklists$checklist_id),]
    o2 <- res_obj$base.mse$checklists[order(res_obj$base.mse$checklists$checklist_id),]
    stats_list <- calcStats(o1$site, o2$site)
    row_df <- as.data.frame(stats_list, row.names = test_names[j])
    
    if(j == 1){
      comp_to_ground_truth_df <- row_df
    } else {
      comp_to_ground_truth_df <- rbind(comp_to_ground_truth_df, row_df)
    }
    j <- j + 1
  }
  return(comp_to_ground_truth_df)
}

# calc similarity btwn each input cluster and 
# the agglom cluster
makeCLUSTER_COMP.DF <- function(res.obj, test_names){
  agglom_df <- res.obj$mse.list[[2]]
  i <- 1
  for(df in res.obj$mse.list){
    
    og <- agglom_df$checklists[order(agglom_df$checklists$checklist_id),]
    pred <- df$checklists[order(df$checklists$checklist_id),]
    stats_list <- mclustcomp(og$site, pred$site, 
                             types = c("overlap", "nmi1", "f"))
   
    
    ARI <- adjustedRandIndex(og$site, pred$site)
    p <- ClusterPurity(as.factor(og$site), as.factor(pred$site))
    
    m_df <- data.frame(t(unlist(c(stats_list$scores, ARI, p))), row.names = test_names[[i]])
    colnames(m_df) <- c("f", "nmi1", "overlap", "ARI", "purity")
    if(i == 1){
      simil_input_method_df <- m_df
    } else {
      simil_input_method_df <- rbind(simil_input_method_df, m_df)
    }
    i <- i + 1
  }
  return(simil_input_method_df)
}

# put mse values into df
makeMSE.DF <- function(res.obj, test_names){
  
  mse_df <- data.frame(res.obj$base.mse$mse, row.names = "base mse")
  colnames(mse_df) <- c("occ", "det")
  k <- 1
  for(df in res.obj$mse.list){
    mse_row_df <- data.frame(df$mse, row.names = test_names[[k]])
    colnames(mse_row_df) <- c("occ", "det")
    mse_df <- rbind(mse_df, mse_row_df)
    k <- k + 1
  }
  return(mse_df)
}



# function to select a subset of correct sites
partialCorrectSites <- function(df, start, int){
  df[-c(start:(start+int)),]$site <- seq(1:nrow(df[-c(start:(start+int)),]))
  return(df)
}



# st_drop_geometry <- function(x) {
#   if(inherits(x,"sf")) {
#     x <- st_set_geometry(x, NULL)
#     class(x) <- 'data.frame'
#   }
#   return(x)
# }



