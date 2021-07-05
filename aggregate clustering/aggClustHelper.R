# AGG Clustering Helper Fcn

library(gtools)
library(mclustcomp)
library(cluster)
library(aricode)

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
    new_CLUSTERS_obj <- mergeClusters(n_CLUSTERS, dMat, n_CLUSTERS[[idx]]@name, n_CLUSTERS[[idx]]@closest, cluster_name)
    new_CLUSTERS <- new_CLUSTERS_obj$clusters
    new_clust <- new_CLUSTERS_obj$mCluster
    cluster_name <- cluster_name + 1
    stopifnot(length(new_CLUSTERS) < length(n_CLUSTERS))
    
    
    # the code below needs to be refactored, it is the
    # bottleneck for the rest of the algorithm it
    # takes about 1s to run and there are many runs (~100)
    # when # unique locations is large
    
    # this hacky-workaround is not 100p correct;
    # it only recalculates for clusters that had 
    # a closest neighbor that was merged 
    
    # recalculate the min distances
    new_CLUSTERS2 <- list()
    for(v in new_CLUSTERS){
      if(run_modification){
        
        # if v's objects and the newly combined objects
        # have a non-0 intersection, find new_cluster_min
        if(is.null(new_CLUSTERS[[as.character(v@closest)]])){
          min_obj <- find_min_cluster(v, new_CLUSTERS, dMat) 
          v@closest <- min_obj$clust
          v@min_dist <- min_obj$dist
          # new_CLUSTERS2 <- append(new_CLUSTERS2, v)
          # new_CLUSTERS2[[as.character(v@name)]] <- v
        } else {
          # only need to check if the distance from 
          # the current cluster, v, to the merged cluster
          # is less than what is currently in min_dist
          
          # in worst case, O(n^2) still but in practice O(n)?
          
          # need to check if v == new_clust
          if(v@name != new_clust@name){
            new_dist <- calcMeanDist(v, new_clust, dMat)
            if(v@min_dist > new_dist){
              v@min_dist <- new_dist
              v@closest <- new_clust@name
            }  
          }
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
  return(list(clusters=clusters, mCluster=new_clst))
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
# TODO: order matters ... but it shouldn't
appendSites <- function(tests, WETA_sites, og_data, covObj=NA, truth_df=NA){
  
  need_L_join <- FALSE
  
  df_to_join <- list()
  
  if(length(tests$agnes) > 0){
    for(i in 1:length(tests$agnes)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs, "latitude", "longitude"))
      weta_clust <- agnes(x, method = "ward") # I think complete linkage makes more sense in the context of sites
      agnes.clusters <- cutree(as.hclust(weta_clust), k = tests$agnes[[i]])
      cpy$site <- as.numeric(agnes.clusters)
      # df_to_join <- append(df_to_join, list(cpy))
      df_to_join[[paste0("agnes-", as.character(tests$agnes[[i]]))]] <- cpy
    }
  }
  
  if(length(tests$clustGeo) > 0){
    for(i in tests$clustGeo){
      disp(i)
      clustGeo_df_i <- clustGeoSites(alpha = i[1], og_data, covObj, num_sites = i[2])
      # df_to_join <- append(df_to_join, list(clustGeo_df_i))
      df_to_join[[paste0("clustGeo-", as.character(i[1]), "-", as.character(i[2]))]] <- clustGeo_df_i
    }
  }
  
  if(!is.na(tests$DBSC)){
    DBSC_df <- runDBSC(og_data, covObj)  
    # df_to_join <- append(df_to_join, list(DBSC_df))
    df_to_join[["DBSC"]] <- DBSC_df
  }
  
  if("all_svs" %in% names(tests)){
    all_svs_df <- og_data
    all_svs_df$site <- seq(1:nrow(all_svs_df))
    # df_to_join <- append(df_to_join, list(DBSC_df))
    df_to_join[["all_svs"]] <- all_svs_df
  }
  
  # if(!is.na(tests$eBird)){
  #   eBird_df <- filter_repeat_visits(
  #     og_data,
  #     min_obs = 2,
  #     max_obs = 10,
  #     annual_closure = TRUE,
  #     date_var = "formatted_date",
  #     site_vars = c("locality_id", "observer_id")
  #   )
  #   df_to_join <- append(df_to_join, list(eBird_df))
  #   need_L_join <- TRUE
  # }
  
  if(!is.na(tests$eBird_simple)){
    eBird_simple_df <- filter_repeat_visits(
      og_data,
      min_obs = 1,
      max_obs = 1000000,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id")
    )
    # df_to_join <- append(df_to_join, list(eBird_simple_df))
    df_to_join[["eBird_simple"]] <- eBird_simple_df
  }
  
  
  if(length(tests$kmeans) > 0){
    for(i in 1:length(tests$kmeans)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs, "latitude", "longitude"))
      kmean_res <- kmeans(x, tests$kmeans[[i]])
      cpy$site <- as.numeric(kmean_res$cluster)
      # df_to_join <- append(df_to_join, list(cpy))
      df_to_join[[paste0("kmeans-", as.character(tests$agnes[[i]]))]] <- cpy
    }
  }
  
  if(length(tests$kmSq) > 0){
    for(i in 1:length(tests$kmSq)){
      kmSq.df <- kmsq.Sites(og_data, rad_m = tests$kmSq[[i]])
      # df_to_join <- append(df_to_join, list(kmSq.df))
      df_to_join[[paste0("kmSq-", as.character(tests$kmSq[[i]]))]] <- kmSq.df
    }
  }
  
  
  if(length(tests$rounded) > 0){
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
      # df_to_join <- append(df_to_join, list(eBird_rounded_df))
      df_to_join[[paste0("rounded-", as.character(i))]] <- eBird_rounded_df
    }
  }
  
  if(!is.na(tests$noisy_gt)){
    ######
    # prob <- tests$noisy_gt
    # nrows <- nrow(truth_df)
    # for(row in 1:nrows){
    #   if(runif(1) <= prob){
    #     # neighbor <- as.integer(runif(1, -10, 10))
    #     # if((neighbor + row) > nrows){
    #     #   neighbor <- -1
    #     # }
    #     # if((neighbor + row) < 1){
    #     #   neighbor <- 1
    #     # }
    #     rand_site <- sample(seq(1:nrows), 1)
    #     truth_df[row,]$site <- truth_df[rand_site,]$site
    #   }  
    # }
    # df_to_join <- append(df_to_join, list(truth_df))
    ######
    stopifnot(1==0)
  }

  
  # need: overlap, interval
  if(length(tests$local) > 0){
    #####
    # overlap_p <- tests$local$overlap
    # interval <- tests$local$interval
    # overlap <- as.integer(overlap_p * interval)
    # stopifnot(overlap <= interval)
    # for(i in seq(from=1, to=nrow(og_data), by=interval)){
    #   if(i - overlap > 0){
    #     start <- i - overlap
    #   } else {
    #     start <- 1
    #   }
    # #   disp(as.character(start))
    # # }
    #   df <- partialCorrectSites(truth_df, start, interval)
    #   df_to_join <- append(df_to_join, list(df))  
    # }
    #####
    stopifnot(1==0)
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
    #changing join on checklist_id to vertex because checklist_id may no longer be unique (boot exp)
    WETA_sites <- inner_join(WETA_sites, df[c("vertex", "site")], by="vertex", suffix=c("", as.character(i)))  
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
      # disp(paste0("# of rows v_by_sites.DF: ", nrow(vertex_by_sites.DF)))
      # disp(paste0("# of rows group_df: ", nrow(group_df)))
      # disp(paste0("# of group_checklists: ", length(group_checklists)))
      # if(length(group_checklists) != nrow(vertex_by_sites.DF)){
      #   disp("STOP HERE!!!!")
      #   disp("STOP HERE!!!!")
      #   disp("STOP HERE!!!!")
      # }
      
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
      
      # if run_mod == TRUE, we make a simplifying assumption that does not hold
      # TODO: refactor and make quicker
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


  og_data <- inner_join(og_data, v_by_s_df[c("vertex", "site")], by="vertex")
    
  return(og_data)
}


# examines marginal difference from 
calcStats <- function(pred_df, og_df){
  # pred_df_uniq <- sqldf("SELECT * from pred_df GROUP BY latitude, longitude")
  # og_df_uniq <- sqldf("SELECT * from og_df GROUP BY latitude, longitude")
  # 
  # pred_df_uniq <- pred_df_uniq[order(pred_df_uniq$checklist_id),]
  # og_df_uniq <- og_df_uniq[order(og_df_uniq$checklist_id),]
  
  pred_df <- pred_df[order(pred_df$checklist_id),]
  og_df <- og_df[order(og_df$checklist_id),]
  
  pred_sites <- pred_df$site
  og_sites <- og_df$site
  # cmetrics = c("jaccard", "overlap", "mirkin", "f", "mmm", "nmi1", "nvi")
  ARI <- adjustedRandIndex(og_sites, pred_sites)
  
  # m <- mclustcomp(og_sites, pred_sites, types=cmetrics)
  nmi <- NMI(og_sites, pred_sites, variant = "joint")
  ami <- AMI(og_sites, pred_sites)
  nid <- NID(og_sites, pred_sites)
  nvi <- NVI(og_sites, pred_sites)
  
  p <- ClusterPurity(as.factor(pred_sites), as.factor(og_sites))
  return(list(purity=p, ARI=ARI, NMI=nmi, AMI=ami, NID=nid, NVI=nvi))
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}


# function to calculate the clusterings for each test specified in 
# variable `tests`
runExp <- function(tests, covObj, WETA_sites, og_data, truth_df, proj_cent, comb_df=NA, comb_aggl_df_fast=NA){
  
  sites_obj <- appendSites(tests, WETA_sites, og_data, covObj, truth_df)
  WETA_sites <- sites_obj[[1]]
  dfs_list <- sites_obj[[2]]
  
  # generate sites!
  if(is.na(comb_aggl_df_fast) & is.na(comb_df)){
    #TODO: 4/13 speed up?
    ptm <- proc.time()
    comb_df <- combineMethods(proj_cent, WETA_sites, og_data)
    end <- proc.time() - ptm
    disp("BALLS alg duration: ", as.character(end))
  
    ptm <- proc.time()
    # TODO: fix run_mod
    comb_aggl_df_fast <- combineMethodsAgg(proj_cent, WETA_sites, og_data, run_mod = TRUE)
    end <- proc.time() - ptm
    disp("AGGLOM-UPDATED alg duration: ", as.character(end))
  }
  
  # initialize mse list
  sites_list <- list()
  
  sites_list[["balls"]] <- comb_df
  sites_list[["agglom"]] <- comb_aggl_df_fast
  sites_list[["base"]] <- truth_df
  
  
  # calc mse for each individual method:
  for(i in 1:length(dfs_list)){
    # mse.obj <- calcOccMSE(dfs_list[[i]], covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF, syn_spec = T)
    sites_list[[names(dfs_list[i])]] <- dfs_list[[i]]
  }
  
  return(list(
    # base.mse=base_mse,
    sites_list=sites_list, 
    pop.WETA.Sites=WETA_sites, 
    balls_df=comb_df, 
    aggl_df=comb_aggl_df_fast)
    )
}


# calc similarity btwn each input cluster and the ground truth
makeSIMIL_TO_GT.DF <- function(res_obj, test_names){
  j <- 1
  # o2 <- res_obj$base.mse$checklists[order(res_obj$base.mse$checklists$checklist_id),]
  for(df in res_obj$mse.list){
    # o1 <- df$checklists[order(df$checklists$checklist_id),]
    stats_list <- calcStats(df$checklists, res_obj$base.mse$checklists)
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


makeINPUT_SIMUL.DF <- function(res_obj){
  # TODO: HARDCODED skipping agglom and balls alg results
  num_algs <- length(res_obj$mse.list) - 2
  pairwise_comb <- combinations(num_algs, 2)
  for(i in 1:nrow(pairwise_comb)){
    first_alg <- pairwise_comb[i,][[1]]
    second_alg <- pairwise_comb[i,][[2]]
  
    stats_list <- calcStats(res_obj$mse.list[[first_alg]]$checklists, res_obj$mse.list[[second_alg]]$checklists)
    if(i == 1){
      df <- as.data.frame(stats_list)  
    } else {
      df <- df + stats_list
    }
  }
  df <- df/num_algs
  return(df)
}



# calc similarity btwn each input cluster and 
# the agglom cluster (DEPRECATED ... I THINK)
makeCLUSTER_COMP.DF <- function(res.obj, test_names){
  # TODO: HARDCODED agglom-fast df
  agglom_df <- res.obj$mse.list[[2]]
  i <- 1
  for(df in res.obj$mse.list){
    # og <- agglom_df$checklists[order(agglom_df$checklists$checklist_id),]
    # pred <- df$checklists[order(df$checklists$checklist_id),]
    m_df <- calcStats(pred, og)
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


# helper function to link Checklists given a reduced
# df (red_df) and a full df
linkChecklists <- function(red_df, full_df){
  #######
  # for(row in 1:nrow(full_df)){
  #   # disp(row)
  #   long <- full_df[row, ]$longitude
  #   lat <- full_df[row, ]$latitude
  #   site_obj <- red_df[red_df$longitude == long & red_df$latitude == lat,]
  #   full_df[row,]$site <- as.character(site_obj$site)
  # }
  #######
  for(row in 1:nrow(red_df)){
    # disp(row)
    long <- red_df[row, ]$longitude
    lat <- red_df[row, ]$latitude
    full_df[full_df$longitude == long & full_df$latitude == lat,]$site <- red_df[row,]$site
  }
  
  return(full_df)
}


# fcn to make test flow smoother
# test_names is a list of strings specifying
# the tests you want to run. the parameters
# of each test are separated with a '-'
#
# for example: clustGeo-.8-850, kmSq-1000, rounded-4
genTests <- function(test_names){
  tests <- list(rounded=list(), 
                base=NA,
                eBird=NA, 
                eBird_simple=NA,
                eBird_lower=NA,
                eBird_upper=NA,
                kmSq=list(),
                DBSC=NA,
                GT=NA,
                noisy_gt=NA,
                clustGeo=list(),
                agnes=list(),
                kmeans=list()
                # local=list()
  )
  for(i in 3:length(test_names)){
    t_name <- strsplit(test_names[[i]], "-")[[1]]
    
    # Need to update if we use local
    if(t_name[1] %in% c("clustGeo", "agnes", "kmeans", "local", "kmSq", "rounded")){
      
      len <- length(tests[[t_name[1]]]) + 1
      tests[[t_name[1]]][len] <- list(as.double(t_name[2:length(t_name)]))
    } else {
      tests[[t_name[[1]][1]]] <- T  
    }
  }
  
  return(tests)
}


# Clustering object that holds the test_name
# the mse values, the checklist dataframe
# and the similarity to GT data frame 
# (generated by clusterStats)
setClass("Clustering", 
         slots = list(
           name="character", 
           occ.mse="list", 
           det.mse="list",
           prob.occ.mse="list", 
           prob.det.mse="list",
           svs="numeric",
           acc.svs="numeric",
           checklist_df="data.frame",
           simil.GT="list")
)


runStabilityExp <- function(exp, truth_df, WETA_2017, WETA_filtered, proj_cent, tests, test_names, deterministic = FALSE){
  
  # list_of_res <- vector("list", length(test_names))
  
  truth_df_list <- list(truth_df)
  
  if(exp == "boot"){
    sampleSizes <- c(nrow(truth_df))
    # sampleSizes <- c(200)
    # numExps <- 2
    numExps <- 10
    replacement <- TRUE
    gen.new.df <- FALSE
  } else if(exp == "location"){
    WETA_uniq <- sqldf("SELECT * FROM WETA_filtered GROUP BY latitude, longitude")
    sampleSizes <- c(200, 400, 600, nrow(WETA_uniq))
    numExps <- 10
    # sampleSizes <- c(100, 120)
    # numExps <- 1
    replacement <- FALSE
    gen.new.df <- TRUE
  } else if(exp == "normal"){
    sampleSizes <- c(round(nrow(truth_df)/2), round(nrow(truth_df)/3), round(nrow(truth_df)/4))
    numExps <- 10
    # sampleSizes <- c(100, 120)
    # numExps <- 2
    replacement <- FALSE
    gen.new.df <- TRUE
  }
  
  for(sampSize in sampleSizes){
    aggl_df <- NA
    balls_df <- NA
    
    results <- list()
    for(i in 1:numExps){
      disp(paste0("running ", as.character(exp), " experiment; sample size: ", as.character(sampSize), " and exp #", as.character(i)))
      
      # use the same underlying dfs for each sample size
      if(length(truth_df_list) > 1){
        # head(truth_df_list[[i]])
        truth_df <- truth_df_list[[i]]
      } else {
        truth_df <- truth_df_list[[1]]
      }
      
      og_data <- subset(WETA_filtered_uniq, select = -c(site))
      
      og_data <- og_data[order(og_data$checklist_id),]
      truth_df <- truth_df[order(truth_df$checklist_id),]
      
      og_data$species_observed_syn <- truth_df$species_observed_syn
      
      WETA_sites <- subset(og_data, select = c(checklist_id, vertex))
      
      if(exp == "location"){
        # NOTE: each exp may have different total number of checklists because we are only sampling locations
        samp <- sample(seq(1:nrow(WETA_uniq)), sampSize, replace = replacement)
        W_uniq <- WETA_uniq[samp,]
        W_full <- linkChecklists(W_uniq, WETA_2017)
        W_full <- W_full[!(W_full$site == -1),]
        W_sites <- subset(W_full, select = c(checklist_id, vertex))
        
        og_d <- og_data[og_data$checklist_id %in% W_sites$checklist_id,]
        t_df <- truth_df[truth_df$checklist_id %in% W_sites$checklist_id,]
        p_cent <- proj_cent[proj_cent$checklist_id %in% W_sites$checklist_id,]
      } else {
        samp <- sample(seq(1:nrow(WETA_2017)), sampSize, replace = replacement)
        W_sites <- WETA_sites[samp,]
        og_d <- og_data[samp,]
        t_df <- truth_df[samp,]
        p_cent <- proj_cent[samp,]

        # disp(og_d$checklist_id == W_sites$checklist_id)
        
        W_sites <- W_sites[order(W_sites$checklist_id),]
        og_d <- og_d[order(og_d$checklist_id),]
        t_df <- t_df[order(t_df$checklist_id),]
        p_cent <- p_cent[order(p_cent$checklist_id),]
        
        
        if(replacement){
          id <- 1
          # get the checklist_ids of the duplicated rows
          chklsts <- unique(W_sites[ave(W_sites$vertex, W_sites$checklist_id, FUN = length) > 1, ]$checklist_id)
          to_change <- list()
          
          # needed to run w/o error on the HPC server
          W_sites$checklist_id <- as.character(W_sites$checklist_id)
          og_d$checklist_id <- as.character(og_d$checklist_id)
          W_sites <- W_sites[order(W_sites$checklist_id),]
          og_d <- og_d[order(og_d$checklist_id),]
          t_df <- t_df[order(t_df$checklist_id),]
          p_cent <- p_cent[order(p_cent$checklist_id),]
          
          for(k in 1:nrow(W_sites)){
            if(W_sites[k,]$checklist_id %in% chklsts){
              if(W_sites[k,]$checklist_id %in% to_change){
                new_vert <- max(W_sites$vertex) + 1
                # new_cid <- as.character(paste0("C", as.character(id)))
                # disp(paste0("new c_id: ", new_cid, " new_vert ", new_vert))
                disp(paste0("current c_id W_sites: ", W_sites[k,]$checklist_id))
                
                W_sites[k,]$checklist_id <- as.character(paste0("C", as.character(id)))
                W_sites[k,]$vertex <- new_vert
                
                og_d[k,]$checklist_id <- as.character(paste0("C", as.character(id)))
                og_d[k,]$vertex <- new_vert
                
                t_df[k,]$checklist_id <- as.character(paste0("C", as.character(id)))
                
                p_cent[k,]$checklist_id <- as.character(paste0("C", as.character(id)))
                id <- id + 1

              } else {
                to_change <- append(to_change, as.character(W_sites[k,]$checklist_id))
              }
            }
          }
          ######
          # for(c_id in chklsts){
          #   # rows <- W_sites[W_sites$checklist_id == c_id,]
          #   for(j in 1:(nrow(W_sites[W_sites$checklist_id == c_id,]) - 1)){
          #     new_vert <- max(W_sites$vertex) + 1
          #     W_sites[W_sites$checklist_id == c_id,][j,]$checklist_id <- as.character(paste0("C", as.character(id)))
          #     W_sites[W_sites$checklist_id == c_id,][j,]$vertex <- new_vert
          #     
          #     og_d[og_d$checklist_id == c_id,][j,]$checklist_id <- as.character(paste0("C", as.character(id)))
          #     og_d[og_d$checklist_id == c_id,][j,]$vertex <- new_vert
          #     
          #     p_cent[p_cent$checklist_id == c_id,][j,]$checklist_id <- as.character(paste0("C", as.character(id)))
          #     
          #     id <- id + 1
          #   }
          # }
        }
        
      }
      
      # W_sites$vertex <- seq(1:nrow(W_sites))
      
      stopifnot(og_d$checklist_id == t_df$checklist_id)
      # we can feed in consensus clusterings because all of the algs are deterministic (should feed this in as a param)
      res_obj <- runExp(tests, covObj, W_sites, og_d, t_df, p_cent, comb_df = balls_df, comb_aggl_df_fast = aggl_df)
      
      if(deterministic){
        aggl_df <- res_obj$aggl_df
        balls_df <- res_obj$balls_df
      }
      
      for(j in 1:length(res_obj$sites_list)){
        clustering <- res_obj$sites_list[j][[1]]
        t_name <- names(res_obj$sites_list)[[j]]
        
        if(!is.na(clustering)){
          res <- clusterStats(clustering, t_df, t_name, covObj, TRUE_OCC_COEFF, TRUE_DET_COEFF)
          if(i == 1){
            results[[t_name]] <- new("Clustering",
                                       name=t_name,
                                       occ.mse=res$mse$occ,
                                       det.mse=res$mse$det,
                                       svs=res$svs,
                                       acc.svs=res$acc.svs,
                                       checklist_df=clustering,
                                       simil.GT=res$simil.GT)
          } else {
            results[[t_name]]@occ.mse <- results[[t_name]]@occ.mse + res$mse$occ
            results[[t_name]]@det.mse <- results[[t_name]]@det.mse + res$mse$det
            results[[t_name]]@svs <- results[[t_name]]@svs + res$svs
            results[[t_name]]@acc.svs <- results[[t_name]]@acc.svs + res$acc.svs
            results[[t_name]]@simil.GT <- results[[t_name]]@simil.GT + res$simil.GT
          }
        }
      }
      
      if(gen.new.df){
        TRUE_OCC_COEFF <- runif(6, -1.5, 1.5)
        TRUE_DET_COEFF <- runif(6, -1.5, 1.5)
        truth_df <- populateDF(WETA_filtered, covObj$siteCovs, covObj$obsCovs, unique(WETA_filtered$site), TRUE_OCC_COEFF, TRUE_DET_COEFF)
        truth_df_list <- append(truth_df_list, list(truth_df))
      }
      
    }
    
    z <- combineIntoDF(results, numExps)
    write.csv(z, paste0("stability-", exp, "-", as.character(sampSize), "-", as.character(numExps), "runs.csv"))
    
    # each sample size uses same dfs
    # (stored in truth_df_list)
    gen.new.df <- FALSE
  }
  
}






# experiment to test similarity and mse of sp. clustering algs
# against the ground truth (no consensus clustering)
#
# this runs the clustering aspect of most of the algorithms.
# evaluation occurs at a later step
baselineExp <- function(tests, og_data, covObj, truth_df){
  
  results <- list()
  results[["base"]] <- truth_df
  
  if(length(tests$rounded) > 0){
    WETA_2017_i <- roundLatLong(og_data, 4)
    eBird_rounded_df <- filter_repeat_visits(
      WETA_2017_i,
      min_obs = 1,
      max_obs = 1000000,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("rounded_locality_id")
    )
    results[[paste0("rounded-", as.character(4))]] <- eBird_rounded_df
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
    results[["eBird"]] <- eBird_df
  }
  
  if(!is.na(tests$eBird_lower)){
    eBird_df <- filter_repeat_visits(
      og_data,
      min_obs = 2,
      max_obs = 1000000,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id", "observer_id")
    )
    results[["eBird_lower"]] <- eBird_df
  }
  
  if(!is.na(tests$eBird_upper)){
    eBird_df <- filter_repeat_visits(
      og_data,
      min_obs = 1,
      max_obs = 10,
      annual_closure = TRUE,
      date_var = "formatted_date",
      site_vars = c("locality_id", "observer_id")
    )
    results[["eBird_upper"]] <- eBird_df
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
    results[["eBird_simple"]] <- eBird_simple_df
  }
  

  if(!is.na(tests$DBSC)){
    DBSC_df <- runDBSC(og_data, covObj)
  
    results[["DBSC"]] <- DBSC_df
  }
  
  if(length(tests$clustGeo) > 0){
    for(i in tests$clustGeo){
      clustGeo_df_i <- clustGeoSites(alpha = i[1], og_data, covObj, num_sites = i[2])
      results[[paste0("clustGeo-", as.character(i[1]), "-", as.character(i[2]))]] <- clustGeo_df_i
    }
  }
  
  if(length(tests$agnes) > 0){
    for(i in 1:length(tests$agnes)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs, "latitude", "longitude"))
      weta_clust <- agnes(x, method = "ward")
      agnes.clusters <- cutree(as.hclust(weta_clust), k = tests$agnes[[i]])
      cpy$site <- as.numeric(agnes.clusters)
      results[[paste0("agnes-", as.character(tests$agnes[[i]]))]] <- cpy
    }
  }
  
  if(length(tests$kmeans) > 0){
    for(i in 1:length(tests$kmeans)){
      cpy <- og_data
      x <- subset(cpy, select = c(covObj$siteCovs, covObj$obsCovs, "latitude", "longitude"))
      kmean_res <- kmeans(x, tests$kmeans[[i]])
      cpy$site <- as.numeric(kmean_res$cluster)
      results[[paste0("kmeans-", as.character(tests$kmeans[[i]]))]] <- cpy
    }
  }
  
  if(length(tests$kmSq) > 0){
    for(i in 1:length(tests$kmSq)){
      kmSq.df <- kmsq.Sites(og_data, rad_m = tests$kmSq[[i]])
      results[[paste0("kmSq-", as.character(tests$kmSq[[i]]))]] <- kmSq.df
    }
  }
  return(results)
}


calcAccSingleVisitSites <- function(checklist_df, truth_df){
  visit.freq <- data.frame(table(checklist_df$site))
  single_visit_sites <- visit.freq[visit.freq$Freq == 1,]$Var1
  chklsts <- checklist_df[checklist_df$site %in% single_visit_sites,]$checklist_id
  
  true.visit.freq <- data.frame(table(truth_df$site))
  true.single_visit_sites <- true.visit.freq[true.visit.freq$Freq == 1,]$Var1
  true.chklsts <- truth_df[truth_df$site %in% true.single_visit_sites,]$checklist_id
  
  acc.svs <- intersect(chklsts, true.chklsts)
  
  return(length(acc.svs))
}


# prob_MSE is MSE calculated as the average difference between true probability and predicted probability
clusterStats <- function(checklist_df, truth_df, og_data, test_name, covObj, occ_coeff, det_coeff, full_df=data.frame(), prob_MSE=FALSE){
  
  stopifnot(og_data$checklist_id == truth_df$checklist_id)
  # stopifnot(checklist_df$checklist_id == truth_df$checklist_id)
  checklist_df <- subset(checklist_df, select = -c(species_observed_syn))
  checklist_df <- sqldf("SELECT c.*, o.species_observed_syn FROM checklist_df c 
                       INNER JOIN og_data o on c.checklist_id == o.checklist_id")

  # checklist_df$species_observed_syn <- og_data$species_observed_syn
  
  # if(nrow(full_df) > 0){
  #   checklist_df <- linkChecklists(checklist_df, full_df)
  #   truth_df <- linkChecklists(truth_df, full_df)
  # }
  if(prob_MSE){
    t_df <- truth_df
  } else {
    t_df <- data.frame()
  }
  mse.obj <- calcOccMSE(
    checklist_df,
    covObj,
    true_occ_coefficients = occ_coeff,
    true_det_coefficients = det_coeff,
    syn_spec = T,
    truth_df = t_df)
  
  # mse.obj <- list()
  # mse.obj$mse <- list(occ=0,det=0)
  
  svs <- data.frame(table(table(checklist_df$site)))["1",]$Freq
  if(as.character(data.frame(table(table(checklist_df$site)))["1",]$Var1) != "1"){
    svs <- 0
  }
  
  acc.svs <- calcAccSingleVisitSites(checklist_df, truth_df)
  
  if(length(checklist_df$site) != length(truth_df$site)){
    # TODO: hard coded CLUSTER METRICS HERE
    simil.GT <- base::as.data.frame(matrix(data=0, nrow=1,ncol=6, dimnames = list(test_name, c("purity", "ARI", "NMI", "AMI", "NID", "NVI"))))
  } else {
    simil.GT <- calcStats(checklist_df, truth_df)
    simil.GT <- as.data.frame(simil.GT, row.names = test_name)  
  }
  
  return(list(mse=mse.obj$mse, prob_mse=mse.obj$prob_mse, simil.GT=simil.GT, svs=svs, acc.svs=acc.svs))
}


combineIntoDF <- function(results.obj, numExps){
  res <- NA
  for(i in 1:length(results.obj)){
    cl <- results.obj[[i]]
    if(!is.na(cl)){
      row_df <- cbind(cl@simil.GT, occ=cl@occ.mse, det=cl@det.mse, svs=cl@svs, acc_svs=cl@acc.svs)
      if(is.na(res)){
        res <- row_df
      } else {
        res <- rbind(res, row_df)
      }
    }
  }
  z <- res/numExps
  return(z)
}



# given a df with some svs, filter
# out all checklists with the same lat/long as these svs
filterOutSVS <- function(df){
  
  # GET THE SINGLE VISIT SITES
  WETA_sites <- sqldf("SELECT site, count(site) visits from df GROUP BY site")
  WETA_svs <- WETA_sites[WETA_sites$visits == 1,]
  # disp("number svs: ", nrow(WETA_svs))
  WETA_svs_updated <- sqldf("SELECT latlong, WETA_svs.* FROM WETA_svs JOIN df on df.site == WETA_svs.site")
  
  # SELECT ONLY THE CHECKLISTS W/ a LATLONG DIFFERENT FROM THE SET OF SVS
  WETA_filter <- df[!(df$latlong %in% WETA_svs_updated$latlong),]
    
  # sqldf("SELECT * FROM df WHERE latlong NOT IN (SELECT latlong FROM WETA_svs_updated)")
  # nrow(WETA_filter)
  return(WETA_filter)
}








