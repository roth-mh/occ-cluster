# AGG Clustering Helper Fcn

##########
# checklist similarity
#  - returns 1 if two checklists are in the same site
#     in one of the clusterings but in different sites
#     in the other clustering
#  - o/w returns 0
##########
checklistSimil <- function(chklst1, chklst2){
  if(as.character(chcklst1$site1) == as.character(chcklst2$site1) & 
     as.character(chcklst1$site2) != as.character(chcklst2$site2)){
    return(1)
  }
  if(as.character(chcklst1$site1) != as.character(chcklst2$site1) &
     as.character(chcklst1$site2) == as.character(chcklst2$site2)){
    return(1)
  }
  return(0)
}

##########
# cluster similarity
#  - counts the number of pairwise checklist disagreements
#     between two clusterings
#  - not scalable, O(n^2)
##############################
# clust.DF contains each object and 2 clusterings in
# the columns site1 and site2
##########
clusterSimil <- function(clust.DF){
  comb <- combinations(nrow(clust.DF))
  for(i in 1:nrow(comb)){
    c1 <- comb[i,][[1]]
    c2 <- comb[i,][[2]]
    
    checklistSimil(clust.DF[c1,], clust.DF[c2,])
  }
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