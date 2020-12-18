
#######
# set up the DS for holding occ/det values
#######
lst_clust <- c()
names_clust <- c()

lst_skate <- c()
names_skate <- c()
i <- 0
num_sites_seq <- seq(200,400,20)
for(n_s in num_sites_seq){
  i <- i + 1
  lst_clust[[i]] <- list(occ=0, det=0)
  lst_skate[[i]] <- list(occ=0, det=0)
  names_clust <- append(names_clust, as.character(n_s))
  names_skate <- append(names_skate, as.character(n_s))
}


#######
# run MSE on clust and SKATER with variable number of sites
#######
numIters <- 10
for(iter in 1:numIters){
  disp("running ITERATION: ", iter)
  i <- 0
  for(n_s in num_sites_seq){
    disp("now running: ", n_s)
    i <- i + 1
    skate.MSE <- calcSKATER.MSE(WETA_2017, covObj, numSites = n_s)
    clust.MSE <- calcClustGeoMSE(.8, WETA_2017, covObj, num_sites = n_s)
    
    lst_clust[[i]]$det <- lst_clust[[i]]$det + clust.MSE$mse$det
    lst_clust[[i]]$occ <- lst_clust[[i]]$occ + clust.MSE$mse$occ
    
    lst_skate[[i]]$det <- lst_skate[[i]]$det + skate.MSE$mse$det
    lst_skate[[i]]$occ <- lst_skate[[i]]$occ + skate.MSE$mse$occ
  }
}

#######
# misc analysis
#######

k<-0
for(m in num_sites_seq){
  k <- k+1
  lst_clust[[k]]$det <- lst_clust[[k]]$det/numIters
  lst_clust[[k]]$occ <- lst_clust[[k]]$occ/numIters
  
  lst_skate[[k]]$det <- lst_skate[[k]]$det/numIters
  lst_skate[[k]]$occ <- lst_skate[[k]]$occ/numIters
}

res1 <- combineDF(lst_clust, names_clust)
res2 <- combineDF(lst_skate, names_skate)
plot(x =num_sites_seq, y=res1$occ, main = "clustGeo avg occ X num sites")
plot(x =num_sites_seq, y=res1$det, main = "clustGeo avg det X num sites")
plot(x =num_sites_seq, y=res2$occ, main = "SKATER avg occ X num sites")
plot(x =num_sites_seq, y=res2$det, main = "SKATER avg det X num sites")
res1
res2

write.csv(res1, file = "eBird/site def/ClustXnumSites.csv")
write.csv(res2, file = "eBird/site def/SKATERxNumSites.csv")
