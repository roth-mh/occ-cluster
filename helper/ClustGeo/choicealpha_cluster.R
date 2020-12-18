inertdiss <- function (D, indices=NULL, wt = NULL) {
  n <- as.integer(attr(D, "Size"))
  if (is.null(indices)) indices <- 1:n
  if (is.null(wt)) wt <- rep(1/n,n)
  D <- as.matrix(D)
  if (length(indices) > 1) {
    subD <- D[indices,indices]
    subw <- wt[indices]
    mu <- sum(wt[indices])
    inert <-sweep(subD^2, 1, FUN="*", STATS=subw)
    inert <-sweep(inert, 2, FUN="*", STATS=subw)
    inert <- sum(inert/(2*mu))  
  }
  else {
    inert <- 0
  }
  return(inert)
}

withindiss <- function (D, part, wt = NULL) {
  n <- as.integer(attr(D, "Size"))
  if (is.null(wt)) wt <- rep(1/n,n)
  k <-length(unique(part))
  W <- 0
  for (i in 1:k)
  {
    A <- which(part==i)
    W <- W + inertdiss(D,A,wt)
  }
  return(W)
}

choicealpha <- function(D0, D1, range.alpha,K, wt=NULL,scale=TRUE,graph=TRUE) {
  
  if (is.null(D0)) stop("D0 must be an argument",call.=FALSE)
  if (is.null(D1)) stop("D1 must be an argument",call.=FALSE)
  if (class(D0)!="dist")
    stop("DO must be of class dist (use as.dist)",call.=FALSE)
  if (class(D1)!="dist")
    stop("D1 must be of class dist (use as.dist)",call.=FALSE)
  
  n.alpha<-length(range.alpha)
  n <- as.integer(attr(D1, "Size"))
  if (is.null(n)) 
    stop("invalid dissimilarities",call.=FALSE)
  if (is.na(n) || n > 65536L) 
    stop("size cannot be NA nor exceed 65536",call.=FALSE)
  if (n < 2) 
    stop("must have n >= 2 objects to cluster",call.=FALSE)
  if (!is.null(D1) && length(D0) != length(D1))
    stop("the two dissimilarity structures must have the same size",call.=FALSE)
  if ((max(range.alpha) >1) || (max(range.alpha) < 0))
    stop("Values range.alpha must be in [0,1]",call. = FALSE)
  
  if (scale==TRUE) {
    #scaled dissimilarity matrices 
    D0 <- D0/max(D0)
    D1 <- D1/max(D1)
  }
  
  if (is.null(wt)) 
    wt <- rep(1/n, n)
  
  # within-cluster inertia obtained either with D0 or D1
  W <- matrix(0,length(range.alpha),2)
  rownames(W)  <- paste("alpha=", range.alpha, sep="")
  colnames(W) <- c("W0","W1")
  for (i in 1:length(range.alpha)) {
    tree <- hclustgeo(D0,D1,range.alpha[i],scale=scale,wt=wt)
    part <- cutree(tree,k=K)
    W[i,1] <- withindiss(D0,part,wt)
    W[i,2] <- withindiss(D1,part,wt)
  }
  
  # total inertia obtained with either with D0 or D1
  T0 <-  inertdiss(D0,wt=wt)
  T1 <-  inertdiss(D1,wt=wt) 
  
  # proportion of explained inertia obtained either with D0 orD1 
  Q <- matrix(0,length(range.alpha),2)
  rownames(Q)  <- rownames(W)
  colnames(Q) <- c("Q0","Q1")
  Q[,1] <- 1-W[,1]/T0
  Q[,2] <- 1-W[,2]/T1
  
  # normalized proportion of explained inertia 
  Qnorm <- matrix(0,length(range.alpha),2)
  rownames(Qnorm)  <- rownames(W)
  colnames(Qnorm) <- c("Q0norm","Q1norm")
  Qnorm[,1] <- Q[,1]/Q[1,1]
  Qnorm[,2] <- Q[,2]/Q[length(range.alpha),2]
  
  retlist <- list(Q=Q,Qnorm=Qnorm,range.alpha=range.alpha,K=K)
  class(retlist) <- "choicealpha"
  return(retlist)
}





loadCovariates <- function(file){
  
  covars <- read.delim(file, header=FALSE, sep = ",", strip.white = TRUE)
  
  occ_cov <- lapply(covars[5,], as.character)
  det_cov <- lapply(covars[6,], as.character)
  
  occCov1 <- as.character(occ_cov[1])
  occCov2 <- as.character(occ_cov[2])
  occCov3 <- as.character(occ_cov[3])
  occCov4 <- as.character(occ_cov[4])
  occCov5 <- as.character(occ_cov[5])
  
  detCov1 <- as.character(det_cov[1])
  detCov2 <- as.character(det_cov[2])
  detCov3 <- as.character(det_cov[3])
  detCov4 <- as.character(det_cov[4])
  detCov5 <- as.character(det_cov[5])
  
  siteCovs <- as.character(occ_cov)
  obsCovs <- as.character(det_cov)
  
  return(list(occCov1=occCov1, occCov2=occCov2, occCov3=occCov3, occCov4=occCov4, occCov5=occCov5,
              detCov1=detCov1, detCov2=detCov2, detCov3=detCov3, detCov4=detCov4, detCov5=detCov5,
              siteCovs=siteCovs, obsCovs=obsCovs, det_cov=det_cov, occ_cov=occ_cov))
}



f_in_WETA <- "/nfs/stak/users/rothmark/Documents/siteDef/NORM_weta_data_baseline_merged_2017.csv"
WETA_2017 <- read.delim(f_in_WETA, header=TRUE, sep = ",")

f_in_syn_spec_form <- "/nfs/stak/users/rothmark/Documents/siteDef/syn_species_0_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)


WETA_env_data <- dist(subset(WETA_2017, select = unlist(covObj$occ_cov)))
WETA_geo_data <- dist(subset(WETA_2017, select = c("latitude", "longitude")))
range.alpha <- c(.2, .4, .6, .8)
K <- 400
cr <- choicealpha(WETA_env_data, WETA_geo_data, range.alpha, K, graph=FALSE)

save(cr, file = "/nfs/stak/users/rothmark/Documents/siteDef/cr.Rdata")



