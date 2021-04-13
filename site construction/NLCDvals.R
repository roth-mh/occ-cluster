library(dplyr)
################################
# NLCD Land Classification Types
# https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend
################################

############
# Calculates conductance by finding
#   & dividing by the max resistance
#   of two cells
############
conductanceVals <- function(x){
  num_layers <- length(x)/2
  # THUS FIRST LAYER MUST BE NLCD
  tot_res <- mean(c(getResistVal(x[1]), getResistVal(x[2])))
  for(i in 2:num_layers){
    j <- 2 * i
    tot_res <- tot_res + mean(x[(j - 1)], x[j])
  }
  return(1/tot_res)
}
# should this be an additive or multiplicative relationship?
# additive is probably fine

# both hwy (time of day)
# and river (time of year) have temporal
# aspects that can be incorporated


# resistance values for highway and river buffers
# hwy 50:     3
# hwy 100:    2
# hwy 200:    1

# L riv 50:   3
# L riv 100:  2
# M riv 50:   2.5
# M riv 100:  1.5

# resistance values for NLCD
getResistVal <- function(val){
  if(val == 0){
    # disp("ERROR, no value for cell")
    return(100)
  }
  if(val == 11){
    return(1)
  }
  if(val == 12){
    return(1)
  }
  if(val == 21){
    return(1.5)
  }
  if(val == 22){
    return(3)
  }
  if(val == 23){
    return(4)
  }
  if(val == 24){
    return(5.5)
  }
  if(val == 31){
    return(1)
  }
  if(val == 41){
    return(5)
  }
  if(val == 42){
    return(5)
  }
  if(val == 43){
    return(5)
  }
  if(val == 51){
    return(2.5)
  }
  if(val == 52){
    return(2.5)
  }
  if(val == 71){
    return(2)
  }
  if(val == 72){
    return(2)
  }
  if(val == 73){
    return(1)
  }
  if(val == 74){
    return(1)
  }
  if(val == 81){
    return(2)
  }
  if(val == 82){
    return(2)
  }
  if(val == 90){
    return(3)
  }
  if(val == 95){
    return(2)
  }
}


############################################
# ? # x == 0 ~  "no value?",
# 1 # x == 11 ~ "water",
# 1 # x == 12 ~ "snow",
# 1.5 # x == 21 ~ "developed/open space",
# 3 # x == 22 ~ "developed, low intensity",
# 4 # x == 23 ~ "developed, medium intensity",
# 5.5 # x == 24 ~ "developed, high intensity",
# 1 # x == 31 ~ "barren land (rock/sand/clay)",
# 5 # x == 41 ~ "deciduous forest",
# 5 # x == 42 ~ "evergreen forest",
# 5 # x == 43 ~ "mixed forest",
# 2.5 # x == 51 ~ "dwarf shrub",
# 2.5 # x == 52 ~ "shrub/scrub",
# 2 # x == 71 ~ "grassland/herbaceous",
# 2 # x == 72 ~ "sedge/herbaceous",
# 1 # x == 73 ~ "lichens",
# 1 # x == 74 ~ "moss",
# 2 # x == 81 ~ "pasture/hay",
# 2 # x == 82 ~ "cultivated crops",
# 3 # x == 90 ~ "woody wetlands",
# 2 # x == 95 ~ "emergent herbaceous wetlands"
############################################

# resistanceVals(c(z,95))

# colnames(WETA_2017_region)
