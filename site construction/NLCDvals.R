library(dplyr)
################################
# NLCD Land Classification Types
# https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend
################################

resistanceVals <- function(x){
  return(max(getResistVal(x[1]), getResistVal(x[2])))
}

getResistVal <- function(val){
  if(val == 0){
    disp("ERROR, no value for cell")
    return(-1)
  }
  if(val == 11){
    return(3)
  }
  if(val == 12){
    return(1)
  }
  if(val == 21){
    return(2)
  }
  if(val == 22){
    return(3)
  }
  if(val == 23){
    return(4)
  }
  if(val == 24){
    return(5)
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
    return(3)
  }
  if(val == 52){
    return(3)
  }
  if(val == 71){
    return(3)
  }
  if(val == 72){
    return(3)
  }
  if(val == 73){
    return(2)
  }
  if(val == 74){
    return(2)
  }
  if(val == 81){
    return(1)
  }
  if(val == 82){
    return(1)
  }
  if(val == 90){
    return(3)
  }
  if(val == 95){
    return(3)
  }
}


############################################
# ? # x == 0 ~  "no value?",
# 3 # x == 11 ~ "water",
# 1 # x == 12 ~ "snow",
# 2 # x == 21 ~ "developed/open space",
# 3 # x == 22 ~ "developed, low intensity",
# 4 # x == 23 ~ "developed, medium intensity",
# 5 # x == 24 ~ "developed, high intensity",
# 1 # x == 31 ~ "barren land (rock/sand/clay)",
# 5 # x == 41 ~ "deciduous forest",
# 5 # x == 42 ~ "evergreen forest",
# 5 # x == 43 ~ "mixed forest",
# 3 # x == 51 ~ "dwarf shrub",
# 3 # x == 52 ~ "shrub/scrub",
# 3 # x == 71 ~ "grassland/herbaceous",
# 3 # x == 72 ~ "sedge/herbaceous",
# 2 # x == 73 ~ "lichens",
# 2 # x == 74 ~ "moss",
# 1 # x == 81 ~ "pasture/hay",
# 1 # x == 82 ~ "cultivated crops",
# 3 # x == 90 ~ "woody wetlands",
# 3 # x == 95 ~ "emergent herbaceous wetlands"
############################################

# resistanceVals(c(z,95))

# colnames(WETA_2017_region)
