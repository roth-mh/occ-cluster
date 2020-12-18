library(dplyr)
################################
# NLCD Land Classification Types
# https://www.mrlc.gov/data/legends/national-land-cover-database-2016-nlcd2016-legend
################################

resistanceVals <- function(x){
  sorted_x <- sort(x)
  if(round(x[1]/10) == round(x[2]/10)){ # come from the same 'CATEGORY'
    return(0) #NOT RIGHT!!!!!!! agh
  }
  
  if(round(x[1]/10) == 1){ # come from water
    if(round(x[2]/10) == 3){ # come from developed
     return() 
    }
  }
  
}


############################################
# x == 0 ~  "no value?",
# x == 11 ~ "water",
# x == 12 ~ "snow",
# x == 21 ~ "developed/open space",
# x == 22 ~ "developed, low intensity",
# x == 22 ~ "developed, medium intensity",
# x == 22 ~ "developed, high intensity",
# x == 31 ~ "barren land (rock/sand/clay)",
# x == 41 ~ "deciduous forest",
# x == 42 ~ "evergreen forest",
# x == 43 ~ "mixed forest",
# x == 51 ~ "dwarf shrub",
# x == 52 ~ "shrub/scrub",
# x == 71 ~ "grassland/herbaceous",
# x == 72 ~ "sedge/herbaceous",
# x == 73 ~ "lichens",
# x == 74 ~ "moss",
# x == 81 ~ "pasture/hay",
# x == 82 ~ "cultivated crops",
# x == 90 ~ "woody wetlands",
# x == 95 ~ "emergent herbaceous wetlands"
############################################

# resistanceVals(c(z,95))

colnames(WETA_2017_region)
