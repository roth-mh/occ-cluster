########
# file to construct
# SDM from OR2020
# observations
########

sampling_region <- readRDS("~/Documents/Oregon State/Research/eBird/occ and grouping checklists/TassledCapOR/sampling_region")
us <- raster::getData('GADM', country = 'US', level = 1)
oregon <- us[us$NAME_1 == "Oregon",]

# process to obtain this:
#   1. load ebd_weta_breeding_or_zf.csv
#   2. load ebd_weta_breeding_or_zf_no2020users.csv
#   3. find (1) - (2)
#   4. find all individual contributors in (3)
OR.2020.obs <- c("obs458115", "obs503692", "obs416347", "obs29674", "obs285064", "obs131077", "obs299653")


z <- WETA_2017[WETA_2017$observer_id %in% OR.2020.obs, c("latitude", "longitude", "observer_id", "duration_minutes")]

######
or2020.maybe <- read.delim("~/Desktop/oregon2020_full_ebirdformat_zerofilled_082720.csv", header = T, sep = ",")
or2020.maybe.WETA <- or2020.maybe[or2020.maybe$Common_Name.x == "Western Tanager",]
or2020.maybe.WETA.2017 <- or2020.maybe.WETA[or2020.maybe.WETA$Year == "2017",]
# or2020.2017 <- or2020[or2020$Year == '2017',]
# or2020.2017.lat.long <- or2020.2017[or2020.2017$Latitude %in% z$latitude & or2020.2017$Longitude %in% z$longitude,]
######


ebd <- read.delim("~/Desktop/ebd_weta_breeding_or_zf.csv", header = T, sep = ",")
ebd.no.or2020 <- read.delim("~/Desktop/ebd_weta_breeding_or_zf_no2020users.csv", header = T, sep = ",")



ebd.2017 <- ebd[ebd$year=="2017",]
ebd.no.or2020.2017 <- ebd.no.or2020[ebd.no.or2020$year=="2017",]


ebd.2017_region <- subset(ebd.2017, ebd.2017$latitude <= 44.5)
ebd.2017_region <- subset(ebd.2017_region, ebd.2017_region$longitude <= -123)

ebd.no.or2020.2017_region <- subset(ebd.no.or2020.2017, ebd.no.or2020.2017$latitude <= 44.5)
ebd.no.or2020.2017_region <- subset(ebd.no.or2020.2017_region, ebd.no.or2020.2017_region$longitude <= -123)

or2020.2017_region <- ebd.2017_region[!(ebd.2017_region$checklist_id %in% ebd.no.or2020.2017_region$checklist_id),]


plot(sampling_region)
points(ebd.no.or2020.2017_region$longitude, ebd.no.or2020.2017_region$latitude)
points(or2020.2017_region$longitude, or2020.2017_region$latitude, col="red")



# write.csv(ebd.no.or2020.2017_region, file = "Documents/Oregon State/Research/eBird/occ and grouping checklists/checklist data/no_OR2020_observations_2017_region.csv")
# write.csv(or2020.2017_region, file = "Documents/Oregon State/Research/eBird/occ and grouping checklists/checklist data/OR2020_observations_2017_region.csv")



nrow(or2020.2017_region)

