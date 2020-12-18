library(ggmap)
library(ggsn)
library(ggspatial)


# WETA_2017_narrow <- subset(WETA_2017, WETA_2017$latitude <= 42.4)
# WETA_2017_narrow <- subset(WETA_2017_narrow, WETA_2017_narrow$longitude >= -119)
# WETA_2017_narrow <- subset(WETA_2017_narrow, WETA_2017_narrow$longitude <= -118.2)

p <- ggmap(get_googlemap(center = c(lon = -123.700000, lat = 43.4),
                         zoom = 9, scale = 2,
                         maptype ='terrain',
                         color = 'color'))

p


get_googlemap(
  center = c(lon = -95.3632715, lat = 29.7632836), zoom = 10,
  size = c(640, 640), scale = 2, 
  maptype = "terrain"
  )

register_google(key="AIzaSyDCUBgYeqcfEhQvpDLF1A9fxXqilVv6-hw")

p <- ggmap(get_googlemap(center = c(lon = -121.300000, lat = 42.150000), zoom = 11))
p_sate <- ggmap(get_googlemap(center = c(lon = -118.535000, lat = 42.220732), zoom = 11, maptype = "satellite"))
p_road <- ggmap(get_googlemap(center = c(lon = -118.535000, lat = 42.220732), zoom = 11, maptype = "roadmap"))
p_hybr <- ggmap(get_googlemap(center = c(lon = -118.535000, lat = 42.220732), zoom = 11, maptype = "hybrid"))
p_sate
p_road
p_hybr

WETA_2017_narrow$lon <- WETA_2017_narrow$longitude
WETA_2017_narrow$lat <- WETA_2017_narrow$latitude

wonk$lon <- wonk$longitude
wonk$lat <- wonk$latitude
# 
clust_wonk$lon <- clust_wonk$longitude
clust_wonk$lat <- clust_wonk$latitude

# round$lon <- round$longitude
# round$lat <- round$latitude


x_WETA <- WETA_2017_narrow[c("lat", "lon")]
x_wonk <- wonk[c("lat", "lon")]
x_clust_wonk <- clust_wonk[c("lat", "lon")]
# x_round <- round[c("lat", "lon")]


x_WETA <- x_WETA[x_WETA$lat >= 42,]
x_WETA <- x_WETA[x_WETA$lat <= 42.3,]

x_WETA <- x_WETA[x_WETA$lon >= -121.5,]
x_WETA <- x_WETA[x_WETA$lon <= -121.09,]

x_WETA$type <- "weta"
x_WETA$long <- x_WETA$lon

x_clust_wonk$type <- "clustGeo"
x_clust_wonk$long <- x_clust_wonk$lon

x_wonk$type <- "SKATER"
x_wonk$long <- x_wonk$lon

# x_clust_wonk$type <- "clustGeo"
# x_round$type <- "round"


df <- rbind(x_WETA, x_wonk, x_clust_wonk)#, x_round) 
df <- subset(df, select = -c(lon))
plot.new()
# p_hybr + geom_point(data = df, aes(x = lat, y = lon, size = 3, colour = c("black", "red", "blue", "yellow")))
p + geom_point(data = df, size = 3,
               aes(x=long, y=lat, color = type)) +
  ggsn::scalebar(data=df, location="topright", dist = 5, model = 'WGS84', 
                   transform = T, dist_unit = "km", st.bottom = T, st.size = 3, border.size = .4)
                   # x.min = -121.5194, y.max = 42.31244, x.max= -121.0799, y.min = 41.98663)

p$coordinates$transform("WGS84")
p$data


p + ggsn::scalebar(location="bottomleft", dist = 5, model = 'WGS84', 
                   transform = T, dist_unit = "km", st.bottom = F, st.size = 3,
                   x.min = -121.5194, y.max = 42.31244, x.max= -121.0799, y.min = 41.98663)
# ggsn::scalebar(df, location="bottomleft", dist = 1000, dist_unit = "km", 
#                st.size=5, model = 'WGS84', transform = T, st.color="red",
#                st.bottom = F)


p +
  geom_point(aes(long, lat), data=x_WETA, color="black", size = 3) +
  ggsn::scalebar(x_WETA, location="bottomleft", dist = 1000, dist_unit = "km", 
                 st.size=5, height=.02, model = 'WGS84', transform = T, st.color="red", 
                 x.min = -118.9, x.max=-118.2, y.min = 42, y.max=42.4)




p + geom_point(data = x_WETA, size = 3) + geom_point(data = x_wonk, colour="red", size = 3) + geom_point(data = x_clust_wonk, colour="blue", size = 3) + geom_point(data = x_round, colour="yellow", size = 3)
legend("bottomright", legend=c("SKATER", "ClustGeo", "rounded"), col=c("red", "blue", "yellow"), pch=1)

plot.new()
plot(p+plot(x=WETA_2017_narrow$longitude, y=WETA_2017_narrow$latitude))
points(x=wonk$longitude, y=wonk$latitude, col="red")
points(x=round$longitude, y=round$latitude, col="green")
points(x=clust_wonk$longitude, y=clust_wonk$latitude, col="blue")


























##### centering clustGeo

p_clust <- ggmap(get_googlemap(center = c(lon = -120.742100, lat = 42.830000), zoom = 11 ))

W <- WETA_2017_region

x_WETA <- W[W$longitude >= min(p_clust$data$lon),]
x_WETA <- x_WETA[x_WETA$longitude <= max(p_clust$data$lon),]

x_WETA <- x_WETA[x_WETA$latitude >= min(p_clust$data$lat),]
x_WETA <- x_WETA[x_WETA$latitude <= max(p_clust$data$lat),]

wonk <- wonk[wonk$long >= min(p_clust$data$lon),]
wonk <- wonk[wonk$long <= max(p_clust$data$lon),]

wonk <- wonk[wonk$lat >= min(p_clust$data$lat),]
wonk <- wonk[wonk$lat <= max(p_clust$data$lat),]




clust_wonk$long <- clust_wonk$longitude
clust_wonk$lat <- clust_wonk$latitude
clust_wonk$algorithm <- "clustGeo site"

x_WETA$long <- x_WETA$longitude
x_WETA$lat <- x_WETA$latitude
x_WETA$algorithm <- "all checklists"

wonk$long <- wonk$longitude
wonk$lat <- wonk$latitude
wonk$algorithm <- "SKATER site"

x_WETA <- subset(x_WETA, select = c(covObj$siteCovs, "long", "lat", "algorithm"))
mean_att_df <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(mean_att_df) <- covObj$siteCovs
for(cov_att in covObj$siteCovs){
  cov_val <- mean(x_WETA[[cov_att]])
  mean_att_df[[cov_att]] <- cov_val
}

x_WETA$dist <- apply(x_WETA, 1, function(x) dist(c(x[covObj$siteCovs], mean_att_df))[1])
wonk$dist <- c(1.2464303, 1.9173307)
clust_wonk$dist <- c(0.8079995, 0.8191715)


x_WETA <- subset(x_WETA, select = c(long, lat, algorithm, dist))
clust_wonk <- subset(clust_wonk, select = c(long, lat, algorithm, dist))
wonk <- subset(wonk, select = c(long, lat, algorithm, dist))
# x_WETA <- subset(x_WETA, select = c(long, lat, algorithm))
# clust_wonk <- subset(clust_wonk, select = c(long, lat, algorithm))
# wonk <- subset(wonk, select = c(long, lat, algorithm))
wonk <- wonk[wonk$lat < 43,]

p_clust + geom_point(data = rbind(x_WETA,clust_wonk, wonk), size = 3, aes(x=long, y=lat, color=algorithm)) +
  ggsn::scalebar(data=rbind(x_WETA,clust_wonk, wonk), location="topright", dist = 2.5, model = 'WGS84',
                 transform = T, dist_unit = "km", st.bottom = T, st.size = 3, border.size = .4,
                 anchor = c(
                   x = -120.570000, 
                   y = 42.950000
                 ))



# idea: fill based on distance from some average midpoint
p_clust + geom_point(data = rbind(x_WETA,clust_wonk, wonk), size = c(7,7,7,7,7,7,7,7,7,7,7,7,7, 4.5, 4.5, 4.5, 4.5), aes(x=long, y=lat, color=algorithm, shape=algorithm, fill =algorithm)) +
  ggsn::scalebar(data=rbind(x_WETA,clust_wonk, wonk), dist = 2.5, model = 'WGS84', st.dist = .028,
                 transform = T, dist_unit = "km", st.bottom = T, st.size = 5, border.size = .2,
                 anchor = c(
                   x = -120.865000, 
                   y = 42.965000
                 )) + scale_shape_manual(values=c(16,17,25)) + 
  scale_fill_manual(values=c("#000000", "#ffa31a", "#ff00ff")) +
  guides(shape = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=13)) +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15,face="bold")) + 
  scale_color_manual(values=c("#000000", "#ffa31a", "#ff00ff"))

#













# p_clust + geom_point(data = rbind(x_WETA,clust_wonk, wonk), size = 4, aes(x=long, y=lat, color=algorithm, shape=algorithm, fill = dist)) +
#   ggsn::scalebar(data=rbind(x_WETA,clust_wonk, wonk), location="topright", dist = 2.5, model = 'WGS84',
#                  transform = T, dist_unit = "km", st.bottom = F, st.size = 3, border.size = .4,
#                  anchor = c(
#                    x = -120.865000, 
#                    y = 42.965000
#                  )) + scale_shape_manual(values=c(21,22,23)) +
#   scale_color_manual(values=c("#ffff00", "#9900ff", "#ff0000"))
# 
# 
# 



