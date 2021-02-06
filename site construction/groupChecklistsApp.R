library(leaflet)
library(shiny)
library(geosphere)
library(sf)

# setwd("/Users/MarkRoth/Documents/Oregon State/Research/")
source("helper/helpers.R")
f_name <- "../../../clusteredSites_2020-12-26_.csv"
clusted_sites <- read.delim(f_name, sep=",")

# 
# f_in_WETA <- "eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_2017.csv"
# WETA_2017 <- read.delim(f_in_WETA, header=TRUE, sep = ",")
# WETA_2017_region <- subset(WETA_2017, WETA_2017$latitude <= 44.5)
# WETA_2017_region <- subset(WETA_2017_region, WETA_2017_region$longitude <= -123)
# 
f_in_syn_spec_form <- "../../../eBird/Class Imbalance/generate syn spec/data/linear/syn_species_1_formula.txt"
covObj <- loadCovariates(f_in_syn_spec_form)
# MIN_OBS <- 1
# MAX_OBS <- 100000
# #
# # plot(x = WETA_2017$longitude, y = WETA_2017$latitude, main = "2017 eBird Checklists of Oregon")
# # points(x=WETA_2017_region$longitude, y=WETA_2017_region$latitude, col="blue")
# # legend("bottomright", legend=c("all checklists", "checklists used in calculations"), col=c("black", "blue"), pch=1)
# 
# unique_location_checklists <- sqldf("SELECT * from WETA_2017_region group by latitude, longitude")
# 
# unique_location_checklists <- filter_repeat_visits(
#   unique_location_checklists,
#   min_obs = MIN_OBS,
#   max_obs = MAX_OBS,
#   annual_closure = TRUE,
#   date_var = "formatted_date",
#   site_vars = c("locality_id")
# )
# 
# # x <- unique_location_checklists[1:5,]
# 
# # unique_location_checklists$site_std <- 0
# unique_location_checklists$site_area <- 0
# unique_location_checklists$site_max_dist <- 0

# unique_location_checklists <- clusted_sites
# unique_location_checklists$site <- as.character(unique_location_checklists$site)


color_li <- colorNumeric(palette = c("orange"), domain=seq(length(unique(unique_location_checklists$site))))
# cf <- colorFactor(topo.colors(length(clusted_sites$site) - length(unique(clusted_sites$site))), clusted_sites$site)

c_func <- function(sites){
  t.f <- sites %in% seq(1:100000)
  cols <- ifelse(t.f == TRUE, "orange", "blue")
  return(cols)
}

t.f <- clusted_sites$site %in% seq(1:100000)
max_site <- max(as.double(as.character(clusted_sites[t.f,]$site)))


# RUN APP
runMap <- function(){
  r_colors <- rgb(t(col2rgb(colors()) / 255))
  names(r_colors) <- colors()

  ui <- fluidPage(
    leafletOutput("mymap"),
    p(),
    actionButton("createSite", "Create Site"),
    actionButton("clearCluster", "Clear Cluster"),
    actionButton("endSession", "Return Checklists"),
    verbatimTextOutput('summary')
  )

  server <- function(input, output, session) {
    x <- reactiveValues()
    y <- reactiveValues()
    pvsChecklist <- reactiveValues()
    pvsChecklist$id <- 0
    
    x$site <- max_site + 1
    y$cluster <- list()
    df_var <- reactiveValues(df=unique_location_checklists)

    output$mymap <- renderLeaflet({
      leaflet() %>% addTiles() %>%
        addCircleMarkers(
          fillColor = c_func(unique_location_checklists$site),
          fillOpacity = .7,
          lng=unique_location_checklists$longitude,
          lat=unique_location_checklists$latitude,
          # data=unique_location_checklists$checklist_id,
          popup = paste0(
            "<b>summer_nbr_TCB_mean_2400: </b>",
            round(unique_location_checklists$summer_nbr_TCB_mean_2400,4),
            "<br>",
            "<b>spring_nbr_TCG_mean_75: </b>",
            round(unique_location_checklists$spring_nbr_TCG_mean_75,4),
            "<br>",
            "<b>fall_nbr_B2_stdDev_300: </b>",
            round(unique_location_checklists$fall_nbr_B2_stdDev_300,4),
            "<br>",
            "<b>slope_mean_300: </b>",
            round(unique_location_checklists$slope_mean_300,4),
            "<br>",
            "<b>summer_b5_B4_stdDev_1200: </b>",
            round(unique_location_checklists$summer_b5_B4_stdDev_1200,4),
            "<br>"),
          options = popupOptions(closeOnClick = FALSE)
        ) %>% addScaleBar(position = "bottomright") %>%
        addProviderTiles('Esri.WorldImagery')
    })
    
    observe({
      click<-input$mymap_marker_click
      if(!is.null(click)){
        checklist <- unique_location_checklists[
          click$lat == unique_location_checklists$latitude & click$lng == unique_location_checklists$longitude,
          ]
        
        c_id <- as.character(checklist$checklist_id)
        if(!(c_id %in% y$cluster) & c_id != pvsChecklist$id){
          y$cluster <- append(y$cluster, c_id)
          pvsChecklist$id <- c_id
          
          
          leafletProxy("mymap") %>%
            addCircleMarkers(
              lng=checklist$longitude,
              lat=checklist$latitude,
              fillColor = color_li(x$site),
              fillOpacity = .4,
              # layerId = x$site,
              # data=unique_location_checklists$checklist_id,
              popup = paste0(
                "<b>summer_nbr_TCB_mean_2400: </b>",
                round(unique_location_checklists$summer_nbr_TCB_mean_2400,4),
                "<br>",
                "<b>spring_nbr_TCG_mean_75: </b>",
                round(unique_location_checklists$spring_nbr_TCG_mean_75,4),
                "<br>",
                "<b>fall_nbr_B2_stdDev_300: </b>",
                round(unique_location_checklists$fall_nbr_B2_stdDev_300,4),
                "<br>",
                "<b>slope_mean_300: </b>",
                round(unique_location_checklists$slope_mean_300,4),
                "<br>",
                "<b>summer_b5_B4_stdDev_1200: </b>",
                round(unique_location_checklists$summer_b5_B4_stdDev_1200,4),
                "<br>"),
              options = popupOptions(closeOnClick = FALSE)
            ) %>% addScaleBar(position = "bottomright") %>%
            addProviderTiles('Esri.WorldImagery')
        }
      }
    })
    
    observeEvent(input$endSession, {
      stopApp(df_var$df)
    })
    
    observeEvent(input$clearCluster, {
      y$cluster <- list()
    })
    
    observeEvent(input$createSite, {
      
      if(length(y$cluster) == 0){
        return()
      }

      for(c_id in y$cluster){
        df_var$df[df_var$df$checklist_id == c_id,]$site <- as.character(x$site)
        # disp(df_var$df[df_var$df$checklist_id == c_id,]$site)
      }

      set_of_pts <- df_var$df[df_var$df$site == x$site,]
      # disp(x$site)
      ch_open <- chull(set_of_pts[c('latitude', 'longitude')])
      ch <- c(ch_open, ch_open[1])
      edge_pts <- set_of_pts[ch,]

      coords <- edge_pts[c('latitude', 'longitude')]
      coordinates(coords) <- ~longitude+latitude
      proj4string(coords)=CRS("+init=epsg:4326")
      coords_trans <- spTransform(coords, CRS("+proj=aea +lat_1=34 +lat_2=47 +lat_0=43 +lon_0=-120 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
      area <- polyarea(coords_trans@coords[,1], coords_trans@coords[,2])
      site_std <- list()
      for(cov in covObj$siteCovs){
        site_std <- append(site_std, sd(df_var$df[df_var$df$site == x$site,][[cov]]))
      }

      max_km <- poly_length(coords_trans@coords[,1], coords_trans@coords[,2])/1000

      df_var$df[df_var$df$site == x$site,]$site_area <- area
      df_var$df[df_var$df$site == x$site,]$site_max_dist <- max_km

      leafletProxy("mymap") %>%
        addPolygons(
          lng=edge_pts$longitude,
          lat=edge_pts$latitude,
          fillColor = color_li(x$site),
          fillOpacity = .75
        ) %>%
        addCircleMarkers(
          lng=set_of_pts$longitude,
          lat=set_of_pts$latitude,
          popup = paste0(
            "<b>SITE AREA (km^2): </b>",
            round(sqrt(area/(1000*1000)),3),
            "<br>",
            "<b>Max Dist Between 2 pts (km): </b>",
            round(max_km,2),
            "<br>",
            "<b>SD of summer_nbr_TCB_mean_2400: </b>",
            round(site_std[[1]],3),
            "<br>",
            "<b>SD spring_nbr_TCG_mean_75: </b>",
            round(site_std[[2]],3),
            "<br>",
            "<b>SD fall_nbr_B2_stdDev_300: </b>",
            round(site_std[[3]],3),
            "<br>",
            "<b>SD slope_mean_300: </b>",
            round(site_std[[4]],3),
            "<br>",
            "<b>SD summer_b5_B4_stdDev_1200: </b>",
            round(site_std[[5]],3),
            "<br>")
        )

      y$cluster <- list()
      x$site <- x$site + 1
      
    })
    #####
    
    output$summary <- renderPrint({
      sprintf("current cluster %s", paste0(y$cluster, collapse=", "))
    })

    # output$summary <- renderPrint({
    #   sprintf("site #: %s with %s checklists and area: %s km^2", z()$site, length(z()$c_ids), round(sqrt(z()$area/(1000*1000)),3))
    #   # sprintf("cluster with following ids: %s", paste0(x$cluster, collapse = ", "))
    # })
  }
  shinyApp(ui, server)
}
df_pls <- runApp(runMap())
# write.csv(df_pls, file=paste("clusteredSites", Sys.Date(), ".csv", sep="_"))





#####
# visualize the sites
#####
sites_to_viz <- DBSC.sites

l <- leaflet() %>% addTiles() %>%
  addCircleMarkers(
    fillColor = c_func(sites_to_viz$site),
    fillOpacity = .7,
    lng=sites_to_viz$longitude,
    lat=sites_to_viz$latitude,
    # data=unique_location_checklists$checklist_id,
    popup = paste0(as.character(sites_to_viz$site)),
    options = popupOptions(closeOnClick = FALSE)
  ) %>% addScaleBar(position = "bottomright") %>%
  addProviderTiles('Esri.WorldImagery') 

for(s in unique(sites_to_viz$site)){
  set_of_pts <- sites_to_viz[sites_to_viz$site == s,]
  # disp(x$site)
  ch_open <- chull(set_of_pts[c('latitude', 'longitude')])
  ch <- c(ch_open, ch_open[1])
  edge_pts <- set_of_pts[ch,]
  l <- l %>% addPolygons(lng=edge_pts$longitude, lat=edge_pts$latitude, fillColor = 'red', fillOpacity = .75)
}
l

