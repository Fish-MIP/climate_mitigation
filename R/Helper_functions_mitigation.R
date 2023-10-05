
# CN 24/05/2023 

# function to be applied
extract_global_outputs<-function(netcdf, file = "new"){
  
  # # trial
  # a<-combinations %>% filter(identifier == "boats_ipsl_historical")
  # netcdf = a$netcdf_name
  # file = "new"
   
  if(file.exists(file.path(dir, netcdf))){
    
    ######### extract info from netcdf name and print warnings ----
    model = sub("\\_.*", "", netcdf)
    
    if(str_detect(netcdf, "gfdl", negate = FALSE)){
      esm = "gfdl-esm4"
    }else if (str_detect(netcdf, "ipsl", negate = FALSE)){
      esm = "ipsl-cm6a-lr"
    }
    
    if(str_detect(netcdf, "monthly", negate = FALSE)){
      time_step = "monthly"
    }else if (str_detect(netcdf, "annual", negate = FALSE)){
      time_step = "annual"
    }
    
    if(str_detect(netcdf, "historical", negate = FALSE)){
      scenario = "historical"
    }else if (str_detect(netcdf, "ssp126", negate = FALSE)){
      scenario = "ssp1"
    }else if (str_detect(netcdf, "ssp585", negate = FALSE)){
      scenario = "ssp5"
    }else if (str_detect(netcdf, "picontrol|2100", negate = FALSE)) {
      scenario = "picontrol_fut"
    } else if (str_detect(netcdf, "picontrol|2014", negate = FALSE)) {
      scenario = "picontrol_hist"}
    
    # extract info from netcdf description: 
    nc_data <- nc_open(file.path(dir, netcdf))
    
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat", verbose = F)
    
    # check function below re time vector
    t_units<-ncatt_get(nc_data, "time", "units")$value
    
    # get time vector based on file name
    if(str_detect(netcdf, "1850", negate = FALSE)){
      stTime = "1850-1-1"
      enTime = "2014-12-31"
    }else if (str_detect(netcdf, "1950", negate = FALSE)){
      stTime = "1950-1-1"
      enTime = "2014-12-31"
    }else if (str_detect(netcdf, "2015", negate = FALSE)){
      stTime = "2015-1-1"
      enTime = "2100-12-31"
    }
    
    if(time_step == "monthly"){
      time_step_vector = "month"
    }else if(time_step == "annual"){
      time_step_vector = "year"
    }
    
    library(lubridate) # https://data.library.virginia.edu/working-with-dates-and-time-in-r-using-the-lubridate-package/
    t1<- as.character(seq(ymd(stTime), ymd(enTime), by = time_step_vector))
    print(paste("Start time from file name ",t1[1], sep = ""))
    print(paste("End time from file name ",t1[length(t1)], sep = ""))
    
    # get time vector based on built in function 
    t <- as.character(nc.get.time.series(nc_data))
    print(paste("Start time with built in function ",t[1], sep = ""))
    print(paste("End time with built in function ",t[length(t)], sep = ""))
    
    if((t1[1] != t[1]) | (t1[length(t1)] != t[length(t)])){
      warning(paste(model, esm, scenario, "incorrect time vector", sep = " "), immediate. = TRUE)
      ## trust the vector from file name (this function does not work with inputs)
      t<-t1
    }
    
    # this is only to FIX zoom size bins names 
    if(model != "zoomss" & file == "new"){
      bins<-ncvar_get(nc_data, "bins")
    }else if (model != "zoomss" & file == "old"){ # this is only to check DBPM old files (not in DKRZ)
      bins<-ncvar_get(nc_data, "size")
    }else if (model == "zoomss"){
      bins<-c(1:6)}
    
    t_units<-ncatt_get(nc_data, "time", "units")$value
    b_units<-ncatt_get(nc_data, "tcblog10", "units")$value
    
    missingValues<-ncatt_get(nc_data, "tcblog10", "missing_value")$value
    
    nc_close(nc_data)
    
    # print warnings 
    
    stLon<-lon[1]
    enLon<-lon[length(lon)]
    stLat<-lat[1]
    enLat<-lat[length(lat)]
    stTime<-t[1]
    enTime<-t[length(t)]
    
    if(stLon != -179.5){
      warning(paste(model, esm, scenario, "incorrect starting Lon", sep = " "), immediate. = TRUE)
    }
    if(enLon != 179.5){
      warning(paste(model, esm, scenario, "incorrect ending Lon", sep = " "), immediate. = TRUE)
    }
    ## dbpm and zoom IPSL incorrect ln fixed below 
    
    if(stLat != 89.5){
      warning(paste(model, esm, scenario, "incorrect starting Lat", sep = " "), immediate. = TRUE)
    }
    if(enLat != -89.5){
      warning(paste(model, esm, scenario, "incorrect ending Lat", sep = " "), immediate. = TRUE)
    }
    if(scenario == "historical" & !stTime %in% c("1950-01-01","1850-01-01")){ # some model include 100 years more 
      warning(paste(model, esm, scenario, "incorrect starting time", sep = " "), immediate. = TRUE)
    }
    if(scenario != "historical" & !stTime %in% c("2015-01-01")){
      warning(paste(model, esm, scenario, "incorrect starting time", sep = " "), immediate. = TRUE)
    }
    if(scenario == "historical" & !enTime %in% c("2014-12-01", "2014-01-01")){ # models can be monthly or annual 
      warning(paste(model, esm, scenario, "incorrect ending time", sep = " "), immediate. = TRUE)
    }
    if(scenario != "historical" & !enTime %in% c("2100-12-01", "2100-01-01")){ # models can be monthly or annual 
      warning(paste(model, esm, scenario, "incorrect ending time", sep = " "), immediate. = TRUE)
    }
    if(t_units != "days since 1601-1-1 00:00:00"){
      warning(paste(model, esm, scenario, "incorrect time units", sep = " "), immediate. = TRUE)
    }
    if(bins[1] != 1 & bins[6] != 6){
      warning(paste(model, esm, scenario, "incorrect bins names", sep = " "), immediate. = TRUE)
    }
    if(bins[length(bins)] != 6){
      warning(paste(model, esm, scenario, "incorrect bins dimension", sep = " "), immediate. = TRUE)
    }
    ## macroecological 5 bins and boats 2-4 bins fixed below 
    
    if(b_units != "g m-2"){
      warning(paste(model, esm, scenario, "incorrect biomass units", sep = " "), immediate. = TRUE)
    }
    
    # extract data as raster object: 
    brick_data<-list()
    
    for (i in 1:length(bins)){
      brick_data[[i]]<-brick(file.path(dir, netcdf), level = i) # level: integer > 0 (default=1). To select the 'level' (4th dimension variable) to use, if the file has 4 dimensions, e.g. to create a RasterBrick of weather over time at a certain height.
      print(dim(brick_data[[i]]))
    }
    
    # # check 
    # plot(brick_data[[1]][[1]])
    # mean(brick_data[[1]]) # BOATS size class 1 and 6 is a matrix of NAs
    
    # check if land is specified as missingValues or as NA
    # options(scipen=999)
    trial<-brick_data[[1]][[1]] # take first layer and first year
    if(is.null(trial[trial == missingValues]) == FALSE){ # if there are missingValues cells 
      warning(paste(model, esm, scenario, "missing values not as NAs", sep = " "), immediate. = TRUE)
      for (i in 1:length(brick_data)){
        brick_data[[i]][brick_data[[i]] == missingValues]<-NA
      }
    }
    # # back check 
    # plot(trial)
    # trial[is.na(trial)] <- missingValues
    # plot(trial)
    
    #### IF JB and DT approach - function end here and output is:
    # return(brick_data)
    
    ######### remove marginal seas using the land-sea IPSL mask provided by Matthias
    
    if(esm == "ipsl-cm6a-lr" & model %in% c("dbpm", "macroecological", "feisty", "zoomss")){
      
      # read in the mask
      dir_land<-"/rd/gem/private/users/camillan"
      brick_land<-brick(file.path(dir_land, "IPSL-CM6A-LR_lsm_nolakes.nc"))
      
      # CHECK 
      # plot(brick_land)
      # plot(brick_data[[5]][[dim(brick_data[[2]])[3]]])
      
      # # check extent
      # extent(brick_land)
      # extent(brick_data[[1]])
      
      # deal with different extent for Zoom and DBPM IPSL - from SO code
      if(stLon != -179.5){
        bb <- extent(-180, 180, -90, 90)
        for(i in 1:6){
          brick_data[[i]]<-setExtent(brick_data[[i]], bb, keepres=FALSE)
        }
      }
      
      # multiply data X mask with na.rm = F to delete lakes
      brick_data[[1]]<-brick_data[[1]]*brick_land
      # plot(brick_data[[1]][[1]])
      # plot(trial[[1]])
      brick_data[[2]]<-brick_data[[2]]*brick_land
      brick_data[[3]]<-brick_data[[3]]*brick_land
      brick_data[[4]]<-brick_data[[4]]*brick_land
      brick_data[[5]]<-brick_data[[5]]*brick_land
      
      if(model != "macroecological"){
        brick_data[[6]]<-brick_data[[6]]*brick_land
      }
      
    }
    
    ######### calculate total annual sums (OR weighted annual means) ----
    
    # STEP 1 - remove 1850-1950 as not all models have them
    indices<-t
    
    if(scenario %in% c("historical","picontrol_hist")){
      indices_subset<-indices[indices>="1950-01-01"]
      indices_position<-match(indices_subset,indices)
      brick_data_subset<-lapply(brick_data, FUN = function(x) raster::subset(x, indices_position))
    }else if (scenario %in% c("ssp1","ssp5","picontrol_fut")){
      brick_data_subset<-brick_data
    } else if(scenario == "picontrol_whole"){
      indices_subset<-indices[indices>="2015-01-01"]
      indices_position<-match(indices_subset,indices)
      brick_data_subset<-lapply(brick_data, FUN = function(x) raster::subset(x, indices_position))
    }
    
    # # CHECK 
    # plot(brick_data_subset[[1]][[1]])
    
    # STEP 2 - calculate annual means
    # https://gis.stackexchange.com/questions/257090/calculating-and-displaying-mean-annual-precipitation-from-cru-data
    # create vector to serve as index
    
    if(scenario %in% c("historical","picontrol_hist", "picontrol_whole")){
      indices2<-as.Date(indices_subset)
    }else if(scenario %in% c("ssp1", "ssp5", "picontrol_fut")){
      indices2<-as.Date(t)}
    
    indices2<-format(indices2, format = "%Y")
    indices2<-as.numeric(indices2)
    tic()
    brick_data_annual<-lapply(brick_data_subset, FUN = function(x) stackApply(x, indices=indices2, fun=mean))
    toc()
    
    # # CHECK
    # dim(brick_data_annual[[1]])
    # plot(brick_data_annual[[1]][[15]])
    
    if(scenario == "picontrol_whole"){scenario = "picontrol_fut"}
    
    # STEP 3 - extract global sums 
    
    # create a raster with area cell values (outside loop for efficiency)
    # crs(brick_data_annual[[1]]) # lon/lat
    # plot(brick_data_annual[[2]][[1]])
    w <- area(brick_data_annual[[1]])*1e6 
    w2<-area(brick_data_annual[[2]][[1]], na.rm = TRUE)*1e6 # size class 2 (available in all models) and year 1
    # plot(w)
    # plot(w2) # same but without land 
    # If x is a Raster* object: RasterLayer or RasterBrick. Cell values represent the size of the cell in km2, 
    # or the relative size if weights=TRUE. 
    # If the CRS is not longitude/latitude the values returned are the product of the cell resolution 
    # (typically in square meter).
     
    weighted_mean_lat_ls<-list()
    
    for(i in 1: length(brick_data_annual)){ # for each size bin
      
      # i =2
      # code from below on weighted means and adjusted
      # https://stackoverflow.com/questions/55230510/calculating-weighted-spatial-global-annual-averages-across-grid-cells-using-netc
      
      # multiply areas with values
      x <- brick_data_annual[[i]] * w2
      
      # CHECK 
      # plot(brick_data_annual[[i]][[1]])
      # plot(x[[i]])
      # plot(w2)
      
      # compute sum
      weighted_mean_lat<-cellStats(x, sum, na.rm = TRUE) 
      weighted_mean_lat_ls[[i]]<-data.frame(Year = unique(indices2), weighted_mean = weighted_mean_lat) %>%
        mutate(
          Year = as.numeric(Year),
          file = netcdf,
          mem = model,
          esm = esm,
          scenario = scenario,
          bin = bins[[i]])
      rownames(weighted_mean_lat_ls[[i]])<-NULL
      
    }
    
    weighted_mean_lat_df<-do.call(rbind, weighted_mean_lat_ls)
    
    ######### sum biomass across common bins ----
    
    weighted_mean_allBio<-weighted_mean_lat_df %>%
      filter(bin %in% c(2:5)) %>% 
      group_by(Year, mem, esm, scenario, file) %>%
      summarise(sum_allBio = sum(weighted_mean)) %>% 
      ungroup()

    ## CHECK     
    # head(weighted_mean_allBio)
    # a<-filter(weighted_mean_lat_df, Year == 1950)
    # sum(a$weighted_mean)
    # filter(weighted_mean_allBio, Year == 1950)
    
    # ### this needs to be done outside the function as future scenarios will not have the ref year
    # ######### calculate delta change ----
    # 
    # # as per protocol 𝑑𝑦(𝑇)=𝑦(𝑇)−𝑦(0)
    # # with ref 1995-2014 inclusive
    # 
    # reference <- weighted_mean_allBio %>% 
    #   filter(Year >= 1995 , Year <= 2014) %>% 
    #   group_by(mem, esm, scenario, file) %>% 
    #   summarise(ref = mean(sum_allBio, na.rm = TRUE))
    # 
    # reference<-reference$ref
    # 
    # weighted_mean_allBio<-weighted_mean_allBio %>% 
    #   mutate(delta = sum_allBio - reference)
    
    #### rename objects:
    brick_data_annual<-brick_data_annual
    SumsBySizeBins_df<-weighted_mean_lat_df
    SumsAllBio_df<-weighted_mean_allBio
    
    return(list(brick_data_annual = brick_data_annual,
                SumsBySizeBins_df = SumsBySizeBins_df, # this is the only data provided as csv
                SumsAllBio_df = SumsAllBio_df))
    
    rm(brick_data, brick_data_subset, indices, indices2, w, x, weighted_sum, 
       brick_data_annual, SumsBySizeBins_df, weighted_mean_lat_df, SumsAllBio_df, weighted_mean_allBio)

  } # end of if file exists 
  
} # end of function  