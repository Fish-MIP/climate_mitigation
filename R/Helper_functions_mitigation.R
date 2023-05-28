
# CN 24/05/2023 

# function to be applied
extract_global_outputs<-function(netcdf, file = "new"){
  
  # # trial
  # netcdf = "ecotroph_gfdl-esm4_nobasd_historical_nat_default_tcblog10_global_annual_1950_2014.nc"
  # netcdf = "apecosm_ipsl-cm6a-lr_nobasd_historical_nat_default_tcblog10_global_monthly_1850_2014.nc"
  # netcdf = "macroecological_ipsl-cm6a-lr_nobasd_ssp585_nat_default_tcblog10_global_annual_2015_2100.nc"
  # netcdf = "boats_gfdl-esm4_nobasd_historical_nat_default_tcblog10_global_monthly_1950_2014.nc"
  # file = "new"
  
  if(file.exists(file.path(dir, netcdf))){
    
    ######### extract info from netcdf name and print warnings ----
    model = sub("\\_.*", "", netcdf)
    
    if(str_detect(netcdf, "gfdl", negate = FALSE)){
      esm = "gfdl-esm4"
    }else if (str_detect(netcdf, "ipsl", negate = FALSE)){
      esm = "ipsl-cm6a-lr"
    }
    
    # WARNING - add in EC? 
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
    t <- as.character(nc.get.time.series(nc_data))
    
    # this is only to FIX zoom size bins names 
    if(model != "zoomss" & file == "new"){
      bins<-ncvar_get(nc_data, "bins")
    }else if (model != "zoomss" & file == "old"){ # this is only to check DBPM old files (not in DKRZ)
      bins<-ncvar_get(nc_data, "size")
    }else if (model == "zoomss"){
      bins<-c(1:6)}
    
    t_units<-ncatt_get(nc_data, "time", "units")$value
    b_units<-ncatt_get(nc_data, "tcblog10", "units")$value
    
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
    if(b_units != "g m-2"){
      warning(paste(model, esm, scenario, "incorrect biomass units", sep = " "), immediate. = TRUE)
    }
    
    # extract data as raster object: 
    brick_data<-list()
    
    for (i in 1:length(bins)){
      brick_data[[i]]<-brick(file.path(dir, netcdf), level = i)
      print(dim(brick_data[[i]]))
    }
    
    # # check 
    # plot(brick_data[[1]][[1]])
    # mean(brick_data[[1]]) # BOATS size class 1 and 6 is a matrix of NAs
    
    #### IF JB and DT approach - function end here and output is:
    # return(brick_data)
    
    ######### WARINING - remove marginal seas using the land-sea IPSL mask provided by Matthias TO DO ----
    
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
    
    # STEP 2 - calculate annual means and sums of values across grid cells
    # https://gis.stackexchange.com/questions/257090/calculating-and-displaying-mean-annual-precipitation-from-cru-data
    # create vector to serve as index
    
    if(scenario %in% c("historical","picontrol_hist", "picontrol_whole")){
      indices2<-as.Date(indices_subset)
    }else if(scenario %in% c("ssp1", "ssp5", "picontrol_fut")){
      indices2<-as.Date(t)}
    
    indices2<-format(indices2, format = "%Y")
    indices2<-as.numeric(indices2)
    ### WARNING: if fun=sum, results is a raster of 0 instead of NA for e.g. BOATS size class 1 (empty and not 0!)
    brick_data_annual<-lapply(brick_data_subset, FUN = function(x) stackApply(x, indices=indices2, fun=mean))
    
    # # CHECK
    # dim(brick_data_annual[[1]])
    # plot(brick_data_annual[[1]][[15]])
    
    if(scenario == "picontrol_whole"){scenario = "picontrol_fut"}
    
    # STEP 3 - extract global sums 
    
    # create a raster with area cell values (outside loop for efficiency)
    w <- area(brick_data_annual[[1]])*1e6 # WARNING - the code on area weighed means in EC and FishingEffort 09 does /10000, why? these functions are not used but they should be corrected. 
    
    weighted_mean_lat_ls<-list()
    
    for(i in 1: length(brick_data_annual)){ # for each size bin
      
      # i =1
      # code from below on weighted means and adjusted
      # https://stackoverflow.com/questions/55230510/calculating-weighted-spatial-global-annual-averages-across-grid-cells-using-netc
      
      # multiply areas with values
      x <- brick_data_annual[[i]] * w
      # plot(x)
      # plot(w)
      
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
      summarise(weighted_mean_allBio = sum(weighted_mean)) %>%
      ungroup()
    
    #### rename objects with more approriate names:
    brick_data_annual<-brick_data_annual
    SumsBySizeBins_df<-weighted_mean_lat_df
    SumsAllBio_df<-weighted_mean_allBio
    
    return(list(brick_data_annual = brick_data_annual,
                SumsBySizeBins_df = SumsBySizeBins_df,
                SumsAllBio_df = SumsAllBio_df))
    
    rm(brick_data, brick_data_subset, indices, indices2, w, x, weighted_sum, 
       brick_data_annual, SumsBySizeBins_df, weighted_mean_lat_df, SumsAllBio_df, weighted_mean_allBio)

  } # end of if file exists 
  
} # end of function  