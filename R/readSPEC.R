
demon_readSPEC <- function(path, metadata, state, format) { 
  #if ( ! ("spectrolab" %in% installed.packages())) {install.packages("spectrolab", dependencies = T)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)}

  
  # Load required packages
  require(spectrolab)
  require(dplyr)
  require(tidyr)
  
  # Filter information form the metadata to be joined with the spectra data
  meta_state <- metadata %>% 
    filter(COUNTRY_STATE_CODE == state) %>% 
    select(SPECIES, DATE, YEAR, LEAF.SPEC.ID, 
           COUNTRY_STATE, SPCODE, fulcrum_id, LEAF_REPLICATE)
  
  # Read spectra data
  spec <- read_spectra(
    path = path, 
    format = format, 
    #recursive = TRUE, 
    exclude_if_matches = c("BAD", "WR"))
  
  # Transform spectra to dataframe
  tm <- as.data.frame(spec, fix_names = "none", metadata = TRUE)
  # Join metadata with spectra
  spec_data <- full_join(meta_state, tm, by = c("LEAF.SPEC.ID" = "sample_name"))
  # Transform to spectra format
  spect <- as_spectra(spec_data, name_idx = 1, meta_idxs = c(2:8))
  # Transform to dataframe format
  spect_df <- as.data.frame(spect, fix_names = "none", metadata = TRUE)
  # Save both objects 
  obj <- list(spectra = spect, spectra_df = spect_df)
  
  return(obj)
}

#state <- "OH"
#ruta = "/Users/jesusnpl/Dropbox/Oak.Project.SVC data/*Ohio/" 
#metaDT <- read.csv("/Users/jesusnpl/Dropbox/Oak.Project.SVC data/2019.META.DATA/2019.SVC.DATA.LABELS.csv")

#bbb <- demon_read_by_state(path = ruta, 
 #                          metadata = metaDT, 
  #                         state = state, 
   #                        format = "sig")

#bbb[[1]]

#bbb[[2]][1:10, 1:10]

#x <- bbb[[2]]
