library(tidyverse)
library(rtrim)
library(progress)

options(dpyr.summarise.inform=FALSE)
options(readr.show_progress=FALSE)

data_schemes<-read_csv2("data/data_schemes.csv")

spschemesfinal<-read_csv2("data/data-bmp-whichstratum.csv") %>% select(-check,-`Latin Name`,-stcp) %>%
  pivot_longer(2:ncol(.),values_to = "scheme") %>%
  filter(!is.na(scheme) & !scheme==0) %>%
  select(-scheme) %>%
  rename(
    Scheme=name,
    Species_nr=EURING,
  ) %>%
  #mutate(
  #  Changepoints="all"
  #)
  left_join(
    read_csv2("data/data-bmp-whichstratum.csv") %>% select(EURING,stcp) %>% rename(Species_nr=EURING),
    by="Species_nr"
  ) %>%
  rename(Changepoints=stcp)
  

schemespecies <- spschemesfinal %>% select(Species_nr) %>% distinct() %>% unlist() %>% unname()

prefix<-"BMP"
country_code<-6

n_iter <- length(schemespecies)
pb <- progress_bar$new(format = "  töötlen [:bar] :percent eta: :eta", total = n_iter, clear = FALSE, width= 60)
for(i in 1:length(schemespecies)) {
  # i<-4
  # i<-which(schemespecies==11870)
  # i<-which(schemespecies==4330)
  # i<-which(schemespecies==1540)
  # i<-which(schemespecies==110)
  spp<-schemespecies[i]
  spsch<-spschemesfinal %>% filter(Species_nr==spp) %>% select(Scheme) %>% unlist() %>% unname()
  spcp<-spschemesfinal %>% filter(Species_nr==spp) %>% select(Changepoints) %>% unlist() %>% unname() %>% unique()
  countdata<-data_schemes %>%
    filter(year>=2013) %>%
    filter(euring==spp)
  #countdata<-spschemesfinal %>% filter(Species_nr==spp) %>% select(-Changepoints) %>%
  #  rename(scheme=Scheme, euring=Species_nr) %>%
  #  left_join(countdatast, by=c("euring","scheme"))
  
  if(nrow(countdata)==0) {
    next;
  }
  
  sppyears<-range(countdata$year)
  
  countsfile<-paste0(prefix,"_",spp,"_1_",country_code,"_counts.csv")
  
  countdata %>%
    mutate(weights=1) %>%
    select(site,year,count,weights) %>%
    write_excel_csv2(paste0("ST_inputs/",countsfile), na="NA")
  
  # 0-base year prevention
  positivecountyears<-countdata %>%
    group_by(year) %>%
    summarise(observations=n(), positive=length(count[!is.na(count)&count>0]) ) %>%
    slice(which(!is.na(positive) & positive>0)) %>%
    select(year) %>%
    unlist() %>% unname() %>% range()
  
  arginputstratumfile<-paste0(prefix,"_",spp,"_1_",country_code,"_arg_input_stratum.csv")
  stratum<-paste0(prefix,"_",spp,"_1_",country_code)
  
  tibble(
    File=stratum,
    Base_year_first_year=positivecountyears[1],
    Base_year_last_year=positivecountyears[1],
    Changepoints=spcp,
    Serial_correlation=ifelse("PR0065" %in% spsch,TRUE,FALSE),
    Overdispersion=TRUE,
    Presence_weights=TRUE,
    Presence_monthfactors=FALSE,
    #Year_from=ifelse(positivecountyears[1]<=2019,2019,positivecountyears[1]),
    Year_from=positivecountyears[1],
    Save_fitted_values=TRUE
  ) %>%
    write_excel_csv2(paste0("ST_inputs/",arginputstratumfile), na="NA" )
  
  pb$tick()
}


rtrim_folder<-getwd()
working_directory<-getwd()
source("R/rtrim-shell-mod.R")
listSpeciesStratumCombinations <<- dir(paste0(working_directory,"//ST_Inputs"), pattern = "arg_input_stratum.csv")
numberSpeciesStratumCombinations <<- length (listSpeciesStratumCombinations)
overview <<- makeOverview(listSpeciesStratumCombinations)
#rtrim.strata("ST_inputs/")


inputfilespath <- "ST_inputs/"
n_iter <- length(schemespecies)
pb <- progress_bar$new(format = "  töötlen [:bar] :percent eta: :eta", total = n_iter, clear = FALSE, width= 60)
for(j in 1:numberSpeciesStratumCombinations) {
  
  arguments <- read_delim(paste0(inputfilespath,listSpeciesStratumCombinations[j]), delim=";", show_col_types = FALSE)
  arguments<-as.list(arguments)                                                  
  if(!arguments$Changepoints%in%c("all","auto")) {                               
    if(grepl("-",arguments$Changepoints)) {
      Changepoints<-as.integer(unlist(strsplit(arguments$Changepoints,"-")))      
    } else {
      Changepoints<-try(as.integer(arguments$Changepoints))
      if(inherits(Changepoints,"try-error") | is.na(Changepoints)) {
        # in case of failure pick the first year.. can change this to "all"
        Changepoints<-1
      }
    }
    arguments$Changepoints<-Changepoints                                         
  }

  counts <- read_delim(paste0(inputfilespath, arguments$File, "_counts.csv"), delim=";", show_col_types = FALSE)#, stringsAsFactors = FALSE, header = TRUE, dec = ".",sep=";")   
  overview$first_year[overview$ss_combinations == arguments$File] <- min(counts$year[counts$count > 0], na.rm = TRUE) 
  overview$last_year[overview$ss_combinations == arguments$File]  <- max(counts$year[counts$count > 0], na.rm = TRUE)
  
  if(Changepoints==1) {
    # just linear model
    result <- tryCatch(
      {
        trim(count ~ site + year, data = counts, model = 1, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
      }, 
      error = warning)
    if(class(result) == "trim") {
      save(x = result,  file = paste0(inputfilespath, arguments$File, ".RData"))
      overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
      overview$success[overview$ss_combinations == arguments$File] <- "yes"
    } else {
      if(arguments$Presence_monthfactors == FALSE) {
        overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
        overview$error_1[overview$ss_combinations == arguments$File] <- result
      } 
    }
  } else {
    # as usual
    
    # Start with the most elaborate model and switch automatically to a more simple model when needed.
    
    result <- tryCatch(
      {
        # Attempt 1
        #cat(" Attempt 1\n") #	19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
        if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
          trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
          
        } else {
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
            if(weight=="eurostat"){
              print("entro")
              trim(count ~ site + year+covariate, data = counts, weights = "weights" ,model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
              
            }else{
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
            }          
          } else{
            if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
              trim(count ~ site + (year + month), data = counts,            model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
            } else{
              trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
            }
          }
        }
      }
      , error = warning)
    if (class(result) == "trim") {
      
      save(x = result,  file = paste0(inputfilespath, arguments$File, ".RData"))
      overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
      overview$success[overview$ss_combinations == arguments$File] <- "yes"
      
    } else {
      
      # First attempt failed, try a less elaborate model by setting serial correlation off. 
      # Also, when month factors are available, no changepoints are estimated in the next model.
      
      overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
      overview$error_1[overview$ss_combinations == arguments$File] <- result
      
      result <- tryCatch(
        {
          # attempt 2 
          #cat(" Attempt 2\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
            trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
          } else {
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              if(weight=="eurostat"){
                print("entro")
                trim(count ~ site + year+covariate, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
              }else{
                trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
              }
            } else{
              if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
                trim(count ~ site + (year + month), data = counts,                      model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
              } else{
                trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
              }
            }
          }
        }
        , error = warning)
      
      if (class(result) == "trim") {
        
        save(x = result,  file = paste0(inputfilespath, arguments$File, ".RData"))
        overview$attempt_2[overview$ss_combinations == arguments$File] <- "success"
        overview$success[overview$ss_combinations == arguments$File] <- "yes"
        
      } else {
        #cat(" Attempt 3\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
        
        # Second attempt also failed. Now try an even more simple model: no month factors, no changepoints, but serial correlation switched on. 
        # Note that no further options are available to include month factors in the model. 
        
        overview$attempt_2[overview$ss_combinations == arguments$File] <- "error"
        overview$error_2[overview$ss_combinations == arguments$File] <- result
        
        if (arguments$Presence_monthfactors == TRUE) {
          
          #cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
          
        }
        
        result <- tryCatch( 
          {
            # attempt 3
            
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              if(weight=="eurostat"){
                print("entro")
                trim(count ~ site + year + covariate, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
                
              }else{
                trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
              }
            } else{
              
              if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
              }
            }
          }
          , error = warning)
        if (class(result) == "trim") {
          
          save(x = result,  file = paste0(inputfilespath, arguments$File, ".RData"))
          overview$attempt_3[overview$ss_combinations == arguments$File] <- "success"
          overview$success[overview$ss_combinations == arguments$File] <- "yes"
          
        } else {
          if (arguments$Presence_monthfactors == FALSE) {
            overview$attempt_3[overview$ss_combinations == arguments$File] <- "error"
            overview$error_3[overview$ss_combinations == arguments$File] <- result
          }
          
          # Third attempt also failed. Now try the most simple model: no changepoints at all and no serial correlation. 
          
          result <- tryCatch(
            {
              # Final attempt
              #cat(" Attempt 4\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
              
              if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
                if(weight=="eurostat"){
                  print("entro")
                  trim(count ~ site + year +covariate, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
                  
                }else{
                  trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)
                }
              } else{
                if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                  trim(count ~ site + year, data = counts,                      model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5, verbose=F)    
                } 
              }
              
            }
            
            , error = warning)
          
          if (class(result) == "trim") {
            
            save(x = result,  file = paste0(inputfilespath, arguments$File, ".RData"))
            overview$attempt_4[overview$ss_combinations == arguments$File] <- "success"
            overview$success[overview$ss_combinations == arguments$File] <- "yes"
            
          } else {
            if (arguments$Presence_monthfactors == FALSE) {
              
              # When this analysis also fails, send a error message to screen.
              
              overview$attempt_4[overview$ss_combinations == arguments$File] <- "error"
              overview$error_4[overview$ss_combinations == arguments$File] <- result
              
              #cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
            }        
          }
        }  
      } 
    }  

    
    
  } # as usual
  
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  write_excel_csv2(overview, paste0(inputfilespath, "overview.csv"), na='') 
  
  pb$tick()
}


#processing.output("ST_inputs/","ST_outputs/")

outputfilespath <- "ST_outputs/"
outputs<-outputfilespath #paste0(rtrim_folder,"03_Outputs\\")
#####################################################################################################################
# SELECTING SUCCESSFUL RUNS.
#####################################################################################################################
# Determines which datafiles (combinations of species and stratum) have been analysed successfully. 
# This information is found in "overview.csv", which has been produced by the script that called rtrim.

overview <- read_csv2(paste0(inputfilespath, "overview.csv"))#,sep=";",dec=".",header = T)   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
listsuccessfulAnalyses <- overview$ss_combinations[overview$success == "yes"]
listSpeciesStratumCombinations <- paste(listsuccessfulAnalyses, "_arg_input_stratum.csv", sep = "")

# Determine how many combinations of species and stratum have been analysed successfully.
numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# File with trends and indices for each combination of species and stratum
all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)

#####################################################################################################################
# PROCESSING OUTPUT OF SUCCESFUL RUNS.
#####################################################################################################################

pb <- progress_bar$new(format = "  töötlen [:bar] :percent eta: :eta", total = numberSpeciesStratumCombinations, clear = FALSE, width= 60)
for (j in 1:numberSpeciesStratumCombinations) {

  arguments <- read_delim(paste0(inputfilespath, listSpeciesStratumCombinations[j]), delim=";", col_types=c("ciicllllil") )#, header = TRUE, stringsAsFactors = FALSE,sep=";",dec=".")
  counts <- read_delim(paste0(inputfilespath, arguments$File, "_counts.csv"), delim=";", col_types=c("cidd") )#, header = TRUE,dec=".",sep=";") 
  if(is.character(counts$count)==T){
    counts <- read_csv2(paste0(inputfilespath, arguments$File, "_counts.csv"), delim=";", col_types=c("cidd") )#, stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
  } 
  load(paste0(inputfilespath, arguments$File, ".RData")) # loads object with name "result"
  
  indices_TT_file <- make_Indices_TT_file(result = result)
  
  arg_output_file <- make_arg_output_file(arguments = arguments)
  
  indices_TT_file <- fill_Indices_TT_file(indices_TT_file = indices_TT_file, result = result, arguments = arguments)
  
  arg_output_file <- fill_arg_output_file(arg_output_file = arg_output_file, result = result, counts = counts)
  arg_output_file$Changepoints<-gsub(",\\ ", "-", arg_output_file$Changepoints) # changes the separator from function "results" from comma to dash # Javier Rivas Salvador 9/2/2023 - v1.6 & sub changed to gsub by Javier Rivas 23/3/2023 according to Meelis Leivits
  all_Indices_All_Trends <- fill_All_Indices_All_Trends(all_Indices_All_Trends = all_Indices_All_Trends, result = result, arguments = arguments, j = j, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  covariant_matrix <- vcov(result)
  
  if (arguments$Save_fitted_values){
    
    FI <- results(result)
    
  }
  
  #####################################################################################################################
  # WRITING OUTPUT FILES. 
  #####################################################################################################################
  # Indices and time totals.
  name_Indices_TT_file <- paste0(outputs,arguments$File, "_indices_TT.csv")
  write.table(indices_TT_file, name_Indices_TT_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots  - v1.6
  #cat("Indices_TT generated.","\n")
  ############################################################ 
  # Slopes entire period and arguments.
  name_arg_output_file <- paste0(outputs,arguments$File, "_arg_output.csv") 
  write.table(arg_output_file, name_arg_output_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
  #cat("Arg_output generated","\n")
  ############################################################
  # Covariant matrix.
  name_covariant_matrix <- paste0(outputs,arguments$File, "_ocv.csv")
  write.table(covariant_matrix, name_covariant_matrix, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
  #cat("ocv generated","\n")
  ############################################################
  # File with fitted values.
  if (arguments$Save_fitted_values){
    name_Fitted_Values_File <- paste0(outputs,arguments$File, "_fitted_values.csv")
    write.table(FI, name_Fitted_Values_File, row.names = FALSE,sep=";",dec=".")  #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    #cat("fitted values generated","\n")
  }
  ############################################################
  # 25/01/2022 Javier Rivas: Generating Tables 1-3 for the PECBMS coordinators - v1.4 IMPORTANT!!: version is extracted out of the path indicated in folder, if coordinators change the name it wont work.
  pb$tick()
}

#cat("Summarizing tables successfully created","\n")
all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
write.table(all_Indices_All_Trends, paste0(outputs,"All_Indices_All_Trends.csv"), row.names = FALSE,sep=";",dec=".")      

