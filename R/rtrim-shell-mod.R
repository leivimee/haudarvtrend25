#####################################################################################################################
# RTRIM-shell functions. Version RTRIM-shell_2.0
# Functions to support the  rtrim shell scripts.
# Marnix de Zeeuw rtrim@cbs.nl Statistics Netherlands 2019. 

# Last: modifications J.Rivas 13/03/2024
#####################################################################################################################
# 1/06/2020 Eva Silarova unified upper and lower cases in variable names - v1.1
make_Indices_TT_file <- function(result) {
  # Aim: this function creates an empty dataframe for indices and time totals for a species-stratum combination.
  
  indices_TT_file <- data.frame(Year = integer(result$nyear), 
                                # Years for which indices and time totals have been calculated.
                                
                                TT_model = numeric(result$nyear), 
                                # Time totals based on the estimated model.
                                
                                TT_model_SE = numeric (result$nyear), 
                                # Standard errors of model time totals.
                                
                                TT_imputed = numeric(result$nyear), 
                                # Imputed time totals.
                                
                                TT_imputed_SE = numeric(result$nyear), 
                                # Standard errors of imputed time totals.
                                
                                Index_model = numeric(result$nyear), 
                                # Indices based on the estimated model.
                                
                                Index_model_SE = numeric(result$nyear), 
                                # Standard errors of model indices.
                                
                                Index_imputed = numeric(result$nyear), 
                                # Indices based on the imputed time totals.
                                
                                Index_imputed_SE = numeric(result$nyear))
  # Standard errors of imputed indices.
  
  return(indices_TT_file)
}

#####################################################################################################################
make_arg_output_file <- function(arguments) {
  # Aim: this function creates an empty dataframe for slopes spanning the entire study period and a subperiod for a species-stratum combination.
  
  arg_output_file <- data.frame(N_sites = integer(1), 
                                # Number of unique sites.
                                
                                N_time_values = integer(1), 
                                # Number of unique years.
                                
                                N_observed_zero_counts = integer(1), 
                                # Number of zero counts.
                                
                                N_observed_positive_counts = integer(1), 
                                # Number of positive counts.
                                
                                N_missing_counts = integer(1), 
                                # Number of missing counts.
                                
                                N_counts_total = integer(1), 
                                # Total number of counts.
                                
                                Base_year_first_year = numeric(1),
                                Base_year_last_year = numeric(1),
                                
                                # Calendar year used as base year for indices. 
                                # If Base_year_first_year equals Base_year_last_year a single year is used as base year.
                                # If Base_year_first_year < Base_year_last_year, a period is used as base time.
                                # In the latter case, Base_year_first_year is the first year of the period.
                                
                                Changepoints = numeric(1), 
                                # Changepoints used.
                                
                                Overdispersion = numeric(1),
                                # Estimated overdispersion.
                                
                                Serial_correlation = numeric(1),
                                # Estimated serial correlation.
                                
                                Slope_imputed_mul = numeric(1),
                                # Multiplicative imputed slope for the entire period.
                                
                                Slope_imputed_mul_SE = numeric(1),
                                # Standard error of multiplicative imputed slope for the entire period.
                                
                                Slope_imputed_classification = character(1),
                                # Trend classification of multiplicative imputed slope for the entire period.
                                
                                Year_from = numeric(1),
                                # First year of the subperiod from which a slope has been calculated. 
                                # Last year of the subperiod is equal to the last year of the entire period.
                                # Note that it is assumed that the first year of the subperiod is closer to the present than the first year of the entire period.
                                
                                Slope_from_imputed_mul = numeric(1),
                                # Multiplicative imputed slope for the subperiod.
                                
                                Slope_from_imputed_mul_SE = numeric(1),
                                # Standard error of multiplicative imputed slope for the subperiod.
                                
                                Slope_from_imputed_classification = character(1)
                                # Trend classification of multiplicative imputed slope for the subperiod.
                                
  )
  
  arg_output_file$Year_from <- arguments$Year_from  
  #   Here, a value is already assigned to Year_from, to make the function that fills the arg_output file more simple. 
  #   The variabele 'arguments' is not needed due to this, when calling the function.
  
  arg_output_file$Base_year_first_year <- arguments$Base_year_first_year
  arg_output_file$Base_year_last_year  <- arguments$Base_year_last_year
  
  return(arg_output_file)
}

#####################################################################################################################
fill_Indices_TT_file <- function(indices_TT_file, result, arguments) {
  # Aim: this function fills the dataframe with indices and time totals for a species-stratum combination.    
  
  indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
  time_totals <- totals(result, which = "both")
  
  indices_TT_file$Year <- indices$time
  # The years over which the indices and time totals have been calculated.
  
  indices_TT_file$TT_model <- time_totals$fitted 
  # Time totals based on the estimated model. 
  
  indices_TT_file$TT_model_SE <- time_totals$se_fit
  # Standard errors of model time totals.
  
  indices_TT_file$TT_imputed <- time_totals$imputed
  # Imputed time totals.
  
  indices_TT_file$TT_imputed_SE <- time_totals$se_imp
  # Standard errors of imputed time totals.
  
  indices_TT_file$Index_model <- indices$fitted
  # Indices based on the estimated model.
  
  indices_TT_file$Index_model_SE <- indices$se_fit
  # Standard errors of model indices.
  
  indices_TT_file$Index_imputed <- indices$imputed
  # Indices based on the imputed time totals.
  
  indices_TT_file$Index_imputed_SE <- indices$se_imp
  # Standard errors of imputed indices.
  
  return(indices_TT_file)
}

#####################################################################################################################
fill_arg_output_file <- function(arg_output_file, result, counts) {
  # Aim: this function fills the dataframe with additional information.  
  
  overviewCounts <- count_summary(counts, count_col = "count", site_col = "site", year_col = "time") 
  # count_summary is a RTRIM function.
  
  slopes_imputed <- overall(result, which = "imputed")
  slopes_fitted  <- try( overall(result, which = "fitted") )
  if(inherits(slopes_fitted,"try-error")) {
    slopes_fitted  <- overall(result, which = "imputed") 
  }
  # The slopes over the entire period. Both variables are lists, a specific class of R-objects.
  
  slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = c(arg_output_file$Year_from)) 
  # Year_from is already available in the arg_output_file.
  slopes_subperiod_fitted <- try( overall(result, which = "fitted", changepoints = c(arg_output_file$Year_from)) )
  if(inherits(slopes_subperiod_fitted,"try-error")) {
    slopes_subperiod_fitted <- overall(result, which = "imputed", changepoints = c(arg_output_file$Year_from))
  }
  # The slopes over a subperiod. Both variables are lists.
  
  arg_output_file$N_sites <- result$nsite
  # Number of unique sites in the counts file used. 
  
  arg_output_file$N_time_values <- result$ntime
  # Number of unique years in the counts file used.
  
  arg_output_file$N_observed_zero_counts <- overviewCounts$zero_counts
  # Number of zero counts in the entire period, in the counts file used.
  
  arg_output_file$N_observed_positive_counts <- overviewCounts$positive_counts
  # Number of positive counts in the entire period, in the counts file used.
  
  arg_output_file$N_missing_counts <- overviewCounts$missing_counts
  # Number of missing counts in the entire period, in the counts file used.
  
  arg_output_file$N_counts_total <- overviewCounts$total_counts
  # Total number of counts in the entire period, in the counts file used.
  
  arg_output_file$Changepoints <- paste(result$changepoints, collapse = ", ")
  # Changepoints used to calculate the indices.
  
  arg_output_file$Overdispersion <- overdispersion(result)
  # Estimated overdispersion.
  
  arg_output_file$Serial_correlation <- serial_correlation(result)
  # Estimated serial correlation.
  
  arg_output_file$Slope_imputed_mul <- slopes_imputed$slope$mul
  # Multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_imputed_mul_SE <- slopes_imputed$slope$se_mul
  # Standard error of multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_imputed_classification <- slopes_imputed$slope$meaning
  # Trend classification of multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_from_imputed_mul <- slopes_subperiod_imputed$slope$mul[slopes_subperiod_fitted$slope$from == arg_output_file$Year_from]
  # Multiplicative imputed slope over the subperiod.
  
  arg_output_file$Slope_from_imputed_mul_SE <- slopes_subperiod_imputed$slope$se_mul[slopes_subperiod_fitted$slope$from == arg_output_file$Year_from]
  # Standard error of multiplicative imputed slope over the subperiod.
  
  arg_output_file$Slope_from_imputed_classification <- slopes_subperiod_imputed$slope$meaning[slopes_subperiod_imputed$slope$from == arg_output_file$Year_from]
  # Trend classification of multiplicative imputed slope over the subperiod.
  
  return(arg_output_file)
}

#####################################################################################################################
makeOverview <- function(listSpeciesStratumCombinations) {
  # Aim: this function makes an empty dataframe to list failure or success of rtrim runs.
  species_list_pcbms=read.table(paste0(rtrim_folder,"\\teatmik\\species_list.csv"),sep=";",dec=".",header = T)
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  overview <- data.frame(ss_combinations = character(numberSpeciesStratumCombinations), 
                         # All files of combination of species and stratum present in the working directory.
                         
                         species_group = character(numberSpeciesStratumCombinations),
                         # Short for specific species group.
                         
                         species_number = integer(numberSpeciesStratumCombinations),
                         # Unique species number.
                         
                         species_name = character(numberSpeciesStratumCombinations),
                         #Latin spcies name
                         
                         stratum_type = character(numberSpeciesStratumCombinations),
                         # Stratumtype. 1 = standard stratum (use this), 2 = combination of strata (not relevant here).
                         
                         stratum_number = integer(numberSpeciesStratumCombinations),
                         # Unique number for a stratum.
                         
                         first_year = integer(numberSpeciesStratumCombinations),
                         # The first year in which counts are available for this specific combination of species and stratum.
                         
                         last_year = integer(numberSpeciesStratumCombinations),
                         # The last year in which counts are available for this specific combination of species and stratum.
                         
                         success = character(numberSpeciesStratumCombinations), 
                         # yes if the analysis were successful, no if not.
                         
                         attempt_1 = character(numberSpeciesStratumCombinations), 
                         # Outcome of the first attempt. The attempt may be successful or not.
                         
                         attempt_2 = character (numberSpeciesStratumCombinations), 
                         # Outcome of the second attempt.
                         # Outcome is only available if the first attempt was not successful, else the value will be n.a.
                         
                         attempt_3 = character(numberSpeciesStratumCombinations), 
                         # Outcome of the third attempt.
                         # Outcome is only avaialble if earlier attempts were not succesful, else the value will be n.a.
                         
                         attempt_4 = character(numberSpeciesStratumCombinations),
                         # Outcome of the fourth attempt.
                         # Outcome is only avaialble if earlier attempts were not succesful, else the value will be n.a.
                         
                         error_1 = character(numberSpeciesStratumCombinations),
                         # Error message of the first attempt.
                         
                         error_2 = character(numberSpeciesStratumCombinations),
                         # Error message of the second attempt.
                         
                         error_3 = character(numberSpeciesStratumCombinations),
                         # Error message of the third attempt.
                         
                         error_4 = character(numberSpeciesStratumCombinations))
  # Error message of the last attempt.
  
  overview$ss_combinations <- gsub(listSpeciesStratumCombinations, pattern = "_arg_input_stratum.csv", replacement = "")
  
  overview$success <- "no"
  overview$attempt_1 <- "n.a."                                                                                       
  overview$attempt_2 <- "n.a."
  overview$attempt_3 <- "n.a."
  overview$attempt_4 <- "n.a."
  overview$error_1 <- "n.a."
  overview$error_2 <- "n.a."
  overview$error_3 <- "n.a."
  overview$error_4 <- "n.a."
  
  # retrieve species number, stratumtype and statumnumber from the filename 
  
  overview$species_group <- gsub(listSpeciesStratumCombinations, pattern = "_[0-9]+_[0-9]+_[0-9]+_arg_input_stratum.csv", replacement = "")
  
  without_Species_Group <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  overview$species_number <- as.integer(gsub(without_Species_Group, pattern = "_[0-9]+_[0-9]+_arg_input_stratum.csv", replacement = ""))
  
  without_Species_Group_without_Species_number <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  overview$stratum_type <- gsub(without_Species_Group_without_Species_number, pattern = "_[0-9]_arg_input_stratum.csv", replacement = "")
  
  without_Species_Group_without_Species_number_without_Stratum_type <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_[0-9]+_", replacement = "")
  
  overview$stratum_number <- gsub(without_Species_Group_without_Species_number_without_Stratum_type, pattern = "_arg_input_stratum.csv", replacement = "")
  
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  for(i in 1:nrow(overview)){
    sp_sub<-subset(species_list_pcbms, species_list_pcbms$Euring==overview$species_number[i])
    if(nrow(sp_sub)!=0){
      overview$species_name[i]<-sp_sub$Latin.Name
    }else{
      overview$species_name[i]<-"Not found in species list"
    }
  }
  
  return(overview)
}

#####################################################################################################################
make_All_Indices_All_Trends <- function(overview, listSpeciesStratumCombinations){
  # Aim: this function makes an empty dataframe to list indices and time totals of all species-stratum combinations.    
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  first_year_over_all_species <- min(overview$first_year)
  last_year_over_all_species <- max(overview$last_year)
  
  complete_period <- first_year_over_all_species:last_year_over_all_species
  columnnames <- as.character(complete_period)
  number_of_columns <- length(complete_period)
  number_of_rows <- numberSpeciesStratumCombinations * 4
  # Four rows for each combination of species and stratum (one row for each record type)
  
  temporary_matrix <- matrix(nrow = number_of_rows, ncol = number_of_columns)
  temporary_dataframe <- as.data.frame(temporary_matrix)
  colnames(temporary_dataframe) <- columnnames
  
  all_Indices_All_Trends <- data.frame(Year_of_analysis = integer(numberSpeciesStratumCombinations * 4),
                                       # Dataframe, named all_Indices_All_Trends, is created.
                                       # First column (Year_of_analysis): Most recent year for which counts are available.
                                       
                                       Species_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for a species.
                                       
                                       Stratum_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for a stratum.
                                       
                                       Recordtype_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for record types. There are four record types: 
                                       # 1: indices. 
                                       # 2: standard errors of indices.
                                       # 3: time totals.
                                       # 4: standard errors of time totals.
                                       
                                       Recordtype_name = character(numberSpeciesStratumCombinations),	
                                       # Name of the record type.
                                       
                                       N_sites = integer(numberSpeciesStratumCombinations),	
                                       # Number of unique sites in the counts file used.
                                       
                                       Slope_imputed_mul =	numeric(numberSpeciesStratumCombinations),
                                       # Multiplicative imputed slope over the entire period.
                                       
                                       Slope_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
                                       # Standard error of the multiplicative imputed slope over the entire period.
                                       
                                       Slope_imputed_classification = character(numberSpeciesStratumCombinations),	
                                       # Trend classification of the imputed slope over the entire period.
                                       
                                       Slope_from_imputed_mul = numeric(numberSpeciesStratumCombinations),
                                       # Multiplicative imputed slope over the subperiod.
                                       
                                       Slope_from_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
                                       # Standard error of the multiplicative imputed slope over the subperiod.
                                       
                                       Slope_from_imputed_classification = character(numberSpeciesStratumCombinations),	
                                       # Trend classification of the imputed slope over the subperiod.
                                       
                                       Year_from = integer(numberSpeciesStratumCombinations),	
                                       # First year of the subperiod over which a slope has been calculated. 
                                       
                                       Date_analysis = character(numberSpeciesStratumCombinations))
  # Date at which the analysis was executed.
  
  # retrieve species code, stratumtype and statumnumber from the filename 
  
  without_Species_Group <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  species_Code  <- as.integer(gsub(without_Species_Group, pattern = "_[0-9]_[0-9]+_arg_input_[a-zA-Z]+.csv", replacement = ""))
  
  without_Species_Group_Without_Species_Number <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  stratumtype_and_Stratumnumber <- gsub(without_Species_Group_Without_Species_Number, pattern = "_arg_input_[a-zA-Z]+.csv", replacement = "")
  
  stratumtype <- gsub(stratumtype_and_Stratumnumber, pattern = "_[0-9]+", replacement = "")
  
  stratum_Number <- gsub(stratumtype_and_Stratumnumber, pattern = "[0-9]+_", replacement = "")
  
  all_Indices_All_Trends$Recordtype_number <- rep(1:4, each = numberSpeciesStratumCombinations)
  # Unique number for record types. There are four record types: 
  # 1: indices. 
  # 2: standard errors of indices.
  # 3: time totals.
  # 4: standard errors of time totals.
  number_of_record_types <- length(unique(all_Indices_All_Trends$Recordtype_number))
  
  all_Indices_All_Trends$Species_number <- rep(species_Code, number_of_record_types)
  all_Indices_All_Trends$Stratum_number <- rep(stratum_Number, number_of_record_types)
  
  names_record_types <- c("indices", "se_indices", "time_totals", "se_time_totals")
  all_Indices_All_Trends$Recordtype_name <- rep(names_record_types, each = numberSpeciesStratumCombinations)	
  # Name of the record type.
  
  all_Indices_All_Trends$Slope_imputed_classification <- ""
  all_Indices_All_Trends$Slope_from_imputed_classification <- ""
  all_Indices_All_Trends$Date_analysis <- format(Sys.Date(), "%d-%m-%Y")
  # Date at which the analysis was done.
  
  all_Indices_All_Trends <- cbind(all_Indices_All_Trends, temporary_dataframe)
  
  return(all_Indices_All_Trends)  
}

#####################################################################################################################
fill_All_Indices_All_Trends <- function(result, arguments, j, listSpeciesStratumCombinations, all_Indices_All_Trends){
  # Aim: this function fills the dataframe with indices and time totals of all species-strata combinations.      
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
  time_totals <- totals(result, which = "both")
  slopes_imputed <- overall(result, which = "imputed")
  slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = c(arguments$Year_from))
  
  names_columns <- as.character(indices$time)
  position_colums <- colnames(all_Indices_All_Trends) %in% names_columns
  
  all_Indices_All_Trends$Year_of_analysis <- max(indices$time)
  # Most recent year in which counts are available.
  
  all_Indices_All_Trends$N_sites[j] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 2] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 3] <- result$nsite
  # Number of unique sites in the counts file used to call rtrim.
  
  all_Indices_All_Trends$Slope_imputed_mul[j] <- slopes_imputed$slope$mul
  # Multiplicative imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_imputed_mul_SE[j] <-	slopes_imputed$slope$se_mul
  # Standard error of the multiplicative imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_imputed_classification[j] <-	slopes_imputed$slope$meaning
  # Trend classification of the imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_from_imputed_mul[j] <- slopes_subperiod_imputed$slope$mul[slopes_subperiod_imputed$slope$from == arguments$Year_from]
  # Multiplicative imputed slope over the subperiod.
  
  all_Indices_All_Trends$Slope_from_imputed_mul_SE[j] <- slopes_subperiod_imputed$slope$se_mul[slopes_subperiod_imputed$slope$from == arguments$Year_from]	
  # Standard error of the multiplicative imputed slope over the subperiod.
  
  all_Indices_All_Trends$Slope_from_imputed_classification[j] <- slopes_subperiod_imputed$slope$meaning[slopes_subperiod_imputed$slope$from == arguments$Year_from]	
  # Trend classification of the imputed slope over the subperiod.
  
  all_Indices_All_Trends$Year_from[j] <-	arguments$Year_from
  # First year of the subperiod over which a slope has been calculated. 
  
  all_Indices_All_Trends[j, position_colums] <- round(100 * indices$imputed, 1)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations, position_colums] <- round(100 * indices$se_imp, 2)    #12/6/2020 Eva Silarova added 100 into <- round(indices$se_imp, 2) to have SE corresponding with index in % - v1.2
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 2, position_colums] <- round(time_totals$imputed, 2)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 3, position_colums] <- round(time_totals$se_imp, 2)
  
  return(all_Indices_All_Trends)  
}






#####################################################################################################################
#This function run the RTRIM models for each species and save the result as .R file
rtrim.strata<- function(inputfilespath="./", weight=0) {
  pb <- progress_bar$new(format = "  töötlen [:bar] :percent eta: :eta", total = numberSpeciesStratumCombinations, clear = FALSE, width= 60)
  #cat("##############################################","\n")
  #cat("Starting Rtrim Strata","\n")
  #cat("##############################################","\n")
  for (j in 1:numberSpeciesStratumCombinations) {
    
    #cat("Processing", listSpeciesStratumCombinations[j], "\n") # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for - v1.6
    
    #j<- 1
    
    # The file with arguments contains the information to analyse the counts for a particular combination of species and stratum. 
    # The arguments are used when calling the function 'rtrim'
    arguments <- read_csv2(paste0(inputfilespath,listSpeciesStratumCombinations[j]), show_col_types = FALSE)#, header = TRUE, stringsAsFactors = FALSE)
    
    arguments<-as.list(arguments)                                                  # Enables the script to read the change points properly, 
    if(!arguments$Changepoints%in%c("all","auto")) {                               # when they are specified in the ?_arg_input_stratum.csv? file as comma-separated numbers.
      if(grepl("-",arguments$Changepoints)) {
        Changepoints<-as.integer(unlist(strsplit(arguments$Changepoints,"-")))       # Function was created 20/03/2022 by Dario Massimino -> modified 11/04/2022 by Meelis Leivits -> modified 20/07/2022 by John Kennedy - v1.6
      } else {
        Changepoints<-try(as.integer(arguments$Changepoints))
        if(inherits(Changepoints,"try-error") | is.na(Changepoints)) {
          # in case of failure pick the first year.. can change this to "all"
          Changepoints<-1
        }
      }
      arguments$Changepoints<-Changepoints                                         # Function was adjusted to use of "-" in the changepoints # 9/02/2023 Javier Rivas Salvador - v1.6
    } 
    
    # The file with counts contains the counts for a particular combination of species and stratum.
    # Weights may also be present in this file.
    #13/01/2023 Javier Rivas: I have modified the following lines to ensure that independently of the decimal symbol the counts$count variable is loaded as numeric
    counts <- read_csv2(paste0(inputfilespath, arguments$File, "_counts.csv"), show_col_types = FALSE)#, stringsAsFactors = FALSE, header = TRUE, dec = ".",sep=";")   
    
    if(is.character(counts$count)==T | is.character(counts$weight)==T){
      counts <- read_csv2(paste0(inputfilespath,arguments$File, "_counts.csv"), show_col_types = FALSE)#, stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
    }
    if(weight=="eurostat"){
      counts$covariate<-as.factor(counts$covariate)
    }
    # 13/01/2023 end of modifications 
    # 30/05/22 Alena Jechumtal Skalova & Martin Stjernman corrected mismatches in the species-codes and their first last year.-v1.6 
    #   Moreover, the code evaluates the first year of detection and not the first year of the scheme.-v1.6
    overview$first_year[overview$ss_combinations == arguments$File] <- min(counts$year[counts$count > 0], na.rm = TRUE) 
    overview$last_year[overview$ss_combinations == arguments$File]  <- max(counts$year[counts$count > 0], na.rm = TRUE)
    # Original code
    # overview$first_year[j] <- min(counts$year, na.rm = TRUE)
    # overview$last_year[j]  <- max(counts$year, na.rm = TRUE)
    # End of the editation 30/05/2022 
    
    # First and last year for each combination of stratum and species is stored for later use.
    
    #####################################################################################################################
    # RUNNING RTRIM.
    #####################################################################################################################
    
    # Start with the most elaborate model and switch automatically to a more simple model when needed.
    
    result <- tryCatch(
      {
        # Attempt 1
        #cat(" Attempt 1\n") #	19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
        if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
          trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          
        } else {
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
            if(weight=="eurostat"){
              print("entro")
              trim(count ~ site + year+covariate, data = counts, weights = "weights" ,model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              
            }else{
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
            }          
          } else{
            if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
              trim(count ~ site + (year + month), data = counts,            model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            } else{
              trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
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
            trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          } else {
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              if(weight=="eurostat"){
                print("entro")
                trim(count ~ site + year+covariate, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              }else{
                trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              }
            } else{
              if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
                trim(count ~ site + (year + month), data = counts,                      model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              } else{
                trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
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
                trim(count ~ site + year + covariate, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
                
              }else{
                trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              }
            } else{
              
              if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
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
                  trim(count ~ site + year +covariate, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
                  
                }else{
                  trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
                }
              } else{
                if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                  trim(count ~ site + year, data = counts,                      model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
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
  }
  
  #####################################################################################################################
  # WRITING OVERVIEW OF RTRIM SUCCESSES and FAILURES.  
  #####################################################################################################################
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  write_excel_csv2(overview, paste0(inputfilespath, "overview.csv"), na='')   #13/01/2023 Javier Rivas: Changes to produce outputs with dots
  #cat("Overview table generated successfully","\n")
  #cat("##############################################","\n")
  #cat("End of Rtrim Strata","\n")
  #cat("##############################################","\n")
  pb$tick()
}



#This function creates the outputs of RTRIM
processing.output<-function(inputfilespath,outputfilespath){
  # inputfilespath<-"ST_inputs/";outputfilespath<-"ST_outputs/"
  #cat("##############################################","\n")
  #cat("Start of Processing Outputs","\n")
  #cat("##############################################","\n")
  outputs<-outputfilespath #paste0(rtrim_folder,"03_Outputs\\")
  #####################################################################################################################
  # SELECTING SUCCESSFUL RUNS.
  #####################################################################################################################
  # Determines which datafiles (combinations of species and stratum) have been analysed successfully. 
  # This information is found in "overview.csv", which has been produced by the script that called rtrim.
  
  overview <- read_delim(paste0(inputfilespath, "overview.csv"), delim=";")#,sep=";",dec=".",header = T)   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
  listsuccessfulAnalyses <- overview$ss_combinations[overview$success == "yes"]
  listSpeciesStratumCombinations <- paste(listsuccessfulAnalyses, "_arg_input_stratum.csv", sep = "")
  
  # Determine how many combinations of species and stratum have been analysed successfully.
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  # File with trends and indices for each combination of species and stratum
  all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  #####################################################################################################################
  # PROCESSING OUTPUT OF SUCCESFUL RUNS.
  #####################################################################################################################
  pb <- txtProgressBar(min=1, max=numberSpeciesStratumCombinations, style = 3, char="=")
  for (j in 1:numberSpeciesStratumCombinations) {
    # j<-1
    #cat("Generating outputs for sp: ",listSpeciesStratumCombinations[j],"\n" )
    # The file with arguments contains the arguments to run the analysis for a particular combination of species and stratum (stratum is e.g. a region). 
    # These arguments are used when calling the function 'rtrim'.
    arguments <- read_delim(paste0(inputfilespath, listSpeciesStratumCombinations[j]), delim=";", col_types=c("ciicllllil") )#, header = TRUE, stringsAsFactors = FALSE,sep=";",dec=".")
    counts <- read_delim(paste0(inputfilespath, arguments$File, "_counts.csv"), delim=";", col_types=c("cidd") )#, header = TRUE,dec=".",sep=";") 
    if(is.character(counts$count)==T){
      counts <- read_csv2(paste0(inputfilespath, arguments$File, "_counts.csv"), delim=";", col_types=c("cidd") )#, stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
    } 
    load(paste0(inputfilespath, arguments$File, ".RData")) # loads object with name "result"
    
    #####################################################################################################################
    # CREATING FILES FOR LATER USE. 
    #####################################################################################################################
    
    # Several output files are created (dataframes). Filling the files with output is done at a later stage.
    # To create the file containing indices and time totals (indices_TT_file) requires the file with results.
    # The file arg_output contains the slopes and arguments used to run the rtrim function.
    
    indices_TT_file <- make_Indices_TT_file(result = result)
    
    arg_output_file <- make_arg_output_file(arguments = arguments)
    
    #####################################################################################################################
    # FILLING OUTPUT.
    # Several functions are used to fill the output files.
    # These functions can be found in "functions_Shell_Rtrim.r".
    # This script goes to that file to include them here, which makes it easier to see what is done in this script. 
    #####################################################################################################################
    
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
    setTxtProgressBar(pb, value = j)
  }
  close(pb)
  #cat("Summarizing tables successfully created","\n")
  all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
  write.table(all_Indices_All_Trends, paste0(outputs,"All_Indices_All_Trends.csv"), row.names = FALSE,sep=";",dec=".")           # 09/04/2021 Eva Silarova changed "all" to capital "All": "all_Indices_All_Trends.csv" was changed to "All_Indices_All_Trends.csv" - v1.3
  #cat("All indices all trends successfully created","\n")
  #18/09/2023 Javi: Creation of the Schedule table
  
  
  #cat("##############################################","\n")
  #cat("End of Processing output","\n")
  #cat("##############################################","\n")
}
