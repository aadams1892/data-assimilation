################################################################################
#' Costfunction for optimizing parameter values
#' @description This function is a cost function that evaluates model performance
#' @param normParameterValues Vector with normalized PFT-specific parameter values (dimensionless)
#' @param lowerBound List with lower bound parameter values 
#' @param upperBound List with upper bound parameter values 
#' @param parameterValueLength Vector with the number of parameter values per parameter
#' @param parameterNames List of parameter names
#' @param parameterFile Path to CLASSIC parameter file (run_parameters.txt)
#' @param mod.list List of paths to CLASSIC model output
#' @param ref.list List of paths to reference output
#' @param ref.id.list List with reference IDs
#' @param ref.unit.conv.list List of factors used for converting reference data units
#' @param run_classic_file Path to CLASSIC submission file (run_classic.sh)
#' @param modelOutputFolder Path to model output folder
#' @return Textfiles (daisyOutput_Ref-ID) with scores (columns 1-6) and parameter values (remaining columns)
#' @export

costFunc <- function(identifier, indivFile, argValsFile, argTypesFile, 
                     argLengthsFile, argNamesFile, farmName) {
  
    # identifier - the line of indivFile to read from (the individual whose fitness will be calculated)
    # indivFile - file containing all the individuals
    # argValsFile - file containing the values of the arguments
    # argTypesFile - file containing the data types of the arguments. Format: type of whole argument, type of element(s)
    # argLengthsFile - file containing the lengths of the arguments
    # argNamesFile - file containing the names of the arguments

    # Source scripts
    source('/home/aadams/daisy/daisy/R/editParameterFile.R')
    source('/home/aadams/daisy/daisy/R/getDataModRef.R')
    source('/home/aadams/daisy/daisy/R/getParameterValues.R')
    source('/home/aadams/daisy/daisy/R/intFun.R')
    source('/home/aadams/daisy/daisy/R/postfit.R')
    
    # Read in the data for the arguments.
    argVals <- readLines(argValsFile)
    argTypes <- readLines(argTypesFile)
    argLengths <- readLines(argLengthsFile)
    argNames <- readLines(argNamesFile)
    
    # Format arguments correctly.
    for (index in 1:length(argVals)) {
      
      types <- unlist(strsplit(argTypes[index], " "))
      dataType <- types[1]
      elementType <- types[2]
      
      if (dataType != "list") {
      
        # Split the argument into a list if it has a length > 1.
        if (argLengths[index] > 1) {
          arg <- strsplit(argVals[index], " ")
          
        # Otherwise, set it just as it is.
        } else {
          arg <- argVals[index]
        }
        
        # Format elements.
        # Double
        if (elementType == "double") {
          arg <- lapply(arg, as.double)
        # Integer
        } else if (elementType == "integer") {
          arg <- lapply(arg, as.integer)
        # Character
        } else if (elementType == "character") {
          arg <- lapply(arg, trimws)
        }
        
        # Unlist.
        arg <- unlist(arg)
        
      # Data type is a list.
      } else if (dataType == "list") {
        arg <- unlist(strsplit(argVals[index], " "))
        
        # If element type is double, convert to double.
        if (elementType == "double") {
          arg <- unlist(lapply(arg, as.double))
        # If element type is character, trim white space
        }
      }
        
      # Assign the argument its name.
      assign(argNames[index], arg)
    }
    
    # MANAGE PARAMETERS
    
    # modelOutputFolder
    # Set runFolder
    runFolder <- sprintf("RUN%d", as.integer(identifier))
    filePath <- paste0("/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/", farmName, "/", runFolder)
    # Set modelOutputFolder
    modelOutputFolder <- filePath
    
    # run_classic_file
    # Copy run_classic_file to the RUNX folder and set run_classic_folder.
    system(sprintf("cp %s run_classic.sh", run_classic_file))
    Sys.sleep(1)
    run_classic_file <- paste0(filePath, "/run_classic.sh")
    # Change the outputDir variable in run_classic_file to be modelOutputFolder
    system(paste0("sed -i 's|/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/simulations/$simulationID|", 
                  modelOutputFolder, "/$simulationID|g' ", run_classic_file))
    
    # parameterFile
    # Move the parameterFile to the RUNX folder.
    system(paste0("cp ", parameterFile, " /home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/", 
                  farmName, "/", runFolder))
    # Set parameterFile to this file.              
    parameterFile <- paste0(filePath, "/run_parameters.txt")               
    
    # dir.mod
    # Set dir.mod.
    dir.mod <- paste0(filePath, "/daisyRun/netcdf_files")
    # Set .nc files
    nc.mod01 <- file.path(dir.mod, "gpp_monthly.nc")
    nc.mod02 <- file.path(dir.mod, "lai_monthly.nc")
    nc.mod03 <- file.path(dir.mod, "hfls_monthly.nc")
    nc.mod04 <- file.path(dir.mod, "hfss_monthly.nc")
    nc.mod05 <- file.path(dir.mod, "ts_monthly.nc")
    nc.mod06 <- file.path(dir.mod, "albs_monthly.nc")
    
    # Set mod.list
    mod.list <- list(
      nc.mod01, nc.mod01,
      nc.mod02, nc.mod02, nc.mod02,
      nc.mod03, nc.mod03,
      nc.mod04, nc.mod04,
      nc.mod05,
      nc.mod06, nc.mod06, nc.mod06
    )
    
    Sys.sleep(1)
    
    n <- length(parameterNames)
    
    .libPaths(c('/home/aadams/projects/def-cseiler-ab/aadams/AMBER/renv/', .libPaths()))
    
    library(readr)
    
    # Read in the data for the arguments
    normParameterValues <- read_lines(file = indivFile, skip = as.integer(identifier)-1, n_max = 1)
    # Split the individual up.
    normParameterValues <- unlist(strsplit(normParameterValues, " "))
    # Convert from character to double (requires unlisting again).
    normParameterValues <- unlist(lapply(normParameterValues, as.double))

    # Convert vector to list, where each element represents one parameter
    breakpoints <- rep(seq_along(parameterValueLength), times = parameterValueLength)
    normParameterValuesList <- split(x = normParameterValues, f = breakpoints)
    
    upperBoundList <- split(x = upperBound, f = breakpoints)
    lowerBoundList <- split(x = lowerBound, f = breakpoints)

    # The loop edits the parameter file for each parameter

    daisyOutput <- numeric()

    for (i in 1:n) {

        # Get the parameter names
        parameterName <- parameterNames[[i]]

        # un-normalize parameter values
        parameterValues <- intFun.unnormalize(normParameterValuesList[[i]], upperBoundList[[i]], lowerBoundList[[i]])
        
        daisyOutput <- c(daisyOutput, parameterValues)

        # Get the default parameter values
        defaultParameterValues <- getParameterValues(parameterFile, parameterName)
        
        non_zero_indices <- which(defaultParameterValues != 0, arr.ind = TRUE)
        newValues <- defaultParameterValues
        
        # Replace prior with new parameter values, skipping all instances where
        # the parameter values equal zero in the parameter file
        
        newValues[non_zero_indices] <- defaultParameterValues[non_zero_indices] -
          defaultParameterValues[non_zero_indices] + parameterValues[!is.na(parameterValues)]
        
        # Overwrite the parameter file with the new parameter values
        lines <- editParameterFile(parameterFile, parameterName, parameterValues = newValues)
        writeLines(lines, parameterFile)

    }

    # Run CLASSIC with new parameter file
    Sys.sleep(10)
    system(run_classic_file)

    # Wait until CLASSIC finishes
    while (!file.exists(mod.list[[1]])) {
        Sys.sleep(10)
    }
    Sys.sleep(10)
    print("CLASSIC run completed")

    # Obtain model and reference data.

    rm(i, n)

    # Check whether the length of the model and reference lists are equal
    condition <- length(mod.list) != length(ref.list)
    if (condition) {
        stop("Error: The length of the model and reference lists are not equal")
    }

    n <- length(ref.list)
    scores <- numeric(n)

    # Loop through all reference datasets
    for (i in 1:n) {

        nc.mod <- mod.list[[i]]
        nc.ref <- ref.list[[i]]
        ref.id <- ref.id.list[[i]]
        ref.unit.conv <- ref.unit.conv.list[[i]]

        myLevel <- 1

        mod <- raster::brick(nc.mod, level = myLevel)
        suppressWarnings(mod <- raster::rotate(mod))

        ref <- raster::brick(nc.ref)
        suppressWarnings(ref <- raster::rotate(ref))

        # model data dates
        dates.mod <- intFun.getZ(nc.mod)
        mod <- raster::setZ(mod, dates.mod)
        names(mod) <- dates.mod
        dates.mod <- format(as.Date(dates.mod), "%Y-%m")  # only year and month
        start.date.mod <- min(dates.mod)
        end.date.mod <- max(dates.mod)

        # reference data dates
        dates.ref <- intFun.getZ(nc.ref)
        ref <- raster::setZ(ref, dates.ref)
        names(ref) <- dates.ref
        dates.ref <- format(as.Date(dates.ref), "%Y-%m")  # only year and month
        start.date.ref <- min(dates.ref)
        end.date.ref <- max(dates.ref)

        # find common time period
        start.date <- max(start.date.mod, start.date.ref)
        end.date <- min(end.date.mod, end.date.ref)

        # subset common time period
        mod <- mod[[which(format(as.Date(raster::getZ(mod)), "%Y-%m") >= start.date &
            format(as.Date(raster::getZ(mod)), "%Y-%m") <= end.date)]]
        ref <- ref[[which(format(as.Date(raster::getZ(ref)), "%Y-%m") >= start.date &
            format(as.Date(raster::getZ(ref)), "%Y-%m") <= end.date)]]

        # get layer names
        mod.names <- base::names(mod)
        ref.names <- base::names(ref)

        # unit conversion if appropriate
        ref <- ref * ref.unit.conv

        #---------------------------------------------------------------------------

        # II Statistical analysis

        #---------------------------------------------------------------------------

        # (1) Bias

        #---------------------------------------------------------------------------
        # create a mask to excludes all grid cells that the model and reference
        # data do not have in common.  This mask varies in time.
        mask <- (mod * ref)
        mask <- mask - mask + 1
        mod <- mod * mask
        names(mod) <- mod.names  # this adds the corresponding dates
        ref <- ref * mask
        names(ref) <- ref.names  # this adds the corresponding dates
        # now mod and ref are based on the same grid cells

        mod.mean <- raster::mean(mod, na.rm = TRUE)  # time mean
        ref.mean <- raster::mean(ref, na.rm = TRUE)  # time mean
        bias <- mod.mean - ref.mean  # time mean

        mod.sd <- raster::calc(mod, fun = sd, na.rm = TRUE)  # standard deviation of model data
        ref.sd <- raster::calc(ref, fun = sd, na.rm = TRUE)  # standard deviation of reference data
        epsilon_bias <- abs(bias)/ref.sd
        epsilon_bias[epsilon_bias == Inf] <- NA  # relative error
        bias.score <- exp(-1 * epsilon_bias)  # bias score as a function of space
        S_bias <- mean(raster::getValues(bias.score), na.rm = TRUE)  # scalar score

        #---------------------------------------------------------------------------

        # (2) root mean square error (rmse)

        #---------------------------------------------------------------------------
        ESM.mode <- FALSE
        # Set ESM.mode to FALSE when forcing model with quasi-observed data
        # (e.g.  reanalysis)
        if (ESM.mode == FALSE) {
            rmse <- intFun.rmse(mod, ref)  # rmse
            mod.anom <- mod - mod.mean  # anomaly
            ref.anom <- ref - ref.mean  # anomaly
        }

        # Set ESM.mode to TRUE for online experiments or when forcing model
        # with earth system model data
        if (ESM.mode == TRUE) {
            index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
            index <- as.numeric(index)
            mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
            ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)

            mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
            ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec

            # Ensure that the order is correct
            JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5",
                "index_6", "index_7", "index_8", "index_9", "index_10", "index_11",
                "index_12")
            mod.clim.mly <- mod.clim.mly[[JanToDec]]
            ref.clim.mly <- ref.clim.mly[[JanToDec]]

            rmse <- intFun.rmse(mod.clim.mly, ref.clim.mly)  # rmse

            mod.anom <- mod.clim.mly - mod.mean  # anomaly
            ref.anom <- ref.clim.mly - ref.mean  # anomaly

        }

        crmse <- intFun.crmse(mod.anom, ref.anom)  # centralized rmse

        #-------------------------------------------------------------
        epsilon_rmse <- crmse/ref.sd
        epsilon_rmse[epsilon_rmse == Inf] <- NA  # relative error
        rmse.score <- exp(-1 * epsilon_rmse)  # rmse score as a function of space
        S_rmse <- mean(raster::getValues(rmse.score), na.rm = TRUE)

        #---------------------------------------------------------------------------

        # (3) phase shift

        #---------------------------------------------------------------------------

        index <- format(as.Date(names(ref), format = "X%Y.%m.%d"), format = "%m")
        index <- as.numeric(index)
        mod.clim.mly <- raster::stackApply(mod, index, fun = raster::mean)
        ref.clim.mly <- raster::stackApply(ref, index, fun = raster::mean)

        mod.cycle <- mod.clim.mly  # does not necessarily start in Jan or end in Dec
        ref.cycle <- ref.clim.mly  # does not necessarily start in Jan or end in Dec

        # Ensure that the order is correct
        JanToDec <- c("index_1", "index_2", "index_3", "index_4", "index_5", "index_6",
            "index_7", "index_8", "index_9", "index_10", "index_11", "index_12")
        mod.clim.mly <- mod.clim.mly[[JanToDec]]
        ref.clim.mly <- ref.clim.mly[[JanToDec]]

        bias.clim.mly <- mod.clim.mly - ref.clim.mly

        # find month of seasonal peak

        # In most cases, we are interested in the timing of the seasonal
        # maximum value

        # In some cases, however, the seasonal peak is a minimum, e.g. NEE =
        # RECO - GPP

        phaseMinMax <- "phaseMax"
        #
        if (phaseMinMax == "phaseMax") {
            mod.max.month <- raster::which.max(mod.clim.mly)
            ref.max.month <- raster::which.max(ref.clim.mly)
        }

        if (phaseMinMax == "phaseMin") {
            mod.max.month <- raster::which.min(mod.clim.mly)
            ref.max.month <- raster::which.min(ref.clim.mly)
        }

        # get shortest time distance between these months
        abs.diff <- abs(mod.max.month - ref.max.month)  # absolute difference from 0 to 12 months
        phase <- raster::calc(abs.diff, intFun.theta)  # shortest distance from 0 to 6 months (theta)
        phase.score <- 0.5 * (1 + cos(2 * pi * phase/12))  # score from 0 (6 months) to 1 (0 months)
        S_phase <- mean(raster::getValues(phase.score), na.rm = TRUE)  # scalar score

        #---------------------------------------------------------------------------

        # (4) interannual variability

        #---------------------------------------------------------------------------

        years <- floor(raster::nlayers(mod)/12)  # total number of years
        months <- years * 12  # number of months considering complete years only
        mod.fullyear <- raster::subset(mod, 1:months)
        ref.fullyear <- raster::subset(ref, 1:months)
        c.mod <- raster::calc(mod.cycle, fun = function(x) {
            rep(x, years)
        })  # climatological cycle for all months (mod)
        c.ref <- raster::calc(ref.cycle, fun = function(x) {
            rep(x, years)
        })  # climatological cycle for all months (ref)
        mod.iav <- sqrt(raster::mean((mod.fullyear - c.mod)^2, na.rm = TRUE))  # interannual variability  (mod)
        ref.iav <- sqrt(raster::mean((ref.fullyear - c.ref)^2, na.rm = TRUE))  # interannual variability  (ref)
        # set values close to zero to NA
        # ref.iav.na <- ref.iav
        # ref.iav.na[ref.iav.na < 10^(-5)] <- NA
        # ref.iav.na[ref.iav.na < 10^(-10)] <- NA

        # epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav.na
        epsilon_iav <- abs((mod.iav - ref.iav))/ref.iav
        
        epsilon_iav[epsilon_iav == Inf] <- NA  # I changed Eq. 26 so that epsilon_iav >= 0
        iav.score <- exp(-1 * epsilon_iav)  # iav score as a function of space
        S_iav <- mean(raster::getValues(iav.score), na.rm = TRUE)  # scalar score (not weighted)

        #---------------------------------------------------------------------------

        # (5) dist

        #---------------------------------------------------------------------------

        mod.sigma.scalar <- stats::sd(raster::getValues(mod.mean), na.rm = TRUE)  # standard deviation of period mean data
        ref.sigma.scalar <- stats::sd(raster::getValues(ref.mean), na.rm = TRUE)  # standard deviation of period mean data
        sigma <- mod.sigma.scalar/ref.sigma.scalar
        y <- raster::getValues(mod.mean)
        x <- raster::getValues(ref.mean)
        reg <- stats::lm(y ~ x)
        R <- sqrt(summary(reg)$r.squared)
        S_dist <- 2 * (1 + R)/(sigma + 1/sigma)^2

        #---------------------------------------------------------------------------

        # overall scores

        #---------------------------------------------------------------------------

        S <- mean(c(S_bias, S_rmse, S_phase, S_iav, S_dist), na.rm = TRUE)

        timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        
        # Since the cost.fun file gets executed from a different directory in parallel
        # than sequentially, new score files get created every run. To avoid this,
        # the directory is changed when creating these score files and then reverted
        # to what it was before.
        origDir <- getwd()
        setwd("/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation/")
        Sys.sleep(1)
        # Write current scores and parameter values to text file for analysis
        # purpose
        fileName <- paste("daisyOutput", ref.id, sep = "_")
        write.table(x = matrix(c(S, S_bias, S_rmse, S_phase, S_iav, S_dist, daisyOutput,
            timeStamp), nrow = 1), file = fileName, append = TRUE, row.names = FALSE,
            col.names = FALSE)
        # Revert to original directory.
        setwd(origDir)
        Sys.sleep(1)
        
        scores[i] <- S
    }

    # Move outputs to a new folder so that they are not overwritten

    NewModelOutputFolder <- paste(modelOutputFolder, timeStamp, sep = "_")
    moveModelOutputFolder <- paste("mv", modelOutputFolder, NewModelOutputFolder,
        sep = " ")
    system(moveModelOutputFolder)

    # Calculate the mean score across all evaluations
    S <- mean(scores, na.rm = TRUE)
    return(S)

}
