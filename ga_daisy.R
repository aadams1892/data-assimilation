##############################################################################
#                                                                            #
#                        GENETIC ALGORITHMS in R                             #
#                                                                            #
##############################################################################

ga_daisy <- function(type = c("binary", "real-valued", "permutation"), 
               fitness, ...,
               lower, upper, nBits,
               population = gaControl(type)$population,
               selection = gaControl(type)$selection,
               crossover = gaControl(type)$crossover, 
               mutation = gaControl(type)$mutation,
               popSize, 
               pcrossover = 0.8, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               updatePop = FALSE,
               postFitness = NULL,
               maxiter,
               run,
               maxFitness = 1, # Inf
               parallel,
               jobTime,
               farmName,
               dataAssimPath,
               previousObject = NULL,
               names = NULL,
               suggestions = NULL,
               optim = FALSE,
               optimArgs = list(method = "L-BFGS-B", 
                                poptim = 0.05,
                                pressel = 0.5,
                                control = list(fnscale = -1, maxit = 100)),
               keepBest = FALSE,
               keepOutput = FALSE,
               archiveFolder,
               finalOutputFolder = NULL,
               monitor = if(interactive()) gaMonitor else FALSE,
               seed = NULL) 
{

  call <- match.call()
  
  type <- match.arg(type, choices = eval(formals(ga)$type))
  
  if(!is.function(population)) population <- get(population)
  if(!is.function(selection))  selection  <- get(selection)
  if(!is.function(crossover))  crossover  <- get(crossover)
  if(!is.function(mutation))   mutation   <- get(mutation)
  
  if(missing(fitness))
    { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
    { stop("A fitness function must be provided") }
  if(popSize < 10) 
    { warning("The population size is less than 10.") }
  if(maxiter < 1) 
    { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
    { stop("The elitism cannot be larger that population size.") }
  elitism <- as.integer(elitism)
  if(pcrossover < 0 | pcrossover > 1)
    { stop("Probability of crossover must be between 0 and 1.") }
  if(is.numeric(pmutation))
  { 
    if(pmutation < 0 | pmutation > 1)
      { stop("If numeric probability of mutation must be between 0 and 1.") }
    else if(!is.function(population))
           { stop("pmutation must be a numeric value in (0,1) or a function.") }
  }

  # check for min and max arguments instead of lower and upper
  callArgs <- list(...)
  if(any("min" %in% names(callArgs)))
  {
    lower <- callArgs$min
    callArgs$min <- NULL
    warning("'min' arg is deprecated. Use 'lower' instead.")
  }
  if(any("max" %in% names(callArgs)))
  {
    upper <- callArgs$max
    callArgs$max <- NULL
    warning("'max' arg is deprecated. Use 'upper' instead.")
  }

  if(missing(lower) & missing(upper) & missing(nBits))
    { stop("A lower and upper range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }

  # check GA search type 
  switch(type, 
         "binary"      = { nBits <- as.vector(nBits)[1]
                           lower <- upper <- NA
                           nvars <- nBits 
                           if(is.null(names))
                             names <- paste0("x", 1:nvars)
                         },
         "real-valued" = { lnames <- names(lower)
                           unames <- names(upper)
                           lower <- as.vector(lower)
                           upper <- as.vector(upper)
                           nBits <- NA
                           if(length(lower) != length(upper))
                             stop("lower and upper must be vector of the same length!")
                           nvars <- length(upper)
                           if(is.null(names) & !is.null(lnames))
                             names <- lnames
                           if(is.null(names) & !is.null(unames))
                             names <- unames
                           if(is.null(names))
                             names <- paste0("x", 1:nvars)
                         },
         "permutation" = { lower <- as.vector(lower)[1]
                           upper <- as.vector(upper)[1]
                           nBits <- NA
                           nvars <- length(seq.int(lower,upper))
                           if(is.null(names))
                             names <- paste0("x", 1:nvars)
                         }
        )
  
  # check suggestions
  if(is.null(suggestions))
    { suggestions <- matrix(nrow = 0, ncol = nvars) }
  else
    { if(is.vector(suggestions)) 
        { if(nvars > 1) suggestions <- matrix(suggestions, nrow = 1)
          else          suggestions <- matrix(suggestions, ncol = 1) }
      else
        { suggestions <- as.matrix(suggestions) }
      if(nvars != ncol(suggestions))
        stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
    }

  # check monitor arg
  if(is.logical(monitor))
    { if(monitor) monitor <- gaMonitor }
  if(is.null(monitor)) monitor <- FALSE
  
  # if optim merge provided and default args for optim()
  if(optim)
    { # merge default and provided parameters
      optimArgs.default <- eval(formals(ga)$optimArgs)
      optimArgs.default$control[names(optimArgs$control)] <- optimArgs$control
      optimArgs$control <- NULL
      optimArgs.default[names(optimArgs)] <- optimArgs
      optimArgs <- optimArgs.default; rm(optimArgs.default)
      if(any(optimArgs$method == c("L-BFGS-B", "Brent")))
        { optimArgs$lower <- lower
          optimArgs$upper <- upper }
      else
        { optimArgs$lower <- -Inf
          optimArgs$upper <- Inf }
      optimArgs$poptim <- min(max(0, optimArgs$poptim), 1)
      optimArgs$pressel <- min(max(0, optimArgs$pressel), 1)
      optimArgs$control$maxit <- as.integer(optimArgs$control$maxit)
      # ensure that optim maximise the fitness
      if(is.null(optimArgs$control$fnscale))
        optimArgs$control$fnscale <- -1
      if(optimArgs$control$fnscale > 0)
        optimArgs$control$fnscale <- -1*optimArgs$control$fnscale
  }

  fitnessSummary <- matrix(as.double(NA), nrow = maxiter, ncol = 6)
  colnames(fitnessSummary) <- names(gaSummary(rnorm(10)))
  bestSol <- if(keepBest) vector(mode = "list", length = maxiter)
             else         list()
  Fitness <- rep(NA, popSize)

  object <- new("ga", 
                call = call, 
                type = type,
                lower = lower, 
                upper = upper, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1,
                maxiter = maxiter,
                suggestions = suggestions,
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = if(is.numeric(pmutation)) pmutation else NA,
                optim = optim,
                fitness = Fitness, 
                summary = fitnessSummary,
                bestSol = bestSol)
                  
  if(maxiter == 0)
    return(object)
  
  # generate beginning population
  Pop <- matrix(as.double(NA), nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if(ng > 0) # use suggestion if provided
  { Pop[1:ng,] <- suggestions }
  
  # fill the rest with a random population
  if(popSize > ng)
    { Pop[(ng+1):popSize,] <- population(object)[1:(popSize-ng),] }
  object@population <- Pop
  
      # If a previous object exists, read it in. - aadams
      start <- 1 # cseiler: start of iteration
      # cseiler: in case the run was interrupted, read the object from your last iteration:
      if(!is.null(previousObject)) {
        print("I will take the most recent object.rds file!")
        object <- readRDS(previousObject)
        fitnessSummary <- object@summary
        Fitness <- object@fitness
        Pop <- object@population
        popSize <- object@popSize
        start <- object@iter + 1
        seed <- 1
        run <- object@maxiter
        } # cseiler

  # start iterations
  # for(iter in seq_len(maxiter))
  for(iter in seq(start, maxiter, 1)) # cseiler: in case of interruption, continue with last iterration
     {
      object@iter <- iter
  
      # Get individuals we need to calculate the fitness for.
      getFitFor <- list()
      for (i in 1:popSize) {
        if (is.na(Fitness[i])) {
          getFitFor <- append(getFitFor, list(Pop[i,]))
        }
      }
        
      # File names.
      argValsFile <- "argValues.txt"
      argTypesFile <- "argTypes.txt"
      argLensFile <- "argLengths.txt"
      argNamesFile <- "argNames.txt"
      names <- names(callArgs)
        
      if (!file.exists(paste0(dataAssimPath, "/",  farmName))) {
          
        # Create the farm and remove the default table.dat file.
        system(sprintf("farm_init.run %s", farmName))
        setwd(farmName)
        Sys.sleep(1)
        system("rm table.dat")
          
        # Edit the SBATCH parameters for the job_script and resubmit_script files.
        # Convert 0-00:10 to the specified jobTime and Your-account-name to def-cseiler-ab.
        # I will keep memory at 4GB
        system(sprintf("sed -i 's|0-00:10|%s|' job_script.sh", jobTime))
        system(sprintf("sed -i 's|Your_account_name|def-cseiler-ab|' job_script.sh"))
        system(sprintf("sed -i 's|Your_account_name|def-cseiler-ab|' resubmit_script.sh"))
          
        Sys.sleep(1)
          
      # If the farm directory already exists, then it means we have previously
      # created an iteration so we need to clean the folder.
      } else {
        system("echo yes | clean.run")
        Sys.sleep(5)
        system("rm table.dat")
      }
      
      # Write argument information to the arg files.
      for (arg in callArgs) {
        cat(as.character(arg), "\n", file = argValsFile, append = TRUE)
        types <- c(typeof(arg), typeof(arg[[1]]))
        cat(types, "\n", file = argTypesFile, append = TRUE)
        write(length(arg), argLensFile, append = TRUE)
      }
      
      write(names(callArgs), argNamesFile)
        
      # File to store all the individuals in.
      popFile <- "population.txt"

      # Check for the existence of the indivs file.
      if (file.exists(popFile)) {system(sprintf("rm %s", popFile))}
        
      Sys.sleep(1)
        
      # Write the individuals to a separate file. The cost function will read in
      # a specific individual to calculate the fitness of using an identifier.
      for (indiv in getFitFor) {
        cat(unlist(indiv), "\n", file = popFile, append = TRUE)
      }
        
      # Create the table.dat file
      dataAssimPath <- "/home/aadams/projects/def-cseiler-ab/aadams/data-assimilation"
      popFile <- paste0(dataAssimPath, "/", farmName, "/", popFile)
      # Each case has its own val file that will be created in setArgs.R.
      valFile <- paste0(dataAssimPath, "/", farmName, "/", argValsFile)
      typeFile <- paste0(dataAssimPath, "/", farmName, "/", argTypesFile)
      lenFile <- paste0(dataAssimPath, "/", farmName, "/", argLensFile)
      nameFile <- paste0(dataAssimPath, "/", farmName, "/", argNamesFile)
        
      argsToChange <- list("parameterFile, run_classic_file")
        
      for (id in 1:length(getFitFor)) {
        # Since each case is a separate instance, I need to add the daisy package to libPaths() each time.
        appendLibPaths <- "R -e \".libPaths(c(.libPaths(), '/home/aadams/daisy/renv'))"
        loadDaisy <- "library(daisy)"
        callCostFunc <- paste0("cost.fun('", id, "', '", popFile, "', '", 
                          valFile, "', '", typeFile, "', '", lenFile, "', '", 
                          nameFile, "', '", dataAssimPath, "', '", 
                          argsToChange, "', '", farmName, "', '", keepOutput, "', '",
                          finalOutputFolder, "')")
        line <- paste0(appendLibPaths, "; ", loadDaisy, "; ", callCostFunc, "\"")
        write(line, "table.dat", append = TRUE)
      }
        
      # Submit meta-jobs.
      system(sprintf("submit.run %d", parallel))
        
      # Wait until the simInfo.txt file exists for every RUN folder.
      Sys.sleep(2)
      runs <- 1:length(getFitFor)
      while (length(runs)) {
        for (runNum in runs) {
          f <- paste0("RUN", runNum)
          # If the simInfo.txt file exists in the run, remove the run's number
          # from the list of runs to check.
          if (file.exists(paste0(f, "/simInfo.txt"))) {
            ind <- match(runNum, runs) # Get index of runNum
            runs <- runs[- ind] # Remove runNum from the list
            # NOTE: doing runs <- runs[- runNum] would remove the value at
            # the runNum-th index, not the value runNum.
          }
        }
        Sys.sleep(10)
      }
        
      # Read in every individual and get their fitness value.
        
      library(readr)

      getFitForInd <- 1
      for (instance in 1:popSize) {
          
        # Info only needs to be read in if the individual's fitness was NA before.
        if (is.na(Fitness[instance])) {
            
          # Move into the RUN folder
          setwd(paste0("RUN", getFitForInd))
            
          # If the output folders are not kept, the score will be on the first line of the
          # simInfo.txt file. Otherwise, it will be on the second line, with the first line
          # being the name of the time-stamped output folder.
          if (keepOutput) {
            score <- read_lines("simInfo.txt", skip = 1, n_max = 1)
          } else {
            score <- readLines("simInfo.txt")
          }
            
          # Return to the farm directory
          setwd(paste0(dataAssimPath, "/", farmName))
            
          Fitness[instance] <- as.numeric(score)
          Pop[instance,] <- unlist(getFitFor[getFitForInd])
          getFitForInd <- getFitForInd + 1
        }
      }
        
      # Move all RUNX folders and population file to the GEN folder
      genFolder <- paste0("GEN", iter)
      # Make the GENX folder.
      system(paste("mkdir", genFolder))
      Sys.sleep(1)
      # Move all RUNX folders to GEN folder.
      for (x in 1:length(getFitFor)) {
        runFolder <- paste0("RUN", x)
        system(paste("mv", runFolder, genFolder))
        Sys.sleep(1)
      }
      
      system(paste("mv", popFile, genFolder))
        
      # if(iter == 81) {print(Pop)} # cseiler Pop OK here
      # if(iter == 81) {print(Fitness)} # cseiler Fitness OK here
      
      # update object
      object@population <- Pop
      object@fitness <- Fitness
      
      # Local search optimisation
      if(optim & (type == "real-valued"))
      {
        if(optimArgs$poptim > runif(1))
        { # perform local search from random selected solution
          # with prob proportional to fitness
          i <- sample(1:popSize, size = 1, 
                      prob = optimProbsel(Fitness, q = optimArgs$pressel))
          # run local search
          opt <- try(suppressWarnings(
                      do.call(stats::optim, 
                              c(list(fn = fitness,
                                     par = Pop[i,],
                                     method = optimArgs$method,
                                     lower = optimArgs$lower,
                                     upper = optimArgs$upper,
                                     control = optimArgs$control), 
                                callArgs))
                      ), silent = TRUE)
          if(is.function(monitor))
            { if(!inherits(opt, "try-error"))
                cat("\b | Local search =", 
                    format(opt$value, digits = getOption("digits")))
              else cat("\b |", opt[1])
              cat("\n")
            }
          if(!inherits(opt, "try-error"))
            { Pop[i,] <- opt$par
              Fitness[i] <- opt$value }
          # update object

      # if(iter == 81) {print(Pop)} # cseiler not being called
      # if(iter == 81) {print(Fitness)} # cseiler not being called

          object@population <- Pop
          object@fitness <- Fitness
          # update iterations summary
          fitnessSummary[iter,] <- gaSummary(object@fitness)
          object@summary <- fitnessSummary
        }
      }
      
      if(keepBest) 
      { 
        object@bestSol[[iter]] <- unique(Pop[Fitness == max(Fitness, na.rm = TRUE),, drop=FALSE]) 
      }

      # apply a user's defined function to update the GA object
      if(is.function(postFitness))
      { 
        object <- do.call(postFitness, c(object, callArgs))
        Fitness <- object@fitness
        Pop <- object@population
      }

      # update iterations summary
      fitnessSummary[iter,] <- gaSummary(object@fitness)
      object@summary <- fitnessSummary

      if(is.function(monitor)) 
        { monitor(object) }
      
      # check stopping criteria
      if(iter > 1)
        object@run <- garun(fitnessSummary[seq(iter),1])
      
      if(object@run >= run) break  
      if(max(Fitness, na.rm = TRUE) >= maxFitness) break
      if(object@iter == maxiter) break  

      ord <- order(Fitness, decreasing = TRUE)
      PopSorted <- Pop[ord,,drop=FALSE]
      FitnessSorted <- Fitness[ord]
        
      # selection
      if(is.function(selection))
        { 
          # set.seed(iter) # cseiler: important, otherwise iterruption has impact on result
          sel <- selection(object)
          # sel <- do.call(selection, c(object, callArgs))
          Pop <- sel$population
          Fitness <- sel$fitness

    #  if(iter == 81) {print(Pop)} # cseiler 
    #  if(iter == 81) {print(Fitness)} # cseiler
        }
      else
        { sel <- sample(1:popSize, size = popSize, replace = TRUE)
          Pop <- object@population[sel,]
          Fitness <- object@fitness[sel]
        }
      object@population <- Pop
      object@fitness <- Fitness
    
      # crossover
      if(is.function(crossover) & pcrossover > 0)
        { nmating <- floor(popSize/2)
          mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
          for(i in seq_len(nmating))
            
             {
                if(pcrossover > runif(1))
                 { parents <- mating[i,]
                   Crossover <- crossover(object, parents)
                   Pop[parents,] <- Crossover$children
                   Fitness[parents] <- Crossover$fitness
                 }
             }             
          object@population <- Pop
          object@fitness <- Fitness
      }

      # mutation
      pm <- if(is.function(pmutation)) pmutation(object) else pmutation
      if(is.function(mutation) & pm > 0)
        { for(i in seq_len(popSize)) 
             {
                if(pm > runif(1)) 
                 { Mutation <- mutation(object, i)
                   Pop[i,] <- Mutation
                   Fitness[i] <- NA
                 }
             }
          object@population <- Pop
          object@fitness <- Fitness
        }

      # elitism
      if(elitism > 0) # (elitism > 0 & iter > 1) 
        { ord <- order(object@fitness, na.last = TRUE)
          u <- which(!duplicated(PopSorted, margin = 1))
          Pop[ord[1:elitism],] <- PopSorted[u[1:elitism],]
          Fitness[ord[1:elitism]] <- FitnessSorted[u[1:elitism]]
          object@population <- Pop
          object@fitness <- Fitness
      } 

      # if(iter == 81) {print(Pop)} # cseiler Pop differs 
      # if(iter == 81) {print(Fitness)} # cseiler Fitness differs

    # cseiler start
    # Sys.sleep(2)
    saveRDS(object, "object.rds") # cseiler
    objectID <- paste("object_iteration_", iter, ".rds", sep = "")
    copyObject <- paste("cp", "object.rds", objectID, sep = " ")
    system(copyObject)
    #  if(iter == 80){stop(iter)}
    # cseiler end
    
    # Move the generation output to the archive folder.
    if (keepOutput) {
      # Move the genFolder and object iteration file to the ARCHIVE folder.
      system(paste0("mv ", genFolder, " ", dataAssimPath, "/", archiveFolder))
      
    } else {
      if (iter > 1) {
        # Remove the previous object iteration file.
        system(paste0("rm ", dataAssimPath, "/", archiveFolder, "/object_iteration_", iter-1, ".rds"))
      }
    }
    
    # Save the most recent object iteration file.
    system(paste0("mv ", "object_iteration_", iter, ".rds ", 
                  dataAssimPath, "/", archiveFolder))
    Sys.sleep(1)
    
    # End of iteration
  }

  # Upon completion of all iterations, copy the GEN folders and object iteration
  # files into the farm directory.
  if (keepOutput) {
    system(paste0("cp ", dataAssimPath, "/", archiveFolder, "/GEN* ", 
                  dataAssimPath, "/", farmName))
  }
  system(paste0("cp ", dataAssimPath, "/", archiveFolder, "/object* ", 
                dataAssimPath, "/", farmName))
  Sys.sleep(1)
  
  # Move the GEN folders and object iteration files into a timestamped folder.
  timeStamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # The name of the folder for the run in the archive folder.
  simArchiveFolder <- paste("CLASSIC", timeStamp, sep = "_")
  system(paste0("mkdir ", archiveFolder, "/", simArchiveFolder))
  Sys.sleep(2)
  
  # Move GEN folders and object iteration files into the simulation archive folder.
  if (keepOutput) {
    system(paste0("mv ", dataAssimPath, "/", archiveFolder, "/GEN* ", 
                  dataAssimPath, "/", archiveFolder, "/", simArchiveFolder))
  }
  system(paste0("mv ", dataAssimPath, "/", archiveFolder, "/object* ", 
                dataAssimPath, "/", archiveFolder, "/", simArchiveFolder))
  # Wait for all GEN folders and object iteration files to be moved.
  while (file.exists(paste0(archiveFolder, "/GEN*")) || file.exists(paste0(archiveFolder, "/object*"))) {
    Sys.sleep(5)
  }
  
  # Remove the final generation folder from the farm directory if output is not being kept.
  if (!keepOutput) {
    system(paste0("rm -rf ", dataAssimPath, "/", farmName, "/GEN*"))
    Sys.sleep(2)
  }
  
  ### cseiler: end of for loop
  # if optim is required perform a local search from the best 
  # solution at the end of GA iterations
  if(optim & (type == "real-valued"))
    { 
      optimArgs$control$maxit <- rev(optimArgs$control$maxit)[1]
      i <- which.max(object@fitness)
      # if not provided suggest approx parscale
      # if(is.null(optimArgs$control$parscale))
      #   optimArgs$control$parscale <- 10^round(log10(abs(object@population[i,])+1))
      # run local search
      opt <- try(suppressWarnings(
                 do.call(stats::optim, 
                         c(list(fn = fitness,
                                par = object@population[i,],
                                method = optimArgs$method,
                                lower = optimArgs$lower,
                                upper = optimArgs$upper,
                                control = optimArgs$control), 
                           callArgs))
                   ), silent = TRUE)
      if(is.function(monitor))
        { if(!inherits(opt, "try-error"))
            cat("\b | Final local search =",
                format(opt$value, digits = getOption("digits")))
          else cat("\b |", opt[1])
        }
      if(!inherits(opt, "try-error"))
        { object@population[i,] <- opt$par
          object@fitness[i] <- opt$value }
  }

  # if(is.function(monitor)) 
  #   { cat("\n"); flush.console() }

  # in case of premature convergence remove NAs from summary 
  # fitness evalutations
  object@summary <- na.exclude(object@summary)
  attr(object@summary, "na.action") <- NULL

  # get solution(s)
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt,,drop=FALSE]
  if(nrow(solution) > 1)
    { # find unique solutions to precision given by default tolerance
      eps <- gaControl("eps")
      solution <- unique(round(solution/eps)*eps, margin = 1)
    }
  colnames(solution) <- parNames(object)
  object@solution <- solution
  if(keepBest)
    object@bestSol <- object@bestSol[!sapply(object@bestSol, is.null)]  
  
  # return an object of class 'ga'
  saveRDS(object, "objectFinal.rds") # cseiler
  return(object)
}

setClassUnion("numericOrNA", members = c("numeric", "logical"))

setClass(Class = "ga", 
         representation(call = "language",
                        type = "character",
                        lower = "numericOrNA", 
                        upper = "numericOrNA", 
                        nBits = "numericOrNA", 
                        names = "character",
                        popSize = "numeric",
                        iter = "numeric", 
                        run = "numeric", 
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = "matrix",
                        elitism = "numeric", 
                        pcrossover = "numeric", 
                        pmutation = "numericOrNA",
                        optim = "logical",
                        fitness = "numericOrNA",
                        summary = "matrix",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix"
                      ),
         package = "GA" 
) 

setMethod("print", "ga", function(x, ...) str(x))

setMethod("show", "ga",
function(object)
 { cat("An object of class \"ga\"\n")
   cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
   cat("Available slots:\n")
   print(slotNames(object))
}) 

summary.ga <- function(object, ...)
{
  nvars <- ncol(object@population)
  varnames <- parNames(object)
  domain <- NULL
  if(object@type == "real-valued")
    { domain <- rbind(object@lower, object@upper)
      rownames(domain) <- c("lower", "upper")
      if(ncol(domain) == nvars) 
         colnames(domain) <- varnames
    }
  suggestions <- NULL
  if(nrow(object@suggestions) > 0) 
    { suggestions <- object@suggestions
      dimnames(suggestions) <- list(1:nrow(suggestions), varnames) 
    }
  
  out <- list(type = object@type,
              popSize = object@popSize,
              maxiter = object@maxiter,
              elitism = object@elitism,
              pcrossover = object@pcrossover,
              pmutation = object@pmutation,
              domain = domain,
              suggestions = suggestions,
              iter = object@iter,
              fitness = object@fitnessValue,
              solution = object@solution)
  class(out) <- "summary.ga"
  return(out)
}

setMethod("summary", "ga", summary.ga)

print.summary.ga <- function(x, digits = getOption("digits"), ...)
{
  dotargs <- list(...)
  if(is.null(dotargs$head)) dotargs$head <- 10
  if(is.null(dotargs$tail)) dotargs$tail <- 2
  if(is.null(dotargs$chead)) dotargs$chead <- 10
  if(is.null(dotargs$ctail)) dotargs$ctail <- 2
  
  cat(cli::rule(left = crayon::bold("Genetic Algorithm"), 
                width = min(getOption("width"),40)), "\n\n")
  # cat("+-----------------------------------+\n")
  # cat("|         Genetic Algorithm         |\n")
  # cat("+-----------------------------------+\n\n")

  cat("GA settings: \n")
  cat(paste("Type                  = ", x$type, "\n"))
  cat(paste("Population size       = ", x$popSize, "\n"))
  cat(paste("Number of generations = ", x$maxiter, "\n"))
  cat(paste("Elitism               = ", x$elitism, "\n"))
  cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
  cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))
  #
  if(x$type == "real-valued")
    { cat(paste("Search domain = \n"))
      do.call(".printShortMatrix", 
              c(list(x$domain, digits = digits), 
                dotargs[c("head", "tail", "chead", "ctail")]))
    }
  #
  if(!is.null(x$suggestions))
    { cat(paste("Suggestions =", "\n"))
      do.call(".printShortMatrix", 
              c(list(x$suggestions, digits = digits), 
                dotargs[c("head", "tail", "chead", "ctail")]))
    }
  #
  cat("\nGA results: \n")
  cat(paste("Iterations             =", format(x$iter, digits = digits), "\n"))
  cat(paste("Fitness function value =", format(x$fitness, digits = digits), "\n"))
  if(nrow(x$solution) > 1) 
    { cat(paste("Solutions = \n")) }
  else
    { cat(paste("Solution = \n")) }
  do.call(".printShortMatrix", 
          c(list(x$solution, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))
  #
  invisible()
}


plot.ga <- function(x, y, ylim, cex.points = 0.7,
                    col = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1)),
                    pch = c(16, 1), lty = c(1,2), legend = TRUE,
                    grid = graphics:::grid, ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  is.final <- !(any(is.na(object@summary[,1])))
  iters <- if(is.final) 1:object@iter else 1:object@maxiter
  summary <- object@summary
  if(missing(ylim)) 
    { ylim <- c(max(apply(summary[,c(2,4)], 2, 
                          function(x) min(range(x, na.rm = TRUE, finite = TRUE)))),
                max(range(summary[,1], na.rm = TRUE, finite = TRUE))) 
  }
  
  plot(iters, summary[,1], type = "n", ylim = ylim, 
       xlab = "Generation", ylab = "Fitness value", ...)
  if(is.final & is.function(grid)) 
    { grid(equilogs=FALSE) }
  points(iters, summary[,1], type = ifelse(is.final, "o", "p"),
         pch = pch[1], lty = lty[1], col = col[1], cex = cex.points)
  points(iters, summary[,2], type = ifelse(is.final, "o", "p"),
         pch = pch[2], lty = lty[2], col = col[2], cex = cex.points)
  if(is.final)
    { polygon(c(iters, rev(iters)), 
              c(summary[,4], rev(summary[,1])), 
              border = FALSE, col = col[3]) }
  else
    { title(paste("Iteration", object@iter), font.main = 1) }
  if(is.final & legend)
    { inc <- !is.na(col)
      legend("bottomright", 
             legend = c("Best", "Mean", "Median")[inc], 
             col = col[inc], pch = c(pch,NA)[inc], 
             lty = c(lty,1)[inc], lwd = c(1,1,10)[inc], 
             pt.cex = c(rep(cex.points,2), 2)[inc], 
             inset = 0.02) }
  
  out <- data.frame(iter = iters, summary)
  invisible(out)
}

setMethod("plot", "ga", plot.ga)
          
setGeneric(name = "parNames", 
           def = function(object, ...) { standardGeneric("parNames") }
          )

setMethod("parNames", "ga",
function(object, ...)
{ 
  names <- object@names
  nvars <- ncol(object@population)
  if(length(names) == 0)
    { names <- paste("x", 1:nvars, sep = "") }
  return(names)
})

gaSummary <- function(x, ...)
{
  # compute summary for each step
  x <- na.exclude(as.vector(x))
  q <- fivenum(x)
  c(max = q[5], mean = mean(x), q3 = q[4], median = q[3], q1 = q[2], min = q[1])
}
