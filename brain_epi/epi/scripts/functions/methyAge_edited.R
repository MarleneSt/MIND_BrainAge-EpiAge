methyAge_edited <- function(betas, clock='HorvathS2013', age_info=NA, fit_method='Linear', 
                     do_plot=TRUE, inputation=TRUE, simple_mode=FALSE, 
                     species='Homo sapiens', MM_array=FALSE, use_cores=detectCores()){
  ## prepare clock coefficients
  usable_clocks <- suppressMessages(availableClock())
  if (!(clock %in% usable_clocks)){
    stop(message(paste0("Unavailable for the input clock: ", clock,
                        ". Please choose one of the available clocks:\n\n", paste(usable_clocks, collapse=", "), '\n')))
  } else {
    if(grepl('^PC[A-Z]', clock)){
      is_PCclock <- TRUE
      data(list='PC-clocks', envir=environment())
      coefs$Coefficient <- coefs[[clock]]
    }else{
      data(list=clock, envir=environment())
    }
    
    is_beta <- TRUE
    plot_simple <- FALSE
    x_lim = y_lim = c(0, 100)
    if (clock == 'HorvathS2013' & !simple_mode){
      betas <- horvathPreprocess(betas, normalizeData=TRUE, use_cores=use_cores)
    } else if (clock == 'ZhangQ2019'){
      betas <- preprocessZhangQ2019(betas)
    } else if (clock == 'YangZ2016'){
      mAge <- EstEpiTOC(betas, epiTOCcpgs, mode="raw", ref.idx=NULL)
      m_age <- as.matrix(mAge)
      is_beta <- FALSE
    } else if(clock == 'epiTOC2'){
      m_age <- epiTOC2(betas, coefs, full_model=!simple_mode)
      m_age <- as.matrix(m_age)
      is_beta <- FALSE
    } else if(clock == 'DunedinPACE'){
      betas <- preprocessDunedinPACE(betas, ref_means=gold_standard_means)
      y_lim <- c(0.5, 2)
    } else if(clock == 'BernabeuE2023c'){
      coefs$Probe <- sub('_2', '', coefs$Probe_2)  
    } else if(clock %in% c('LuA2023p1', 'LuA2023p2', 'LuA2023p3')){
      if(!MM_array){
        betas <- arrayConverter_EPICtoMM(betas, coefs$Probe[-1])
        inputation <- FALSE
      }
    }
    ## Free the Y limits in plotting
    if(clock %in% c('YangZ2016', 'ZhangY2017', 'LuA2019')){
      y_lim = c(NULL, NULL)
      plot_simple <- TRUE
    }
    
  }  
  
  if(is_beta){
    r_coefs <- coefs
    #data('HorvathS2013')
    coefs <- setNames(coefs$Coefficient, coefs$Probe)
    ## add intercept
    betas <- rbind(betas, Intercept=1)
    
    ## identify missing probes, set their beta values as 0.5
    betas <- betas[rownames(betas) %in% names(coefs), ]
    missing_probe <- setdiff(names(coefs), rownames(betas))
    if(length(missing_probe) > 0){
      warning(paste(c("Found ", length(missing_probe), "out of", length(coefs),
                      "probes missing! They will be assigned with mean values from reference dataset, missing probes are:\n ", missing_probe), collapse=" "))
    }
    if (inputation){
      ## Mean imputation
      data(list='golden_ref', envir=environment())
      ref_mean <- setNames(golden_ref$Mean, rownames(golden_ref))
      ref_mean <- ref_mean[names(ref_mean) %in% names(coefs)]
      betas <- meanImputation(mt=betas, ref=ref_mean, only_ref_rows=FALSE)
    }else{
      betas[is.na(betas)] <- 0
    }
    
    if(clock %in% c('BernabeuE2023c')){
      ## Some CpGs are used as quadratic 
      probes_2 <- r_coefs$Probe[grepl('_2', r_coefs$Probe_2)]
      betas_2 <- betas[rownames(betas) %in% probes_2, ]^2
      rownames(betas_2) <- paste0(rownames(betas_2), '_2')
      betas <- betas[rownames(betas) %in% r_coefs$Probe_2, ]
      betas <- rbind(betas, betas_2)
      coefs <- setNames(r_coefs$Coefficient, r_coefs$Probe_2) 
      coefs_L <- setNames(r_coefs$Coefficient_L, r_coefs$Probe_2)
    }
    
    ## matrix multiplication
    m_age <- t(betas) %*% matrix(data=coefs[rownames(betas)])
    
    ## post transformation
    if(clock %in% c('HorvathS2013', 'ShirebyG2020', 'HorvathS2018', 
                    'McEwenL2019', 'PCHorvathS2013', 'PCHorvathS2018')){
      HorvathS2013_transform <- function(x){
        if (x > 0){
          x <- x *(20 + 1) + 20
        }else{
          x <- exp(x + log(20 + 1)) - 1
        }
      }
      m_age[,1] <- sapply(m_age[,1], HorvathS2013_transform)
    } else if(clock %in% c('CBL_specific', 'CBL_common', 'Cortex_common')){
      m_age[,1] <- exp(m_age[,1])
    } else if(clock %in% c('BernabeuE2023c')){
      ## log predictor for samples younger than 20
      less_20 <- m_age[,1] < 20
      # if(sum(less_20) > 0){
      #   # m_age_L <- t(betas[less_20, ]) %*% matrix(data=coefs_L[rownames(m_age)[less_20]]) # VB commented this out as caused error with dimensions
      #   m_age_L <- t(betas[, less_20]) %*% matrix(data=coefs_L) # VB replaced above line with this line
      #   m_age[less_20, 1]  <- exp(m_age_L[, 1]) 
      # }
      if(sum(less_20) > 0){ # version 2
        sample_names <- rownames(m_age)[less_20]  # the sample IDs (should match colnames in betas)
        betas_L <- betas[, sample_names, drop = FALSE]  # get DNAm data for those samples
        
        # coefs_L must be ordered the same as the rows (CpGs) in betas
        coefs_vec <- coefs_L[rownames(betas)]  # ensure alignment
        coefs_mat <- matrix(data = coefs_vec, nrow = length(coefs_vec))  # convert to column matrix
        
        m_age_L <- t(betas_L) %*% coefs_mat  # matrix multiplication
        
        m_age[less_20, 1] <- exp(m_age_L[, 1])  # back-transform from log
      }
      
    } else if(clock %in% c('LuA2023p1')){
      m_age[,1] <- exp(m_age[,1]) - 2
    } else if(clock %in% c('LuA2023p2')){
      data(list='SpeciesAgeInfo', envir=environment())
      my_max <- ifelse(species %in% c('Homo sapiens', 'Mus musculus'), 1, 1.3)
      y_maxAge <- SpeciesAgeInfo[species, 'maxAge'] * my_max
      y_gestation <- SpeciesAgeInfo[species, 'GestationTimeInYears']
      m_age[,1] <- exp(-exp(-1 * m_age[,1]))
      m_age[,1] <- m_age[,1]*(y_maxAge + y_gestation) - y_gestation
    } else if(clock %in% c('LuA2023p3')){
      data(list='SpeciesAgeInfo', envir=environment())
      y_gestation <- SpeciesAgeInfo[species, 'GestationTimeInYears']
      y_ASM <- SpeciesAgeInfo[species, 'averagedMaturity.yrs']
      F2_revtrsf_clock3 <- function(y.pred, m1){
        ifelse(y.pred<0, (exp(y.pred)-1)*m1 + m1, y.pred*m1+m1)  #wyc
      }
      a_Logli <- 5 * (y_gestation / y_ASM)^0.38
      r_adult_age <- F2_revtrsf_clock3(m_age[,1], a_Logli)
      m_age[,1] <- r_adult_age *(y_ASM + y_gestation) - y_gestation
    }
    
  }
  ## save results
  m_age <- data.frame(Sample=rownames(m_age), mAge=as.numeric(m_age))
  
  ## calculate age acceleration
  warning_message <- "\n'age_info' should be a dataframe which contains sample ID and age information, like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge acceleration will not be calculated."
  if (class(age_info) == "data.frame") {
    if (all(c('Sample', 'Age') %in% colnames(age_info))){
      m_age <- merge(age_info, m_age, by='Sample', sort=FALSE)
      if (nrow(m_age) < 1){
        stop(message("Colnames of the input beta dataframe do not match any of the values of the 'Sample' column in age_info!"))
      }
      if(clock == 'PCGrimAge'){
        if('Sex' %in% colnames(age_info)){
          #m_age$is_Female <- NA
          m_age$is_Female <- gsub('^F.*', 1, m_age$Sex)
          m_age$is_Female <- gsub('^M.*', 0, m_age$is_Female)
          m_age$is_Female <- as.numeric(m_age$is_Female)
          
          m_age$mAge <- m_age$mAge + as.matrix(m_age[, c('is_Female', 'Age')]) %*% PCGrimAge_agesex$PCGrimAge
        }else{
          stop(message("\nTo calculate 'PCGrimage', 'age_info' should include a 'Sex' column that contains binary sex annotation, i.e. either Female or Male."))
        }
      }
      if("Color" %in% colnames(age_info)){
        point_color <- m_age$Color
      }else{
        point_color <- NA
      }
      if("Shape" %in% colnames(age_info)){
        point_pch <- m_age$Shape
      }else{
        point_pch <- NA
      }
      m_age$Age_Acceleration <- NA
      m_age$Age_Acceleration[!is.na(m_age$Age)] <- getAccel(m_age$Age[!is.na(m_age$Age)], m_age$mAge[!is.na(m_age$Age)],, method=fit_method, title=clock, do_plot=do_plot, point_color=point_color, point_shape=point_pch, simple=plot_simple, x_lim=x_lim, y_lim=y_lim)
    }else{
      warning(message("\nThe colnames of age_info should include both 'Sample' and 'Age', like:\nSample\tAge\nname1\t30\nname2\t60\nname3\t40\nAge\nAge acceleration will not be calculated."))
    }
  } else if (is.na(age_info[1])){
    if(clock == 'PCGrimAge'){
      stop(message("\nTo calculate 'PCGrimage': \n'age_info' should be a dataframe which contains sample ID, age, sex information, like:\nSample\tAge\tSex\nname1\t30\tFemale\nname2\t60\tMale\nname3\t40\tFemale\n"))
    }
  } else {
    warning(message(warning_message))
  }
  return(m_age)
}