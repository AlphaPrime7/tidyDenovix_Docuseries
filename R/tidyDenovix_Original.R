tidyDenovix = function(dfile, file_type= c('csv','excel','txt'), sample_type = c('RNA','DNA'), check_level = c('strict','lax'), qc_omit = NULL, normalized = c('yes','no'), fun = NA){

  suppressWarnings({

    #PRE-TRANSPOSE:workflow before transposing

    nofun <- is.na(fun)

    xdf_qc = read_denovix_data(dfile=dfile, file_type = file_type)
    xdf = read_denovix_data(dfile=dfile, file_type = file_type)

    #get the qc df
    qc_df = qc_attributes(dfile = dfile,xdf)

    #get the lambda df
    lambda_df = extract_wavelength(xdf)
    #lambda_df = make_wavelength()
    lambda_df_special = lambda_df[complete.cases(lambda_df), ]
    #lambda_df_special = make_wavelength
    #lambda_df_special = lambda_df[complete.cases(lambda_df), ]

    #initial col_names for qc crosscheck
    xdfc = janitor::clean_names(xdf)
    esn = extract_sample_names(dfile=dfile)
    sample_names = xdfc[, c("sample_name")]
    sample_names_wl = append('wave_length', sample_names)

    #TRANSPOSE the df-for workflow after transposing

    xdf = data.table::transpose(l = xdf)

    #sample_col_names<- vector("list")
    sample_col_names<- c()
    for(j in xdf[2,]){
      if(is.na(j) != nofun){
        sample_col_names <- c(sample_col_names,j)
      }
    }
    sample_col_names = append('wave_length', sample_col_names)

    #column names quality control
    qc_colnames = as.numeric((sample_names_wl == sample_col_names))

    if (sum(qc_colnames) == length(sample_names_wl)){
      col_names = sample_names
    } else {
      col_names = c(1:length(sample_names))
    }

    #wrangle df
    cutoff = which( colnames(xdf)=="Exposure" )
    aftercut = as.numeric(cutoff + 1)
    xdf = xdf[aftercut:nrow(xdf),]

    #bind qc_df with wrangled df for QC
    xdf =rbind(qc_df, xdf)
    colnames(xdf) = esn

    #QC now
    xdf = lambda_check(xdf,sample_type = sample_type, check_level = check_level)
    samp_names = lambda_check_source(xdf_qc,sample_type = sample_type, check_level = check_level)
    samp_names_wl = append('wave_length', samp_names)

    if( is.null(qc_omit) || qc_omit == 'yes'){

      if(is.null(normalized) || 'no' %in% normalized){
        xdf = xdf[6:nrow(xdf),]
        xdf = cbind(lambda_df, xdf)
        xdf = xdf[complete.cases(xdf), ]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        #xdf = xdf %>% drop_na()
        #xdf = na.omit(xdf)
        rownames(xdf) = c(1:nrow(xdf))
        colnames(xdf) = samp_names_wl

        return(xdf)

      } else {
        xdf = xdf[6:nrow(xdf),]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        xdf = xdf[complete.cases(xdf), ]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], min_max_norm))
        xdf = cbind(lambda_df_special, xdf)
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        #xdf = xdf %>% drop_na()
        #xdf = na.omit(xdf)
        #rownames(xdf) = c(1:nrow(xdf))
        colnames(xdf) = samp_names_wl

        return(xdf)

      }

    } else {

      if(is.null(normalized) || 'no' %in% normalized){
        n <- 'qc_df'
        assign( paste0(n, '_check'), as.data.frame(xdf), envir = parent.frame())
        xdf = xdf[6:nrow(xdf),]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        xdf = cbind(lambda_df, xdf)
        xdf = xdf[complete.cases(xdf), ]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        #xdf = xdf %>% drop_na()
        #xdf = na.omit(xdf)
        rownames(xdf) = c(1:nrow(xdf))
        colnames(xdf) = samp_names_wl

        return(xdf)

      } else {
        n <- 'qc_df'
        assign( paste0(n, '_check'), as.data.frame(xdf), envir = parent.frame())
        xdf = xdf[6:nrow(xdf),]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        xdf = xdf[complete.cases(xdf), ]
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], min_max_norm))
        xdf = cbind(lambda_df_special, xdf)
        xdf = as.data.frame(lapply(xdf[1:ncol(xdf)], as.numeric))
        #xdf = xdf %>% drop_na()
        #xdf = na.omit(xdf)
        #rownames(xdf) = c(1:nrow(xdf))
        colnames(xdf) = samp_names_wl

        return(xdf)

      }

    }

  })

}
