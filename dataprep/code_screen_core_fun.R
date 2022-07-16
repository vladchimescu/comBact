##==========================SCREEN FUNCTIONS==========================
##-------to extract plates--------
Ext_384plates_Eli <- function(file_path, output_dir, barcode = TRUE){
  
  header <- c("Time", "Barcode",  "PlateID", "Family", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20", "A21", "A22", "A23", "A24", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18", "B19", "B20", "B21", "B22", "B23", "B24", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12", "D13", "D14", "D15", "D16", "D17", "D18", "D19", "D20", "D21", "D22", "D23", "D24", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", "E20", "E21", "E22", "E23", "E24", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19", "F20", "F21", "F22", "F23", "F24", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12", "G13", "G14", "G15", "G16", "G17", "G18", "G19", "G20", "G21", "G22", "G23", "G24", "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", "H20", "H21", "H22", "H23", "H24", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10", "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", "I20", "I21", "I22", "I23", "I24", "J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8", "J9", "J10", "J11", "J12", "J13", "J14", "J15", "J16", "J17", "J18", "J19", "J20", "J21", "J22", "J23", "J24", "K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10", "K11", "K12", "K13", "K14", "K15", "K16", "K17", "K18", "K19", "K20", "K21", "K22", "K23", "K24", "L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L20", "L21", "L22", "L23", "L24", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "M20", "M21", "M22", "M23", "M24", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19", "N20", "N21", "N22", "N23", "N24", "O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8", "O9", "O10", "O11", "O12", "O13", "O14", "O15", "O16", "O17", "O18", "O19", "O20", "O21", "O22", "O23", "O24", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P20", "P21", "P22", "P23", "P24")
  dataset <- read.table(file_path, header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, skip = 1, col.names=header)
  
  plates_vec <- unique(dataset$PlateID)
  #Convert Date and time columns to time-intervals in hours and write the files
  for (plate in plates_vec){
    plate_data <- subset(dataset, PlateID == plate)
    Date_time <- as.character(as.vector(plate_data[[1]]))
    dtm <- strptime(Date_time, format = "%m/%d/%Y %H:%M:%S", tz = "CET")
    dtm[is.na(dtm)] = strptime(Date_time[is.na(dtm)], format = "%m/%d/%Y %H:%M", tz = "CET")
    
    Time_h = vector(mode = "numeric", length = length(dtm))  
    for (i in 2:length(Time_h))
    {
      j = i-1
      diff = dtm[i]-dtm[j]
      Time_h[i] = Time_h[j] + as.numeric(diff, units="hours")
    }
    
    Abs = as.matrix(plate_data[,5:length(plate_data)])
    new_table = cbind(Time_h, Abs)
    PlateID = as.vector(plate_data[[3]])
    Barcode = as.vector(plate_data[[2]])
    if(barcode == TRUE)
    {out_file = paste0(output_dir,Barcode [1],".txt",collapse=NULL)} else {
      out_file = paste0(output_dir,PlateID [1],".txt",collapse=NULL)
    }
    write.table(new_table,file = out_file, append = FALSE, sep = "\t",row.names = FALSE)
    
  }
  
} 

#Ext_384plates_Eli("/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/screenK1_biomek_files/to_extract/run_K_13.txt", "/Users/Elisabetta/Documents/EMBL/combinatorials/screen_K/screenK1_extract_reads/test/", barcode = FALSE) 


##these work for list of mapped plates

#add time points
add_tps <- function (df) {
  v <- c(1:8)  
  df$Time_points <- rep (v)
  return(df)}

#annotated <- lapply (annotated, add_tps)

#replace function 
torep = ","
toberep = "."

repl <-  function (df, torep, toberep) {
  
  df[, c("Concentration", "Flag_Conc", "Drug_Flag_Conc", "Full.name")] <- 
    lapply(df[, c("Concentration", "Flag_Conc", "Drug_Flag_Conc", "Full.name")], 
           function(x) {
             x <- gsub(x, pattern = torep, replacement= toberep)
             return(as.numeric(x))
           })
  return(df)
  
}



#=======do bg subtraction keeping the list=======

#just first tp
bg_sub <- function (df) {
  
  #subtract 1st time point
  firsttp <- subset (df,  Time_points == 1) %>% dplyr::group_by (Well) %>%
    dplyr::summarise(mean_1sttp = mean (OD))
  
  test <- dplyr::inner_join (df, firsttp) %>%
    plyr::mutate(test, OD_no1st = OD - mean_1sttp)
  
  return(test)
}


#first tp or any minimum if lower than 1st (to avoid negative values)
bg_sub_mod <- function (df) {
  #get 1st time point value
  firsttp <- subset (df,  Time_points == 1) %>% 
    dplyr::group_by (Well) %>%
    dplyr::summarise(mean_1sttp = mean (OD))
  
  test <- dplyr::inner_join (df, firsttp)
  final <- dplyr::rename(test, OD_to_subtract = mean_1sttp) 
  
  # add this for plates that have missing values for the first time point
  bg_plate = median(final$OD_to_subtract, na.rm = T)
  final = dplyr::mutate(final, 
                        OD_to_subtract = ifelse(!is.na(OD_to_subtract),
                                                OD_to_subtract,
                                                bg_plate))
  final = dplyr::mutate(final, OD_no1st = OD - OD_to_subtract)
  final = dplyr::mutate(final, OD_no1st = ifelse(OD_no1st > 0, OD_no1st, 0))
  
  return(final)
}

bgsub_estimate <- function(df) {
  df_out = df
  
  non.spiking = dplyr::group_by(df, Well) %>%
    dplyr::filter(OD[Time_points == 2] > OD[Time_points == 1]) %>%
    dplyr::filter(OD[Time_points == 3] > OD[Time_points == 1]) %>%
    dplyr::filter(OD[Time_points == 4] > OD[Time_points == 1]) 
  
  if(length(unique(non.spiking$Well)) > 5) {
    # median of non-spiking wells as background estimate for spiking wells
    bg_df = dplyr::ungroup(non.spiking) %>%
      dplyr::filter(Time_points == 1)
    bg_est = median(bg_df$OD)
    
    spiking = dplyr::anti_join(df, non.spiking, by="Well")
    if(!is.na(bg_est)) {
      spiking = dplyr::mutate(spiking, OD_no1st = ifelse(OD > bg_est, OD - bg_est, 0))
      
      df_out = dplyr::bind_rows(ungroup(non.spiking),
                                spiking)
    }
    
  }
  df_out
  
}

# smooth the spike in OD between time points 4 and 5
smooth_OD <- function(df) {
  # wells to smooth are those in which tp4 > tp5
  to_smooth = dplyr::group_by(df, Well) %>%
    filter(OD[Time_points == 5] < OD[Time_points == 4]) %>%
    distinct(Well)
  # for these wells to be smoothed: tp4 = tp5
  df = group_by(df, Well) %>%
    mutate(OD = ifelse(Well %in% to_smooth$Well & Time_points == 4, OD[Time_points == 5], OD))
  df
}


#this subtracts whatever is lower among first 3 time points
bg_sub_among3tp <- function (df) {
  #get 1st time point value
  firsttp <- subset (df,  Time_points == 1) %>% dplyr::group_by (Well) %>%
    dplyr::summarise(mean_1sttp = mean (OD))
  
  test <- dplyr::inner_join (df, firsttp)
  
  #get 2nd time point value
  sectp <- subset (df,  Time_points == 2) %>% dplyr::group_by (Well) %>%
    dplyr::summarise(mean_2ndtp = mean (OD))
  
  test2 <- dplyr::inner_join (test, sectp)
  
  #get 3rd time point value
  thirdtp <- subset (df,  Time_points == 3) %>% dplyr::group_by (Well) %>%
    dplyr::summarise(mean_3rdtp = mean (OD))
  
  test3 <- dplyr::inner_join (test2, thirdtp) %>%
    plyr::mutate(test3, OD_to_subtract = pmin(mean_1sttp,mean_2ndtp,mean_3rdtp)) %>%
    plyr::mutate(test3, OD_no1st = OD - OD_to_subtract)
  
  #test3$OD_no1st = OD - OD_to_subtract
  return(test3)
}

#this does key stats on each plate
stats_whole <- function(df, output_dir) {
  all_stats <- dplyr::group_by (df, Time_points, plate_name, run, donor_conc, plate_tp) %>%
    dplyr::summarise (med_plate = median(OD_no1st), sd=sd(OD_no1st), se = sd(OD_no1st)/sqrt(length(OD_no1st)),
                      mean_plate = mean(OD_no1st), min_plate = min(OD_no1st), max_plate = max(OD_no1st))
  
  write.table(all_stats, file = paste0(output_dir, "all_stats_screen.tsv", collapse=NULL), sep = "\t",row.names = FALSE)
  return(all_stats)
  
}

#this extracts the bug wells and does stats on them
extract_bug <- function(df) {
  bug <- df[df$Drug == "BUG", ]
  bug <- dplyr::group_by (bug, Drug, Time_points, plate_name, run, donor_conc, plate_tp) %>%
    dplyr::summarise (med_plate = median(OD_no1st), sd=sd(OD_no1st), se = sd(OD_no1st)/sqrt(length(OD_no1st)),
                      mean_plate = mean(OD_no1st), min_plate = min(OD_no1st), max_plate = max(OD_no1st))
  
  return(bug)
}


###-----function to annotates plates NOT TESTED YET--------
read_annotate <- function (read_dir, map_dir, output_dir) {
  
  read <- read.table(read_dir, header= T, sep="\t")
  out <- melt (read, id.vars = "Time_h", variable.name = "wells", value.name= "OD")
  
  #fix names
  colnames = c("Time_h","Well","OD")
  colnames(out) = colnames
  
  #load map
  map <- read.table (map_dir, header = TRUE, sep = "\t")
  
  #===fix the map for subsequent plotting====
  #add leading 0 before digits 1-9 in Conc_Flag (otherwise it thinks 2 > 10)
  map$flag_fixed <- sprintf("%02d", as.numeric(map$Conc.flag))
  #put the fixed column back in the map
  map <- cbind(map, flag_fixed)
  #remove old flag
  map$Conc_flag = NULL
  #round decimals in concentrations to avoid long labels
  map$Concentration <- round(map$Concentration, digits = 4)
  #create columns "Flag + real conc" and "Drug_flag_conc" with this digit fixed
  map$Flag_Conc <- paste0(map$flag_fixed,"_",map$Concentration)
  map$Drug_Flag_Conc <- paste0(map$Drug,"_",map$Flag_Conc)
  
  #==============match reads with the map=============
  annotated_read <- inner_join (out, map, by = "Well")
  
  #add time points
  v <- c(1:25)
  annotated_read$Time_points <- rep (v)
  
  write.csv(annotated_read, paste0(output_dir, file_path_sans_ext(read_dir), ".csv"), row.names=FALSE)
  return (annotated_read)
  
}


#function to plot all curves of a plate at once (TO TEST)

#plot all curves by well
#annotated read has to be the unlisted version of the dataset
plot_curves <- function (annotated_read, plot_dir) {
  
  title_plot = as.character(unique(annotated_read$plate_name))
  pdf (paste0(plot_dir, title_plot, "_plate_curves.pdf"), width=13, height=10) #onefile=TRUE
  pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)
  #if you want 96 panels per page (12x8: good format for titration plates) use the following
  #pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)
  
  
  for (i in seq_len(pages)){
    page <- ggplot(annotated_read, aes (Time_h, OD, color=Replicate)) +
      geom_point(shape = 1, size = 0.5, show.legend = F) +
      #here the no of cols and rows suits well the titration format: 12 concentrations x 8 drugs per page)
      #facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 12, nrow = 8, page = i) +
      facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 9, nrow = 6, page = i) +
      theme(axis.text.x = element_text(size=6, angle=60, hjust=1)) +
      theme(axis.text.y = element_text(size=6)) +
      theme(strip.text.x = element_text(size=6)) +
      labs(title=title_plot, x="Time", y="OD") +
      theme(plot.title = element_text(hjust = 0.5,size=18)) +
      scale_y_continuous(limits=c(0.0, 0.8)) +
      #scale_y_continuous(limits=c(min(annotated_read$OD), max(annotated_read$OD))) +
      theme(panel.spacing = unit(1, "lines"))
    print(page)
    
  }
  dev.off()
  
}


plot_curves_titration <- function (df, plot_dir) {
  
  title_plot = as.character(unique(df$plate_name))
  pdf(paste0(plot_dir, title_plot, "_plate_curves.pdf"), width=13, height=10) #onefile=TRUE
  #pages = ceiling(length(unique(df$Drug)) / 16)
  #if you want 96 panels per page (12x8: good format for titration plates) use the following
  #pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)
  
  p <- ggplot(subset(df, Time_points == 16), aes(x= OD_minsub, y=as.factor(Concentration), colour = Replicate)) +
    geom_point(shape = 1, size = 0.5, show.legend = F) +
    geom_path(aes(group = Replicate)) +
    theme(axis.text.x = element_text(size=6, angle=60, hjust=1)) +
    theme(axis.text.y = element_text(size=6)) +
    theme(strip.text.x = element_text(size=6)) +
    #here the no of cols and rows suits well the titration format: 12 concentrations x 8 drugs per page)
    #facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 12, nrow = 8, page = i) +
    facet_rep_wrap( ~ Drug, nrow = 4, ncol = 4, scales = 'free_x', repeat.tick.labels = 'all') +
    labs(title=title_plot, x="OD at 6th tp", y="Concentration") +
    theme(plot.title = element_text(hjust = 0.5,size=18)) +
    coord_flip() +
    #ylim(c(min(df$OD_minsub), max(df$OD_minsub))) +
    theme(panel.spacing = unit(1, "lines")) +
    #scale_y_continuous(trans = "reverse") +
    # scale_x_continuous(breaks = pretty(inter_scores$fit_donor, n = 10)) +
    # scale_y_continuous(breaks = pretty(inter_scores$fit_comb, n = 10)) +
    theme(plot.title = element_text(hjust = 0.5, size =28))
  
  print(p)
  
  dev.off()
  
}



plot_curves_titration_AUC <- function (df, plot_dir) {
  
  title_plot = as.character(unique(df$plate_name))
  pdf(paste0(plot_dir, title_plot, "_plate_curves_AUC.pdf"), width=13, height=10) #onefile=TRUE
  #pages = ceiling(length(unique(df$Drug)) / 16)
  #if you want 96 panels per page (12x8: good format for titration plates) use the following
  #pages = ceiling(length(unique(annotated_read$Drug_Flag_Conc)) / 54)
  
  p <- ggplot(subset(df, Time_points == 16), aes(x= AUC_well_ODno1st, y=as.factor(Concentration), colour = Replicate)) +
    geom_point(shape = 1, size = 0.5, show.legend = F) +
    geom_path(aes(group = Replicate)) +
    theme(axis.text.x = element_text(size=6, angle=60, hjust=1)) +
    theme(axis.text.y = element_text(size=6)) +
    theme(strip.text.x = element_text(size=6)) +
    #here the no of cols and rows suits well the titration format: 12 concentrations x 8 drugs per page)
    #facet_wrap_paginate(Drug ~ Flag_Conc, ncol= 12, nrow = 8, page = i) +
    facet_rep_wrap( ~ Drug, nrow = 4, ncol = 4, scales = 'free_x', repeat.tick.labels = 'all') +
    labs(title=title_plot, x="AUC", y="Concentration") +
    theme(plot.title = element_text(hjust = 0.5,size=18)) +
    coord_flip() +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.title = element_text(hjust = 0.5, size =28))
  
  print(p)
  
  dev.off()
  
}


##-----ggplot extension----

# for printing the geom_smooth statistics 
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 3), 
                                                  b = format(coef(m)[2], digits = 3), 
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)



##================================NORMALISATION FUNCTIONS===============================

#plate based normalisation
#Z score: subtracts the mean and divides by sd
Zscore <- function (df){
  plate_mean <- dplyr::group_by(df, Time_points) %>%
    dplyr::summarize(mean_plate = mean(OD_no1st), med_plate = median(OD_no1st), sd_plate = sd(OD_no1st), se = sd(OD_no1st)/sqrt(length(OD_no1st)),
                     min_plate= min(OD_no1st), max_plate = max (OD_no1st))
  df = inner_join(df, plate_mean, by = "Time_points") %>% mutate(z_score = (OD_no1st - mean_plate) / sd_plate)
  return(df)                                                               
  
}

#list_all <- lapply (list_all, Zscore)

#control based normalisation
#POC (% of controls): I can use as positive control (unperturbed measurements) BUG wells (just 6) or more robust measurement (e.g. top 5% of plate = 20 wells)
#this might be refined --> for some plates overkilling is widespread so 20 wells would shift the correction way down

#POC with bugs
POC_bugwells <- function (df) {
  plate_mean_bug <- dplyr::group_by(df, Time_points) %>%
    dplyr::filter(Drug == "BUG") %>%
    dplyr::summarize(mean_plate_bug = mean(OD_no1st), 
                     med_plate_bug = median(OD_no1st),
                     sd_plate_bug = sd(OD_no1st), 
                     se_bug = sd(OD_no1st)/sqrt(length(OD_no1st)),
                     min_plate_bug= min(OD_no1st),
                     max_plate_bug = max (OD_no1st))
  
  df = inner_join(df, plate_mean_bug, by = "Time_points") %>%
    mutate(POC_bug = ((OD_no1st/mean_plate_bug)*100))
  return(df)
}

# Control-well robust mean normalization of AUC_well_OD_no1st
POC_bugwells_AUC <- function (df) {
  # first check if there are at least 3 BUG wells
  if(length(unique(dplyr::filter(df, Drug == "BUG")$Well)) >= 3) {
    plate_mean_bug <- dplyr::filter(df, Drug == "BUG") %>%
      dplyr::group_by(Time_points) %>%
      dplyr::summarise(mean_AUC_bug = mean(AUC_well_ODno1st),
                       med_AUC_bug = median(AUC_well_ODno1st),
                       mean_AUC_bug_rob = robustMean(AUC_well_ODno1st))
  } else {
    # if not, use 8 wells with highest fitness in the plate
    plate_mean_bug <- dplyr::filter(df, AUC_well_ODno1st >= quantile(AUC_well_ODno1st, 0.98)) %>%
      dplyr::group_by(Time_points) %>%
      dplyr::summarise(mean_AUC_bug = mean(AUC_well_ODno1st),
                       med_AUC_bug = median(AUC_well_ODno1st),
                       mean_AUC_bug_rob = robustMean(AUC_well_ODno1st))
  }

  df = inner_join(df, plate_mean_bug, by = "Time_points") %>% 
    mutate(AUCnorm_bug_mean = AUC_well_ODno1st/mean_AUC_bug,
           AUCnorm_bug_rob = AUC_well_ODno1st/mean_AUC_bug_rob,
           AUCnorm_bug_med = AUC_well_ODno1st/med_AUC_bug)
  return(df)
}


# Control-well robust mean normalization of OD_well_OD_no1st
POC_bugwells_OD <- function (df) {
  # first check if there are at least 3 BUG wells
  if(length(unique(dplyr::filter(df, Drug == "BUG")$Well)) >= 3) {
    plate_mean_bug <- dplyr::filter(df, Drug == "BUG") %>%
      group_by(Time_points) %>%
      dplyr::summarise(mean_OD_bug = mean(OD_no1st),
                       med_OD_bug = median(OD_no1st),
                       mean_OD_bug_rob = robustMean(OD_no1st))
  } else {
    # if not, use 8 wells with highest fitness in the plate
    plate_mean_bug <- dplyr::filter(df, OD_no1st >= quantile(OD_no1st, 0.98)) %>%
      group_by(Time_points) %>%
      dplyr::summarise(mean_OD_bug = mean(OD_no1st),
                       med_OD_bug = median(OD_no1st),
                       mean_OD_bug_rob = robustMean(OD_no1st))
  }
  
  df = inner_join(df, plate_mean_bug) %>% 
    mutate(ODnorm_bug_mean = OD_no1st/mean_OD_bug,
           ODnorm_bug_rob = OD_no1st/mean_OD_bug_rob,
           ODnorm_bug_med = OD_no1st/med_OD_bug)
  return(df)
}


#POC with bugs (customised for titration, so if no BUG wells it picks the null concentration wells - so no drug in there)
POC_bugwells_mod <- function (df) {
  if ( "BUG" %in% unique(df$Drug)) {
    plate_mean_bug <- dplyr::group_by(df, Time_points) %>%
      dplyr::filter(Drug == "BUG") %>%
      dplyr::summarize(mean_plate_bug = mean(OD_no1st), med_plate_bug = median(OD_no1st), sd_plate_bug = sd(OD_no1st), se_bug = sd(OD_no1st)/sqrt(length(OD_no1st)),
                       min_plate_bug= min(OD_no1st), max_plate_bug = max (OD_no1st))
    
    df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(POC_bug = ((OD_no1st/mean_plate_bug)*100))
    return(df)
  }
  
  else {
    plate_mean_bug <- dplyr::group_by(df, Time_points) %>%
      dplyr::filter(Concentration == 0.000) %>%
      dplyr::summarize(mean_plate_bug = mean(OD_no1st), med_plate_bug = median(OD_no1st), sd_plate_bug = sd(OD_no1st), se_bug = sd(OD_no1st)/sqrt(length(OD_no1st)),
                       min_plate_bug= min(OD_no1st), max_plate_bug = max (OD_no1st))
    
    df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(POC_bug = ((OD_no1st/mean_plate_bug)*100))
    return(df)
  }
}

#POC bugwells with scaled OD
POC_bugwells_scaled = function (df) {
  plate_mean_bug <- dplyr::group_by(df, Time_points) %>%
    dplyr::filter(Drug == "BUG") %>%
    dplyr::summarize(mean_plate_bug_scaled = mean(OD_scaled), med_plate_bug_scaled = median(OD_scaled), sd_plate_bug_scaled = sd(OD_scaled), se_bug_scaled = sd(OD_scaled)/sqrt(length(OD_scaled)),
                     min_plate_bug_scaled = min(OD_scaled), max_plate_bug_scaled = max (OD_scaled))
  
  df = inner_join(df, plate_mean_bug, by = "Time_points") %>% mutate(POC_bug_scaled = ((OD_scaled/mean_plate_bug_scaled)*100))
  return(df)
}

#basic function to calculate AUC for one well
AUC_calc_ODno1st <- function (x) {
  AUCs <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_well_ODno1st = DescTools::AUC(Time_h, OD_no1st))
  x= inner_join(x, AUCs, by = "wellID")
  return(x)
}

#list_all <- lapply (list_all, AUC_calc_ODno1st)


#proper ranking function (120 is 15*8 time points)
POC_top4 <- function (df) {
  ranked_AUC <- df %>% mutate(AUC_rank = BiocGenerics::rank(AUC_well_ODno1st, ties.method = "min")) %>%
    arrange(desc(AUC_rank))
  top15 <- ranked_AUC[1:120, ]
  
  plate_mean_top15 <- dplyr::group_by(top15, Time_points) %>%
    dplyr::summarize(mean_plate_top = mean(OD_no1st), med_plate_top = median(OD_no1st), sd_plate_top = sd(OD_no1st), se_top = sd(OD_no1st)/sqrt(length(OD_no1st)),
                     min_plate_top= min(OD_no1st), max_plate_top = max (OD_no1st))
  
  df = inner_join(df, plate_mean_top15, by = "Time_points") %>% mutate(POC_top = ((OD_no1st/mean_plate_top)*100))
  return(df)
  
}

#fix Nan
#make is.nan work for dfs (normally is doesnt)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

fix_nan <- function (x) {
  x[is.nan(x)] <- 0
  return(x)
}

#list_all <- lapply (list_all, fix_nan)

#recompute AUC using POCs instead of OD_no1st
AUC_raw <- function (x) {
  AUCs <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_well = DescTools::AUC(Time_h, OD))
  x= inner_join(x, AUCs, by = "wellID")
  return(x)
}


AUC_calc_poc_top <- function (x) {
  AUCs <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_well_top = DescTools::AUC(Time_h, POC_top))
  x= inner_join(x, AUCs, by = "wellID")
  return(x)
}

AUC_calc_poc_bug <- function (x) {
  AUCs <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_well_bug = DescTools::AUC(Time_h, POC_bug))
  x= inner_join(x, AUCs, by = "wellID")
  return(x)
}

#PROPER FUNCTION: NORMALISE AUC FOR THE MEDIAN AUC OF BUG WELLS

#always run POC_bugwells before this!!!
AUC_med_bug <- function (x) {
  test <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_bug = DescTools::AUC(Time_h, med_plate_bug))
  x= dplyr::inner_join(x, test, by = "wellID")%>%
    plyr::mutate(x, AUC_norm_bug = AUC_well_ODno1st/AUC_bug)
  return(x)
}

#screen_50_AUCnorm <- lapply(screen_50, AUC_med_bug)

#alternative: normalise for the AUC of the top4% wells
AUC_med_top <- function (x) {
  test <- dplyr::group_by(x,wellID) %>%
    dplyr::summarise(AUC_top = DescTools::AUC(Time_h, med_plate_top))
  x= dplyr::inner_join(x, test, by = "wellID")%>%
    plyr::mutate(x, AUC_norm_top = AUC_well_ODno1st/AUC_top)
  return(x)
}
#screen_50_AUCnorm <- lapply(screen_50, AUC_med_top)


#B score calculation (normalised value according to median polish)
Bscore <- function(df){
  tps =  unique (df$Time_points)
  B_list = list()
  for (t in tps) {
    plate <- filter(df, Time_points == t)
    #casting the signal to a matrix
    values <- acast(plate, rowID ~ colID, value.var = "OD_no1st")
    #Tukey's two way median polish
    med.pol <- medpolish(values, trace.iter=F, na.rm=T)
    
    #extract values
    estPlateAverage <- med.pol$overall
    rowEffects <- med.pol$row
    colEffects <- med.pol$col
    
    #computation function
    calcB <- function(x, mu_p, rE, cE){ return(x[["OD_no1st"]] - (mu_p + rE[x$rowID] + cE[x$colID]))}
    #calculate Bscore for each well
    result <- adply(plate,1, calcB, mu_p=estPlateAverage, rE=rowEffects, cE=colEffects)
    colnames(result)[ncol(result)] <- "Bscore"      
    result$Bscore <- result$Bscore / mad(result[["OD_no1st"]], na.rm=T)
    B_list[[i]] = result
  }
  B_unl <- plyr::ldply(B_list, 'data.frame')
  return(B_unl)
}

#list_all_B <- lapply(list_all, Bscore)


extract_bug <- function(df) {
  bug <- df[df$Drug == "BUG", ]
  bug <- dplyr::group_by (bug, Drug, Time_points, plate_name, run, donor_conc, plate_tp) %>%
    dplyr::summarise (med_plate = median(OD_no1st), sd=sd(OD_no1st), se = sd(OD_no1st)/sqrt(length(OD_no1st)),
                      mean_plate = mean(OD_no1st), min_plate = min(OD_no1st), max_plate = max(OD_no1st))
  return(bug)
}



#rescaling/shifting OD to get AUCs with values between 0 and 1

rescale_OD <- function(df) {
  OD_scale_shift <- dplyr::group_by(df,Well) %>%
    dplyr::mutate( OD_scaled = (OD595 - min(OD595)) / abs(1 - min(OD595)),
                   OD_shifted = OD595 - min(OD595) ) 
  df = merge(df,OD_scale_shift)
  return(df)
}

#recalculate all AUCs based on the scaled and shifted value
AUC_calc_ODscaled = function (x) {
  AUCs <- dplyr::group_by(x,Well) %>%
    dplyr::summarise(AUC_well_ODscaled = DescTools::AUC(Time_h, OD_scaled))
  x= inner_join(x, AUCs, by = "Well")
  return(x)
}


AUC_calc_ODshifted <- function (x) {
  AUCs <- dplyr::group_by(x,Well) %>%
    dplyr::summarise(AUC_well_ODshifted = DescTools::AUC(Time_h, OD_shifted))
  x= inner_join(x, AUCs, by = "Well")
  return(x)
}


#this normalises AUCs again against BUG wells but using scaled values
AUC_med_bug_scaled <- function (x) {
  test <- dplyr::group_by(x,Well) %>%
    dplyr::summarise(AUC_bug_scaled = DescTools::AUC(Time_h, med_plate_bug_scaled))
  x= dplyr::inner_join(x, test, by = "Well")%>%
    plyr::mutate(x, AUC_scaled_norm_bug = AUC_well_ODscaled / AUC_bug_scaled)
  return(x)
}


#basic AUC function in case DescTools doesnt work
calcAUC <- function(df) {
  N <- length(df$OD_no1st)
  if (N<2) return(NA)
  s1 <- 1:(N-1)
  s2 <- s1 + 1
  sum( (df$OD_no1st[s1] + df$OD_no1st[s2])/2 * (df$Time_h[s2]-df$Time_h[s1]) )
}

#calculate robust mean of BUG wells and normalise for that
robustMean <- function(x){
  if(length(x) == 1){return(x)}
  else if (unique(x) == 0) {return (mean(x))}
  else{
    return(smoothmest::smhuber(x)$mu)
  } 
}

#write pearson corr coeff in plots
corr_eqn = function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}


#convert factor in numeric (fast)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#cool function to get sd ,se and ci 95%
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
