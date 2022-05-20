
server <- function(input, output, session) {
  
  rvalues <- reactiveValues(df_input=NULL, df_compounds=NULL, df_heatmap=NULL, mat_normalized=NULL,
                            plot_boxplots=NULL, compounds_list=NULL, plot_ht=NULL)
  
  shinyjs::hide("Button_generate_boxplots")
  shinyjs::hide("Button_generate_heatmap")
  
  # Refresh csv files list
  observeEvent(input$Button_refresh_csv, ignoreInit = T, ignoreNULL = T, {
    updateSelectInput(session, 'filename', choices = rev(list.files(wddir, 
                                                                    pattern = ".*bile.*csv|.*PFBBr.*csv|.*Indole.*csv|.*Tryptophan.*csv", 
                                                                    ignore.case = T)))
  })
  
  
  # 1. INPUT DATA ######################################################################################################
  
  observeEvent(input$Button_upload_csv, ignoreInit = T, ignoreNULL = T, {
    
    shinyjs::show("Button_generate_boxplots")
    shinyjs::show("Button_generate_heatmap")
    
    ## 1.1 Get input data  =============================================================================================
    
    filename <- file.path(wddir,input$filename)
    
    # Outputs
    output$Textout_panel <- renderText({paste("Panel: ", rvalues$panel)}) # Display identified panel
    output$Textout_filename <- renderText({paste("Uploaded file: ", input$filename)}) # Display filename
    
    # Identify panel and read equivalent quant compounds csv
    rvalues$panel <- case_when(
      grepl("Bile", input$filename, ignore.case = T) ~ "BileAcids",
      grepl("PFBBr", input$filename, ignore.case = T) ~ "PFBBr",
      grepl("Tryptophan|Indole", input$filename, ignore.case = T) ~ "Tryptophan",
      grepl("TMS", input$filename, ignore.case = T) ~ "TMS",
      TRUE ~ "NOT IDENTIFIED... CHECK FILENAME")
    
    # Read and clean input data
    if (rvalues$panel == "BileAcids") {
    rvalues$df_input <- Function_readin_csv_1(filename, zero_threshold) 
    } else {
      rvalues$df_input <- Function_readin_csv_2(filename, zero_threshold) 
    }
  
    
    ## 1.2 Average and median of ITSD (PANEL SPECIFIC)  ================================================================
    
    if (rvalues$panel == "BileAcids") { # For bile acids panel
      rvalues$df_itsd_stats <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid, letter) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
    } else {  # For non bile acids panel
      
      rvalues$df_itsd_stats <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
    }

    
    ## 1.3 Input for rhandsontable (concentration selection & grouping) ================================================

    compounds_categories <- c("1","2","3","4","5", bile_categories) 
    
    rvalues$df_compounds <- rvalues$df_input %>% 
      filter(is.na(itsd)) %>%
      select(compound_name, conc) %>% 
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(compound_name) %>%
      group_by(compound_name) %>% 
      mutate(found_dil = ifelse(any(grepl("diluted", conc)), "Yes", "No"),
             found_conc = ifelse(any(grepl("concentrated", conc)), "Yes", "No")) %>% 
      mutate(conc = factor(conc, levels = sort(unique(conc)))) %>% 
      distinct(compound_name, .keep_all = T) %>% 
      ungroup() %>% 
      mutate(group_compounds = factor(1, levels = compounds_categories, ordered = TRUE))
    
    
    ## 1.4 List of compounds ===========================================================================================
    compounds_list <- rvalues$df_compounds %>% 
      distinct(compound_name) %>%
      `$`(compound_name)
    
    
  })
    
  
  
  # 2. RAW DATA ########################################################################################################
  
  
  observeEvent(input$Button_generate_boxplots, ignoreInit = T, ignoreNULL = T, {

    output$Table_compounds_settings <- renderRHandsontable({
      if (!is.null(rvalues$df_compounds))
        rhandsontable(rvalues$df_compounds, rowHeaders = NULL,overflow = 'visible') %>%
        hot_col("compound_name", readOnly = T) %>%
        hot_col("conc", allowInvalid = F) %>% 
        hot_col("found_conc", readOnly = T, halign = "htCenter") %>% 
        hot_col("found_dil", readOnly = T, halign = "htCenter") %>% 
        hot_col("group_compounds", allowInvalid = F) 
    })
    
    # Table for grouping samples (rHandsontable)
    
    rvalues$df_samples <- rvalues$df_input %>%
      filter(!grepl("MB",sampleid, ignore.case = T),
             !grepl("Pooled",sampleid, ignore.case = T),
             !grepl("Plasma",sampleid, ignore.case = T),
             !grepl("Standard",sampleid, ignore.case = T),
             !grepl("CC[0-9]+", sampleid)) %>%
      distinct(sampleid) %>%
      select(sampleid) %>%
      arrange(sampleid) %>%
      mutate(group_samples = factor(1, levels = c(1,2,3,4,5), ordered = TRUE))

    output$Table_samples_settings <- renderRHandsontable({
      if (!is.null(rvalues$df_samples))
        rhandsontable(rvalues$df_samples, rowHeaders = NULL, overflow = 'visible') %>%
        hot_col("sampleid", readOnly = T) %>%
        hot_col("group_samples", allowInvalid = F)
    })


    # Boxplots
    
    output$Plot_boxplots <- renderPlot({
      
      rvalues$plot_boxplots <- rvalues$df_input %>%
        filter(is.na(itsd)) %>%
        ggplot(aes(y=factor(compound_name, levels = rev(levels(factor(compound_name)))), x=peakarea)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width=0.1,height=0.1,alpha=0.3) +
        theme_bw() +
        theme(strip.text.x=element_text(size=15),
              axis.text.y=element_text(size=13)) +
        ylab("") +
        xlab("peak area") +
        facet_grid(itsd ~ conc,scales="free_y",space="free") +
        scale_x_continuous(trans=log_epsilon_trans(epsilon=1000))
      
      rvalues$plot_boxplots
    }, height=nrow(rvalues$df_compounds)*24.5)
    
    
  })
  
 
  # 3. HEATMAP #########################################################################################################
 
  observeEvent(input$Button_generate_heatmap, ignoreInit = T, ignoreNULL = T, {

     #saveRDS(hot_to_r(input$Table_compounds_settings), "Table_compounds_settings.rds")
     #saveRDS(hot_to_r(input$Table_samples_settings), "Table_samples_settings.rds")
     
    #Table_compounds_settings <- readRDS("Table_compounds_settings.rds")
    #Table_samples_settings <- readRDS("Table_samples_settings.rds")
    
    
    ## 3.1 Get concentration selection information from Table_compounds_settings =======================================
    compounds_dil <- hot_to_r(input$Table_compounds_settings) %>%
      filter(conc=="diluted") %>%
      distinct(compound_name) %>%
      pull(compound_name)


    compounds_conc <- hot_to_r(input$Table_compounds_settings) %>%
      filter(conc=="concentrated") %>%
      distinct(compound_name) %>%
      pull(compound_name)

    if( length(compounds_conc) >0){
      df_conc <- data.frame(compounds_conc, conc="concentrated") %>% dplyr::rename(compound_name=compounds_conc)
    }else{
      df_conc <- data.frame(compound_name= character(0), conc= character(0))
    }

    if( length(compounds_dil) >0){
      df_dil <- data.frame(compounds_dil, conc="diluted") %>% dplyr::rename(compound_name=compounds_dil)
    }else{
      df_dil <- data.frame(compound_name= character(0), conc= character(0))
    }

    rvalues$df_conc_type <- rbind(df_conc, df_dil)


    ## 3.2 Normalized dataframe ========================================================================================

    if (rvalues$panel == "BileAcids") {
      
      rvalues$df_normalized <- rvalues$df_input %>%
        filter(is.na(itsd)) %>%
        inner_join(rvalues$df_conc_type, by=c("compound_name", "conc")) %>%
        inner_join(rvalues$df_itsd_stats, by=c("sampleid", "letter")) %>%
        mutate(norm_peak = peakarea / avg) %>% 
        mutate(norm_peak = ifelse(peakarea < zero_threshold, NA, norm_peak)) %>% 
        mutate(norm_peak = ifelse(is.infinite(norm_peak), 0, norm_peak))
    } else {
      
      rvalues$df_normalized <- rvalues$df_input %>%
        filter(is.na(itsd)) %>%
        inner_join(rvalues$df_conc_type, by=c("compound_name", "conc")) %>%
        inner_join(rvalues$df_itsd_stats, by="sampleid") %>%
        mutate(norm_peak = peakarea / avg) %>% 
        mutate(norm_peak = ifelse(peakarea < zero_threshold, NA, norm_peak)) %>% 
        mutate(norm_peak = ifelse(is.infinite(norm_peak), 0, norm_peak))
    }
    
    # Method blanks mean dataframe
    rvalues$df_MB_mean <- rvalues$df_normalized %>%
      filter(grepl("MB",sampleid)) %>% 
      group_by(compound_name) %>% 
      summarize(mean_mb = mean(norm_peak)) %>% 
      mutate(mean_mb = ifelse(is.na(mean_mb), 0, mean_mb))
 
    # Subtract mean MB if checkbox is selected
    rvalues$df_normalized <- rvalues$df_normalized %>%
      left_join(rvalues$df_MB_mean, by="compound_name") %>%
      rowwise() %>%
      mutate(norm_peak = ifelse(input$Checkbox_subtract_MB==T, max(norm_peak-mean_mb,0), norm_peak)) %>% 
      mutate(norm_peak = ifelse(norm_peak==0, NA, norm_peak))
    
    
    # df_normalized <- df_normalized %>%
    #   left_join(df_MB_mean, by="compound_name") %>%
    #   rowwise() %>%
    #   mutate(norm_peak = max(norm_peak-mean_mb,0)) %>% 
    #   mutate(norm_peak = ifelse(norm_peak==0, NA, norm_peak))
    
    
    ## 3.3 Heatmap dataframe ===========================================================================================

    rvalues$df_heatmap <- rvalues$df_normalized %>%
      filter(!grepl("MB|Pooled|Plasma|Standard|CC|Spiked",sampleid, ignore.case = T)) %>% 
      group_by(compound_name) %>%
      mutate(compound_med = median(norm_peak, na.rm=T)) %>%
      ungroup() %>%
      mutate(heat_val = log((norm_peak / compound_med), base = 2)) %>% 
      left_join(hot_to_r(input$Table_compounds_settings), by=c("compound_name","conc")) %>%
      left_join(hot_to_r(input$Table_samples_settings), by="sampleid") 
    

    ## 3.4 Plot heatmap ================================================================================================

    rvalues$mat_normalized <- rvalues$df_heatmap %>%
      select(sampleid, compound_name, heat_val) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
      column_to_rownames(var = "sampleid") %>%
      t()

    # Get compounds with all NAs
    rvalues$nd_compounds <- rvalues$mat_normalized[rowSums(is.na(rvalues$mat_normalized)) == ncol(rvalues$mat_normalized), ]
    rvalues$nd_compounds <- rownames(rvalues$nd_compounds) %>% as.data.frame()
    if(nrow(rvalues$nd_compounds)!=0) {colnames(rvalues$nd_compounds) <- c("ND_compounds")}
    
    # Filter ND compounds (across all samples)
    rvalues$mat_normalized <- rvalues$mat_normalized[rowSums(is.na(rvalues$mat_normalized)) != ncol(rvalues$mat_normalized), ]

    # Replace NA with a very high number
    mat_normalized_na_removed <- rvalues$mat_normalized 
    mat_normalized_na_removed[is.na(mat_normalized_na_removed)] <- 99999

    df_column_split <- rvalues$df_heatmap %>% 
      filter(compound_name %in% rownames(mat_normalized_na_removed)) %>% 
      select(sampleid, group_samples) %>% 
      distinct()
    
    df_row_split <- rvalues$df_heatmap %>% 
      filter(compound_name %in% rownames(mat_normalized_na_removed)) %>%
      select(compound_name, group_compounds) %>% 
      distinct()
    
    hm_min <- min(rvalues$df_heatmap$heat_val, na.rm=T)
    hm_max <- max(rvalues$df_heatmap$heat_val, na.rm=T)
    hm_lim <- max(c((-1*hm_min),hm_max))

    color_min <- ceiling(hm_min)
    color_max <- ceiling(hm_max)
    color_lim <- ceiling(hm_lim)

    col_fun1 <- colorRamp2(c(-hm_lim, 0, hm_lim, 99999), c("blue", "white", "red", "grey"))
    
    
    

    
    heatmap_width <- nrow(rvalues$df_samples)*0.6
    

    output$Plot_heatmap <- renderPlot({

      rvalues$plot_ht <- Heatmap(as.matrix(mat_normalized_na_removed), col = col_fun1, width = unit(heatmap_width, "cm"), 
                                 row_names_max_width = unit(8, "cm"),
                                 heatmap_legend_param = list(title = 'Log2FoldChange', title_position = 'topcenter', direction="horizontal",
                                                             at = c(-color_lim, 0, color_lim)),
                                 row_names_side = c("left"),
                                 row_gap = unit(5, 'mm'),
                                 border_gp = gpar(col = "black", lty = 1),
                                 rect_gp = gpar(col = "black", lwd = 0.2),
                                 row_dend_side = "right",
                                 row_dend_width = unit(10, "cm"),
                                 column_dend_height = unit(4, "cm"),
                                 row_split = df_row_split$group_compounds,
                                 column_split = df_column_split$group_samples,
                                 cluster_rows = input$Checkbox_cluster_compounds,
                                 cluster_columns = input$Checkbox_cluster_samples,
                                 show_heatmap_legend = TRUE)
      

      draw(rvalues$plot_ht, heatmap_legend_side="top")
      
      
      # lgd1 = Legend(col_fun = col_fun1, title = 'Log2FoldChange',
      #               title_position = 'topcenter',
      #               legend_height = unit(6, "cm"),
      #               at = c(-color_lim, 0, color_lim))
      # 
      # lgd2 = Legend(title = "Not Detected", labels = "", at = 1:1, legend_gp = gpar(fill = 8:9), title_position = c("topcenter"))
      # 
      # pd = packLegend(lgd1, lgd2, direction = "horizontal", column_gap = unit(20, "mm"))
      # pushViewport(viewport(width = 1, height = 1))
      # grid.rect(gp = gpar(col = "white"))
      # 
      # draw(pd, x = unit(0.9, "npc"), y = unit(0.9, "npc"), just = c("right", "top"))
      # 
      # popViewport()
      
      
      

    }, height=max(nrow(rvalues$mat_normalized)*14.5, 500)
    )

 
  })
  
  
  # 4. SAVE RESULTS ####################################################################################################
  

  # Save normalized csv
  output$Button_download_normalized_csv <- downloadHandler(
    
    filename = function(){
      paste0("normalized_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      
      rvalues$df_normalized %>% 
        select(sampleid, compound_name, norm_peak) %>% 
        mutate(sampleid = ifelse(grepl("MB|Pooled|Plasma|CC|Standard|Spiked", sampleid, ignore.case = T),
                                 sampleid,
                                 gsub("^[0-9]{3}_", "", sampleid))) %>% 
        pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
        write.csv(file=file, row.names=F,quote=F)
    })
  
  
  
  # Save normalized csv (no QCs)
  output$Button_download_normalized_csv_no_qc <- downloadHandler(
    
    filename = function(){
      paste0("removed_qcs_normalized_results_",
             gsub("\\.csv","",input$filename),"_",
             gsub("\\-","",Sys.Date()),".csv")
    },
    
    content = function(file) {
      
      rvalues$df_normalized %>% 
        select(sampleid, compound_name, norm_peak) %>% 
        filter(!grepl("MB|Pooled|Plasma|CC|Standard|Spiked",sampleid, ignore.case = T)) %>% 
        mutate(sampleid = gsub("^[0-9]{3}_", "", sampleid)) %>% 
        pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
        write.csv(file=file, row.names = F, quote=F)
    })
  
  
  
  # Download heatmap pdf
  output$Button_download_heatmap <- downloadHandler(
    
    filename = function(){
      paste0(rvalues$panel, "_Heatmap_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },
    
    content = function(file) {
      pdf(file, height = nrow(rvalues$mat_normalized)/2, width = max(ncol(rvalues$mat_normalized)/2, 20))
      draw(rvalues$plot_ht, heatmap_legend_side = "top")
      plot.new()
      print( gridExtra::grid.arrange(gridExtra::tableGrob(rvalues$nd_compounds, rows = NULL)) )
      dev.off()
    })
  
  
}