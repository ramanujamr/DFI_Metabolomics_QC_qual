# Server - qual

server <- function(input, output, session) {
  
  rvalues <- reactiveValues(df_input=NULL, df_compounds=NULL, df_heatmap=NULL, mat_normalized=NULL,
                            plot_boxplots=NULL, compounds_list=NULL, plot_ht=NULL)
  
  shinyjs::hide("Button_generate_boxplots")
  shinyjs::hide("Button_generate_heatmap")
  shinyjs::hide("Button_new_file")
  shinyjs::hide("Button_itsd_stats")
  
  # Refresh csv files list
  observeEvent(input$Button_refresh_csv, ignoreInit = T, ignoreNULL = T, {
    updateSelectInput(session, 'filename', choices = rev(list.files(wddir, 
                                                                    pattern = ".*bile.*csv|.*PFBBr.*csv|.*Indole.*csv|.*Tryptophan.*csv|.*TMS.*csv", 
                                                                    ignore.case = T)))
  })
  
  
  ## 1.2 Upload new file ===============================================================================================
  observeEvent(input$Button_new_file, ignoreInit = T, ignoreNULL = T, {
    session$reload()
  })
  
  
  # 1. INPUT DATA ######################################################################################################
  
  observeEvent(input$Button_upload_csv, ignoreInit = T, ignoreNULL = T, {
    
    
    shinyjs::hide("filename")
    shinyjs::hide("Button_refresh_csv")
    shinyjs::hide("Button_upload_csv")
    
    shinyjs::show("Button_new_file")
    shinyjs::show("Button_itsd_stats")
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
    
    rvalues$zero_threshold <- ifelse(rvalues$panel=="Tryptophan", 100, 1000)

    # Read and clean input data
    if (rvalues$panel == "BileAcids") {
    rvalues$df_input <- Function_readin_csv_1(filename, rvalues$zero_threshold)
    } else {
      rvalues$df_input <- Function_readin_csv_2(filename, rvalues$zero_threshold)
    }
  
    
    ## 1.2 Average and median of ITSD (by sampleid for normalization)  =================================================
    
    if (rvalues$panel == "BileAcids") { # For bile acids panel
      rvalues$df_itsd_samples <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid, letter) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
    } else if(rvalues$panel == "TMS") {  # For TMS
      
      rvalues$df_itsd_samples <- rvalues$df_input %>%
        filter(itsd=="ITSD") %>%
        group_by(sampleid) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
    } else if(rvalues$panel == "PFBBr") {
      
      temp1 <- rvalues$df_input %>%
        filter(itsd == "ITSD") %>%
        filter(compound_name == "phenol" | compound_name == "proline_d7") %>% 
        filter(conc == "concentrated") %>% 
        group_by(sampleid, conc) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
      
      temp2 <- rvalues$df_input %>%
        filter(itsd == "ITSD") %>%
        filter(compound_name == "valerate" | compound_name == "valine_d8") %>% 
        filter(conc == "diluted") %>% 
        group_by(sampleid, conc) %>%
        summarize(avg = mean(peakarea),
                  med = median(peakarea))
      
      rvalues$df_itsd_samples = rbind(temp1, temp2)
    
    } else if(rvalues$panel == "Tryptophan") {
      
      rvalues$df_itsd_samples <- rvalues$df_input %>%
        filter(itsd == "ITSD") %>%
        filter(grepl("serotonin|melatonin", compound_name)) %>% 
        group_by(sampleid, conc) %>%
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
    
    if(rvalues$panel=="PFBBr") {
      rvalues$df_compounds <- rvalues$df_compounds %>% 
      mutate(conc = "concentrated")
    }
    
    if(rvalues$panel=="BileAcids") {
      rvalues$df_compounds <- rvalues$df_compounds %>% 
        select(-group_compounds) %>% 
        left_join(bile_acids_groups, by="compound_name") %>% 
        mutate(group_compounds = replace_na(group_compounds, "1")) %>% 
        mutate(group_compounds = factor(group_compounds, levels=compounds_categories, ordered = T))
        
    }
    
    
    
    ## 1.4 List of compounds ===========================================================================================
    rvalues$compounds_list <- rvalues$df_compounds %>% 
      distinct(compound_name) %>%
      `$`(compound_name)

    
    ## 1.2 ITSD Stat by compound (for QC)  =============================================================================
    
    rvalues$df_itsd_stats <- rvalues$df_input %>%
      filter(itsd=="ITSD") %>%
      filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= rvalues$zero_threshold, 0, peakarea)) %>% 
      group_by(batch, compound_name) %>%
      summarise(stdev = sd(peakarea),
                average = mean(peakarea),
                middle = median(peakarea),
                cv = stdev / average,
                cv_med = stdev / median(peakarea)) # Don't turn into % here since it will be applied in the y-axis scale
    
    rvalues$df_itsd <- rvalues$df_input %>% 
      filter(itsd == "ITSD") %>% 
      filter(!grepl("MB|Pooled|Plasma|Standard|Spiked",sampleid, ignore.case = T)) %>% 
      mutate(peakarea = ifelse(peakarea <= rvalues$zero_threshold, 0, peakarea),
             num = str_extract(sampleid, "[0-9][0-9][0-9]"),
             num = as.numeric(num),
             cc_shape = ifelse(grepl("CC[0-9]+", sampleid), "CC Sample", "ITSD")) %>% 
      left_join(rvalues$df_itsd_stats) %>% 
      mutate(num = as.numeric(num),
             flag = ifelse(peakarea > average + (1.5 * stdev) | peakarea < average - (1.5 * stdev), 
                           paste(num,batch,sampleid, sep = "_"), NA))

    
  })
  
  # ITSD stats (bsModal)
  output$Table_ITSD_stats <- DT::renderDataTable({
    rvalues$df_itsd_stats %>% 
      dplyr::rename(Batch = batch,
                    `Internal Standard` = compound_name,
                    StDev = stdev,
                    Mean = average,
                    Median = middle,
                    CV = cv,
                    `CV Median` = cv_med) %>%
      mutate(StDev = round(StDev, digits = 0),
             Mean = round(Mean, digits = 0),
             Median = round(Median, digits = 0),
             `CV (%)` = round(CV * 100, digits = 1),
             `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
      select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`) %>% 
      datatable(options = list(columnDefs = list(list(className='dt-center', targets="_all"))))
    
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

    # saveRDS(hot_to_r(input$Table_compounds_settings), "Table_compounds_settings.rds")
    # saveRDS(hot_to_r(input$Table_samples_settings), "Table_samples_settings.rds")
      
    # Table_compounds_settings <- readRDS("Table_compounds_settings.rds")
    # Table_samples_settings <- readRDS("Table_samples_settings.rds")
    
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
        inner_join(rvalues$df_itsd_samples, by=c("sampleid", "letter")) %>%
        mutate(norm_peak = peakarea / avg) %>% 
        mutate(norm_peak = ifelse(peakarea < rvalues$zero_threshold, NA, norm_peak)) %>% 
        mutate(norm_peak = ifelse(is.infinite(norm_peak), 0, norm_peak))
    } else {
      
      rvalues$df_normalized <- rvalues$df_input %>%
        filter(is.na(itsd)) %>%
        inner_join(rvalues$df_conc_type, by=c("compound_name", "conc")) %>%
        inner_join(rvalues$df_itsd_samples) %>%
        mutate(norm_peak = peakarea / avg) %>% 
        mutate(norm_peak = ifelse(peakarea < rvalues$zero_threshold, NA, norm_peak)) %>% 
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
    #   left_join(rvalues$df_MB_mean, by="compound_name") %>%
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
    
    
    

    
    heatmap_width <- nrow(rvalues$df_samples)*0.7
    

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
                                 row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                                 row_title_rot = 0,
                                 row_split = df_row_split$group_compounds,
                                 column_split = df_column_split$group_samples,
                                 cluster_rows = input$Checkbox_cluster_compounds,
                                 cluster_columns = input$Checkbox_cluster_samples,
                                 show_heatmap_legend = TRUE)
      


      rvalues$lgd1 = Legend(title = "Not Detected", labels = "", at = 1:1, legend_gp = gpar(fill = 8:9), title_position = c("topcenter"))
      draw(rvalues$lgd1, x = unit(0.7, "npc"), y = unit(0.99, "npc"), just = c("right", "top"))
      draw(rvalues$plot_ht, heatmap_legend_side="top")
      


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
        dplyr::select(sampleid, compound_name, norm_peak) %>% 
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
      paste0(rvalues$panel, "_Heatmap_",unique(rvalues$df_input$batch)[1],"_",Sys.Date(),".pdf")
    },
    
    content = function(file) {
      pdf(file, height = nrow(rvalues$mat_normalized)/2, width = max(ncol(rvalues$mat_normalized)/2, 30))
      draw(rvalues$plot_ht, heatmap_legend_side = "top")
      draw(rvalues$lgd1, x = unit(0.7, "npc"), y = unit(0.99, "npc"), just = c("right", "top"))
      plot.new()
      if(nrow(rvalues$nd_compounds) > 0)
      {
        print( gridExtra::grid.arrange(gridExtra::tableGrob(rvalues$nd_compounds, rows = NULL)) )
      }
      dev.off()
    })
  
  
  
  ## 4.6 QC report =====================================================================================================
  
  output$Button_download_qc_report <- downloadHandler(
    
    filename = function(){
      paste0(rvalues$panel, "_QC_Report_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },
    
    content = function(file) {
      
      plot_itsd <- rvalues$df_itsd %>%
        ggplot(., aes(x = num, y = peakarea, group = compound_name, label = flag)) +
        geom_point(aes(color = compound_name, shape = cc_shape), fill = "black", alpha = 0.6, size = 2) +
        geom_line(aes(y = average, color = compound_name)) +
        geom_line(aes(y = average + 1.5*stdev, color = compound_name), linetype = "dashed") +
        geom_line(aes(y = average - 1.5*stdev, color = compound_name), linetype = "dashed") +
        geom_ribbon(aes(ymin = average - stdev, ymax = average + stdev,
                        fill = compound_name), alpha=0.2) +
        ggrepel::geom_label_repel(size = 1.2, max.overlaps = Inf,
                                  min.segment.length = 0.1, label.padding = 0.1, na.rm=T) +
        theme_bw()+
        theme(panel.grid.minor= element_blank(),
              panel.grid.major.x = element_blank(),
              legend.text = element_text(color = "black", size = 8),
              legend.title = element_text(color = "black", size = 10),
              legend.position = "top",
              strip.text=element_text(color = "black", size=5),
              axis.text =element_text(color = "black", size=8),
              axis.title = element_text(color = "black", size = 10),
              plot.margin = margin(1,0.5,0,0.6, unit = 'cm')) +
        scale_fill_manual(values = c(paletteer::paletteer_d("ggsci::default_igv", length(rvalues$df_itsd_stats$compound_name)))) +
        guides(color = guide_legend(title = "Internal Standard Compound",
                                    override.aes = list(size = 2.5), nrow = 2,
                                    title.position="top", title.hjust = 0.5,
                                    label.position = "right"), fill = "none",
               shape = guide_legend(title = "",
                                    override.aes = list(size = 2.5), nrow = 2,
                                    title.position="top", title.hjust = 0.5,
                                    label.position = "right")) +
        scale_shape_manual(values = c(24,16))+
        scale_y_continuous(label = scales::scientific) +
        scale_x_continuous(breaks = seq(0,150,25)) +
        ylab("Raw Peak Area\n") +
        xlab("\nInjection Number")+
        facet_wrap_paginate(~ compound_name, ncol = 2, nrow = 3) +
        ggtitle(paste0(rvalues$panel, " Normalized QC Report - ", unique(rvalues$df_input$batch)[1]))
      
      
      
      
      
      pdf(file, onefile=T, height = 11, width = 12)

      for(i in 1:n_pages(plot_itsd))
      {
        print(plot_itsd + facet_wrap_paginate(~ compound_name, ncol = 2, nrow = 3, page = i))
      }
      
      temp <- rvalues$df_itsd_stats %>% 
        dplyr::rename(Batch = batch,
                      `Internal Standard` = compound_name,
                      StDev = stdev,
                      Mean = average,
                      Median = middle,
                      CV = cv,
                      `CV Median` = cv_med) %>%
        mutate(StDev = round(StDev, digits = 0),
               Mean = round(Mean, digits = 0),
               Median = round(Median, digits = 0),
               `CV (%)` = round(CV * 100, digits = 1),
               `CV Median (%)` = round(`CV Median` * 100, digits = 1)) %>%
        select(Batch,`Internal Standard`,StDev,Mean,Median,`CV (%)`,`CV Median (%)`)
      
      
      print( gridExtra::grid.arrange(gridExtra::tableGrob(temp, rows = NULL)) )
      
      invisible(dev.off())
    })
  
  
  
}