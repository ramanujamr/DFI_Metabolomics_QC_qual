
server <- function(input, output, session) {
  
  rvalues <- reactiveValues(df_input=NULL, df_compounds=NULL, plot_boxplots=NULL, compounds_list=NULL, plot_ht=NULL)
  
  shinyjs::hide("Button_generate_boxplots")
  shinyjs::hide("Button_generate_heatmap")
  
  # Refresh csv files list
  observeEvent(input$Button_refresh_csv, ignoreInit = T, ignoreNULL = T, {
    updateSelectInput(session, 'filename', choices = list.files(wddir, pattern = "csv$"))
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
    
    rvalues$df_compounds <- rvalues$df_input %>% 
      filter(is.na(itsd)) %>%
      select(compound_name, conc) %>% 
      mutate(compound_name=as.character(compound_name)) %>%
      arrange(compound_name) %>%
      group_by(compound_name) %>% 
      mutate(found_dil = ifelse(any(grepl("diluted", conc)), "Yes", "No"),
             found_conc = ifelse(any(grepl("concentrated", conc)), "Yes", "No")) %>% 
      distinct() %>% 
      mutate(conc = factor(conc, levels = sort(unique(conc)))) %>% 
      ungroup() %>% 
      mutate(group_compounds = factor(1, levels = c("1","2","3","4","5"), ordered = TRUE))
    
    
    ## 1.4 List of compounds ===========================================================================================
    rvalues$compounds_list <- rvalues$df_compounds %>% 
      distinct(compound_name) %>%
      `$`(compound_name)
    
    
  })
  
  
  
  # 2. RAW DATA ########################################################################################################
  
  
  observeEvent(input$Button_generate_boxplots, ignoreInit = T, ignoreNULL = T, {

    # Table for conc/dil selection (rHandsontable)
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
    
    df_samples <- rvalues$df_input %>%
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
      if (!is.null(df_samples))
        rhandsontable(df_samples, rowHeaders = NULL, overflow = 'visible') %>%
        hot_col("sampleid", readOnly = T) %>%
        hot_col("group_samples", allowInvalid = F)
    })


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
        left_join(rvalues$df_itsd_stats, by=c("sampleid", "letter")) %>%
        mutate(norm_peak = peakarea / avg)
      
    } else {
      
      rvalues$df_normalized <- rvalues$df_input %>%
        filter(is.na(itsd)) %>%
        inner_join(rvalues$df_conc_type, by=c("compound_name", "conc")) %>%
        left_join(rvalues$df_itsd_stats, by="sampleid") %>%
        mutate(norm_peak = peakarea / avg)
    }
 


    ## 3.3 Heatmap dataframe ===========================================================================================

    df_heatmap <- rvalues$df_normalized %>%
      group_by(compound_name) %>%
      mutate(compound_med = median(norm_peak)) %>%
      ungroup() %>%
      mutate(heat_val = log((norm_peak / (compound_med +1))+1, base = 2)) %>%
      filter(!grepl("MB",sampleid, ignore.case = T),
             !grepl("Pooled",sampleid, ignore.case = T),
             !grepl("Plasma",sampleid, ignore.case = T),
             !grepl("Standard",sampleid, ignore.case = T),
             !grepl("CC[0-9]+", sampleid)) %>%
      mutate(heat_val = ifelse(is.infinite(heat_val), 0, heat_val)) %>%
      left_join(hot_to_r(input$Table_compounds_settings), by=c("compound_name","conc")) %>%
      left_join(hot_to_r(input$Table_samples_settings), by="sampleid")
    
    df_column_split <- df_heatmap %>% 
      select(sampleid, group_samples) %>% 
      distinct()
    
    df_row_split <- df_heatmap %>% 
      select(compound_name, group_compounds) %>% 
      distinct()
    

    ## 3.4 Plot heatmap ================================================================================================

    rvalues$mat_normalized <- df_heatmap %>%
      mutate(heat_val = ifelse(peakarea<zero_threshold, 0, heat_val)) %>%
      select(sampleid, compound_name, heat_val) %>%
      pivot_wider(names_from = compound_name, values_from = heat_val, values_fill = NA) %>%
      column_to_rownames(var = "sampleid") %>%
      janitor::remove_empty(which = "cols") %>%
      t()

    hm_min <- min(df_heatmap$heat_val)
    hm_max <- max(df_heatmap$heat_val)
    hm_lim <- max(c((-1*hm_min),hm_max))

    color_min <- ceiling((hm_min - 1))
    color_max <- ceiling(hm_max)
    color_lim <- ceiling(hm_lim)

    col_fun1 <- colorRamp2(c(-hm_lim, 0, hm_lim), c("blue", "white", "red"))

    output$Plot_heatmap <- renderPlot({

      rvalues$plot_ht <- Heatmap(as.matrix(rvalues$mat_normalized), col = col_fun1, width = unit(30, "cm"),
                    row_names_max_width = unit(8, "cm"),
                    heatmap_legend_param = list(title = 'Log2FoldChange', title_position = 'topcenter', direction="horizontal",
                                                legend_height = unit(6, "cm"),
                                                at = c(-color_lim,(ceiling(color_lim - color_lim/2) * -1),
                                                       0,
                                                       (ceiling(color_lim - color_lim /2)), color_lim)),

                    row_split = df_row_split$group_compounds,
                    column_split = df_column_split$group_samples,
                    row_names_side = c("left"),
                    row_gap = unit(5, 'mm'),
                    border_gp = gpar(col = "black", lty = 1),
                    rect_gp = gpar(col = "black", lwd = 0.2),
                    row_dend_width = unit(10, "cm"),
                    column_dend_height = unit(4, "cm"),
                    cluster_rows = input$Checkbox_cluster_samples,
                    cluster_columns = input$Checkbox_cluster_compounds,
                    show_heatmap_legend = TRUE)

      draw(rvalues$plot_ht, heatmap_legend_side = "top")

    })



    # Get compounds that are not detected in any of the samples
    nd_compounds <- rvalues$compounds_list[ !(rvalues$compounds_list %in% colnames(rvalues$mat_normalized)) ]


  })
  
  
  # 4. SAVE RESULTS ####################################################################################################
  
  # Save normalized csv
  observeEvent(input$Button_download_normalized_csv, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_normalized %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      mutate(sampleid = ifelse(grepl("MB|Pooled|Plasma|CC|Standard", sampleid, ignore.case = T),
                               sampleid,
                               gsub("^[0-9]{3}_", "", sampleid))) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/normalized_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
  })
  

  # Save normalized csv (no QCs)
  observeEvent(input$Button_download_normalized_csv_no_qc, ignoreInit = T, ignoreNULL = T, {
    
    rvalues$df_normalized %>% 
      select(sampleid, compound_name, norm_peak) %>% 
      filter(!grepl("MB|Pooled|Plasma|CC|Standard",sampleid, ignore.case = T)) %>% 
      mutate(sampleid = gsub("^[0-9]{3}_", "", sampleid)) %>% 
      pivot_wider(names_from = compound_name, values_from = norm_peak, values_fill = NA) %>% 
      write.csv(paste0("~/Downloads/removed_qcs_normalized_results_",
                       gsub("\\.csv","",input$filename),"_",
                       gsub("\\-","",Sys.Date()),".csv"), row.names=F,quote=F)
    
    shinyalert(title = "File saved to ~/Downloads", type = "success")
    
  })
  
  
  # Download heatmap pdf
  output$Button_download_heatmap <- downloadHandler(
    
    filename = function(){
      paste0(rvalues$panel, "_Heatmap_",unique(rvalues$df_input$batch),"_",Sys.Date(),".pdf")
    },
    
    content = function(file) {
      pdf(file, height = nrow(rvalues$mat_normalized)/2, width = ncol(rvalues$mat_normalized)/2)
      draw(rvalues$plot_ht, heatmap_legend_side = "top")
      dev.off()
    })
  
  
  
}