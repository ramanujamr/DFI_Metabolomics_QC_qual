
library(shiny)
library(shinythemes)
library(shinyjs)
library(shinyalert)
library(shinyBS)

library(tidyverse)
library(reshape2)
library(DT)

library(broom)
library(ggsci)
library(gridExtra)

library(rhandsontable)
library(yingtools2)
library(ComplexHeatmap)
library(circlize)
library(janitor)



# 0. Initializations ###################################################################################################

zero_threshold=1000


# 1. set up directory and files ########################################################################################

#wddir <- "/Volumes/chaubard-lab/shiny_workspace/csvs/"
wddir <- "/Users/ramanujam/GitHub/dfimmfshiny_test/test_files"


# 2. FUNCTIONS #########################################################################################################

## 2.1 Read and clean data from input file =============================================================================

# For bile acids panel
Function_readin_csv_1 <- function(filename, zero_threshold, recursive=F)
{
  df_input <- read.csv(file = filename) %>%
    select(-CAS.ID) %>%
    reshape2::melt(id.vars=c("Compound.Name", "Formula", "Mass", "RT")) %>%
    replace_na(list(value = 0)) %>%
    mutate(itsd = str_extract(Compound.Name,pattern="ITSD"),
           com = gsub("\\_ITSD","",Compound.Name),
           Data.File=variable) %>%
    separate(com,into=c("num","compound_name","letter"),sep="\\_") %>%
    separate(variable,into=c("num2","date_run","batch","sampleid","conc"),sep="\\_\\_") %>%
    mutate(num2 = gsub("[Xx]", "", num2),
           sampleid = paste(num2, sampleid, sep = "_")) %>% 
    select(Data.File, sampleid, date_run, Compound.Name, compound_name,
           batch, letter, itsd, conc, peakarea=value) %>%
    mutate(filename = filename,
           compound_name = gsub("D[0-9]+\\-","",compound_name),
           compound_name = tolower(compound_name),
           conc = ifelse(grepl("dil",conc),"diluted","concentrated"),
           peakarea = as.numeric(peakarea),
           peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
           cc = str_extract(sampleid, pattern=".*CC.*"),
           cc = sub(".*_", "", cc))
  
  return(df_input)
}

# For non bile acids panels
Function_readin_csv_2 <- function(filename, zero_threshold, recursive=F){
  
  df_input_raw <- read.csv(file=filename, check.names=F)
  colnames(df_input_raw)[2] <- "garbage"
  colnames(df_input_raw)[3] <- "sampleid"
  colnames(df_input_raw)[4] <- "Data.File"
  colnames(df_input_raw)[5] <- "Type"
  colnames(df_input_raw)[6] <- "Level"
  colnames(df_input_raw)[7] <- "Acq.Date.Time"
  colnames(df_input_raw) <- gsub(" Results","",colnames(df_input_raw))
  df_input_raw = df_input_raw[-1,]
  
  df_input <- df_input_raw %>% 
    select(-Sample, -garbage,-Type,-Level,-Acq.Date.Time, -Data.File) %>% 
    slice(-1) %>% 
    reshape2::melt(id.vars="sampleid") %>% 
    mutate(value = as.numeric(value)) %>% 
    replace_na(list(value = 0)) %>% 
    dplyr::rename(compound_name = variable) %>% 
    mutate(itsd = str_extract(compound_name, pattern="ITSD"),
           compound_name = tolower(gsub("\\_ITSD","",compound_name))) %>% 
    separate(sampleid, into=c("inj_num","date_run","batch","sampleid","conc"),sep="\\_\\_") %>% 
    select(sampleid, date_run, batch, compound_name, itsd, conc, peakarea=value) %>% 
    mutate(conc = ifelse(grepl("dil",conc),"diluted","concentrated"),
           peakarea = ifelse(peakarea <= zero_threshold, 0, peakarea),
           cc = str_extract(sampleid, pattern=".*CC.*"),
           cc = sub(".*_", "", cc))
    
  return(df_input)
}






# 2. Execution #########################################################################################################

source('ui.R', local=TRUE)
source('server.R', local=TRUE)
shinyApp(ui=ui, server=server)
#runApp(list(ui=ui, server=server), host="0.0.0.0",port=1000)