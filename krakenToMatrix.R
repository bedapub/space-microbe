#' convert Kraken2 output to sparse matrix (input for Seurat)
#' 
#' @param filePath A string. Path to input file.
#' @param counts A string. Either "umi_counts" or "read_counts". Default "umi_counts".
#' @param tax_level A string. Either "species", "genus" or "family". Default "genus".
#' @param remove_singletons A logical. Only applicable if counts=="umi_counts". Default TRUE.
#' @returns A list of $taxid_list (list with taxids ordered by abundance) and $matrix (matrix with rows=taxids, cols=barcodes).
#' @examples


krakenToMatrix <- function(filePath, counts="umi_counts", tax_level="genus", remove_singletons=TRUE){
  
  ##check parameters
  if (is.null(counts)){
    counts <- "umi_counts"
    warning("counts set to default: umi_counts")
  } else if (!(counts %in% c("read_counts", "umi_counts"))){
    stop("counts has to be either read_counts or umi_counts")
  }
  
  if (is.null(tax_level)){
    tax_level <- "genus"
    warning("tax_level set to default: genus")
  } else if (!(tax_level %in% c("species", "genus", "family"))){
    stop("tax_level has to be either species, genus or family")
  }
  
  if (is.null(remove_singletons)){
    remove_singletons <- TRUE
    warning("remove_singletons set to default: TRUE")
  } else if (!is.logical(remove_singletons)){
    stop("remove_singletons has to be a logical (TRUE or FALSE)")
  }
  
  if(!file.exists(filePath)){
    stop("file in filePath doesn't exist")
  }
  
  ##load libraries
  library(tidyverse)
  library(plyr)
  library(Matrix)
  library(taxonomizr, lib.loc="/home/anzboecs/R/x86_64-pc-linux-gnu-library/4.0.5-foss")
  
  ## define variables
  DIR<-"/projects/site/pred/microbiome/database/taxonomizr_DB/"
  FILE<- paste0(DIR,"nameNode.sqlite")
  whitelistPath <- "/projects/site/pred/microbiome/projects/ST_microbiome/umi_pipeline/05_kraken2/whitelist.tsv"
  
  ## import modified kraken2 output file
  data <- read.table(file=filePath, sep="\t", header=FALSE)
  colnames(data) <- c("barcode", "umi", "taxid")
  #trim whitespaces
  data$barcode <- trimws(data$barcode)
  data$umi <- trimws(data$umi)
  data$taxid <- trimws(data$taxid)
  
  ## extract valid barcodes from whitelist
  whitelist <- read.table(file=whitelistPath, header=FALSE)
  colnames(whitelist) <- c("barcode")
  data <- subset(data, data$barcode %in% whitelist$barcode)
  data <- tidyr::as_tibble(data)
  
  ##remove taxid 131567: cellular organisms and 1:root and 2:Bacteria and 2759:Eukaryota
  data <- data %>% filter(!taxid %in% c(131567, 1, 2, 2759))
  
  ########################## counts==read_counts #####################################################
  
  if (counts=="read_counts"){
    
    data <- data %>% select(barcode, taxid)
    
    if (tax_level=="genus"){
      ## remove all taxids > genus-level
      genus_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("genus"))
      data <- data[!is.na(genus_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("genus"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    }
    else if (tax_level=="species"){
      ## remove all taxids > species-level
      species_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("species"))
      data <- data[!is.na(species_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("species"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    } 
    else if (tax_level=="family"){
      ## remove all taxids > species-level
      family_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("family"))
      data <- data[!is.na(family_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("family"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    }
    
    ## collapse all reads with same barcode & taxid
    data <- plyr::ddply(data,.(barcode, taxid),nrow)
    colnames(data) <- c("barcode", "taxid", "counts")
    
    ## report most abundant taxa
    taxa <- data %>% select(taxid, counts)
    taxa_abundance <- aggregate(counts ~ taxid, data=taxa, FUN=sum) %>% arrange(desc(counts))
    print(paste0(sum(taxa_abundance$counts), " total read counts at ", tax_level, " level"))
    taxa_list <- taxa_abundance$taxid
    
    ## convert dataframe to sparse matrix 
    data$counts <- as.numeric(data$counts)
    data$taxid <- as.character(data$taxid)
    matrix <- tidyr::pivot_wider(data, names_from=taxid, values_from = counts, values_fill=0)
    #complete barcodes from whitelist with 0s
    matrix <- left_join(whitelist, matrix)
    matrix[is.na(matrix)] <- 0
    #add -1 to barcode column
    matrix$barcode <- paste0(matrix$barcode, "-1")
    matrix <- as.matrix(matrix)
    rownames(matrix) <- NULL
    rownames(matrix) <- matrix[,1]
    matrix <- matrix[,2:ncol(matrix)]
    #transpose matrix and convert matrix to dgCMatrix class 
    matrix <- t(matrix)
    colnames <- colnames(matrix)
    rownames <- rownames(matrix)
    matrix <- Matrix::Matrix(as.numeric(matrix), nrow=length(rownames), ncol=length(colnames), sparse=TRUE, dimnames=list(rownames, colnames))
    
    ## return
    return(list("taxid_list"=taxa_list, "matrix"=matrix))
  }
  
  
  ########################## counts==umi_counts #####################################################
  
  if (counts=="umi_counts"){
    ## get umi counts: collapse all reads with same taxid
    data <- plyr::ddply(data,.(barcode, umi, taxid),nrow)
    data <- data %>% select(barcode, umi, taxid)
    
    if (tax_level=="genus"){
      ## remove all taxids > genus-level
      genus_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("genus"))
      data <- data[!is.na(genus_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("genus"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    }
    else if (tax_level=="species"){
      ## remove all taxids > species-level
      species_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("species"))
      data <- data[!is.na(species_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("species"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    } 
    else if (tax_level=="family"){
      ## remove all taxids > species-level
      family_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("family"))
      data <- data[!is.na(family_level),]
      ## assign genus_taxid to all reads
      data$taxid <- taxonomizr::getId(taxa=getTaxonomy(ids=data$taxid,sqlFile=FILE, desiredTaxa=c("family"))[,1],sqlFile=FILE)
      #select first taxid in case of multiple taxids (usually bacterial)
      data$taxid <- gsub("\\,.*","",data$taxid)
    }

    ## transform tibble into a matrix with feature=taxid as row and barcode as column
    #drop umi and count column
    data <- plyr::ddply(data,.(barcode, taxid),nrow)
    colnames(data) <- c("barcode", "taxid", "counts")
    data$counts <- as.numeric(data$counts)
    data$taxid <- as.character(data$taxid)
    
    if(remove_singletons==TRUE){
      data <- data %>% filter(counts!=1)
    }
    
    ## report most abundant taxa
    taxa <- data %>% select(taxid, counts)
    taxa_abundance <- aggregate(counts ~ taxid, data=taxa, FUN=sum) %>% arrange(desc(counts))
    print(paste0(sum(taxa_abundance$counts), " total UMI counts at ", tax_level, " level"))
    taxa_list <- taxa_abundance$taxid
    
    matrix <- tidyr::pivot_wider(data, names_from=taxid, values_from = counts, values_fill=0)
    #complete barcodes from whitelist with 0s
    matrix <- left_join(whitelist, matrix)
    matrix[is.na(matrix)] <- 0
    #add -1 to barcode column
    matrix$barcode <- paste0(matrix$barcode, "-1")
    matrix <- as.matrix(matrix)
    rownames(matrix) <- NULL
    rownames(matrix) <- matrix[,1]
    matrix <- matrix[,2:ncol(matrix)]
    #transpose matrix and convert matrix to dgCMatrix class 
    matrix <- t(matrix)
    colnames <- colnames(matrix)
    rownames <- rownames(matrix)
    matrix <- Matrix::Matrix(as.numeric(matrix), nrow=length(rownames), ncol=length(colnames), sparse=TRUE, dimnames=list(rownames, colnames))
    
    #return
    return(list("taxid_list"=taxa_list, "matrix"=matrix))
  }
  
}
