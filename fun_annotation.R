#!/usr/bin/env Rscript
#@author: Thomas NEFF

library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(dplyr)

getAnnotations <- function(specie.db=c("Homo sapiens", "EnsDb"),vars=c("gene_id", "gene_name", "gene_biotype", "seq_name", "description", "entrezid")){
  #' @param specie.db select specie's database of interest in Annotationhub db
  #' @param vars select the informations of interest
  #' 
  #' @return dataframe with annotation for each gene
  
  #call clone database, and update the copy cache on computer if necessary
  ah <- AnnotationHub()
  #select the last update
  ah.data <- ah[[tail(rownames(mcols(query(ah,pattern = specie.db, ignore.case = TRUE))),n=1)]]
  #get data with gene function
  ah.data.frame <- genes(ah.data,return.type = "data.frame")
  #select vars
  ah.data.frame <- ah.data.frame[,vars]
  
  #put the gene id in rownames, prepare for merge
  rownames(ah.data.frame) <- ah.data.frame[["gene_id"]]
  
  #use function in AnnotationHub() package
  #res <- genes(ah[[tail(rownames(mcols(query(ah,pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE))),n=1)]],return.type = "data.frame")[,c("gene_id", "gene_name", "gene_biotype", "seq_name", "description", "entrezid")]
    
  return(ah.data.frame)
}

mergeAnnotation <- function(annotations.obj,results.obj,default.col="id"){
  #' @param annotation.oobj dataframe with gene's information
  #' @param results.obj data.frame have to be annotated
  #' @param defaut.col default colname to merge datasets between annotation.obj and results.obj
  #' 
  #'@return res merge data.frame between annotation.obj and results.obj
  
  #define name on rownames column (it's safer because the rownames can miss after save or change var type)
  #Add robsutness to the code
  annotations.obj[[default.col]] <- rownames(annotations.obj)
  results.obj[[default.col]] <- rownames(results.obj)
  
  #check the transcript in DEGs results occording to getAnnotation() results
  #and indicate if intersections getAnnotation() results and DEGs results seems to be a problem
  if (all(annotations.obj$defaut.col %in% results.obj$default.col)) {
    #merge DEG result dataframe with annotation dataframe
    res <- merge(x = results.obj,y = annotations.obj, by= default.col)
  } else {
    print("Impossible to merge because se rownames are not similar (null intersection)" )
  }
  
  return(res)
}

addAnnotations <- function(results.obj,annotations.obj=getAnnotations(),merge.value = "id",gene_idx=FALSE){
  #' @param results.obj data to be annotated
  #' @param reference data for annotation
  #' @param merge.value name of column to be merged between results.obj and annotation obj
  #' @param gene_idx if annotation.obj are indexing with gene's symbol 
  #' 
  #' @return data annoteted with annotation extract to AnnotationHub()
  
  if (gene_idx) {
    #delete the duplicates lines in annotation dataframe
    annotations.obj <- annotations.obj %>% group_by(gene_name) %>% dplyr::filter(!duplicated(gene_name)) %>% tibble::column_to_rownames(var = "gene_name")
  }
  
  annotations.obj <- as.data.frame(annotations.obj)
  results.obj <- as.data.frame(results.obj)
  
  #check if there are a missing value between the annotation dataframe and the DEG result dataframe
  if (all(rownames(results.obj) %in% rownames(annotations.obj))) {
    #select and sort the rownames
    annotations.obj <- annotations.obj[rownames(results.obj),]
    #merge the dataframe
    res <- mergeAnnotation(annotations.obj,results.obj,default.col = merge.value)
    
    print("All dataset have been merged")
  } else {
    #select and sort the rownames
    annotations.obj <- annotations.obj[rownames(results.obj),]
    #merge the dataframe
    res <- mergeAnnotation(annotations.obj,results.obj,default.col = merge.value)
    
    print("Missing values in annotations object")
  }
  
  return(res)
}

convertCountsByGenes <- function(counts){
  #' @param counts counts matrix with with EnsID transcript indexing 
  #' 
  #' @return counts matrix with gene symbol indexing 

  #merge counts matrix with dataframe annotation
  tmp <- addAnnotations(results.obj = counts,gene_idx = FALSE)
  #select the transcripts associated with protein coding
  tmp <- tmp[which(tmp[,"gene_biotype"] == "protein_coding"),]
  #select colnames of genes names in dataframe for future rownames
  tmp <- tmp[c(colnames(counts),"gene_name")]
  #delete the missing data
  tmp <- tmp[!tmp$gene_name=="",]
  #summarize the counts matrix according to genes names
  new.counts <- tmp %>% group_by(gene_name) %>% summarise(across(everything(), sum)) %>% tibble::column_to_rownames(var = "gene_name")
  
  return(new.counts)
}

correct_ENSid <- function(counts){
  #' @param counts counts matrix idexing with EnsID transcript
  #' 
  #' @return counts matrix with EnsID corrected
  #' 
  #' warning: The Ens.ID are updated regurlarly and add to the following to cannonical transcript id
  #' example: cannonical transcript id.nb add -> ENSGXXXXX.XX, the second part are the number of the addition
  
  #delete the second part of ENS.ID
  counts["tmp"] <- sapply(rownames(counts), function(x) strsplit(x=as.character(x),split = "\\.")[[1]][1])
  #summarize the counts because of  transcipt with multiple names 
  new.counts <- counts %>% group_by(tmp) %>% summarise(across(everything(), sum)) %>% tibble::column_to_rownames(var = "tmp")
  
  return(new.counts)
}




