# pyclone pipeline

## Params
seed=123123
restarts=10
chrprefix="Y"

## Get args
args <- commandArgs(trailingOnly=TRUE)
file_list <- args[1]

## Set outdir if needed
if(!is.na(args[2])){
  out_dir <- args[2]
} else {
  out_dir <- ""
}

## arg check
if(!file.exists(args[1])){
	stop("input file not found")
}

if(out_dir != ""){
	if(!dir.exists(out_dir)){
		stop("output dir not found")
	}
}	

## libs
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

## Functions
files <- read.table(file_list,header=T,sep="\t")

pat_split <- split(files,f=files$patient_id)
	
read_cn_data <- function(x){
  grange_tables <- lapply(x,function(y){
    if(nrow(y) == 1){
      tabs <- read.table(y$cn_file[1],header=T,sep="\t")
      tabs$sample <- y$sample_id
      tabs <- list(makeGRangesFromDataFrame(tabs,keep.extra.columns=T,start.field="startpos",end.field="endpos"))
    } else {
      for(i in 1:nrow(y)){
        if(i == 1){
          tabs <- read.table(y$cn_file[i],header=T,sep="\t")
          tabs$sample <- y$sample_id[i]
          tabs <- list(makeGRangesFromDataFrame(tabs,keep.extra.columns=T,start.field="startpos",end.field="endpos"))
        } else {
          tab_new <- read.table(y$cn_file[i],header=T,sep="\t")
          tab_new$sample <- y$sample_id[i]
          tab_new <- makeGRangesFromDataFrame(tab_new,keep.extra.columns=T,start.field="startpos",end.field="endpos")
          tabs <- append(tabs,list(tab_new))
        }
      }
    }
    names(tabs) <- y$sample_id
    return(tabs)
  })
  return(grange_tables)
}

read_snv_data <- function(x){
  snv_tables <- lapply(x,function(y){
    if(nrow(y) == 1){
      tabs <- read.csv(y$snv_file[1],header=T,sep=",") %>%
            select(Chromosome,Position,Ref,Alt,Depth_tumour,Depth_normal,Tumour)
      tabs$Tumour <- y$sample_id
      tabs <- list(makeGRangesFromDataFrame(tabs,keep.extra.columns=T,
                                                         start.field="position",
                                                         end.field="position"))
    } else {
      for(i in 1:nrow(y)){
        if(i == 1){
          tabs <- read.csv(y$snv_file[i],header=T,sep=",") %>%
            select(Chromosome,Position,Ref,Alt,Depth_tumour,Depth_normal,Tumour)
          tabs$Tumour <- y$sample_id[i]
          tabs <- list(makeGRangesFromDataFrame(tabs,keep.extra.columns=T,
                                                           start.field="position",
                                                           end.field="position"))
        } else {
          tab_new <- read.csv(y$snv_file[i],header=T,sep=",") %>%
            select(Chromosome,Position,Ref,Alt,Depth_tumour,Depth_normal,Tumour)
          tab_new$Tumour <- y$sample_id[i]
          tab_new <- makeGRangesFromDataFrame(tab_new,keep.extra.columns=T,
                                              start.field="position",
                                              end.field="position")
          tabs <- append(tabs,list(tab_new))
        }
      }
    }
    names(tabs) <- y$sample_id
    return(tabs)
  })
  return(snv_tables)
}

findCNoverlaps <- function(snv_grange,cn_grange){
  overlaps <- GenomicRanges::findOverlaps(snv_grange,cn_grange)
  ## select overlapping snvs
  snv_grange_filt <- snv_grange[queryHits(overlaps),]
  ## Recompute overlaps
  overlaps <- findOverlaps(snv_grange_filt,cn_grange)
  mcols(snv_grange_filt)$major_cn <- mcols(cn_grange[subjectHits(overlaps)])$nMajor
  mcols(snv_grange_filt)$minor_cn <- mcols(cn_grange[subjectHits(overlaps)])$nMinor
  
  snv_grange_filt.df <- as.data.frame(snv_grange_filt) %>%
    mutate(normal_cn = ifelse(seqnames != chrprefix,2,1)) %>%
    mutate(mutation_id = paste(seqnames,start,Ref,Alt,sep="_")) %>%
    dplyr::rename("ref_counts" = "Depth_normal","alt_counts" = "Depth_tumour","sample_id" = "Tumour") %>%
    select(c("mutation_id","sample_id","ref_counts","alt_counts","major_cn","minor_cn","normal_cn")) %>%
    mutate(tumour_content = getPurityVal(.))
  
  sample_name <- unique(snv_grange_filt.df$sample_id)
  
  return(snv_grange_filt.df)
}

getPurityVal <- function(x){
  pu <- files$purity[match(x$sample_id,files$sample_id)]
  return(pu)
}

getCNoverlaps <- function(snv,cn){
  ov <- mapply(snv,cn,SIMPLIFY = F,FUN = function(a,b){
    if(any(names(a) != names(b))){
      stop("lists out of order")
    }
    t <- do.call(rbind,mapply(a,b,USE.NAMES = T,SIMPLIFY = F,FUN = function(x,y){
      findCNoverlaps(x,y)
    }))
  })
  lapply(names(ov),function(x){
    sample_name <- x
    write.table(x = ov[[x]],
                file = paste0(out_dir,sample_name,"_pyclone_input.tsv"),
                quote = F,append = F,sep = "\t",row.names = F,col.names = T)
  })
  return(ov)
}

callPyclone <- function(x){
  lapply(names(x),function(y){
    sample_name <- y
    pyclone_fit <- paste0("pyclone-vi fit -i ",
                        paste0(out_dir,sample_name,"_pyclone_input.tsv")," -o ",
                        paste0(out_dir,sample_name,"_pyclone_fit.h5d")," --seed ",seed," -r ",restarts)
    system(pyclone_fit)
    pyclone_write <- paste0("pyclone-vi write-results-file -i ",
                          paste0(out_dir,sample_name,"_pyclone_fit.h5d")," -o ",
                          paste0(out_dir,sample_name,"_pyclone_results.tsv"))
    system(pyclone_write)
  })
  return(cat("pyclone finished"))
}

formatOutput <- function(x){
  lapply(names(x),function(y){
    sample_name <- y
    pyclone_results <- read.table(paste0(out_dir,sample_name,"_pyclone_results.tsv"),header = T,sep = "\t")
    pyclone_results <- pyclone_results %>%
                        pivot_wider(id_cols = mutation_id,
                                    names_from = sample_id,
                                    values_from = -c(1,2))
    snv_list_file <- pat_split[[y]]$snv_file[1]
    snv_list <- read.table(snv_list_file,header=T,sep=",",row.names=1)
    snv_list$mutation_id <- paste(snv_list$Chromosome,snv_list$Position,snv_list$Ref,snv_list$Alt,sep="_")
    joined <- dplyr::left_join(pyclone_results,snv_list,by="mutation_id")
    write.table(x = joined,file = paste0(out_dir,sample_name,"_pyclone_results_final.tsv"),quote = F,append = F,sep = "\t",row.names = F,col.names = T)
  })
}

## Run pipeline
cn <- read_cn_data(pat_split)
snvs <- read_snv_data(pat_split)
pyclone_input <- getCNoverlaps(snvs,cn)
callPyclone(pyclone_input)
formatOutput(pyclone_input)

# END
