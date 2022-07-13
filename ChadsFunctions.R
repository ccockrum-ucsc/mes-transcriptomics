# Importing counts

importCounts<- function(dir, filterAttribute = F, attribute = 'protein_coding', attribute_name = 'gene_biotype') {
  files <- list.files(dir)
  tmp <- lapply(files, function(x) {
    
    name <- sub(".txt", "", x) %>%
      sub("-", ".", .)
    data <- read.table(paste0(dir,'/',x), header = T, skip = 1)
    
      data <- data %>% select(Geneid, !!name := 8, {{attribute_name}})
      if(filterAttribute == T) {
        print(attribute)
        data %<>% filter(gene_biotype == 'protein_coding') 
      }
      data %<>% select(-{{attribute_name}})
    })
  all.counts <- tmp %>%
    purrr::reduce(left_join,by='Geneid') %>%
    select(order(colnames(.)))
  return(all.counts)
}
  
#### convert counts to TPM ####

countsToTPM <- function(x, lengths) {
  transcripts <- x/( {{ lengths }}/1000)
  libsize <- sum(transcripts, na.rm = T)/1e6
  TPM <- transcripts/libsize
  return(TPM)
}
  
##### Make BED file #####

makeBEDfile <- function(my.list, dir, annoTable = master.wide.table, chrName = TRUE) {
  bed.list <- lapply(seq_along(my.list), function(i, ...) {
    name <- names(my.list)[[i]]
    data <- my.list[[i]]
    data <- annoTable %>% filter(Gene.WB.ID %in% data)
    bed_file <- data.frame(
      'chrom'= data$Chr,
      'chromStart'= ifelse(startsWith(data$Operon,"C"),data$Operon.Start, data$Start), 
      'chromEnd'=ifelse(startsWith(data$Operon,"C"),data$Operon.End, data$End),
      'name'=data$Gene.WB.ID,
      'score' = '0',
      'strand'=data$Strand
    )
    bed_file %<>% distinct(.keep_all = TRUE) %>% na.omit()
    if(chrName == FALSE) { 
          bed_file <-  mutate(bed_file, chrom = gsub('chr','', bed_file$chrom))
    }

    write.table(bed_file, paste0(dir,name,".BED",sep=""),row.names = FALSE, quote = FALSE, col.names = FALSE,sep='\t')
    return(bed_file)
  })
  names(bed.list) <- names(my.list)
  return(bed.list)
}
  
#### convert BED to GRanges ####

bedToGRanges <- function(bedFile) {
  bedFile %<>% distinct() %>% dplyr::select(seqnames = 1, start = 2, end = 3, strand = 6)
  GenomicRanges::makeGRangesFromDataFrame(bedFile, starts.in.df.are.0based = TRUE)
}

##### Extract wormbase IDs from table(s) for GO analysis #####

saveId = function(x, y, dir, id) {
  result <- master.wide.table %>% filter(Gene.WB.ID %in% x) %>% select({{id}})
  write.table(result, paste0(dir,y,".txt",sep=""),row.names = FALSE, quote = FALSE, col.names = FALSE,sep='\t')
  return(result)
  
}

##### DESeq2 functions #####

deseq2_pair <- function(sce, level_list, lfcThreshold = 0, alpha = 0.05, cooksCutoff = F, independentFiltering, minmu = 0.5, minReplicatesForReplace = Inf) {
  
  result_list <- lapply(seq_along(level_list), function(i) {
    name <- names(level_list[i])
    print(name)
    sce_data = sce[,sce$condition %in% level_list[[i]]]
    test_lvl <- level_list[[i]][1]
    control_lvl <- level_list[[i]][2]
    sce_data$condition <- droplevels(sce_data$condition)
    sce_data$batch <- droplevels(sce_data$batch)
    sce_data$condition <- relevel(sce_data$condition, control_lvl)
    
    dds <- convertTo(sce_data, type = 'DESeq2')
    design(dds) <- ~1 + condition
    
    dds <- DESeq(dds, fitType='parametric', useT = F, minReplicatesForReplace = Inf)
    print(resultsNames(dds))
    res <- results(dds, lfcThreshold = lfcThreshold, 
                   contrast = c("condition", test_lvl, control_lvl),
                   alpha = alpha, cooksCutoff = F, minmu = minmu,
                   independentFiltering = independentFiltering)
    x <- list('.dds' = dds, 
              '.res' = res)
  })
  
  names(result_list) <- names(level_list)
  return(result_list)
}
  
combineResults <- function(dds.list) {
  res.list <- dds.list %>% purrr::map(2)
  tmp.list <- lapply(seq_along(res.list), function(i) {
    basename <- names(res.list[i])
    res <- data.frame(res.list[[i]])
    res <- res %>%
      dplyr::rename(!!paste0(basename, '.log2FoldChange') := log2FoldChange,
                    !!paste0(basename, '.padj') := padj,
                    !!paste0(basename, '.baseMean') := baseMean) %>%
      mutate(Gene.Sequence.Name = rownames(res)) %>%
      select(Gene.Sequence.Name, everything(), -lfcSE, -pvalue)
    print(res)
    return(res)
  })
  purrr::reduce(tmp.list, full_join, by = 'Gene.Sequence.Name')
}

  
