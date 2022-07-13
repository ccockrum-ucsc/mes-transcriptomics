##### This script uses the ChIPSeeker package to assign TF ChIP peaks (from modENCODE and modERN) to genes

library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(org.Ce.eg.db)
library(tidyverse)
library(magrittr)
library(biomaRt)

select <- dplyr::select

# make ce10 Txdb object

ce10_ensembl_txdb <- makeTxDbFromEnsembl(organism = 'Caenorhabditis elegans', release = 66)
ce10_ensembl_txdb <- renameSeqlevels(ce10_ensembl_txdb, c('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX', 'chrM'))
save(ce10_ensembl_txdb, file = '~/stromelab/bioinformatics/omes/ws220/transcriptome/ce10_ensembl_txdb.rda')
seqlevels(ce10_ensembl_txdb)

##### modENCODE data #####

library(rtracklayer)

files = list.files('~/stromelab/bioinformatics/chipData/LIN-15B/peaks/')
gRangesList_modEncode <- lapply(files, function(x) {
  print(x)
  gff <- readGFF(filepath = paste0('~/stromelab/bioinformatics/chipData/LIN-15B/peaks/', x))
  gff %<>% mutate(seqid = paste0('chr',.$seqid)) 
  gff <- makeGRangesFromDataFrame(gff, starts.in.df.are.0based = TRUE)
  return(gff)
  
})
newName <-gsub("#.*","",files)
names(gRangesList_modEncode) <- files
names(gRangesList_modEncode) <- newName

#### load chIP peak bedfiles from modERN ####

files = list.files('~/stromelab/bioinformatics/chipData/modERN/modErn_chipPeakCalls_ce10/beds/', pattern = '.Peak' )
gRangesList <- lapply(files, function(x) {
  name <- sub('.gz','', x)
  dt <- read.table(paste0('~/stromelab/bioinformatics/chipData/modERN/modErn_chipPeakCalls_ce10/beds', '/', x))
  dt %<>% distinct() %>% dplyr::select(seqnames = 1, start = 2, end = 3, strand = 6)
  makeGRangesFromDataFrame(dt, starts.in.df.are.0based = TRUE)
})

fileNames <- files %>% map(~paste0(.x %>% str_extract('optimal.*?_') %>% str_remove('optimal.'),
                                   .x %>% str_remove('_IP.*') %>% str_remove('.*_')))
names(gRangesList) <- fileNames

gRangesList <- c(gRangesList, gRangesList_modEncode)

#### chipseeker to annotate genes ####

annotatedPeaks <- sapply(names(gRangesList), USE.NAMES = TRUE, simplify = FALSE, function(x) {
  dt <- gRangesList[[x]]
  annotatePeak(peak = dt, 
               TxDb = ce10_ensembl_txdb , verbose = F, level = 'transcript', overlap = 'TSS', tssRegion = c(-500, 0))
})

load('~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_wide.rda')

genesWithPeaks <- sapply(names(annotatedPeaks), simplify = FALSE, USE.NAMES = TRUE, function(x) {
  dt <- as.data.frame(annotatedPeaks[[x]])
  genes <- dt %>% 
    dplyr::select(geneId)
  genes <- pgc_wide %>% filter(Gene.Sequence.Name %in% genes$geneId | Gene.Public.Name %in% genes$geneId) %>% pluck('Gene.WB.ID')
  write.table(genes, paste0('~/stromelab/bioinformatics/chipData/modERN/modErn_chipPeakCalls_ce10/genes/', x, '.txt'),
              sep = '/t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(genes)
})


genesWithPeaks <- genesWithPeaks[!duplicated(genesWithPeaks)]
genesWithPeaks <- 
  genesWithPeaks %>% 
  map(tibble) %>% 
  map2(names(genesWithPeaks), ~ mutate(.x, 'tf' = .y)) %>%
  bind_rows()%>% 
  separate(tf, into = c('tf', 'stage'), sep = '_') %>% 
  select(-stage)%>% 
  group_by(tf) 

genesWithPeaks <- split(genesWithPeaks$`<chr>`, genesWithPeaks$tf)
lin15bGenes <- genesWithPeaks[grep('LIN-15B', names(genesWithPeaks))] %>% flatten_chr() %>% unique()

save(lin15bGenes, file ='~/stromelab/pgcAnalysis/lin15bGenes.rda')

