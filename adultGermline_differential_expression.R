library(tidyverse)
library(magrittr)
library(dplyr)
library(DESeq2)

source('~/stromelab/bioinformatics/R/ChadsFunctions.R')

setwd('~/stromelab/mesAdultGermlineSeq')
counts <- importCounts('~/stromelab/mesAdultGermlineSeq/counts')
counts %<>% column_to_rownames('Geneid')
counts %<>% mutate(across(.fns = as.numeric))

anno <- read.table('~/stromelab/mesAdultGermlineSeq/20201223_colData.txt', header = T, sep = '\t', row.names = 1, stringsAsFactors = T)

##### deseq2 #####
 
 deseq2_pair_adults <- function(counts, colData, level.list, design = ~condition, lfcThreshold = 0, alpha = 0.05, cooksCutoff = F, independentFiltering = T, minmu = 0.5, useT = FALSE, minReplicatesForReplace = Inf) {

  result.list <- lapply(seq_along(level.list), function(i) {
    name <- names(level.list[i])
    print(name)
    test.lvl <- level.list[[i]][1]
    print(test.lvl)
    control.lvl <- level.list[[i]][2]
    print(control.lvl)

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design =  design)
    dds = dds[,dds$condition %in% level.list[[i]]]
    dds$condition <- droplevels(dds$condition)
    dds$condition <- relevel(dds$condition, control.lvl)

    dds <- DESeq(dds, useT = useT, minReplicatesForReplace = minReplicatesForReplace)
    res <- results(dds, lfcThreshold = lfcThreshold,
                   contrast = c("condition", test.lvl, control.lvl),
                   alpha = alpha, cooksCutoff = cooksCutoff, minmu = minmu,
                   independentFiltering = independentFiltering)

    res <- lfcShrink(dds, res = res, coef = 2, type = 'ashr', lfcThreshold = lfcThreshold)
    print(summary(res))
    x <- list('.dds' = dds,
              '.res' = res)
  })
  names(result.list) <- names(level.list)
  return(result.list)
}

level.list <- list(
  'mes4AdHerm_wtAdHerm' = c('mes4_herm','wt_herm'),
  'mes4AdMale_wtAdMale' = c('mes4_male','wt_male'),
  'mes4lin15bAdHerm_wtAdHerm' = c('mes4lin15b_herm','wt_herm_set2'),
  'mes4lin15bAdHerm_mes4AdHerm' = c('mes4lin15b_herm', 'mes4_herm'),
  'wtAdHerm2_wtAdHerm1' = c('wt_herm_set2', 'wt_herm'),
  'drh3AdHerm_wtAdHerm' = c('drh3_herm','wt_herm')
  )
  


deseq_results <- deseq2_pair_adults(counts = counts, colData = anno, level.list = level.list, lfcThreshold = 0)

dds_list <- lapply(deseq_results, "[[", 1)
res_list <- lapply(deseq_results, "[[", 2)

lapply(res_list, function(x) {
  plotMA(x,
         ylim = c(-5,5))
})

lapply(res_list, function(x) {
  hist(x$pvalue, breaks = 20)
})

lapply(dds_list, function(x) {
  DESeq2::plotDispEsts(x)
})

lapply(res_list_noX, function(x) {
  summary(x)
})

combinedDf <- combineResults(dds.list =  deseq_results)


##### get normalized counts #####

dds.bulk <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = anno, design =  ~condition)
dds.bulk <- dds.bulk[, dds.bulk$condition %in% c('wt_herm', 'wt_herm_set2','mes4_herm', 'mes4lin15b_herm')]
dds.bulk$condition <- relevel(dds.bulk$condition, 'wt_herm')
dds.bulk$condition <- droplevels(dds.bulk$condition)
dds.bulk <- DESeq(dds.bulk, minReplicatesForReplace = Inf)

normCounts <- data.frame(counts(dds.bulk, normalized = TRUE))
normCounts %<>% rename_all(~paste0(dds.bulk$condition, '_',dds.bulk$replicate))

y <- levels(colData(dds.bulk)$condition) %>% as.character()
for(level in y) {
  normCounts <- normCounts %>% mutate(!!paste0({{level}},'.average.counts') := rowMeans(dplyr::select(normCounts,starts_with({{level}}))))
}

normCounts <- normCounts %>% mutate('Gene.Sequence.Name' = rownames(rowData(dds.bulk))) %>% select(Gene.Sequence.Name, ends_with('counts'))


#### write table of raw counts ####

rawCounts <- data.frame(counts(dds.bulk, normalized = FALSE)) %>% rownames_to_column('gene')
rawCounts %<>% select(
                    gene,
                    'wt_adgerm_batch01_rep01' = CC1,
                    'wt_adgerm_batch01_rep02' = CC8,
                    'wt_adgerm_batch01_rep03' = CC15,
                    'mes4_adgerm_batch01_rep01' = CC3,
                    'mes4_adgerm_batch01_rep02' = CC10,
                    'mes4_adgerm_batch01_rep03' = CC17,
                     'mes4lin15b_adgerm_batch02_rep01' = mes4lin15b_adultGerm_herm_set2_rep01,
                     'mes4lin15b_adgerm_batch02_rep02' = mes4lin15b_adultGerm_herm_set2_rep02,
                     'mes4lin15b_adgerm_batch02_rep03' = mes4lin15b_adultGerm_herm_set2_rep03,
                     'wt_adgerm_batch02_rep01' = wt_adultGerm_herm_set2_rep01,
                     'wt_adgerm_batch02_rep02' = wt_adultGerm_herm_set2_rep02,
                     'wt_adgerm_batch02_rep03' = wt_adultGerm_herm_set2_rep03)
write.table(rawCounts, file = "~/stromelab/adgerm_rawcounts.txt", row.names = F, quote = F, sep = '\t')
pgcRawCount <- read.table('~/stromelab/pgcs_egcs_rawcounts.txt', sep = '\t', header = T)
combinedCount <- pgcRawCount %>% left_join(rawCounts, by = 'gene')
write.table(combinedCount, file = "~/stromelab/allsamples_rawcounts_matrix.txt", row.names = F, quote = F, sep = '\t')

##### add normalized counts to deseq2_resultsTable #####

ad_deseq2_resultsTable <- combinedDf %>% left_join(normCounts, by = 'Gene.Sequence.Name')
save(ad_deseq2_resultsTable, file = '~/stromelab/mesAdultGermlineSeq/ad_deseq2_resultsTable.rda')


##### Merge table to master table to make wide table ######

load(file = '~/stromelab/bioinformatics/geneSets/Chad_new.adapted.tables/master.wide.table.rda')
load(file = '~/stromelab/mesAdultGermlineSeq/ad_deseq2_resultsTable.rda')

ad_wide <- master.wide.table %>% left_join(ad_deseq2_resultsTable, by = 'Gene.Sequence.Name')
save(ad_wide, file = '~/stromelab/mesAdultGermlineSeq/ad_wide.rda')


#### make tidy table ####
  
ad_tidy <- ad_wide %>%
  mutate(identifier = 1:nrow(.)) %>%
  gather("comparison", "value", 
         contains(names(level.list)),
         
          factor_key = TRUE, convert = TRUE) %>%
   separate("comparison", c('comparison','metric'), sep = "\\.", convert = TRUE) %>%
   spread('metric', value, convert = TRUE) %>%
  mutate('change' = case_when(
    log2FoldChange > 0 & padj < .05 ~ 'up',
    log2FoldChange < 0 & padj < .05 ~ 'down',
    padj >= 0.05 ~ 'noChange',
    is.na(padj) ~ 'noTest'
  )) %>% 
  mutate(change = factor(change, levels = c('noChange', 'noTest' ,'down','up')))
save(ad_tidy, file = '~/stromelab/mesAdultGermlineSeq/ad_tidy.rda')


#### plot PCA ####

rld <- vst(dds.bulk)
pcaData <- DESeq2::plotPCA(rld, intgroup = 'condition', ntop = 10000, returnData = T)
assay(rld)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = condition), size = 2) +
  coord_equal() + theme_bw() +
  labs(x = paste0('PC1: ',percentVar[1],'%'),
       y = paste0('PC2: ',percentVar[2],'%'))

res <- results(dds.bulk, contrast = c('condition','mes4_maleXsp', 'wt_maleXsp'), alpha = 0.05, lfcThreshold = 0)
summary(res)


#### make gene sets ####

tpms <- ad_wide %>% 
  summarize(Gene.WB.ID, across('wt_herm.average.counts', ~ countsToTPM(.x, lengths = cur_data()$exonig.length))) %>%
  filter(wt_herm.average.counts > 5)


ad_tidy %<>% filter(Transcript.Type == 'coding')

ad_bkgList <-  ad_tidy %>% list(
  'X' = filter(., comparison == 'mes4AdHerm_wtAdHerm', X.vs.Autosome == 'X'),
  'Autosome'  = filter(., comparison == 'mes4AdHerm_wtAdHerm', X.vs.Autosome == 'Autosome'),
  'adGermExp' = tpms
  
) %>% map(~pluck(.x, 'Gene.WB.ID'))

ad_deList <- ad_tidy %>% list(
  'mes4Herm_wtHerm_up' = filter(., comparison == 'mes4AdHerm_wtAdHerm', change == 'up'),
  'mes4Male_wtMale_up' = filter(., comparison == 'mes4AdMale_wtAdMale', change == 'up'),
  'mes4lin15bHerm_wtHerm_up' = filter(., comparison == 'mes4lin15bAdHerm_wtAdHerm', change == 'up'),
  'mes4lin15bHerm_mes4Herm_up' = filter(., comparison == 'mes4lin15bAdHerm_mes4AdHerm', change == 'up'),
  'drh3AdHerm_wtAdHerm_up' = filter(., comparison == 'drh3AdHerm_wtAdHerm', change == 'up'),
  
  'mes4Herm_wtHerm_down' = filter(., comparison == 'mes4AdHerm_wtAdHerm', change == 'down'),
  'mes4lin15bHerm_wtHerm_down' = filter(., comparison == 'mes4lin15bAdHerm_wtAdHerm', change == 'down'),
  'mes4lin15bHerm_mes4Herm_down' = filter(., comparison == 'mes4lin15bAdHerm_mes4AdHerm', change == 'down'),
  'drh3AdHerm_wtAdHerm_down' = filter(., comparison == 'drh3AdHerm_wtAdHerm', change == 'down')
  
)%>% map(~pluck(.x, 'Gene.WB.ID'))

save(ad_deList, file = '~/stromelab/mesAdultGermlineSeq/ad_deList.rda' )
save(ad_bkgList, file = '~/stromelab/mesAdultGermlineSeq/ad_bkgList.rda' )


