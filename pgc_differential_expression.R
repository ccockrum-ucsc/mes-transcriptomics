library(scran)
library(scater)
library(SingleCellExperiment)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(purrr)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(ggpubr)
library(colorspace)


load(file = '~/stromelab/bioinformatics/geneSets/Chad_new.adapted.tables/master.wide.table.rda')

select <- dplyr::select

source('~/stromelab/bioinformatics/R/ChadsFunctions.R')
dir <- "~/stromelab/illuminaSequencing/SSCC005/counts/"
dir2 <- "~/stromelab/illuminaSequencing/SSCC004/counts2/"

##### import counts #####

sscc005_counts <- importCounts(dir)
sscc004_counts <- importCounts(dir2)

all_counts <- left_join(sscc005_counts, sscc004_counts, by = "Geneid") %>%
  select(order(colnames(.)))

all_counts %<>% column_to_rownames('Geneid')
all_counts %<>% 
  mutate(across(everything(), as.numeric)) 

##### prepare annotation data #####

anno <- read.table("~/stromelab/pgcAnalysis/pgcProfiling_R/combined_anno.txt", sep="\t",header=TRUE)
anno <- anno %>% dplyr::mutate(batch = str_c('batch',batch)) %>%
  dplyr::mutate(batch = as.factor(batch), condition = as.factor(condition), platform = as.factor(platform)) %>%
  dplyr::arrange(cell)
rownames(anno) <- anno$cell


### Create Single Cell Experiment Object

sce <- SingleCellExperiment(assays =  list(counts = all_counts),
                            colData = anno)

spikes <- grepl("^ERCC", rownames(sce))
feat_type <- ifelse(grepl("^ERCC", rownames(sce)), 'spikes', 'endogenous')
splitAltExps(sce, feat_type)

### Sample QC stats ###

per_cell <- perCellQCMetrics(sce, subsets = list(ERCC=grep("^ERCC", rownames(sce))))
per_feat <- perFeatureQCMetrics(sce)

### Calculate logcounts_raw slot with just counts transformed and the logcounts slot that is adjusted for library size.

logcounts(sce) <- log2(calculateCPM(sce) + 1)
assay(sce, "logcounts_raw") <- log2(counts(sce) + 1)

### Sample QC stats ###

qc_stats <- quickPerCellQC(per_cell, percent_subsets="subsets_ERCC_percent")

### Filter low quality samples ###

sce <- sce[,!qc_stats$discard]

### Filter mrg1_l1pgcs_rep05 because it has wild-type reads and therefore isn't mrg-1 null mutant ###

sce$use <- (
  sce$cell != "mrg1_l1pgcs_rep05" 
)

### make table of rawcounts ###

rawcounts <- assays(sce)$counts %>% as.data.frame() %>% rownames_to_column('gene')
write.table(rawcounts, file = "~/stromelab/pgcs_egcs_rawcounts.txt", row.names = F, quote = F, sep = '\t')


##### Make filtered sce #####

sce_filt <- sce[, colData(sce)$use]
keep_feature <- nexprs(sce_filt, byrow=TRUE, detection_limit = 0) > 0
sce_filt <- sce_filt[keep_feature,]

##### calculate normalization scaling factors using scran #####

clust <- quickCluster(sce_filt, min.size = 1)
table(clust)
sce_norm <- scran::computeSumFactors(sce_filt, clusters = clust)

##### DESeq2 #####

level_list <- list(
  'mrg1L1_wtL1' = c('mrg1','wt'),
  'mes4L1_wtL1' = c('mes4','wt'),
  'mes3L1_wtL1' = c('mes3','wt'),
  'met1L1_wtL1' = c('met1','wt'),
  'drh3L1_wtL1' = c('drh3','wt'),
  'mes4L2_wtL2' = c('mes4_l2germ','wt_l2germ'),
  'mes3L2_wtL2' = c('mes3_l2germ','wt_l2germ'),
  'wtL2_wtL1' = c('wt_l2germ','wt'),
  'mes4L2_mes4L1' = c('mes4_l2germ','mes4'),
  'mes3L2_mes3L1' = c('mes3_l2germ','mes3'),
  'mes4L2_mes3L2' = c('mes4_l2germ','mes3_l2germ'),
  'mes4L1_mes3L1' = c('mes4','mes3')
)

deseq_results <- deseq2_pair(sce = sce_norm, level_list = level_list, cooksCutoff = F, independentFiltering = T, minmu = 1e-6)

deseq_results <- map(deseq_results, function(x) {
  res  <- lfcShrink(dds = x[[1]], coef = 2, res = x[[2]], type = 'ashr')
  return(list('.dds' = x[[1]], '.res' = res))
})

dds_list <- lapply(deseq_results, "[[", 1)
res_list <- lapply(deseq_results, "[[", 2)

lapply(res_list, function(x) {
  plotMA(x,
         ylim = c(-15,15))
  hist(x$pvalue, breaks = 20)
  summary(x)
})

lapply(dds_list, function(x) {
  DESeq2::plotDispEsts(x)
})

pgc_deseq2_resultsTable <- combineResults(deseq_results)
save(pgc_deseq2_resultsTable, file = '~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_deseq2_resultsTable.rda')


##### normalize sce #####

##### vst transform dds.bulk that includes all samples.

dds.bulk <- convertTo(sce_norm, type = 'DESeq2')
rld <- vst(dds.bulk, blind = T)
plotPCA(rld)
save(dds.bulk, file = "~/stromelab/pgcAnalysis/pgcProfiling_R/dds.bulk_R.rda")

pgc_normCounts <- counts(dds.bulk, normalized = TRUE) %>% data.frame()

y <- colnames(pgc_normCounts) %>% as.character()
y <- str_remove(y, '_rep..') %>% unique()
for(level in y) {
  pgc_normCounts <- pgc_normCounts %>% mutate(!!paste0({{level}},'.average.counts') := rowMeans(select(pgc_normCounts,starts_with({{level}}))))
}
pgc_normCounts <- pgc_normCounts %>% mutate('Gene.Sequence.Name' = rownames(sce_norm)) %>% select(Gene.Sequence.Name, ends_with('counts'))
save(pgc_normCounts, file = '~/stromeLab/pgcAnalysis/pgcProfiling_R/pgc_normCounts.rda')


##### add normalized counts to deseq2_resultsTable #####

pgc_deseq2_resultsTable %<>% right_join(pgc_normCounts, by = 'Gene.Sequence.Name')
save(pgc_deseq2_resultsTable, file = '~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_deseq2_resultsTable.rda')


##### Merge table to master table ######

load(file = '~/stromelab/bioinformatics/geneSets/Chad_new.adapted.tables/master.wide.table.rda')
load(file = '~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_deseq2_resultsTable.rda')

pgc_wide <- master.wide.table %>% left_join(pgc_deseq2_resultsTable, by = 'Gene.Sequence.Name')
pgc_wide <- pgc_wide %>% filter(Transcript.Type == 'coding')  

save(pgc_wide, file = '~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_wide.rda')


##### make tidy table #####

pgc_tidy <- pgc_wide %>%
  mutate(identifier = 1:nrow(.)) %>%
  gather("comparison", "Value", 
         contains(names(level_list)),
         factor_key = TRUE, convert = TRUE) %>%
  separate("comparison", c('comparison','Metric'), sep = "\\.", convert = TRUE) %>%
  spread('Metric', Value, convert = TRUE) %>%
  mutate('change' = case_when(
    log2FoldChange > 0 & padj < .05 ~ 'up',
    log2FoldChange < 0 & padj < .05 ~ 'down',
    padj >= 0.05 ~ 'noChange',
    is.na(padj) ~ 'noTest'
  )) %>% 
  mutate(change = factor(change, levels = c('noChange', 'noTest' ,'down','up')))
save(pgc_tidy, file = '~/stromelab/pgcAnalysis/pgcProfiling_R/pgc_tidy.rda')


##### plotPCA #####

qualcol <- qualitative_hcl(n = 6, c = 100, l = 50)

rld <- vst(dds.bulk[,!dds.bulk$condition %in% c('drh3')])
pcaDt <- plotPCA(rld, returnData = T, intgroup = c('condition', 'genotype', 'stage'))
pc1 <- attributes(pcaDt)$percentVar[1]*100
pc2 <- attributes(pcaDt)$percentVar[2]*100

pcaPl <-
  pcaDt %>% ggplot(aes(x = PC1, y = PC2, color = genotype, shape = stage)) + geom_point(size = 2) + coord_equal() + theme_pubr(legend = 'top', border = T) +
  scale_x_continuous(name = paste0('PC1: ', sprintf("%.0f", pc1), '% variance')) +
  scale_y_continuous(name = paste0('PC2: ', sprintf("%.0f", pc2), '% variance')) + 
  scale_color_brewer(palette = 'Dark2' )
pcaPl
ggsave(plot = pcaPl, filename = '~/stromelab/pgcAnalysis/pca/pgc_all_pcaPlot.tiff', width = 8, height = 8, units = 'cm')


#### make gene sets ####

# load RNA-seq files from adult germline sequencing

load('~/stromelab/mesAdultGermlineSeq/ad_bkgList.rda')
load('~/stromelab/pgcAnalysis/lin15bGenes.rda')

pgc_bkgList <- list(
  "allGenes" = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1'),
  "X" = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', X.vs.Autosome == 'X'),
  "Autosome" = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', X.vs.Autosome == 'Autosome'),
  'ubiq' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', gs.ubiq.SAGE == 1),
  'germSpec' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', gs.germl.spec == 1),
  'germEnr' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', gs.germl.herm == 1),
  'somaSpec' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', gs.soma.spec.any == 1, !Gene.WB.ID %in% ad_bkgList$adGermExp),
  'oogenicEnr_ortiz' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', Ortiz.Oogenic == 1),
  'spermEnr_ortiz' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', Ortiz.Spermatogenic == 1),
  'genderNeutral_ortiz' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', Ortiz.GenderNeutral == 1),
  'l1Exp' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', wt_l1pgcs.average.counts >= 1 ),
  'l2Exp' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', wt_l2germ.average.counts >= 1 ),
  'l1l2Exp' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', wt_l2germ.average.counts >= 1 | wt_l1pgcs.average.counts >= 1),
  'lin15b' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', lin15bGenes == 1),
  'mes4_bound' = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', mes4Bound == 1)
) %>% purrr::map(select, 'Gene.WB.ID') %>% purrr::map(flatten_chr)

pgc_deList <- list(
  "mes4L1_wtL1_up" = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', change == 'up'), 
  "mes3L1_wtL1_up" = pgc_tidy %>% filter(comparison == 'mes3L1_wtL1', change == 'up'), 
  "mrg1L1_wtL1_up" = pgc_tidy %>% filter(comparison == 'mrg1L1_wtL1', change == 'up'),
  "met1L1_wtL1_up" = pgc_tidy %>% filter(comparison == 'met1L1_wtL1', change == 'up'),
  "drh3L1_wtL1_up" = pgc_tidy %>% filter(comparison == 'drh3L1_wtL1', change == 'up'),
  "mes4L2_wtL2_up" = pgc_tidy %>% filter(comparison == 'mes4L2_wtL2', change == 'up'), 
  "mes3L2_wtL2_up" = pgc_tidy %>% filter(comparison == 'mes3L2_wtL2', change == 'up'), 
  "wtL2_wtL1_up" = pgc_tidy %>% filter(comparison == 'wtL2_wtL1', change == 'up'), 
  "mes4L1_mes3L1_up" = pgc_tidy %>% filter(comparison == 'mes4L1_mes3L1', change == 'up'), 
  "mes4L2_mes3L2_up" = pgc_tidy %>% filter(comparison == 'mes4L2_mes3L2', change == 'up'), 
  "mes4L1_wtL1_down" = pgc_tidy %>% filter(comparison == 'mes4L1_wtL1', change == 'down'),
  "mes3L1_wtL1_down" = pgc_tidy %>% filter(comparison == 'mes3L1_wtL1', change == 'down'),
  "mrg1L1_wtL1_down" = pgc_tidy %>% filter(comparison == 'mrg1L1_wtL1', change == 'down'),
  "met1L1_wtL1_down" = pgc_tidy %>% filter(comparison == 'met1L1_wtL1', change == 'down'),
  "drh3L1_wtL1_down" = pgc_tidy %>% filter(comparison == 'drh3L1_wtL1', change == 'down'),
  "mes4L2_wtL2_down" = pgc_tidy %>% filter(comparison == 'mes4L2_wtL2', change == 'down'),
  "mes3L2_wtL2_down" = pgc_tidy %>% filter(comparison == 'mes3L2_wtL2', change == 'down'),
  "wtL2_wtL1_down" = pgc_tidy %>% filter(comparison == 'wtL2_wtL1', change == 'down')
) %>% purrr::map(select, 'Gene.WB.ID') %>% purrr::map(flatten_chr)

pgc_deList <- c(
  pgc_deList, 
  list(
    'mes3_up' = union(pgc_deList$mes3L1_wtL1_up, pgc_deList$mes3L2_wtL2_up),
    'mes4_up' = union(pgc_deList$mes4L1_wtL1_up, pgc_deList$mes4L2_wtL2_up),
    'mes3_down' = union(pgc_deList$mes3L1_wtL1_down, pgc_deList$mes3L1_wtL1_down),
    'mes4_down' = union(pgc_deList$mes4L1_wtL1_down, pgc_deList$mes4L2_wtL2_down)))


save(pgc_deList, file = '~/stromelab/pgcAnalysis/diffExpGenes/pgc_deList.rda' )
save(pgc_bkgList, file = '~/stromelab/pgcAnalysis/diffExpGenes/pgc_bkgList.rda' )

makeBEDfile(pgc_bkgList, dir = '~/stromelab/pgcAnalysis/diffExpGenes/beds/wholeGene/', chrName = T)
makeBEDfile(pgc_deList, dir = '~/stromelab/pgcAnalysis/diffExpGenes/beds/wholeGene/', chrName = T)

map2(pgc_deList, names(pgc_deList), ~saveId(.x, .y, dir = '~/stromelab/pgcAnalysis/diffExpGenes/wormbaseIds/', id = 'Gene.WB.ID'))
map2(pgc_bkgList, names(pgc_bkgList), ~saveId(.x, .y, dir = '~/stromelab/pgcAnalysis/diffExpGenes/wormbaseIds/', id = 'Gene.WB.ID'))







