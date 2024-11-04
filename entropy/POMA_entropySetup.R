library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(MASS)
library(LEA)

setwd('~/gh/POMA/entropy/')

makePopId <- function(fileIndv){
  PopIDdf = read.table(fileIndv, sep="\t") %>%
    as.data.frame() %>%
    rename(All = V1) %>%
    mutate(Population = str_split_i(All, '_', 2),
           ID = str_split_i(All, '_', 3))
  return(PopIDdf)
}

PCA_entropy <- function(g){
  colmean = apply(g, 2, mean, na.rm = T)
  normalize = matrix(nrow = nrow(g), ncol = ncol(g))
  af = colmean/2
  for (m in 1:length(af)){
    nr = g[,m]-colmean[m]
    dn = sqrt(af[m]*(1-af[m]))
    normalize[,m] = nr/dn
  }
  normalize[is.na(normalize)] = 0
  method1 = prcomp(normalize, scale. = F,center = F)
  pca_df = method1$x[,1:27]
  return(pca_df)
}

PopID <- makePopId('POMA.entropy1.012.indv')

g <- read.table('pntest_mean_POMA.bithin.q30.i70.maf02.a90.recode.txt', header = F)

pca_df <- PCA_entropy(t(g)) %>%
  .[,1:10] %>%
  cbind(PopID)

writeLDAfile <- function(pcaDF, k){
  kCluster = kmeans(pcaDF[,1:5], k, iter.max = 10, nstart = 10, algorithm = 'Hartigan-Wong')
  ldaOut = lda(x = pcaDF[,1:5], grouping = kCluster$cluster, CV = T)
  write.table(round(ldaOut$posterior, 5),
              file = paste('ldak', as.character(k), '.txt', sep = ''),
              quote = F, row.names = F, col.names = F)
}

writeLDAfile(pca_df, 2)
writeLDAfile(pca_df, 3)
writeLDAfile(pca_df, 4)
writeLDAfile(pca_df, 5)
writeLDAfile(pca_df, 6)
writeLDAfile(pca_df, 7)

PopID_list <- paste(PopID$Pop, PopID$ID, sep = '_')

header <- data.frame(dims = NA, PopID_list)

df <- t(header)
dims <- paste(dim(g)[2], dim(g)[1], sep = " ")

df[1,1] <- dims

write.table(df, 'entropy_header.txt',
            sep = " ", na = "",
            quote = FALSE, row.names = FALSE, col.names = FALSE)