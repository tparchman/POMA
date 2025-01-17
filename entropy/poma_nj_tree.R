######################################
## poma_nj_tree.R by JPJ 17i25
######################################

## PURPOSE: to create a neighbor joining tree plot
## USAGE: Rscript poma_nj_tree.R


## load data

gprobs <- read.csv("gprob2.txt", header=TRUE)

gprobs_noname <- gprobs[,-1]

ids <- read.csv("POMA_PopID.csv", header=TRUE)


## calculate allele frequencies

uniq_pops <- sort(unique(ids[,3]))
afreqs <- matrix(0, length(uniq_pops), dim(gprobs_noname)[2])
for (i in 1:length(uniq_pops))
	{
	sub_dat <- subset(gprobs_noname, ids[,3]==as.character(uniq_pops[i]))
	for (j in 1:dim(gprobs_noname)[2])
		{
		af <- mean(sub_dat[,j]) / 2
		afreqs[i,j] <- af
		}
	}



## calculate fst

hudsonFst2<-function(p1=NA, p2=NA, n1=NA, n2=NA)
	{
    numerator<- p1 * (1 - p1) + p2 * (1 - p2)
    denominator<- p1 * (1 - p2) + p2 * (1 - p1)
    fst<- 1 - numerator/denominator
    out <- cbind(numerator, denominator, fst)
    return(out)
	}

fst_mean_out <- matrix(0, length(uniq_pops), length(uniq_pops))

for (i in 1:length(uniq_pops))
	{
	for (j in 1:length(uniq_pops))
		{
		if (i==j) { fst_mean_out[i,j] <- 0 }
		else if (i>j)
			{
			locus_vector <- vector()
			for (k in 1:dim(gprobs_noname)[2])
				{
				if 		(afreqs[i,k]==0 && afreqs[j,k]==0) { locus_vector <- append(locus_vector, 0) }
				else if (afreqs[i,k]==1 && afreqs[j,k]==1) { locus_vector <- append(locus_vector, 0) }
				else
					{
					loc_fst <- hudsonFst2(p1=afreqs[i,k], p2=afreqs[j,k])
					locus_vector <- append(locus_vector, loc_fst[3])
					}
				}
			mean_fst <- mean(locus_vector)
			fst_mean_out[i,j] <- mean_fst
			fst_mean_out[j,i] <- mean_fst
			}
		}
	}


## nj tree

library(ape)
library(MetBrewer)
nj_cols <- met.brewer("Archambault", 7)
tr <- nj(as.dist(fst_mean_out))
pdf("poma_nj_tree.pdf", height=6, width=6)
#quartz(height=6, width=6)
par(mar=c(0,0,0,0))
plot(tr, type="unrooted", use.edge.length=TRUE, edge.width=2, show.tip.label=FALSE)
tiplabels(bg=nj_cols, cex=2.5, pch=21)
legend("topright", legend=uniq_pops, pt.bg=nj_cols, pch=22, pt.cex=3, cex=1.5)
dev.off()












