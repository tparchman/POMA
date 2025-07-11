################################################
## poma_div.R by JPJ 11 vii 25
################################################

## load files
gprobs <- read.csv("gprob2.txt", header=TRUE)
	dim(gprobs)
	gprobs[1:10,1:10]
	
gprobs_noname <- gprobs[,-1]
	dim(gprobs_noname)
	gprobs_noname[1:10,1:10]

ids <- read.csv("POMA_PopID.csv", header=TRUE)
	dim(ids)
	head(ids)

pop_list <- sort(unique(ids[,3]))


## functions

calc_he <- function(input_matrix=NA){
	he_vect <- vector(length=dim(input_matrix)[2])
	for (i in 1:dim(input_matrix)[2]){
		afreq <- mean(input_matrix[,i]) / 2		## calculate the locus allele frequency
		he_vect[i] <- 2 * afreq * (1-afreq)		## calculate expected heterozygosity from the allele frequency
	}
	return(mean(he_vect))
}

calc_ho <- function(input_matrix=NA){
	ho_vect <- vector(length=dim(input_matrix)[2])
	for (i in 1:dim(input_matrix)[2]){
		ho_vect[i] <- mean(input_matrix[,i] > 0.9 & input_matrix[,i] < 1.1)
	}
	return(mean(ho_vect))
}


## analyses

div_mat <- matrix(NA, length(pop_list), 4)
colnames(div_mat) <- c("pop", "he", "ho", "fis")
for (i in 1:length(pop_list)){
	div_mat[i,1] <- pop_list[i]
	pop_sub <- subset(gprobs_noname, ids[,3]==pop_list[i])
	div_mat[i,2] <- calc_he(pop_sub)
	div_mat[i,3] <- calc_ho(pop_sub)
	div_mat[i,4] <- 1 - (as.numeric(div_mat[i,3])/as.numeric(div_mat[i,2]))
	print(i)
}


div_mat
write.table(div_mat, file="poma_div_mat.txt", row.names=F, col.names=F, quote=F)

