div <- read.csv("poma_div_mat.csv")
div

# paired t test since these samples aren't independent
t.test(div$He,div$Ho,alternative="two.sided",paired=TRUE)

#### based on my excel hand calcs
0.02485099/(0.000196656*sqrt(7)) ### ~an order of magnitude off...?
