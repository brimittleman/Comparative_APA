load("../data/mediation_DF/negative_mediation.rda")
library(MASS)

# columns: genes, rows: random sample
mc_de <- list()
for (g in 1:length(is_de_neg)) {
  mu <- c(fit_de_neg$a[g], fit_de_neg$b[g])
  acov <- cbind(c(fit_de_neg$sigma2_alpha[g],0), c(0, fit_de_neg$sigma2_beta[g]))
  rnorm_de <- mvrnorm(1000, mu, acov, empirical=FALSE)
  ab_de <- rnorm_de[,1]*rnorm_de[,2]
  mc_de <- cbind(mc_de, ab_de)
}
names(mc_de) <- gvec_neg[is_de_neg]

saveRDS(mc_de, file = "../data/mediation_DF/mc_de_negative.rds")

mc_node <- list()
for (g in 1:length(isnot_de_gen)) {
  mu <- c(fit_node_neg$a[g], fit_node_neg$b[g])
  acov <- cbind(c(fit_node_neg$sigma2_alpha[g],0), c(0, fit_node_neg$sigma2_beta[g]))
  rnorm_node <- mvrnorm(1000, mu, acov, empirical=FALSE)
  ab_de <- rnorm_node[,1]*rnorm_node[,2]
  mc_node <- cbind(mc_node, ab_de)
}
names(mc_node) <- gvec_neg[isnot_de_gen]


saveRDS(mc_node, file = "../data/mediation_DF/mc_node_negative.rds")
