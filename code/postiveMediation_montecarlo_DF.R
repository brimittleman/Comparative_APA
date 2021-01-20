load("../data/mediation_DF/positive_mediation.rda")
library(MASS)

# columns: genes, rows: random sample
mc_de <- list()
for (g in 1:length(is_de)) {
  mu <- c(fit_de_pos$a[g], fit_de_pos$b[g])
  acov <- cbind(c(fit_de_pos$sigma2_alpha[g],0), c(0, fit_de_pos$sigma2_beta[g]))
  rnorm_de <- mvrnorm(1000, mu, acov, empirical=FALSE)
  ab_de <- rnorm_de[,1]*rnorm_de[,2]
  mc_de <- cbind(mc_de, ab_de)
}
names(mc_de) <- gvec[is_de]

saveRDS(mc_de, file = "../data/mediation_DF/mc_de_postive.rds")

mc_node <- list()
for (g in 1:length(isnot_de)) {
  mu <- c(fit_node_pos$a[g], fit_node_pos$b[g])
  acov <- cbind(c(fit_node_pos$sigma2_alpha[g],0), c(0, fit_node_pos$sigma2_beta[g]))
  rnorm_node <- mvrnorm(1000, mu, acov, empirical=FALSE)
  ab_de <- rnorm_node[,1]*rnorm_node[,2]
  mc_node <- cbind(mc_node, ab_de)
}
names(mc_node) <- gvec[isnot_de]


saveRDS(mc_node, file = "../data/mediation_DF/mc_node_postive.rds")
