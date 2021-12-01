# Parameters:
# nchains, outname, dir_output, dir_work, out_prefix, data_stan

### Read results from csv into R ###
cat("Reading draws from Stan csv output into R (large files take time)...")
nchains <- threads[[1]]
csv_name <- do.call(c, lapply(1:nchains, function(i) paste0(outname,"-",i,".csv")))
draws_beta <- read_cmdstan_csv(file.path(dir_work, csv_name), variables = "beta_ind", format = "draws_list", sampler_diagnostics = "accept_stat__")
cat("DONE")
                              
### Save output files and check convergence ###
draws_name <- paste0(out_prefix,"_draws_beta.rds")
util_name <- paste0(out_prefix,"_utilities_r.csv")
pdf_name <- paste0(out_prefix,"_trace_plots.pdf")
fit_name <-  paste0(out_prefix,"_fit_stats.csv")
message(paste0(
  "\nSaving post warm-up files for:\n",
  " respondent point estimates:    ", util_name,"\n",  
  " draws of utilities as R list:  ", draws_name,"\n",
  " convergence stats of mean:     ", fit_name, "\n",
  " PDF of detailed traceplots:    ", pdf_name,"\n",
  "\nShowing post warm-up:\n",
    " Acceptance rate across iterations (histogram)\n",
    " Traceplots of all mean utilities together (Sawtooth chart)"
))

hist(do.call(rbind,draws_beta$post_warmup_sampler_diagnostics)$accept_stat__, breaks = 30, main = "Acceptance Rate - Sampling", xlab = "", xlim = c(0,1))
saveRDS(modifyList(draws_beta,list(warmup_draws = NULL)), file.path(dir_work, draws_name)) # drop warmup
utilities <- matrix(
            Reduce("+",lapply(draws_beta$post_warmup_draws, colMeans))/nchains,
            data_stan$I, data_stan$P,
            byrow = TRUE) # First P entries are respondent 1, next P are resp 2
utilities_r <- utilities %*% t(data_stan$code_master)
write.table(cbind(id = data_stan$resp_id, utilities_r), file = file.path(dir_work, util_name), sep = ",", na = ".", row.names = FALSE)
                              
# Convergence charts saved as pdf and in fit_stats
fit_stats <- data.frame(
  variable = colnames(data_stan$ind),
  mean = NA,
  rhat = NA,
  ESS = NA
)
ndraws <- nrow(draws_beta$post_warmup_draws[[1]])
draws_beta_mu <- list() # Creates the mean of respondent utilities for each iteration, like alpha
for (chain_i in (1:nchains)){
  draws_beta_list <- as.matrix(draws_beta$post_warmup_draws[[chain_i]])
  draws_beta_mu[[chain_i]] <- t(sapply(1:ndraws, function(draw){
    beta_mu <- colMeans(matrix(draws_beta_list[draw,],
                               data_stan$I,data_stan$P, byrow = TRUE))
  }))
  matplot(1:nrow(draws_beta_mu[[chain_i]]), draws_beta_mu[[chain_i]],
          type = "l" , lty = 1, lwd = 1, main = paste0("Chain ", chain_i), xlab = "Iteration", ylab = "Mean Beta")   
} 

pdf(file = file.path(dir_work, pdf_name),   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches
hist(do.call(rbind,draws_beta$post_warmup_sampler_diagnostics)$accept_stat__, breaks = 30, main = "Acceptance Rate - Sampling", xlab = "", xlim = c(0,1))
for (chain_i in (1:nchains)){
  matplot(1:nrow(draws_beta_mu[[chain_i]]), draws_beta_mu[[chain_i]],
          type = "l" , lty = 1, lwd = 1, main = paste0("Chain ", chain_i), xlab = "Iteration", ylab = "Mean Beta")   
}  
chain_cols <- c("red","blue","green","black")
for (i in 1:ncol(draws_beta_mu[[1]])){
  x <- sapply(1:length(draws_beta_mu), function(chain){
    draws_beta_mu[[chain]][,i]     
  }) # x is set of column i across draws_beta_mu
  fit_stats$mean[i] <- round(mean(x), 2)
  fit_stats$rhat[i] <- round(rhat(x),2)
  fit_stats$ESS[i] <- round(ess_basic(x),1)
  plot(x[,1], type = "l", col = chain_cols[1], ylim = c(min(x), max(x)),
                          xlab = "Sample Iteration", ylab = "Mean Beta",
                          main = paste(colnames(data_stan$ind)[i],
                                                    "| rhat = ", round(rhat(x),2),
                                                    "| ESS = ", round(ess_basic(x),1)
  ))
  for (chain in 2:nchains){
    lines(x[,2], type = "l", col = chain_cols[chain])
  }
}
dev.off()
write.table(fit_stats, file = file.path(dir_work, paste0(out_prefix,"_fit_stats.csv")), sep = ",", na = ".", row.names = FALSE)
