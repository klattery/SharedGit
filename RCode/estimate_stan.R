outname <- paste0(out_prefix, "_StanOut_", 
                  format(Sys.time(), '%Y%m%d-%H%M%S')) # edit as desired
message("ESTIMATING...")
message(paste0("Optional lines to run in terminal to check progress:\n",
    "cd ", dir_work, "   # Change to your working directory and then:\n",
    "  awk 'END { print NR - 45 } ' '",outname,"-1.csv'", "                # Count lines in output\n",
    "  tail -n +45 '",outname,"-1.csv'  | cut -d, -f 1-300 > temp.csv", "  # Create temp.csv with first 300 columns\n"))

#####  Run Stan Model ###############
HB_model$sample(modifyList(data_stan, data_model),
                iter_warmup = data_model$iter_warmup,
                iter_sampling = data_model$iter_sampling,
                output_dir = dir_work,
                output_basename = outname,
                chains = threads[[1]],
                parallel_chains = threads[[2]],
                threads_per_chain = threads[[3]],
                save_warmup = TRUE,
                refresh = data_model$refresh,
                seed = 271,
                init = .1,
                adapt_delta = data_stan$adapt_delta,
                show_messages = FALSE,
                validate_csv = FALSE
)
cat("### ESTIMATION DONE ###\n")
