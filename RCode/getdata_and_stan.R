message("Working directory for output set to: ", dir_work)
if (!file.exists(file.path(dir_data, my_csv_file))){
  message(paste0("Cannot find your data file: ", my_csv_file))    
} else {
  data_in <- read.csv(file.path(dir_data, my_csv_file), as.is=TRUE)
  message("Read data into R file: data_in")
  str(data_in)
}
if (!file.exists(file.path(dir_stan, my_stan_code))){
  message(paste0("Cannot find Stan Code in: ", dir_stan))    
} else {
  message("Compiling Stan code...")
  HB_model <- cmdstan_model(file.path(dir_stan,my_stan_code), quiet = TRUE, cpp_options = list(stan_threads = TRUE))
  message("DONE")
}
