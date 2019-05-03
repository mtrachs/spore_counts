require(rjags)


wd.path <- '~/spore_counts/'
setwd(wd.path)

data.loc <- 'data'
code.loc <- 'R/'
source(paste0(code.loc,'multinomial_real_data.R'))

model.path <-  'jags_code'

data <- read.csv(paste0(data.loc,'/simulated_data_v1.csv'),header=TRUE)


nr_tracer <- data$X..tracers.added..true.
Total_count <- data$X..counts
spore_counted <- data$X..spores.counted
tracer_counted <- data$X..tracers.counted


out1 <- 
  lapply(1:nrow(data),function(x){  
    binom_model(N_trials = 1,
            N=Total_count[x],
            z=tracer_counted[x],
            nr.tracer = nr_tracer[x],
            model.path=model.path,
            model_binom ="/binom.txt",
            n_iter1 = 5000,
            nchain1 = 3,
            thin1 = 10)
  })  
