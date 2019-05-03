require(rjags)

wd.path <- '~/spore_counts/'
setwd(wd.path)

data.loc <- 'data'
code.loc <- 'R/'
output.loc <- 'output/'
source(paste0(code.loc,'multinomial_real_data.R'))

model.path <-  'jags_code'

data <- read.csv(paste0(data.loc,'/simulated_data_v1.csv'),header=TRUE)

nr_tracer <- data$X..tracers.added..true.
Total_count <- data$X..counts
spore_counted <- data$X..spores.counted
tracer_counted <- data$X..tracers.counted

start <- Sys.time()
out1<- 
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
end <- Sys.time()
end-start

saveRDS(out1,paste0(output.loc,'estimations_simulated_data_v1.RDS'))


#########################################################################################################################
# negative binomial model
#########################################################################################################################
data.neg.bin <- read.csv(paste0(data.loc,'/simulated_data_count_tracers_v1.csv'),header=TRUE)

start <- Sys.time()
neg_bin_out <- 
  lapply(1:nrow(data.neg.bin),function(x){
neg_binom_model(N_trials = 1,
            r=data.neg.bin$X..spores.counted[x],
            z=data.neg.bin$X..tracers.counted[x],
            nr.tracer = data.neg.bin$X..tracers.added..true.[x],
            model.path=model.path,
            model_neg_binom ="/neg_binom.txt",
            n_iter1 = 5000,
            nchain1 = 3,
            thin1 = 10)

  })
end <- Sys.time()
end-start
 

saveRDS(neg_bin_out,paste0(output.loc,'simulated_data_count_tracers_v1_neg_binom.RDS'))




start <- Sys.time()
neg_bin_out_use_bin <- 
  lapply(1:nrow(data.neg.bin),function(x){
    binom_model(N_trials = 1,
                    N=data.neg.bin$X..counts[x],
                    z=data.neg.bin$X..tracers.counted[x],
                    nr.tracer = data.neg.bin$X..tracers.added..true.[x],
                    model.path=model.path,
                    model_binom ="/binom.txt",
                    n_iter1 = 5000,
                    nchain1 = 3,
                    thin1 = 10)
    
  })
end <- Sys.time()
end-start


saveRDS(neg_bin_out_use_bin,paste0(output.loc,'simulated_data_count_tracers_v1_neg_binom_use_binom.RDS'))
##