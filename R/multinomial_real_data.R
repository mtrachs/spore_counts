require(rjags)


wd.path <- '~/spore_counts/'
setwd(wd.path)

data.loc <- 'data'

model.path <-  'jags_code'

data <- read.csv(paste0(data.loc,'/spore_counts.csv'),header=TRUE)



##################################################################################################################
# simple model to make inference on theta (probability of success) for data that follows a binomial distribution
# model can be used to make predictions from posterior
# model has a uniform prior on theta
##################################################################################################################

sink(paste(model.path,"/binom.txt",sep=''))
cat('
    model {
	        for (i in 1:N_trials) {
            z[i] ~ dbin(theta , N[i])
	        }
	        theta  ~ dbeta( 1 , 1 )
          }
    ',fill = TRUE)
sink()


sink(paste(model.path,"/multinom.txt",sep=''))
cat('
    model {
	          for (i in 1:N_sample) {
                y[i,] ~ dmulti(q[], N_multinom[i])  
            }
	          q[1:K] ~ ddirch(alpha[])
            for (k in 1:K)  {
              alpha[k] <- 1
            }
          }
    ',fill = TRUE)
sink()





##################################################################################################################
#start bringing data into form that JAGS will need
##################################################################################################################
nr_tracer <- data$Tracer.Added
data.use <- data[,-c(1:3,ncol(data))] 
data.use <- data.use[,colSums(data.use)>0]



#################################################################################################################
# binomial model 
#################################################################################################################
binom_model <- function(N_trials=1, #number of samples from binomial distribution to consider currently set to 1
                N, #total count for binomial distribution
                z, # number of tracers counted
                nr.tracer, # total number of tracers added to the sample
                model.path, # where is the JAGS model code saved
                model_binom, # name of the binomial model code
                n_iter1, # number of iterations for binomial model
                thin1, #thinning for binomial model
                nchain1){ #number of chains for binomial model
                

  jags.data.binom = list(N_trials = N_trials,N = N, z = z)
  modjags.binom <-jags.model(file = paste(model.path,model_binom,sep=''),data=jags.data.binom,n.chains = nchain1)
  samples.binom <- coda.samples(modjags.binom,variable.names=c('theta'),n.iter = n_iter1,thin=thin1)
  plot(samples.binom)
  
  prob.tracer <- unlist(samples.binom) 
  pred.spore <- round(nr.tracer*(1-prob.tracer)/prob.tracer)
  par(mfrow=c(1,1))
  hist(pred.spore,nclass=15)
  
  output <- data.frame(probability_tracer =prob.tracer,
                 pred_nr_spores = pred.spore)
}
##################################################################################################################
# combination of binomial and multinomial model
##################################################################################################################

bin_multinom <- function(N_trials=1, #number of samples from binomial distribution to consider currently set to 1
                         N, #total count for binomial distribution
                         z, # number of tracers counted
                         nr.tracer, # total number of tracers added to the sample
                         N_sample=1, # number of samples from multinomial distribution to consider
                         y, #vector of counted spores
                         K, #number of taxa counted
                         N_multinom, # total spores counted
                         model.path, # where is the JAGS model code saved
                         model_binom, # name of the binomial model code
                         model_multinom, # name of the multinomial model code
                         n_iter1, # number of iterations for binomial model
                         n_iter2,
                         thin1, #thinning for binomial model
                         thin2,
                         nchain1, #number of chains for binomial model
                         nchain2){




  
    
#################################################################################################################
# binomial model
# both models have the same structure: 
#  1)prepare data for model, variable names have to match variable names in model specification
#  2) pass data to JAGS
#  3) ask JAGS to sample from the posterior 
#  4) post-process output   
#################################################################################################################

  
  
  jags.data.binom = list(N_trials = N_trials,N = N, z = z)
  modjags.binom <-jags.model(file = paste(model.path,model_binom,sep=''),data=jags.data.binom,n.chains = nchain1)
  samples.binom <- coda.samples(modjags.binom,variable.names=c('theta'),n.iter = n_iter1,thin=thin1)
  plot(samples.binom)

  prob.tracer <- unlist(samples.binom) 
  pred.spore <- round(nr.tracer*(1-prob.tracer)/prob.tracer)
  par(mfrow=c(1,1))
  hist(pred.spore,nclass=15)








#################################################################################################################
# multinomial model estimate probabilities for spores
#################################################################################################################
  jags.data.multinom = list(N_sample = N_sample,y = y, K = K,N_multinom = N_multinom)

    modjags.multinom <-jags.model(file = paste(model.path, model_multinom,sep=''),data=jags.data.multinom,n.chains = nchain2)
# sample from the posterior
  samples.multinom <- coda.samples(modjags.multinom,variable.names=c('q'),n.iter = n_iter2,thin=thin2)
  plot(samples.multinom)
  
  if(length(samples.multinom)>1){
    prob_multinom <- rep(NA,ncol(samples.multinom[[1]])) 
    
    for(i in 1:length(samples.multinom)){
      prob_multinom <- rbind(prob_multinom,samples.multinom[[i]])
    }
    
     prob_multinom <- prob_multinom[-1,] 
  }else{
    prob_multinom <- samples.multinom[[1]]
  }

  #make output
  output <- list(probability_tracer =prob.tracer,
                 pred_nr_spores = pred.spore,
                 prob_mult_spores =  prob_multinom)
}




test <- bin_multinom(N_trials = 1,
         N = sum(data.use[3,]),
         z = data.use[3,ncol(data.use)],
         nr.tracer = nr_tracer[3],
         N_sample=1,
         y=data.use[3,-ncol(data.use)],
         K=ncol(data.use)-1,
         N_multinom= sum(data.use[3,-ncol(data.use)]),
         model.path=model.path,
         model_binom ="/binom.txt",
         model_multinom = "/multinom.txt",
         n_iter1=5000,
         n_iter2=10000,
         thin1=5,
         thin2=10,
         nchain1=3,
         nchain2=3)













