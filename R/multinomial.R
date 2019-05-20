require(rjags)
setwd('~/spore_counts/')
plot.loc <- 'plots/'

source(paste0('R/multinomial_real_data.R'))


# binom_result <- 
# binom_model(N_trials=1,
#             N = 78,
#             z = 48,
#             nr.tracer = 10000,
#             model.path = model.path,
#             model_binom = 'binom.txt',
#             n_iter1 = 100000,
#             thin=1,
#             nchain=10)
# 
# 
# dens <- density(binom_result$pred_nr_spores)
# dens1 <- density(1-binom_result$probability_tracer)
# 
# plot(dens)
# abline(v=6250)
# 
# 
# plot(dens1)
# abline(v=30/78)
# 
# jags.data.binom = list(N_trials = N_trials,N = N, z = z)
# modjags.binom <-jags.model(file = paste(model.path,model_binom,sep=''),data=jags.data.binom,n.chains = nchain1)
# samples.binom <- coda.samples(modjags.binom,variable.names=c('theta'),n.iter = n_iter1,thin=thin1)
# plot(samples.binom)
# 
# prob.tracer <- unlist(samples.binom) 
# pred.spore <- round(nr.tracer*(1-prob.tracer)/prob.tracer)


counts <- as.data.frame(matrix(ncol=3,
                               data=c(30,48,0,30,48,10,30,48,50,30,48,100),
                               byrow=TRUE))

multinom_prob <- lapply(1:4,function(bb){ 
multinom(N_sample=1, # number of samples from multinomial distribution to consider
                     y=counts[bb,], #vector of counted spores
                     K=length(counts[bb,]), #number of taxa counted
                     N_multinom=sum(counts[bb,]), # total spores counted
                     model.path=model.path, # where is the JAGS model code saved
                     model_multinom='multinom.txt', # name of the multinomial model code
                     n_iter2=1000000,
                     thin2=1,
                     nchain2=10)
})

n_tracer <- 10000



# n_spores1 <- n_tracer/multinom_prob[[1]][,2]*multinom_prob[[1]][,1] 
# n_spores2 <- n_tracer/multinom_prob[[2]][,2]*multinom_prob[[2]][,1]
# n_spores3 <- n_tracer/multinom_prob[[3]][,2]*multinom_prob[[3]][,1]
# n_spores4 <- n_tracer/multinom_prob[[4]][,2]*multinom_prob[[4]][,1]
# 
# dens1 <- density(n_spores1)
# dens2 <- density(n_spores2)
# dens3 <- density(n_spores3)
# dens4 <- density(n_spores4)
# 
# plot(dens1)
# lines(dens2,col=2)
# lines(dens3,col=3)
# lines(dens4,col=3)
# abline(v=6250)


####################################################################################################################
#
####################################################################################################################

d1 <- density(multinom_prob[[1]][,1]/(multinom_prob[[1]][,1]+multinom_prob[[1]][,2]))
d2 <- density(multinom_prob[[2]][,1]/(multinom_prob[[2]][,1]+multinom_prob[[2]][,2]))
d3 <- density(multinom_prob[[3]][,1]/(multinom_prob[[3]][,1]+multinom_prob[[3]][,2]))
d4 <- density(multinom_prob[[4]][,1]/(multinom_prob[[4]][,1]+multinom_prob[[4]][,2]))



pdf(paste0(plot.loc,'bivariate_probability.pdf'),height=5,width=6)
  plot(d1,xlab='',main='')
  lines(d2,col=2,lty=2)
  lines(d3,col=3,lty=3)
  lines(d3,col=4,lty=4)
  abline(v=30/78,lty=2)
  mtext(side=1,line=2.2,'p')
  mtext(side=3,line=0,font=2,'Posterior Distribution')
  legend('topright',lty=1:4,col=1:4,paste(c(0,10,50,100),'pollen grains'))
dev.off()

test <- multinom_prob[[1]][,1]/(multinom_prob[[1]][,1]+multinom_prob[[1]][,2])
mean(test)
