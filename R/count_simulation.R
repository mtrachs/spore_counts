####
setwd('~/spore_counts/')
plot.loc <- 'plots/'
n_tracer <- 10000
c_spore <- 30
c_tracer <- 48
c_pollen <- 10



count_simulation <- 
function(n_tracer,
         c_spore,
         c_tracer,
         c_pollen,
         cutoff_criterion,
         n_replicates,
         tot_count){

    n_spore  <- round(c_spore/c_tracer*n_tracer)
    n_pollen <- round(c_pollen/(c_spore+c_tracer)*(n_spore+n_tracer))
    total_counts <- c(rep(1,n_tracer),rep(2,n_spore),rep(3,n_pollen))  
  
    
    sim_count <- 
      replicate(n_replicates,{
        sim_counts <- sample(total_counts,length(total_counts),replace=FALSE)
        if(cutoff_criterion=='spore') cutoff_point <- which(sim_counts==2)[c_spore]
        if(cutoff_criterion=='tracer') cutoff_point <- which(sim_counts==1)[c_tracer]
        if(cutoff_criterion=='pollen') cutoff_point <- which(sim_counts==3)[c_pollen]
        if(cutoff_criterion=='total_counts') cutoff_point <- tot_count
        
        sim <- c(sim_counts[1:cutoff_point],1:3)
        table(sim)-1
      })
    
}


observed_counts <- matrix(ncol=3, data=c(rep(48,4),rep(30,4),c(0,10,50,100)))

start <- Sys.time()
simulated_counts <- lapply(1:4,function(x){     
  count_simulation(n_tracer=10000,
                   c_spore=observed_counts[x,2],
                   c_tracer=observed_counts[x,1],
                   c_pollen=observed_counts[x,3],
                   cutoff_criterion = 'tracer',
                   n_replicates = 10000000)

})  
end <- Sys.time()  

end - start


n_tracer <- 10000

pdf(paste0(plot.loc,'simulated_counts_cutoff_tracer'),height=5,width=6)
sapply(1:4,function(x){
  cc <- simulated_counts[[x]]
  dd <- table(cc[2,]/cc[1,]*n_tracer)
  if(x==1){
    plot(as.numeric(names(dd)),dd/sum(dd),xlim=c(0,20000),xlab='',ylab='')
    axis(2)
    abline(v=6250,lty=2,col='brown')
    legend('topright',legend=paste(c(0,10,50,100),'pollen grains'),pch = 1:4, col=1:4)
    mtext(side=3,line=0,font=2,'Cutoff: 48 tracers')
    mtext(side=1,line=2.2,font=2,'# tracers')
    mtext(side=2,line=2.2,font=2,'Probability')
  }else{
    points(as.numeric(names(dd)),dd/sum(dd),pch=x,col=x)
  }
})
dev.off()




#########################################################################################################################
#
######################################################################################################################
observed_counts <- matrix(ncol=3, data=c(rep(48,3),rep(30,3),c(10,50,100)))

start <- Sys.time()
simulated_counts_pollen_cutoff <- lapply(1:3,function(x){     
  count_simulation(n_tracer=10000,
                   c_spore=observed_counts[x,2],
                   c_tracer=observed_counts[x,1],
                   c_pollen=observed_counts[x,3],
                   cutoff_criterion = 'pollen',
                   n_replicates = 5000000)
  
})  
end <- Sys.time()  

end - start

pdf(paste0(plot.loc,'simulated_counts_cutoff_pollen'),height=5,width=6)
sapply(1:3,function(x){
  cc <- simulated_counts_pollen_cutoff[[x]]
  dd <- round(cc[2,]/cc[1,],2)*n_tracer
  #dd <- dd[is.finite(dd)]
  if(x==1){
    ee <- hist(dd,breaks =seq(0,100000,500),plot=FALSE)
    plot(ee$breaks[-1]+250,ee$counts/sum(ee$counts),xlim=c(0,20000),xlab='',ylab='')
    abline(v=6250,lty=2,col=4)
    mtext(side=1,line=2.2,font=2,'# tracers')
    mtext(side=2,line=2.2,font=2,'Probability')
    mtext(side=3,line=0,font=2,'Cutoff: Pollen')
    legend('topright',legend=paste(c(10,50,100),'pollen grains'),pch = 1:3, col=1:3)
  }else{
    ee <- hist(dd,breaks =seq(0,100000,500),plot=FALSE)
    points(ee$breaks[-1]+250,ee$counts/sum(ee$counts),xlim=c(0,20000),col=x,pch=x)
  }
})
dev.off()






########################################################################################################################
# simulate spore cutoff
#########################################################################################################################
observed_counts <- matrix(ncol=3, data=c(rep(48,4),rep(30,4),c(0,10,50,100)))

start <- Sys.time()
simulated_counts_spore_cutoff <- lapply(1:4,function(x){     
  count_simulation(n_tracer=10000,
                   c_spore=observed_counts[x,2],
                   c_tracer=observed_counts[x,1],
                   c_pollen=observed_counts[x,3],
                   cutoff_criterion = 'spore',
                   n_replicates = 1000000)
  
})  
end <- Sys.time()  

end - start

pdf(paste0(plot.loc,'simulated_counts_cutoff_spores'),height=5,width=6)
sapply(1:4,function(x){
  cc <- simulated_counts_spore_cutoff[[x]]
  dd <- table(cc[2,]/cc[1,]*n_tracer)
  if(x==1){
    plot(as.numeric(names(dd)),dd/sum(dd),xlab='',ylab='',xlim=c(0,20000))
    mtext(side=1,line=2.2,font=2,'# tracers')
    mtext(side=2,line=0,font=2,'Probability')
    mtext(side=3,line=2.2,font=2,'Cutoff: 30 spores')
    abline(v=6250,col='orange',lty=2)
    legend('topright',legend=paste(c(0,10,50,100),'pollen grains'),pch = 1:4, col=1:4)
  }else{
    points(as.numeric(names(dd)),dd/sum(dd),pch=x,col=x)
  }
})
dev.off()



########################################################################################################################
# simulate total count cutoff
#########################################################################################################################
observed_counts <- matrix(ncol=3, data=c(rep(48,4),rep(30,4),c(0,10,50,100)))

start <- Sys.time()
simulated_counts_total_cutoff <- lapply(1:4,function(x){     
  count_simulation(n_tracer=10000,
                   c_spore=observed_counts[x,2],
                   c_tracer=observed_counts[x,1],
                   c_pollen=observed_counts[x,3],
                   cutoff_criterion = 'total_counts',
                   tot_count = sum(observed_counts[x,]),
                   n_replicates = 10000000)
  
})  
end <- Sys.time()  

end - start

pdf(paste0(plot.loc,'simulated_counts_cutoff_total counts'),height=5,width=6)
sapply(1,function(x){
  cc <- simulated_counts_total_cutoff[[x]]
  dd <- round(cc[2,]/cc[1,],2)*n_tracer
  #dd <- dd[is.finite(dd)]
  if(x==1){
    ee <- hist(dd,breaks =seq(0,100000,500),plot=FALSE)
    plot(ee$breaks[-1]+250,ee$counts/sum(ee$counts),xlim=c(0,20000),xlab='',ylab='')
    abline(v=6250,lty=2,col=4)
    mtext(side=1,line=2.2,font=2,'# tracers')
    mtext(side=2,line=2.2,font=2,'Probability')
    mtext(side=3,line=0,font=2,'Cutoff: Pollen')
    legend('topright',legend=paste(c(10,50,100),'pollen grains'),pch = 1:3, col=1:3)
  }else{
    ee <- hist(dd,breaks =seq(0,100000,500),plot=FALSE)
    points(ee$breaks[-1]+250,ee$counts/sum(ee$counts),xlim=c(0,20000),col=x,pch=x)
  }
})
dev.off()
