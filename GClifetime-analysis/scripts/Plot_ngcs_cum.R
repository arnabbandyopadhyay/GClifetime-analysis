# script to plot cumulative number of initialized GCs

library(ggplot2)


times <- seq(0,1000,by=2)

cum_gcs_count_mean <- NULL
cum_gcs_count_sd <- NULL

cum_gcs <- NULL

ntimepoints <- length(times)

for(i in 1:nsims){

    # read file with initiation time of GCs
    read <- read.delim(paste0('ini_time_',i),header=FALSE,sep=" ")
    read_ngcs <- read.delim(paste0('nGCs_',i),header=FALSE,sep=" ")

    ngcs=nrow(read)
    
    # calculated number of GCs initialized at or before time point
    
    tmp <- rep(0,ntimepoints)
    
    for(j in 1:ntimepoints){
    
        for(k in 1:ngcs){
        
            if(read[k,1]<=times[j]){
            
               tmp[j]=tmp[j]+read_ngcs[k,1];

            
            }
        
        
        }
        
    
    }

    cum_gcs <- cbind(cum_gcs,tmp)
   
}

times <- times/24

# calculate mean and SD

if(nsims==1){

cum_gcs_count_mean=cbind(times,cum_gcs)

cum_gcs_count_sd=cbind(times,rep(0,nrow(cum_gcs)))

}else{

cum_gcs_count_mean=cbind(times,rowMeans(cum_gcs))

cum_gcs_count_sd=cbind(times,apply(cum_gcs,1,sd))

}

cum_gcs_count_mean <- data.frame(cum_gcs_count_mean)
cum_gcs_count_sd <- data.frame(cum_gcs_count_sd)

write.table(cum_gcs_count_mean,file="cum_gcs_count_mean", sep = " ", row.names = FALSE, col.names = FALSE)
write.table(cum_gcs_count_sd,file="cum_gcs_count_sd", sep = " ", row.names = FALSE, col.names = FALSE)


png(file="totalgcs_cumulative.png",width=4, height=4, units="in", res=300)

pngcs_new <- ggplot() + 
   xlab("Time (days)")+
   ylab("new GCs")+
   theme_bw()+theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=22),axis.text.x  = element_text(size=19),legend.position ="none",legend.title=element_blank())+
   geom_line(data = cum_gcs_count_mean, aes(x=cum_gcs_count_mean[,1], y=cum_gcs_count_mean[,2]),color="magenta")+
   geom_ribbon(data=cum_gcs_count_sd,aes(x=cum_gcs_count_sd[,1],ymin=cum_gcs_count_mean[,2]-cum_gcs_count_sd[,2],ymax=cum_gcs_count_mean[,2]+cum_gcs_count_sd[,2]),fill = "magenta",alpha = 0.15,colour=NA,show.legend=FALSE)
                 
print(pngcs_new)                 
dev.off()



