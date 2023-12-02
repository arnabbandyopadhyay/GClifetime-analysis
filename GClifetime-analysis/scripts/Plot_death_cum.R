# script to plot cumulative number of terminated GCs

library(ggplot2)

times <- seq(0,1000,by=2)

times <- times/24

cum_gcs_death_mean <- NULL
cum_gcs_death_sd <- NULL

cum_gcs_death <- NULL

ntimepoints <- length(times)

# read file with shutdown timing of GCs
read <- read.delim("death_time.out",header=FALSE,sep="")

read_ngcs <- read.delim(paste0('nGCs_',1),header=FALSE,sep=" ")

# calculate number of GCs terminated before every time point

nsims=ncol(read)

for(i in 1:nsims){

    ngcs=nrow(read)
    
    tmp <- rep(0,ntimepoints)
    
    for(j in 1:ntimepoints){
    
        for(k in 1:ngcs){
        
            if(read[k,i]<=times[j]){
            
               tmp[j]=tmp[j]+read_ngcs[k,1];

            
            }
        
        
        }
        
    
    }

    cum_gcs_death <- cbind(cum_gcs_death,tmp)
   
}

# calculate mean and SD

if(nsims==1){

cum_gcs_death_mean=cbind(times,cum_gcs_death)

cum_gcs_death_sd=cbind(times,rep(0,nrow(cum_gcs_death)))

}else{

cum_gcs_death_mean=cbind(times,rowMeans(cum_gcs_death))

cum_gcs_death_sd=cbind(times,apply(cum_gcs_death,1,sd))

}

cum_gcs_death_mean <- data.frame(cum_gcs_death_mean)
cum_gcs_death_sd <- data.frame(cum_gcs_death_sd)

write.table(cum_gcs_death_mean,file="cum_gcs_death_mean", sep = " ", row.names = FALSE, col.names = FALSE)
write.table(cum_gcs_death_sd,file="cum_gcs_death_sd", sep = " ", row.names = FALSE, col.names = FALSE)


png(file="totalgcs_cumulative_death.png",width=4, height=4, units="in", res=300)

pngcs_new <- ggplot() + 
   xlab("Time (days)")+
   ylab("terminated GCs")+
   theme_bw()+theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=22),axis.text.x  = element_text(size=19),legend.position ="none",legend.title=element_blank())+
   geom_line(data = cum_gcs_death_mean, aes(x=cum_gcs_death_mean[,1], y=cum_gcs_death_mean[,2]),color="magenta")+
   geom_ribbon(data=cum_gcs_death_sd,aes(x=cum_gcs_death_sd[,1],ymin=cum_gcs_death_mean[,2]-cum_gcs_death_sd[,2],ymax=cum_gcs_death_mean[,2]+cum_gcs_death_sd[,2]),fill = "magenta",alpha = 0.15,colour=NA,show.legend=FALSE)
                 
print(pngcs_new)                 
dev.off()



