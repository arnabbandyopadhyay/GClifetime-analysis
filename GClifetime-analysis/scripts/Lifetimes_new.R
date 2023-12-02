# script to calculate lifetime of GCs, uses threshold set in Classify_new.R

print(threshold)

death_time <- NULL

lifetime <- NULL
lifetime_mean <- NULL
lifetime_sd <- NULL

for(i in 1:nsims){

    read_s_gcvol <- read.delim(paste0('GCvolume_',i),header=FALSE,sep="")
    
    # for each representative GC, calculate the lifetime
    
    n_rep_gcs <- ncol(read_s_gcvol)-1
    
    gcnumbers <- c(1:n_rep_gcs)
    how_long <- c(rep(0,n_rep_gcs))
    
    death_time_tmp <- c(rep(200,n_rep_gcs))
    
    for(i in 1:n_rep_gcs){

        birth=0
        death=0
        
        done=0
        
        for(j in 2:(nrow(read_s_gcvol)-1)){
        
        if(done==0 & read_s_gcvol[j-1,i+1]<=threshold & read_s_gcvol[j,i+1]>=threshold){
            
            birth=read_s_gcvol[j,1]
            done=1
        }
        
        if(read_s_gcvol[j,i+1]>=threshold & read_s_gcvol[j+1,i+1]<=threshold){
            
            death=read_s_gcvol[j,1]
            
        
        }
        
        
        }
        print(birth)
        print(death)
        
        how_long[i]=round((death-birth),2)
        
        print(how_long[i])
        
        if(how_long[i]<0){
          how_long[i]=NA
        }
        
        if(death>0){
           death_time_tmp[i]=death
        }
        
        
    }
    
    # shutdown timing
  death_time <- cbind(death_time,death_time_tmp)

  # lifetime 
  lifetime <- cbind(lifetime,how_long)

}

if(nsims==1){

    lifetime_mean <- cbind(gcnumbers,lifetime)
    lifetime_sd <- cbind(gcnumbers,c(rep(0,length(gcnumbers))))

}else{

    lifetime_mean <- cbind(gcnumbers,rowMeans(lifetime))
    lifetime_sd <- cbind(gcnumbers,apply(lifetime,1,sd))

}


lifetime_mean <- data.frame(lifetime_mean)
lifetime_sd <- data.frame(lifetime_sd)

write.table(lifetime_mean, file = "lifetime_mean.out", sep = " ", row.names = FALSE, col.names = FALSE)
write.table(lifetime_sd, file = "lifetime_sd.out", sep = " ", row.names = FALSE, col.names = FALSE)

write.table(death_time, file = "death_time.out", sep = " ", row.names = FALSE, col.names = FALSE)


png(file="lifetimes_newcalc.png",width=4, height=4, units="in", res=300)

pl_new <- ggplot(data=lifetime_mean, aes(x=as.factor(lifetime_mean[,1]), y=lifetime_mean[,2],fill=as.factor(lifetime_mean[,1]),color=as.factor(lifetime_mean[,1]))) + 
  geom_bar(stat = "identity",width=0.5)+ 
   xlab("Simulated GC")+
   ylab("Lifetime (days)")+
   geom_errorbar(data=lifetime_sd,aes(x=as.factor(lifetime_sd[,1]),ymin=lifetime_mean[,2]-lifetime_sd[,2], ymax=lifetime_mean[,2]+lifetime_sd[,2]), width=0.1,position=position_dodge(0.05))+
  theme_bw()+theme(axis.title.y = element_text(size=19),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=19),axis.text.x  = element_text(angle=90,size=10),legend.position = "none",legend.title=element_blank())

                 
print(pl_new)                 
dev.off()
