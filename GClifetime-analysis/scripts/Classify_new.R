# script to plot total number of GCs


# cutoff of number of B cells to classify gcs # only threshold is used for GClifetime paper

# < threshold - not visible
# threshold - small -> small
# small - medium -> medium
# medium - large -> large
# total - small+medium+large

threshold = 100

large = 4000

small = 800

medium = 2500

# data from Rao et al, JI:

time_rao <- c(4,6,8,10,12,14,16,20,23,26,30,40)
ngcs_rao <- c(4.39,67.36,66.66,81.40,87.53,78.76,67.19,56.31,51.21,40.33,39.45,14.53) # mean
Ebu_rao <- c(NA,71.75,68.95,84.56,92.10,84.03,70.88,59.47,58.25,43.50,49.98,16.48) # upper error bar
SD_rao <- Ebu_rao-ngcs_rao # standard deviation

NGCs_rao <- cbind(time_rao,ngcs_rao,SD_rao)
NGCs_rao <- data.frame(NGCs_rao)

npcgg_data <- read.delim("NP-CGG.txt",header=TRUE,sep="") # Al-Qahtani data

# count number of GCs

new_gcs_count_mean <- NULL
new_gcs_count_sd <- NULL


calc_s_total <- NULL
calc_s_small <- NULL
calc_s_medium <- NULL
calc_s_large <- NULL

for(i in 1:nsims){


# Read the GCvolume files

read_s_gcvol <- read.delim(paste0('GCvolume_',i),header=FALSE,sep=" ")

# Total number of GCs simulated

ngcs_sim <- read.delim("nGCs_1",header=FALSE,sep=" ")
max_gcs <- sum(ngcs_sim)


ntimepoints <- nrow(read_s_gcvol)

calc_s_tmp <- c(rep(0,ntimepoints))

calc_s_tmp_small <- c(rep(0,ntimepoints))
calc_s_tmp_medium <- c(rep(0,ntimepoints))
calc_s_tmp_large <- c(rep(0,ntimepoints))



for(j in 1:ntimepoints){

    for(k in 1:nrow(ngcs_sim)){
    
        if(read_s_gcvol[j,k+1] >=threshold){
        
            # count total number of GCs
            
            calc_s_tmp[j] <- calc_s_tmp[j]+ngcs_sim[k,1]
        
        }
        
        # count small, medium and large GCs
    
         if(read_s_gcvol[j,k+1]>=threshold && read_s_gcvol[j,k+1]<small){
            
            calc_s_tmp_small[j] <- calc_s_tmp_small[j]+ngcs_sim[k,1]
        
        }
        
        if(read_s_gcvol[j,k+1]>=small && read_s_gcvol[j,k+1]<medium){
            
            calc_s_tmp_medium[j] <- calc_s_tmp_medium[j]+ngcs_sim[k,1]
        
        }
        
        if(read_s_gcvol[j,k+1]>=medium){
            
            calc_s_tmp_large[j] <- calc_s_tmp_large[j]+ngcs_sim[k,1]
        
        }
    
    }

}

calc_s_total <- cbind(calc_s_total,calc_s_tmp)

calc_s_small <- cbind(calc_s_small,calc_s_tmp_small)
calc_s_medium <- cbind(calc_s_medium,calc_s_tmp_medium)
calc_s_large <- cbind(calc_s_large,calc_s_tmp_large)

}


# calculate mean and sd 

for(i in 1:ntimepoints){

    if(nsims==1){

        tmp_count <- c(read_s_gcvol[i,1],calc_s_total[i,1],calc_s_small[i,1],calc_s_medium[i,1],calc_s_large[i,1])
        tmp_count_sd <- c(read_s_gcvol[i,1],0,0,0,0)

    }else{
    
        tmp_count <- c(read_s_gcvol[i,1],mean(calc_s_total[i,1:nsims]),mean(calc_s_small[i,1:nsims]),mean(calc_s_medium[i,1:nsims]),mean(calc_s_large[i,1:nsims]))
        tmp_count_sd <- c(read_s_gcvol[i,1],sd(calc_s_total[i,1:nsims]),sd(calc_s_small[i,1:nsims]),sd(calc_s_medium[i,1:nsims]),sd(calc_s_large[i,1:nsims]))    
    
    }

    new_gcs_count_mean <-  rbind(new_gcs_count_mean,tmp_count)
    new_gcs_count_sd <- rbind(new_gcs_count_sd,tmp_count_sd)

}

new_gcs_count_mean <- data.frame(new_gcs_count_mean)
new_gcs_count_sd <- data.frame(new_gcs_count_sd)

write.table(new_gcs_count_mean,file="new_gcs_count_mean", sep = " ", row.names = FALSE, col.names = FALSE)
write.table(new_gcs_count_sd,file="new_gcs_count_sd", sep = " ", row.names = FALSE, col.names = FALSE)


# extract simulation results corresponding to timepoints in experimental data

data_tmp_ngcs <- NULL
time_tmp_ngcs <- NULL
data_ngcs <- NULL

for(loop in 1:nrow(NGCs_rao)){

    index_extract<-which(abs(new_gcs_count_mean[,1]-NGCs_rao[loop,1])==min(abs(new_gcs_count_mean[,1]-NGCs_rao[loop,1])))
    data_tmp_ngcs <- c(data_tmp_ngcs,new_gcs_count_mean[index_extract,2])
    time_tmp_ngcs <- c(time_tmp_ngcs,new_gcs_count_mean[index_extract,1])

}

data_ngcs <- cbind(NGCs_rao[,1:2],time_tmp_ngcs,data_tmp_ngcs)

write.table(data_ngcs,file="../data_ngcs", sep = " ", row.names = FALSE, col.names = FALSE)




# plot with SRBC data

png(file="totalgcs_newcalc.png",width=4, height=4, units="in", res=300)

pngcs_new <- ggplot() + 
   geom_point(data=NGCs_rao, aes(x=NGCs_rao[,1], y=NGCs_rao[,2],color="black"))+
   xlab("Time (days)")+
   ylab("No. of GCs")+
   geom_errorbar(data=NGCs_rao,aes(x=NGCs_rao[,1],ymin=NGCs_rao[,2]-NGCs_rao[,3], ymax=NGCs_rao[,2]+NGCs_rao[,3],color="black"), width=1.0,position=position_dodge(0.05))+theme_bw()+theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=22),axis.text.x  = element_text(size=19),legend.position ="none",legend.title=element_blank())+
   geom_line(data = new_gcs_count_mean, aes(x=new_gcs_count_mean[,1], y=new_gcs_count_mean[,2],color="red"))+
   geom_ribbon(data=new_gcs_count_sd,aes(x=new_gcs_count_sd[,1],ymin=new_gcs_count_mean[,2]-new_gcs_count_sd[,2],ymax=new_gcs_count_mean[,2]+new_gcs_count_sd[,2]),fill = "red",alpha = 0.15,colour=NA,show.legend=FALSE)+
       scale_color_manual(values = c("black","red"),labels=c("Data","Simulation"))+guides(col = guide_legend())+coord_cartesian(ylim = c(0, 100),xlim=c(0,40))
                 
print(pngcs_new)                 
dev.off()

# plot with NP-CGG data 

png(file="totalgcs_npcgg_sim_newcalc.png",width=4, height=4, units="in", res=300)

pngcs <- ggplot() + 
   geom_point(data=npcgg_data, aes(x=npcgg_data[,1], y=npcgg_data[,2],color="black"))+
   xlab("Time (days)")+
   ylab("No. of GCs")+
   geom_errorbar(data=npcgg_data,aes(x=npcgg_data[,1],ymin=npcgg_data[,2]-npcgg_data[,3], ymax=npcgg_data[,2]+npcgg_data[,3],color="black"), width=1.0,position=position_dodge(0.05))+theme_bw()+theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=22),axis.text.x  = element_text(size=19),legend.position = "top",legend.title=element_blank())+
   geom_line(data = new_gcs_count_mean , aes(x=new_gcs_count_mean[,1], y=new_gcs_count_mean[,2],color="red"))+
   geom_ribbon(data=new_gcs_count_sd,aes(x=new_gcs_count_sd[,1],ymin=new_gcs_count_mean[,2]-new_gcs_count_sd[,2],ymax=new_gcs_count_mean[,2]+new_gcs_count_sd[,2]),fill = "red",alpha = 0.3,colour=NA,show.legend=FALSE)+
   scale_color_manual(values = c("black","red"),labels=c("Simulation","Data"))+coord_cartesian(xlim=c(0,40))
                 
print(pngcs)                 
dev.off()

