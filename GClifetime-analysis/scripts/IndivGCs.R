# this script plots individual GC volume and plasma cell output

gcvol_int <- read.delim("GCvolume_mean",header=FALSE,sep=" ")
pcn_int <- read.delim("Plasmacells_mean",header=FALSE,sep=" ")

gcvol_int_sd <- read.delim("GCvolume_sd",header=FALSE,sep=" ")
pcn_int_sd <- read.delim("Plasmacells_sd",header=FALSE,sep=" ")


final_int <- NULL # gcvol
final_int_sd <- NULL 

finalpcn_int <- NULL # no. of plasmacells
finalpcn_int_sd <- NULL 

gcs_int_tmp <- NULL 
gcs_int_tmp_sd <- NULL 
pcn_int_tmp <- NULL 
pcn_int_tmp_sd <- NULL

# reformatting the results

for(i in 1:ngcs){
  for(j in 1:nrow(gcvol_int)){
        
        gcs_int_tmp <- c(gcvol_int[j,1],gcvol_int[j,i+1],i)
        pcn_int_tmp <- c(pcn_int[j,1],pcn_int[j,i+1],i)

        gcs_int_tmp_sd <- c(gcvol_int_sd[j,1],gcvol_int_sd[j,i+1],i)
        pcn_int_tmp_sd <- c(pcn_int_sd[j,1],pcn_int_sd[j,i+1],i)
        
        final_int <- rbind(final_int,gcs_int_tmp)
        finalpcn_int <- rbind(finalpcn_int,pcn_int_tmp)

        final_int_sd <- rbind(final_int_sd,gcs_int_tmp_sd)
        finalpcn_int_sd <- rbind(finalpcn_int_sd,pcn_int_tmp_sd)
        
  }

}


final_int<- data.frame(final_int)
finalpcn_int <- data.frame(finalpcn_int)

final_int_sd <- data.frame(final_int_sd)
finalpcn_int_sd <- data.frame(finalpcn_int_sd)

# plot the gcvolume of simulated GCs

png(file="indiv_gcvolume.png",width=4,height=4,units="in",res=200)

c1=ggplot() +
    geom_line(data = final_int, mapping = aes(x = final_int[,1], y = final_int[,2], color = as.factor(final_int[,3])),show.legend=FALSE)+
    geom_ribbon(data=final_int_sd,aes(x=final_int_sd[,1],ymin=final_int[,2]-final_int_sd[,2],ymax=final_int[,2]+final_int_sd[,2],fill = as.factor(final_int[,3])),alpha = 0.3,colour=NA,show.legend=FALSE)+
       xlab("Time (days)")+ 
    ylab("GC volume")+theme_bw()+
   theme(axis.title.x = element_text(size=22),axis.text.x = element_text(size=19),axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),legend.position = "none")+coord_cartesian(xlim=c(0,40))

print(c1)
dev.off()

# plot the plasma cell output from simulated GCs

png(file="indiv_npcs.png",width=4,height=4,units="in",res=200)

c3=ggplot(data = finalpcn_int, mapping = aes(x = finalpcn_int[,1], y = finalpcn_int[,2], color = as.factor(finalpcn_int[,3])),show.legend=FALSE) +
geom_ribbon(data=finalpcn_int_sd,aes(x=finalpcn_int_sd[,1],ymin=finalpcn_int[,2]-finalpcn_int_sd[,2],ymax=finalpcn_int[,2]+finalpcn_int_sd[,2],fill = as.factor(finalpcn_int[,3])),alpha = 0.3,colour=NA,show.legend=FALSE)+
    geom_line()+
       xlab("Time (days) ")+ 
    ylab("No. of PCs")+theme_bw()+
   theme(axis.title.x = element_text(size=22),axis.text.x = element_text(size=19),axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),legend.position = "none")+coord_cartesian(xlim=c(0,40))

print(c3)
dev.off()

