# script to plot overall GC volume


# data from Hollowood, Eur J Immunol 1992: (i.p, SRBC - 2*10^8, suspended in 0.25 ml sterile saline): readout: GCvolume in mm^3 
# data from Wang, Immunity 2005: (i.p, SRBC-0.1 ml):readout:% PNA+ Fas+ cells

# read simulation data
tt_gcvol <- read.delim("tt_gcvol_mean",header=FALSE,sep=" ")
tt_gcvol[,2]=tt_gcvol[,2]/(10^5)
tt_gcvol <- data.frame(tt_gcvol)

tt_gcvol_sd <- read.delim("tt_gcvol_sd",header=FALSE,sep=" ")
tt_gcvol_sd[,2]=tt_gcvol_sd[,2]/(10^5)

tt_gcvol_sd <- data.frame(tt_gcvol_sd)


# data from hollowood
data_hollowood_gcvol <- read.delim("Hollowood_gcvolume.txt",header=TRUE,sep="")
zero_hollowood <- 1.16

# data from Wang et al, Immunity 2005  

data_wang_gcvol <- read.delim("Wang_gcvolume.txt",header=TRUE,sep="")                                                                                                                                                                  
zero_wang <- 1.03

# rescale the data points to simulation
# rescale to the 7th day in the simulation

index_resc <- which(tt_gcvol[,1]<7.05 & tt_gcvol[,1]>6.95)[[1]]
value_resc <- tt_gcvol[index_resc,2]

hollowood_gcvol_mean_res <- ((data_hollowood_gcvol$mean_gcvol-zero_hollowood)*value_resc)/(3.01-zero_hollowood) 
hollowood_gcvol_ueb_res <- ((data_hollowood_gcvol$upperEB_gcvol-zero_hollowood)*value_resc)/(3.01-zero_hollowood)
hollowood_gcvol_sd_res <- hollowood_gcvol_ueb_res-hollowood_gcvol_mean_res

# rescale Wang data
# rescale 9.3 to value_resc

#gcvol_wang_res <- (gcvol_wang*value_resc)/9.3 
#EBu_wang_res <- (EBu_wang*value_resc)/9.3

wang_gcvol_mean_res <- ((data_wang_gcvol$mean_gcvol-zero_wang)*value_resc)/(9.3-zero_wang) 
wang_gcvol_ebu_res <- ((data_wang_gcvol$upperEB_gcvol-zero_wang)*value_resc)/(9.3-zero_wang)
wang_gcvol_sd_res <- wang_gcvol_ebu_res-wang_gcvol_mean_res

# combine the data of Hollowood and wang into a dataframe

finaldata_gcvol_time <- c(data_hollowood_gcvol$time,data_wang_gcvol$time)
finaldata_gcvol_mean <- c(hollowood_gcvol_mean_res,wang_gcvol_mean_res)
finaldata_gcvol_sd <- c(hollowood_gcvol_sd_res,wang_gcvol_sd_res)
finaldata_gcvol_classes <- c(rep(1,length(data_hollowood_gcvol$time)),rep(2,length(data_wang_gcvol$time)))

finaldata_gcvol <- cbind(finaldata_gcvol_time,finaldata_gcvol_mean,finaldata_gcvol_sd,finaldata_gcvol_classes)
finaldata_gcvol <- data.frame(finaldata_gcvol)


# extract simulation results corresponding to timepoints in experimental data

data_tmp <- NULL
time_tmp <- NULL
data_gcvol <- NULL

for(loop in 1:nrow(finaldata_gcvol)){

    index_extract<-which(abs(tt_gcvol[,1]-finaldata_gcvol[loop,1])==min(abs(tt_gcvol[,1]-finaldata_gcvol[loop,1])))
    data_tmp <- c(data_tmp,tt_gcvol[index_extract,2])
    time_tmp <- c(time_tmp,tt_gcvol[index_extract,1])

}

data_gcvol <- cbind(finaldata_gcvol[,1:2],time_tmp,data_tmp)

write.table(data_gcvol,file="../data_gcvol", sep = " ", row.names = FALSE, col.names = FALSE)


# plot overall GC volume with data

png(file="total_gcvolume.png",width=4, height=4, units="in", res=200)

p <- ggplot() +
   geom_point(data=finaldata_gcvol, aes(x=finaldata_gcvol[,1], y=finaldata_gcvol[,2],color=as.factor(finaldata_gcvol[,4])))+
   xlab("Time (days)")+
      ylab(bquote(.("GC volume")~"(x"*10^-5*")"))+
   #ylab(" GC volume (x 10) ")+
   geom_errorbar(data=finaldata_gcvol,aes(x=finaldata_gcvol[,1],ymin=finaldata_gcvol[,2]-finaldata_gcvol[,3], ymax=finaldata_gcvol[,2]+finaldata_gcvol[,3],color=as.factor(finaldata_gcvol[,4])),width=1.0,position=position_dodge(0.005))+ theme_bw() +
   theme(axis.title.y = element_text(size=22),axis.text.y  = element_text(size=19),axis.title.x = element_text(size=22),axis.text.x  = element_text(size=19),legend.position = "none",legend.title=element_blank())+geom_line(data = tt_gcvol ,aes(x=tt_gcvol[,1], y=tt_gcvol[,2],color="red"))+geom_ribbon(data=tt_gcvol_sd,aes(x=tt_gcvol[,1],ymin=tt_gcvol[,2]-tt_gcvol_sd[,2],ymax=tt_gcvol[,2]+tt_gcvol_sd[,2]),alpha = 0.15,fill = "red",colour=NA,show.legend=FALSE)+ scale_color_manual(values = c("black", "blue", "red"),labels=c("Hollowood data","Wang data","simulation"))+coord_cartesian(xlim=c(0,40))#+guides(shape = guide_legend())
  
print(p)                 
dev.off()


                 
