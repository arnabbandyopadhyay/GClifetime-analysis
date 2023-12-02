# main script used for GC lifetime paper

library(ggplot2)


nsims=50 # number of simulation repeats 


# output files are collected in the results folder

system('mkdir results')
system('cd results')
setwd(paste0(getwd(),'/results'))

files=c("tt_gcvol","tt_pcs","tt_PCaff","GCvolume","nGCs","Plasmacells")

Agn_files=c("PCaffinity_Ag")

for(j in 1:nsims){

   system(paste0('cp ../',sprintf("%03d", j-1),'/ini_time.out',' ','ini_time','_',j))

    
    for(f in 1:length(files)){
      
        system(paste0('cp ../',sprintf("%03d", j-1),'/',files[f],'.out',' ',files[f],'_',j))

    }
}


# get the number of simulated GCs and number of antigens 

ngcs_read <- read.delim("nGCs_1",header=FALSE,sep="")
nags_read <- read.delim("tt_PCaff_1",header=FALSE,sep="")

ngcs = nrow(ngcs_read)
nags = ncol(nags_read)-1

for(j in 1:nsims){     
        
        for(f in 1:length(Agn_files)){
        
          for(agn in 1:nags){
              system(paste0('cp ../',sprintf("%03d", j-1),'/',Agn_files[f],agn-1,'.out',' ',Agn_files[f],agn-1,'_',j))
          }
        
          
        }    

}

allfiles = c("tt_gcvol","tt_pcs","tt_PCaff","GCvolume","Plasmacells")

for(i in 1:length(Agn_files)){
  for(j in 1:nags){
      allfiles=c(allfiles,paste0(Agn_files[i],j-1))
  }
}

print(allfiles)

# calculate mean and SD of simulation repeats

for(f1 in 1:length(allfiles)){

    calc_mean <- NULL
    calc_sd <- NULL
    
    read <- read.delim(paste0(allfiles[f1],'_',1),header=FALSE,sep="")
    
    calc_time <- read[,1]
    calc_mean <- cbind(calc_mean,calc_time)
    calc_sd <- cbind(calc_sd,calc_time)
    
    for(j in 2:ncol(read)){
        
        tmp <- NULL
    
        for(i in 1:nsims){
        
        read <- read.delim(paste0(allfiles[f1],'_',i),header=FALSE,sep="")
        
        tmp <- cbind(tmp,read[,j])
        
        }
                
        if(nsims==1){
        
            tmp_mean <- tmp
            tmp_sd <- c(rep(0,length(tmp)))
            
        }else{
        
            tmp_mean<- rowMeans(tmp)
            tmp_sd <- apply(tmp[,1:nsims],1,sd)
        
        }
        
        calc_mean <- cbind(calc_mean,tmp_mean)
        calc_sd <- cbind(calc_sd,tmp_sd)
  
    
    }
    
    write.table(data.frame(calc_mean),file=paste0(allfiles[f1],'_mean'), sep = " ", row.names = FALSE, col.names = FALSE)
    write.table(data.frame(calc_sd),file=paste0(allfiles[f1],'_sd'), sep = " ", row.names = FALSE, col.names = FALSE)

}


system('cp ../*.R .')
system('cp ../*.txt .')


# Overall GC volume
source('GCvolume.R')

# count number of GCs
source('Classify_new.R')

# calculate GC lifetime
source('Lifetimes_new.R')

# cumulative number of GCs initialized
source('Plot_ngcs_cum.R')

# cumulative number of GCs terminated
source('Plot_death_cum.R')

# Individual germinal centre readouts
source('IndivGCs.R')



