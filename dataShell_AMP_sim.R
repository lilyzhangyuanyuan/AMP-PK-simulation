##################################################################################
#The code of this paper includes three files:
# (1)"dataShell_AMP_sim.R" : simulating the infusion and non-infusion visits 
#                            according to the AMP protocol specification with 
#                            different adherences and sample sizes.
# (2)"AMP_ConcentrationSim_hvtn104IVFinal.ctl": NONMEM control file simulating concentrations of AMP 
#                                               participants by using the master
#                                               popPK model, which is the model 
#                                               that best describes the HVTN104
#                                               data.
# (3)"AMP_estimate_PK_parameters": NONMEM control file estimating popPK model parameters based on the 
#                                  simulated concentration of AMP participants.
#
##################################################################################

rm(list=ls())
library(plyr)
library(reshape2)
library(parallel)

###################
#simulation function
##################
#' Simulation
#' @param n_infu  : number of infusions
#' @param B  : number of datasets simulated 
#' @param n_ppt  : number of participants
#' @param prob_misInfu   : probability of an independently missed single infusion
#' @param prob_misVis : probability of an independently missed non-infusion visit
#' @param r_infuDiscont : cumulative probability of permanent infusion discontinuation
#' @param r_dropout  : annual dropout rate



sim <- function(n_infu=10
                ,B=1000
                ,n_ppt
                ,prob_misInfu
                ,prob_misVis
                ,r_infuDiscont
                ,r_dropout){
  #### redefine some parameters #####
  n_male <- n_ppt/2
  n_female <- n_ppt/2
  prob_misInfu <- rep(prob_misInfu,n_infu)
  prob_outWin <- rep(0.2,n_infu)
  prob_disContInfu_Win <- c(0.2/2,1-0.2,0.2/2)
  #directory to save the current scenario data
  folder <- unique(paste0("m_",n_ppt,"_p1_",prob_misInfu,"_p2_",prob_misVis,"_r1_",r_infuDiscont,"_r2_",r_dropout))
  dir <- file.path("/fh_fast/AMP_PK_Sim_new/nonmem/dataShell",folder)
  
  if(file.exists(dir)==FALSE){
    dir.create(dir)
  }
  ##### start the simulation #####
  lapply(1:B,function(b){
    #-------- infusion visits -------------------
    dat <- data.frame(ID=1:n_ppt)
    dat$infu1 <- 1
    for(i in 2:n_infu){
      var_name <- paste0("infu",i)
      dat[,var_name] <- rbinom(n=n_ppt,size=1,prob=1-prob_misInfu[i])
    }
    #long format
    dat <- melt(dat,id.vars="ID")
    dat <- dat[order(dat$ID),]
    dat <- rename(dat,c("value"="f_infu","variable"="infusion"))
    dat$infusion <- as.numeric(gsub("infu","",dat$infusion))
    #out of target window but within allowable window
    dat <- ddply(dat, .(infusion), function(df){
      if(unique(df$infusion)==1){
        df$out_win <- NA
      }else{
        infu <- unique(df$infusion)
        dat_infu <- subset(df,f_infu==1)
        dat_infu$out_win <- rbinom(n=nrow(dat_infu),size=1,prob=prob_outWin[infu])
        dat_misInfu <- subset(df,f_infu==0)
        if(nrow(dat_misInfu)==0){
          df <- dat_infu
        }else{
          dat_misInfu$out_win <- NA
          df <- rbind(dat_infu,dat_misInfu)
        }
      }
      return(df)
    })
    #random draw from windows
    dat_inWin <- subset(dat,out_win==0)
    dat_inWin$window <- round(runif(n=nrow(dat_inWin),min=-7,max=7),0)
    dat_outWin <- subset(dat,out_win==1)
    dat_outWin$window <- round(runif(n=nrow(dat_outWin),min=8,max=48),0)
    dat_naWin <- subset(dat,is.na(out_win))
    dat_naWin$window <- 0
    dat <- rbind(dat_inWin,dat_outWin,dat_naWin)
    dat <- dat[order(dat$ID,dat$infusion),]
    #infusion day
    dat <- ddply(dat,.(ID),function(df){
      df$day_infu[df$infusion==1] <- 0
      for(i in 2:n_infu){
        df$day_infu[df$infusion==i] <- df$day_infu[df$infusion==(i-1)]+56+df$window[df$infusion==i]
      }
      return(df)
    })
    dat$day_infu <- ifelse(dat$f_infu==0,NA,dat$day_infu)
    dat$out_win <- NULL
    dat$window <- NULL
    #------------ post-infusion visits ------------------
    dat$window <- round(runif(n=nrow(dat),min=-7,max=7),0)
    dat$day_4wkPostInfu <- dat$day_infu+28+dat$window
    dat$day_8wkPostInfu10 <- dat$day_infu+56+dat$window
    dat$window_5dayPostInfu2 <- round(runif(n=nrow(dat),min=-2,max=2),0)
    dat$day_5dayPostInfu2 <- dat$day_infu+5+dat$window_5dayPostInfu2
    #missing visit 
    dat <- ddply(dat,.(infusion),function(df){
      nppt_infu <- length(df$f_infu[df$f_infu==1])
      df$f_4wkPostInfuVis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)
      if(unique(df$infusion)==2){
        df$f_5dayPostInfu2Vis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)
      }else{
        df$f_5dayPostInfu2Vis <- 0
      }
      if(unique(df$infusion)==10){
        df$f_8wkPostInfu10Vis <- rbinom(n=nrow(df),size=1,prob=1-prob_misVis)
      }else{
        df$f_8wkPostInfu10Vis <- 0
      }
      return(df)
      
    })
    dat$day_4wkPostInfu <- ifelse(dat$f_4wkPostInfuVis==0,NA,dat$day_4wkPostInfu)
    dat$day_8wkPostInfu10 <- ifelse(dat$f_8wkPostInfu10Vis==0,NA,dat$day_8wkPostInfu10)
    dat$day_5dayPostInfu2 <- ifelse(dat$f_5dayPostInfu2Vis==0,NA,dat$day_5dayPostInfu2)
    dat$window <- NULL
    dat$window_5dayPostInfu2 <- NULL
    #------------ change data to long format for NONMAM ----------------
    #infusion visit
    dat_infu <- subset(dat,select=c(ID,infusion,f_infu,day_infu))
    dat_infu <- rename(dat_infu,c("infusion"="visit_n","f_infu"="f_vis","day_infu"="day_vis"))
    dat_infu$visit_n <- ifelse(dat_infu$visit_n%in%c(1,2),dat_infu$visit_n*2,dat_infu$visit_n*2+1)
    dat_infu$EVID <- 1
    #trough measurement
    dat_trgh <- subset(dat_infu, visit_n!=2)
    dat_trgh$EVID <- 0 #(0 for observation event, 1 for dose event)
    #4 week post infusion
    dat_4wkPostInfu <- subset(dat,select=c(ID,infusion,f_4wkPostInfuVis,day_4wkPostInfu))
    dat_4wkPostInfu <- rename(dat_4wkPostInfu,c("infusion"="visit_n"
                                                ,"f_4wkPostInfuVis"="f_vis"
                                                ,"day_4wkPostInfu"="day_vis"))
    dat_4wkPostInfu$visit_n <- ifelse(dat_4wkPostInfu$visit_n==1,3,dat_4wkPostInfu$visit_n*2+2)
    dat_4wkPostInfu$EVID <- 0
    #8 week post 10th infusion
    dat_8wkPostInfu10 <- subset(dat
                                ,infusion==10
                                ,select=c(ID,infusion,f_8wkPostInfu10Vis,day_8wkPostInfu10))
    dat_8wkPostInfu10 <- plyr::rename(dat_8wkPostInfu10,c("infusion"="visit_n"
                                                          ,"f_8wkPostInfu10Vis"="f_vis"
                                                          ,"day_8wkPostInfu10"="day_vis"))
    dat_8wkPostInfu10$visit_n <- 23
    dat_8wkPostInfu10$EVID <- 0
    # 5 days post 2nd infusion 
    dat_5dayPostInfu2 <- subset(dat
                                ,infusion==2
                                ,select=c(ID,infusion,f_5dayPostInfu2Vis,day_5dayPostInfu2))
    dat_5dayPostInfu2 <- rename(dat_5dayPostInfu2,c("infusion"="visit_n"
                                                    ,"f_5dayPostInfu2Vis"="f_vis"
                                                    ,"day_5dayPostInfu2"="day_vis"))
    dat_5dayPostInfu2$visit_n <- 5
    dat_5dayPostInfu2$EVID <- 0
    #combine
    dat_long <- rbind(dat_infu,dat_trgh,dat_5dayPostInfu2,dat_4wkPostInfu,dat_8wkPostInfu10)
    dat_long <- subset(dat_long,f_vis==1)
    dat_long <- dat_long[order(dat_long$ID,dat_long$visit_n),]
    #----------- permanent infusion discontinuation -------------------------------------------------
    #time to discontinuation of infusion
    if(r_infuDiscont==0){
      r_infuDiscont <- 0.00000001
    }
    time_discontInfu <- data.frame(ID=1:n_ppt,day_discontInfu=round(rexp(n=n_ppt,rate=r_infuDiscont/365),0))
    #don't allow discontiune at time 0
    time_discontInfu$day_discontInfu[time_discontInfu$day_discontInfu==0] <- 1
    dat_long <- merge(dat_long,time_discontInfu,by="ID")
    dat_long$f_discontInfu <- ifelse(dat_long$day_vis<dat_long$day_discontInfu,0,1)
    dat_long$f_vis <- ifelse(dat_long$f_discontInfu==1,0,dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$f_discontInfu==1, NA,dat_long$day_vis)
    # visits after the ppts discontinue infusions
    id_discontInfu <- unique(dat_long$ID[dat_long$f_discontInfu==1])
    dat_discont <- unique(dat_long[dat_long$ID%in%id_discontInfu,c("ID","day_discontInfu")])
    dat_discont <- data.frame(ID=rep(dat_discont$ID,each=6)
                              ,day_discontInfu=rep(dat_discont$day_discontInfu,each=6)
                              ,wk_post_enroll = rep(c("20wk","32wk","44wk","56wk","68wk","80wk"),nrow(dat_discont))
    ) 
    dat_discont$lw[dat_discont$wk_post_enroll=="20wk"] <- 99
    dat_discont$lw[dat_discont$wk_post_enroll=="32wk"] <- 183
    dat_discont$lw[dat_discont$wk_post_enroll=="44wk"] <- 267
    dat_discont$lw[dat_discont$wk_post_enroll=="56wk"] <- 351
    dat_discont$lw[dat_discont$wk_post_enroll=="68wk"] <- 435
    dat_discont$lw[dat_discont$wk_post_enroll=="80wk"] <- 519
    
    dat_discont$day_tar[dat_discont$wk_post_enroll=="20wk"] <- 140
    dat_discont$day_tar[dat_discont$wk_post_enroll=="32wk"] <- 224
    dat_discont$day_tar[dat_discont$wk_post_enroll=="44wk"] <- 308
    dat_discont$day_tar[dat_discont$wk_post_enroll=="56wk"] <- 392
    dat_discont$day_tar[dat_discont$wk_post_enroll=="68wk"] <- 476
    dat_discont$day_tar[dat_discont$wk_post_enroll=="80wk"] <- 560
    
    dat_discont$f_discont <- ifelse(dat_discont$day_discontInfu<dat_discont$lw,1,0)
    dat_discont <- subset(dat_discont,f_discont==1)
    
    win_type <- data.frame(t(rmultinom(nrow(dat_discont), size = 1, prob=prob_disContInfu_Win)))
    names(win_type) <- c("win_type1","win_type2","win_type3")
    dat_discont <- cbind(dat_discont,win_type)
    dat_discont$win1 <- round(runif(n=nrow(dat_discont),min=-41,max=-14),0)
    dat_discont$win2 <- round(runif(n=nrow(dat_discont),min=-14,max=14),0)
    dat_discont$win3 <- round(runif(n=nrow(dat_discont),min=14,max=42),0)
    dat_discont$win[dat_discont$win_type1==1] <- dat_discont$win1[dat_discont$win_type1==1]
    dat_discont$win[dat_discont$win_type2==1] <- dat_discont$win2[dat_discont$win_type2==1]
    dat_discont$win[dat_discont$win_type3==1] <- dat_discont$win3[dat_discont$win_type3==1]
    dat_discont$day_vis <- dat_discont$day_tar+dat_discont$win
    #missing visit
    dat_discont$f_vis <- rbinom(n=nrow(dat_discont),size=1,prob=1-prob_misVis)
    dat_discont$day_vis <- ifelse(dat_discont$f_vis==0,NA,dat_discont$day_vis)
    #visits for ppt who discontinue infusions 
    dat_discont <- subset(dat_discont,select=c(ID,wk_post_enroll,f_vis,day_vis))
    dat_discont <- rename(dat_discont,c("wk_post_enroll"="visit_n"))
    dat_discont$visit_n  <- as.character(dat_discont$visit_n)
    dat_discont$visit_n[dat_discont$visit_n=="20wk"] <- 72
    dat_discont$visit_n[dat_discont$visit_n=="32wk"] <- 73
    dat_discont$visit_n[dat_discont$visit_n=="44wk"] <- 74
    dat_discont$visit_n[dat_discont$visit_n=="56wk"] <- 75
    dat_discont$visit_n[dat_discont$visit_n=="68wk"] <- 76
    dat_discont$visit_n[dat_discont$visit_n=="80wk"] <- 77
    dat_discont$EVID[!is.na(dat_discont$ID)] <- 0 #for situation that no ppts who discontinue infusions
    #combine
    dat_long$day_discontInfu <- NULL
    dat_long$f_discontInfu <- NULL
    dat_long <- rbind(dat_long,dat_discont)
    dat_long <- subset(dat_long,f_vis==1)
    #------------------- drop out ----------------------
    if(r_dropout==0){
      r_dropout <- 0.000000001
    }
    time_dropout <- data.frame(ID=1:n_ppt, time_dropout=rexp(n=n_ppt,rate = r_dropout/365))
    dat_long <- merge(dat_long,time_dropout)
    dat_long$f_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,0, dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,NA, dat_long$day_vis)
    dat_long <- dat_long[order(as.numeric(dat_long$visit_n),dat_long$EVID),]
    #final data
    dat_long <- subset(dat_long,f_vis==1)
    dat_long$time_dropout <- NULL
    
    #-------- NONMAM datashell -------------------------
    #get the distribution of age from hvtn503 female and hvtn502 male 
    m_age_male <- 30.52342
    sd_age_male <- 7.860361
    m_age_female <- 30.52342
    sd_age_female <- 4.591396
    #WT data
    dataDirWT <- "/trials/vaccine/AMP_PK/AMP_PK_sim/data"
    dat_502_WT <- read.csv(file.path(dataDirWT,"v502_wt_list.csv"))
    dat_502_WT <- unique(dat_502_WT)
    dat_503_WT <- read.csv(file.path(dataDirWT,"503_blweight_women.csv"))
    #WT, age and dose
    dat_m <- data.frame(ID=1:n_male
                        ,WT=sample(dat_502_WT$WT_KLG,n_male,replace=T)
                        ,AGE = round(rnorm(n=n_male,mean=m_age_male,sd=sd_age_male),0)
                        ,sex=1
                        ,dose=c(rep(10,floor(n_male/2)),rep(30,ceiling(n_male/2)))
                        )
    
    dat_m$AGE <- ifelse(dat_m$AGE>50,50,dat_m$AGE)
    dat_m$AGE <- ifelse(dat_m$AGE<18,18,dat_m$AGE)
    dat_f <- data.frame(ID=(n_male+1):(n_female+n_male)
                        ,WT=sample(dat_503_WT$weight,n_female,replace=T)
                        ,AGE = round(rnorm(n=n_female,mean=m_age_female,sd=sd_age_female),0)
                        ,sex=0
                        ,dose=c(rep(10,floor(n_female/2)),rep(30,ceiling(n_female/2)))
                        )
   
    dat_f$AGE <- ifelse(dat_f$AGE>40,40,dat_f$AGE)
    dat_f$AGE <- ifelse(dat_f$AGE<18,18,dat_f$AGE)
    dat_cov <- rbind(dat_m,dat_f)
    #complete data with TIME on each day before last visit
    day_max <- ddply(dat_long,.(ID),summarise,day_max=max(day_vis))
    dat_com <- ddply(day_max,.(ID),function(df){
      dat <- data.frame(ID=df$ID,TIME=1:df$day_max)
      return(dat)
    })
    #merge visit information
    dat_shell <- subset(dat_long,select=c(ID,day_vis,EVID,f_vis,visit_n))
    dat_shell <- rename(dat_shell,c("day_vis"="TIME"))
    dat_shell_infu <- subset(dat_shell,EVID==1)
    dat_shell_sample <- subset(dat_shell,EVID==0)
    dat_shell_sample <- merge(dat_com,dat_shell_sample,by=c("ID","TIME"),all.x=T)
    dat_shell_sample$EVID <- 0
    dat_shell_sample$f_vis[is.na(dat_shell_sample$f_vis)] <- 0
    dat_shell <- rbind(dat_shell_infu,dat_shell_sample)
    dat_shell <- merge(dat_shell,dat_cov,by="ID")
    dat_shell$AMT <- ifelse(dat_shell$EVID==1,dat_shell$WT*dat_shell$dose,0)
    dat_shell$DV <- ifelse(dat_shell$EVID==0,1,0)
    dat_shell$RATE <- dat_shell$AMT*24*(60/45)
    dat_shell$visit_n <- as.numeric(dat_shell$visit_n)
    dat_shell <- dat_shell[order(dat_shell$ID,dat_shell$TIME,dat_shell$EVID),]
    dat_shell <- dat_shell[,c("ID","TIME","AMT","RATE","DV","WT","AGE","sex","f_vis","dose","visit_n")]
    dat_shell$boot <- b
    #save data
    csv.file <- paste0(b,".csv")
    write.table(dat_shell,file.path(dir,csv.file)
                ,na=".",row.names = FALSE
                ,col.names = FALSE
                ,sep=",")
    
  })
}

###############################################
#start simulations under different scenarios
###############################################
n_ppt <- list(30,60,120,240)
scenario <- list(c(0,0,0,0)
                 ,c(0.05,0.1,0.1,0.1)
                 ,c(0.1,0.15,0.15,0.15)
                 ,c(0.15,0.2,0.2,0.2)
                 )
scenario <- expand.grid(n_ppt,scenario)
scenario <- apply(scenario,2,function(x)do.call("rbind",x))
scenario <- do.call("cbind",scenario)
scenario <- data.frame(scenario)
names(scenario) <- c("n_ppt","prob_misInfu","prob_misVis","r_infuDiscont","r_dropout")

set.seed(111)
for(i in 1:nrow(scenario)){
  sim(n_ppt=scenario$n_ppt[i]
      ,prob_misInfu=scenario$prob_misInfu[i]
      ,prob_misVis=scenario$prob_misVis[i]
      ,r_infuDiscont=scenario$r_infuDiscont[i]
      ,r_dropout=scenario$r_dropout[i])
}


