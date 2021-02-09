
rm(list = ls())
options(scipen=7)
require('deSolve')
library(ggplot2)
library(openxlsx)
require(gridExtra)

## simulated singapore imported cases
simulate_input <- read.xlsx('sim_singapore_20200525.xlsx')


## SEIR function
SEIR_func=function(t, x, vparameters){
  S = x[1]  
  E = x[2]  
  I = x[3] 
  R = x[4]
  N=S+E+I+R
  with(as.list(vparameters),{
    #npop = S+E+I+R
    dS = (-R0/DI * S/N * I)
    dE = (R0/DI * S/N * I) - (1/DE * E)
    dI = (1/DE * E) - (1/DI * I)
    dR = (1/DI * I)
    vout = c(dS,dE,dI,dR)
    list(vout)
  })
}

## set date
date_real <- "2020/02/25"
simulate_input$date <- seq.Date(from=as.Date(date_real,format = "%Y/%m/%d"), by = "day",length.out = nrow(simulate_input)) 


date0 <- "2020/03/04"  ## start of simulation period
row0<-which(simulate_input$date==date0)
tmin = 1
tmax=nrow(simulate_input)-row0+1
vt = seq(tmin,tmax)
date1 <- "2020/04/08" ## quarantine policy changed 
row1<-which(simulate_input$date==date1)
t1 <- row1-row0+1 

simulate_input <- simulate_input[c(row0:nrow(simulate_input)),]


# initials
npop = 5700000  # population
R_0=79 # removed cases
I_0<-0 # infected local cases
# parameters
DI <- 2.3 # infectious period
DE <- 5.2 # incubation period
prop_a <- 0.31 # asymptomatic fraction
prop_e<-0.643 # presymptomatic fraction

# Rt
Rt_df <- read.csv('estimate_rt_sg.csv',as.is = T)

## inbound traveller reduction
k1<-0.988
k2<-0.554
k3<-0.305

## scenario simulation

## 1.three scenarios of quarantine policies with the estimated time-varying Rt and actual inbound traveller reduction


# quarantine all inbound travellers
inputcase<-simulate_input$simulate

input1<-list(input = data.frame(new_input=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,
                                nfree_I=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,  ## not free
                                free_Is=0,   ## symptomatic free of quarantine
                                free_E=0,    ## presymptomatic free of quarantine
                                free_Ia=0),  ## asymptomatic free of quarantine
             
             I_0 = I_0,
             R_0=R_0)
input1[['E_0']]<-input1[['I_0']]*(1-prop_a)*prop_e
input1[['S_0']] = npop-input1[['I_0']]-input1[['R_0']]-input1[['E_0']]
# quarantine symptomatic travellers
input2<-list(input = data.frame(new_input=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,
                                nfree_I=inputcase*(1-prop_a)*(1-prop_e),  ## not free
                                free_Is=0,                                ## symptomatic free of quarantine
                                free_E=inputcase*(1-prop_a)*prop_e,       ## presymptomatic free of quarantine
                                free_Ia=inputcase*prop_a),                ## asymptomatic free of quarantine
             
             I_0=I_0,
             R_0=R_0)
input2[['E_0']]<-input2[['I_0']]*(1-prop_a)*prop_e
input2[['S_0']] = npop-input2[['I_0']]-input2[['R_0']]-input2[['E_0']]


# no quarantine required
input3<-list(input = data.frame(new_input=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,
                                nfree_I=0,  ##not free
                                free_Is=inputcase*(1-prop_a)*(1-prop_e),  ## symptomatic free of quarantine
                                free_Ia=inputcase*(prop_a),               ## asymptomatic free of quarantine
                                free_E=inputcase*(1-prop_a)*prop_e),      ## presymptomatic free of quarantine
             
             I_0=I_0,
             R_0=R_0)
input3[['E_0']]<-input3[['I_0']]*(1-prop_a)*prop_e
input3[['S_0']] = npop-input3[['I_0']]-input3[['R_0']]-input3[['E_0']]

input_all <- list(input1,input2,input3)



# simulation
results <- matrix(NA,nrow = length(vt),ncol = 7)
colnames(results) <- c("S","E","I","R","Rt","new_input","new_case")
results_all <-list()

for (i in c(1:length(input_all))){
  results[1,c(1:7)] <- c(input_all[[i]]$S_0,input_all[[i]]$E_0,input_all[[i]]$I_0,R_0,Rt_df$rt[1],input_all[[i]]$input$new_input[1],input_all[[i]]$input$new_input[1])
  results <- data.frame(results,stringsAsFactors = F)
  
  for (t in c(2:tmax)){
    results[t,"new_input"] <- input_all[[i]]$input$new_input[t]
    vparameters = c(R0=Rt_df$rt[t],DI=DI,DE=DE)
    inits <- c(S=results$S[t-1],
               E=results$E[t-1]+input_all[[i]]$input$free_E[t],
               I=results$I[t-1]+input_all[[i]]$input$free_Is[t]+input_all[[i]]$input$free_Ia[t],
               R=results$R[t-1])
    
    solved_model = as.data.frame(lsoda(inits, seq(1,tmax), SEIR_func, vparameters))
    
    results[t,"S"] <- solved_model$S[2]
    results[t,"E"] <- solved_model$E[2]
    results[t,"I"] <- solved_model$I[2]
    results[t,"R"] <- solved_model$R[2]
    results[t,'Rt'] <- Rt_df$rt[t]
    results[t,"new_case"] <- results[t,"I"]+results[t,"R"]-
      (results[t-1,"I"]+results[t-1,"R"])+input_all[[i]]$input$nfree_I[t]
    results[t,'under_case_all'] <- ((results[t,'new_case'])^2)/(1+results[t-1,'new_case'])            ## under report cases
    results[t,'under_case_import'] <- ((results[t,'new_input'])^2)/(1+results[t-1,'new_input'])       ## under report cases
    results[t,"all_new_case"] <- results[t,"new_case"]+results[t,"under_case_all"]
    results[t,"all_new_case_import"] <- results[t,"new_input"]+results[t,"under_case_import"]
    results[t,"all_new_case_local"] <- results[t,"all_new_case"]-results[t,"all_new_case_import"]
    
  }
  
  results$Date <- seq.Date(from = as.Date(date0,format = "%Y/%m/%d"), by = "day", length.out = tmax)
  results1 <- results[,c(13,1:12)]
  name <- paste('results',i,sep='_')
  results_all[[name]]<-results1
}


## 2.three scenarios of intensity under the estimated time-varying Rt and actual quarantine policy

I_0=I_0
R_0=R_0
E_0<-I_0*(1-prop_a)*prop_e
S_0 = npop-I_0-R_0-E_0

for(k in c(k1,k2,k3)){
  inputcase<-simulate_input$`simulate-no-restriction`*(1-k)
  results <- matrix(NA,nrow = length(vt),ncol = 7)
  colnames(results) <- c("S","E","I","R","Rt","new_input","new_case")
  results[1,c(1:7)] <- c(S_0,E_0,I_0,R_0,Rt_df$rt[1],inputcase[1],inputcase[1])
  results <- data.frame(results,stringsAsFactors = F)
  ##quarantine symptomatic
  for (t in c(2:t1)){
    input = data.frame(new_input=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,
                       nfree_I=inputcase*(1-prop_a)*(1-prop_e),  ##not free
                       free_Is=0,
                       free_E=inputcase*(1-prop_a)*prop_e,
                       free_Ia=inputcase*prop_a)
    results[t,"new_input"] <- input$new_input[t]
    vparameters = c(R0=Rt_df$rt[t],DI=DI,DE=DE)
    inits <- c(S=results$S[t-1],
               E=results$E[t-1]+input$free_E[t],
               I=results$I[t-1]+input$free_Is[t]+input$free_Ia[t],
               R=results$R[t-1])
    
    solved_model = as.data.frame(lsoda(inits, seq(1,t1), SEIR_func, vparameters))
    
    results[t,"S"] <- solved_model$S[2]
    results[t,"E"] <- solved_model$E[2]
    results[t,"I"] <- solved_model$I[2]
    results[t,"R"] <- solved_model$R[2]
    results[t,'Rt'] <- Rt_df$rt[t]
    results[t,"new_case"] <- results[t,"I"]+results[t,"R"]-
      (results[t-1,"I"]+results[t-1,"R"])+input$nfree_I[t]
    results[t,'under_case_all'] <- ((results[t,'new_case'])^2)/(1+results[t-1,'new_case'])
    results[t,'under_case_import'] <- ((results[t,'new_input'])^2)/(1+results[t-1,'new_input'])
    results[t,"all_new_case"] <- results[t,"new_case"]+results[t,"under_case_all"]
    results[t,"all_new_case_import"] <- results[t,"new_input"]+results[t,"under_case_import"]
    results[t,"all_new_case_local"] <- results[t,"all_new_case"]-results[t,"all_new_case_import"]
  }
  ## quarantine all inbound travellers
  for (t in c((t1+1):tmax)){
    input = data.frame(new_input=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,
                       nfree_I=inputcase*(1-prop_a)*(1-prop_e)+inputcase*prop_a,  ##not free
                       free_Is=0,
                       free_E=0,
                       free_Ia=0)
    results[t,"new_input"] <- input$new_input[t]
    vparameters = c(R0=Rt_df$rt[t],DI=DI,DE=DE)
    inits <- c(S=results$S[t-1],
               E=results$E[t-1]+input$free_E[t],
               I=results$I[t-1]+input$free_Is[t]+input$free_Ia[t],
               R=results$R[t-1])
    
    solved_model = as.data.frame(lsoda(inits, seq(1,tmax-t1+1), SEIR_func, vparameters))
    
    results[t,"S"] <- solved_model$S[2]
    results[t,"E"] <- solved_model$E[2]
    results[t,"I"] <- solved_model$I[2]
    results[t,"R"] <- solved_model$R[2]
    results[t,'Rt'] <- Rt_df$rt[t]
    results[t,"new_case"] <- results[t,"I"]+results[t,"R"]-
      (results[t-1,"I"]+results[t-1,"R"])+input$nfree_I[t]
    results[t,'under_case_all'] <- ((results[t,'new_case'])^2)/(1+results[t-1,'new_case'])
    results[t,'under_case_import'] <- ((results[t,'new_input'])^2)/(1+results[t-1,'new_input'])
    results[t,"all_new_case"] <- results[t,"new_case"]+results[t,"under_case_all"]
    results[t,"all_new_case_import"] <- results[t,"new_input"]+results[t,"under_case_import"]
    results[t,"all_new_case_local"] <- results[t,"all_new_case"]-results[t,"all_new_case_import"]
    
  }
  
  results$Date <- seq.Date(from = as.Date(date0,format = "%Y/%m/%d"), by = "day", length.out = tmax)
  results1 <- results[,c(13,1:12)]
  name <- paste('results',k,sep='_')
  results_all[[name]]<-results1

}

##cumulative infectious cases
cumulative_all <- matrix(NA,nrow = length(results_all),ncol = 1)
cumulative_all <- data.frame(cumulative_all,stringsAsFactors = F)

for (i in c(1:length(results_all))){
  cumulative_all[i,1] <- sum(results_all[[i]][,11],na.rm = T)
}
