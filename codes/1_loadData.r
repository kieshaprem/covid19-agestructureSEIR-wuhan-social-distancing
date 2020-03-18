require(data.table)
require(Matrix)
require(matrixcalc)

## load data: Wuhan's population age structure and Contact matrices

loadPopData = TRUE
loadContactMatrices = TRUE
loadCaseData =TRUE
loadR0posterior =TRUE


## load data: Population age structure (currently using China's) and Contact matrices

# 1) population data
if(loadPopData) 
{ 
  wuhanpop = read.csv('data/wuhanpop.csv',as.is = TRUE)
}

# 2) (projected) contact matrices 
### Acknowlegments: codes from Petra Klepac (petrakle) 
normalize.contact.matrices <- function(C, popv, make.sym = F){
  # FUN normalize them so that
  # the matrix with all the contacts is normalized so that its dominant eigenvalue is 1 
  # and other matrices keep their contributions 
  if (make.sym){
    Csym <- lapply(C, function(x, popv) (x + t(x)*((popv)%*%t(1/popv)))/2, popv) # make sure contacts are reciprocal
  } else {
    Csym <- C # if make.sym = F leave it as is
  }
  eig1 <- Re(eigen(Csym["all"]$all)$value[1])  # save dominant eigenvalue of the matrix with all contacts
  
  # divide all matrices by the real part of the dominant matrix of all of the contacts
  # that way all the rescaled matrices still sum to C_all = C_work + C_home + C_school + C_other
  Cnorm <- lapply(Csym, function(x,eig1) x/eig1, eig1)
  
  return(Cnorm)
}

# Synthetic contact matrices for China from Prem, Cook and Jit (2017) https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697   
if(loadContactMatrices)
{
  load(paste0('data/contacts_china.rdata'))
  contacts <- contacts_china # normalize.contact.matrices(contacts_china,wuhanpop$popage, make.sym=T)
  rm(contacts_china)
}

# case age distribution
if(loadCaseData)
{
  wuhancaseraw = read.csv('data/wuhan_pop_case_dist.csv',as.is = TRUE)
  caseage = rep(wuhancaseraw$wuhan,each=2)/2
  wuhancase = data.frame(agegroup = 1:16, caseage = c(caseage[1:15],sum(caseage[16:20])))
  rm(wuhancaseraw,caseage)
}

if(loadR0posterior)
{
  # --- read in R0 posterior
  R0_plot <-read.csv(paste0("data/out_R0.csv"))
  R0_dates <- read.csv(paste0('data/out_date.csv'))
  start_date <- as.Date(R0_dates[1,1]) # first case
  end_date <- as.Date(R0_dates[nrow(R0_dates),1]) # period to forecast ahead
  date_range <- seq.Date(start_date,end_date,1)
  
  # extract all estimates from 01.01.2020 - 23.01.2020
  R0_posterior <- R0_plot[which(date_range == as.Date("2020-01-01") ):which(date_range == as.Date("2020-01-23")),]
  range(R0_posterior)
  r0posterior = as.vector((unlist(R0_posterior)))
  par(mfrow=c(2,1))
  R0_dense = (density((r0posterior)))
  plot(x = R0_dense$x,y=R0_dense$y,type='l',xlab='R0',ylab='Density',lwd=2)
  R0_dense = (density(log(r0posterior)))
  plot(x = R0_dense$x,y=R0_dense$y,type='l',xlab='ln(R0)',ylab='Density',lwd=2)
  
  
  rm(R0_dense,R0_plot,R0_posterior,date_range,end_date,start_date)
  
}


rm(loadContactMatrices,loadPopData,loadR0posterior,loadCaseData)
