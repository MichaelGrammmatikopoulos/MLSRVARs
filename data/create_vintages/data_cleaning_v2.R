# Forecasting with Machine Learning Shadow-Rate VARs, Michael Grammatikopoulos (2023)
# 
# Data cleaning code for the ALFRED data vintages
#
# This code comes without technical support of any kind.  It is expected to
# reproduce the data-sets that are used in the paper. Under no circumstances will
# the author be held responsible for any use (or misuse) of this code in
# any way.

# WARNING: ALL VINTAGE DATASETS SHOULD START WITH  19XX-01-01 to avoid misalignment of data

rm(list=ls())
path_link<-"T:/PROJECTS/_MG_docs/severity/Michael/Final_Program_Round_1/matlab/MLSVARs/ALFRED_raw"
#Install some packages
source(paste0(path_link,"/00_install_packages.R",sep=""))

load("T:/PROJECTS/_MG_docs/severity/Michael/Final_Program_Round_1/matlab/MLSVARs/ALFRED_raw/raw_data.RData")

###########################################################################################################

gdp_data<-merge(gdp_raw1,gdp_raw2)
flbr_data<-merge(merge(flbr_raw1,flbr_raw2),flbr_raw3)
cpi_data<-merge(merge(cpi_raw1,cpi_raw2),cpi_raw3)
fedfunds_data<-fedfunds_raw1
payrolls_data<-merge(payrolls_raw1,payrolls_raw2)
houst_data<-merge(merge(houst_raw1,houst_raw2),houst_raw3)
fip_data<-merge(merge(fip_raw3,fip_raw4),fip_raw5)
fcu_data<-merge(fcu_raw1,fcu_raw2)
ftfx_data<-merge(merge(ftfx_raw1,ftfx_raw2),ftfx_raw3)
m2_data<-merge(merge(merge(m2_raw1,m2_raw2),m2_raw3),m2_raw4)
pincome_data<-merge(pincome_raw1,pincome_raw2)
pcomp_data<-pcomp_raw1
psavert_data<-psavert_raw1
ppiaco_data<-ppiaco_raw1
frgt3M_data<-merge(frgt3M_raw1,frgt3M_raw2)
frgt6M_data<-frgt6M_raw1
frgt1y_data<-frgt1y_raw1
frgt3y_data<-merge(frgt3y_raw1,frgt3y_raw2)
frgt5y_data<-frgt5y_raw1
frgt10y_data<-frgt10y_raw1
baamoodys_data<-merge(baamoodys_raw1,baamoodys_raw2)

###########################################################################################################

rm(list=ls(pattern="*_raw*"))  
raw_data_list<-list(    gdp_data,  flbr_data,   cpi_data,  fedfunds_data,  payrolls_data,  houst_data,  fip_data,  fcu_data, 
                        ftfx_data,   m2_data,  ppiaco_data, frgt3M_data, frgt6M_data, frgt1y_data, frgt3y_data,  frgt5y_data,  frgt10y_data, baamoodys_data)
names(raw_data_list)<-c("gdp_data","flbr_data","cpi_data","fedfunds_data","payrolls_data","houst_data","fip_data","fcu_data",
                        "ftfx_data","m2_data","ppiaco_data", "frgt3M_data", "frgt6M_data","frgt1y_data", "frgt3y_data","frgt5y_data","frgt10y_data","baamoodys_data")

for(name_i in names(raw_data_list)){
  if(name_i!="gdp_data"){
    
    data<-as.data.frame(raw_data_list[name_i])
    colnames_backup<-colnames(data)
    vintage_dates<-as.numeric(substr(colnames_backup,nchar(colnames_backup)-7,nchar(colnames_backup)))
      
    colnames(data)<-substr(str_replace_all(colnames(data),paste0(name_i,".",sep=""),""),1,nchar(colnames(data))-2)
    colnames(data)<-substr(colnames(data),1,nchar(colnames(data))-2)
    colnames(data)<-paste0(substr(colnames(data),1,nchar(colnames(data))-2),"m",substr(colnames(data),nchar(colnames(data))-1,nchar(colnames(data))),sep="")
    colnames(data)[1]<-"observation_date"
    ALFRED_mnemonic<-substr(colnames(data[-1][1]),1,nchar(colnames(data[-1][1]))-7)
    
    clean_columns_data<-as.data.frame(data[,which(colnames(data)=="observation_date")])
    colnames(clean_columns_data)[1]<-"observation_date"
    #ncol(fip_data)
    counter<-1
    for(j in 2:ncol(data)){
      if(str_detect(colnames(data)[j],"m01|m04|m07|m10")==TRUE){
        clean_columns_data_i<-as.data.frame(data[,j])
        colnames(clean_columns_data_i) <- colnames(data[j])
        clean_columns_data<-cbind(clean_columns_data,clean_columns_data_i)
        counter<-counter+1
        if(colnames(clean_columns_data)[counter]==colnames(clean_columns_data)[counter-1] & vintage_dates[j]>vintage_dates[j-1]){
          clean_columns_data<-clean_columns_data[,-(counter-1)]
          counter<-counter-1
        }
        
      }
    }
    monthly_data<-ts(clean_columns_data[,-1],start=c(as.numeric(substr(min(clean_columns_data$observation_date),1,4)),1),frequency=12)
    quarterly_data<-as.data.frame(aggregate(monthly_data,nfrequency = 4)/3)
    start_quarter<-paste0(ALFRED_mnemonic,paste0(substr(as.character(as.Date(as.yearqtr(as.Date(start_quarter_backup))+0.25)),1,4),"m",substr(as.character(as.Date(as.yearqtr(as.Date(start_quarter_backup))+0.25)),6,7),sep=""),sep="")
    start_quarter_position<-which(colnames(quarterly_data)==start_quarter) #starting vintage of estimation is 2008q4
    
    first_prints<-quarterly_data[,start_quarter_position:ncol(quarterly_data)]

    first_prints$observation_date<-as.Date(ts(1:nrow(quarterly_data), start = c(as.numeric(substr(min(clean_columns_data$observation_date),1,4)), as.numeric(substr(min(clean_columns_data$observation_date),6,7))), frequency = 4))
    first_prints$observation_date<-as.yearqtr(as.Date(first_prints$observation_date))
    
    assign(paste0(name_i,"_first_prints",sep=""),first_prints)
    
  }else{
    quarter_dates<-gdp_data$observation_date
    start_quarter<-"2005-01-01"
    start_quarter_backup<-start_quarter
    start_quarter_position<-which(quarter_dates==start_quarter) #starting vintage of estimation is 2008q4
    colnames(gdp_data)<-substr(colnames(gdp_data),1,nchar(colnames(gdp_data))-2)
    colnames(gdp_data)<-paste0(substr(colnames(gdp_data),1,nchar(colnames(gdp_data))-2),"m",substr(colnames(gdp_data),nchar(colnames(gdp_data))-1,nchar(colnames(gdp_data))))
    colnames(gdp_data)[which(colnames(gdp_data)=="observation_mda")]<-"observation_date"
    colnames(gdp_data)[which(!is.na(gdp_data[which(gdp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]] #the 2 in the bracket skips the observatin date
    
    gdp_data_first_prints_i<-as.data.frame(gdp_data[,colnames(gdp_data)[which(!is.na(gdp_data[which(gdp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]])
    colnames(gdp_data_first_prints_i)[1]<-colnames(gdp_data)[which(!is.na(gdp_data[which(gdp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]
    gdp_data_first_prints<-as.data.frame(gdp_data$observation_date)
    colnames(gdp_data_first_prints)[1]<-"observation_date"
    gdp_data_first_prints<-cbind(gdp_data_first_prints,gdp_data_first_prints_i)
    
    for(i in (start_quarter_position+1):nrow(gdp_data)){
      gdp_data_first_prints_i<-as.data.frame(gdp_data[,colnames(gdp_data)[which(!is.na(gdp_data[which(gdp_data$observation_date==quarter_dates[i]),]))[2]]])
      colnames(gdp_data_first_prints_i)[1]<-colnames(gdp_data)[which(!is.na(gdp_data[which(gdp_data$observation_date==quarter_dates[i]),]))[2]]
      gdp_data_first_prints<-cbind(gdp_data_first_prints,gdp_data_first_prints_i)
    }
    gdp_data_first_prints$observation_date<-as.yearqtr(as.Date(gdp_data_first_prints$observation_date))
  }
}

###########################################################################################################

quarter_dates<-pincome_data$observation_date
start_quarter<-"2005-01-01"
start_quarter_backup<-start_quarter
start_quarter_position<-which(quarter_dates==start_quarter) #starting vintage of estimation is 2008q4
colnames(pincome_data)<-substr(colnames(pincome_data),1,nchar(colnames(pincome_data))-2)
colnames(pincome_data)<-paste0(substr(colnames(pincome_data),1,nchar(colnames(pincome_data))-2),"m",substr(colnames(pincome_data),nchar(colnames(pincome_data))-1,nchar(colnames(pincome_data))))
colnames(pincome_data)[which(colnames(pincome_data)=="observation_mda")]<-"observation_date"
colnames(pincome_data)[which(!is.na(pincome_data[which(pincome_data$observation_date==quarter_dates[start_quarter_position]),]))[2]] #the 2 in the bracket skips the observatin date

pincome_data_first_prints_i<-as.data.frame(pincome_data[,colnames(pincome_data)[which(!is.na(pincome_data[which(pincome_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]])
colnames(pincome_data_first_prints_i)[1]<-colnames(pincome_data)[which(!is.na(pincome_data[which(pincome_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]
pincome_data_first_prints<-as.data.frame(pincome_data$observation_date)
colnames(pincome_data_first_prints)[1]<-"observation_date"
pincome_data_first_prints<-cbind(pincome_data_first_prints,pincome_data_first_prints_i)

for(i in (start_quarter_position+1):nrow(pincome_data)){
  pincome_data_first_prints_i<-as.data.frame(pincome_data[,colnames(pincome_data)[which(!is.na(pincome_data[which(pincome_data$observation_date==quarter_dates[i]),]))[2]]])
  colnames(pincome_data_first_prints_i)[1]<-colnames(pincome_data)[which(!is.na(pincome_data[which(pincome_data$observation_date==quarter_dates[i]),]))[2]]
  pincome_data_first_prints<-cbind(pincome_data_first_prints,pincome_data_first_prints_i)
}
pincome_data_first_prints$observation_date<-as.yearqtr(as.Date(pincome_data_first_prints$observation_date))

###########################################################################################################

quarter_dates<-pcomp_data$observation_date
start_quarter<-"2005-01-01"
start_quarter_backup<-start_quarter
start_quarter_position<-which(quarter_dates==start_quarter) #starting vintage of estimation is 2008q4
colnames(pcomp_data)<-substr(colnames(pcomp_data),1,nchar(colnames(pcomp_data))-2)
colnames(pcomp_data)<-paste0(substr(colnames(pcomp_data),1,nchar(colnames(pcomp_data))-2),"m",substr(colnames(pcomp_data),nchar(colnames(pcomp_data))-1,nchar(colnames(pcomp_data))))
colnames(pcomp_data)[which(colnames(pcomp_data)=="observation_mda")]<-"observation_date"
colnames(pcomp_data)[which(!is.na(pcomp_data[which(pcomp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]] #the 2 in the bracket skips the observatin date

pcomp_data_first_prints_i<-as.data.frame(pcomp_data[,colnames(pcomp_data)[which(!is.na(pcomp_data[which(pcomp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]])
colnames(pcomp_data_first_prints_i)[1]<-colnames(pcomp_data)[which(!is.na(pcomp_data[which(pcomp_data$observation_date==quarter_dates[start_quarter_position]),]))[2]]
pcomp_data_first_prints<-as.data.frame(pcomp_data$observation_date)
colnames(pcomp_data_first_prints)[1]<-"observation_date"
pcomp_data_first_prints<-cbind(pcomp_data_first_prints,pcomp_data_first_prints_i)

for(i in (start_quarter_position+1):nrow(pcomp_data)){
  pcomp_data_first_prints_i<-as.data.frame(pcomp_data[,colnames(pcomp_data)[which(!is.na(pcomp_data[which(pcomp_data$observation_date==quarter_dates[i]),]))[2]]])
  colnames(pcomp_data_first_prints_i)[1]<-colnames(pcomp_data)[which(!is.na(pcomp_data[which(pcomp_data$observation_date==quarter_dates[i]),]))[2]]
  pcomp_data_first_prints<-cbind(pcomp_data_first_prints,pcomp_data_first_prints_i)
}
pcomp_data_first_prints$observation_date<-as.yearqtr(as.Date(pcomp_data_first_prints$observation_date))

# merge data
minimum_dates_set<-c(min(gdp_data$observation_date),min(cpi_data$observation_date),
                     min(flbr_data$observation_date),min(fedfunds_data$observation_date),
                     min(houst_data$observation_date),min(payrolls_data$observation_date),
                     min(fip_data$observation_date),min(fcu_data$observation_date),
                     min(ftfx_data$observation_date),min(frgt3M_data$observation_date),
                     min(frgt6M_data$observation_date),min(frgt1y_data$observation_date),
                     min(frgt3y_data$observation_date),min(frgt5y_data$observation_date),
                     min(frgt10y_data$observation_date),min(m2_data$observation_date),
                     min(pincome_data$observation_date),min(pincome_data$observation_date),
                     min(pcomp_data$observation_date),min(baamoodys_data$observation_date),
                     min(psavert_data$observation_date),min(ppiaco_data$observation_date))

max_min_date<-max(minimum_dates_set)

if(length(unique(minimum_dates_set))>1){
  gdp_data_first_prints<-gdp_data_first_prints[which(as.Date(gdp_data_first_prints$observation_date)==max_min_date):nrow(gdp_data_first_prints),]
  flbr_data_first_prints<-flbr_data_first_prints[which(as.Date(gdp_data_first_prints$observation_date)==max_min_date):nrow(flbr_data_first_prints),]
  cpi_data_first_prints<-cpi_data_first_prints[which(as.Date(cpi_data_first_prints$observation_date)==max_min_date):nrow(cpi_data_first_prints),]
  fedfunds_data_first_prints<-fedfunds_data_first_prints[which(as.Date(fedfunds_data_first_prints$observation_date)==max_min_date):nrow(fedfunds_data_first_prints),]
  payrolls_data_first_prints<-payrolls_data_first_prints[which(as.Date(payrolls_data_first_prints$observation_date)==max_min_date):nrow(payrolls_data_first_prints),]
  houst_data_first_prints<-houst_data_first_prints[which(as.Date(houst_data_first_prints$observation_date)==max_min_date):nrow(houst_data_first_prints),]
  fip_data_first_prints<-fip_data_first_prints[which(as.Date(fip_data_first_prints$observation_date)==max_min_date):nrow(fip_data_first_prints),]
  fcu_data_first_prints<-fcu_data_first_prints[which(as.Date(fcu_data_first_prints$observation_date)==max_min_date):nrow(fcu_data_first_prints),]
  ftfx_data_first_prints<-ftfx_data_first_prints[which(as.Date(ftfx_data_first_prints$observation_date)==max_min_date):nrow(ftfx_data_first_prints),]
  frgt3M_data_first_prints<-frgt3M_data_first_prints[which(as.Date(frgt3M_data_first_prints$observation_date)==max_min_date):nrow(frgt3M_data_first_prints),]
  frgt6M_data_first_prints<-frgt6M_data_first_prints[which(as.Date(frgt6M_data_first_prints$observation_date)==max_min_date):nrow(frgt6M_data_first_prints),]
  frgt1y_data_first_prints<-frgt1y_data_first_prints[which(as.Date(frgt1y_data_first_prints$observation_date)==max_min_date):nrow(frgt1y_data_first_prints),]
  frgt3y_data_first_prints<-frgt3y_data_first_prints[which(as.Date(frgt3y_data_first_prints$observation_date)==max_min_date):nrow(frgt3y_data_first_prints),]
  frgt5y_data_first_prints<-frgt5y_data_first_prints[which(as.Date(frgt5y_data_first_prints$observation_date)==max_min_date):nrow(frgt5y_data_first_prints),]
  frgt10y_data_first_prints<-frgt10y_data_first_prints[which(as.Date(frgt10y_data_first_prints$observation_date)==max_min_date):nrow(frgt10y_data_first_prints),]
  m2_data_first_prints<-m2_data_first_prints[which(as.Date(m2_data_first_prints$observation_date)==max_min_date):nrow(m2_data_first_prints),]
  pincome_data_first_prints<-pincome_data_first_prints[which(as.Date(pincome_data_first_prints$observation_date)==max_min_date):nrow(pincome_data_first_prints),]
  pcomp_data_first_prints<-pcomp_data_first_prints[which(as.Date(pcomp_data_first_prints$observation_date)==max_min_date):nrow(pcomp_data_first_prints),]
  baamoodys_data_first_prints<-baamoodys_data_first_prints[which(as.Date(baamoodys_data_first_prints$observation_date)==max_min_date):nrow(baamoodys_data_first_prints),]
  ppiaco_data_first_prints<-ppiaco_data_first_prints[which(as.Date(ppiaco_data_first_prints$observation_date)==max_min_date):nrow(ppiaco_data_first_prints),]
}

############################################################################
# correct data irregularities
############################################################################

find_irregularities<-grep("11",substr(names(gdp_data_first_prints),nchar(names(gdp_data_first_prints))-1,nchar(names(gdp_data_first_prints))))
names(gdp_data_first_prints)[find_irregularities]<-"GDP_2013m10"
find_irregularities<-grep("02",substr(names(gdp_data_first_prints),nchar(names(gdp_data_first_prints))-1,nchar(names(gdp_data_first_prints))))
names(gdp_data_first_prints)[find_irregularities]<-"GDP_2019m01"


find_irregularities<-grep("11",substr(names(pcomp_data_first_prints),nchar(names(pcomp_data_first_prints))-1,nchar(names(pcomp_data_first_prints))))
names(pcomp_data_first_prints)[find_irregularities]<-"PCECC96_2013m10"
find_irregularities<-grep("02",substr(names(pcomp_data_first_prints),nchar(names(pcomp_data_first_prints))-1,nchar(names(pcomp_data_first_prints))))
names(pcomp_data_first_prints)[find_irregularities]<-"PCECC96_2019m01"

houst_data_first_prints$HOUST_2013m10<-houst_data_first_prints$HOUST_2014m01
houst_data_first_prints$HOUST_2019m01<-houst_data_first_prints$HOUST_2019m04

pincome_data_first_prints$PINCOME_2009m07<-pincome_data_first_prints$PINCOME_2009m10
find_irregularities<-grep("11",substr(names(pincome_data_first_prints),nchar(names(pincome_data_first_prints))-1,nchar(names(pincome_data_first_prints))))
names(pincome_data_first_prints)[find_irregularities]<-"PINCOME_2013m10"
find_irregularities<-grep("02",substr(names(pincome_data_first_prints),nchar(names(pincome_data_first_prints))-1,nchar(names(pincome_data_first_prints))))
names(pincome_data_first_prints)[find_irregularities]<-"PINCOME_2019m01"


data<-psavert_data
colnames_backup<-colnames(data)
vintage_dates<-as.numeric(substr(colnames_backup,nchar(colnames_backup)-7,nchar(colnames_backup)))

colnames(data)<-substr(str_replace_all(colnames(data),paste0(name_i,".",sep=""),""),1,nchar(colnames(data))-2)
colnames(data)<-substr(colnames(data),1,nchar(colnames(data)))
colnames(data)<-paste0(substr(colnames(data),1,nchar(colnames(data))-2),"m",substr(colnames(data),nchar(colnames(data))-1,nchar(colnames(data))),sep="")
colnames(data)[1]<-"observation_date"
ALFRED_mnemonic<-substr(colnames(data[-1][1]),1,nchar(colnames(data[-1][1]))-7)

clean_columns_data<-as.data.frame(data[,which(colnames(data)=="observation_date")])
colnames(clean_columns_data)[1]<-"observation_date"
#ncol(fip_data)
counter<-1
for(j in 2:ncol(data)){
  if(str_detect(colnames(data)[j],"m01|m04|m07|m10")==TRUE){
    clean_columns_data_i<-as.data.frame(data[,j])
    colnames(clean_columns_data_i) <- colnames(data[j])
    clean_columns_data<-cbind(clean_columns_data,clean_columns_data_i)
    counter<-counter+1
    if(colnames(clean_columns_data)[counter]==colnames(clean_columns_data)[counter-1] & vintage_dates[j]>vintage_dates[j-1]){
      clean_columns_data<-clean_columns_data[,-(counter-1)]
      counter<-counter-1
    }
    
  }
}

clean_columns_data$PSAVERT_2005m07<-psavert_data$PSAVERT_20050802
clean_columns_data$PSAVERT_2006m04<-psavert_data$PSAVERT_20060501
clean_columns_data$PSAVERT_2006m07<-psavert_data$PSAVERT_20060801
clean_columns_data$PSAVERT_2007m01<-psavert_data$PSAVERT_20061222
clean_columns_data$PSAVERT_2007m10<-psavert_data$PSAVERT_20071101
clean_columns_data$PSAVERT_2008m04<-psavert_data$PSAVERT_20080501
clean_columns_data$PSAVERT_2008m07<-psavert_data$PSAVERT_20080804
clean_columns_data$PSAVERT_2009m01<-psavert_data$PSAVERT_20090202
clean_columns_data$PSAVERT_2009m07<-psavert_data$PSAVERT_20080804
clean_columns_data$PSAVERT_2010m01<-psavert_data$PSAVERT_20091223
clean_columns_data$PSAVERT_2010m04<-psavert_data$PSAVERT_20100503
clean_columns_data$PSAVERT_2010m07<-psavert_data$PSAVERT_20100803
clean_columns_data$PSAVERT_2011m07<-psavert_data$PSAVERT_20110802
clean_columns_data$PSAVERT_2013m07<-psavert_data$PSAVERT_20130802
clean_columns_data$PSAVERT_2013m10<-psavert_data$PSAVERT_20131223
clean_columns_data$PSAVERT_2014m04<-psavert_data$PSAVERT_20140501
clean_columns_data$PSAVERT_2014m07<-psavert_data$PSAVERT_20140801
clean_columns_data$PSAVERT_2015m01<-psavert_data$PSAVERT_20150202
clean_columns_data$PSAVERT_2015m07<-psavert_data$PSAVERT_20150803
clean_columns_data$PSAVERT_2016m01<-psavert_data$PSAVERT_20151125
clean_columns_data$PSAVERT_2016m07<-psavert_data$PSAVERT_20160802
clean_columns_data$PSAVERT_2017m04<-psavert_data$PSAVERT_20170501
clean_columns_data$PSAVERT_2017m07<-psavert_data$PSAVERT_20170801
clean_columns_data$PSAVERT_2019m01<-psavert_data$PSAVERT_20190301

monthly_data<-ts(clean_columns_data[,-1],start=c(as.numeric(substr(min(clean_columns_data$observation_date),1,4)),1),frequency=12)
quarterly_data<-as.data.frame(aggregate(monthly_data,nfrequency = 4)/3)
start_quarter<-paste0(ALFRED_mnemonic,paste0(substr(as.character(as.Date(as.yearqtr(as.Date(start_quarter_backup))+0.25)),1,4),"m",substr(as.character(as.Date(as.yearqtr(as.Date(start_quarter_backup))+0.25)),6,7),sep=""),sep="")
start_quarter_position<-which(colnames(quarterly_data)==start_quarter) #starting vintage of estimation is 2008q4

first_prints<-quarterly_data[,start_quarter_position:ncol(quarterly_data)]
first_prints$observation_date<-as.Date(ts(1:nrow(quarterly_data), start = c(as.numeric(substr(min(clean_columns_data$observation_date),1,4)), as.numeric(substr(min(clean_columns_data$observation_date),6,7))), frequency = 4))
first_prints$observation_date<-as.yearqtr(as.Date(first_prints$observation_date))

assign(paste0("psavert_data","_first_prints",sep=""),first_prints)
psavert_data_first_prints<-psavert_data_first_prints[which(as.Date(psavert_data_first_prints$observation_date)==max_min_date):nrow(psavert_data_first_prints),]

baamoodys_data_first_prints$BAA_2017m01<-baamoodys_data_first_prints$BAA_2017m04

full_data_first_prints<-merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(
 gdp_data_first_prints
,flbr_data_first_prints,by='observation_date')
,cpi_data_first_prints,by='observation_date')
,fedfunds_data_first_prints,by='observation_date')
,payrolls_data_first_prints,by='observation_date')
,houst_data_first_prints,by='observation_date')
,fip_data_first_prints,by='observation_date')
,fcu_data_first_prints,by='observation_date')
,ftfx_data_first_prints,by='observation_date')
,m2_data_first_prints,by='observation_date')
,pincome_data_first_prints,by='observation_date')
,pcomp_data_first_prints,by='observation_date')
,ppiaco_data_first_prints,by='observation_date')
,frgt3M_data_first_prints,by='observation_date')
,frgt6M_data_first_prints,by='observation_date')
,frgt1y_data_first_prints,by='observation_date')
,frgt3y_data_first_prints,by='observation_date')
,frgt5y_data_first_prints,by='observation_date')
,frgt10y_data_first_prints,by='observation_date')
,baamoodys_data_first_prints,by='observation_date')


year = 2007:2010
min_y = min(year)
max_y = max(year)
month = c(1,4,7,10)
all_vintages<-list()
for(year_i in year){
  for(month_i in month){
    if(month_i == 10){
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m",as.character(month_i)))]))
    }else{
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m0",as.character(month_i)))]))
    }
  }
}

require(xlsx)
for(i_ in 1:length(all_vintages)){
  vintage_i<-as.data.frame(all_vintages[[i_]])
  new_colnames<-as.data.table(str_split(colnames(vintage_i),"_"))[1,]
  colnames(vintage_i)<-new_colnames[1,]
  #xlsx::write.xlsx(, file = names(all_vintages[i_]), sheetName = names(all_vintages[i_]), append = FALSE, row.names = FALSE)
  xlsx::write.xlsx(vintage_i, file = paste0(path_link,"/vintages_",paste0(min_y,"_",max_y,".xlsx",sep=""),sep=""), sheetName = names(all_vintages[i_]), append = TRUE, row.names = FALSE)
  Sys.sleep(1)
}

year = 2011:2014 # Choose the forecast evaluation vintages. 2009q1 vintage means data available till 2008q4
min_y = min(year)
max_y = max(year)
month = c(1,4,7,10)
all_vintages<-list()
for(year_i in year){
  for(month_i in month){
    if(month_i == 10){
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m",as.character(month_i)))]))
    }else{
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m0",as.character(month_i)))]))
    }
  }
}

require(xlsx)
for(i_ in 1:length(all_vintages)){
  vintage_i<-as.data.frame(all_vintages[[i_]])
  new_colnames<-as.data.table(str_split(colnames(vintage_i),"_"))[1,]
  colnames(vintage_i)<-new_colnames[1,]
  #xlsx::write.xlsx(, file = names(all_vintages[i_]), sheetName = names(all_vintages[i_]), append = FALSE, row.names = FALSE)
  xlsx::write.xlsx(vintage_i, file = paste0(path_link,"/vintages_",paste0(min_y,"_",max_y,".xlsx",sep=""),sep=""), sheetName = names(all_vintages[i_]), append = TRUE, row.names = FALSE)
  Sys.sleep(1)
}

year = 2015:2018 # Choose the forecast evaluation vintages. 2009q1 vintage means data available till 2008q4
min_y = min(year)
max_y = max(year)
month = c(1,4,7,10)
all_vintages<-list()
for(year_i in year){
  for(month_i in month){
    if(month_i == 10){
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m",as.character(month_i)))]))
    }else{
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m0",as.character(month_i)))]))
    }
  }
}

require(xlsx)
for(i_ in 1:length(all_vintages)){
  vintage_i<-as.data.frame(all_vintages[[i_]])
  new_colnames<-as.data.table(str_split(colnames(vintage_i),"_"))[1,]
  colnames(vintage_i)<-new_colnames[1,]
  #xlsx::write.xlsx(, file = names(all_vintages[i_]), sheetName = names(all_vintages[i_]), append = FALSE, row.names = FALSE)
  xlsx::write.xlsx(vintage_i, file = paste0(path_link,"/vintages_",paste0(min_y,"_",max_y,".xlsx",sep=""),sep=""), sheetName = names(all_vintages[i_]), append = TRUE, row.names = FALSE)
  Sys.sleep(1)
}


year = 2019:2022 # Choose the forecast evaluation vintages. 2009q1 vintage means data available till 2008q4
min_y = min(year)
max_y = max(year)
month = c(1,4,7,10)
all_vintages<-list()
for(year_i in year){
  for(month_i in month){
    if(month_i == 10){
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m",as.character(month_i)))]))
    }else{
      all_vintages[[paste0(as.character(year_i),"m",as.character(month_i),sep="")]]<-as.data.frame(na.omit(full_data_first_prints[str_detect(colnames(full_data_first_prints),paste0("obs|",as.character(year_i),"m0",as.character(month_i)))]))
    }
  }
}

require(xlsx)
for(i_ in 1:length(all_vintages)){
  vintage_i<-as.data.frame(all_vintages[[i_]])
  new_colnames<-as.data.table(str_split(colnames(vintage_i),"_"))[1,]
  colnames(vintage_i)<-new_colnames[1,]
  #xlsx::write.xlsx(, file = names(all_vintages[i_]), sheetName = names(all_vintages[i_]), append = FALSE, row.names = FALSE)
  xlsx::write.xlsx(vintage_i, file = paste0(path_link,"/vintages_",paste0(min_y,"_",max_y,".xlsx",sep=""),sep=""), sheetName = names(all_vintages[i_]), append = TRUE, row.names = FALSE)
  Sys.sleep(1)
}
