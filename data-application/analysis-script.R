# This script is developed by Lingxuan Shao #

######wash data
library("ggplot2")
library("rgl") 
library("matrixcalc")
library("pracma")

# load data that have been prepared by other scripts; see README in https://github.com/linulysses/iRRDA-sparse
load('hippo.RData')

dti <- dti * 1000 # amplify each DTI by 1000 for numerical stability

d=3;

# remove some outliers
checkdti=array(0,nrow(metainfo));
for(l in 1:nrow(metainfo)){
    checkdti[l]=sum(abs(dti[,,l]));
}

co = setdiff(1:nrow(metainfo),which(checkdti>5))
newdti=dti[,,co];
newmetainfo=metainfo[co,]

iname=names(table(newmetainfo$subject_id));
n=length(iname);
gname=names(table(newmetainfo$group));
#make Y
Y=array(0,c(n,20,d,d));
t=array(-10,c(n,20));
gin=array(0,n);
m=array(0,n);
for(l in 1:nrow(newmetainfo)){
    thisid=newmetainfo$subject_id[l];
    thisi=0;
    for(i in 1:n){
        if(thisid==iname[i]){
            thisi=i;
            break;
        }
    }
    if(thisi==0){
        print(l);
        break;
    }
    m[thisi]=m[thisi]+1;
    Y[thisi,m[thisi],,]=newdti[,,l];
    t[thisi,m[thisi]]=newmetainfo$age[l];
    gin[thisi]=as.numeric(newmetainfo$group[l]=="CN");
}

#record
orin=n;
oriY=Y;
orit=t;
orim=m;


#CN/AZ;
CNAZ=0; #1 for CN and 0 for AZ
source('analysis.R')
save(file='AZ.RData',meancurve,matrixeigen,hmu,hc,mu_leave_one_out_error,hmulist,
     c_leave_one_out_error,hclist,pcaresult_values,tanchor)


CNAZ=1;
source('analysis.R')
save(file='CN.RData',meancurve,matrixeigen,hmu,hc,mu_leave_one_out_error,hmulist,
     c_leave_one_out_error,hclist,pcaresult_values,tanchor)
