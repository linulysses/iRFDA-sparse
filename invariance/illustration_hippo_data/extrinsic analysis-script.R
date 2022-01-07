# This script is developed by Lingxuan Shao #

######wash data
library("ggplot2")
library("rgl") 
library("matrixcalc")
library("pracma")

load('hippo.RData')

d=3;
iname=names(table(metainfo[,4]));
n=length(iname);
gname=names(table(metainfo[,5]));
#make Y
Y=array(0,c(n,20,d,d));
t=array(-10,c(n,20));
gin=array(0,n);
m=array(0,n);
for(l in 1:977){
    thisid=metainfo[l,4];
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
    Y[thisi,m[thisi],,]=dti[,,l];
    t[thisi,m[thisi]]=metainfo[l,7];
    gin[thisi]=as.numeric(metainfo[l,5]=="CN");
}

#record
orin=n;rm(n);
oriY=Y;rm(Y);
orit=t;rm(t);
orim=m;rm(m);

#consider AZ;
CNAZ=0; #1 for CN and 0 for AZ

#no rotation for the frame
rotation=0;
source('extrinsic analysis.R')
save(file='extrinsic_AZ_rotation0.RData',meancurve,matrixeigen,hmu,hc,pcaresult_values,tanchor)

#some rotation for the frame
rotation=1;
source('extrinsic analysis.R')
save(file='extrinsic_AZ_rotation1.RData',meancurve,matrixeigen,hmu,hc,pcaresult_values,tanchor)
