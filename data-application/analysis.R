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
CNAZ=0;#1 for CN and 0 for AZ
n=length(which(gin==CNAZ));
Y=oriY[which(gin==CNAZ),,,];
t=orit[which(gin==CNAZ),];
m=orim[which(gin==CNAZ)];
maxm=max(m);


###use log-cholesky norm
cut=function(A,B){
  newA=A;
  na=length(A);
  nb=length(B);
  nc=na-nb;
  C=array(0,nc);
  ib=1;
  for(ia in 1:na){
    if(A[ia]==B[ib]){
      newA[ia]=0;
      ib=ib+1;
    }
    if(ib>nb){
      break;
    }
  }
  ic=1;
  for(ia in 1:na){
    if(newA[ia]>0){
      C[ic]=newA[ia];
      ic=ic+1;
    }
    if(ic>nc){
      break;
    }
  }
  return(C);
}
#trans, make matrix to vector
trans=function(A){
  RA=array(0,6);
  k=0;
  for(i in 1:d){
    for(j in 1:i){
      k=k+1;
      RA[k]=A[i,j];
      if(i==j){
        RA[k]=log(RA[k]);
      }
    }
  }
  return(RA);
}
newy=array(0,c(n,maxm,6));
for(i in 1:n){
  for(j in 1:m[i]){
    newy[i,j,]=trans(t(chol(Y[i,j,,])));
  }
}


#Kernel function K=C(1-|t|^3)^3
Ker=function(t,h){
  if((t-h)>0){
    return(0);
  }
  else if((t+h)<0){
    return(0);
  }
  else{
    return(((1-(abs(t/h))**3)**3)*70/(h*81));
  }
}


#hatmu is a function
hatmu=function(s,inputnewn,inputnewm,inputt,inputy,hmu){
  hatu1=0;
  hatu2=0;
  for(i in 1:inputnewn){
    for(j in 1:inputnewm[i]){
      hatu1=hatu1+1000*Ker(inputt[i,j]-s,hmu)*(inputt[i,j]-s);
      hatu2=hatu2+1000*Ker(inputt[i,j]-s,hmu)*(inputt[i,j]-s)*(inputt[i,j]-s);
    }
  }
  up=array(0,6);
  down=0;
  for(i in 1:inputnewn){
    for(j in 1:inputnewm[i]){
      if(Ker(inputt[i,j]-s,hmu)!=0){
        up=up+1000*Ker(inputt[i,j]-s,hmu)*(hatu2-hatu1*(inputt[i,j]-s))*inputy[i,j,];
        down=down+1000*Ker(inputt[i,j]-s,hmu)*(hatu2-hatu1*(inputt[i,j]-s));
      }
    }
  }
  if(down==0){
    print("miss")
    print(s);
    return(array(-100,6));
  }else{
    return((1/down)*up);
  }
}

#hatc is a function
hatc=function(s1,s2,n,m,t,logY,hc){
  S00=0;
  S01=0;
  S10=0;
  S11=0;
  S20=0;
  S02=0;
  R00=array(0,c(6,6));
  R01=array(0,c(6,6));
  R10=array(0,c(6,6));
  for(i in 1:n){
    for(j1 in 1:m[i]){
      for(j2 in 1:m[i]){
        if((abs(t[i,j1]-s1)<hc)&&(abs(t[i,j2]-s2)<hc)&&(j1!=j2)){
          translogYlogY=logY[i,j1,]%*%t(logY[i,j2,]);
          S00=S00+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc);
          S10=S10+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc);
          S01=S01+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc);
          S11=S11+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*((t[i,j2]-s2)/hc);
          S20=S20+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*((t[i,j1]-s1)/hc);
          S02=S02+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc)*((t[i,j2]-s2)/hc);
          R00=R00+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*translogYlogY;
          R10=R10+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j1]-s1)/hc)*translogYlogY;
          R01=R01+Ker(t[i,j1]-s1,hc)*Ker(t[i,j2]-s2,hc)*((t[i,j2]-s2)/hc)*translogYlogY;
        }
      }
    }
  }
  S20021111=S20*S02-S11*S11;
  S10020111=S10*S02-S01*S11;
  S10110120=S10*S11-S01*S20;
  if((S20021111*S00-S10020111*S10+S10110120*S01)==0){
    return(array(-100,c(6,6)));
    print("wrong");
  }
  else{
    return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
  }
}

#################estimation

##leave-one-out to select hmu
hmulong=10;
hmulist=exp(c(0:(hmulong-1))*log(3)/(30-1)+log(3));
mu_leave_one_out_error=array(0,hmulong);
for(hmuit in 1:hmulong){
  print(hmuit);
  hmu=hmulist[hmuit];
  for(testindex in 1:n){
    trainindex=cut(c(1:n),c(testindex));
    newn=n-1;
    newm=m[trainindex];
    newt=t[trainindex,];
    newnewy=newy[trainindex,,];
    temp=0;
    for(j in 1:m[testindex]){
      temp=temp+sum((hatmu(t[testindex,j],newn,newm,newt,newnewy,hmu)-newy[testindex,j,])**2);
    }
    temp=temp/m[testindex];
    mu_leave_one_out_error[hmuit]=mu_leave_one_out_error[hmuit]+temp;
  }
}

#estimate hatc
hmu=4.218843;
logY=array(0,c(n,maxm,6));
for(i in 1:n){
  for(j in 1:m[i]){
    logY[i,j,]=newy[i,j,]-hatmu(t[i,j],n,m,t,newy,hmu);
  }
}

##leave-one-out to select hc
hclong=13;
hclist=exp(c(0:(hclong-1))*log(3)/(13-1)+log(3));
c_leave_one_out_error=array(0,hclong);
for(hcit in 1:hclong){
  print(hcit);
  hc=hclist[hcit];
  for(testindex in 1:n){
    trainindex=cut(c(1:n),c(testindex));
    newn=n-1;
    newm=m[trainindex];
    newt=t[trainindex,];
    newlogY=logY[trainindex,,];
    temp=0;
    for(j1 in 1:m[testindex]){
      for(j2 in 1:m[testindex]){
        if(j1!=j2){
          temp=temp+sum((hatc(t[testindex,j1],t[testindex,j2],newn,newm,newt,newlogY,hc)-logY[testindex,j1,]%*%t(logY[testindex,j2,]))**2);
        }
      }
    }
    temp=temp/(m[testindex]*(m[testindex]-1));
    c_leave_one_out_error[hcit]=c_leave_one_out_error[hcit]+temp;
  }
}


####estimate C in a grid

inv=function(RA){
  kk=0;
  A=array(0,c(3,3));
  for(i in 1:3){
    for(j in 1:i){
      kk=kk+1;
      if(i==j){
        A[i,j]=exp(RA[kk]);
      }
      if(i!=j){
        A[i,j]=RA[kk];
      }
    }
  }
  AA=A%*%t(A);
  return(AA);
}


hc=4.726621;
l=99;
tbegin=55.2082;
tend=93.5500;
tanchor=c(0:(l-1))*(tend-tbegin)/(l-1)+tbegin;
meancurve=array(0,c(l,d,d));
for(j in 1:l){
  meancurve[j,,]=inv(hatmu(tanchor[j],n,m,t,newy,hmu));
}
matrixhatc=array(0,c(l,l,6,6));
for(j1 in 1:l){
  print(j1);
  for(j2 in 1:l){
    matrixhatc[j1,j2,,]=hatc(tanchor[j1],tanchor[j2],n,m,t,logY,hc);
  }
}
#pca
bigmatrix=array(0,c(6*l,6*l));
for(j1 in 1:l){
  for(j2 in 1:l){
    for(j3 in 1:6){
      for(j4 in 1:6){
        bigmatrix[(j1-1)*6+j3,(j2-1)*6+j4]=matrixhatc[j1,j2,j3,j4];
      }
    }
  }
}
for(j1 in 1:(6*l)){
  for(j2 in 1:(6*l)){
    temp=0.5*(bigmatrix[j1,j2]+bigmatrix[j2,j1]);
    bigmatrix[j1,j2]=temp;
    bigmatrix[j2,j1]=temp;
  }
}
pcaresult=eigen(bigmatrix);
pcaresult_values=pcaresult$values/l;
ggplot(data.frame(cbind(c(1:30),pcaresult$values[1:30])))+geom_path(aes(x=X1,y=X2),size=1)+
  xlab(label="")+ylab(label="eigenvalue");
#eigenfunction
eigen=array(0,c(l*6,l,6));
for(k in 1:(l*6)){
  for(j1 in 1:l){
    for(j2 in 1:6){
      eigen[k,j1,j2]=pcaresult$vectors[(j1-1)*6+j2,k]*sqrt(l/(tend-tbegin));
    }
  }
}

matrixeigen=array(0,c(l*6,l,3,3));
for(k in 1:(l*6)){
  print(k);
  for(j in 1:l){
    matrixeigen[k,j,,]=inv(eigen[k,j,]+hatmu(tanchor[j],n,m,t,newy,hmu));
  }
}

save(file='cc_AZ.RData',meancurve,matrixeigen,hmu,hc,mu_leave_one_out_error,hmulist,
     c_leave_one_out_error,hclist,pcaresult_values,tanchor)
