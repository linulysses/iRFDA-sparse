# Simulation: Affine-Invariant
# Developed by Lingxuan Shao

# l time points of the time interval are focused on in the computation
l=20;
# calculate the estimates on tanchor
tanchor=c(1:l)/l;

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

#tensor product of vectors A,B
trans=function(A,B){
  newc=array(0,c(2,2));
  newc[1,1]=A[1]*B[1];
  newc[1,2]=A[1]*B[2];
  newc[2,1]=A[2]*B[1];
  newc[2,2]=A[2]*B[2];
  return(newc);
}

#calculate \hat{\mu} on time s 
hatmu=function(s,inputnewn,inputnewm,inputt,inputy,hmu){
  hatu1=0;
  hatu2=0;
  for(i in 1:inputnewn){
    for(j in 1:inputnewm[i]){
      hatu1=hatu1+1000*Ker(inputt[i,j]-s,hmu)*(inputt[i,j]-s);
      hatu2=hatu2+1000*Ker(inputt[i,j]-s,hmu)*(inputt[i,j]-s)*(inputt[i,j]-s);
    }
  }
  up=array(0,2);
  down=0;
  for(i in 1:inputnewn){
    for(j in 1:inputnewm[i]){
      up=up+1000*Ker(inputt[i,j]-s,hmu)*(hatu2-hatu1*(inputt[i,j]-s))*inputy[i,j,];
      down=down+1000*Ker(inputt[i,j]-s,hmu)*(hatu2-hatu1*(inputt[i,j]-s));
    }
  }
  if(down==0){
    print("miss")
    print(s);
    return(array(0,2));
  }else{
    return((1/down)*up);
  }
}


# calculate \hat{C} on time (s1,s2)
hatc=function(s1,s2,n,m,t,logY,hc){
  S00=0;
  S01=0;
  S10=0;
  S11=0;
  S20=0;
  S02=0;
  R00=array(0,c(2,2));
  R01=array(0,c(2,2));
  R10=array(0,c(2,2));
  for(i in 1:n){
    for(j1 in 1:m[i]){
      for(j2 in 1:m[i]){
        if((abs(t[i,j1]-s1)<hc)&&(abs(t[i,j2]-s2)<hc)&&(j1!=j2)){
          translogYlogY=trans(logY[i,j1,],logY[i,j2,]);
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
    return(array(0,c(2,2)));
  }
  else{
    return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
  }
}

#real vales of C(s1,s2)
realc=function(s1,s2){
  newc=array(0,c(2,2));
  newc[1,1]=s1*s2/300;
  newc[2,2]=s1*s2/300;
  return(newc);
}

#real vales of C on the grid of tanchor*tanchor
realcdata=array(0,c(l,l,2,2));
for(k1 in 1:l){
  for(k2 in 1:l){
    realcdata[k1,k2,,]=realc(tanchor[k1],tanchor[k2]);
  }
}


################
#Mentor Carlo for itermax times
itermax=100;
#record the sup and l2 norms of difference of our proposed estimator
supnormlist=array(0,itermax);
l2normlist=array(0,itermax);
for(it in 1:itermax){
  #generate m,t
  print(it);
  m=rpois(n,mp)+2;
  mmax=max(m);
  N=sum(m);
  t=array(0,c(n,mmax));
  for(i in 1:n){
    t[i,]=runif(mmax,0,1);
  }
  #generate data Y 
  #Y is of the form Y=P^{T}diag{exp(t+tZ1),exp(t+tZ2)}P
  #M is isometric to PMP^{T}
  #For diagonal matrix, affine-invariant = log-Euclidean
  #transfer Y to R^{2}
  Z1=runif(n,-0.1,0.1);
  Z2=runif(n,-0.1,0.1);
  Y=array(0,c(n,mmax,2));
  for(i in 1:n){
    for(j in 1:m[i]){
      vare=runif(1,-0.025,0.025);
      Y[i,j,1]=t[i,j]+t[i,j]*Z1[i]+vare;
      Y[i,j,2]=t[i,j]+t[i,j]*Z2[i]+vare;
    }
  }
  #decompose the data into foldnum=2 parts
  #first part for estimation and second part for validation, and plus
  #second part for estimation and first part for validation
  foldnum=2;
  validationindex=array(0,c(foldnum,n/foldnum));
  estimaindex=array(0,c(foldnum,(n-n/foldnum)));
  for(fold in 1:foldnum){
    validationindex[fold,]=(fold-1)*n/foldnum+c(1:(n/foldnum));
    estimaindex[fold,]=c(1:n)[-validationindex[fold,]];
  }
  #select hmu with minimum valierror, hmulist are the candidates of hmu
  hmulong=8;
  hmulist=exp(c(0:(hmulong-1))*log(10)/(hmulong-1)+log(0.05));
  valierror=array(0,hmulong);#record the validation error for each hmu
  for(hmuit in 1:hmulong){
    hmu=hmulist[hmuit];
    for(fold in 1:foldnum){
      #use this fold for validation
      for(i in validationindex[fold,]){
        for(j in 1:m[i]){
          valierror[hmuit]=valierror[hmuit]+sum((hatmu(t[i,j],length(estimaindex[fold,]),m[estimaindex[fold,]],t[estimaindex[fold,],],Y[estimaindex[fold,],,],hmu)-Y[i,j,])**2)/m[i];
        }
      }
    }
  }
  hmu=hmulist[which.min(valierror)];
  #compute hat{mu}(t[i,j])
  hatmuofY=array(0,c(n,mmax,2));
  for(i in 1:n){
    for(j in 1:m[i]){
      hatmuofY[i,j,]=hatmu(t[i,j],n,m,t,Y,hmu);
    }
  }
  #compute raw-covariances Y(t[i,j])-hat{mu}(t[i,j])
  loghatmuY=array(0,c(n,mmax,2));
  for(i in 1:n){
    for(j in 1:m[i]){
      loghatmuY[i,j,]=Y[i,j,]-hatmuofY[i,j,];
    }
  }
  #select hc with minimum valierror, hclist are the candidates of hc
  hclong=8;
  hclist=exp(c(0:(hclong-1))*log(10)/(hclong-1)+log(0.05));
  valierror=array(0,hclong);#record the validation error for each hc
  for(hcit in 1:hclong){
    hc=hclist[hcit];
    for(fold in 1:foldnum){
      #use this fold for validation
      for(i in validationindex[fold,]){
        if(m[i]==1){next;}
        for(j1 in 1:m[i]){
        for(j2 in 1:m[i]){
        if(j1<j2){
          valierror[hcit]=valierror[hcit]+sum((hatc(t[i,j1],t[i,j2],length(estimaindex[fold,]),m[estimaindex[fold,]],t[estimaindex[fold,],],loghatmuY[estimaindex[fold,],,],hc)-trans(loghatmuY[i,j1,],loghatmuY[i,j2,]))**2)/(m[i]*(m[i]-1));
        }}}
      }
    }
  }
  hc=hclist[which.min(valierror)];
  #record the supnorm and l2 norm
  supnorm=0;
  l2norm=0;
  for(k1 in 1:l){
    for(k2 in 1:l){
      temp=hatc(tanchor[k1],tanchor[k2],n,m,t,loghatmuY,hc);
      supnorm=max(supnorm,sqrt(sum((realcdata[k1,k2,,]-temp)**2)));
      l2norm=l2norm+sum((realcdata[k1,k2,,]-temp)**2);
    }
  }
  supnormlist[it]=supnorm;
  l2normlist[it]=sqrt(l2norm/(l*l));
}


#the supnorm and l2 norm of C
realsupnorm=0;
reall2norm=0;
for(k1 in 1:l){
  for(k2 in 1:l){
    realsupnorm=max(realsupnorm,sqrt(sum(realcdata[k1,k2,,]**2)));
    reall2norm=reall2norm+sum(realcdata[k1,k2,,]**2);
  }
}
reall2norm=sqrt(reall2norm/(l*l));

#relative error
relasup=mean(supnormlist)/realsupnorm*100;
relal2=mean(l2normlist)/reall2norm*100;
relasupsd=sd(supnormlist)/realsupnorm*100;
relal2sd=sd(l2normlist)/reall2norm*100;
print("Answer is:")
print("relasup");
print(relasup);
print("relal2");
print(relal2);
print("relasupsd");
print(relasupsd);
print("relal2sd");
print(relal2sd);

#output
text=data.frame(c(n,mp,relasup,relal2,relasupsd,relal2sd));
name=paste("affine_n",n,"_mp",mp,".txt",sep="");
write.table(text,file=name);
