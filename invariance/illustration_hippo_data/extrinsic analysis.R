#set basic parameter
n=length(which(gin==CNAZ));
Y=oriY[which(gin==CNAZ),,,];
t=orit[which(gin==CNAZ),];
m=orim[which(gin==CNAZ)];
maxm=max(m);


#scatter plot of pairT
dataT=c();
for(i in 1:n){
  for(j1 in 1:m[i]){
    for(j2 in 1:m[i]){
      if(j1!=j2){
        temp=c(t[i,j1],t[i,j2]);
        dataT=rbind(dataT,temp);
      }
    }
  }
}
dataT=as.data.frame(dataT);
ggplot(data=dataT, aes(x=V1, y=V2)) + geom_point() + xlab("age") + xlab("age")


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
#inverse trans, make vector to matrix
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
#trans y into newy by trans
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
    print("error");
    print(s);
    return(array(0,6));
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
    print("error");
    return(array(0,c(6,6)));
  }
  else{
    return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
  }
}

Rotate=function(s){
  P=array(0,c(6,6));
  P[1,1]=cos(4*pi*s);
  P[2,2]=cos(4*pi*s);
  P[3,3]=cos(4*pi*s);
  P[4,4]=cos(4*pi*s);
  P[5,5]=cos(4*pi*s);
  P[6,6]=cos(4*pi*s);
  P[1,2]=sin(4*pi*s);
  P[2,1]=-sin(4*pi*s);
  P[3,4]=sin(4*pi*s);
  P[4,3]=-sin(4*pi*s);
  P[5,6]=sin(4*pi*s);
  P[6,5]=-sin(4*pi*s);
  return(P);
}

#estimate hatc
hmu=10;
logY=array(0,c(n,maxm,6));
for(i in 1:n){
  for(j in 1:m[i]){
    logY[i,j,]=newy[i,j,]-hatmu(t[i,j],n,m,t,newy,hmu);
  }
}
if(rotation==1){#rotate the frame
  for(i in 1:n){
    for(j in 1:m[i]){
      logY[i,j,]=Rotate(t[i,j])%*%logY[i,j,];
    }
  }
}




####estimate C in a grid
hc=20;
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
  print("matrixhatc")
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
  print("matrixeigen")
  print(k);
  for(j in 1:l){
    matrixeigen[k,j,,]=inv(eigen[k,j,]+hatmu(tanchor[j],n,m,t,newy,hmu));
  }
}
