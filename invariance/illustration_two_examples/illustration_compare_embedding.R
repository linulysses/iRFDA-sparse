# Simulation: sphere case
# n samples, mi satisfies Possion(mp)+2
mp=10;
n=100;
#pre-deterministic hc
hccandidatelen=5;
hccandidate=exp(seq(log(0.1),log(0.5),length.out = hccandidatelen));
# estimation over a grid
l=10;
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
  newc=array(0,c(3,3));
  newc[1,1]=A[1]*B[1];
  newc[1,2]=A[1]*B[2];
  newc[1,3]=A[1]*B[3];
  newc[2,1]=A[2]*B[1];
  newc[2,2]=A[2]*B[2];
  newc[2,3]=A[2]*B[3];
  newc[3,1]=A[3]*B[1];
  newc[3,2]=A[3]*B[2];
  newc[3,3]=A[3]*B[3];
  return(newc);
}

#hatc
hatc=function(s1,s2,n,m,t,logY,hc){
  S00=0;
  S01=0;
  S10=0;
  S11=0;
  S20=0;
  S02=0;
  R00=array(0,c(3,3));
  R01=array(0,c(3,3));
  R10=array(0,c(3,3));
  for(i in 1:n){
    for(j1 in 1:m[i]){
      for(j2 in 1:m[i]){
        if((abs(t[i,j1]-s1)<hc)&&(abs(t[i,j2]-s2)<hc)&&(j1!=j2)){
          #logY[i,j,] is a vector of dimension 2 with respect to coordinate B1,B2
          #logY[i,j,] is on T_{hat{mu}(t[i,j])}, which is Log_{hat{mu}(t[i,j])}Y(t[i,j])
          #parallel transport logY[i,j1,] to T_{hat{mu}(s1)} and then T_{mu(s1)}
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
    return(array(0,c(3,3)));
  }
  else{
    return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
  }
}

#real C under differen frames
realc1=array(0,c(3,3));
realc1[1,1]=1/300;
realc1[2,2]=1/300;
realc1[3,3]=0;

realc2=array(0,c(3,3));
realc2[1,1]=0;
realc2[2,2]=1/300;
realc2[3,3]=1/300;

realc3=array(0,c(3,3));
realc3[1,1]=1/300;
realc3[2,2]=0;
realc3[3,3]=1/300;

realnorm=sum(realc1**2);

#estimate C
#Mentor Carlo for itermax=100 times
itermax=100;
normlist=array(0,c(itermax,hccandidatelen,3));
for(it in 41:itermax){
  cat("it:",it,"\n")
  #generate m,t
  m=rpois(n,mp)+2;
  mmax=max(m);
  N=sum(m);
  t=array(0,c(n,mmax));
  for(i in 1:n){
    t[i,]=runif(mmax,0,1);
  }
  #generate data logmuY
  #logmuY[i,j,] is the vector of Log_{mu(t[i,j])}Y_{i,j} with respect to B1,B2
  #Y_{i,j} is restored in the form of Log_{mu(t[i,j])}Y_{i,j} rather than a point in R^{3}
  Z1=runif(n,-0.1,0.1);
  Z2=runif(n,-0.1,0.1);
  manicoor=array(0,c(n,mmax,2));
  for(i in 1:n){
    for(j in 1:m[i]){
      manicoor[i,j,1]=0.25+0.5*t[i,j]+Z1[i];
      manicoor[i,j,2]=0.5+Z2[i];
    }
  }
  #asuume mu is known, set hatmu=mu
  #the extrinsic method under three embeddings
  newcoor=array(0,c(3,n,mmax,3));
  #embedding 1
  for(i in 1:n){
    for(j in 1:m[i]){
      newcoor[1,i,j,1]=Z1[i];
      newcoor[1,i,j,2]=Z2[i];
      newcoor[1,i,j,3]=0;
    }
  }
  #embedding 2
  for(i in 1:n){
    for(j in 1:m[i]){
      newcoor[2,i,j,1]=Z1[i]*cos(pi*(0.25+0.5*t[i,j]));
      newcoor[2,i,j,2]=-Z1[i]*sin(pi*(0.25+0.5*t[i,j]));
      newcoor[2,i,j,3]=Z2[i];
    }
  }
  #embedding 3
  for(i in 1:n){
    for(j in 1:m[i]){
      newcoor[3,i,j,1]=Z1[i]*cos(2*pi*(0.25+0.5*t[i,j]));
      newcoor[3,i,j,2]=-Z1[i]*sin(2*pi*(0.25+0.5*t[i,j]));
      newcoor[3,i,j,3]=Z2[i];
    }
  }
  #for different hc
  for(hcit in 1:hccandidatelen){
    cat("hcit:",hcit,"\n")
    hc=hccandidate[hcit];
    #record the l2 norm
    for(k1 in 1:l){
      for(k2 in 1:l){
        s1=tanchor[k1];
        s2=tanchor[k2];
        normlist[it,hcit,1]=normlist[it,hcit,1]+sum((hatc(s1,s2,n,m,t,newcoor[1,,,],hc)-realc1)**2);
        normlist[it,hcit,2]=normlist[it,hcit,2]+sum((hatc(s1,s2,n,m,t,newcoor[2,,,],hc)-realc2)**2);
        normlist[it,hcit,3]=normlist[it,hcit,3]+sum((hatc(s1,s2,n,m,t,newcoor[3,,,],hc)-realc3)**2);
      }
    }
    normlist[it,hcit,1]=normlist[it,hcit,1]/(l*l);
    normlist[it,hcit,2]=normlist[it,hcit,2]/(l*l);
    normlist[it,hcit,3]=normlist[it,hcit,3]/(l*l);
  }
}
  
#relative error
errortable=array(0,c(hccandidatelen,3,2));
for(hcit in 1:hccandidatelen){
  for(me in 1:3){
    errortable[hcit,me,1]=sqrt(mean(normlist[,hcit,me]/realnorm));
    errortable[hcit,me,2]=sqrt(sd(normlist[,hcit,me]/realnorm));
  }
}
print(errortable[,,1])
print(errortable[,,2])










