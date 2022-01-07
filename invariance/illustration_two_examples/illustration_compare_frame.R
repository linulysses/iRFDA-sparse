# Simulation: sphere case
# n samples, mi satisfies Possion(mp)+2
mp=10;
n=100;
#pre-deterministic hc
hccandidatelen=5;
hccandidate=exp(seq(log(0.2),log(0.8),length.out = hccandidatelen));
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
  newc=array(0,c(2,2));
  newc[1,1]=A[1]*B[1];
  newc[1,2]=A[1]*B[2];
  newc[2,1]=A[2]*B[1];
  newc[2,2]=A[2]*B[2];
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
  R00=array(0,c(2,2));
  R01=array(0,c(2,2));
  R10=array(0,c(2,2));
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
    return(array(0,c(2,2)));
  }
  else{
    return(((S20021111*R00-S10020111*R10+S10110120*R01)/(S20021111*S00-S10020111*S10+S10110120*S01)));
  }
}

#real vales of C(s1,s2)
#all the frames share the following same matrix representation
realc=function(s1,s2){
  newc=array(0,c(2,2));
  newc[1,1]=s1*s2/300;
  newc[2,2]=s1*s2/300;
  return(newc);
}
realnorm=0; 
for(k1 in 1:l){
  for(k2 in 1:l){
    s1=tanchor[k1];
    s2=tanchor[k2];
    realnorm=realnorm+sum((realc(s1,s2))**2)
  }
}
realnorm=realnorm/(l*l);

#estimate C
#Mentor Carlo for itermax=100 times
itermax=100;
normlist=array(0,c(itermax,hccandidatelen,3));
for(it in 1:itermax){
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
  logmuY=array(0,c(n,mmax,2));
  for(i in 1:n){
    for(j in 1:m[i]){
      logmuY[i,j,1]=t[i,j]*Z1[i];
      logmuY[i,j,2]=t[i,j]*Z2[i];
    }
  }
  #assume mu is known, set hatmu=mu
  #the extrinsic method under rotation 4pikt
  loghatmuY=array(0,c(3,n,mmax,2));
  #frame 1
  for(i in 1:n){
    for(j in 1:m[i]){
      loghatmuY[1,i,j,1]=logmuY[i,j,1];
      loghatmuY[1,i,j,2]=logmuY[i,j,2];
    }
  }
  #frame 2
  for(i in 1:n){
    for(j in 1:m[i]){
      loghatmuY[2,i,j,1]=cos(4*pi*t[i,j])*logmuY[i,j,1]+sin(4*pi*t[i,j])*logmuY[i,j,2];
      loghatmuY[2,i,j,2]=sin(4*pi*t[i,j])*logmuY[i,j,1]+cos(4*pi*t[i,j])*logmuY[i,j,2];
    }
  }
  #frame 3
  for(i in 1:n){
    for(j in 1:m[i]){
      loghatmuY[3,i,j,1]=sin(0.5*pi*t[i,j])*logmuY[i,j,1]+cos(0.5*pi*t[i,j])*logmuY[i,j,2];
      loghatmuY[3,i,j,2]=cos(0.5*pi*t[i,j])*logmuY[i,j,1]-sin(0.5*pi*t[i,j])*logmuY[i,j,2];
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
        normlist[it,hcit,1]=normlist[it,hcit,1]+sum((hatc(s1,s2,n,m,t,loghatmuY[1,,,],hc)-realc(s1,s2))**2);
        normlist[it,hcit,2]=normlist[it,hcit,2]+sum((hatc(s1,s2,n,m,t,loghatmuY[2,,,],hc)-realc(s1,s2))**2);
        normlist[it,hcit,3]=normlist[it,hcit,3]+sum((hatc(s1,s2,n,m,t,loghatmuY[3,,,],hc)-realc(s1,s2))**2);
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










