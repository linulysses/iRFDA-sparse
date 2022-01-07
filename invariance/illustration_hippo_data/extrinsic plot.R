library(rgl) 
library(matrixcalc)
library(pracma)

FA <- function(A)
{
    v <- eigen(A,symmetric=T)$values
    vbar <- mean(v)
    sqrt(3/2)*sqrt(sum((v-vbar)^2)) / sqrt(sum(v^2))
}


my.rgl.ellipsoid <- function(x=0, y=0, z=0, 
                             a = 1, b=1, c=1, 
                             rotM=diag(3),
                             subdivide = 3, 
                             smooth = TRUE,
                             col = col)
{
    
    sphere <- rgl::subdivision3d(rgl::cube3d(col=col), subdivide)
    class(sphere) <- c("mesh3d","shape3d")
    
    norm <- sqrt(sphere$vb[1, ]^2 + 
                     sphere$vb[2, ]^2 + 
                     sphere$vb[3, ]^2 )
    for (i in 1:3) sphere$vb[i, ] <- sphere$vb[i, ]/norm
    sphere$vb[4, ] <- 1
    sphere$normals <- sphere$vb
    result <- rgl::scale3d(sphere, a,b,c)
    #rotM <- cpp_euler_passive(phi,theta,psi)
    result <- rgl::rotate3d(result,matrix=rotM)
    result <- rgl::translate3d(result, x,y,z)
    invisible(result)
}




draw <- function(M.list){
    
    N <- length(M.list)
    colours <- rep(0,N)
    
    minf <- 0.01
    maxf <- 0.5
    
    for(i in 1:length(M.list))
    {
        EV1 <- eigen(M.list[[i]])$vectors[,1]
        f <- FA(M.list[[i]])
        
        #EV1[1] <- abs(EV1[1])
        #EV1[2] <- (1+EV1[2])/2
        #EV1[3] <- (1+EV1[3])/2
        #EV1 <- abs(EV1)
        colours[i] <- rgb((maxf-f)/(maxf-minf),0,(f-minf)/(maxf-minf)) #rgb(EV1[1],EV1[2],EV1[3])
    }
    
    ll <- lapply(seq(1,N), function(ii){
        A <- M.list[[ii]]
        s <- svd(A)
        rotM <- s$u
        if(det(rotM) < 0)  rotM[,1] <- -rotM[,1]
        abc <- s$d
        my.rgl.ellipsoid(3.5*ii,(N-ii)*0.4,0,abc[1],abc[2],abc[3],
                         rotM=rotM, col = colours[ii])})
    
    rgl::shapelist3d(ll)
    
}

draw.ellipsoid <- function(A)
{
    s <- svd(A)
    rotM <- s$u
    abc <- s$d
    my.rgl.ellipsoid(0,0,0,abc[1],abc[2],abc[3],
                     rotM=rotM, col = 'red')
}


get.list <- function(X,k)
{
    L <- list()
    p <- dim(X)[2]
    for(i in 1:p)
    {
        L[[i]] <- X[k,i,,]
    }
    L
}

## Find ranage of FA
load('extrinsic_AZ_rotation0.RData')
Z <- matrixeigen[1:3,seq(1,99,7),,]
fresult <- NULL
for(k in 1:3)
{
    M.list <- get.list(Z,k)
    fresult <- c(fresult,sapply(M.list,FA))
}
fresult <- c(fresult,sapply(M.list,FA))
load('extrinsic_AZ_rotation1.RData')
Z <- matrixeigen[1:3,seq(1,99,7),,]
for(k in 1:3)
{
    M.list <- get.list(Z,k)
    fresult <- c(fresult,sapply(M.list,FA))
}
fresult <- c(fresult,sapply(M.list,FA))
range(fresult)


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
    
    dev.new(width=0.2, height=5)
    par(mar=c(0,4,0,0))
    plot(c(0,1.5), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,0.8,y+1/scale, col=lut[i], border=NA)
    }
}

color.bar(colorRampPalette(c("red",'blue'))(100), min=0.01,max=0.5,nticks=6)



## Plot eigenfunctions

load('3d.par.2.RData')
load('extrinsic_AZ_rotation0.RData')
rotation0EF <- matrixeigen[1:3,seq(1,99,7),,]
load('extrinsic_AZ_rotation1.RData')
rotation1EF <- matrixeigen[1:3,seq(1,99,7),,]

for(k in 1:3)
{
    open3d()
    par3d(windowRect = c(20, 30, 800, 120))
    M.list <- get.list(rotation0EF,k)
    draw(M.list)
    par3d(my.par)
    rgl.postscript(paste0('extrinsic_AZ_rotation0',k,'.eps'), fmt = "eps")
    rgl.snapshot(paste0('extrinsic_AZ_rotation0',k,'.png'))
    
    open3d()
    par3d(windowRect = c(20, 30, 800, 120))
    M.list <- get.list(rotation1EF,k)
    draw(M.list)
    par3d(my.par)
    rgl.postscript(paste0('extrinsic_AZ_rotation1',k,'.eps'), fmt = "eps")
    rgl.snapshot(paste0('extrinsic_AZ_rotation1',k,'.png'))
}    
