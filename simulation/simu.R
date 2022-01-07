# Perform simulation for all settings

for(n in c(100,200,400))
{
    for(mp in c(5,10,20,30))
    {
        source('affine.R')
        source('logcholesky.R')
        source('sphere.R')
    }
}
