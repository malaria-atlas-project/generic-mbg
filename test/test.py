from __future__ import division
import pymc as pm
import numpy as np
import pylab as pl

def gencirc(x,r):
    out = np.zeros(x.shape[0])
    out[np.where(x[:,0]**2+x[:,1]**2<r**2)]=1
    return out
    
N = 300
r1 = .4
r2 = .6
c1 = -2
c2 = 1

if __name__ == '__main__':
    
    # Make covariate asciis
    lon_ax = np.linspace(-1.5,1.5,500)
    lat_ax = np.linspace(-1.5,1.5,500)
    lon,lat = np.meshgrid(lon_ax, lat_ax)
    x = np.dstack((lon,lat)).reshape((-1,2))
    cov1 = gencirc(x,r1).reshape((500,500))
    cov2 = gencirc(x,r2).reshape((500,500))    
    
    import map_utils
    map_utils.exportAscii2(lon_ax, lat_ax, np.ma.masked_array(cov1,mask=np.zeros((500,500))), 'cov1.asc')
    map_utils.exportAscii2(lon_ax, lat_ax, np.ma.masked_array(cov2,mask=np.zeros((500,500))), 'cov2.asc')    
            
    # Make fake data
    lon = np.random.uniform(-1,1,size=N)
    lat = np.random.uniform(-1,1,size=N)
    t = np.zeros(N)
    
    x = np.vstack((lon,lat,t)).T
    
    cov1 = gencirc(x,r1)
    cov2 = gencirc(x,r2)
    
    M = c1*cov1+c2*cov2
    S = pm.gp.FullRankCovariance(pm.gp.cov_funs.exponential.aniso_geo_rad, amp=.5, scale=.08, inc=.5, ecc=.5).cholesky(x[:,:2])
    
    y = pm.rmv_normal_chol(M,S.T)+np.random.normal(N)*.1
    z = pm.flib.invlogit(y)
    
    lo_age = np.ones(N)*2
    up_age = np.ones(N)*10
    n = np.random.randint(10,500,size=N)
    pos = pm.rbinomial(n, z)
    neg = n-pos
    
    data_file = np.rec.fromarrays([pos,neg,lo_age,up_age,lon,lat,t,cov1,cov2],names='pos,neg,lo_age,up_age,lon,lat,t,cov1,cov2')


# where_0 = np.where(M==0)
# where_1 = np.where(M==1)
# where_n1 = np.where(M==-1)
# 
# pl.figure(1)
# pl.clf()
# pl.hist(z[where_0])
# 
# pl.figure(2)
# pl.clf()
# pl.hist(z[where_n1])
# 
# pl.figure(3)
# pl.clf()
# pl.hist(z[where_1])
# 
# pl.figure(4)
# pl.clf()
# pl.plot(x[where_0,0],x[where_0,1],'b.')