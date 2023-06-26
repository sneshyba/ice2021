# Code for roughness statistics
import numpy as np
import copy
import imagestuff as ims
from scipy.interpolate import griddata
#import importlib; importlib.reload(ims)



def pWeibull(r, sigma, eta):
    ''' Weibull function '''
    from numpy import exp
    mu = 1-r
    ret = 2*eta/sigma**2/mu**3 * \
        (((mu**(-2)-1)/sigma**2)**(eta-1)) * \
        exp(-((mu**(-2)-1)/sigma**2)**eta)
    return ret

def pWeibullr(r, sigma, eta):
    ''' Weibull function times r '''
    return pWeibull(r, sigma, eta)*r

def pGaussian(r, sigma):
    ''' Gaussian function '''
    return pWeibull(r, sigma, 1)    

def pGaussianr(r, sigma):
    ''' Gaussian function times r '''
    return pWeibullr(r, sigma, 1)

def bimodal(r, sigma1, sigma2, N):
    ''' Bimodal Gaussian function '''
    pdf1 = pWeibull(r,sigma1,1.0)
    pdf2 = pWeibull(r,sigma2,1.0)
    return (1-N)*pdf1 + N*pdf2 

def bimodalr(r, sigma1, sigma2, N):
    ''' Bimodal Gaussian function times r'''
    pdf1 = pWeibullr(r,sigma1,1.0)
    pdf2 = pWeibullr(r,sigma2,1.0)
    return (1-N)*pdf1 + N*pdf2 


def bimodalfunc(r, sigma1, sigma2, N):
    ''' Bimodal Gaussian function '''
    pdf1 = pWeibullr(r,sigma1,1.0)
    pdf2 = pWeibullr(r,sigma2,1.0)
    return (1-N)*pdf1 + N*pdf2 


def sigma2meanr(sigma):
    ''' Converting sigma to <r> 
        Usage: 
        
        sigmalist = np.linspace(.01,.9)
        meanr = sigma2meanr(sigmalist)
        plt.figure()
        plt.plot(sigmalist,meanr,'o')
        plt.grid(True)        
    '''
    p = np.array([ 4.57899291e-01, -2.27236062e+00,  4.72080621e+00, -5.09338608e+00,
        2.57626515e+00,  1.77811151e-01, -8.38705493e-01,  1.49765369e-02,
        4.98822342e-01,  3.87112620e-05, -3.41914362e-07])
    meanr = np.polyval(p,sigma)
    return meanr

def R_squar(y,yfit):
    SS_res = np.sum((y-yfit)**2)
    SS_tot = np.sum((y-np.mean(y))**2)
    return 1-SS_res/SS_tot

# def makehistogram(\
#                   nsegments,nx1list,nx2list,ny1list,ny2list,dx,dy,solution,\
#                   accumlist, newrbins):

#     # Setting up null arrays
#     hbins_accum = []
#     meanrsub_accum = []
#     zsigma_accum = []
#     Z2_accum = []
#     Zsquared_accum = []
#     rsub_accum = []
#     meanrsublist = []
#     Zsigmalist = []
#     Z2list = []

#     # Now, to evaluate the roughness ... First step is to flatten each panel via rotation
#     # Here we explicitly flip the y-coordinate (to make it a right-handed system) so we don't have to invert on the fly

#     for isegment in range(0,nsegments):

#         # Extract this segment
#         nx1=nx1list[isegment]; nx2=nx2list[isegment]; nxsegment = nx2-nx1+1
#         ny1=ny1list[isegment]; ny2=ny2list[isegment]; nysegment = ny2-ny1+1
#         surf_xseg = np.linspace(0,(nxsegment-1)*dx,nxsegment); 
#         surf_yseg = np.linspace(0,(nysegment-1)*dy,nysegment); 
#         surf_xseggrid, surf_yseggrid = np.meshgrid(surf_xseg,surf_yseg) # 1st index is y, 2nd is x
#         surf_zseggrid = copy.copy(np.flipud(solution[ny1:ny2+1,nx1:nx2+1])) # This flips the y-coordinate

#         # Fit a plane to the data and adjust data to start at the origin
#         m = ims.polyfit2d(\
#                       surf_xseggrid.reshape(nysegment*nxsegment), \
#                       surf_yseggrid.reshape(nysegment*nxsegment), \
#                       surf_zseggrid.reshape(nysegment*nxsegment), \
#                       linear=True,order=1)

#         # Get the angles of the plane
#         dzdy = m[1]; thetay = np.arctan(dzdy)*180/np.pi; #print 'y:', thetay

#         # Get rotation matrix & flatten in one direction
#         Roty = ims.myrotation_matrix([1,0,0], -thetay)
#         surf_xseggridp, surf_yseggridp, surf_zseggridp = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid, Roty)

#         # Fit a plane to the data and adjust data to start at the origin
#         mp = ims.polyfit2d(\
#                       surf_xseggridp.reshape(nysegment*nxsegment), \
#                       surf_yseggridp.reshape(nysegment*nxsegment), \
#                       surf_zseggridp.reshape(nysegment*nxsegment), \
#                       linear=True,order=1)

#         # Get the angle of the plane in another direction
#         dzdx = mp[2]; thetaxp = np.arctan(dzdx)*180/np.pi; #print 'x:', thetaxp

#         # Get rotation matrix & flatten in another direction
#         Rotxp = ims.myrotation_matrix([0,1,0], thetaxp)
#         surf_xseggridpp, surf_yseggridpp, surf_zseggridpp = \
#             ims.flatten(surf_xseggridp, surf_yseggridp, surf_zseggridp, Rotxp)


#         # Trying out the polyval2d
#         surf_zseggrid_theory_long = ims.polyval2d(\
#                       surf_xseggrid.reshape(nysegment*nxsegment), \
#                       surf_yseggrid.reshape(nysegment*nxsegment), \
#                       m)
#         surf_zseggrid_theory = surf_zseggrid_theory_long.reshape(nysegment,nxsegment)
#         surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid_theory, Roty)
#         surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory = \
#             ims.flatten(surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory, Rotxp)

#         # Now rotate
#         deltay = surf_yseggridpp_theory[0,-1]-surf_yseggridpp_theory[0,0]
#         deltax = surf_xseggridpp_theory[0,-1]-surf_xseggridpp_theory[0,0]
#         thetazpp = -np.arctan(deltay/deltax)*180/np.pi;
#         Rotzpp = ims.myrotation_matrix([0,0,1], thetazpp)
#         surf_xseggridppp, surf_yseggridppp, surf_zseggridppp = \
#             ims.flatten(surf_xseggridpp, surf_yseggridpp, surf_zseggridpp, Rotzpp)
#         surf_xseggridppp_theory, surf_yseggridppp_theory, surf_zseggridppp_theory = \
#             ims.flatten(surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory, Rotzpp)

#         # Now we have to extract an orthogonal subset
#         dxsub = dysub = dx
#         xsubstart = np.max(surf_xseggridppp_theory[[0,-1],0])+dxsub*2
#         xsubstop = np.min(surf_xseggridppp_theory[[0,-1],-1])-dxsub*2
#         ysubstart = np.max(surf_yseggridppp_theory[0,[0,-1]])+dysub*2
#         ysubstop = np.min(surf_yseggridppp_theory[-1,[0,-1]])-dysub*2
#         xsub = np.arange(xsubstart,xsubstop,dxsub)
#         ysub = np.arange(ysubstart,ysubstop,dysub)
#         sub_xseggrid, sub_yseggrid = np.meshgrid(xsub,ysub) # 1st index is y, 2nd is x
#         nsuby, nsubx = np.shape(sub_xseggrid)
#         surf_xseggridppp_theory_long = np.reshape(surf_xseggridppp_theory,nysegment*nxsegment)
#         surf_yseggridppp_theory_long = np.reshape(surf_yseggridppp_theory,nysegment*nxsegment)
#         points = np.vstack((surf_xseggridppp_theory_long,surf_yseggridppp_theory_long)).T # rows are x,y pairs
#         values = np.reshape(surf_zseggridppp,nysegment*nxsegment)
#         sub_zseggrid_long = griddata(points, values, (sub_xseggrid, sub_yseggrid), method='cubic')
#         sub_zseggrid = np.reshape(sub_zseggrid_long,(nsuby, nsubx))

#         # Now we get the roughness
#         dzsub_dx = np.diff(sub_zseggrid,axis=1)/np.diff(sub_xseggrid,axis=1)
#         dzsub_dy = np.diff(sub_zseggrid,axis=0)/np.diff(sub_yseggrid,axis=0)
#         Zsquared = dzsub_dx[1:,:]**2+dzsub_dy[:,1:]**2
#         rsub = 1.0 - 1/np.sqrt(1+Zsquared)
#         mu = 1-rsub
#         phi = np.arccos(mu)
#         Zplus = Zsquared**.5
#         Z = np.hstack((Zplus,-Zplus)) # Need +/- to generate a two-sided distribution
#         thismeanrsub = np.round(np.mean(rsub)*1000)/1000; meanrsublist.append(thismeanrsub)
#         thissigma = np.round(np.std(Z)*100)/100; Zsigmalist.append(thissigma)
#         thismeanZ2 = np.mean(Zsquared); Z2list.append(thismeanZ2)

#         # Numerical distribution functions
#         rsub_long = np.reshape(rsub,np.size(rsub))
#         hist = np.histogram(rsub_long,bins=newrbins)
#         rbins = hist[1][0:-1]
#         rbins1 = hist[1][1:]
#         hbins = hist[0] 

#         # Normalizing ... this might be wrong
#         norm = np.trapz(hbins,rbins)
#         hbins = hbins/norm

#         # Defining the analytical distribution function bins
#         rwidth = rbins1-rbins
#         rbinsW = (rbins+rwidth/2.0)        

#         # Accumulate the binned data
#         if isegment in accumlist:
#             hbins_accum.append(hbins)
#             print ('Accumulating ...', isegment+1, 'out of', len(accumlist))

#     # Combine the histograms of individual segments
#     hbins_total = np.sum((hbins_accum),axis=0)/len(accumlist)
#     norm = np.trapz(hbins_total,np.log(rbinsW)); print('Norm =', norm)
#     hbins_total = hbins_total/norm

#     # Get out
#     return hbins_total, rbinsW


# def makehistogram_heights(\
#                   nsegments,nx1list,nx2list,ny1list,ny2list,dx,dy,solution,\
#                   accumlist, newzbins):

#     # Setting up null arrays
#     hbins_accum = []
#     z_accum = []
#     zlist = []

#     # Looping over all the segments
#     for isegment in range(0,nsegments):

#         # Extract this segment
#         nx1=nx1list[isegment]; nx2=nx2list[isegment]; nxsegment = nx2-nx1+1
#         ny1=ny1list[isegment]; ny2=ny2list[isegment]; nysegment = ny2-ny1+1
#         surf_xseg = np.linspace(0,(nxsegment-1)*dx,nxsegment); 
#         surf_yseg = np.linspace(0,(nysegment-1)*dy,nysegment); 
#         surf_xseggrid, surf_yseggrid = np.meshgrid(surf_xseg,surf_yseg) # 1st index is y, 2nd is x
#         surf_zseggrid = copy.copy(np.flipud(solution[ny1:ny2+1,nx1:nx2+1])) # This flips the y-coordinate

#         # Fit a plane to the data and adjust data to start at the origin
#         m = ims.polyfit2d(\
#                       surf_xseggrid.reshape(nysegment*nxsegment), \
#                       surf_yseggrid.reshape(nysegment*nxsegment), \
#                       surf_zseggrid.reshape(nysegment*nxsegment), \
#                       linear=True,order=1)

#         # Get the angles of the plane
#         dzdy = m[1]; thetay = np.arctan(dzdy)*180/np.pi; #print 'y:', thetay

#         # Get rotation matrix & flatten in one direction
#         Roty = ims.myrotation_matrix([1,0,0], -thetay)
#         surf_xseggridp, surf_yseggridp, surf_zseggridp = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid, Roty)

#         # Fit a plane to the data and adjust data to start at the origin
#         mp = ims.polyfit2d(\
#                       surf_xseggridp.reshape(nysegment*nxsegment), \
#                       surf_yseggridp.reshape(nysegment*nxsegment), \
#                       surf_zseggridp.reshape(nysegment*nxsegment), \
#                       linear=True,order=1)

#         # Get the angle of the plane in another direction
#         dzdx = mp[2]; thetaxp = np.arctan(dzdx)*180/np.pi; #print 'x:', thetaxp

#         # Get rotation matrix & flatten in another direction
#         Rotxp = ims.myrotation_matrix([0,1,0], thetaxp)
#         surf_xseggridpp, surf_yseggridpp, surf_zseggridpp = \
#             ims.flatten(surf_xseggridp, surf_yseggridp, surf_zseggridp, Rotxp)


#         # Trying out the polyval2d
#         surf_zseggrid_theory_long = ims.polyval2d(\
#                       surf_xseggrid.reshape(nysegment*nxsegment), \
#                       surf_yseggrid.reshape(nysegment*nxsegment), \
#                       m)
#         surf_zseggrid_theory = surf_zseggrid_theory_long.reshape(nysegment,nxsegment)
#         #surf_zseggrid_theory -= z0
#         surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid_theory, Roty)
#         surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory = \
#             ims.flatten(surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory, Rotxp)

#         # Now rotate
#         deltay = surf_yseggridpp_theory[0,-1]-surf_yseggridpp_theory[0,0]
#         deltax = surf_xseggridpp_theory[0,-1]-surf_xseggridpp_theory[0,0]
#         thetazpp = -np.arctan(deltay/deltax)*180/np.pi;
#         Rotzpp = ims.myrotation_matrix([0,0,1], thetazpp)
#         surf_xseggridppp, surf_yseggridppp, surf_zseggridppp = \
#             ims.flatten(surf_xseggridpp, surf_yseggridpp, surf_zseggridpp, Rotzpp)
#         surf_xseggridppp_theory, surf_yseggridppp_theory, surf_zseggridppp_theory = \
#             ims.flatten(surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory, Rotzpp)

#         # Now we have to extract an orthogonal subset
#         dxsub = dysub = dx
#         xsubstart = np.max(surf_xseggridppp_theory[[0,-1],0])+dxsub*2
#         xsubstop = np.min(surf_xseggridppp_theory[[0,-1],-1])-dxsub*2
#         ysubstart = np.max(surf_yseggridppp_theory[0,[0,-1]])+dysub*2
#         ysubstop = np.min(surf_yseggridppp_theory[-1,[0,-1]])-dysub*2
#         xsub = np.arange(xsubstart,xsubstop,dxsub)
#         ysub = np.arange(ysubstart,ysubstop,dysub)
#         sub_xseggrid, sub_yseggrid = np.meshgrid(xsub,ysub) # 1st index is y, 2nd is x
#         nsuby, nsubx = np.shape(sub_xseggrid)
#         surf_xseggridppp_theory_long = np.reshape(surf_xseggridppp_theory,nysegment*nxsegment)
#         surf_yseggridppp_theory_long = np.reshape(surf_yseggridppp_theory,nysegment*nxsegment)
#         points = np.vstack((surf_xseggridppp_theory_long,surf_yseggridppp_theory_long)).T # rows are x,y pairs
#         values = np.reshape(surf_zseggridppp,nysegment*nxsegment)
#         sub_zseggrid_long = griddata(points, values, (sub_xseggrid, sub_yseggrid), method='cubic')
#         sub_zseggrid = np.reshape(sub_zseggrid_long,(nsuby, nsubx))

#         # Now we get the heights relative to a reference
#         zreference = np.median(sub_zseggrid)
        
#         # Numerical distribution functions
#         hist = np.histogram(sub_zseggrid_long-zreference,bins=newzbins)
#         zbins = hist[1][0:-1]
#         zbins1 = hist[1][1:]
#         hbins = hist[0] 
#         norm = np.trapz(hbins,zbins)
#         hbins = hbins/norm

#         # Defining the analytical distribution function bins
#         zwidth = zbins1-zbins
#         zbinsW = (zbins+zwidth/2.0)        

#         # Accumulate the binned data
#         if isegment in accumlist:
#             hbins_accum.append(hbins)
#             print ('Accumulating ...', isegment+1, 'out of', len(accumlist))

#     # Combine the histograms of individual segments
#     hbins_total = np.sum((hbins_accum),axis=0)/len(accumlist)
#     norm = np.trapz(hbins_total,zbinsW); print('Norm =', norm)
#     hbins_total = hbins_total/norm

#     # Get out
#     return hbins_total, zbinsW


# def get_heights(nsegments,nx1list,nx2list,ny1list,ny2list,dx,dy,solution,isegment):

#         # Extract this segment
#         nx1=nx1list[isegment]; nx2=nx2list[isegment]; nxsegment = nx2-nx1+1
#         ny1=ny1list[isegment]; ny2=ny2list[isegment]; nysegment = ny2-ny1+1
#         surf_xseg = np.linspace(0,(nxsegment-1)*dx,nxsegment); 
#         surf_yseg = np.linspace(0,(nysegment-1)*dy,nysegment); 
#         surf_xseggrid, surf_yseggrid = np.meshgrid(surf_xseg,surf_yseg) # 1st index is y, 2nd is x
#         surf_zseggrid = copy.copy(np.flipud(solution[ny1:ny2+1,nx1:nx2+1])) # This flips the y-coordinate

#         # Fit a plane to the data and adjust data to start at the origin
#         m = ims.polyfit2d(surf_xseggrid.reshape(nysegment*nxsegment), \
#                           surf_yseggrid.reshape(nysegment*nxsegment), \
#                           surf_zseggrid.reshape(nysegment*nxsegment), \
#                           linear=True,order=1)

#         # Get the angles of the plane
#         dzdy = m[1]; thetay = np.arctan(dzdy)*180/np.pi; #print 'y:', thetay

#         # Get rotation matrix & flatten in one direction
#         Roty = ims.myrotation_matrix([1,0,0], -thetay)
#         surf_xseggridp, surf_yseggridp, surf_zseggridp = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid, Roty)

#         # Fit a plane to the data and adjust data to start at the origin
#         mp = ims.polyfit2d(surf_xseggridp.reshape(nysegment*nxsegment), \
#                            surf_yseggridp.reshape(nysegment*nxsegment), \
#                            surf_zseggridp.reshape(nysegment*nxsegment), \
#                            linear=True,order=1)

#         # Get the angle of the plane in another direction
#         dzdx = mp[2]; thetaxp = np.arctan(dzdx)*180/np.pi; #print 'x:', thetaxp

#         # Get rotation matrix & flatten in another direction
#         Rotxp = ims.myrotation_matrix([0,1,0], thetaxp)
#         surf_xseggridpp, surf_yseggridpp, surf_zseggridpp = \
#             ims.flatten(surf_xseggridp, surf_yseggridp, surf_zseggridp, Rotxp)


#         # Trying out the polyval2d
#         surf_zseggrid_theory_long = ims.polyval2d( \
#                                         surf_xseggrid.reshape(nysegment*nxsegment), \
#                                         surf_yseggrid.reshape(nysegment*nxsegment), \
#                                         m)
#         surf_zseggrid_theory = surf_zseggrid_theory_long.reshape(nysegment,nxsegment)
#         #surf_zseggrid_theory -= z0
#         surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory = \
#             ims.flatten(surf_xseggrid, surf_yseggrid, surf_zseggrid_theory, Roty)
#         surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory = \
#             ims.flatten(surf_xseggridp_theory, surf_yseggridp_theory, surf_zseggridp_theory, Rotxp)

#         # Now rotate
#         deltay = surf_yseggridpp_theory[0,-1]-surf_yseggridpp_theory[0,0]
#         deltax = surf_xseggridpp_theory[0,-1]-surf_xseggridpp_theory[0,0]
#         thetazpp = -np.arctan(deltay/deltax)*180/np.pi;
#         Rotzpp = ims.myrotation_matrix([0,0,1], thetazpp)
#         surf_xseggridppp, surf_yseggridppp, surf_zseggridppp = \
#             ims.flatten(surf_xseggridpp, surf_yseggridpp, surf_zseggridpp, Rotzpp)
#         surf_xseggridppp_theory, surf_yseggridppp_theory, surf_zseggridppp_theory = \
#             ims.flatten(surf_xseggridpp_theory, surf_yseggridpp_theory, surf_zseggridpp_theory, Rotzpp)

#         # Now we have to extract an orthogonal subset
#         dxsub = dysub = dx
#         xsubstart = np.max(surf_xseggridppp_theory[[0,-1],0])+dxsub*2
#         xsubstop = np.min(surf_xseggridppp_theory[[0,-1],-1])-dxsub*2
#         ysubstart = np.max(surf_yseggridppp_theory[0,[0,-1]])+dysub*2
#         ysubstop = np.min(surf_yseggridppp_theory[-1,[0,-1]])-dysub*2
#         xsub = np.arange(xsubstart,xsubstop,dxsub)
#         ysub = np.arange(ysubstart,ysubstop,dysub)
#         sub_xseggrid, sub_yseggrid = np.meshgrid(xsub,ysub) # 1st index is y, 2nd is x
#         nsuby, nsubx = np.shape(sub_xseggrid)
#         surf_xseggridppp_theory_long = np.reshape(surf_xseggridppp_theory,nysegment*nxsegment)
#         surf_yseggridppp_theory_long = np.reshape(surf_yseggridppp_theory,nysegment*nxsegment)
#         points = np.vstack((surf_xseggridppp_theory_long,surf_yseggridppp_theory_long)).T # rows are x,y pairs
#         values = np.reshape(surf_zseggridppp,nysegment*nxsegment)
#         sub_zseggrid_long = griddata(points, values, (sub_xseggrid, sub_yseggrid), method='cubic')
#         sub_zseggrid = np.reshape(sub_zseggrid_long,(nsuby, nsubx))

#         # Now we get the heights relative to a reference
#         zreference = np.median(sub_zseggrid)
        
# #         # Accumulate the binned data
# #         if isegment in accumlist:
            

#         # Get out
#         return sub_zseggrid

    
# def getrhoofz2(sollast_in,dx,dy,nbins=10,transposeflag=False):
    
#     # Transpose, if flagged
#     if transposeflag:
#         sollast = sollast_in.T
#     else:
#         sollast = sollast_in
    
#     # Dimensions 
#     Nx, Ny = np.shape(sollast)
    
#     # Calculate the gradient squared (Z2)
#     dzdx = np.diff(sollast, axis=0)/dx
#     dzdy = np.diff(sollast, axis = 1)/dy #we are not sure which axis is which
#     Z2 = dzdx[:, 1:]**2+dzdy[1:, :]**2
    
#     # Get the probability distribution
#     Z2flat = np.reshape(Z2, (Nx-1)*(Ny-1))
#     counts, bins = np.histogram(Z2flat,bins=nbins)
#     counts = counts/np.sum(counts)
#     newbins = bins[1:]
# #     subset = np.array([i for i in range(4,len(bins))])-1
# #     logcounts = np.log(counts[subset])

#     #plt.semilogy(newbins, counts, 'o', label='Numerical result')
#     return counts, newbins

def randomcorrelation(nacross,ndown,nsamples=1000):
    d1d2_list = []
    for i in range(nsamples):
        d1 = np.random.randn(nacross,ndown); #print('d1=',d1)
        d2 = np.random.randn(nacross,ndown); #print('d2=',d2)
        d1d2 = np.mean(d1*d2); #print('random prediction = ', d1d2*percent)
        d1d2_list.append(d1d2)
    return np.std(d1d2_list)

def getinfomatrix(cseg):
    size = len(cseg)
    csegrel = []
    for i in range(len(cseg)):
        thisone = (cseg[i]-np.mean(cseg[i])) / np.std(cseg[i])
        csegrel.append(thisone)
    infomatrix = np.matrix(np.zeros(shape=(size,size)))
    for i in range(len(csegrel)):
        for j in range(i,len(csegrel)):
            product = np.array(csegrel[i])*np.array(csegrel[j])
            infomatrix[i,j] = np.mean(product)
    return infomatrix

def getinfoscore(cseg):
    info = getinfomatrix(cseg)
    score = 0; count = 0
    for i in range(len(cseg)):
        for j in range(i+1,len(cseg)):
            count += 1
            nextone = info[i,j]
            print(i,j,nextone*100)
            score += nextone**2
    return score**.5*100
            
#     cAseg = cA[ny1:ny2,nx1:nx2]; cAsegmean = np.mean(cAseg); print('<cA> =',cAsegmean)
#     cBseg = cB[ny1:ny2,nx1:nx2]; cBsegmean = np.mean(cBseg); print('<cB> =',cBsegmean)
#     cCseg = cC[ny1:ny2,nx1:nx2]; cCsegmean = np.mean(cCseg); print('<cC> =',cCsegmean)
#     cDseg = cD[ny1:ny2,nx1:nx2]; cDsegmean = np.mean(cDseg); print('<cD> =',cDsegmean)
#     cArel = cAseg - cAsegmean; cArel = cArel/np.std(cArel)
#     cBrel = cBseg - cBsegmean; cBrel = cBrel/np.std(cBrel)
#     cCrel = cCseg - cCsegmean; cCrel = cCrel/np.std(cCrel)
#     cDrel = cDseg - cDsegmean; cDrel = cDrel/np.std(cDrel)
#     infomatrix = np.matrix(np.zeros(shape=(4,4)))
#     infomatrix[0,0] = np.mean(cArel*cArel)
#     infomatrix[0,1] = np.mean(cArel*cBrel)
#     infomatrix[0,2] = np.mean(cArel*cCrel)
#     infomatrix[0,3] = np.mean(cArel*cDrel)
#     infomatrix[1,1] = np.mean(cBrel*cBrel)
#     infomatrix[1,2] = np.mean(cBrel*cCrel)
#     infomatrix[1,3] = np.mean(cBrel*cDrel)
#     infomatrix[2,2] = np.mean(cCrel*cCrel)
#     infomatrix[2,3] = np.mean(cCrel*cDrel)
#     infomatrix[3,3] = np.mean(cDrel*cDrel)

def Weibull(Z2,sigma2W,etaW):
    # Getting the Weibull distribution
    term1 = etaW/(sigma2W)
    term2 = (Z2/sigma2W)**(etaW-1)
    term3 = np.exp(-(Z2/sigma2W)**etaW)
    rhoW = term1*term2*term3
    #rhoW = etaW/(sigma2W)*(Z2/sigma2W)**(etaW-1)*np.exp(-(Z2/sigma2W)**etaW)
    return rhoW

def logWeibull(Z2,sigma2W,etaW):
    temp1 = Weibull(Z2,sigma2W,etaW)
    temp2 = np.log(temp1)
    return temp2

def Gaussian(Z2,sigma2G):
    term1 = 1/(sigma2G)
    term2 = np.exp(-(Z2/sigma2G))
    rhoW = term1*term2
    return rhoW
    
def logGaussian(Z2,sigma2G):
    temp1 = Gaussian(Z2,sigma2G)
    temp2 = np.log(temp1)
    return temp2
