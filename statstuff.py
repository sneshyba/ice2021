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

def quadfit_to_Weibull_parameters(b,c):
    etaW = np.sqrt(-b**2/c - b*np.sqrt(b**2 - 8*c)/c + 8)/2
    sigmaW = 2**(3/4)*((1 - etaW**2)/c)**(1/4)/2
    return sigmaW, etaW
