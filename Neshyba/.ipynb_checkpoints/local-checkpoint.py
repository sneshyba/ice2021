import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import imagestuff as ims
import retrievestuff as rs
import gradstuff as gds
import copy
import f90nml

def response_function(navec,nbvec,ndvec,cA,cB,cC,cD,nx1list,nx2list,ny1list,ny2list,imageroot,graphics=False, verbose=False):        
        ABCDangle_deg = 16
        theta = ABCDangle_deg*np.pi/180
        nboxes = len(nx1list)
        if verbose: print ("nboxes =", nboxes)
        
        # For each detector, get s-values
        sAa, sBa, sCa, sDa = ims.mygets(navec,theta)
        sAb, sBb, sCb, sDb = ims.mygets(nbvec,theta)
        sAd, sBd, sCd, sDd = ims.mygets(ndvec,theta)
        slistA = np.squeeze(np.array([sAa, sAb, sAd]))
        slistB = np.squeeze(np.array([sBa, sBb, sBd]))
        slistC = np.squeeze(np.array([sCa, sCb, sCd]))
        slistD = np.squeeze(np.array([sDa, sDb, sDd]))
        
        # Extract the observed intensities
        cA_obs = []
        cB_obs = []
        cC_obs = []
        cD_obs = []
        for isegment in range(nboxes):
            nx1=nx1list[isegment]; nx2=nx2list[isegment] 
            ny1=ny1list[isegment]; ny2=ny2list[isegment]
            cA_obs.append(np.mean(cA[ny1:ny2,nx1:nx2].astype('float')))
            cB_obs.append(np.mean(cB[ny1:ny2,nx1:nx2].astype('float')))
            cC_obs.append(np.mean(cC[ny1:ny2,nx1:nx2].astype('float')))
            cD_obs.append(np.mean(cD[ny1:ny2,nx1:nx2].astype('float')))
        
        # Defines the range of the s-parameter
        srange = [-.6,.6]

        # See what the A-D detector parameters look like graphically
        if graphics:
            plt.figure()
            markersize = 10
            plt.plot(slistA,cA_obs,'or',markersize=15)
            plt.plot(slistB,cB_obs,'<b',markersize=15)
            plt.plot(slistC,cC_obs,'sy',markersize=15)
            plt.plot(slistD,cD_obs,'>g',markersize=15)
            plt.legend(['A', 'B', 'C', 'D'],loc='upper left')
            plt.plot(slistA[0],cA_obs[0],'k*')
            plt.plot(slistA[1],cA_obs[1],'kx')
            plt.plot(slistA[2],cA_obs[2],'k+')
            plt.plot(slistB[0],cB_obs[0],'k*')
            plt.plot(slistB[1],cB_obs[1],'kx')
            plt.plot(slistB[2],cB_obs[2],'k+')
            plt.plot(slistC[0],cC_obs[0],'k*')
            plt.plot(slistC[1],cC_obs[1],'kx')
            plt.plot(slistC[2],cC_obs[2],'k+')
            plt.plot(slistD[0],cD_obs[0],'k*')
            plt.plot(slistD[1],cD_obs[1],'kx')
            plt.plot(slistD[2],cD_obs[2],'k+')
            plt.grid()
            plt.xlim(srange)
            plt.xlabel('$s$')
            plt.ylabel('$c_{obs}$')
        
        # Fitting
        maxorder = 1
        order = min(len(slistA)-1,maxorder)
        pA = np.polyfit(slistA,cA_obs,order)
        if verbose: print('pA =', pA[0], ',', pA[1])
        pB = np.polyfit(slistB,cB_obs,order)
        if verbose: print('pB =', pB[0], ',', pB[1])
        pC = np.polyfit(slistC,cC_obs,order)
        if verbose: print('pC =', pC[0], ',', pC[1])
        pD = np.polyfit(slistD,cD_obs,order)
        if verbose: print('pD =', pD[0], ',', pD[1])
        s_theory = np.linspace(srange[0],srange[1])
        cA_theory = np.polyval(pA,s_theory)
        cB_theory = np.polyval(pB,s_theory)
        cC_theory = np.polyval(pC,s_theory)
        cD_theory = np.polyval(pD,s_theory)
        plt.plot(s_theory,cA_theory,'-r',linewidth=2)
        plt.plot(s_theory,cB_theory,'--b',linewidth=2)
        plt.plot(s_theory,cC_theory,'-.y',linewidth=2)
        plt.plot(s_theory,cD_theory,':g',linewidth=2)
        plt.title(imageroot)
        
        # Error in the fit
        errorA = (np.polyval(pA,slistA)-cA_obs)**2
        errorB = (np.polyval(pB,slistB)-cB_obs)**2
        errorC = (np.polyval(pC,slistC)-cC_obs)**2
        errorD = (np.polyval(pD,slistD)-cD_obs)**2
        error_tot = errorA.sum()+errorB.sum()+errorC.sum()+errorD.sum()

        return error_tot, pA, pB, pC, pD

def get_nvecs(angleManager, draw, cvecdir, boxa, boxb, nx1list, ny1list, verbose=False):
        avec = angleManager.avec
        bvec = angleManager.bvec
        cvec = angleManager.cvec
        xorigin = angleManager.xorigin
        yorigin = angleManager.yorigin
        
        # Draw the crystal orientation vectors
        filla = 150
        fillb = 70
        fillc = 0
        scale = 150
        linewidth = 10
        linea = [xorigin,yorigin,xorigin+avec[0]*scale,yorigin+avec[1]*scale]
        lineb = [xorigin,yorigin,xorigin+bvec[0]*scale,yorigin+bvec[1]*scale]
        linec = [xorigin,yorigin,xorigin+cvec[0]*scale,yorigin+cvec[1]*scale]
        draw.line(linea, fill=filla,width=linewidth)
        draw.line(lineb, fill=fillb,width=linewidth)
        draw.line(linec, fill=fillc,width=linewidth)
        
        # The basal normal
        # cvecdir = 'pointingup'
        cvecdir = 'pointingdown'
        if cvecdir== 'pointingup':
            ndvec = cvec
            if verbose: 
                print ('unit normal basal facet =\n',ndvec)
                print ('this likely means use +28 for avec pyramidal and -28 for bvec pyramidal')
        elif cvecdir== 'pointingdown':
            ndvec = -cvec; 
            if verbose:
                print ('unit normal basal facet =\n',ndvec)
                print ('this likely means use -28 for avec pyramidal and +28 for bvec pyramidal')
        else: 
            print('WARNING! SPECIFY C VEC DIRECTION')
        dboxcenterx = nx1list[2]
        dboxcentery = ny1list[2]
        scaledcenterxterm = (dboxcenterx+ndvec[0]*scale).item()
        scaledcenteryterm = (dboxcentery+ndvec[1]*scale).item()
        lined_disp = list(np.squeeze([dboxcenterx,dboxcentery,scaledcenterxterm,scaledcenteryterm]).astype(int))
        draw.line(lined_disp, fill=fillc,width=linewidth)
        
        # Normal associated with vector a
        #boxa= 'prismatic'
        boxa= 'pyramidal'
        navec = np.matrix(np.cross(cvec.T,avec.T)).T
        if verbose:
            print ('unit normal a-prismatic facet =\n',navec)
        aboxcenterx = nx1list[0]
        aboxcentery = ny1list[0]
        if boxa== 'prismatic':
            scaledcenterxterm = (aboxcenterx+navec[0]*scale).item()
            scaledcenteryterm = (aboxcentery+navec[1]*scale).item()
            linea_disp = list(np.squeeze([aboxcenterx,aboxcentery,scaledcenterxterm,scaledcenteryterm]).astype(int))
            draw.line(linea_disp, fill=filla,width=linewidth)
        elif boxa== 'pyramidal':
            Rot28 = ims.myrotation_matrix(avec,-28) #this could be +/- 28, check
            navec = Rot28*navec
            if verbose:
                print ('unit normal a-pyramidal facet =\n',navec)
            scaledcenterxterm = (aboxcenterx+navec[0]*scale).item()
            scaledcenteryterm = (aboxcentery+navec[1]*scale).item()
            linea_disp = list(np.squeeze([aboxcenterx,aboxcentery,scaledcenterxterm,scaledcenteryterm]).astype(int))
            draw.line(linea_disp, fill=filla,width=linewidth)
        else: 
            print('WARNING! SPECIFY BOX A FACET')
        
        # Normal associated with vector b
        #boxb='prismatic'
        boxb='pyramidal'
        nbvec = np.matrix(np.cross(bvec.T,cvec.T)).T
        if verbose:
            print ('unit normal b-prismatic facet =\n',nbvec)
        bboxcenterx = nx1list[1]
        bboxcentery = ny1list[1]
        if boxb=='prismatic':
            scaledbvecxterm = (bboxcenterx+nbvec[0]*scale).item()
            scaledbvecyterm = (bboxcentery+nbvec[1]*scale).item()
            lineb_disp = list(np.squeeze([bboxcenterx,bboxcentery,scaledbvecxterm,scaledbvecyterm]).astype(int))
            draw.line(lineb_disp, fill=fillb,width=linewidth)
        elif boxb=='pyramidal':
            Rot28 = ims.myrotation_matrix(bvec,28) #this could be +/- 28, check
            nbvec = Rot28*nbvec
            if verbose:
                print ('unit normal b-pyramidal facet =\n',nbvec)
            bboxcenterx = nx1list[1]
            bboxcentery = ny1list[1]
            scaledbvecxterm = (bboxcenterx+nbvec[0]*scale).item()
            scaledbvecyterm = (bboxcentery+nbvec[1]*scale).item()
            lineb_disp = list(np.squeeze([bboxcenterx,bboxcentery,scaledbvecxterm,scaledbvecyterm]).astype(int))
            draw.line(lineb_disp, fill=fillb,width=linewidth)
        else:
            print('WARNING! SPECIFY BOX B FACET')

        return navec, nbvec, ndvec

def get_box_lists(Boxesfile,draw):
        
        Boxes=f90nml.read(Boxesfile) # Reads the locations of the boxes
        nx1list=Boxes['Boxes']['nx1list']
        ny1list=Boxes['Boxes']['ny1list']
        labellist=Boxes['Boxes']['labellist']
        boxsize=Boxes['Boxes']['boxsize']; #print (boxsize)
        nboxes = len(nx1list); print ("nboxes =", nboxes)
        nx2list = np.array(nx1list)+boxsize; #print(nx2list)
        ny2list = np.array(ny1list)+boxsize; #print(ny2list)
        
        # Show the boxes
        for i in range(nboxes):
            nx1 = nx1list[i]
            nx2 = nx2list[i]
            ny1 = ny1list[i]
            ny2 = ny2list[i]
            ims.myrectanglelabel(draw,(nx1,ny1),(nx2,ny2),labellist[i])

        return nx1list,nx2list,ny1list,ny2list

def retrieve_segments(\
    pA, pB, pC, pD, cA, cB, cC, cD, nx1list, nx2list, ny1list, ny2list, imageroot, dx=1, dy=1, \
    overlapping=True, nptsx=103, nptsy=101, nstride=3):
    
    # Number of segments
    nsegments = len(nx1list); print('nsegments ', nsegments)
    
    # Set up sub-grids in case of multiple segments
    nyxgrid = []
    for i in range(nsegments):
        nyxgridi = [ (y, x) for y in range(ny1list[i], ny2list[i]+1) for x in range(nx1list[i], nx2list[i]+1) ]
        nyxgrid.append(nyxgridi)
    
    # Set up a grid of surface normal vectors and the backscatter response on them
    nxmid = int(nptsx/2); #print nxmid
    nymid = int(nptsy/2); #print nymid
    nmax = 5
    nxi = np.linspace(-nmax,nmax,nptsx); dnx = nxi[1]-nxi[0]
    nyi = np.linspace(-nmax,nmax,nptsy); dny = nyi[1]-nyi[0]
    nxigrid,nyigrid = np.meshgrid(nxi,nyi)
    theta = 15*np.pi/180
    sA = (-nxigrid*np.sin(theta)+np.cos(theta)-1)/(1+nxigrid**2+nyigrid**2)**.5
    sB = (-nyigrid*np.sin(theta)+np.cos(theta)-1)/(1+nxigrid**2+nyigrid**2)**.5
    sC = (+nxigrid*np.sin(theta)+np.cos(theta)-1)/(1+nxigrid**2+nyigrid**2)**.5
    sD = (+nyigrid*np.sin(theta)+np.cos(theta)-1)/(1+nxigrid**2+nyigrid**2)**.5
    
    # Set up the grids     
    BSgridA = np.polyval(pA,sA)
    BSgridB = np.polyval(pB,sB)
    BSgridC = np.polyval(pC,sC)
    BSgridD = np.polyval(pD,sD)

    # Generating the response function for each detector
    BSgridN = [BSgridA, BSgridB, BSgridC, BSgridD]
    BSgridL = ['A', 'B', 'C', 'D']
    BSmax = 150 # this for display purposes

    # Set up interpolators for detector responses
    Arule, Brule, Crule, Drule, \
    KAxrule, KAyrule, KBxrule, KByrule, KCxrule, KCyrule, KDxrule, KDyrule =\
    rs.setupdetectorresponse3(BSgridA, BSgridB, BSgridC, BSgridD, nxi, nyi, dnx, dny)

    # Create a blank slate
    solution = np.zeros(cA.shape)

    # Generic retrieval parameters
    maxiter = 5
    tolerance = 10

    # Define the variance in the observations (BS units^2)
    noiseamp = 25.0
    print('Std deviation in input signal is', noiseamp**.5)
    
    # Define parameters determining the variance in the a priori (microns^2)
    apriorivar0 = 225.0
    
    print('Std deviation in a priori is', apriorivar0**.5)
    
    # Create the initial a priori variance
    apriorivar = np.ones(cA.shape)*apriorivar0
    
    # Create the initial a priori set
    Normal_list = []
    aprioriset = np.zeros(cA.shape)

    # Loop to retrieve each segment
    for isegment in range(nsegments):
        
        # Choose the particular location of the dataset to analyze
        nx1=nx1list[isegment]; nx2=nx2list[isegment]; nx = nx2-nx1+1
        ny1=ny1list[isegment]; ny2=ny2list[isegment]; ny = ny2-ny1+1
    
        # Construct gradients
        Ny_unscaled, Nx_unscaled = gds.makeNxNy(ny,nx)
        Ny = Ny_unscaled/dy
        Nx = -Nx_unscaled/dx #fixing x inversion
        
        # Number of observations
        nobs = (nx-1)*(ny-1)*4
    
        # Number of desired points (heights)
        nzpts = ny*nx-1
        
        # Extract the a priori variance
        vartemp = apriorivar[ny1:ny2+1,nx1:nx2+1]
        #vartemp_long = np.reshape(vartemp,nzpts+1,0) # This appears to have been deprecated
        vartemp_long = vartemp.flatten(order='C')
        Sa = np.diag(vartemp_long[:-1]); #print "apriorivar", shape(Sa)
        
        # Extract the starting z
        settemp = solution[ny1:ny2+1,nx1:nx2+1]
        #settemp_long = np.reshape(settemp,nzpts+1,0)
        settemp_long = settemp.flatten(order='C')
        settemp_longminus1 = settemp_long[:-1]
        z_start = np.matrix(settemp_longminus1).T; #print "aprioriset", shape(z_start)
        z_start = z_start*0.0; #print "aprioriset", shape(z_start)
    
        # Construct the variance in observation + model
        Se = np.matrix(np.eye(nobs))*noiseamp # Variance in observation + model (c)
    
        # Do the retrieval
        print('')
        print("Segment:", isegment, '(', isegment+1, "of", nsegments, ')')
        print("for", nx1, ny1)
        z_retrieved = rs.retrievesegment(\
            nx1,ny1,nx2,ny2,cA,cB,cC,cD,\
            Sa,Se,z_start,maxiter,tolerance,\
            Nx,Ny,\
            Arule, Brule, Crule, Drule,\
            KAxrule, KAyrule, \
            KBxrule, KByrule, \
            KCxrule, KCyrule, \
            KDxrule, KDyrule)
        
        if isegment == 0:
            solution[ny1:ny2+1,nx1:nx2+1] = copy.copy(z_retrieved)
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.invert_yaxis() # invert y axis (this fixes the right-hand-oriented vs left-hand-oriented system)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            ax.set_title(imageroot)
            ax.view_init(azim=-120,elev=22)

        else:
            nextsolution = np.zeros(cA.shape)
            nextsolution[ny1:ny2+1,nx1:nx2+1] = copy.copy(z_retrieved)
            overlap = []
            for i in range(isegment):
                nextoverlap = list( set(nyxgrid[i])&set(nyxgrid[isegment]) )
                overlap = overlap + nextoverlap
                Noverlap = len(overlap); 
            print("Noverlap =", Noverlap)
            if overlapping:
                diff = 0.0
                for j in range(Noverlap):
                    diff += nextsolution[overlap[j]] - solution[overlap[j]]
                diffavg = diff/Noverlap
                z_retrieved -= diffavg
            solution[ny1:ny2+1,nx1:nx2+1] = copy.copy(z_retrieved)

        # Plotting
        surf_y = np.arange(ny1list[isegment],ny2list[isegment]+1,1)*dy
        surf_x = np.arange(nx1list[isegment],nx2list[isegment]+1,1)*dx
        surf_xgrid, surf_ygrid = np.meshgrid(surf_x,surf_y)
        ax.plot_surface(surf_xgrid, surf_ygrid, z_retrieved, rstride=nstride ,cstride=nstride)

        # Fitting
        A = surf_xgrid.reshape(nzpts+1)
        B = surf_ygrid.reshape(nzpts+1)
        C = z_retrieved.reshape(nzpts+1)
        C = np.array(C).reshape(-1)
        m = ims.polyfit2d(A,B,C,linear=True,order=1)

        # Find the normal vectors to each best-fit plane segment
        Normal = np.matrix([-m[2],-m[1],1])
        Normal /= np.linalg.norm(Normal)
        Normal_list.append(Normal)

    # Report angle between each pair of planes
    print('Angles between segments')
    for i in range(0,nsegments):
        for j in range(i+1,nsegments):
            dotprod = np.ndarray.item(np.dot(Normal_list[i],Normal_list[j].T))
            print(i,j,np.arccos(dotprod)/np.pi*180, ' degrees')

    # Package up the whole surface
    nx1tot = min(nx1list)
    nx2tot = max(nx2list)
    ny1tot = min(ny1list)
    ny2tot = max(ny2list)
    nxtot = nx2tot-nx1tot; print (nxtot)
    nytot = ny2tot-ny1tot; print (nytot)
    ymaxtot = (nytot-1)*dy; xmaxtot = (nxtot-1)*dx
    
    surf_ytot = np.linspace(0,ymaxtot,nytot); #print surf_ytot[1]-surf_ytot[0]; 
    surf_xtot = np.linspace(0,xmaxtot,nxtot); #print surf_xtot[1]-surf_xtot[0]; 
    surf_xgridtot, surf_ygridtot = np.meshgrid(surf_xtot,surf_ytot)
    settemp = solution[ny1tot:ny2tot,nx1tot:nx2tot]


    # Return
    return [surf_xgridtot, surf_ygridtot, settemp], [sA, sB, sC, sD]
