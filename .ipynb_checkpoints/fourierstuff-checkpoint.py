
# Does 2d Fourier transform (to/from) 

import numpy as np
def FT(grid,x,y):

    # Fourier transform it
    grid_FT = np.fft.fft2(grid)
    grid_FTshift = np.fft.fftshift(grid_FT)

    # And get the frequencies
    Ny, Nx = np.shape(grid)
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    kx = np.fft.fftfreq(Nx)*2*np.pi/dx
    ky = np.fft.fftfreq(Ny)*2*np.pi/dy
    kxshift = np.fft.fftshift(kx)
    kyshift = np.fft.fftshift(ky)
    return grid_FTshift, kxshift, kyshift

def IFT(sollast_FTshift_filtered):

    # Un-shift it
    sollast_FT_filtered = np.fft.ifftshift(sollast_FTshift_filtered)

    # Inverse Fourier transform
    sollast_FT_filtered_IFT = np.fft.ifft2(sollast_FT_filtered)

    return sollast_FT_filtered_IFT

def filterseg(kmax,xgrid,ygrid,sollast):
#     xgrid = xseggrid[i]
#     ygrid = yseggrid[i]
#     sollast = zseggrid[i]
    x = xgrid[0,:]; #print(x)
    y = ygrid[:,0]; #print(y)
    Ny, Nx = np.shape(sollast)

    # Fourier transform it
    sollast_offset = np.mean(sollast)
    sollast_atzero = sollast-sollast_offset
    #sollast_FTshift,kxshift,kyshift = fs.FT(sollast_atzero,x,y)
    sollast_FTshift,kxshift,kyshift = FT(sollast_atzero,x,y)

    # This makes gridded versions of the k-arrays
    kxshiftgrid,kyshiftgrid = np.meshgrid(kxshift,kyshift);

    # Find the absolute square
    sollast_FTshift_square = np.real(sollast_FTshift)**2 +  np.imag(sollast_FTshift)**2 

    # Make a copy of the Fourier representation
    sollast_FTshift_filtered = sollast_FTshift*1

    # Apply a mask
    for ix in range(Nx):
        for iy in range(Ny):
            ktest = np.sqrt(kxshiftgrid[iy,ix]**2+kyshiftgrid[iy,ix]**2)
            if(ktest>kmax):
                sollast_FTshift_filtered[iy,ix]=0

    # Inverse FT
    #sollast_FT_filtered_IFT = fs.IFT(sollast_FTshift_filtered)+sollast_offset
    sollast_FT_filtered_IFT = IFT(sollast_FTshift_filtered)+sollast_offset
    sollast_FT_filtered_IFT_real = np.real(sollast_FT_filtered_IFT)
    
    # Get out
    return (sollast_FT_filtered_IFT_real, kxshiftgrid, kyshiftgrid, sollast_FTshift_square)

from itertools import product
from numpy import empty, roll
def autocorrelate(xgrid,ygrid,x):
    """
    Compute the multidimensional autocorrelation of an nd array.
    input: an nd array of floats
    output: an nd array of autocorrelations
    This is from https://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
    """
    
    # grids
    xgrid_new = xgrid[1:-1,1:-1]
    xgrid_mean = np.mean(xgrid)
    xgrid_new = xgrid_new - xgrid_mean
    ygrid_new = ygrid[1:-1,1:-1]
    ygrid_mean = np.mean(ygrid)
    ygrid_new = ygrid_new - ygrid_mean

    # used for transposes
    t = roll(range(x.ndim), 1)

    # pairs of indexes
    # the first is for the autocorrelation array
    # the second is the shift
    ii = [list(enumerate(range(1, s - 1))) for s in x.shape]

    # initialize the resulting autocorrelation array
    acor = empty(shape=[len(s0) for s0 in ii])

    # iterate over all combinations of directional shifts
    for i in product(*ii):
        # extract the indexes for
        # the autocorrelation array 
        # and original array respectively
        i1, i2 = np.asarray(i).T

        x1 = x.copy()
        x2 = x.copy()

        for i0 in i2:
            # clip the unshifted array at the end
            x1 = x1[:-i0]
            # and the shifted array at the beginning
            x2 = x2[i0:]

            # prepare to do the same for 
            # the next axis
            x1 = x1.transpose(t)
            x2 = x2.transpose(t)

        # normalize shifted and unshifted arrays
        x1 -= x1.mean()
        x1 /= x1.std()
        x2 -= x2.mean()
        x2 /= x2.std()

        # compute the autocorrelation directly
        # from the definition
        acor[tuple(i1)] = (x1 * x2).mean()

    return xgrid_new, ygrid_new, acor

def autocovariance(xgrid,ygrid,Z,debug=False):
    # Trying to implement what Jawaid & Seewig, 2023, seem to be saying.
    Zmean = np.mean(Z); print('Z mean ', Zmean)
    Zstd = np.std(Z); print('Z std dev ',Zstd)
    Znew = Z-Zmean

    # Decide how far out to take the displacement
    ny,nx = np.shape(Znew)
    print('ny,nx ', ny,nx)
    fraction_to_show = 0.25
    
    # Get the autocovariance, crudely
    fullACVF = np.zeros((2*ny-1, 2*nx-1))
    ny_fullACVF,nx_fullACVF = np.shape(fullACVF)
    if debug: print('ny_full, nx_full ', ny_fullACVF,nx_fullACVF)
    for ix in range(nx_fullACVF):
        for iy in range(ny_fullACVF):
            
            iy1_m1 = ny-iy-1
            iy2_m1 = ny
            ix1_m1 = nx-ix-1
            ix2_m1 = nx
            
            if iy1_m1 < 0: 
                iy2_m1 += iy1_m1
                iy1_m1 = 0
            if ix1_m1 < 0: 
                ix2_m1 += ix1_m1
                ix1_m1 = 0
            
            matrix1 = Znew[iy1_m1:iy2_m1,ix1_m1:ix2_m1]
            
            iy1_m2 = 0
            iy2_m2 = iy+1
            ix1_m2 = 0
            ix2_m2 = ix+1
            
            if iy2_m2 > ny:
                iy1_m2 = iy2_m2 - ny
                iy2_m2 = ny
            if ix2_m2 > nx:
                ix1_m2 = ix2_m2 - nx
                ix2_m2 = nx
                
            matrix2 = Znew[iy1_m2:iy2_m2,ix1_m2:ix2_m2]
            
            ny_matrix1,nx_matrix1 = np.shape(matrix1)
            ny_matrix2,nx_matrix2 = np.shape(matrix2)
            if ny_matrix1==ny_matrix2 and nx_matrix1==nx_matrix2:
                E = np.mean(matrix1*matrix2)
            else:
                print('We have a problem at ...')
                print('iy,ix ', iy,ix)
                print('ny_matrix1,ny_matrix2', ny_matrix1,ny_matrix2)
                print('nx_matrix1,nx_matrix2', nx_matrix1,nx_matrix2)
                break
            
            if (ix<4 and iy<4 and debug):
                print('For iy,ix ', iy,ix, np.shape(matrix1), np.shape(matrix2))
                print(' m1 ', matrix1)
                print(' m2 ', matrix2)
                print(' E ', E)
                print('\n')
                
            fullACVF[iy,ix] = E

    # Pare it down to a more manageable size
    iy_center = int(ny_fullACVF/2)
    ix_center = int(nx_fullACVF/2)
    ny_to_show = int(ny_fullACVF*fraction_to_show)
    nx_to_show = int(nx_fullACVF*fraction_to_show)
    ACVF = fullACVF[\
            iy_center-ny_to_show:iy_center+ny_to_show, \
            ix_center-nx_to_show:ix_center+nx_to_show]
    ny_ACVF,nx_ACVF = np.shape(ACVF)
    print('Size before paring ', ny_fullACVF,nx_fullACVF)
    print('Size after paring ', ny_ACVF,nx_ACVF)
    
    # grids
    ygrid_new = ygrid[0:ny_ACVF,0:nx_ACVF]
    xgrid_new = xgrid[0:ny_ACVF,0:nx_ACVF] 
    
    ygrid_new_mean = np.mean(ygrid_new)
    xgrid_new_mean = np.mean(xgrid_new)
    
    ygrid_new = ygrid_new - ygrid_new_mean
    xgrid_new = xgrid_new - xgrid_new_mean

            
    return xgrid_new, ygrid_new, ACVF