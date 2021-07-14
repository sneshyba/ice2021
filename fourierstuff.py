
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

