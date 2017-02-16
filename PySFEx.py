import numpy as np
from astropy.io import fits
import copy as cp

def interpsfex(dotpsfpath, x, y):
    """Use PSFEx generated model to perform spatial PSF interpolation.

        Parameters
        ----------
        dotpsfpath : string
            Path to .psf file (PSFEx output).

        x : float
            Position along the x-axis in the field of view.

        y : float
            Position along the y-axis in the field of view.

        Returns
        -------
        PSF : ndarray
            PSF imagette.
        """
    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]
    PSF_basis = np.array(PSF_model.data)[0][0]
    deg = PSF_model.header['POLDEG1']
    x_interp, x_scale = PSF_model.header['POLZERO1'], PSF_model.header['POLSCAL1']
    y_interp, y_scale = PSF_model.header['POLZERO2'], PSF_model.header['POLSCAL2']
    x, y = (x-x_interp)/x_scale, (y-y_interp)/y_scale
    
    # compute polynomial coefficients
    nb_monomials = int((deg+1)*(1+ float(deg)/2))
    coeff_vect = np.empty(nb_monomials)
    for i in xrange(0,deg+1):
        coeff_vect[i] = x**i
    count = deg+1
    for i in xrange(1,deg+1):
        for j in xrange(0,deg-i+1):
            coeff_vect[count] = (x**j)*(y**i)
            count +=1
    # compute interpolated PSF
    PSF = np.zeros((PSF_basis.shape[1], PSF_basis.shape[2]))
    for i in range(0,nb_monomials):
        PSF += coeff_vect[i]*PSF_basis[i,:,:]
    return PSF
        
        
def create_seximage(PSF_stack, fov, filename='im.fits', size=[10000,10000], noise_level=.001):
    # determine field of view
    posmin, posmax = np.min(fov,axis=0), np.max(fov,axis=0)
    span = posmax-posmin
    posmin -= .1*span
    posmax += .1*span
    span *= 1.2
    # initialize image with gaussian noise
    size = np.array([size[1],size[0]])
    im = np.random.normal(scale=noise_level,size=size)
    # determine central pixel positions
    pixel_pos = ((fov-posmin)/span * size).astype(int)
    shap = PSF_stack.shape
    if len(shap) == 2:
        PSFs = PSF_stack.reshape(shap[0], int(shap[1]**.5), int(shap[1]**.5))
    else:
        PSFs = cp.copy(PSF_stack)
    half_size = int(PSFs.shape[1]/2)
    for psf, pos in zip(PSFs,pixel_pos):
        im[size[0]-pos[1]-half_size:size[0]
-pos[1]+half_size+1, 
           pos[0]-half_size:pos[0]+half_size+1] = psf
    hdu = fits.PrimaryHDU(im)
    hdu.writeto(filename, clobber=True)
    return im
    
    
    
    
    
    
    
    
        
        
        