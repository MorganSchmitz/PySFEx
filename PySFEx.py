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
    try:
        deg = PSF_model.header['POLDEG1']
    except KeyError:
        return PSF_basis[0,:,:]
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
        
        
def create_seximage(PSF_stack, fov, filename='im.fits', size=[10000,10000], noise_level=1e-12):
    """Create a dummy astronomical .fits image, from stack of theoretical PSFs, to be fed to SExtractor.

        Parameters
        ----------
        PSF_stack : np.ndarray
            Stack of PSFs. Each row should be a single PSF, either flattened or as a square array.

        fov : np.ndarray
            Positions of each PSF in the field of view.

        filename : string, optional
            Name of the .fits output file. Default saves to 'im.fits' in the current working directory.
            
        size : array-like, optional
            Size of the dummy image, in pixels per pixels.
            
        noise_level : float, optional
            Standard deviation of the Gaussian background noise.

        Returns
        -------
        im : ndarray
            Numpy array version of the dummy image..
    """
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
    # change shape if PSF stack is flattened
    if len(shap) == 2:
        PSFs = PSF_stack.reshape(shap[0], int(shap[1]**.5), int(shap[1]**.5))
    else:
        PSFs = cp.copy(PSF_stack)
    half_size = int(PSFs.shape[1]/2)
    for psf, pos in zip(PSFs,pixel_pos):
        im[size[0]-pos[1]-half_size:size[0]-pos[1]+half_size+1, 
           pos[0]-half_size:pos[0]+half_size+1] = psf
    hdu = fits.PrimaryHDU(im)
    hdu.writeto(filename, clobber=True)
    return im
    
    
def convert_coords(catname,fov_train,fov_test):
    """Convert true field of view positions to SExtractor default convention. Assumes SExtractor 
    was ran on an headerless .fits image created by create_seximage.

        Parameters
        ----------
        catname : string
            Path to the .cat file (SExtractor output)
            
        fov_train : np.ndarray
            Positions in the field of view used to create dummy image with create_seximage.        
        
        fov_test : np.ndarray
            Positions in the field of view to be converted. Assumes same convention as fov_train.

        Returns
        -------
        newcoords : np.ndarray
            Positions in SExtractor/PSFEx convention, to be fed to interpsfex.
    """
    # read coordinates from SExtractor cat
    cat = fits.open(catname)
    objs = cat[2].data
    xs, ys = np.array([obj[-8] for obj in objs]), np.array([obj[-7] for obj in objs])
    # determine extreme coords
    xmin, ymin, xmax, ymax = np.min(xs), np.min(ys), np.max(xs), np.max(ys)
    fovmin = np.min(fov_train,axis=0)
    fovmax = np.max(fov_train,axis=0)
    fovspan = fovmax - fovmin
    # convert to PSFEx positions
    newx = np.array([xmin + (fov_pos - fovmin[0])/fovspan[0] * (xmax-xmin) 
                            for fov_pos in fov_test[:,0]])
    newy = np.array([(ymax*fovmax[1] - ymin*fovmin[1] - fov_pos*(ymax-ymin))/(fovspan[1])
                            for fov_pos in fov_test[:,1]])      
    return np.vstack((newx,newy)).T  
    
    
    
    
    
        
        
        