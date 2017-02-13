import numpy as np
from astropy.io import fits


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
            PSF postage stamp.
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


def interpsfex_v2(dotpsfpath, x, y):
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
            PSF postage stamp.
        """

    # read PSF model and extract basis and polynomial degree and scale position
    PSF_model = fits.open(dotpsfpath)[1]
    PSF_basis = np.array(PSF_model.data)[0][0]
    deg = PSF_model.header['POLDEG1']
    x = (x - PSF_model.header['POLZERO1']) / PSF_model.header['POLSCAL1']
    y = (y - PSF_model.header['POLZERO2']) / PSF_model.header['POLSCAL2']

    # compute polynomial coefficients
    coeff_vect = [x ** i for i in xrange(deg + 1)]
    [[coeff_vect.append((x ** j) * (y ** i)) for j in xrange(deg - i + 1)]
     for i in xrange(1, deg + 1)]

    # compute interpolated PSF
    return np.sum((cv * pb for cv, pb in zip(coeff_vect, PSF_basis)))
