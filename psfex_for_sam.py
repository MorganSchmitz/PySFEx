from astropy.io import fits
import subprocess
import os
from numpy import zeros, random


def randfilename(ext='.txt'):
    filename = str(random.randint(0, 1000000))+ext
    return filename

def map_to_cube(map,n1,n2):
    shap = map.shape
    ax1 = shap[0]/n1
    ax2 = shap[1]/n2
    nb_slice = ax1*ax2
    cube = zeros((n1,n2,nb_slice))
    for i in range(0,ax1):
        for j in range(0,ax2):
            cube[:,:,i+j*ax1] = map[n1*i:n1*(i+1),n2*j:n2*(j+1)]
    return cube


def poly_val(x,y,deg):
    nb_monomials = int((deg+1)*(1+ float(deg)/2))
    coeff_vect = zeros((nb_monomials,))

    for i in range(0,deg+1):
        coeff_vect[i] = x**i
    count = deg+1
    for i in range(1,deg+1):
        for j in range(0,deg-i+1):
            coeff_vect[count] = (x**j)*(y**i)
            count +=1
    return coeff_vect,nb_monomials

# ----------------------------------------------------------------------------------------------------------------------------------------------

def sextractor(in_filename,detect_thresh,PSF_SIZE,out_cat_name=None,out_check_name=None,check_type='BACKGROUND'):
    """
	|	Builts a catalogue of point sources from an astronomical image
	|
	|	Parameters
	|	----------
	|	in_filename : string
	|		Astronomical image file name (FITS format).
	|
    |   detect_thresh : scalar
	|		Pixels detection threshold : pixels with intensities below detect_thresh x the background standard deviation are discarded.
	|
    |	PSF_SIZE : array_like
	|		Point sources stamp sizes.
    |
	|	out_cat_name : string
	|		SEXTRACTOR catalogues's name (with no extension).
    |
    |	out_check_name : string
	|		Diagnostic output image file name.
    |
    |	check_type : string
	|		Diagnostic output image type.
    |
    |
	|	Return
	|	------
	|	out_cat_name : string
	|		SEXTRACTOR catalogues's name (with no extension).
    |
    |	out_check_name : string
	|		Diagnostic output image file name.
    |
	"""

    if out_cat_name is None:
       out_cat_name = randfilename(ext='')
    if out_check_name is None:
        out_check_name = randfilename(ext='.check')

    cat_fits = fits.getdata(in_filename)
    siz_cat = cat_fits.shape

    param = 'temp.param'
    subprocess.call(['cp']+['./def.param']+[param])
    f = open(param, 'a')
    f.write('VIGNET('+str(PSF_SIZE[0])+','+str(PSF_SIZE[1])+')')
    f.close()
    filt_filename = 'gauss_2.0_5x5.conv'

    com_sext = ['sex',in_filename,'-CATALOG_TYPE','FITS_LDAC','-CATALOG_NAME',out_cat_name+'.cat','-CHECKIMAGE_NAME',out_check_name\
                ,'-DETECT_MINAREA','4','-DETECT_THRESH',str(detect_thresh),'-FILTER','N','-FILTER_NAME',filt_filename,\
                '-DEBLEND_MINCONT','1','-WEIGHT_TYPE','NONE','-SATUR_LEVEL','1','-GAIN','0','-VERBOSE_TYPE','FULL','-PARAMETERS_NAME',param]
    subprocess.call(com_sext)
    #out = fits.getdata(out_check_name)
    os.remove(param)
    return out_cat_name,out_check_name

def psfextractor(PSF_SIZE,PSF_SAMPLING,cat_name,var_deg=3):
    """
	|	Estimates a PSF model from a catalogue of point sources outputed by SEXTRACTOR
	|
	|	Parameters
	|	----------
	|	PSF_SIZE : array_like
	|		Point sources stamp sizes.
	|
    |   PSF_SAMPLING : scalar
	|		Upsampling factor.
	|
	|	cat_name : string
	|		SEXTRACTOR catalogues's name (with no extension).
    |
    |   var_deg : scalar
    |       PSFEx's polynomial model's degree (might be decreased by the program).
	|
	|	Return
	|	------
	|	src : array_like
	|		"Eigen" PSFs.
    |
    |   deg : scalar
    |       Final PSFEx's polynomial model's degree
    |
    |   offsetx, scalex, offsety, scaley : scalars
    |       PSFEX's model normalization parameters
    |
    |   chi, moff, resi, samp, snap_cube, sub_sym, sub_moff : arrays_like
    |       Diagnostic outputs and parametric models.
	|
	|
	"""

    samp = 1./PSF_SAMPLING
    com_psfsex = ['psfex',cat_name+'.cat','-PSF_SAMPLING',str(samp),\
              '-PSF_SIZE',str(PSF_SIZE[0])+','+str(PSF_SIZE[1]),'-PSFVAR_DEGREES',\
              str(var_deg),'-SAMPLE_AUTOSELECT','N','-BASIS_TYPE','PIXEL','-SAMPLE_VARIABILITY',\
              '1','-SAMPLE_MAXELLIP','1','-SAMPLE_FWHMRANGE','0.00001,900','-VERBOSE_TYPE','FULL','-SAMPLE_MINSN','0.001']
    subprocess.call(com_psfsex)

    chi = fits.getdata('chi_'+cat_name+'.fits')
    #moff = fits.getdata('moffat_'+cat_name+'.fits')
    resi = fits.getdata('resi_'+cat_name+'.fits')
    samp = fits.getdata('samp_'+cat_name+'.fits')
    snap = fits.getdata('snap_'+cat_name+'.fits')
    snap_cube = map_to_cube(snap,PSF_SIZE[0],PSF_SIZE[1])
    #sub_sym = fits.getdata('subsym_'+cat_name+'.fits')
    #sub_moff = fits.getdata('submoffat_'+cat_name+'.fits')

    #proto_map = fits.getdata('proto_'+cat_name+'.fits')
    hdu = fits.open(cat_name+'.psf')
    offsetx = 0
    scalex = 1
    offsety = 0
    scaley = 1
    src = None
    '''if nb_comp==1:
        src = proto_map
    else:'''
    psf_mod = hdu[1].data
    nb_comp = psf_mod[0][0].shape[0]
    src = zeros((PSF_SIZE[0],PSF_SIZE[1],nb_comp))
    for i in range(0,nb_comp):
        src[:,:,i] = psf_mod[0][0][i,:,:]
    offsetx = hdu[1].header['POLZERO1']
    offsety = hdu[1].header['POLZERO2']
    scalex = hdu[1].header['POLSCAL1']
    scaley = hdu[1].header['POLSCAL2']
    deg = hdu[1].header['POLDEG1']
    hdu.close()

    return src,deg,offsetx,scalex,offsety,scaley,chi,resi,samp,snap_cube
#src,deg,offsetx,scalex,offsety,scaley,chi,moff,resi,samp,snap_cube,sub_sym,sub_moff


def psfex_est(posi,src,deg,offx=0,scalex=1,offy=0,scaley=1,dir = True):
    """
	|	Estimates a PSF at the given location in the FOV based on PSFEx polynomial model
	|
	|	Parameters
	|	----------
	|	posi : array_like
	|		Position in the FOV in pixels (in the model resolution) of the PSF to be estimated.
	|
    |	src : array_like
	|		"Eigen" PSFs.
    |
    |   deg : scalar
    |       Final PSFEx's polynomial model's degree
    |
    |   offsetx, scalex, offsety, scaley : scalars
    |       PSFEX's model normalization parameters
    |
    |   dir : boolean
    |       Coordinates convention variable
    |
	|	Return
	|	------
	|	est : array_like
	|		Estimated PSF.
    |
	|
	"""

    est = None
    coeff_vect = None
    nb_monomials = None
    if dir:
        coeff_vect,nb_monomials = poly_val((posi[0]-offx)/scalex,(posi[1]-offy)/scaley,deg)
    else:
        coeff_vect,nb_monomials = poly_val((posi[1]-offy)/scaley,(posi[0]-offx)/scalex,deg)
    shap = src.shape
    est = zeros((shap[0],shap[1]))
    if nb_monomials>1:
        for i in range(0,nb_monomials):
            est += coeff_vect[i]*src[:,:,i]
    else:
        est+=coeff_vect[0]*src
    return est
