import numpy as np
from astropy import wcs
from astropy.io import fits
from sip_transform import *

fname = "gal_test_nn_F200W_481_001.slp.fits"

hdu = fits.open(fname)

header = hdu[0].header
del header['CD3_3']
del header['CRPIX3']
del header['CTYPE3']
del header['CUNIT3']
del header['CRVAL3']
del header['V3I_YANG']
del header['VPARITY']
del header['NAXIS3']

header['NAXIS'] = 2
header['WCSAXES'] = 2
#print(header['A_ORDER'])


print('A02 = ',header['A_0_2'])



#CDi_j = CDij_sip(header)
#print(CDi_j[0,0])
#print(CDi_j[0,1])
#print(CDi_j[1,0])
#print(CDi_j[1,1])

header['CRVAL1'] = 53.181671732945411 
header['CRVAL2'] = -27.823842206235913 
CRPIX1 = header['CRPIX1']
CRPIX2 = header['CRPIX2']
CRVAL1 = header['CRVAL1']
CRVAL2 = header['CRVAL2']
 

CRPIX = np.array([CRPIX1,CRPIX2],dtype=np.float64)
print("CRPIX = ",CRPIX1, CRPIX2)
print("CRVAL = ",CRVAL1, CRVAL2)


A_p_q = Apq_sip(header)
B_p_q = Bpq_sip(header)
print("A10, A01 = ",A_p_q[1,0],A_p_q[0,1])
print("B10, B01 = ",B_p_q[1,0],B_p_q[0,1])

xy_off = [header['CRVAL1'],header['CRVAL2']]
#xy = xy_from_sip(0.,0.,header) + xy_off
#print("For uv=[0,0], xy=",xy,xy_off)

xy_det = np.array([1,1],dtype=np.float64)
uv = xy_det.copy() - CRPIX
#print("First, xy_det = ",xy_det)
#print("First, uv = ",uv)
#xy = xy_from_sip(uv[0],uv[1], header) + xy_off
#print("For uv=[-1023.5,-1023.5], xy=",xy)

xy = xy_from_sip(1.,1.,header) + xy_off
s = "For uv=[%f,%f], xy= [%13.11f, %13.11f]" % (uv[0],uv[1],xy[0],xy[1])
print(s)

#print(header)

xy_det = [[1.0,1.0],[CRPIX1+1,CRPIX2+1]]
print("Check xy_det = ",xy_det)
pixcrd = np.array(xy_det,dtype=np.float64)
w = wcs.WCS(header)
world = w.all_pix2world(pixcrd, 0)
print(world)





