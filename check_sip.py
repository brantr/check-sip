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

A_p_q = Apq_sip(header)
B_p_q = Bpq_sip(header)

#print(A_p_q[1,0],A_p_q[0,1],A_p_q[5,0])
#print(B_p_q[1,0],B_p_q[0,1])

CDi_j = CDij_sip(header)
#print(CDi_j[0,0])
#print(CDi_j[0,1])
#print(CDi_j[1,0])
#print(CDi_j[1,1])
CRPIX1 = header['CRPIX1']
CRPIX2 = header['CRPIX2']
CRVAL1 = header['CRVAL1']
CRVAL2 = header['CRVAL2']

print("CRPIX = ",CRPIX1, CRPIX2)
print("CRVAL = ",CRVAL1, CRVAL2)

xy_off = [header['CRVAL1'],header['CRVAL2']]
xy = xy_from_sip(0.,0.,header) + xy_off
print("For uv=[0,0], xy=",xy,xy_off)
xy = xy_from_sip(1.,1.,header) + xy_off
print("For uv=[1,1], xy=",xy)

#print(header)

pixcrd = np.array([[CRPIX1-0.5,CRPIX2-0.5],[CRPIX1+0.5,CRPIX2+0.5]])
w = wcs.WCS(header)
world = w.all_pix2world(pixcrd, 0)[0]
print(world)





