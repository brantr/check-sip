import numpy as np
from astropy import wcs
from astropy.io import fits
from sip_transform import *
from convert_siaf_to_sip import *


#set CRPIX
CRPIX1  = np.double(1024.5)
CRPIX2  = np.double(1024.5)
print("CRPIX1, CRPIX2 = ",CRPIX1,CRPIX2)

fname_siaf_sca = ["NRCA1_FULL.ascii", "NRCA2_FULL.ascii", "NRCA3_FULL.ascii","NRCA4_FULL.ascii","NRCA5_FULL.ascii"]

siaf, Sci2IdlCoeffX, Sci2IdlCoeffY, nidlc = siaf_read_parameters(fname_siaf_sca[0])

siaf["x_sci_ref"] = 1024.5
siaf["y_sci_ref"] = 1024.5
siaf["x_det_ref"] = 1024.5
siaf["y_det_ref"] = 1024.5
print("XSciRef = ",siaf["x_sci_ref"])
print("YSciRef = ",siaf["y_sci_ref"])
print("V3SciXAng = ",siaf["v3_sci_x_angle"])
print("V3SciYAng = ",siaf["v3_sci_y_angle"])
print("V3IdlYAng = ",siaf["v3_idl_yang"])
print("VSciParity = ",siaf["v_idl_parity"])
print("SciIdlCoeffX[0,0] = ",Sci2IdlCoeffX[0,0])
print("SciIdlCoeffY[0,0] = ",Sci2IdlCoeffY[0,0])
print("SciIdlCoeffX[1,0] = ",Sci2IdlCoeffX[1,0])
print("SciIdlCoeffY[1,1] = ",Sci2IdlCoeffY[1,1])
print("SciIdlCoeffX[5,5] = ",Sci2IdlCoeffX[5,5])
print("SciIdlCoeffY[5,5] = ",Sci2IdlCoeffY[5,5])


u = np.double(1.0)
v = np.double(1.0)

#set the detector coordinates
xydet = np.zeros(2)
xydet[0] = u + CRPIX1
xydet[1] = v + CRPIX2

#get the science coordinates from the 
#detector coordinates
xysci = siaf_det_to_sci(xydet, siaf)
s = "For xydet=(%f,%f), xysci=(%f,%f)" % (xydet[0],xydet[1],xysci[0],xysci[1])
print(s)

#get the ideal coordinates from the
#science coordinates
xyidl = siaf_sci2idl(xysci, siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)
s = "For xysci=(%f, %f), xyidl = (%14.12e, %14.12e)" % (xysci[0],xysci[1],xyidl[0],xyidl[1])
print(s)

#get the sip coefficients from the Sci2IdlCoeffs

#Apq, Bpq = sip_from_siaf(siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)

exit() 

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





