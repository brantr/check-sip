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
print("SciIdlCoeffX[0,0] = ",Sci2IdlCoeffX[0,0])
print("SciIdlCoeffX[1,0] = ",Sci2IdlCoeffX[1,0])
print("SciIdlCoeffX[1,1] = ",Sci2IdlCoeffX[1,1])
print("SciIdlCoeffX[2,0] = ",Sci2IdlCoeffX[2,0])
print("SciIdlCoeffX[2,1] = ",Sci2IdlCoeffX[2,1])
print("SciIdlCoeffX[2,0] = ",Sci2IdlCoeffX[2,2])

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

CDELT1, CDELT2, Apq, Bpq, A_ORDER, B_ORDER = siaf_to_sip(siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)

verbose = False
if(verbose):
	for p in range(Apq.shape[0]):
		for q in range(Apq.shape[1]-p):
			s = "p %d q %d i+j %d j %d CX % 14.13e CY % 14.13e " % (p,q,p+q,q,Apq[p,q],Bpq[p,q])
			print(s)

	for p in range(Apq.shape[0]):
		for q in range(Apq.shape[1]-p):
			s = "p %d q %d Apq % 14.13e Bpq % 14.13e " % (p,q,Apq[p,q],Bpq[p,q])
			print(s)

#set CDij
CDij = ((CDELT1, 0), (0, CDELT2))

#compute SIP expansion
xy = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
s = "For u,v=(%f,%f), xy=(%14.12e, %14.12e)" % (u,v,xy[0],xy[1])
print(s)

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


#i 1 j 0 CX  3.1139e-02 CY  2.5754e-05 xi -3.11392706e-02 yi -2.57544769e-05
#i 1 j 1 CX  0.0000e+00 CY  3.1322e-02 xi -3.11392706e-02 yi  3.12965629e-02
#i 2 j 0 CX  2.6479e-09 CY  6.6304e-08 xi -3.11392679e-02 yi  3.12966292e-02
#i 2 j 1 CX -2.0728e-07 CY  3.9248e-08 xi -3.11390606e-02 yi  3.12965899e-02
#i 2 j 2 CX -3.4585e-08 CY -1.4301e-07 xi -3.11390952e-02 yi  3.12964469e-02
#i 3 j 0 CX  9.9814e-12 CY  5.6684e-13 xi -3.11390952e-02 yi  3.12964469e-02
#i 3 j 1 CX  1.9623e-12 CY  9.3897e-12 xi -3.11390952e-02 yi  3.12964469e-02
#i 3 j 2 CX  1.1383e-11 CY  1.5175e-12 xi -3.11390952e-02 yi  3.12964469e-02
#i 3 j 3 CX  9.9310e-13 CY  1.1021e-11 xi -3.11390952e-02 yi  3.12964469e-02
#i 4 j 0 CX -7.5499e-16 CY -1.2449e-16 xi -3.11390952e-02 yi  3.12964469e-02
#i 4 j 1 CX -8.4298e-16 CY -6.7492e-16 xi -3.11390952e-02 yi  3.12964469e-02
#i 4 j 2 CX -9.6680e-16 CY -1.0749e-15 xi -3.11390952e-02 yi  3.12964469e-02
#i 4 j 3 CX -8.3657e-16 CY -6.4470e-16 xi -3.11390952e-02 yi  3.12964469e-02
#i 4 j 4 CX -1.7431e-16 CY -8.8158e-16 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 0 CX  1.4128e-19 CY  3.7366e-21 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 1 CX  9.6091e-21 CY  1.5857e-19 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 2 CX  2.9035e-19 CY  1.6127e-20 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 3 CX  7.8866e-21 CY  2.9722e-19 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 4 CX  1.4538e-19 CY  1.1329e-20 xi -3.11390952e-02 yi  3.12964469e-02
#i 5 j 5 CX -4.3172e-22 CY  1.3700e-19 xi -3.11390952e-02 yi  3.12964469e-02

#p 0 q 0 i+j 0 j 0 CX -0.0000000000000e+00 CY  0.0000000000000e+00
#p 0 q 1 i+j 1 j 1 CX  0.0000000000000e+00 CY  3.1322317344598e-02
#p 0 q 2 i+j 2 j 2 CX -3.4585136135670e-08 CY -1.4301268526321e-07
#p 0 q 3 i+j 3 j 3 CX  9.9309668917683e-13 CY  1.1021330139259e-11
#p 0 q 4 i+j 4 j 4 CX -1.7431228696624e-16 CY -8.8158086621863e-16
#p 0 q 5 i+j 5 j 5 CX -4.3171561141436e-22 CY  1.3699613093905e-19
#p 1 q 0 i+j 1 j 0 CX  3.1139270573204e-02 CY  2.5754476873595e-05
#p 1 q 1 i+j 2 j 1 CX -2.0728141170157e-07 CY  3.9247807945874e-08
#p 1 q 2 i+j 3 j 2 CX  1.1382764595793e-11 CY  1.5174647869935e-12
#p 1 q 3 i+j 4 j 3 CX -8.3656634546455e-16 CY -6.4470479152582e-16
#p 1 q 4 i+j 5 j 4 CX  1.4538020584747e-19 CY  1.1329298384821e-20
#p 2 q 0 i+j 2 j 0 CX  2.6478689607189e-09 CY  6.6303854242098e-08
#p 2 q 1 i+j 3 j 1 CX  1.9623269772875e-12 CY  9.3896655770785e-12
#p 2 q 2 i+j 4 j 2 CX -9.6680295804739e-16 CY -1.0749273161663e-15
#p 2 q 3 i+j 5 j 3 CX  7.8866015622962e-21 CY  2.9721796374688e-19
#p 3 q 0 i+j 3 j 0 CX  9.9813648422923e-12 CY  5.6683700303698e-13
#p 3 q 1 i+j 4 j 1 CX -8.4297715210866e-16 CY -6.7491979672650e-16
#p 3 q 2 i+j 5 j 2 CX  2.9034833520774e-19 CY  1.6126926350745e-20
#p 4 q 0 i+j 4 j 0 CX -7.5498954275526e-16 CY -1.2448546551204e-16
#p 4 q 1 i+j 5 j 1 CX  9.6091055290891e-21 CY  1.5856757785470e-19
#p 5 q 0 i+j 5 j 0 CX  1.4128314056249e-19 CY  3.7365835042945e-21
