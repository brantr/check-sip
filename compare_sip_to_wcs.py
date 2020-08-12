import numpy as np
import astropy.wcs as wcs
from astropy.io import fits
from sip_transform import *
from convert_siaf_to_sip import *
import sys


# list of NIRCAM SCAs
fname_siaf_sca = ["NRCA1_FULL.ascii", "NRCA2_FULL.ascii", "NRCA3_FULL.ascii","NRCA4_FULL.ascii","NRCA5_FULL.ascii","NRCB1_FULL.ascii", "NRCB2_FULL.ascii", "NRCB3_FULL.ascii","NRCB4_FULL.ascii","NRCB5_FULL.ascii"]

# select at sca
sca = 0
if(len(sys.argv)>1): #optionally use a command line arg to specify a different sca
	sca = np.int(sys.argv[1])

#print the info about the test to the screen
s = "Testing the WCS vs. SIP transforms for %s" % fname_siaf_sca[sca]
print(s)

#set CRPIX
CRPIX1  = np.double(1024.5)
CRPIX2  = np.double(1024.5)

# read in SIAF info and coefficients
siaf, Sci2IdlCoeffX, Sci2IdlCoeffY, nidlc = siaf_read_parameters(fname_siaf_sca[sca])

#reset XSciRef, YSciRef, XDetRef, YDetRef
siaf["x_sci_ref"] = CRPIX1
siaf["y_sci_ref"] = CRPIX2
siaf["x_det_ref"] = CRPIX1
siaf["y_det_ref"] = CRPIX2

#reset XSciRef, YSciRef, XDetRef, YDetRef
siaf["x_sci_ref"] = CRPIX1
siaf["y_sci_ref"] = CRPIX2
siaf["x_det_ref"] = CRPIX1
siaf["y_det_ref"] = CRPIX2

#print some info if desired
verbose = False
if(verbose):
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


#get the sip coefficients from the Sci2IdlCoeffs
CDELT1, CDELT2, NC, Apq, Bpq, A_ORDER, B_ORDER = siaf_to_sip(siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)

#if desired, print the SIP coefficients to screen
verbose = False
if(verbose):
	sip_print_coeffs(Apq,Bpq)

#set CDij
#here NC = ((v_idl_parity*cos(v3_idl_yang,0)),(0,cos(v3_idl_yang)))
CDij = np.zeros((2,2))
CDij[0,0] = CDELT1*NC[0,0]
CDij[1,1] = CDELT2*NC[1,1]



print("*******************")
print("Comparing WCS vs. SIP")


fname = "gal_test_nn_F200W_481_001.slp.fits"

hdu = fits.open(fname)



header = hdu[0].header

#CORRECT THE SIP EXPANSION IN THE HEADER
header['WCSAXES'] = 2

verbose = False
for p in range(Apq.shape[0]):
	for q in range(Apq.shape[1]-p):
		A_term = "A_%d_%d" % (p,q)
		try:
			del header[A_term]
		except KeyError:
			s=''
		AP_term = "AP_%d_%d" % (p,q)
		try:
			del header[AP_term]
		except KeyError:
			s=''
		B_term = "B_%d_%d" % (p,q)
		try:
			del header[B_term]
		except KeyError:
			s=''
		BP_term = "BP_%d_%d" % (p,q)
		try:
			del header[BP_term]
		except KeyError:
			s = ''

		if(np.abs(Apq[p,q])>1.0e-20):
			s = "A_%d_%d" % (p,q)
			header[s] = Apq[p,q]
		if(np.abs(Bpq[p,q])>1.0e-20):
			s = "B_%d_%d" % (p,q)
			header[s] = Bpq[p,q]

		#if((p==1)&(q==0)):
		#	Bpq[p,q]=0
		if(verbose):
			s = "p %d q %d Apq % 14.13e Bpq % 14.13e " % (p,q,Apq[p,q],Bpq[p,q])
			print(s)

#remove other 3D properties that prevent WCS from using SIP
del header['CRPIX3']
del header['CRVAL3']
del header['CTYPE3']
del header['CUNIT3']
del header['CD3_3']
header['NAXIS'] = 2
del header['NAXIS3']

#add in CD in arcsec
CDij/=3600.
header['CD1_1'] = CDij[0,0]
header['CD1_2'] = CDij[0,1]
header['CD2_1'] = CDij[1,0]
header['CD2_2'] = CDij[1,1]


#set CRVAL
CRVAL1 = header['CRVAL1']
CRVAL2 = header['CRVAL2']

#print(header)
print("CRPIX1, CRPIX2 = ",CRPIX1,CRPIX2)
print("CRVAL1,2 = ",header['CRVAL1'], header['CRVAL2'])
print("CD1_1, CD1_2 = ",header['CD1_1'], header['CD1_2'])
print("CD2_1, CD2_2 = ",header['CD2_1'], header['CD2_2'])


#get the WCS from the header
W = wcs.WCS(header)


# do a test with u=1, v=1
u = 1.
v = 1.
ra, dec = W.all_pix2world([u+CRPIX1],[v+CRPIX2],1)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("********")
s = "u, v = %16.15e %16.15e" % (u,v)
print(s)
u_wcs, v_wcs = W.all_world2pix(ra,dec,1)
u_wcs -= CRPIX1
v_wcs -= CRPIX2
print("INVERSE ERROR ",u_wcs,v_wcs)
s = "INVERSE ERROR % 16.5e, %16.5e " % (u_wcs-u,v_wcs-v)
print(s)
s = "FROM WCS = %16.15e %16.15e" % (ra,dec)
print(s)
dec_sip = xysip[1]+CRVAL2
ra_sip = xysip[0]/np.cos(dec_sip * np.pi/180.)+CRVAL1
s = "FROM SIP = %16.15e %16.15e" % (ra_sip, dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [degrees]" % (ra[0]-ra_sip, dec[0]-dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [arcsec]" % ( (ra[0]-ra_sip)*3600., (dec[0]-dec_sip)*3600.0 )
print(s)

# do a test with u=10, v=10
u = 10.
v = 10.
ra, dec = W.all_pix2world([u+CRPIX1],[v+CRPIX2],1)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("********")
s = "u, v = %16.15e %16.15e" % (u,v)
print(s)
u_wcs, v_wcs = W.all_world2pix(ra,dec,1)
u_wcs -= CRPIX1
v_wcs -= CRPIX2
print("INVERSE ERROR ",u_wcs,v_wcs)
s = "INVERSE ERROR % 16.5e, %16.5e " % (u_wcs-u,v_wcs-v)
print(s)
s = "FROM WCS = %16.15e %16.15e" % (ra,dec)
print(s)
dec_sip = xysip[1]+CRVAL2
ra_sip = xysip[0]/np.cos(dec_sip * np.pi/180.)+CRVAL1
s = "FROM SIP = %16.15e %16.15e" % (ra_sip, dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [degrees]" % (ra[0]-ra_sip, dec[0]-dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [arcsec]" % ( (ra[0]-ra_sip)*3600., (dec[0]-dec_sip)*3600.0 )
print(s)

# do a test with u=100, v=100
u = 100.
v = 100.
ra, dec = W.all_pix2world([u+CRPIX1],[v+CRPIX2],1)
u_wcs, v_wcs = W.all_world2pix(ra,dec,1)
u_wcs -= CRPIX1
v_wcs -= CRPIX2
print("INVERSE ",u_wcs,v_wcs)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("********")
s = "u, v = %16.15e %16.15e" % (u,v)
print(s)
u_wcs, v_wcs = W.all_world2pix(ra,dec,1)
u_wcs -= CRPIX1
v_wcs -= CRPIX2
print("INVERSE ERROR ",u_wcs,v_wcs)
s = "INVERSE ERROR % 16.5e, %16.5e " % (u_wcs-u,v_wcs-v)
print(s)
s = "FROM WCS = %16.15e %16.15e" % (ra,dec)
print(s)
dec_sip = xysip[1]+CRVAL2
ra_sip = xysip[0]/np.cos(dec_sip * np.pi/180.)+CRVAL1
s = "FROM SIP = %16.15e %16.15e" % (ra_sip, dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [degrees]" % (ra[0]-ra_sip, dec[0]-dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [arcsec]" % ( (ra[0]-ra_sip)*3600., (dec[0]-dec_sip)*3600.0 )
print(s)

# do a test with u=1000, v=1000
u = 1000.
v = 1000.
ra, dec = W.all_pix2world([u+CRPIX1],[v+CRPIX2],1)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("********")
s = "u, v = %16.15e %16.15e" % (u,v)
print(s)
u_wcs, v_wcs = W.all_world2pix(ra,dec,1)
u_wcs -= CRPIX1
v_wcs -= CRPIX2
print("INVERSE ERROR ",u_wcs,v_wcs)
s = "INVERSE ERROR % 16.5e, %16.5e " % (u_wcs-u,v_wcs-v)
print(s)
s = "FROM WCS = %16.15e %16.15e" % (ra,dec)
print(s)
dec_sip = xysip[1]+CRVAL2
ra_sip = xysip[0]/np.cos(dec_sip * np.pi/180.)+CRVAL1
s = "FROM SIP = %16.15e %16.15e" % (ra_sip, dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [degrees]" % (ra[0]-ra_sip, dec[0]-dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [arcsec]" % ( (ra[0]-ra_sip)*3600., (dec[0]-dec_sip)*3600.0 )
print(s)

# do a test with u=10000, v=10000
u = 10000.
v = 10000.
ra, dec = W.all_pix2world([u+CRPIX1],[v+CRPIX2],1)
print("********")
s = "u, v = %16.15e %16.15e" % (u,v)
print(s)
s = "FROM WCS = %16.15e %16.15e" % (ra,dec)
print(s)
dec_sip = xysip[1]+CRVAL2
ra_sip = xysip[0]/np.cos(dec_sip * np.pi/180.)+CRVAL1
s = "FROM SIP = %16.15e %16.15e" % (ra_sip, dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [degrees]" % (ra[0]-ra_sip, dec[0]-dec_sip)
print(s)
s = "DIFFERENCE  = % 16.15e, % 16.15e [arcsec]" % ( (ra[0]-ra_sip)*3600., (dec[0]-dec_sip)*3600.0 )
print(s)

