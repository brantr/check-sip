import numpy as np
from astropy import wcs
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
s = "Testing the conversion from SIAF to SIP for %s" % fname_siaf_sca[sca]
print(s)

#set CRPIX
CRPIX1  = np.double(1024.5)
CRPIX2  = np.double(1024.5)
print("CRPIX1, CRPIX2 = ",CRPIX1,CRPIX2)

# read in SIAF info and coefficients
siaf, Sci2IdlCoeffX, Sci2IdlCoeffY, nidlc = siaf_read_parameters(fname_siaf_sca[sca])

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
verbose = True
if(verbose):
	sip_print_coeffs(Apq,Bpq)

#set CDij
#here NC = ((v_idl_parity*cos(v3_idl_yang,0)),(0,cos(v3_idl_yang)))
CDij = np.zeros((2,2))
CDij[0,0] = CDELT1*NC[0,0]
CDij[1,1] = CDELT2*NC[1,1]



#Let's do some tests!
print("*************************")

#set a pixel location to compute the transform at
u = np.double(1.0)
v = np.double(1.0)

#set the detector coordinates
xydet = np.zeros(2)
xydet[0] = u + CRPIX1
xydet[1] = v + CRPIX2

#get the science coordinates from the 
#detector coordinates
xysci = siaf_det_to_sci(xydet, siaf)
verbose = False
if(verbose):
	s = "For xydet=(%f,%f), xysci=(%f,%f)" % (xydet[0],xydet[1],xysci[0],xysci[1])
	print(s)

#get the ideal coordinates from the
#science coordinates
xyidl = siaf_sci2idl(xysci, siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)
verbose = False
if(verbose):
	s = "For xysci=(%f, %f), xyidl = (%14.12e, %14.12e)" % (xysci[0],xysci[1],xyidl[0],xyidl[1])
	print(s)

#compute SIP expansion
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
verbose = False
if(verbose):
	s = "For u,v=(%f,%f), xysip=(%14.12e, %14.12e)" % (u,v,xysip[0],xysip[1])
	print(s)

s = "For u,v=(%f,%f), xyidl = (%14.12e, %14.12e)" % (u,v,xyidl[0],xyidl[1])
print(s)
s = "For u,v=(%f,%f), xysip = (%14.12e, %14.12e)" % (u,v,xysip[0],xysip[1])
print(s)
s = "Error in SIP vs. SIAF (for u=%f,v=%f):" % (u,v)
print(s)
s = "Delta x = % 15.14e " %(xysip[0]-xyidl[0])
print(s)
s = "Delta y = % 15.14e " %(xysip[1]-xyidl[1])
print(s)

#Do some more tests
u = np.double(100.0)
v = np.double(100.0)
xydet[0] = u + CRPIX1
xydet[1] = v + CRPIX2
xysci = siaf_det_to_sci(xydet, siaf)
xyidl = siaf_sci2idl(xysci, siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("*************************")
s = "For u,v=(%f,%f), xyidl = (%14.12e, %14.12e)" % (u,v,xyidl[0],xyidl[1])
print(s)
s = "For u,v=(%f,%f), xysip = (%14.12e, %14.12e)" % (u,v,xysip[0],xysip[1])
print(s)
#Compute Error
s = "Error in SIP vs. SIAF (for u=%f,v=%f):" % (u,v)
print(s)
s = "Delta x = % 15.14e " %(xysip[0]-xyidl[0])
print(s)
s = "Delta y = % 15.14e " %(xysip[1]-xyidl[1])
print(s)

#Do some more tests
u = np.double(1000.0)
v = np.double(1000.0)
xydet[0] = u + CRPIX1
xydet[1] = v + CRPIX2
xysci = siaf_det_to_sci(xydet, siaf)
xyidl = siaf_sci2idl(xysci, siaf, Sci2IdlCoeffX, Sci2IdlCoeffY)
xysip = xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER)
print("*************************")
s = "For u,v=(%f,%f), xyidl = (%14.12e, %14.12e)" % (u,v,xyidl[0],xyidl[1])
print(s)
s = "For u,v=(%f,%f), xysip = (%14.12e, %14.12e)" % (u,v,xysip[0],xysip[1])
print(s)
#Compute Error
s = "Error in SIP vs. SIAF (for u=%f,v=%f):" % (u,v)
print(s)
s = "Delta x = % 15.14e " %(xysip[0]-xyidl[0])
print(s)
s = "Delta y = % 15.14e " %(xysip[1]-xyidl[1])
print(s)


