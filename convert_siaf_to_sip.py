import numpy as np

def siaf_read_parameters(siaf_file):
	fp = open(siaf_file,"r")
	fl = fp.readlines()
	fp.close()

	siaf = {"x_det_ref":0.0, "y_det_ref":0.0, "det_sci_yangle":0.0, \
	"det_sci_parity":0.0, "x_sci_ref":0.0, "y_sci_ref":0.0, \
	"x_sci_scale":0.0, "x_sci_scale":0.0, \
	"v3_sci_x_angle":0.0, "v3_sci_y_angle":0.0, "v3_idl_yang":0.0, \
	"v_idl_parity":0.0, "v2_ref":0.0, "v3_ref":0.0, "sci_to_idl_degree":0, \
	"nterms":0}

	siaf["x_det_ref"]         = np.double(fl[0])
	siaf["y_det_ref"]         = np.double(fl[1])
	siaf["det_sci_yangle"]    = np.double(fl[2])
	siaf["det_sci_parity"]    = np.double(fl[3])
	siaf["x_sci_ref"]         = np.double(fl[4])
	siaf["y_sci_ref"]         = np.double(fl[5])
	siaf["x_sci_scale"]       = np.double(fl[6])
	siaf["y_sci_scale"]       = np.double(fl[7])
	siaf["v3_sci_x_angle"]    = np.double(fl[8])
	siaf["v3_sci_y_angle"]    = np.double(fl[9])
	siaf["v3_idl_yang"]       = np.double(fl[10])
	siaf["v_idl_parity"]      = np.double(fl[11])
	siaf["v2_ref"]            = np.double(fl[12])
	siaf["v3_ref"]            = np.double(fl[13])
	siaf["sci_to_idl_degree"] = np.int(fl[14].split()[0])


	verbose = False
	if(verbose):
		print("x_sci_scale = ",siaf["x_sci_scale"])
		print("y_sci_scale = ",siaf["y_sci_scale"])
		print("sci_to_idl_degree = ",siaf["sci_to_idl_degree"])

	nx = siaf["sci_to_idl_degree"]
	nterms = nx*(nx+1)/2
	siaf["nterms"] = nterms

	Sci2IdlCoefX = np.zeros((nx,nx))
	Sci2IdlCoefY = np.zeros((nx,nx))

	#current line number in the file
	nl = 15

	#read in Sci2IdlCoefX
	for i in range(nx):
		for j in range(i+1):
			Sci2IdlCoefX[i,j] = np.double(fl[nl].split()[2])
			nl+=1	#iterate line number
	
	#skip Sci2IdlY header
	nl+=1

	#read in Sci2IdlCoefY
	for i in range(nx):
		for j in range(i+1):
			Sci2IdlCoefY[i,j] = np.double(fl[nl].split()[2])
			nl+=1	#iterate line number

	return siaf, Sci2IdlCoefX, Sci2IdlCoefY, siaf["sci_to_idl_degree"]


def siaf_sci2idl(xysci, siaf, Sci2IdlCoefX, Sci2IdlCoefY):

	verbose = False

	xyidl = np.zeros(2)
	degree = siaf["sci_to_idl_degree"]
	dxsci = xysci[0] - siaf['x_sci_ref']
	dysci = xysci[1] - siaf['y_sci_ref']

	if(verbose):
		print("dxsci, dysci = ",dxsci,dysci)

	if(verbose):
		print("degree = ",degree)
	for i in range(1,degree):
		for j in range(i+1):
			xyidl[0] += Sci2IdlCoefX[i,j]*(dxsci**(i-j))*(dysci**j)
			xyidl[1] += Sci2IdlCoefY[i,j]*(dxsci**(i-j))*(dysci**j)
			if(verbose):
				s = "i %d j %d CX % 5.4e CY % 5.4e xi % 9.8e yi % 9.8e" % (i,j,Sci2IdlCoefX[i,j],Sci2IdlCoefY[i,j],xyidl[0],xyidl[1])
				print(s)

	return xyidl

def siaf_det_to_sci(xydet,siaf):
	#ported from Christopher Willmer
	#Convert from detector to "science" coordinates. Uses expressions
	#in JWST-STScI-001550, SM-12, section 4.1

	xysci = np.zeros(2)

	deg_to_rad = np.pi/180.0

	angle = siaf["det_sci_yangle"] * deg_to_rad
	cosa  = np.cos(angle)
	sina  = np.sin(angle)

	#get siaf info
	x_det_ref      = siaf["x_det_ref"]
	y_det_ref      = siaf["y_det_ref"]
	x_sci_ref      = siaf["x_sci_ref"]
	y_sci_ref      = siaf["y_sci_ref"]
	det_sci_parity = siaf["det_sci_parity"]

	#from section 4.1 of JWST-STScI-001550
	t1       = (xydet[0]- x_det_ref)*cosa + (xydet[1]-y_det_ref)*sina
	xysci[0] = x_sci_ref + det_sci_parity * t1

	t2 = -1.*(xydet[0]-x_det_ref)*sina + (xydet[1]-y_det_ref)*cosa
	xysci[1] = y_sci_ref + t2

	return xysci


def siaf_to_sip(siaf,Sci2IdlCoefX, Sci2IdlCoefY):

	#convert siaf expansion coefficients
	#to SIP expansion coefficients

	#*************************************************************************
	# NC11 = DetSciParity * cos(DetSciYAngle)
	# NC12 = DetSciParity * sin(DetSciYAngle) = 0
	# NC21 = -sin(DetSciYAngle) = 0
	# NC22 = cos(DetSciYAngle)
	# A_i_j = CX_(i+j)_j * (PC11)**(i+1) * (PC22)**j
	# B_i_j = CY_(i+j)_j * (PC11)**i     * (PC22)**(j+1)
	# CDELT1 = Sci2IdlX10
	# CDELT2 = Sci2IdlY11 
	#*************************************************************************
	nx  = siaf["sci_to_idl_degree"]
	Apq = np.zeros((nx,nx))
	Bpq = np.zeros((nx,nx))

	deg_to_rad = np.pi/180.0
	angle = siaf["det_sci_yangle"] * deg_to_rad

	NC11 = siaf["det_sci_parity"] * np.cos(angle)
	NC22 = np.cos(angle)

	NC = np.zeros((2,2))
	NC[0,0] = NC11
	NC[1,1] = NC22

	print("NC11, NC22 = ",NC11,NC22)

	CDELT1 = np.abs(Sci2IdlCoefX[1,0])
	CDELT2 = np.abs(Sci2IdlCoefY[1,1])

	#print("CDELT1, CDELT2 = ",CDELT1, CDELT2)

	degree = nx
	for i in range(0,degree):
		for j in range(0,degree-i):
			#print(i,j,nx)
			Apq[i,j] = (Sci2IdlCoefX[i+j,j]/CDELT1) *(NC11**(i+1))*(NC22**j)
			Bpq[i,j] = (Sci2IdlCoefY[i+j,j]/CDELT2) *(NC11**i)*(NC22**(j+1))

	#remove the constant and linear terms
	#in sip, the (xx,yy) linear terms are modeled by 
	#the rotation matrix, while the (yy,xx) linear
	#terms are "bonus" and may be non zero
	#note -- will leave extra linear term Bpq[1,0]
	#print("Apq 00", Apq[0,0], Sci2IdlCoefX[0,0])
	#print("Apq 10", Apq[1,0], Sci2IdlCoefX[1,0])
	#print("Bpq 00", Bpq[0,0], Sci2IdlCoefX[0,0])
	#print("Bpq 01", Bpq[0,1], Sci2IdlCoefX[1,0])

	Apq[0,0] = 0
	Apq[1,0] = 0
	Bpq[0,0] = 0
	Bpq[0,1] = 0


	#store the expansion order
	A_ORDER = degree
	B_ORDER = degree

	return CDELT1, CDELT2, NC, Apq, Bpq, A_ORDER, B_ORDER
