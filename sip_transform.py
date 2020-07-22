import numpy as np
from astropy.io import fits

def xy_from_sip(u, v, header):

	#define the vector (u,v)
	uv = np.array([u,v],dtype=np.float64)
	print("uu,vv = ",uv)

	#get the A_p_q from the header
	A_p_q = Apq_sip(header)
	print(A_p_q)

	#get the B_p_q from the header
	B_p_q = Bpq_sip(header)

	#get the CDi_j from the header
	CDi_j = CDij_sip(header)
	print("CDij = ",CDi_j)

	#add the f(u,v) and g(u,v) functions
	fuv = fuv_sip(u,v,A_p_q)
	guv = guv_sip(u,v,B_p_q)
	print("xsum, ysum = ",fuv,guv)

	uv[0] += fuv # u + f(u,v)
	uv[1] += guv # v + g(u,v)

	#evaluate the transform
	print("xx = uu + xsum ",uv[0])
	print("yy = vv + ysum ",uv[1])

	xy = np.dot(CDi_j, uv)



	#xy[0]/=np.cos(xy[0]*np.pi/180.0)
	cosdec = np.cos(header['CRVAL2']*np.pi/np.double(180.0))
	print(cosdec)
	xy[0]/=cosdec
	print(type(xy[0]))
	s = "%13.10e %13.10e" % (xy[0],xy[1])
	print(s)

	#print((CDi_j[0,0]*uv[0] + CDi_j[0,1]*uv[1])/np.cos(header['CRVAL1']*np.pi/180.0) )
	#print(CDi_j[1,0]*uv[0] + CDi_j[1,1]*uv[1])

	return xy

def Apq_sip(header):

	A_ORDER = header['A_ORDER']

	A_p_q = np.zeros((A_ORDER+1,A_ORDER+1),dtype=np.float64)

	for p in range(A_ORDER+1):
		for q in range(A_ORDER+1):
			s = 'A_%d_%d' % (p,q)
			try:
				A_p_q[p,q] = header[s]
			except KeyError:
				s = ""
				#print("No "+s)

	return A_p_q

def Bpq_sip(header):

	B_ORDER = header['B_ORDER']

	B_p_q = np.zeros((B_ORDER+1,B_ORDER+1),dtype=np.float64)

	for p in range(B_ORDER+1):
		for q in range(B_ORDER+1):
			s = 'B_%d_%d' % (p,q)
			try:
				B_p_q[p,q] = header[s]
			except KeyError:
				s=""
				#print("No "+s)

	return B_p_q


def CDij_sip(header):
	CDi_j = np.zeros((2,2),dtype=np.float64)
	#CDi_j[0,0] = np.double(header['CD1_1'])
	#CDi_j[0,1] = np.double(header['CD1_2'])
	#CDi_j[1,0] = np.double(header['CD2_1'])
	#CDi_j[1,1] = np.double(header['CD2_2'])
	CDi_j[0,0] = np.double(-2.114214532367e-4)
	CDi_j[0,1] = np.double(1.801706494922e-4)
	CDi_j[1,0] = np.double(1.801706494922e-4)
	CDi_j[1,1] = np.double(2.114214532367e-4)
	#print(CDi_j)

	#print(type(CDi_j[0,0]))

	return CDi_j

def fuv_sip(u,v,A_p_q):

	fuv = np.double(0.0)

	print(A_p_q.shape[0],A_p_q.shape[1])

	for p in range(A_p_q.shape[0]):
		for q in range(A_p_q.shape[1]-p):
			fuv += A_p_q[p,q]*(u**p)*(v**q)

	return fuv

def guv_sip(u,v,B_p_q):

	guv = np.double(0.0)
	print(B_p_q.shape[0],B_p_q.shape[1],u,v)

	for p in range(B_p_q.shape[0]):
		for q in range(B_p_q.shape[1]-p):
			guv += B_p_q[p,q]*(u**p)*(v**q)

	return guv

