import numpy as np
from astropy.io import fits

def xy_from_sip(u, v, CDij, Apq, Bpq, A_ORDER, B_ORDER):

	#define the vector (u,v)
	uv = np.array([u,v],dtype=np.float64)
	#print("uu,vv = ",uv)

	#add the f(u,v) and g(u,v) functions
	fuv = sip_fuv(u,v,Apq)
	guv = sip_guv(u,v,Bpq)
	#print("xsum, ysum = ",fuv,guv)

	uv[0] += fuv # u + f(u,v)
	uv[1] += guv # v + g(u,v)

	#evaluate the transform
	#print("xx = uu + xsum ",uv[0])
	#print("yy = vv + ysum ",uv[1])

	# (x,y) = ((CD11, CD12), (CD21, CD22)) * (u+f, v+g)
	xy = np.dot(CDij, uv)

	return xy

def sip_Apq_from_header(header):

	A_ORDER = header['A_ORDER']

	A_p_q = np.zeros((A_ORDER+1,A_ORDER+1),dtype=np.float64)

	for p in range(A_ORDER+1):
		for q in range(A_ORDER+1):
			s = 'A_%d_%d' % (p,q)
			try:
				A_p_q[p,q] = header[s]
			except KeyError:
				#s = ""
				print("No "+s)

	return A_p_q

def sip_Bpq_from_header(header):

	B_ORDER = header['B_ORDER']

	B_p_q = np.zeros((B_ORDER+1,B_ORDER+1),dtype=np.float64)

	for p in range(B_ORDER+1):
		for q in range(B_ORDER+1):
			s = 'B_%d_%d' % (p,q)
			try:
				B_p_q[p,q] = header[s]
			except KeyError:
				#s=""
				print("No "+s)

	return B_p_q


def sip_fuv(u,v,A_p_q):

	fuv = np.double(0.0)

	#print(A_p_q.shape[0],A_p_q.shape[1])

	for p in range(A_p_q.shape[0]):
		for q in range(A_p_q.shape[1]-p):
			fuv += A_p_q[p,q]*(u**p)*(v**q)

	return fuv

def sip_guv(u,v,B_p_q):

	guv = np.double(0.0)
	#print(B_p_q.shape[0],B_p_q.shape[1],u,v)

	for p in range(B_p_q.shape[0]):
		for q in range(B_p_q.shape[1]-p):
			guv += B_p_q[p,q]*(u**p)*(v**q)

	return guv

def sip_print_coeffs(Apq, Bpq):
	for p in range(Apq.shape[0]):
		for q in range(Apq.shape[1]-p):
			s = "p %d q %d Apq % 14.13e Bpq % 14.13e " % (p,q,Apq[p,q],Bpq[p,q])
			print(s)
