#!/usr/bin/env python

import pysiaf
import numpy
print(pysiaf.JWST_PRD_VERSION)

instrument = "NIRCAM"

siaf = pysiaf.Siaf(instrument)
print(siaf)
aperture = ['NRCA1_FULL','NRCA2_FULL','NRCA3_FULL','NRCA4_FULL','NRCA5_FULL',
            'NRCB1_FULL','NRCB2_FULL','NRCB3_FULL','NRCB4_FULL','NRCB5_FULL']

for ap in aperture:
	nrc = siaf[ap]
	print("****")
	print("Aperture = ",ap)
	print("x_det_ref, y_det_ref = ",nrc.reference_point('det'))
	print("x_sci_ref, y_sci_ref = ",nrc.reference_point('sci'))
