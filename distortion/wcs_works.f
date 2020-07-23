c     fake WCS using the STScI coefficients 
c
c     Calculate x_det, y_det for an object at ra, dec for PA_V3
c     requires
c     v2_ref (= v2_offset
c     v3_ref (= v3_offset)
c     ra_ref, dec_ref (corresponding to v2_ref, v3_ref)
c
      implicit none
c
      double precision x_det, y_det, x_det_ref, y_det_ref, 
     &     x_det_size, y_det_size, x_det_org, y_det_org
      double precision x_sci, y_sci, x_sci_ref, y_sci_ref, 
     &     x_sci_size, y_sci_size,det_sci_yangle
      double precision x_ideal, y_ideal
      double precision 
     &     sci_to_ideal_x,sci_to_ideal_y,
     &     ideal_to_sci_x, ideal_to_sci_y,
     &     x_sci_scale, y_sci_scale,
     &     v3_sci_x_angle,v3_sci_y_angle,
     &     v3_idl_yang

      double precision xterm, yterm, dx, dy, xsum, ysum
      double precision zp_idl_sci_x, zp_idl_sci_y
      double precision v2_offset, v3_offset
      integer ideal_to_sci_degree, v_idl_parity,
     &     sci_to_ideal_degree, det_sci_parity
      character filename*180, coeff_name*40, files*80,
     &     new_name*10
      dimension 
     &     sci_to_ideal_x(6,6), sci_to_ideal_y(6,6), 
     &     ideal_to_sci_x(6,6), ideal_to_sci_y(6,6)
c
      double precision c11, c12, c21, c22,cosdec
      double precision new_coeff, det_sign
      double precision aa, bb, ap, bp, xx, yy,
     &     fake_x, fake_y
c
      double precision crpix1, crpix2, crval1, crval2,
     &     cd1_1, cd1_2, cd2_1, cd2_2, 
     &     crpix3, crval3, cd3_3, cdelt3, equinox, roll_ref
      double precision inv1_1, inv1_2, inv2_1, inv2_2, detcd
      integer wcsaxes
      character ctype1*15, ctype2*15, ctype3*15, read_patt*10
      character cunit1*40,cunit2*40,cunit3*40
c
      double precision attitude_dir, attitude_inv
      integer sca, verbose, init, ii, jj, kk, ll, nterms, nx, ny,
     &     ix, iy, precise
      double precision v2_ref, v3_ref, ra_ref, dec_ref
      double precision v2, v3, ra, dec, pa_v3, q,
     &     ra_rad, dec_rad, pa_rad, radius, angle,cosx,sinx,
     &     v2_rad, v3_rad, v2_arcsec,v3_arcsec
      double precision ra_sca, dec_sca, ra_nrc, dec_nrc
      double precision rra, ddec, xlin, ylin, mag,
     &     mstar, z, semi_a, semi_b, pa, sersic_n
      double precision nrcall_v2, nrcall_v3, nrcall_v3idlyangle
      double precision uu, vv, ww, zz, d_a,d_d, u_sum, v_sum
      integer id, iexp, jexp, a_order, b_order, ap_order, bp_order
      integer nt_full, nt_truncated
c
      character subarray*8
      integer colcornr, rowcornr, naxis1, naxis2
c
      double precision zbqlnor
      integer seed
c
      dimension attitude_dir(3,3), attitude_inv(3,3)
      dimension aa(9,9), bb(9,9), ap(9,9), bp(9,9)
      dimension files(10)
c
      real image
      dimension image(2048, 2048)
c
      data nrcall_v2,nrcall_v3,nrcall_v3idlyangle/-0.3276288867102224d0,
     &     -492.96588887005885d0, -0.11255124074467837d0/
      data aa, bb, ap, bp/81*0.0d0,81*0.0d0,81*0.0d0,81*0.0d0/
c
      double precision x_sci_scale_norm, y_sci_scale_norm
c

      files(1)   = './data/NRCA1_FULL.ascii'
      files(2)   = './data/NRCA2_FULL.ascii'
      files(3)   = './data/NRCA3_FULL.ascii'
      files(4)   = './data/NRCA4_FULL.ascii'
      files(5)   = './data/NRCA5_FULL.ascii'
c
      files(6)   = './data/NRCB1_FULL.ascii'
      files(7)   = './data/NRCB2_FULL.ascii'
      files(8)   = './data/NRCB3_FULL.ascii'
      files(9)   = './data/NRCB4_FULL.ascii'
      files(10)  = './data/NRCB5_FULL.ascii'
c
      q          = dacos(-1.0d0)/180.0d0
      verbose    = 1
      precise    = 0
      subarray   = 'FULL'
      colcornr   = 1
      rowcornr   = 1
      naxis1     = 2048
      naxis2     = 2048
      wcsaxes    = 2
      equinox    = 2000.d0
      a_order  = 5
      b_order  = 5
      ap_order = 5
      bp_order = 5
c
c     Position of NIRCam (NRCALL)
c
c      ra_ref  =  53.157734d0
c      dec_ref = -27.807316d0
c      pa_v3   =  30.d0
      ra_nrc  =  53.160149d0
      dec_nrc = -27.794571d0
      pa_v3   =  41.d0
      roll_ref = pa_v3

      x_det = 1024.50d0
      y_det = 1024.50d0
      x_det = 1.d0
      y_det = 1.d0
c      x_det = 820.528 
c      y_det = 847.386
c      x_det = 255.88304037368289
c      y_det = 1963.3211995498791
cc      
      x_det = 1025.50d0
      y_det = 1025.50d0
      x_det_org = x_det
      y_det_org = y_det
c
c     Process NRCLONGB (B5)
c
      sca      = 1
      call read_siaf_parameters(sca,
     &     sci_to_ideal_x, sci_to_ideal_y, sci_to_ideal_degree,
     &     ideal_to_sci_x, ideal_to_sci_y, ideal_to_sci_degree,
     &     x_det_ref, y_det_ref, x_sci_ref, y_sci_ref,
     &     x_sci_scale, y_sci_scale,
     &     det_sci_yangle, det_sci_parity,
     &     v3_idl_yang, v_idl_parity, v2_ref, v3_ref, verbose)
      
      call prep_wcs(
     &     sci_to_ideal_x, sci_to_ideal_y, sci_to_ideal_degree,
     &     ideal_to_sci_x, ideal_to_sci_y, ideal_to_sci_degree,
     &     x_det_ref, y_det_ref, x_sci_ref, y_sci_ref,
     &     det_sci_yangle, det_sci_parity,
     &     v3_idl_yang, v_idl_parity, v2_ref, v3_ref,
     &     crpix1, crpix2,
     &     crval1, crval2,
     &     ctype1, ctype2,
     &     cunit1, cunit2,
     &     cd1_1, cd1_2, cd2_1, cd2_2,
     &     ra_nrc, dec_nrc, pa_v3,
     &     a_order, aa, b_order, bb,
     &     ap_order, ap, bp_order, bp,
     &     attitude_dir, attitude_inv, 
     &     ra_sca, dec_sca,
     &     verbose)
      dec_rad = dec_sca *q
      cosdec = dcos(dec_rad)
c
c     inverse  (dra, ddec) -> (x_idl, y_idl)
c     (not saved, but necessary for testing)
c
      detcd = (cd1_1*cd2_2) - (cd1_2*cd2_1)
      inv1_1 =  cd2_2/detcd
      inv1_2 = -cd1_2/detcd
      inv2_1 = -cd2_1/detcd
      inv2_2 =  cd1_1/detcd
c
      uu      = x_det - crpix1
      vv      = y_det - crpix2
 20   format('ORG_',i1,'_',i1)
 30   format('AP_',i1,'_',i1)
c
c-----------------------------------------------------------------------
c
      call det_to_sci(x_det, y_det, x_sci, y_sci, 
     &     x_det_ref, y_det_ref, x_sci_ref, y_sci_ref, 
     &     det_sci_yangle, det_sci_parity)
      call sci_to_ideal(x_sci, y_sci, x_sci_ref, y_sci_ref,
     &     x_ideal, y_ideal,  sci_to_ideal_x, sci_to_ideal_y, 
     &     sci_to_ideal_degree,verbose)
      call det_to_v2v3(
     &     x_det_ref, y_det_ref,
     &     x_sci_ref, y_sci_ref,
     &     sci_to_ideal_x,sci_to_ideal_y,sci_to_ideal_degree,
     &     ideal_to_sci_x, ideal_to_sci_y,ideal_to_sci_degree,
     &     v3_sci_x_angle,v3_sci_y_angle,
     &     v3_idl_yang,v_idl_parity,
     &     det_sci_yangle, det_sci_parity,
     &     v2_ref, v3_ref,
     &     v2_arcsec, v3_arcsec, x_det, y_det,
     &     precise,verbose)
      v2_rad = q*v2_arcsec/3600.d0
      v3_rad = q*v3_arcsec/3600.d0
      call rot_coords(attitude_dir, v2_rad, v3_rad, ra_rad, dec_rad)
      call coords_to_ra_dec(ra_rad, dec_rad, ra_ref,dec_ref)
      print *, '"truth":'
      print *, 'ra_nrc, dec_nrc     ', ra_nrc, dec_nrc
      print *, 'ra_sca, dec_sca, sca', ra_sca, dec_sca, sca
      print *, 'v2_sca, v3_sca      ', v2_ref, v3_ref
      print *, 'x_det, y_det        ', x_det, y_det
      print *, 'x_sci, y_sci        ', x_sci, y_sci
      print *, 'x_sci_scale, y_sci_scale     ', x_sci_scale, y_sci_scale
      x_sci_scale_norm = x_sci_scale/3600.
      y_sci_scale_norm = y_sci_scale/3600.
      print *, 'normed sci_scale', x_sci_scale_norm, y_sci_scale_norm
      print *, 'x_ideal, y_ideal    ', x_ideal, y_ideal
      print *, 'v2_arcsec, v3_arcsec', v2_arcsec, v3_arcsec
      print *, 'ra_rad, dec_rad     ', ra_rad, dec_rad
      print *, 'ra_matr, dec_matr   ', ra_ref, dec_ref

 22   format('ORG_',i1,'_',i1)
 23   format('BP_',i1,'_',i1)
 24   format(2(2x,a10,2x,e18.11),2(2x,i2))
 60   format('BP_',i1,'_',i1)
c
c     Using TAN-SIP coefficients
c
      uu = x_det - crpix1
      vv = y_det - crpix2 
      xsum  = 0.0d0
      ysum  = 0.0d0
      print 99
      do ii = 1, a_order + 1
         iexp = ii - 1
           do jj = 1, a_order + 1
            jexp = jj - 1
            if(iexp+jexp.le. a_order) then
               write(coeff_name, 100) iexp, jexp, iexp, jexp
               if(verbose.gt.0) then
                  print 110, coeff_name, aa(ii,jj), bb(ii,jj),
     &                 ii, jj, iexp, jexp
               end if
               xterm =  aa(ii, jj) * (uu**iexp) * (vv**jexp)
               xsum = xsum + xterm
               print 115, ii-1, jj-1, aa(ii,jj),uu, vv, xterm, xsum
            end if
         end do
      end do
c
      print *,' '
      print 99
      do ii = 1, a_order + 1
         iexp = ii - 1
           do jj = 1, a_order + 1
            jexp = jj - 1
            if(iexp+jexp.le. a_order) then
               write(coeff_name, 100) iexp, jexp, iexp, jexp
               if(verbose.gt.0) then
                  print 110, coeff_name, aa(ii,jj), bb(ii,jj),
     &                 ii, jj, iexp, jexp
               end if
               yterm = bb(ii, jj) * (uu**iexp) * (vv**jexp)
               ysum = ysum +  yterm
               print 116, ii-1, jj-1, bb(ii,jj),uu, vv, yterm, ysum
            end if
         end do
      end do
 99   format('aa or bb',3x,'value', 11x,'uu',9x,'vv',
     &     8x,'term',10x,'sum')
 100           format('A_',i1,'_',i1,2x,'B_',i1,'_',i1)
 110              format(A15,2(2x, e18.11),4(i6))
 115           format('aa_',i1,'_',i1,2(1x,e13.6,2(1x,f10.4)))
 116           format('bb_',i1,'_',i1,2(1x,e13.6,2(1x,f10.4)))
c
 117  format(' cd1_1',2x,'cd1_2',13x,2(1x,e13.6))
 118  format(' cd2_1',2x,'cd2_2',13x,2(1x,e13.6))
c
      xx   =  uu + xsum
      yy   =  vv + ysum
      dy   =  (cd2_1 * xx + cd2_2 * yy)
      ddec = crval2 + dy
      dx   =  (cd1_1 * xx + cd1_2 * yy)/cosdec
      rra  = crval1 +  dx
      print *,'using TAN-SIP coefficients:'
      print *, 'x_det, y_det            ', x_det, y_det
      print *, 'uu, vv                  ', uu, vv
      print *, 'xsum,  ysum             ', xsum, ysum
      print *, 'xx = uu + xsum          ', xx
      print *, 'yy = vv + ysum          ', yy
      print 117, cd1_1, cd1_2
      print 118, cd2_1, cd2_2
      print *, '(cd1_1 * xx + cd1_2 * yy)/cosdec', dx
      print *, '(cd2_1 * xx + cd2_2 * yy)       ', dy
      print *, 'rra, ddec               ', rra, ddec
      print *, 'using full matrices     ', ra_ref, dec_ref
      print *, 'difference (arc sec)    ', 
     &     (ra_ref-rra)*3600.d0*cosdec, (dec_ref-ddec)*3600.d0
      print *, 'x_ideal-xx, y_ideal-yy  ',x_ideal-xx,
     &     y_ideal-yy
 130     format('dir ',8(1x,f16.10))
 150     format(i8,14(2x,f12.7))
 140  format('inv ',8(1x,f16.10))
 300  continue
      stop
c####################
      print *,' '
      print *,'inverse'
c
c--------------------------------------------------
c
      call rot_coords(attitude_inv, ra_rad, dec_rad, 
     &     v2_rad, v3_rad)
      call coords_to_v2v3(v2_rad, v3_rad,v2_arcsec,v3_arcsec)
      call v2v3_to_det(
     &     x_det_ref, y_det_ref, 
     &     x_sci_ref, y_sci_ref,
     &     sci_to_ideal_x,sci_to_ideal_y,sci_to_ideal_degree,
     &     ideal_to_sci_x,ideal_to_sci_y,ideal_to_sci_degree,
     &     v3_sci_x_angle,v3_sci_y_angle,
     &     v3_idl_yang, v_idl_parity,
     &     det_sci_yangle,det_sci_parity,
     &     v2_ref, v3_ref,
     &     v2_arcsec, v3_arcsec, x_det, y_det,
     &     precise,2)

      print *, '"truth"'
      print *, 'ra_rad, dec_rad        ', ra_rad, dec_rad
      print *, 'ra_ref, dec_ref        ', ra_ref, dec_ref
      print *, 'v2_arcsec, v3_arcsec   ', v2_arcsec, v3_arcsec
      print *, 'x_ideal, y_ideal (rep) ', x_ideal, y_ideal
      print *, 'x_det, y_det (input)   ', x_det_org, y_det_org
      print *, 'x_det, y_det (recov)   ', x_det, y_det
      print *, 'diff_x, diff_y (pixel) ', x_det_org-x_det, 
     &     y_det_org-y_det
      print *,' '
c
c     starting from "true" position
c
c      d_a = (ra_ref  - crval1) * cosdec
c      d_d = (dec_ref - crval2) 
c      print *,'using ww, zz ("ideal")  '
c
c      print *, 'starting from rra, ddec', rra, ddec
      d_a = (rra  - crval1) * cosdec
      d_d = (ddec - crval2) 
c     calculate the "ideal" coordinates:
      ww  =  inv1_1 * d_a + inv1_2 * d_d
      zz  =  inv2_1 * d_a + inv2_2 * d_d
c     
      print *,'ww, zz                  ', ww, zz
      print *,'x_ideal-ww, y_ideal-zz  ', x_ideal-ww,
     &     y_ideal-zz
c
      u_sum = 0.d0
      do ii = 1, ap_order+1
         iexp = ii -1
         do jj = 1, ap_order+1
            jexp = jj - 1
            if(iexp+jexp .le. ap_order) then
               u_sum = u_sum + ap(ii,jj)*(ww**(ii-1))*(zz**(jj-1))
            end if
         end do
      end do
c
      v_sum = 0.d0
      do ii = 1, bp_order+1
         iexp = ii -1
         do jj = 1, bp_order+1 
            jexp = jj - 1
            if(iexp+jexp .le. bp_order) then
               if(verbose.gt.0) then
                  write(coeff_name, 100) iexp, jexp, iexp, jexp
                  print 110, coeff_name, ap(ii,jj), bp(ii,jj),
     &                 ii, jj, iexp, jexp
               end if
               v_sum = v_sum + 
     &              bp(ii,jj)*(ww**(ii-1))*(zz**(jj-1))
            end if
         end do
      end do
c
      det_sign  = dcos(q*det_sci_yangle)

      xx = crpix1 + ww + u_sum
      yy = crpix2 + zz + v_sum
c
      print *,' SCA, v_idl_parity     ', sca, det_sign,det_sci_yangle
      print *,'crpix1, crpix2         ', crpix1, crpix2
      print *,'u_sum, v_sum           ', u_sum, v_sum
      print *,'ww+u_sum, zz+v_sum     ',ww+ u_sum, zz+v_sum
      print *,'x_det, ydet (TAN-SIP)  ', xx, yy
      print *,'x_det_org, y_det_org   ',x_det_org, y_det_org
      print *,'diff                   ', x_det_org-xx, 
     &     y_det_org-yy
c
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
c---------------------------------------------------
c
c      print 10, ctype1, ctype2, crval1, crval2, crpix1,crpix2,
c     &     cd1_1,cd1_2,cd2_1,cd2_2
 10   format(
     &     'CTYPE1  = ',A15,/,
     &     'CTYPE2  = ',A15,/,
     &     'CRVAL1  = ',f15.8,/,
     &     'CRPIX2  = ',f15.8,/,
     &     'CRPIX1  = ',f15.3,/,
     &     'CRPIX2  = ',f15.3,/,
     &     'CD1_1   = ',e15.8,/,
     &     'CD1_1   = ',e15.8,/,
     &     'CD2_1   = ',e15.8,/,
     &     'CD2_2   = ',e15.8)
c
      open(2,file = 'tan_sip3.dat')
      write(2, 310)  crval1, crval2, crpix1,crpix2,
     &     cd1_1,cd1_2,cd2_1,cd2_2
 310  format(8(e15.8,/))
      write(2,320) a_order, b_order
 320  format(2(2x,i6))

      do ii = 1, a_order + 1
         iexp = ii -1
         do jj = 1, a_order + 1
            jexp = jj-1
            if((iexp+ jexp) .le. a_order) then
               write(coeff_name, 30) iexp, jexp
               write(2, 330) ii, jj, aa(ii,jj), coeff_name
 330           format(2(2x,i2),2x,e15.8,2x,a10)
            end if
         end do
      end do
c
      do ii = 1, b_order + 1
         iexp = ii -1
         do jj = 1, b_order + 1
            jexp = jj-1
            if((iexp+ jexp) .le. b_order) then
               write(coeff_name, 60) iexp, jexp
               write(2, 330) ii, jj, bb(ii,jj),coeff_name
            end if
         end do
      end do
      write(2,320) ap_order, bp_order
      do ii = 1, ap_order + 1
         iexp = ii -1
         do jj = 1, ap_order + 1
            jexp = jj-1
            if((iexp+ jexp) .le. ap_order) then
               write(coeff_name, 410) iexp, jexp
 410           format('AP_',i1,'_',i1)
               write(2, 330) ii, jj, ap(ii,jj),coeff_name
            end if
         end do
      end do
c
      do ii = 1, bp_order + 1
         iexp = ii -1
         do jj = 1, bp_order + 1
            jexp = jj-1
            if((iexp+ jexp) .le. bp_order) then
               write(coeff_name, 410) iexp, jexp
 420           format('BP_',i1,'_',i1)
               write(2, 330) ii, jj, bp(ii,jj),coeff_name
            end if
         end do
      end do
      close(2)
c
c--------------------------------------------------
c
      seed = 0
      call zbqlini(seed)
      do jj = 1, 2048
         do ii = 1, 2048
               image(ii,jj) = zbqlnor(10.d0, 5.d0)
         end do
      end do
c
      filename = 'test3.fits'
      call write_float_2d_image(filename, image, naxis1, naxis2,
     &     subarray, colcornr, rowcornr,
     &     wcsaxes, 
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ctype1, ctype2, ctype3,
     &     cunit1, cunit2, cunit3,
     &     cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     cdelt3,
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, read_patt,
     &     v3_idl_yang, v_idl_parity,
     &     a_order, aa, b_order, bb,
     &     ap_order, ap, bp_order, bp)
c      close(7)
c      close(8)
      stop
      end
c
c----------------------------------------------------------------------     
c
      subroutine write_float_2d_image(filename, image, naxis1, naxis2,
     &     subarray, colcornr, rowcornr,
     &     wcsaxes, 
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ctype1, ctype2, ctype3,
     &     cunit1, cunit2, cunit3,
     &     cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     cdelt3,
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, read_patt,
     &     v3i_yang, vparity, 
     &     a_order, aa, b_order, bb,
     &     ap_order, ap, bp_order, bp)
c
      implicit none
      character telescop*20, instrume*20, filter*20,
     *     module*20, partname*4,comment*45, object*20      
      real image
      double precision tframe, tgroup
      integer status, bitpix, naxes, naxis, pcount, gcount, block,
     *     groupgap, group, sca_id, nframe, ngroup,
     *     nnn, job, iunit
      logical simple,extend
      character subarray*8
      integer colcornr, rowcornr, naxis1, naxis2
      character filename*(*)

      double precision cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, cdelt3
      double precision v3i_yang, aa, bb, ap, bp
      integer wcsaxes, vparity
      integer a_order, b_order, ap_order, bp_order
      character read_patt*10
      character card*80
      character ctype1*15, ctype2*15, ctype3*15
      character cunit1*40, cunit2*40, cunit3*40
c
      dimension aa(9,9), bb(9,9), ap(9,9), bp(9,9)
c      dimension aa(a_order+1, a_order+1), bb(b_order+1, b_order+1),
c     &     ap(ap_order+1, ap_order+1), bp(bp_order+1, bp_order+1)
c
      parameter (nnn=2048)
c
      dimension image(nnn,nnn),naxes(2)
c
c     define the required primary array keywords
c
      bitpix   = -32
      simple   = .true.
      extend   = .true.
      naxis    =  2
      naxes(1) = naxis1
      naxes(2) = naxis2
c
      status = 0
      iunit = 91
      call ftgiou(iunit, status)
      if (status .gt. 0) then 
         call printerror(status)
         print *,'iunit ',iunit, status
         stop
      end if
c     
c     delete previous version of the file, if it exists
c
      call ftopen(iunit, filename, 1, block, status)
      if (status .eq. 0) then
         call ftdelt(iunit, status)
      else
c     clear the error message stack
         call ftcmsg
      end if
      status = 0
c
c     create the fits file
c
      call ftinit(iunit, filename, 1, status)
      if (status .gt. 0) then 
         call printerror(status)
         print *,'iunit ',iunit
         print *,'pause: enter return to continue'
         read(*,'(A)')
      endif
      status =  0
c     
      call ftphpr(iunit,simple, bitpix, naxis, naxes, 
     & 0,1,extend,status)
      if (status .gt. 0) then
         print *, 'simple,bitpix,naxis'
         call printerror(status)
      end if
      status =  0
c
c     write more header keywords
c
      
      call  fake_wcs(
     &     iunit, wcsaxes, 
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ctype1, ctype2, ctype3,
     &     cunit1, cunit2, cunit3,
     &     cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     cdelt3,
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, read_patt,
     &     v3i_yang, vparity, 
     &     a_order, aa, b_order, bb,
     &     ap_order, ap, bp_order, bp)
c
c      call write_nircam_keywords(iunit,nframe, 
c     *     real(tframe), groupgap, real(tgroup), 
c     *     ngroup, object, partname, sca_id,module, filter, 
c     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
c     Write out image
c
      print *,'write float matrix'
      group=0
      call ftp2de(iunit,group,nnn, naxes(1),naxes(2),
     *     image,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'writing images'
      end if
c     
      call closefits(iunit)
      return
      end
c
c----------------------------------------------------------------------
c-------------------------------------------------------------------
c
      subroutine  fake_wcs(
     &     iunit, wcsaxes, 
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ctype1, ctype2, ctype3,
     &     cunit1, cunit2, cunit3,
     &     cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     cdelt3,
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, read_patt,
     &     v3i_yang, v_idl_parity,
     &     a_order, aa, b_order, bb,
     &     ap_order, ap, bp_order, bp)
      implicit none
      integer iunit, status
      double precision cd1_1, cd1_2, cd2_1, cd2_2, cd3_3,
     &     crpix1, crpix2, crpix3, 
     &     crval1, crval2, crval3, 
     &     ra_ref, dec_ref, roll_ref, pa_v3,
     &     equinox, cdelt3
      double precision v3i_yang, aa, bb, ap, bp
      integer wcsaxes, v_idl_parity
      integer a_order, b_order, ap_order, bp_order
      character read_patt*10
      character card*80, comment*40
      character ctype1*15, ctype2*15, ctype3*15
      character cunit1*40, cunit2*40, cunit3*40
c
      integer iexp, jexp, ii, jj, kk, ll
      character coeff_name*10
c
      dimension aa(9,9), bb(9,9), ap(9,9), bp(9,9)
c      dimension aa(a_order+1, a_order+1), bb(b_order+1, b_order+1),
c     &     ap(ap_order+1, ap_order+1), bp(bp_order+1, bp_order+1)
c
c     fake WCS keywords
c 
      card = '                                              '
      call ftprec(iunit,card, status)
      status =  0
      card = '         WCS information                             '
      call ftprec(iunit,card, status)
      status =  0
      card = '                                              '
      call ftprec(iunit,card, status)
      status =  0
c
      comment='                                       '
      call ftpkyj(iunit,'WCSAXES',wcsaxes,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'WCSAXES'
      end if
      status =  0
c     
      comment='Axis 1 coordinate of reference pixel   '
      call ftpkyd(iunit,'CRPIX1',crpix1,-7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRPIX1'
      end if
      status =  0
c     
      comment='Axis 2 coordinate of reference pixel   '
      call ftpkyd(iunit,'CRPIX2',crpix2,-7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRPIX2'
      end if
      status =  0
c     
      if(wcsaxes.eq.3) then
         comment='Axis 3 coordinate of reference pixel   '
         if(read_patt .eq. 'RAPID') then
            call ftpkyd(iunit,'CRPIX3',crpix3,-7,comment,status)
         else
            call ftpkyd(iunit,'CRPIX3',crpix3/2.d0,-7,comment,status)
         end if
         if (status .gt. 0) then
            call printerror(status)
            print *, 'CRPIX3'
         end if
         status =  0
      end if
c      
      comment='RA at reference pixel (degrees)        '
      call ftpkyd(iunit,'CRVAL1',crval1,-12,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRVAL1'
      end if
      status =  0
c     
      comment='DEC at reference pixel (degrees)       '
      call ftpkyd(iunit,'CRVAL2',crval2,-12,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRVAL2'
       end if
       status =  0
c     
       if(wcsaxes .eq. 3) then
          comment='T at reference pixel (seconds)       '
          call ftpkyd(iunit,'CRVAL3',crval3,-12,comment,status)
          if (status .gt. 0) then
             call printerror(status)
             print *, 'CRVAL3'
          end if
          status =  0
       end if
c     
      comment='Projection type                        '
      call ftpkys(iunit,'CTYPE1',ctype1,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CTYPE1'
      end if
      status =  0
c     
      call ftpkys(iunit,'CTYPE2',ctype2,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CTYPE2'
      end if
      status =  0
c     
      if(wcsaxes .eq. 3) then
         call ftpkys(iunit,'CTYPE3','',comment,status)
         if (status .gt. 0) then
            call printerror(status)
            print *, 'CTYPE3'
         end if
         status =  0
      end if
c
      cunit1 = 'deg'
      comment='First axis units                       '
      call ftpkys(iunit,'CUNIT1',cunit1,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CUNIT1'
      end if
      status =  0
c     
      cunit2 = 'deg'
      comment='Second axis units                      '
      call ftpkys(iunit,'CUNIT2',cunit2,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CUNIT2'
      end if
      status =  0
c     
      if(wcsaxes .eq. 3) then
         cunit3 = 'sec'
         comment='Third axis units                      '
         call ftpkys(iunit,'CUNIT3',cunit3,comment,status)
         if (status .gt. 0) then
            call printerror(status)
            print *, 'CUNIT3'
         end if
         status =  0
      end if
c     
      comment='RA  of the reference point (deg)'
      call ftpkyd(iunit,'RA_REF',ra_ref,-16,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'ra_ref'
      end if
      status =  0
c     
      comment='Dec of the reference point (deg)'
      call ftpkyd(iunit,'DEC_REF',dec_ref,-16,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'dec_ref'
      end if
      status =  0
c     
      comment='Telescope roll angle of V3 at ref point'
      call ftpkyd(iunit,'ROLL_REF',roll_ref,-8,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'roll_ref'
      end if
      status =  0
c     
      comment='Relative sense of rotation between Ideal xy and V2V3 '
      call ftpkyj(iunit,'V_IDL_PARITY',v_idl_parity,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'v_idl_parity'
      end if
      status =  0
c     
      comment='Angle from V3 axis to Ideal y axis (deg)'
      call ftpkyd(iunit,'V3I_YANG',v3i_yang,4,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'V3I_YANG'
      end if
      status =  0
c     
c     These are non-STScI standard
c     
c     
c      if(naxis .eq.4) then
c         comment='Axis 4 coordinate of reference pixel   '
c         call ftpkyd(iunit,'CRPIX4',1.d0,-7,comment,status)
c         if (status .gt. 0) then
c            call printerror(status)
c            print *, 'CRPIX4'
c         end if
c         status =  0
c         comment='4th dimension at reference pixel (pixel)   '
c         call ftpkyd(iunit,'CRVAL4',1.d0,-7,comment,status)
c         if (status .gt. 0) then
c            call printerror(status)
c            print *, 'CRVAL3'
c         end if
c         
c         status =  0
c         comment='Fourth axis increment per pixel         '      
c         call ftpkyd(iunit,'CDELT4',1.d0,12,comment,status)
c         if (status .gt. 0) then
c            call printerror(status)
c            print *, 'CDELT4'
c         end if
c         status =  0
c      end if
c
      comment='                                       '
      call ftpkyd(iunit,'EQUINOX',equinox,-7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'EQUINOX'
      end if
      status =  0
c
      comment='                                       '
      call ftpkyd(iunit,'CD1_1',cd1_1,10,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD1_1'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD1_2',cd1_2,10,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD1_2'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD2_1',cd2_1,10,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD2_1'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD2_2',cd2_2,10,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD2_2'
      end if
      status =  0
c     
      if(wcsaxes .eq. 3) then
         comment='                                       '
         cd3_3  = cdelt3
         call ftpkyd(iunit,'CD3_3',cd3_3,10,comment,status)
         if (status .gt. 0) then
            call printerror(status)
            print *, 'CD3_3'
         end if
         status =  0
      end if
c
      comment='polynomial order axis 1, detector to sky'
      call ftpkyj(iunit,'A_ORDER',a_order,comment,status)
      if (status .gt. 0) then
         print *, 'A_ORDER'
         call printerror(status)
      end if
      status =  0
c     
      comment='distortion coefficient'
      do ii = 1, a_order+1
         iexp = ii - 1
         do jj = 1, a_order+1
            jexp = jj - 1
            if(iexp+jexp.le. a_order.and.aa(ii,jj).ne.0.0d0) then
               write(coeff_name, 100) iexp, jexp
 100           format('A_',i1,'_',i1)
               call ftpkyd(iunit,coeff_name,aa(ii,jj),10,comment,status)
c               print 115, coeff_name, aa(ii,jj)
 115           format(a10,2x,e15.8)
               if (status .gt. 0) then
                  print 130,coeff_name, aa(ii,jj)
 130              format(a8,2x,e15.8)
                  call printerror(status)
               end if
               status =  0
            end if
         end do
      end do
c
      comment='polynomial order axis 2, detector to sky'
      call ftpkyj(iunit,'B_ORDER',b_order,comment,status)
      if (status .gt. 0) then
         print *, 'B_ORDER'
         call printerror(status)
      end if
      status =  0
c     
      comment='distortion coefficient'
      do ii = 1, b_order+1
         iexp = ii - 1
         do jj = 1, b_order+1
            jexp = jj - 1
            if(iexp+jexp.le. b_order .and.bb(ii,jj).ne.0.d0) then
               write(coeff_name, 140) iexp, jexp
 140           format('B_',i1,'_',i1)
               call ftpkyd(iunit,coeff_name,bb(ii,jj),10,comment,status)
c     print 115, coeff_name, bb(ii,jj)
               if (status .gt. 0) then
                  print 130,coeff_name, bb(ii,jj)
                  call printerror(status)
               end if
               status =  0
            end if
         end do
      end do
c
c     inverse 
c
      comment='polynomial order axis 1, sky to detector'
      call ftpkyj(iunit,'AP_ORDER',ap_order,comment,status)
      if (status .gt. 0) then
         print *, 'AP_ORDER'
         call printerror(status)
      end if
      status =  0
c     
      comment='distortion coefficient'
      do ii = 1, ap_order+1
         iexp = ii - 1
         do jj = 1, ap_order+1
            jexp = jj - 1
            if(iexp+jexp.le. ap_order.and.ap(ii,jj).ne.0.0d0) then
               write(coeff_name, 150) iexp, jexp
 150           format('AP_',i1,'_',i1)
               call ftpkyd(iunit,coeff_name,ap(ii,jj),10,comment,status)
c               print 115, coeff_name, ap(ii,jj)
               if (status .gt. 0) then
                  print 130,coeff_name, ap(ii,jj)
                  call printerror(status)
               end if
               status =  0
            end if
         end do
      end do
      
      comment='polynomial order axis 2, sky to detector'
      call ftpkyj(iunit,'BP_ORDER',bp_order,comment,status)
      if (status .gt. 0) then
         print *, 'BP_ORDER'
         call printerror(status)
      end if
      status =  0
c     
      comment='distortion coefficient'
      do ii = 1, bp_order+1
         iexp = ii - 1
         do jj = 1, bp_order+1
            jexp = jj - 1
            if(iexp+jexp.le. bp_order.and.bp(ii,jj).ne.0.0d0) then
               write(coeff_name, 160) iexp, jexp
 160           format('BP_',i1,'_',i1)
               call ftpkyd(iunit,coeff_name,bp(ii,jj),10,comment,status)
c               print 115, coeff_name, bp(ii,jj)
               if (status .gt. 0) then
                  print 130,coeff_name, bp(ii,jj)
                  call printerror(status)
               end if
               status =  0
            end if
         end do
      end do
      
      return
      end
c
c-----------------------------------------------------------------------
      subroutine printerror(status)
C     Print out the FITSIO error messages
      integer status
      character errtext*30,errmessage*80

C     check if status is OK (no error); if so, simply return
      if (status .le. 0)return

C     get the text string which describes the error
 1    call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

C     read and print out all the error messages on the FITSIO stack
 2    call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end
c-----------------------------------------------------------------------
      subroutine closefits(unit)
      integer unit, status
 1    status=0

C     close the file, 
 9    call ftclos(unit, status)
c     release this logical unit
      call ftfiou(unit,status)

 10   if (status .gt. 0) call printerror(status)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine new_coeffs(degree, parity, det_sign,
     &     poly, new, order, verbose)
      implicit none
      integer parity, degree, verbose
      double precision det_sign, poly,temp_coeff,
     &     new
      integer ii, jj, kk, ll, iexp, jexp, order
      character coeff_name*10, new_name*10
c      dimension poly(6,6), new(6,6)
      dimension poly(6,6), new(9,9)
      
      do  ii = 2, degree
         iexp = ii - 1
         do jj = 1, ii
            jexp = jj - 1
            temp_coeff = (parity**(iexp-jexp)) * (det_sign**iexp) 
     &           * poly(ii,jj)
            kk = iexp-jexp + 1
            ll = jexp      + 1
            if(kk+ll-2.le.order) then
               new(kk,ll) = temp_coeff
               if(verbose.gt.0) then 
                  write(coeff_name, 20) iexp, jexp
 20               format('ORG_',i1,'_',i1)
                  write(new_name, 30) kk-1, ll-1
 30               format('NEW_',i1,'_',i1)
                  print  40, coeff_name, poly(ii,jj),
     &                 new_name, temp_coeff, iexp-jexp,jexp
 40               format(2(2x,a10,2x,e18.11),2(2x,i2))
               end if
            end if
         end do
      end do
      return
      end

c
c---------------------------------------------------------------
c
c
c      u_sum = 0.d0
c      fake_x =0.d0
c      nt_full = 0
c      nt_truncated = 0
c      do ii = 2, ideal_to_sci_degree
c         iexp = ii -1
c         do jj = 1, ii
c            jexp = jj - 1
c            u_sum = u_sum + ideal_to_sci_x(ii,jj)
c     &           * (ww**(iexp-jexp)) * (zz**jexp)
c            nt_full = nt_full + 1
c            kk = iexp-jexp + 1
c            ll = jexp      + 1
c            if(kk+ll-2 .le.ap_order) then
cc            if(iexp+jexp .le.ap_order) then
c               ap(kk,ll)   = ideal_to_sci_x(ii,jj)
c               fake_x = fake_x + ap(kk,ll) * (ww**(kk-1)) *
c     &              (zz**(ll-1))
c               nt_truncated = nt_truncated + 1
c            else
c               print *,'skip coeff ii, jj',ii, jj, ideal_to_sci_x(ii,jj)
c            end if
c         end do
c      end do
cc      u_sum = u_sum * det_sci_parity * det_sign
cc
c      v_sum = 0.d0
c      fake_y = 0.0d0
c      do ii = 2, ideal_to_sci_degree
c         iexp = ii -1
c         do jj = 1, ii
c            jexp = jj - 1
c            new_coeff = ideal_to_sci_y(ii,jj) * det_sign
c            v_sum = v_sum + new_coeff
c     &           * (ww**(iexp-jexp)) * (zz**jexp)
c            kk = iexp-jexp + 1
c            ll = jexp      + 1
c            if(kk+ll-2 .le. bp_order) then
cc            if(iexp+jexp .le. bp_order) then
c               bp(kk,ll)   = new_coeff
c               fake_y = fake_y + bp(kk,ll) * (ww**(kk-1)) *
c     &              (zz**(ll-1))
c            else
c               print *,'skip coeff ii, jj',ii, jj, 
c     &              ideal_to_sci_y(ii,jj),
c     &              new_coeff*(ww**(iexp-jexp))*(zz**jexp)
c            end if
c         end do
c      end do
c
c      xx   = u_sum + crpix1 
c      yy   = v_sum + crpix2 
c      fake_x = fake_x + crpix1
c      fake_y = fake_y + crpix2
c      print *,'x_det all coeffs ', xx
c      print *,'y_det all coeffs ', yy, nt_full
c      print *,' diff            ', x_det - xx, y_det-yy
c      print *,' zp ',zp_idl_sci_x,zp_idl_sci_y
c      print *,'x_det truncated  ', fake_x
c      print *,'y_det truncated  ', fake_y, nt_truncated
c      print *,' diff            ', x_det - fake_x, y_det-fake_y
