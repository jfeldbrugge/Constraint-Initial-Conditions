
      PROGRAM HOFRIB



c     --------------------------------------------------------------------
c     This is the 3-D version of the Hofman-Ribak algorithm for 
c     setting up a constrained random field.
c     The program consists of seven basic parts:
c        I.    Set up lattice and power spectrum
c        II.   Set up unconstrained field sample
c        III.  Input of constraints
c        IV.   Compute covariance matrix of constraints
c        V.    Computing either mean field or constrained realization
c        VI.   Checking the values of the constraints for the 
c              constrained realization.
c        VII.  Fouriertransforming the density field and computing 
c              displacement psi-field in ftpsi
c
c     SET n1,n2,n3 = THE NUMBER OF LATTICE POINTS PER DIMENSION.  SET IN
c     main AND constr.  SET ncon IN constr.
c     N.B. THE USER MUST PROVIDE SUBROUTINE constr(con,g,ncon)!!!
c
c     Note: i,i1,i2 label constraints; j labels spatial grid points; 
c     k,k1,k2,k3 label grid points in Fourier space.  There are n1 
c     lattice points in each dimension for a total of ngrid=n1*n1*n1.  
c     In Fourier space there are (n1/2+1)*n1*n1 grid points (k1 < 0 is 
c     redundant because rho is real).
c     -------------------------------------------------------------------

         
      implicit double precision (a-h,o-z)

      double precision   twopi
      parameter          (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter          (twopi=6.283185307179586d0,sqr2=1.41421356)
      parameter          (n0n0=n0*n0,ngrid0=n0n0*n0,ngr0=ngrid0/2)
      parameter          (ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer            nharm
      double precision   sigma,sigmav,chisq
      double precision   q(nconmax,nconmax),qinv(nconmax,nconmax)
      real               g(nconmax),g0(nconmax),ghelp(nconmax)
      complex            conc(ngdim0,nconmax),rhoc(ngdim0)
      character          chlog,filenm*80,chlin1*78,chlin2*78

      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common             /rhoc/rhoc
      common             /chlin1/chlin1/chlin2/chlin2
      common             /conc/conc

c ----------------------------------------------------------------------
        write(*,*)
        write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%'      
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%       HOFFMAN-RIBAK CONSTRAINED GAUSSIAN RANDO',
     2             'M FIELD CODE        %' 
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  Program for setting up an initial Gaussian ra',
     2             'ndom density- and   %'
        write(*,*) '%  velocity field subject to one or more peak co',
     2             'nstraints (at the   %'
        write(*,*) '%  moment maximally 5 peaks).                   ',
     2             '                    %'
        write(*,*) '%  The statistical properties of the field are d',
     2             'etermined by its    %'
        write(*,*) '%  power spectrum P(k). Several choices of power',
     2             'spectrum are        %'
        write(*,*) '%  available, namely CDM, HDM and power law spec',
     2             'tra, P(k) ~ k^n.    %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  Each constrained peak or dip in the field is ',
     2             'characterized by a  %'
        write(*,*) '%  position {vec r} and a (Gaussian smoothing)  ',
     2             'scale R_G.          %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  The following constraints can be put on each ',
     2             'of these peaks:     %'
        write(*,*) '%  (1) height of peak                           ',
     2             '                    %'
        write(*,*) '%  (2) 1st derivatives equal to 0.0             ',
     2             '                    %'
        write(*,*) '%  (3) curvature of peak                        ',
     2             '                    %'
        write(*,*) '%  (4) two axis ratios of peak                  ',
     2             '                    %'
        write(*,*) '%  (5) orientation of peak (3 Euler angles)     ',
     2             '                    %'
        write(*,*) '%  (6) velocity of peak:                        ',
     2             '                    %'
        write(*,*) '%      (a) amplitude of velocity                ',
     2             '                    %'
        write(*,*) '%      (b) direction of velocity                ',
     2             '                    %'
        write(*,*) '%  (7) shear at peak:                           ',
     2             '                    %'
        write(*,*) '%      (a) two independent eigenvalues          ',
     2             '                    %'
        write(*,*) '%      (b) direction principal axes (3 Euler ang',
     2             'les)                %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  After having specified the constraints and th',
     2             'eir values, a       %'
        write(*,*) '%  realization of a random Gaussian field with t',
     2             'he desired power    %'
        write(*,*) '%  spectrum obeying these constraints is generat',
     2             'ed by means of the  %'
        write(*,*) '%  Hoffman-Ribak algorithm (Ap.J.Lett. 380, L5-L',
     2             '8, Oct. 10, 1991).  %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  Program based on original implementation of a',
     2             'lgorithm by         %'
        write(*,*) '%  E. Bertschinger (MIT), 29 July 1991.         ',
     2             '                    %'
        write(*,*) '%  Peak constraint implementation and further ad',
     2             'aptation:           %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%            Rien van de Weygaert,              ',
     2             '                    %'      
        write(*,*) '%            Canadian Institute for Theoretical ',
     2             'Astrophysics,       %'
        write(*,*) '%                     Toronto, Canada           ',
     2             '                    %'
        write(*,*) '%            Sterrewacht Leiden,                ',
     2             '                    %'
        write(*,*) '%                     Leiden, The Netherlands   ',
     2             '                    %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%  1st version:   March 1992                    ',
     2             '                    %'
        write(*,*) '%  2nd version:   April 1993                    ',
     2             '                    %'
        write(*,*) '%                                               ',
     2             '                    %'
        write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%'      
        write(*,*)
c ----------------------------------------------------------------------


        do 10 m=1,78
          chlin1(m:m)='='
          chlin2(m:m)='-'
10      continue

        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)

        write(*,'('' Do you want to have a log file [y/n]: '',$)')
        read(*,'(a)') chlog
        if ((chlog.eq.'y').or.(chlog.eq.'Y')) then
          write(*,'('' Name of log file: '',$)')
          read(*,'(a80)') filenm
          open(unit=1,file=filenm,status='new',err=5000)
        else
          open(unit=1,status='scratch',err=5000)
        endif

        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)

c ----------------------------------------------------------------------
        write(1,*)
        write(1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%'      
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%       HOFFMAN-RIBAK CONSTRAINED GAUSSIAN RANDO',
     2             'M FIELD CODE        %' 
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  Program for setting up an initial Gaussian ra',
     2             'ndom density- and   %'
        write(1,*) '%  velocity field subject to one or more peak co',
     2             'nstraints (at the   %'
        write(1,*) '%  moment maximally 5 peaks).                   ',
     2             '                    %'
        write(1,*) '%  The statistical properties of the field are d',
     2             'etermined by its    %'
        write(1,*) '%  power spectrum P(k). Several choices of power',
     2             'spectrum are        %'
        write(1,*) '%  available, namely CDM, HDM and power law spec',
     2             'tra, P(k) ~ k^n.    %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  Each constrained peak or dip in the field is ',
     2             'characterized by a  %'
        write(1,*) '%  position {vec r} and a (Gaussian smoothing)  ',
     2             'scale R_G.          %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  The following constraints can be put on each ',
     2             'of these peaks:     %'
        write(1,*) '%  (1) height of peak                           ',
     2             '                    %'
        write(1,*) '%  (2) 1st derivatives equal to 0.0             ',
     2             '                    %'
        write(1,*) '%  (3) curvature of peak                        ',
     2             '                    %'
        write(1,*) '%  (4) two axis ratios of peak                  ',
     2             '                    %'
        write(1,*) '%  (5) orientation of peak (3 Euler angles)     ',
     2             '                    %'
        write(1,*) '%  (6) velocity of peak:                        ',
     2             '                    %'
        write(1,*) '%      (a) amplitude of velocity                ',
     2             '                    %'
        write(1,*) '%      (b) direction of velocity                ',
     2             '                    %'
        write(1,*) '%  (7) shear at peak:                           ',
     2             '                    %'
        write(1,*) '%      (a) two independent eigenvalues          ',
     2             '                    %'
        write(1,*) '%      (b) direction principal axes (3 Euler ang',
     2             'les)                %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  After having specified the constraints and th',
     2             'eir values, a       %'
        write(1,*) '%  realization of a random Gaussian field with t',
     2             'he desired power    %'
        write(1,*) '%  spectrum obeying these constraints is generat',
     2             'ed by means of the  %'
        write(1,*) '%  Hoffman-Ribak algorithm (Ap.J.Lett. 380, L5-L',
     2             '8, Oct. 10, 1991).  %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  Program based on original implementation of a',
     2             'lgorithm by         %'
        write(1,*) '%  E. Bertschinger (MIT), 29 July 1991.         ',
     2             '                    %'
        write(1,*) '%  Peak constraint implementation and further ad',
     2             'aptation:           %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%            Rien van de Weygaert,              ',
     2             '                    %'      
        write(1,*) '%            Canadian Institute for Theoretical ',
     2             'Astrophysics,       %'
        write(1,*) '%                     Toronto, Canada           ',
     2             '                    %'
        write(1,*) '%            Sterrewacht Leiden,                ',
     2             '                    %'
        write(1,*) '%                     Leiden, The Netherlands   ',
     2             '                    %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%  1st version:   March 1992                    ',
     2             '                    %'
        write(1,*) '%  2nd version:   April 1993                    ',
     2             '                    %'
        write(1,*) '%                                               ',
     2             '                    %'
        write(1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%'      
        write(1,*)

c ------------------------------------------------------------------------

        write(*,*)
        write(*,'('' What kind of action do you want of me: '')')
        write(*,*)
	write(*,*) '  ido=1:  Compute unconstrained realization'
	write(*,*) '  ido=2:  Compute mean field of constraints'
	write(*,*) '  ido=3:  Compute constrained realization'
        write(*,*)
        write(*,'('' Enter ido: '',$)')
	read(*,*) ido
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'('' What kind of action do you want of me: '')')
        write(1,*)
	write(1,*) '  ido=1:  Compute unconstrained realization'
	write(1,*) '  ido=2:  Compute mean field of constraints'
	write(1,*) '  ido=3:  Compute constrained realization'
        write(1,*)
        write(1,'('' Enter ido: '',$)')
	write(1,'(i2)') ido
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)
	if (ido.lt.1.or.ido.gt.3) stop

c ------------------------------------------------------------------------

c I.    Some general cosmological parameters:

        call cosmin
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c ------------------------------------------------------------------------

c II.   Set up lattice

        call latini
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c ------------------------------------------------------------------------

c III.  Initializing random number generator

        call randini
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c ------------------------------------------------------------------------

c IV.   Set up power spectrum.

        call spectrum
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c ----------------------------------------------------------------------

c V.    Set up unconstrained field sample and print some information
c       on this unconstrained field

        call ucfldg(sigma,sigmav,chisq,nharm)
        call ucfldd(sigma,sigmav,chisq,nharm)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)
        if (ido.eq.1) goto 4000


c ----------------------------------------------------------------------

c VI.   Input constraints

c       ---------------------------------------------------------------------
c       Constraints are of the form g(i)=sum from j=1 to 2*ngr2 of 
c       con(j,i)*rho(j).
c       The constraints are Fourier transformed so that the convolutions with 
c       the covariance matrix required to compute the mean field and the 
c       covariance of the constraints becomes multiplication by the power 
c       spectrum.
c       ---------------------------------------------------------------------

        call constr(g,ncon)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c -----------------------------------------------------------------------

c VII.  Comparison sampled and desired constraints

        icmp=1
        call concmp(g,g0,ncon,icmp)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)
        
c ------------------------------------------------------------------------

c VIII. (1) Compute covariance matrix of constraints (\xi_ij in 
c           Hoffman-Ribak), and invert this matrix.
c       (2) Calculate chi-squared for the constraints and print this 
c           information.

        call covmat(qinv,ncon)
        call conchs(qinv,g,g0,ncon)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c -------------------------------------------------------------------------

c IX.   Computing either mean field or constrained realization.
          
        call fldprp(ido,g0,ncon)
        call fldcal(qinv,g,g0,ncon)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c -------------------------------------------------------------------------

c X.    Check the values of the constraints for the constrained 
c       realization.

        icmp=2
        call concmp(g,ghelp,ncon,icmp)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)

c -----------------------------------------------------------------------

c XI.   (1) Computing displacement psi-field in Fourier-space and 
c           Fouriertransforming psi and rho to spatial domain.
c       (2) Print information on the resulting density and velocity 
c           field.
c       (3) Write the rho and psi field to file.

4000    call psical(sigma,sigmav,chisq,nharm)
        call psiinf(sigma,sigmav,chisq,nharm,ncon,ido)
        write(*,*)
        write(*,'(1x,a78)') chlin1
        write(*,*)
        write(1,*)
        write(1,'(1x,a78)') chlin1
        write(1,*)
        call psifil

c -----------------------------------------------------------------------

c VIII. Closing the program

        call jimcls


5000   stop
      end


c     *************************************************************************
c     *************************************************************************


      subroutine cosmin


c     ------------------------------------
c     Initializing cosmological parameters
c     ------------------------------------

      implicit double precision(a-h,o-z)

      external fpeebl

      double precision  h0,hubbl0,omega0,omegahh,hf0

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0


c    Hubble constant (in units of 100 km/s/Mpc).

        write(*,'('' Cosmological quantities: '')')
        write(*,*)
        write(*,'('' Omega:                         '',$)')
        read(*,*) omega0
        write(*,'('' Hubble parameter h0:           '',$)')
        read(*,*) h0
        write(1,'('' Cosmological quantities: '')')
        write(1,*)
        write(1,'('' Omega:                         '',$)') 
        write(1,'(f10.3)') omega0
        write(1,'('' Hubble parameter h0:           '',$)')
        write(1,'(f10.3)') h0

        hubbl0=h0*100.0d0
        omegahh=omega0*(h0**2)
        hf0=hubbl0*fpeebl(omega0)

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine latini


c     ------------------
c     Setting up lattice
c     ------------------


      integer            n1,n12,n21,n1n1,ngrid,n11,ngr2,
     2                   ngdim2,ngdim
      real               boxlen
      double precision   h0,hubbl0,omega0,omegahh,hf0

      common             /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common             /box/boxlen
      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


        write(*,'('' Box/Lattice information: '')')
        write(*,*)
        write(*,'('' Comoving size of box (h^{-1} Mpc): '',$)')
        read(*,*) boxlen
        write(*,*)
        write(1,'('' Box/Lattice information: '')')
        write(1,*)
        write(1,'('' Comoving size of box (h^{-1} Mpc): '',$)') 
        write(1,'(f10.3)') boxlen
        write(1,*)

c    Convert boxlen to Mpc
        boxlen=boxlen/h0

        write(*,'('' How many grid points in each direction ?: '')')
        write(*,'('' power of 2:                           '',$)')
        read(*,*) n1
        write(1,'('' How many grid points in each direction ?: '')')
        write(1,'('' power of 2:                           '',$)')
        write(1,'(i6)') n1

        n12=n1/2
        n21=n12+1
        n1n1=n1*n1
        ngrid=n1n1*n1
        n11=n1-1
        ntot=n1n1*n12
        ngr2=ngrid/2
        ngdim2=ngr2+n1n1
        ngdim=2*ngdim2
        write(*,*) 'Number of grid points:      ',2*ntot
        write(1,'(1x,a22,i22)') 'Number of grid points:',2*ntot

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine ucfldg(sigma,sigmav,chisq,nharm)


c     -------------------------------------------
c     Generating unconstrained field realization.
c     -------------------------------------------

      implicit double precision (a-h,o-z)

      double precision   twopi
      parameter          (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter          (twopi=6.283185307179586d0,sqr2=1.41421356)
      parameter          (n0n0=n0*n0,ngrid0=n0n0*n0,ngr0=ngrid0/2)
      parameter          (ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer            n1,n12,n21,n1n1,ngrid,n11,ngr2,
     2                   ngdim2,ngdim,nharm
      double precision   sigma,sigmav,chisq
      double precision   power
      real               xr,boxlen,cosarg,sinarg
      complex            z
      complex            rhoc

      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common             /box/boxlen
      common             /power/power(ngdim0)           
      common             /rhoc/rhoc(ngdim0)


        dx=boxlen/n1
        dk=twopi/boxlen
        d3k=dk*dk*dk
        akmax=n1*dk

	sigma=0.0
	sigmav=0.0
	chisq=0.0
	nharm=0

c       ----------------------------------------------------------------
c       Generate unconstrained sample of rho in Fourier transform space.
c       ----------------------------------------------------------------
	do 30 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          ak33=ak3*ak3
	  do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax
            ak23=ak2*ak2+ak33
  	    k23=(k2-1+(k3-1)*n1)*n12

	    do 10 k1=2,n12
c             --------------------------------------
c             Do k1=1 and k1=n12+1 separately below.
c             --------------------------------------
              ak1=(k1-1)*dk
  	      k=k1+k23
              akk=ak1*ak1+ak23
	      sigft=sqrt(power(k))
	      call randa(xr)
	      arg=twopi*xr
              cosarg=dcos(arg)
              sinarg=dsin(arg)
	      z=cmplx(cosarg,sinarg)
5	      call randa(xr)
	      if (xr.eq.0.0) go to 5
c             ----------------------------------------------------------------
c             Reduce z by sqrt(2) for real f since need sum of |z|^2 over all
c             wavenumbers (positive and negative) to be Chi-square for N 
c             degrees of freedom, not 2N as in the complex case. Equivalently,
c             the power is shared with the conjugate harmonics with k1 > n12+1
c             (k1 < 0).
c             ----------------------------------------------------------------
	      z=sqrt(-log(xr))*z
              rhoc(k)=z*sigft
c             -------------------------------------------------
c             Double the contribution to account for modes with 
c             k1 > n12+1 (k1 < 0).
c             -------------------------------------------------
	      sigma=sigma+2.0*power(k)
	      sigmav=sigmav+2.0*power(k)/akk
	      chisq=chisq+2.0*conjg(z)*z
	      nharm=nharm+2
10	    continue

	    k23=k2-1+(k3-1)*n1
	    k23c=n1+1-k2+(n1+1-k3)*n1
	    if (k2.eq.1) k23c=k23c-n1
            if (k3.eq.1) k23c=k23c-n1n1
	    if (k23.gt.k23c) go to 20

	    do 15 i=1,2
c             ---------------------
c             Do k1=1 and k1=n12+1.
c             ---------------------
	      if (i.eq.1) then
                k=1+n12*k23
		kc=1+n12*k23c
		ak1=0.0
	      else
		k=ngr2+1+k23
	  	kc=ngr2+1+k23c
		ak1=-0.5*akmax
	      end if
	      akk=ak1*ak1+ak23
	      sigft=sqrt(power(k))
	      call randa(xr)
	      arg=twopi*xr
	      if (k23.ne.k23c) then
                cosarg=dcos(arg)
                sinarg=dsin(arg)
	        z=cmplx(cosarg,sinarg)
	      else
                cosarg=dcos(arg)
	        z=sqr2*cmplx(cosarg,0.0)
	      end if
12	      call randa(xr)
	      if (xr.eq.0.0) go to 12
	      z=sqrt(-log(xr))*z
	      rhoc(k)=z*sigft
c             --------------------------------------------
c             Determine conjugate harmonic using symmetry.
c             --------------------------------------------
	      rhoc(kc)=conjg(rhoc(k))
	      if (k23.ne.k23c) then
	        sigma=sigma+2.0*power(k)
	        sigmav=sigmav+2.0*power(k)/akk
	        chisq=chisq+2.0*conjg(z)*z
	        nharm=nharm+2
	      else if (akk.ne.0.0) then
	        sigma=sigma+power(k)
	        chisq=chisq+conjg(z)*z
	        nharm=nharm+1
	      end if
15	    continue
20	  continue
30	continue
	rhoc(1)=cmplx(0.0,0.0)

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine ucfldd(sigma,sigmav,chisq,nharm)


c     -----------------------------------------------------
c     Writing some data on unconstrained field realization.
c     -----------------------------------------------------

      implicit double precision (a-h,o-z)

      external           fpeebl

      integer            nharm
      double precision   sigma,sigmav,chisq
      double precision   h0,hubbl0,omega0,omegahh,hf0

      common             /cosmos/h0,hubbl0,omega0,omegahh,hf0

	sigma=sqrt(sigma)
	sigmav=sqrt(sigmav)*hubbl0*fpeebl(omega0)
        write(*,*) 'Density and velocity data of',
     2    ' the unconstrained field: '
        write(*,*)
	write(*,'(1x,a22,f15.6)') ' sigma              = ',real(sigma)
        write(*,'(1x,a22,f15.6,a5)') ' sigma_v            = ',
     2    real(sigmav),' km/s'
	anu=(chisq-nharm)/sqrt(2.0*nharm)
        write(*,'(1x,a22,f15.6)') ' Chi-squared        = ',real(chisq)
        write(*,'(1x,a22,i8)') ' degrees of freedom = ',nharm
        write(*,'(1x,a22,f15.6)') ' nu                 = ',anu
        write(1,*) 'Density and velocity data of',
     2    ' the unconstrained field: '
        write(1,*)
	write(1,'(1x,a22,f15.6)') ' sigma              = ',real(sigma)
        write(1,'(1x,a22,f15.6,a5)') ' sigma_v            = ',
     2    real(sigmav),' km/s'
        write(1,'(1x,a22,f15.6)') ' Chi-squared        = ',real(chisq)
        write(1,'(1x,a22,i8)') ' degrees of freedom = ',nharm
        write(1,'(1x,a22,f15.6)') ' nu                 = ',anu

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine fldprp(ido,g0,ncon)


c     ----------------------------------------------------------------------
c     Field calculation preparation routine. If the full constrained field 
c     needs to be calculated (in fldcal) nothing needs to be done, but if 
c     mean field needs only to be computed the content of the arrays rhoc(k) 
c     and g0(k) should be put to 0.0. 
c     In this way for both cases the routine fldcal can be called without 
c     further specifications.
c     ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter         (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ido,ncon
      real              g0(nconmax)
      complex           rhoc      

      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /rhoc/rhoc(ngdim0)


	if (ido.eq.2) then
c         -------------------------------------------------------------
c         If computing mean field only put:
c         unconstrained realization rhoc(k)=0.0
c         values of constraints for unconstrained realization g0(i)=0.0
c         -------------------------------------------------------------
	  write(*,*) 'Computing mean field only'
	  write(1,*) 'Computing mean field only'
	  do 105 k=1,ngdim2
	    rhoc(k)=cmplx(0.0,0.0)
105       continue
	  do 106 i=1,ncon
	    g0(i)=0.0
106       continue
	else
c         ------------------------------------------
c         If calculating the full constrained field.
c         ------------------------------------------
	  write(*,*) 'Computing a constrained realization'
	  write(1,*) 'Computing a constrained realization'
	endif

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine fldcal(qinv,g,g0,ncon)


c     ---------------------------------------------------------------
c     Compute corrected density field 
c     The computation is the same for both mean field and constrained 
c     realization because rhoc(k) and g0(i) have been set to the 
c     appropriate values in both cases (done in fldprp).
c     This is calculating in Fourier space the expression
c       f(r)=f~(r)+\xi_i(r) \xi_ij^{-1} (c_j-c~_j)
c     in Hoffman-Ribak.
c     ---------------------------------------------------------------

      implicit double precision(a-h,o-z)

      parameter         (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              g(nconmax),g0(nconmax)
      double precision  qinv(nconmax,nconmax)
      double precision  power
      complex           conc,rhoc

      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /rhoc/rhoc(ngdim0)
      common            /conc/conc(ngdim0,nconmax)
      common            /power/power(ngdim0)


	do 130 i2=1,ncon
	  do 120 i1=1,ncon
	    fact=(g(i2)-g0(i2))*qinv(i1,i2)
	    do 110 k=1,ngdim2
	      rhoc(k)=rhoc(k)+fact*conc(k,i1)*power(k)
110	    continue
120	  continue
130	continue

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine psical(sigma,sigmav,chisq,nharm)

c     ----------------------------------------------------------------------

c     ftpsi finishes the job of hofrib by computing the displacement 
c     field and Fourier transforming density and displacement to bring 
c     them into the spatial domain.
c
c     To get starting positions and velocities for an N-body simulation,
c     psi has to be multiplied by the density resp. velocity growth 
c     factor (both functions of the expansion factor a and present 
c     parameter omega of the Universe, program INITEB can carry out this 
c     job), which are just equal to the scale factor a=1/(1+z) if Omega=1.
c     Then the comoving position of particle (i,j,k) is given by:
c         x=(i-1)*dx+grow*psi1(m)
c         y=(j-1)*dx+grow*psi2(m)
c         z=(k-1)*dx+grow*psi3(m)
c     where m=i+n1*(j-1)+n1*n1*(k-1), dx=boxlen/n1, and grow a function 
c     of the expansion factor aexp and omega: grow=grow(aexp,omega).
c     The peculiar velocity is given by H*aexp*grow*f(omega)*(psi,psi2,psi3)
c     where H=(da/dt)/a is the Hubble parameter at aexp(t), grow is the 
c     density growth factor, depending on the omega at aexp(t), and 
c     f(omega) is the velocity growth factor (see Peebles LSS, par. 14), 
c     a function of omega at aexp(t).
c     rhoc=delta-rho/rho and is dimensionless; psi has units of Mpc. Both 
c     are evaluated at aexp=1 (z=0).
c     ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      external  fpeebl 
      
      double precision   twopi
      parameter          (n0=128,twopi=6.283185307179586d0)
      parameter          (n0n0=n0*n0,ngr0=n0n0*n0/2,ngdim0=ngr0+n0n0)

      integer            nharm
      real               boxlen
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma,sigmav,chisq
      double precision   power
      complex            z
      complex            rhoc,psi1,psi2,psi3

      common             /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common             /box/boxlen
      common             /rhoc/rhoc(ngdim0)
      common             /psi/psi1(ngdim0),psi2(ngdim0),psi3(ngdim0)
      common             /power/power(ngdim0)
      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


        dk=twopi/boxlen
        akmax=float(n1)*dk
        n12=n1/2
        ngr2=n1*n1*n12

        sigma=0.0
        sigmav=0.0
        chisq=0.0
        nharm=0

c       ---------------------------
c       Compute displacement field:
c       ---------------------------

	psi1(1)=cmplx(0.0,0.0)
	psi2(1)=cmplx(0.0,0.0)
	psi3(1)=cmplx(0.0,0.0)

	do 190 k3=1,n1
	  ak3=(k3-1)*dk
	  if (k3.gt.n12) ak3=ak3-akmax
	  ak33=ak3*ak3
	  do 180 k2=1,n1
	    ak2=(k2-1)*dk
	    if (k2.gt.n12) ak2=ak2-akmax
	    ak23=ak2*ak2+ak33
	    k23=(k2-1+(k3-1)*n1)*n12

c           --------------------------------------
c           Do k1=1 and k1=n12+1 separately below.
c           --------------------------------------
	    do 170 k1=2,n12
	      ak1=(k1-1)*dk
	      k=k1+k23
	      akk=ak1*ak1+ak23
c             ------------------------------------------
c             Calculate displacement in units of boxlen.
c             ------------------------------------------
	      z=cmplx(0.0,1.0)*rhoc(k)/akk/boxlen
	      psi1(k)=ak1*z
	      psi2(k)=ak2*z
	      psi3(k)=ak3*z
	      if (k2.eq.n21) psi2(k)=cmplx(0.0,0.0)
	      if (k3.eq.n21) psi3(k)=cmplx(0.0,0.0)

	      dsig=2.0*conjg(rhoc(k))*rhoc(k)
	      sigma=sigma+dsig
	      dsigv=conjg(psi1(k))*psi1(k)+conjg(psi2(k))*psi2(k)
     2             +conjg(psi3(k))*psi3(k)
	      sigmav=sigmav+2.0*dsigv
	      chisq=chisq+dsig/power(k)
	      nharm=nharm+2
170	    continue

	    k23=k2-1+(k3-1)*n1
	    k23c=n1+1-k2+(n1+1-k3)*n1
	    if (k2.eq.1) k23c=k23c-n1
	    if (k3.eq.1) k23c=k23c-n1n1
	    if (k23.gt.k23c) go to 180

c           ---------------------
c           Do k1=1 and k1=n12+1.
c           ---------------------
	    do 175 i=1,2
	      if (i.eq.1) then
		k=1+n12*k23
		kc=1+n12*k23c
		ak1=0.0
	      else
		k=ngr2+1+k23
		kc=ngr2+1+k23c
		ak1=-0.5*akmax
	      end if
	      akk=ak1*ak1+ak23
	      if (k23.ne.k23c) then
	        z=cmplx(0.0,1.0)*rhoc(k)/akk/boxlen
	        dsig=2.0*conjg(rhoc(k))*rhoc(k)
	        sigma=sigma+dsig
		chisq=chisq+dsig/power(k)
		nharm=nharm+2
	      else if (akk.ne.0.0) then
	        z=cmplx(0.0,0.0)
		dsig=conjg(rhoc(k))*rhoc(k)
	        sigma=sigma+dsig
		chisq=chisq+dsig/power(k)
		nharm=nharm+1
	      end if
	      psi1(k)=cmplx(0.0,0.0)
	      psi2(k)=ak2*z
	      psi3(k)=ak3*z
	      if (k2.eq.n21) psi2(k)=cmplx(0.0,0.0)
	      if (k3.eq.n21) psi3(k)=cmplx(0.0,0.0)

c             --------------------------------------------
c             Determine conjugate harmonic using symmetry.
c             --------------------------------------------
	      psi1(kc)=conjg(psi1(k))
	      psi2(kc)=conjg(psi2(k))
	      psi3(kc)=conjg(psi3(k))
	      dsigv=psi1(kc)*psi1(k)+psi2(kc)*psi2(k)
     2             +psi3(kc)*psi3(k)
	      sigmav=sigmav+2.0*dsigv
175	    continue

180	  continue
190	continue

	sigma=sqrt(sigma)
	sigmav=sqrt(sigmav)*hubbl0*fpeebl(omega0)*boxlen

c       -------------------------------------------------
c       Fourier transform density and displacement fields 
c       to position space.
c       -------------------------------------------------
	call fft3rinv(rhoc,n1)
	call fft3rinv(psi1,n1)
	call fft3rinv(psi2,n1)
	call fft3rinv(psi3,n1)

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine psiinf(sigma,sigmav,chisq,nharm,ncon,ido)


c     ------------------------------------------------------------------
c     Calculating maximum density and displacement in field and printing 
c     this information, together with some other info.
c     ------------------------------------------------------------------

      implicit double precision(a-h,o-z)

      external fpeebl

      double precision   twopi
      parameter          (n0=128,twopi=6.283185307179586d0)
      parameter          (n0n0=n0*n0,ngr0=n0n0*n0/2,ngdim0=ngr0+n0n0)

      integer            nharm,ncon,ido
      real               boxlen
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma,sigmav,chisq
      complex            rhoc,psi1,psi2,psi3

      common             /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common             /box/boxlen
      common             /rhoc/rhoc(ngdim0)
      common             /psi/psi1(ngdim0),psi2(ngdim0),psi3(ngdim0)
      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


c       -------------------------------------------   
c       Determine maximum density and displacement.
c       -------------------------------------------
        delm=0.0
        dispm=0.0
        do 200 j=1,ngr2
          del=abs(real(rhoc(j)))
          delm=max(delm,del)
          del=abs(aimag(rhoc(j)))
          delm=max(delm,del)
          ddisp=real(psi1(j))**2+real(psi2(j))**2+real(psi3(j))**2
          dispm=max(dispm,ddisp)
          ddisp=aimag(psi1(j))**2+aimag(psi2(j))**2+aimag(psi3(j))**2
          dispm=max(dispm,ddisp)
200     continue
        dispm=sqrt(dispm)*hubbl0*fpeebl(omega0)*boxlen

c       ------------------------------------------------
c       Print information on density and velocity field.
c       ------------------------------------------------
        if (ido.gt.1) then
          if (ido.eq.2) then
            write(*,*) 'Density and velocity data of',
     2        ' the mean constrained field: '
            write(1,*) 'Density and velocity data of',
     2        ' the mean constrained field: '
          else
            write(*,*) 'Density and velocity data of',
     2        ' the constrained field realization: '
            write(1,*) 'Density and velocity data of',
     2        ' the constrained field realization: '
          endif
          write(*,*)
          write(*,'(1x,a24,f15.6)') ' sigma                = ',
     2      real(sigma)
          write(*,'(1x,a24,f15.6)') ' maximum delta        = ',delm
          write(*,'(1x,a24,f15.6,a5)') ' sigma_v              = ',
     2      real(sigmav),' km/s'
          write(*,'(1x,a24,f15.6,a5)') ' maximum displacement = ',
     2      dispm,' km/s'
          write(*,*)
          write(*,'(1x,a24,f15.6)') ' Chi-squared          = ',
     2      real(chisq)
          write(*,'(1x,a24,i8)') ' degrees of freedom   = ',
     2      nharm-ncon
          write(1,*)
          write(1,'(1x,a24,f15.6)') ' sigma                = ',
     2      real(sigma)
          write(1,'(1x,a24,f15.6)') ' maximum delta        = ',delm
          write(1,'(1x,a24,f15.6,a5)') ' sigma_v              = ',
     2      real(sigmav),' km/s'
          write(1,'(1x,a24,f15.6,a5)') ' maximum displacement = ',
     2      dispm,' km/s'
          write(1,*)
          write(1,'(1x,a24,f15.6)') ' Chi-squared          = ',
     2      real(chisq)
          write(1,'(1x,a24,i8)') ' degrees of freedom   = ',
     2      nharm-ncon
        end if

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine psifil


c     ----------------------------------------------
c     Write density and displacement fields to file.
c     ----------------------------------------------

      implicit double precision(a-h,o-z)

      double precision   twopi
      parameter          (n0=128,twopi=6.283185307179586d0)
      parameter          (n0n0=n0*n0,ngr0=n0n0*n0/2,ngdim0=ngr0+n0n0)

      real               boxlen
      double precision   h0,hubbl0,omega0,omegahh,hf0
      complex            rhoc,psi1,psi2,psi3
      character          filename*80

      common             /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common             /box/boxlen
      common             /rhoc/rhoc(ngdim0)
      common             /psi/psi1(ngdim0),psi2(ngdim0),psi3(ngdim0)
      common             /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


c       -------------
c       Opening file:
c       -------------
        write(*,*)
        write(*,'('' Filename for density- and psi-field: '',$)')
        read(*,'(a80)') filename
        write(1,*)
        write(1,'('' Filename for density- and psi-field: '',$)')
        write(1,'(a80)') filename
        open(10,file=filename,status='new',form='unformatted')
        rewind 10

c       --------------------------------------------------------
c       p=a**alpha will be the time variable in the N-body code.
c       --------------------------------------------------------
	alpha=1.0
        aexp=1.0
	npart=n1**3

c       -----------------------
c       Writing header of file:
c       -----------------------
        write(*,*)
        write(*,'(1x,a9,i8)') ' n1:     ',n1
        write(*,'(1x,a9,i8)') ' npart:  ',npart
        write(*,'(1x,a9,f12.3,a5)') ' boxlen: ',REAL(hubbl0*boxlen),
     2                        ' km/s'
        write(*,'(1x,a9,f12.3)') ' Omega:  ',REAL(omega0)
        write(*,'(1x,a9,f12.3)') ' H0:     ',REAL(hubbl0)
        write(*,'(1x,a9,f12.3)') ' alpha:  ',REAL(alpha)
        write(10)  n1,npart,REAL(hubbl0*boxlen),REAL(aexp),REAL(omega0),
     2             REAL(hubbl0),REAL(alpha)      
        write(*,*)
        write(1,*)
        write(1,'(1x,a9,i8)') ' n1:     ',n1
        write(1,'(1x,a9,i8)') ' npart:  ',npart
        write(1,'(1x,a9,f12.3,a5)') ' boxlen: ',REAL(hubbl0*boxlen),
     2                        ' km/s'
        write(1,'(1x,a9,f12.3)') ' Omega:  ',REAL(omega0)
        write(1,'(1x,a9,f12.3)') ' H0:     ',REAL(hubbl0)
        write(1,'(1x,a9,f12.3)') ' alpha:  ',REAL(alpha)
        write(1,*)

c       -----------------------------------------------
c       Writing density and displacement field to file:
c       -----------------------------------------------
        do 210 j=1,ngr2
          write(10) real(rhoc(j)),
     2              real(psi1(j)),real(psi2(j)),real(psi3(j))
          write(10) aimag(rhoc(j)),
     2              aimag(psi1(j)),aimag(psi2(j)),aimag(psi3(j))
210     continue                

c       -------------
c       Closing file.
c       -------------
	close(10)

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine constr(g,ncon)


c     --------------------------------------------------------------------
c     constr3 provides the constraint matrix for hofmrib.f.  
c     con(j1,j2,j3,i)=con((i-1)*n1*n1*n1+j), j=j1+(j2-1)*n1+(j3-1)*n1*n1, 
c     is defined so that the constraints are
c		g(i)=sum from j=1 to n1*n1*n1 of con(j,i)*f(j),
c     for i=1 to ncon.  E.g., if n=64 and there are two constraints
c     f(5)=2, f(7)=3,  then con(5,1)=con(5)=2, con(7,2)=con(71)=3, and all
c     other values of con vanish.
c     --------------------------------------------------------------------


      implicit double precision (a-h,o-z)

      parameter         (npeakm=3,nconmax=18*npeakm)

      real              boxlen
      real              g(nconmax)
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  cenpk(3),reuler(3,3),lambda(3),lamb21,lamb31
      character         chlin1*78,chlin2*78,answer*1
      logical           lgcon(5),lpkqst(npeakm,7),lfsecder(npeakm)


      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /box/boxlen
      common /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim 
      common /chlin1/chlin1/chlin2/chlin2


        write(*,*)
        write(*,'('' CONSTRAINT INFORMATION: '')')
        write(*,'('' ======================'')')
        write(*,*)
        write(*,'('' The constraints concern peaks or dips in a'',$)')
        write(*,'('' Gaussian random field. '')')
        write(*,'('' In addition to the specific properties of'',$)')
        write(*,'('' extremum (peak or dips) '')')
        write(*,'('' the following general information has to be'',$)')
        write(*,'('' provided first: '')')
        write(*,'('' [1] the number of constrained peaks/dips '')')
        write(*,'('' [2] peak or dip '')')
        write(*,'('' [3] the Gaussian smoothing scale on which '',$)')
        write(*,'(''it is defined '')')
        write(*,'('' [4] the position of the peak centre '')')
        write(*,*)
        write(1,*)
        write(1,'('' CONSTRAINT INFORMATION: '')')
        write(1,'('' ======================'')')
        write(1,*)
        write(1,'('' The constraints concern peaks or dips in a'',$)')
        write(1,'('' Gaussian random field. '')')
        write(1,'('' In addition to the specific properties of'',$)')
        write(1,'('' extremum (peak or dips) '')')
        write(1,'('' the following general information has to be'',$)')
        write(1,'('' provided first: '')')
        write(1,'('' [1] the number of constrained peaks/dips '')')
        write(1,'('' [2] peak or dip '')')
        write(1,'('' [3] the Gaussian smoothing scale on which '',$)')
        write(1,'(''it is defined '')')
        write(1,'('' [4] the position of the peak centre '')')
        write(1,*)

        xmax=float(n1)

        write(*,*)
        write(*,'('' Length scales of constraints in units of '')')
        write(*,'(''  [1] boxlength '')')
        write(*,'(''  [2] h^{-1} Mpc '')')  
        write(*,'('' Unit system [1/2]: '',$)')
        read(*,*) iqu
        if (iqu.ne.1) iqu=2
        write(*,*)
        write(1,*)
        write(1,'('' Length scales of constraints in units of '')')
        write(1,'(''  [1] boxlength '')')
        write(1,'(''  [2] h^{-1} Mpc '')')  
        write(1,'('' Unit system [1/2]: '',$)')
        write(1,'(i3)') iqu
        write(1,*)

        write(*,'(1x,a14,f7.2,a13)') 'Boxlength:    ',h0*boxlen,
     2        ' h^{-1} Mpc.'
        write(*,*)
        write(1,'(1x,a14,f7.2,a13)') 'Boxlength:    ',h0*boxlen,
     2        ' h^{-1} Mpc.'
        write(1,*)

        npm=npeakm
        write(*,'('' Number of objects [1-'',i2,'']: '',$)') npm
        read(*,*) npeak
        write(*,*)
        write(1,'('' Number of objects [1-'',i2,'']: '',$)') npm
        write(1,'(i3)') npeak
        write(1,*)

c       --------------------------------------------------------
c       Part for reading in the peak/dip constraint information:
c       --------------------------------------------------------
        do 10 ipk=1,npeak

          write(*,'(1x,a78)') chlin1
          write(*,*)
          write(*,'(1x,a9,2x,i6,a2)') 'Object',ipk,': '
          write(*,'(1x,a9)') '======'
          write(1,'(1x,a78)') chlin1
          write(1,*)
          write(1,'(1x,a9,2x,i6,a2)') 'Object',ipk,': '
          write(1,'(1x,a9)') '======'

          call cspsgn(ipk)
          call cspscl(ipk,iqu)
          call csppos(ipk,iqu)
          call cspqst(ipk,lpkqst)
          call cspspi(ipk)

          lfsecder(ipk)=.false.
          lfsecder(ipk)=lfsecder(ipk).or.lpkqst(ipk,3)
          lfsecder(ipk)=lfsecder(ipk).or.lpkqst(ipk,4)
          lfsecder(ipk)=lfsecder(ipk).or.lpkqst(ipk,5)          
          if (lpkqst(ipk,1)) call csphgt(ipk)
          if (lpkqst(ipk,2)) call csgrad(ipk)
          if (lfsecder(ipk)) call cspcur(ipk,lpkqst)
          if (lfsecder(ipk)) call cspshp(ipk,lpkqst)
          if (lfsecder(ipk)) call csporn(ipk,lpkqst)
          if (lpkqst(ipk,6)) call cspvel(ipk)
          if (lpkqst(ipk,7)) call cspshr(ipk)

10      continue         

c       --------------------------------------------------------
c       Part for determining the values of the constraints g(j):
c       --------------------------------------------------------
        ipc=0

        do 20 ipk=1,npeak

          if (lpkqst(ipk,1)) then
            g(ipc+1)=pnu(ipk)
            ipc=ipc+1
c           -------------------
c           density of peak/dip
c           -------------------
          endif

          if (lpkqst(ipk,2)) then
            do 25 k=1,3
              g(ipc+k)=pgrd(k,ipk)
25          continue
            ipc=ipc+3
c           ---------------------------------
c           partial derivatives near object,
c           in case of peak/dip equal to zero
c           ---------------------------------
          endif

          if (lfsecder(ipk)) then
            lamb21=a12(ipk)*a12(ipk)
            lamb31=a13(ipk)*a13(ipk)
            lambda(1)=-px(ipk)/(1.0+lamb21+lamb31)
            lambda(2)=lamb21*lambda(1)
            lambda(3)=lamb31*lambda(1)
            call euler(reuler,phi(ipk),theta(ipk),psi(ipk))
            do 30 k=1,6
              g(ipc+k)=0.0
30          continue
            do 35 j=1,3
              g(ipc+1)=g(ipc+1)+lambda(j)*reuler(j,1)*reuler(j,1)
              g(ipc+2)=g(ipc+2)+lambda(j)*reuler(j,2)*reuler(j,2)
              g(ipc+3)=g(ipc+3)+lambda(j)*reuler(j,3)*reuler(j,3)
              g(ipc+4)=g(ipc+4)+lambda(j)*reuler(j,1)*reuler(j,2)
              g(ipc+5)=g(ipc+5)+lambda(j)*reuler(j,1)*reuler(j,3)
              g(ipc+6)=g(ipc+6)+lambda(j)*reuler(j,2)*reuler(j,3)      
35          continue
c           ------------------------------------------------------
c           second order partial derivatives, describing shape and 
c           orientation of density field near peak.
c           ------------------------------------------------------
            ipc=ipc+6
          endif 

          if (lpkqst(ipk,6)) then
            g(ipc+1)=vx(ipk)/sigmav(ipk)
            g(ipc+2)=vy(ipk)/sigmav(ipk)
            g(ipc+3)=vz(ipk)/sigmav(ipk)
            ipc=ipc+3
          endif
c         ----------------------------------        
c         the peculiar velocity of the peak.
c         ----------------------------------

          if (lpkqst(ipk,7)) then
            g(ipc+1)=sh1(ipk)/(hf0*sigma0(ipk))
            g(ipc+2)=sh2(ipk)/(hf0*sigma0(ipk))
            g(ipc+3)=sh3(ipk)/(hf0*sigma0(ipk))
            g(ipc+4)=sh4(ipk)/(hf0*sigma0(ipk))
            g(ipc+5)=sh5(ipk)/(hf0*sigma0(ipk))
            ipc=ipc+5
          endif

20      continue

c       -------------------------------------------------
c       Part for initializing constraint matrix con(j,i):
c       -------------------------------------------------
        ncon=0
        do 120 ipk=1,npeak
          do 125 m=1,5
            lgcon(m)=.false.
125       continue
          if (lpkqst(ipk,1)) lgcon(1)=.true.
          if (lpkqst(ipk,2)) lgcon(2)=.true.
          if (lfsecder(ipk)) lgcon(3)=.true.
          if (lpkqst(ipk,6)) lgcon(4)=.true.
          if (lpkqst(ipk,7)) lgcon(5)=.true.
          rgauss=pscale(ipk)
          cenpk(1)=centre(1,ipk)
          cenpk(2)=centre(2,ipk)
          cenpk(3)=centre(3,ipk)
          call conkmn(ncon,lgcon,rgauss,cenpk,ipk)
120     continue


        write(*,'(1x,a78)') chlin2
        write(*,*)
        write(*,*) 'End of constraint specification part.'
        write(1,'(1x,a78)') chlin2
        write(1,*)
        write(1,*) 'End of constraint specification part.'

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspqst(ipk,lpkqst)


      implicit double precision (a-h,o-z)

      parameter (npeakm=3)

      logical         lpkqst(npeakm,7)
      character       answer*1

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    


        do 10 n=1,7
          lpkqst(ipk,n)=.false.
10      continue

        write(*,*)
        write(*,'('' (*) Object parameters to be specified [y/n]:  '')')
        write(*,'(''     --------------------------------- '')')
        write(*,*)
        write(1,*)
        write(1,'('' (*) Object parameters to be specified [y/n]:  '')')
        write(1,'(''     --------------------------------- '')')
        write(1,*)
        write(*,'(''     [1] Object density:      '',$)')
        write(1,'(''     [1] Object density:      '',$)')
        read(*,'(a)') answer
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,1)=.true.
        write(*,'(''     [2] Density gradient:    '',$)')
        write(1,'(''     [2] Density gradient:    '',$)')
        read(*,'(a)') answer
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,2)=.true.
        write(*,'(''     [3] Curvature:           '',$)')
        read(*,'(a)') answer
        write(1,'(''     [3] Curvature:           '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,3)=.true.
        write(*,'(''     [4] Shape:               '',$)')
        read(*,'(a)') answer
        write(1,'(''     [4] Shape:               '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,4)=.true.
        write(*,'(''     [5] Orientation:         '',$)')
        read(*,'(a)') answer
        write(1,'(''     [5] Orientation:         '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,5)=.true.
        write(*,'(''     [6] Peculiar velocity:   '',$)')
        read(*,'(a)') answer
        write(1,'(''     [6] Peculiar velocity:   '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,6)=.true. 
        write(*,'(''     [7] Shear:               '',$)')
        read(*,'(a)') answer
        write(1,'(''     [7] Shear:               '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'y').or.(answer.eq.'Y')) lpkqst(ipk,7)=.true.        
        write(*,*)
        write(1,*)


       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspsgn(ipk)


      implicit double precision (a-h,o-z)

      parameter (npeakm=3)

      character       chpv*1

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    


        psign(ipk)=1.0d0
        write(*,*)
        write(*,'('' (*) Peak or Void ? '')')
        write(*,'(''     ------------'')')
        write(*,*)
        write(*,'(''      [p/v]: '',$)')
        read(*,'(a1)') chpv
        write(*,*)
        write(1,*)
        write(1,'('' (*) Peak or Void ? '')')
        write(1,'(''     ------------'')')
        write(1,*)
        write(1,'(''      [p/v]: '',$)')
        write(1,'(a1)') chpv
        write(1,*)
        if ((chpv.eq.'v').or.(chpv.eq.'V')) psign(ipk)=-1.0d0

       return
      end


C     *************************************************************************
C     *************************************************************************

      
      subroutine cspscl(ipk,iqu)


      implicit double precision (a-h,o-z)

      parameter         (npeakm=3)

      integer           ipk,iqu
      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /box/boxlen


        if (psign(ipk).lt.0.0) then
          write(*,'('' (*) Scale of void: '')')
          write(*,'(''     -------------'')')
          write(1,'('' (*) Scale of void: '')')
          write(1,'(''     -------------'')')
        else
          write(*,'('' (*) Scale of peak: '')') 
          write(*,'(''     -------------'')')
          write(1,'('' (*) Scale of peak: '')') 
          write(1,'(''     -------------'')')
        endif

1110    write(*,*)
        write(1,*)
        if (iqu.eq.1) then
          write(*,'(''      peakscale [0.1-0.4]: '',$)')
          read(*,*) pscale(ipk)
          write(1,'(''      peakscale [0.1-0.4]: '',$)')
          write(1,*) pscale(ipk)
        else
          write(*,'(''      peakscale [h^{-1} Mpc]: '',$)')
          read(*,*) pscale(ipk)
          write(1,'(''      peakscale [h^{-1} Mpc]: '',$)')
          write(1,'(f10.3)') pscale(ipk)
          pscale(ipk)=pscale(ipk)/(h0*boxlen)
          if (pscale(ipk).gt.0.4d0) then
            write(*,'(''     This is too large. '')')
            write(*,'(''     Try new scale: '')')
            write(1,'(''     This is too large. '')')
            write(1,'(''     Try new scale: '')')
            goto 1110
          endif
        endif
        write(*,*)
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspspi(ipk)


      implicit double precision (a-h,o-z)

      parameter         (npeakm=3)

      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  rgauss,s0,s1,s2,rsthm1,rstgs,gmpeak,rspeak

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /box/boxlen


        rgauss=pscale(ipk)*boxlen
        call spparam(rgauss,s0,s1,s2,gmpeak,rspeak)
        sigma0(ipk)=s0
        sigma1(ipk)=s1
        sigma2(ipk)=s2
        gamma(ipk)=gmpeak
        rstar(ipk)=rspeak
        rsthm1=h0*rstar(ipk)
        rstgs=rstar(ipk)/rgauss

        write(*,'('' (*) Spectral parameter information: '')')
        write(*,'(''     ------------------------------'')')
        write(*,*)
        write(*,'(''      sigma0(R_f):           '',$)')
        write(*,'(3x,f15.9)') sigma0(ipk)
        write(*,'(''      sigma1(R_f):           '',$)')
        write(*,'(3x,f15.9,a11)') sigma1(ipk),'   Mpc^{-1}'
        write(*,'(''      sigma2(R_f):           '',$)')       
        write(*,'(3x,f15.9,a11)') sigma2(ipk),'   Mpc^{-2}'
        write(*,*)
        write(*,'(''      Gamma:                 '',$)')
        write(*,'(3x,f15.9)') gamma(ipk)
        write(*,'(''      R_*:                   '',$)')
        write(*,'(3x,f15.9,a11)') rstar(ipk),'        Mpc'
        write(*,'(''                             '',$)')
        write(*,'(3x,f15.9,a11)') rsthm1,' h^{-1} Mpc'
        write(*,'(''      R_*/R_f:               '',3x,f15.9)') rstgs
        write(*,*)
        write(1,'('' (*) Spectral parameter information: '')')
        write(1,'(''     ------------------------------'')')
        write(1,*)
        write(1,'(''      sigma0(R_f):           '',$)')
        write(1,'(3x,f15.9)') sigma0(ipk)
        write(1,'(''      sigma1(R_f):           '',$)')
        write(1,'(3x,f15.9,a11)') sigma1(ipk),'   Mpc^{-1}'
        write(1,'(''      sigma2(R_f):           '',$)')       
        write(1,'(3x,f15.9,a11)') sigma2(ipk),'   Mpc^{-2}'
        write(1,*)
        write(1,'(''      Gamma:                 '',$)')
        write(1,'(3x,f15.9)') gamma(ipk)
        write(1,'(''      R_*:                   '',$)')
        write(1,'(3x,f15.9,a11)') rstar(ipk),'        Mpc'
        write(1,'(''                             '',$)')
        write(1,'(3x,f15.9,a11)') rsthm1,' h^{-1} Mpc'
        write(1,'(''      R_*/R_f:               '',3x,f15.9)') rstgs
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine csppos(ipk,iqu)


      implicit double precision (a-h,o-z)

      parameter         (npeakm=3)

      integer           ipk,iqu
      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0
      character         answer*1

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /box/boxlen


        if (psign(ipk).lt.0.0) then
          write(*,'('' (*) Position of void centre : '')')
          write(*,'(''     -----------------------'')')
          write(1,'('' (*) Position of void centre: '')')
          write(1,'(''     -----------------------'')')
        else
          write(*,'('' (*) Position of peak centre: '')')
          write(*,'(''     -----------------------'')')
          write(1,'('' (*) Position of peak centre: '')')
          write(1,'(''     -----------------------'')')
        endif
        write(*,*)
        write(1,*)
        write(*,'(''      At centre of simulation box [y/N]: '',$)')
        read(*,'(a)') answer
        write(1,'(''      At centre of simulation box [y/N]: '',$)')
        write(1,'(a)') answer
        if ((answer.eq.'n').or.(answer.eq.'N')) then
          write(*,*)
          write(1,*)
          if (iqu.eq.1) then
            write(*,'(''      Position in boxlength units'',$)')
            write(*,'(''   [0.0,1.0]: '')')
            write(1,'(''      Position in boxlength units'',$)')
            write(1,'(''   [0.0,1.0]: '')')
          else
            write(*,'(''      Position in h^{-1} Mpc'',$)')
            write(*,'(''   [0.0,'',f6.2,'']: '')') h0*boxlen
            write(1,'(''      Position in h^{-1} Mpc'',$)')
            write(1,'(''   [0.0,'',f6.2,'']: '')') h0*boxlen
          endif
          write(*,'(''       [X]: '',$)')
          read(*,*) centre(1,ipk)
          write(*,'(''       [Y]: '',$)')
          read(*,*) centre(2,ipk)
          write(*,'(''       [Z]: '',$)')
          read(*,*) centre(3,ipk)
          write(1,'(''      [X]: '',$)') 
          write(1,'(f10.3)') centre(1,ipk)
          write(1,'(''      [Y]: '',$)')
          write(1,'(f10.3)') centre(2,ipk)
          write(1,'(''      [Z]: '',$)')
          write(1,'(f10.3)') centre(3,ipk)
          if (iqu.eq.2) then
            centre(1,ipk)=centre(1,ipk)/(h0*boxlen)
            centre(2,ipk)=centre(2,ipk)/(h0*boxlen)
            centre(3,ipk)=centre(3,ipk)/(h0*boxlen)
          endif
        else
          centre(1,ipk)=0.5
          centre(2,ipk)=0.5
          centre(3,ipk)=0.5
        endif
        write(*,*)
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine csphgt(ipk)


      implicit double precision (a-h,o-z)

      parameter    (npeakm=3)

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)


        if (psign(ipk).ge.0.0) then
          write(*,'('' (*) Peak height:  '')')
          write(*,'(''     -----------'')')
          write(1,'('' (*) Peak height:  '')')
          write(1,'(''     -----------'')')
        endif
        if (psign(ipk).lt.0.0) then
          write(*,'('' (*) Void depth: '')')
          write(*,'(''     ----------'')')
          write(1,'('' (*) Void depth: '')')
          write(1,'(''     ----------'')')
        endif
        write(*,*)
        write(*,'(''      How to specify'',$)')
        write(*,'('' |delta|=|delta rho/rho|: '')')
        write(*,'(''      [1] |delta| '')')
        write(*,'(''      [2] nu=|delta|/sigma0 '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(*,*)
        write(1,*)
        write(1,'(''      How to specify'',$)')
        write(1,'('' |delta|=|delta rho/rho|: '')')
        write(1,'(''      [1] |delta| '')')
        write(1,'(''      [2] nu=|delta|/sigma0 '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic
        write(1,*)

        if (ichoic.eq.1) then
          write(*,'(''      |delta|:        '',$)')
          read(*,*) pnu(ipk)
          write(1,'(''      |delta|:        '',$)')
          write(1,'(f10.3)') pnu(ipk)
          pnu(ipk)=pnu(ipk)/sigma0(ipk)
        endif
        if (ichoic.eq.2) then
          write(*,'(''      nu:             '',$)')
          read(*,*) pnu(ipk)
          write(1,'(''      nu:             '',$)')
          write(1,'(f10.3)') pnu(ipk)
        endif
        if (psign(ipk).lt.0.0) then
          pnu(ipk)=-pnu(ipk)
        endif

        write(*,*)
        write(*,'(''      delta:      '',f12.6)') pnu(ipk)*sigma0(ipk)
        write(*,'(''      nu:         '',f12.6)') pnu(ipk)
        write(*,*)
        write(1,*)
        write(1,'(''      delta:      '',f17.6)') pnu(ipk)*sigma0(ipk)
        write(1,'(''      nu:         '',f17.6)') pnu(ipk)
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine csgrad(ipk)


      implicit double precision (a-h,o-z)

      parameter    (npeakm=3)

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)


        write(*,'('' (*) Gradient at object location:  '')')
        write(*,'(''     ----------------------------'')')
        write(1,'('' (*) Gradient at object location:  '')')
        write(1,'(''     ----------------------------'')')
        write(*,*)
        write(1,*)

c        write(*,'(''      [1]   zero (peak) value: 0.0   '')')
c        write(*,'(''      [2]   non-zero value           '')')
c        write(*,'(''      Choice:  '',$)')
c        read(*,*) ichoice
c        write(1,'(''      [1]   zero (peak) value: 0.0   '')')
c        write(1,'(''      [2]   non-zero value           '')')
c        write(1,'(''      Choice:  '',$)')
c        read(1,*) ichoice
        ichoice=1
        if (ichoice.eq.1) then
          do 10 k=1,3
            pgrd(k,ipk)=0.0
10        continue
        endif
c        if (ichoice.eq.2) then
c          do 20 k=1,3
c            pgrd(k,ipk)=1.0
c 20       continue
c        endif    

        write(*,*)
        write(*,'(''      grad f:     '',3f12.6)') (pgrd(m,ipk),m=1,3)
        write(*,*)
        write(1,*)
        write(1,'(''      grad f:     '',3f12.6)') (pgrd(m,ipk),m=1,3)
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspcur(ipk,lpkqst)


      implicit double precision(a-h,o-z)

      external  curvav

      parameter (npeakm=3,nxmax=5000)

      double precision  pxpnuc 
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,
     2                  gamma,rstar
      real              xr
      logical           lpkqst(npeakm,7)

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /curvat/xlow,xdel,pxpnuc(nxmax),npx


        pnupk=dabs(pnu(ipk))
        gampk=gamma(ipk)    
        xmean=curvav(gampk,pnupk)
        call xfunc(gampk,pnupk)

        if (lpkqst(ipk,3)) then
          write(*,'('' (*) Curvature x=Delta_F/sigma_2: '')')
          write(*,'(''     ---------------------------'')')
          write(*,*)
          write(*,'(''      <x>  (mean curvature): '',$)')
          write(*,'(3x,f15.9)') xmean
          write(*,*)
          write(*,'(''      How to specify x : '')')
          write(*,'(''      [1]  x '')')
          write(*,'(''      [2]  fx=x/<x> '')')
          write(*,'(''      [3]  percentile '')')
          write(*,'(''      [4]  random '')')
          write(*,'(''      choice:  '',$)')
          read(*,*) ichoic
          write(*,*)
          write(1,'('' (*) Curvature x=Delta_F/sigma_2: '')')
          write(1,'(''     ---------------------------'')')
          write(1,*)
          write(1,'(''      <x>  (mean curvature): '',$)')
          write(1,'(3x,f15.9)') xmean
          write(1,*)
          write(1,'(''      How to specify x : '')')
          write(1,'(''      [1]  x '')')
          write(1,'(''      [2]  fx=x/<x> '')')
          write(1,'(''      [3]  percentile '')')
          write(1,'(''      [4]  random '')')
          write(1,'(''      choice:  '',$)')
          write(1,'(i3)') ichoic
          write(1,*)
        else
          ichoic=4
        endif        
  
        if (ichoic.eq.1) then
          write(*,'(''      x:               '',$)')
          read(*,*) px(ipk)
          write(1,'(''      x:                 '',$)')
          write(1,'(f11.6)') px(ipk)
        endif
        if (ichoic.eq.2) then
          write(*,'(''      fx:              '',$)')
          read(*,*) px(ipk)
          write(1,'(''      fx:                '',$)')
          write(1,'(f11.6)') px(ipk)
          px(ipk)=px(ipk)*xmean
        endif
        if (ichoic.le.2) then
          jx=int((px(ipk)-xlow)/xdel)+1
          x0=float(jx-1)*xdel
          xperc=pxpnuc(jx)+(pxpnuc(jx+1)-pxpnuc(jx))*(px(ipk)-x0)/xdel
        endif
        if (ichoic.ge.3) then
          if (ichoic.eq.3) then
            write(*,'(''      percentile:    '',$)')
            read(*,*) xperc
            write(1,'(''      percentile:      '',$)')
            write(1,'(f11.6)') xperc
            xperc=xperc/100.0
          endif         
          if (ichoic.eq.4) then
            call randa(xr)
            xperc=xr
          endif        
          if (xperc.ge.pxpnuc(npx)) then
            px(ipk)=pxpnuc(npx)
          else
            call locate(pxpnuc,npx,xperc,jx)
            px(ipk)=xdel*(xperc-pxpnuc(jx))/(pxpnuc(jx+1)-pxpnuc(jx))
            px(ipk)=float(jx-1)*xdel+px(ipk)
          endif
        endif
        xperc=100.0*xperc
        if (psign(ipk).lt.0.0) then
          px(ipk)=-px(ipk)
        endif

        if (lpkqst(ipk,3)) then
          write(*,*)
          write(*,'(''      x:            '',f11.6)') px(ipk)
          write(*,'(''      percentile:   '',f11.6)') xperc
          write(*,*)
          write(1,*)
          write(1,'(''      x:            '',f16.6)') px(ipk)
          write(1,'(''      percentile:   '',f16.6)') xperc
          write(1,*)
        endif

1000   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspshp(ipk,lpkqst)


      implicit double precision(a-h,o-z)

      parameter (npeakm=3,nemax=257,npmax=257,nepmax=nemax*npmax)

      integer           nep,ihshp
      double precision  pshpcm,hshp 
      logical           lepsen,lpkqst(npeakm,7)

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /shape/pshpcm(nepmax),hshp(nepmax),ihshp(2,nepmax),nep
      common /shapep/emean,pmean,emax,pmax,eppred,pppred


        if (lpkqst(ipk,4)) then
          if (psign(ipk).lt.0.0) then
            write(*,'('' (*) Shape near void centre: '')')
            write(*,'(''     ----------------------'')')
            write(1,'('' (*) Shape near void centre: '')')
            write(1,'(''     ----------------------'')')
          else
            write(*,'('' (*) Shape near peak centre: '')')
            write(*,'(''     ----------------------'')')
            write(1,'('' (*) Shape near peak centre: '')')
            write(1,'(''     ----------------------'')')
          endif
        endif

        xpeak=dabs(px(ipk))
        call shfunc(xpeak)

        if (lpkqst(ipk,4)) then
          call epax(emean,pmean,a12mn,a13mn)
          call epax(emax,pmax,a12mx,a13mx)
          call epax(eppred,pppred,a12pr,a13pr)

          a21mn=1.0/a12mn
          a31mn=1.0/a13mn
          a21mx=1.0/a12mx
          a31mx=1.0/a13mx
          a21pr=1.0/a12pr
          a31pr=1.0/a13pr

          write(*,*)
          write(*,'(''      Average values e,p: '')')
          write(*,'(''      (1) e:       '',f12.6)') emean
          write(*,'(''      (2) p:       '',f12.6)') pmean
          write(*,'(''      Corresponding values axis ratios: '')')
          write(*,'(''      (1) a2/a1:   '',f12.6)') a21mn
          write(*,'(''      (2) a3/a1:   '',f12.6)') a31mn
          write(*,*)
          write(*,'(''      Most probable values e,p'',$)')
          write(*,'('' (+ approximate values from BBKS): '')')
          write(*,'(''      (1) e:       '',2f12.6)') emax,eppred
          write(*,'(''      (2) p:       '',2f12.6)') pmax,pppred
          write(*,'(''      Corresponding values axis ratios: '')')
          write(*,'(''      (1) a2/a1:   '',2f12.6)') a21mx,a21pr
          write(*,'(''      (2) a3/a1:   '',2f12.6)') a31mx,a31pr
          write(1,*)
          write(1,'(''      Average values e,p: '')')
          write(1,'(''      (1) e:       '',f12.6)') emean
          write(1,'(''      (2) p:       '',f12.6)') pmean
          write(1,'(''      Corresponding values axis ratios: '')')
          write(1,'(''      (1) a2/a1:   '',f12.6)') a21mn
          write(1,'(''      (2) a3/a1:   '',f12.6)') a31mn
          write(1,*)
          write(1,'(''      Most probable values e,p '',$)')
          write(1,'('' (+ approximate values from BBKS): '')')
          write(1,'(''      (1) e:       '',2f12.6)') emax,eppred
          write(1,'(''      (2) p:       '',2f12.6)') pmax,pppred
          write(1,'(''      Corresponding values axis ratios: '')')
          write(1,'(''      (1) a2/a1:   '',2f12.6)') a21mx,a21pr
          write(1,'(''      (2) a3/a1:   '',2f12.6)') a31mx,a31pr   

          write(*,*)    
          write(*,'(''      How to specify shape:  '')')
          write(*,'(''      [1] axis ratios a2/a1, a3/a1   '',$)')
          write(*,'('' (a1>a2>a3) '')')
          write(*,'(''      [2] e,p    '')')
          write(*,'(''      [3] random '')')
          write(*,'(''      choice:  '',$)')
          read(*,*) ichoic
          write(*,*)
          write(1,*)    
          write(1,'(''      How to specify shape:  '')')
          write(1,'(''      [1] axis ratios a2/a1, a3/a1   '',$)')
          write(1,'('' (a1>a2>a3) '')')
          write(1,'(''      [2] e,p    '')')
          write(1,'(''      [3] random '')')
          write(1,'(''      choice:  '',$)')
          write(1,'(i3)') ichoic
          write(1,*)
        else
          ichoic=3
        endif

        if (ichoic.eq.1) then
          write(*,'(''      a1>a2>a3:  '')')
          write(*,'(''       a2/a1                    : '',$)')
          read(*,*) a21pk
          write(*,'(''       a3/a1                    : '',$)')
          read(*,*) a31pk
          write(1,'(''      a1>a2>a3:  '')')
          write(1,'(''       a2/a1                    : '',$)')
          write(1,'(f10.3)') a21pk
          write(1,'(''       a3/a1                    : '',$)')
          write(1,'(f10.3)') a31pk
        endif
        if (ichoic.eq.2) then
1120      write(*,'(''      0 < e < 0.25   ---> -e < p < e       '')')
          write(*,'(''      0.25 < e < 0.5 ---> -(1-3e) < p < e  '')')
          write(*,'(''       e                        : '',$)')
          read(*,*) epk
          write(*,'(''       p                        : '',$)')
          read(*,*) ppk
          write(1,'(''      0 < e < 0.25   ---> -e < p < e       '')')
          write(1,'(''      0.25 < e < 0.5 ---> -(1-3e) < p < e  '')')
          write(1,'(''       e                        : '',$)')
          write(1,'(f10.3)') epk
          write(1,'(''       p                        : '',$)')
          write(1,'(f10.3)') ppk
          lepsen=.false.
          if ((epk.ge.0.0).and.(epk.le.0.25)) then
            if ((ppk.ge.-epk).and.(ppk.le.epk)) then
              lepsen=.true.
            endif
          endif
          if ((epk.ge.0.25).and.(epk.le.0.5)) then
            if ((ppk.ge.(3.0*epk-1.0)).and.(ppk.le.epk)) then
              lepsen=.true.
            endif
          endif
          if (.not.lepsen) then
            write(*,'(''      e,p outside allowable range, '',$)')
            write(*,'(''make new choice: '')')
            write(1,'(''      e,p outside allowable range, '',$)')
            write(1,'(''make new choice: '')')
            goto 1120
          endif       
        endif
        if (ichoic.eq.3) then
          call shaprn(epk,ppk,px(ipk))
        endif

        if (ichoic.eq.1) then
          a12pk=1.0/a21pk
          a13pk=1.0/a31pk
          call axep(a12pk,a13pk,epk,ppk)
        endif
        if ((ichoic.eq.2).or.(ichoic.eq.3)) then
          call epax(epk,ppk,a12pk,a13pk)
          a21pk=1.0/a12pk
          a31pk=1.0/a13pk
        endif
        a12(ipk)=1.0/a21pk
        a13(ipk)=1.0/a31pk

        if (lpkqst(ipk,4)) then
          pshape=pepx(epk,ppk,xpeak)
          call locate(hshp,nep,pshape,jx)
          shperc=(pshpcm(jx+1)-pshpcm(jx))*(pshape-hshp(jx))
          shperc=pshpcm(jx)+shperc/(hshp(jx+1)-hshp(jx))
          shperc=100.0*shperc

          write(*,*)
          write(*,'(''      a21, a31 :    '',$)')
          write(*,'(2f11.6)') a21pk,a31pk
          write(*,'(''      e, p:         '',$)')
          write(*,'(2f11.6)') epk,ppk
          write(*,'(''      percentile:   '',f11.6)') shperc
          write(*,*)
          write(*,*)
          write(1,*)
          write(1,'(''      a21, a31 :    '',$)')
          write(1,'(2f11.6)') a21pk,a31pk
          write(1,'(''      e, p:         '',$)')
          write(1,'(2f11.6)') epk,ppk
          write(1,'(''      percentile:   '',f11.6)') shperc
          write(1,*)
          write(1,*)
        endif

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine csporn(ipk,lpkqst)


      implicit double precision(a-h,o-z)

      double precision  pi,raddeg
      parameter (pi=3.141592653d0,raddeg=180.0/pi)
      parameter (npeakm=3)

      real      xr
      logical   lpkqst(npeakm,7)

      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    


        if (lpkqst(ipk,5)) then
          if (psign(ipk).lt.0.0) then
            write(*,'('' (*) Orientation of void at centre: '')')
            write(*,'(''     ----------------------------- '')')
            write(1,'('' (*) Orientation of void at centre: '')')
            write(1,'(''     ----------------------------- '')')
          else
            write(*,'('' (*) Orientation of peak at centre: '')')
            write(*,'(''     ----------------------------- '')')
            write(1,'('' (*) Orientation of peak at centre: '')')
            write(1,'(''     ----------------------------- '')')
          endif

          write(*,*)    
          write(*,'(''      How to specify orientation:  '')')
          write(*,'(''      [1] Euler angles phi, theta, psi'',$)')
          write(*,'('' (degrees):   '')')
          write(*,'(''      [2] random '')')
          write(*,'(''      choice:  '',$)')
          read(*,*) ichoic
          write(*,*)
          write(1,*)    
          write(1,'(''      How to specify orientation:  '')')
          write(1,'(''      [1] Euler angles phi, theta, psi'',$)')
          write(1,'('' (degrees):   '')')
          write(1,'(''      [2] random '')')
          write(1,'(''      choice:  '',$)')
          write(1,'(i3)') ichoic
          write(1,*)
        else
          ichoic=2
        endif

        if (ichoic.eq.1) then 
          write(*,'(''      "phi"    [0.0,180.0]: '',$)')
          read(*,*) phi(ipk)
          write(*,'(''      "theta"   [0.0,90.0]: '',$)')
          read(*,*) theta(ipk)
          write(*,'(''      "psi"    [0.0,180.0]: '',$)')
          read(*,*) psi(ipk)
          write(*,*)
        endif
        if (ichoic.eq.2) then
          call randa(xr)
          phi(ipk)=raddeg*pi*xr
          call randa(xr)
          theta(ipk)=raddeg*acos(1.0-xr)
          call randa(xr)
          psi(ipk)=raddeg*pi*xr
        endif

        if (lpkqst(ipk,5)) then
          if (ichoic.eq.2) then
            write(*,'(''      "phi"    : '',$)')
            write(*,'(f10.3)') phi(ipk)
            write(*,'(''      "theta"  : '',$)')
            write(*,'(f10.3)') theta(ipk)
            write(*,'(''      "psi"    : '',$)')
            write(*,'(f10.3)') psi(ipk)
            write(*,*)
          endif
          write(1,'(''      "phi"    : '',$)')
          write(1,'(f10.3)') phi(ipk)
          write(1,'(''      "theta"  : '',$)')
          write(1,'(f10.3)') theta(ipk)
          write(1,'(''      "psi"    : '',$)')
          write(1,'(f10.3)') psi(ipk)
          write(1,*)
        endif

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspvel(ipk)


      implicit double precision(a-h,o-z)

      external          pgssig

      integer           npeakm
      double precision  twopi
      parameter         (npeakm=3,twopi=6.283185307179586d0)

      double precision   sigmvp,sigmvf
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar,rgauss
      real               boxlen

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common /box/boxlen


        write(*,'('' (*) Velocity of peak: '')')
        write(*,'(''     ---------------- '')')
        write(1,'('' (*) Velocity of peak: '')')
        write(1,'(''     ---------------- '')')

        rgauss=pscale(ipk)*boxlen
        pindex=2.0*(-1.0)
        sigmm1=pgssig(rgauss)
        gammav=sigma0(ipk)*sigma0(ipk)/(sigmm1*sigma1(ipk))
        sigmvp=hubbl0*sigmm1*sqrt(1.0-gammav*gammav)          
        sigmvf=hubbl0*sigmm1
        sigmav(ipk)=sigmvp

        write(*,*)
        write(*,'(''      velocity dispersion peaks,'',$)')
        write(*,'('' sigma_v_pk: '',$)')
        write(*,'(1x,f12.6,a6)') sigmvp,'  km/s'
        write(*,'(''      velocity dispersion field,'',$)')
        write(*,'('' sigma_v_f:  '',$)')
        write(*,'(1x,f12.6,a6)') sigmvf,'  km/s'
        write(*,*)
        write(1,*)
        write(1,'(''      velocity dispersion peaks,'',$)')
        write(1,'('' sigma_v_pk: '',$)')
        write(1,'(1x,f12.6,a6)') sigmvp,'  km/s'
        write(1,'(''      velocity dispersion field,'',$)')
        write(1,'('' sigma_v_f:  '',$)')
        write(1,'(1x,f12.6,a6)') sigmvf,'  km/s'
        write(1,*)

        write(*,'(''      How to specify velocity: '')')
        write(*,'(''      [1] Cartesian components vx,vy,vz '')')
        write(*,'(''      [2] Spherical components v,theta_v,'',$)')
        write(*,'(''phi_v '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,'(''      How to specify velocity: '')')
        write(1,'(''      [1] Cartesian components vx,vy,vz '')')
        write(1,'(''      [2] Spherical components v,theta_v,'',$)')
        write(1,'(''phi_v '')')
        write(1,'(''      choice: '',$)')
        write(1,*) ichoic

        if (ichoic.eq.1) call cspvlc(sigmvp,ipk)
        if (ichoic.eq.2) call cspvls(sigmvp,ipk)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspvlc(sigmvp,ipk)


      implicit double precision(a-h,o-z)

      integer           npeakm,nvmax
      double precision  twopi
      parameter         (npeakm=3,nvmax=5000,twopi=6.283185307179586d0)

      double precision   sigmvp
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar
      double precision   vpeak,vtheta,vphi
      double precision   xr1,xr2
      real               xr

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /pecvel/vlow,vdel,velpkc(nvmax),npv
      common /vpeak/vpeak,vtheta,vphi


        call vcfunc(sigmvp)

        write(*,*)
        write(*,'(''      How to specify velocity: '')')
        write(*,'(''      [1]    v in km/s '')')
        write(*,'(''      [2]    fv=v/sigma_v_pk '')')
        write(*,'(''      [3]    random '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,*)
        write(1,'(''      How to specify velocity: '')')
        write(1,'(''      [1]    v in km/s '')')
        write(1,'(''      [2]    fv=v/sigma_v_pk '')')
        write(1,'(''      [3]    random '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic
                             
        if (ichoic.eq.1) then
          write(*,*)
          write(*,'(''      v_x:           '',$)')
          read(*,*) vx(ipk)
          write(*,'(''      v_y:           '',$)')
          read(*,*) vy(ipk)
          write(*,'(''      v_z:           '',$)')
          read(*,*) vz(ipk)
          write(1,*)
          write(1,'(''      v_x:             '',$)')
          write(1,'(f11.6)') vx(ipk)
          write(1,'(''      v_y:             '',$)')
          write(1,'(f11.6)') vy(ipk)
          write(1,'(''      v_z:             '',$)')
          write(1,'(f11.6)') vz(ipk)
        endif
        if (ichoic.eq.2) then
          write(*,*)
          write(*,'(''      fv_x:              '',$)')
          read(*,*) vx(ipk)
          write(*,'(''      fv_y:              '',$)')
          read(*,*) vy(ipk)
          write(*,'(''      fv_z:              '',$)')
          read(*,*) vz(ipk)
          write(1,*)
          write(1,'(''      fv_x:                '',$)')
          write(1,'(f11.6)') vx(ipk)
          write(1,'(''      fv_y:                '',$)')
          write(1,'(f11.6)') vy(ipk)
          write(1,'(''      fv_z:                '',$)')
          write(1,'(f11.6)') vz(ipk)
          vx(ipk)=vx(ipk)*sigmvp
          vy(ipk)=vy(ipk)*sigmvp
          vz(ipk)=vz(ipk)*sigmvp
        endif
        if (ichoic.eq.3) then
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          vx(ipk)=sigmvp*dsqrt(-2.0*dlog(xr1))*dcos(xr2)/sqrt(3.0)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          vy(ipk)=sigmvp*dsqrt(-2.0*dlog(xr1))*dcos(xr2)/sqrt(3.0)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          vz(ipk)=sigmvp*dsqrt(-2.0*dlog(xr1))*dcos(xr2)/sqrt(3.0)
        endif

        vpeak=dsqrt(vx(ipk)*vx(ipk)+vy(ipk)*vy(ipk)+vz(ipk)*vz(ipk))
        if (vpeak.gt.0.0) then
          vtheta=vz(ipk)/vpeak
          vtheta=dacosd(vtheta)
          if ((vx(ipk).lt.1.0d-10).and.(vy(ipk).lt.1.0d-10)) then
            vphi=0.0
          else
            vphi=vx(ipk)/dsqrt(vx(ipk)*vx(ipk)+vy(ipk)*vy(ipk))
            vphi=dacosd(vphi)
          endif
          if (vy(ipk).lt.0.0) vphi=360.0-vphi
        else
          vtheta=0.0
          vphi=0.0
        endif

        jv=int((vpeak-vlow)/vdel)+1
        v0=float(jv-1)*vdel
        vperc=velpkc(jv)+(velpkc(jv+1)-velpkc(jv))*(vpeak-v0)/vdel
        vperc=100.0*vperc

        write(*,*)
        write(*,'(''      v:            '',f11.6,'' km/s'')') vpeak
        write(*,'(''      theta_v:      '',f11.6)') vtheta
        write(*,'(''      phi_v:        '',f11.6)') vphi
        write(*,*)
        write(*,'(''      v_x:          '',f11.6,'' km/s'')') vx(ipk)
        write(*,'(''      v_y:          '',f11.6,'' km/s'')') vy(ipk)
        write(*,'(''      v_z:          '',f11.6,'' km/s'')') vz(ipk)
        write(*,*)
        write(*,'(''      percentile:   '',f11.6)') vperc
        write(*,*)
        write(1,*)
        write(1,'(''      v:          '',f16.6,'' km/s'')') vpeak
        write(1,'(''      theta_v:         '',f11.6)') vtheta
        write(1,'(''      phi_v:           '',f11.6)') vphi
        write(1,*)
        write(1,'(''      v_x:             '',f11.6,'' km/s'')') vx(ipk)
        write(1,'(''      v_y:             '',f11.6,'' km/s'')') vy(ipk)
        write(1,'(''      v_z:             '',f11.6,'' km/s'')') vz(ipk)
        write(1,*)
        write(1,'(''      percentile: '',f16.6)') vperc
        write(1,*)

1000   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspvls(sigmvp,ipk)


      implicit double precision(a-h,o-z)

      integer           npeakm,nvmax,dim
      double precision  twopi,raddeg
      parameter         (npeakm=3,nvmax=5000,dim=3)
      parameter         (twopi=6.283185307179586d0,raddeg=360.0/twopi)

      double precision   sigmvp
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar
      double precision   vpeak,vtheta,vphi
      double precision   reuliv(3,3),vpk(dim)
      real               xr

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /pecvel/vlow,vdel,velpkc(nvmax),npv
      common /vpeak/vpeak,vtheta,vphi


        call vcfunc(sigmvp)

        write(*,*)
        write(*,'(''      How to specify velocity amplitude: '')')
        write(*,'(''      [1]    v in km/s '')')
        write(*,'(''      [2]    fv=v/sigma_v_pk '')')
        write(*,'(''      [3]    percentile '')')
        write(*,'(''      [4]    random '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,*)
        write(1,'(''      How to specify velocity amplitude: '')')
        write(1,'(''      [1]    v in km/s '')')
        write(1,'(''      [2]    fv=v/sigma_v_pk '')')
        write(1,'(''      [3]    percentile '')')
        write(1,'(''      [4]    random '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic
                             
        if (ichoic.eq.1) then
          write(*,*)
          write(*,'(''      v_pec:           '',$)')
          read(*,*) vpeak
          write(1,*)
          write(1,'(''      v_pec:             '',$)')
          write(1,'(f11.6)') vpeak
        endif
        if (ichoic.eq.2) then
          write(*,*)
          write(*,'(''      fv:              '',$)')
          read(*,*) vpeak
          write(1,*)
          write(1,'(''      fv:                '',$)')
          write(1,'(f11.6)') vpeak
          vpeak=vpeak*sigmvp
        endif
        if (ichoic.le.2) then
          jv=int((vpeak-vlow)/vdel)+1
          v0=float(jv-1)*vdel
          vperc=velpkc(jv)+(velpkc(jv+1)-velpkc(jv))*(vpeak-v0)/vdel
        endif
        if (ichoic.ge.3) then
          if (ichoic.eq.3) then
            write(*,*)
            write(*,'(''      percentile:    '',$)')
            read(*,*) vperc
            write(1,*)
            write(1,'(''      percentile:      '',$)')
            write(1,'(f11.6)') vperc
            vperc=vperc/100.0
          endif         
          if (ichoic.eq.4) then
            call randa(xr)
            vperc=xr
          endif        
          if (vperc.ge.velpkc(npv)) then
            vpeak=velpkc(npv)
          else
            call locate(velpkc,npv,vperc,jv)
            vpeak=vdel*(vperc-velpkc(jv))/(velpkc(jv+1)-velpkc(jv))
            vpeak=float(jv-1)*vdel+vpeak
          endif
        endif

        write(*,*)
        write(*,'(''      How to specify velocity direction: '')')
        write(*,'(''      [1]    with respect to box axes'')')
        write(*,'(''      [2]    with respect to peak axes'')')
        write(*,'(''      [3]    random '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,*)
        write(1,'(''      How to specify velocity direction: '')')
        write(1,'(''      [1]    with respect to box axes'')')
        write(1,'(''      [2]    with respect to peak axes'')')
        write(1,'(''      [3]    random '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic
                             
        write(*,*)
        write(*,*)
        write(1,*)
        write(1,*)
        if (ichoic.le.2) then 
          write(*,'(''      Velocity direction angles (degrees):  '')')
          write(*,'(''      "theta_v"  : '',$)')
          read(*,*) vtheta
          write(*,'(''      "phi_v"    : '',$)')
          read(*,*) vphi
          write(*,*)
          write(1,'(''      Velocity direction angles (degrees):  '')')
          write(1,'(''      "theta_v"  : '',$)')
          write(1,'(f10.3)') vtheta
          write(1,'(''      "phi_v"    : '',$)')
          write(1,'(f10.3)') vphi
          write(1,*)
        endif
        if (ichoic.eq.3) then
          call randa(xr)
          vtheta=raddeg*acos(1.0-2.0*xr)
          call randa(xr)
          vphi=raddeg*twopi*xr
          write(*,'(''      "theta_v"  : '',$)')
          write(*,'(f10.3)') vtheta
          write(*,'(''      "phi_v"    : '',$)')
          write(*,'(f10.3)') vphi
          write(1,'(''      "theta_v"  : '',$)')
          write(1,'(f10.3)') vtheta
          write(1,'(''      "phi_v"    : '',$)')
          write(1,'(f10.3)') vphi
        endif

        vx(ipk)=vpeak*sind(vtheta)*cosd(vphi)
        vy(ipk)=vpeak*sind(vtheta)*sind(vphi)
        vz(ipk)=vpeak*cosd(vtheta)

        if (ichoic.eq.2) then
          vpk(1)=vx(ipk)
          vpk(2)=vy(ipk)
          vpk(3)=vz(ipk)
          vx(ipk)=0.0
          vy(ipk)=0.0
          vz(ipk)=0.0
          call eulinv(reuliv,phi(ipk),theta(ipk),psi(ipk))
          do 10 m=1,dim
            vx(ipk)=vx(ipk)+vpk(m)*reuliv(1,m)
            vy(ipk)=vy(ipk)+vpk(m)*reuliv(2,m)
            vz(ipk)=vz(ipk)+vpk(m)*reuliv(3,m)
10        continue

          vtheta=vz(ipk)/vpeak
          vphi=vx(ipk)/dsqrt(vx(ipk)*vx(ipk)+vy(ipk)*vy(ipk))
          vtheta=dacosd(vtheta)
          vphi=dacosd(vphi)
          if (vy(ipk).lt.0.0) vphi=360.0-vphi
        endif

        vperc=100.0*vperc
        write(*,*)
        write(*,'(''      v:            '',f11.6,'' km/s'')') vpeak
        write(*,'(''      theta_v:      '',f11.6)') vtheta
        write(*,'(''      phi_v:        '',f11.6)') vphi
        write(*,*)
        write(*,'(''      v_x:          '',f11.6,'' km/s'')') vx(ipk)
        write(*,'(''      v_y:          '',f11.6,'' km/s'')') vy(ipk)
        write(*,'(''      v_z:          '',f11.6,'' km/s'')') vz(ipk)
        write(*,*)
        write(*,'(''      percentile:   '',f11.6)') vperc
        write(*,*)
        write(1,*)
        write(1,'(''      v:            '',f16.6,'' km/s'')') vpeak
        write(1,'(''      theta_v:      '',f11.6)') vtheta
        write(1,'(''      phi_v:        '',f11.6)') vphi
        write(1,*)
        write(1,'(''      v_x:          '',f11.6,'' km/s'')') vx(ipk)
        write(1,'(''      v_y:          '',f11.6,'' km/s'')') vy(ipk)
        write(1,'(''      v_z:          '',f11.6,'' km/s'')') vz(ipk)
        write(1,*)
        write(1,'(''      percentile:   '',f16.6)') vperc
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspshr(ipk)


      implicit double precision(a-h,o-z)

      integer           npeakm
      double precision  twopi
      parameter         (npeakm=3,twopi=6.283185307179586d0)

      double precision   sigmvp,sigmvf
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar,rgauss
      double precision   sigms1f,sigms1p,sigms2f,sigms2p
      real               boxlen

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common /box/boxlen


        write(*,'('' (*) Shear of peak: '')')
        write(*,'(''     ------------- '')')
        write(1,'('' (*) Shear of peak: '')')
        write(1,'(''     ------------- '')')

        sigms1f=hf0*sigma0(ipk)/sqrt(5.0)
        sigms2f=hf0*sigma0(ipk)/sqrt(15.0)
        sigms1p=dsqrt(1.0-gamma(ipk)*gamma(ipk))*sigms1f
        sigms2p=dsqrt(1.0-gamma(ipk)*gamma(ipk))*sigms2f
        sigmas(1,ipk)=sigms1p
        sigmas(2,ipk)=sigms2p

        write(*,*)
        write(*,'(''      shear dispersion peaks: '')')
        write(*,'(''         diagonal components:     '',$)') 
        write(*,'(1x,f12.6,a10)') sigms1p,'  km/s/Mpc'
        write(*,'(''         off-diagonal components: '',$)')
        write(*,'(1x,f12.6,a10)') sigms2p,'  km/s/Mpc'
        write(*,*)
        write(*,'(''      shear dispersion field: '')')
        write(*,'(''         diagonal components:     '',$)') 
        write(*,'(1x,f12.6,a10)') sigms1f,'  km/s/Mpc'
        write(*,'(''         off-diagonal components: '',$)')
        write(*,'(1x,f12.6,a10)') sigms2f,'  km/s/Mpc'
        write(*,*)
        write(1,*)
        write(1,'(''      shear dispersion peaks: '')')
        write(1,'(''         diagonal components:     '',$)') 
        write(1,'(1x,f12.6,a10)') sigms1p,'  km/s/Mpc'
        write(1,'(''         off-diagonal components: '',$)')
        write(1,'(1x,f12.6,a10)') sigms2p,'  km/s/Mpc'
        write(1,*)
        write(1,'(''      shear dispersion field: '')')
        write(1,'(''         diagonal components:     '',$)') 
        write(1,'(1x,f12.6,a10)') sigms1f,'  km/s/Mpc'
        write(1,'(''         off-diagonal components: '',$)')
        write(1,'(1x,f12.6,a10)') sigms2f,'  km/s/Mpc'
        write(1,*)

        write(*,'(''      How to specify shear: '')')
        write(*,'(''      [1] 5 Cartesian components'',$)')
        write(*,'('' s11,s22,s12,s13,s23 '')')
        write(*,'(''      [2] Eigenvalues + orientation '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,'(''      How to specify shear: '')')
        write(1,'(''      [1] 5 Cartesian components'',$)')
        write(1,'('' s11,s22,s12,s13,s23 '')')
        write(1,'(''      [2] Eigenvalues + orientation '')')
        write(1,'(''      choice: '',$)')
        write(1,*) ichoic

        if (ichoic.eq.1) call cspshc(sigms1p,sigms2p,ipk)
        if (ichoic.eq.2) call cspshe(sigms1p,sigms2p,ipk)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspshc(sigms1p,sigms2p,ipk)


      implicit double precision(a-h,o-z)

      integer           npeakm,dim
      double precision  pi,twopi,raddeg
      parameter         (dim=3,npeakm=3)
      parameter         (pi=3.141592653d0,raddeg=180.0/pi,
     2                   twopi=2.0*pi)

      double precision   sigms1p,sigms2p
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   phish,thtsh,psish,shegv1,shegv2,psih,
     2                   shear(3,3),amat(3,3),dsh(3),vsh(3,3),
     2                   v1(3),v2(3),v3(3),vnl(3)
      double precision   epeak,ppeak,a12pk,a13pk,wrpk,vrpk
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar,rgauss
      real               boxlen,xr

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /vpeak/vpeak,vtheta,vphi
      common /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common /box/boxlen


        write(*,*) 
        write(*,'(''      How to specify shear components: '')')
        write(*,'(''      [1]    in km/s/Mpc '')')
        write(*,'(''      [2]    fs=s/(Hf sigma_0/sqrt(5.0)), '',$)')
        write(*,'('' resp. '')')
        write(*,'(''               =s/(Hf sigma_0/sqrt(15.0)) '')')
        write(*,'(''      [3]    random '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,*)
        write(1,'(''      How to specify shear eigenvalues: '')')
        write(1,'(''      [1]    s1,s2 in km/s/Mpc '')')
        write(1,'(''      [2]    fs=s/(Hf sigma_0/sqrt(5.0)), '',$)')
        write(1,'('' resp. '')')
        write(1,'(''               =s/(Hf sigma_0/sqrt(15.0)) '')')
        write(1,'(''      [3]    random '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic

        if (ichoic.eq.1) then
          write(*,*)
          write(*,'(''      s_11:           '',$)')
          read(*,*) sh1(ipk)
          write(*,'(''      s_22:           '',$)')
          read(*,*) sh2(ipk)
          write(*,'(''      s_12:           '',$)')
          read(*,*) sh3(ipk)
          write(*,'(''      s_13:           '',$)')
          read(*,*) sh4(ipk)
          write(*,'(''      s_23:           '',$)')
          read(*,*) sh5(ipk)
          write(1,*)
          write(1,'(''      s_11:         '',$)')
          write(1,'(f11.6)') sh1(ipk)
          write(1,'(''      s_22:         '',$)')
          write(1,'(f11.6)') sh2(ipk)
          write(1,'(''      s_12:         '',$)')
          write(1,'(f11.6)') sh3(ipk)
          write(1,'(''      s_13:         '',$)')
          write(1,'(f11.6)') sh4(ipk)
          write(1,'(''      s_23:         '',$)')
          write(1,'(f11.6)') sh5(ipk)
        endif
        if (ichoic.eq.2) then
          write(*,*)
          write(*,'(''      f_s_11:         '',$)')
          read(*,*) sh1(ipk)
          write(*,'(''      f_s_22:         '',$)')
          read(*,*) sh2(ipk)
          write(*,'(''      f_s_12:         '',$)')
          read(*,*) sh3(ipk)
          write(*,'(''      f_s_13:         '',$)')
          read(*,*) sh4(ipk)
          write(*,'(''      f_s_23:         '',$)')
          read(*,*) sh5(ipk)
          write(1,*)
          write(1,'(''      f_s_11:       '',$)')
          write(1,'(f11.6)') sh1(ipk)
          write(1,'(''      f_s_22:       '',$)')
          write(1,'(f11.6)') sh2(ipk)
          write(1,'(''      f_s_12:       '',$)')
          write(1,'(f11.6)') sh3(ipk)
          write(1,'(''      f_s_13:       '',$)')
          write(1,'(f11.6)') sh4(ipk)
          write(1,'(''      f_s_23:       '',$)')
          write(1,'(f11.6)') sh5(ipk)
          sh1(ipk)=sh1(ipk)*sigms1p
          sh2(ipk)=sh2(ipk)*sigms1p
          sh3(ipk)=sh3(ipk)*sigms2p
          sh4(ipk)=sh4(ipk)*sigms2p
          sh5(ipk)=sh5(ipk)*sigms2p
        endif
        if (ichoic.eq.3) then
          a12pk=a12(ipk)
          a13pk=a13(ipk)
          call axep(a12pk,a13pk,epeak,ppeak)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=xr
          wrpk=sigms1p*dsqrt(-2.0*dlog(xr1))*dcos(xr2)
          wrpk=wrpk+gamma(ipk)*px(ipk)*ppeak*hf0*sigma0(ipk)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=xr
          vrpk=sigms2p*dsqrt(-2.0*dlog(xr1))*dcos(xr2)
          vrpk=wrpk+gamma(ipk)*px(ipk)*epeak*hf0*sigma0(ipk)
          sh2(ipk)=2.0*wrpk/3.0
          sh1(ipk)=-vrpk-0.5*sh2(ipk)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          sh3(ipk)=sigms2p*dsqrt(-2.0*dlog(xr1))*dcos(xr2)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          sh4(ipk)=sigms2p*dsqrt(-2.0*dlog(xr1))*dcos(xr2)
          call randa(xr)
          xr1=xr
          call randa(xr)
          xr2=twopi*xr
          sh5(ipk)=sigms2p*dsqrt(-2.0*dlog(xr1))*dcos(xr2)
        endif

        shear(1,1)=sh1(ipk)
        shear(2,2)=sh2(ipk)
        shear(3,3)=-sh1(ipk)-sh2(ipk)
        shear(1,2)=sh3(ipk)
        shear(1,3)=sh4(ipk)
        shear(2,3)=sh5(ipk)
        shear(2,1)=shear(1,2)
        shear(3,1)=shear(1,3)
        shear(3,2)=shear(2,3)

        do 200 i=1,3
          do 210 j=1,3
            amat(i,j)=shear(i,j)
210       continue
200     continue

        na=3
        call jacobi(amat,na,dsh,vsh,nrot)
        call eigsrt(dsh,vsh,na)

        do 300 m=1,3
          v3(m)=vsh(m,3)
          v1(m)=vsh(m,1)
300     continue
        thtsh=dacosd(v3(3))
        if (dabs(v3(1)).lt.1.0d-10) then
          vnl(1)=1.0
        else
          if (dabs(v3(2)).lt.1.0d-10) then
            vnl(2)=1.0
          else
            if (v3(1).lt.0.0) then
              vnl(1)=v3(2)
              vnl(2)=-v3(1)
            else
              vnl(1)=-v3(2)
              vnl(2)=v3(1)
            endif
          endif
        endif
        vnln=dsqrt(vnl(1)*vnl(1)+vnl(2)*vnl(2))
        vnl(1)=vnl(1)/vnln
        vnl(2)=vnl(2)/vnln
        phish=dacosd(vnl(1))
        v2(1)=v1(3)*v3(2)-v1(2)*v3(3)
        v2(2)=v1(1)*v3(3)-v1(3)*v3(1)
        v2(3)=v1(2)*v3(1)-v1(1)*v3(2)
        psish=0.0
        psih=0.0
        do 310 m=1,3
          psish=psish+v1(m)*vnl(m)
          psih=psih+v2(m)*vnl(m)
310     continue
        if (psih.le.0.0) then
          psish=dacosd(psish)
        else
          psish=180.0-dacosd(psish)
        endif

        write(*,*)
        write(*,'(''      s_1:          '',f11.6,$)') dsh(1)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_2:          '',f11.6,$)') dsh(2)
        write(*,'('' km/s/Mpc'')') 
        write(*,'(''      phi_sh:       '',f11.6)') phish
        write(*,'(''      tht_sh:       '',f11.6)') thtsh
        write(*,'(''      psi_sh:       '',f11.6)') psish
        write(*,*)
        write(*,'(''      s_11:         '',f11.6,$)') sh1(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_22:         '',f11.6,$)') sh2(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_33:         '',f11.6,$)') shear(3,3)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_12:         '',f11.6,$)') sh3(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_13:         '',f11.6,$)') sh4(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_23:         '',f11.6,$)') sh5(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,*)
        write(1,*)
        write(1,'(''      s_1:          '',f11.6,$)') dsh(1)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_2:          '',f11.6,$)') dsh(2)
        write(1,'('' km/s/Mpc'')') 
        write(1,'(''      phi_sh:       '',f11.6)') phish
        write(1,'(''      tht_sh:       '',f11.6)') thtsh
        write(1,'(''      psi_sh:       '',f11.6)') psish
        write(1,*)
        write(1,'(''      s_11:         '',f11.6,$)') sh1(ipk)
        write(1,'('' km/s/Mpc'')') 
        write(1,'(''      s_22:         '',f11.6,$)') sh2(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_33:         '',f11.6,$)') shear(3,3)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_12:         '',f11.6,$)') sh3(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_13:         '',f11.6,$)') sh4(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_23:         '',f11.6,$)') sh5(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine cspshe(sigms1p,sigms2p,ipk)


      implicit double precision(a-h,o-z)

      integer           npeakm,nshmax,dim
      double precision  pi,raddeg
      parameter         (dim=3,npeakm=3,nshmax=5000)
      parameter         (pi=3.141592653d0,raddeg=180.0/pi)

      double precision   sigms1p,sigms2p
      double precision   h0,hubbl0,omega0,omegahh,hf0
      double precision   phish,thtsh,psish,shegv1,shegv2,psih
      double precision   vpeak,vtheta,vphi
      double precision   transf(3,3),tranf1(3,3),tranf2(3,3),lambsh(3),
     2                   shear(3,3)    
      double precision   v1(3),v2(3),v3(3),vnl(3)
      double precision   sigma0,sigma1,sigma2,sigmav,sigmas,
     2                   gamma,rstar,rgauss
      real               boxlen,xr

      common /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common /peak/pscale(npeakm),centre(3,npeakm),psign(npeakm),
     2             pnu(npeakm),pgrd(3,npeakm),
     2             px(npeakm),a12(npeakm),a13(npeakm),
     2             phi(npeakm),theta(npeakm),psi(npeakm),
     2             vx(npeakm),vy(npeakm),vz(npeakm),
     2             sh1(npeakm),sh2(npeakm),sh3(npeakm),
     2             sh4(npeakm),sh5(npeakm)    
      common /sigmap/sigma0(npeakm),sigma1(npeakm),sigma2(npeakm),
     2               sigmav(npeakm),sigmas(2,npeakm),
     2               gamma(npeakm),rstar(npeakm)
      common /vpeak/vpeak,vtheta,vphi
      common /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common /box/boxlen


        write(*,*) 
        write(*,'(''      How to specify shear eigenvalues: '')')
        write(*,'(''      [1]    s1,s2 in km/s/Mpc '')')
        write(*,'(''      [2]    fs=s/(Hf sigma_0/sqrt(5.0)) '')')
c        write(*,'(''      [3]    random '')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichoic
        write(1,*)
        write(1,'(''      How to specify shear eigenvalues: '')')
        write(1,'(''      [1]    s1,s2 in km/s/Mpc '')')
        write(1,'(''      [2]    fs=s/(Hf sigma_0/sqrt(5.0)) '')')
c        write(1,'(''      [3]    random '')')
        write(1,'(''      choice: '',$)')
        write(1,'(i3)') ichoic

        if (ichoic.eq.1) then
          write(*,*)
          write(*,'(''      s_1:           '',$)')
          read(*,*) lambsh(1)
          write(*,'(''      s_2:           '',$)')
          read(*,*) lambsh(2)
          write(1,*)
          write(1,'(''      s_1:          '',$)')
          write(1,'(f11.6)') lambsh(1)
          write(1,'(''      s_2:          '',$)')
          write(1,'(f11.6)') lambsh(2)
        endif
        if (ichoic.eq.2) then
          write(*,*)
          write(*,'(''      fs_1:              '',$)')
          read(*,*) lambsh(1)
          write(*,'(''      fs_2:              '',$)')
          read(*,*) lambsh(2)
          write(1,*)
          write(1,'(''      fs_1:         '',$)')
          write(1,'(f11.6)') lambsh(1)
          write(1,'(''      fs_2:         '',$)')
          write(1,'(f11.6)') lambsh(2)
          lambsh(1)=lambsh(1)*sigms1p
          lambsh(2)=lambsh(2)*sigms1p
        endif
c        if (ichoic.eq.3) then
c          call shrrn(lambsh(1),lambsh(2))
c        endif

        lambsh(3)=-lambsh(1)-lambsh(2)

        write(*,*)
        write(*,'(''      Shear orientation specified with'',$)')
        write(*,'('' respect to: '')')
        write(*,'(''      [1] box axes'')')
        write(*,'(''      [2] peak axes'')')
        write(*,'(''      [3] velocity-defined frame'')')
        write(*,'(''      choice: '',$)')
        read(*,*) ichc1
        write(*,*)
        write(*,'(''      How to specify the shear orientation:  '')')
        write(*,'(''      [1] Euler angles phi, theta, psi'',$)')
        write(*,'('' (degrees):   '')')
        write(*,'(''      [2] random '')')
        write(*,'(''      choice:  '',$)')
        read(*,*) ichc2
        write(*,*)
        write(1,*)
        write(1,'(''      Shear orientation specified with'',$)')
        write(1,'('' respect to: '')')
        write(1,'(''      [1] box axes'')')
        write(1,'(''      [2] peak axes'')')
        write(1,'(''      [3] velocity-defined frame'')')
        write(1,'(''      choice: '',$)')
        write(1,*) ichc1
        write(1,*)
        write(1,'(''      How to specify the shear orientation:  '')')
        write(1,'(''      [1] Euler angles phi, theta, psi'',$)')
        write(1,'('' (degrees):   '')')
        write(1,'(''      [2] random '')')
        write(1,'(''      choice:  '',$)')
        write(1,*) ichc2
        write(1,*)

        if (ichc2.eq.1) then 
          write(*,'(''      "phi_sh"    [0.0,180.0]: '',$)')
          read(*,*) phish
          write(*,'(''      "theta_sh"   [0.0,90.0]: '',$)')
          read(*,*) thtsh
          write(*,'(''      "psi_sh"    [0.0,180.0]: '',$)')
          read(*,*) psish
          write(*,*)
        endif
        if (ichc2.eq.2) then
          call randa(xr)
          phish=raddeg*pi*xr
          call randa(xr)
          thtsh=raddeg*acos(1.0-xr)
          call randa(xr)
          psish=raddeg*pi*xr
        endif
        write(1,'(''      "phi_sh"    [0.0,180.0]: '',$)')
        write(1,*) phish
        write(1,'(''      "theta_sh"   [0.0,90.0]: '',$)')
        write(1,*) thtsh
        write(1,'(''      "psi_sh"    [0.0,180.0]: '',$)')
        write(1,*) psish
        write(1,*)

        if (ichc1.eq.1) then
          call euler(transf,phish,thtsh,psish)
        endif
        if (ichc1.ge.2) then
          if (ichc1.eq.2) then
            call euler(tranf1,phi(ipk),theta(ipk),psi(ipk))
          endif
          if (ichc1.eq.3) then
            call vtrans(tranf1,vphi,vtheta)
          endif
          call euler(tranf2,phish,thtsh,psish)
          do 10 i=1,dim
            do 20 j=1,dim
              transf(i,j)=0.0
              do 30 k=1,dim
                transf(i,j)=transf(i,j)+tranf2(i,k)*tranf1(k,j)
30            continue
20          continue
10        continue
        endif

        do 40 i=1,dim
          do 50 j=1,dim
            shear(i,j)=0.0
            do 60 k=1,dim
              shear(i,j)=shear(i,j)+lambsh(k)*transf(k,j)*transf(k,i)
60          continue
50        continue
40      continue
  
        sh1(ipk)=shear(1,1)
        sh2(ipk)=shear(2,2)
        sh3(ipk)=shear(2,1)
        sh4(ipk)=shear(3,1)
        sh5(ipk)=shear(3,2)

        if (ichc1.ne.0) then
          v3(1)=transf(3,1)
          v3(2)=transf(3,2)
          v3(3)=transf(3,3)
          v1(1)=transf(1,1)
          v1(2)=transf(1,2)
          v1(3)=transf(1,3)
          thtsh=dacosd(v3(3))
          vnl(1)=0.0
          vnl(2)=0.0
          vnl(3)=0.0
          if (dabs(v3(1)).lt.1.0d-10) then
            vnl(1)=1.0
          else
            if (dabs(v3(2)).lt.1.0d-10) then
              vnl(2)=1.0
            else
              if (v3(1).lt.0.0) then
                vnl(1)=v3(2)
                vnl(2)=-v3(1)
              else
                vnl(1)=-v3(2)
                vnl(2)=v3(1)
              endif
            endif
          endif
          vnln=dsqrt(vnl(1)*vnl(1)+vnl(2)*vnl(2))
          vnl(1)=vnl(1)/vnln
          vnl(2)=vnl(2)/vnln
          phish=dacosd(vnl(1))
          v2(1)=v1(3)*v3(2)-v1(2)*v3(3)
          v2(2)=v1(1)*v3(3)-v1(3)*v3(1)
          v2(3)=v1(2)*v3(1)-v1(1)*v3(2)
          psish=0.0
          psih=0.0
          do 70 m=1,dim
            psish=psish+v1(m)*vnl(m)
            psih=psih+v2(m)*vnl(m)
70        continue
          if (psih.le.0.0) then
            psish=dacosd(psish)
          else
            psish=180.0-dacosd(psish)
          endif
        endif

        write(*,*)
        write(*,'(''      s_1:          '',f11.6,$)') lambsh(1)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_2:          '',f11.6,$)') lambsh(2)
        write(*,'('' km/s/Mpc'')') 
        write(*,'(''      phi_sh:       '',f11.6)') phish
        write(*,'(''      tht_sh:       '',f11.6)') thtsh
        write(*,'(''      psi_sh:       '',f11.6)') psish
        write(*,*)
        write(*,'(''      s_11:         '',f11.6,$)') sh1(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_22:         '',f11.6,$)') sh2(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_33:         '',f11.6,$)') shear(3,3)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_12:         '',f11.6,$)') sh3(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_13:         '',f11.6,$)') sh4(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,'(''      s_23:         '',f11.6,$)') sh5(ipk)
        write(*,'('' km/s/Mpc'')')
        write(*,*)
        write(1,*)
        write(1,'(''      s_1:          '',f11.6,$)') lambsh(1)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_2:          '',f11.6,$)') lambsh(2)
        write(1,'('' km/s/Mpc'')') 
        write(1,'(''      phi_sh:       '',f11.6)') phish
        write(1,'(''      tht_sh:       '',f11.6)') thtsh
        write(1,'(''      psi_sh:       '',f11.6)') psish
        write(1,*)
        write(1,'(''      s_11:         '',f11.6,$)') sh1(ipk)
        write(1,'('' km/s/Mpc'')') 
        write(1,'(''      s_22:         '',f11.6,$)') sh2(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_33:         '',f11.6,$)') shear(3,3)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_12:         '',f11.6,$)') sh3(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_13:         '',f11.6,$)') sh4(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,'(''      s_23:         '',f11.6,$)') sh5(ipk)
        write(1,'('' km/s/Mpc'')')
        write(1,*)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conkmn(ncon,lgcon,rgauss,centre,ipk)


c     -------------------------------------------------------------------
c     Subroutine `conkmn' is the main Fourier space (k-space) routine for 
c     setting up the constraint matrix conc in Fourier space for a given 
c     peak with position `centre' and a (Gaussian smoothing) scale R_G
c     (rgauss). 
c     Depending on what constraints were decided to put on the peak, 
c     information contained in the logical array lgcon, it will call 
c     several other routines for determining the relevant elements of 
c     the constraint matrix.  
c     First, it determines the common element in each of the constraint 
c     matrix elements, e^(-k^2 R_G^2/2) e^(ikr), determined by the 
c     scale and position of the peak. Then it can call successively the 
c     5 subroutines :
c       (1)  peak height
c       (2)  first derivatives, equal to 0.0 for a peak/dip
c       (3)  second derivatives:
c            curvature, axis ratios and orientation of peak
c       (4)  velocity of peak:
c            amplitude and direction
c       (5)  shear at peak:
c            2 independent eigenvalues and direction principal axes.
c     At every extra imposed constraint the integer inmcon is increased 
c     by 1.
c     ------------------------------------------------------------------


      implicit double precision (a-h,o-z)

      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  centre(dim),rgauss
      real              boxlen
      real              coskr,sinkr,kr,rx,ry,rz
      complex           conc,gskrnc
      logical           lgcon(5)

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common            /conc/conc(ngdim0,nconmax)/gskrnc/gskrnc(ngdim0)


        dk=twopi/boxlen
        akmax=float(n1)*dk

        rgauss=rgauss*boxlen
        rgs2=rgauss*rgauss/2.0

        rx=centre(1)*boxlen
        ry=centre(2)*boxlen
        rz=centre(3)*boxlen
   
c       -----------------------------------------------------------------
c       Setting up general part of constraint kernels in Fourier space, 
c       ie. e^(-k^2 R_G^2/2) e^{ikr), determined by scale (Gaussian, R_G) 
c       of peak and its position r.
c       -----------------------------------------------------------------

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            kr=ak2*ry+ak3*rz
            coskr=cos(kr)
            sinkr=sin(kr)
            gskrnc(k)=cmplx(coskr,-sinkr)*exp(-akk*rgs2)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            kr=ak1*rx+ak2*ry+ak3*rz
            coskr=cos(kr)
            sinkr=sin(kr)
            gskrnc(k)=cmplx(coskr,-sinkr)*exp(-akk*rgs2)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              kr=ak1*rx+ak2*ry+ak3*rz
              coskr=cos(kr)
              sinkr=sin(kr)
              gskrnc(k)=cmplx(coskr,-sinkr)*exp(-akk*rgs2)
30          continue         
20        continue
10      continue


c       -----------------------------------------------------------
c       Now the separate constraint matrix elements are determined, 
c       depending on your choice of constraints, as given by the 
c       elements of the logical array lgcon. 
c       lgcon(1) = .true.            peak height constraint
c       lgcon(2) = .true.            first derivative constraints
c       lgcon(3) = .true.            second derivatives constraint
c       lgcon(4) = .true.            velocity constraint
c       lgcon(5) = .true.            shear constraint
c       -----------------------------------------------------------
        if (lgcon(1)) call conk1(ncon,ipk)
        if (lgcon(2)) call conk2(ncon,ipk)
        if (lgcon(3)) call conk3(ncon,ipk)
        if (lgcon(4)) call conk4(ncon,ipk)
        if (lgcon(5)) call conk5(ncon,ipk)

1000   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conk1(ncon,ipk)


      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              boxlen
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,
     2                  gamma,rstar
      complex           conc,gskrnc

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /conc/conc(ngdim0,nconmax)/gskrnc/gskrnc(ngdim0)
      common            /sigmap/sigma0(npeakm),sigma1(npeakm),
     2                          sigma2(npeakm),sigmav(npeakm),
     2                          sigmas(2,npeakm),
     2                          gamma(npeakm),rstar(npeakm)


        ncon=ncon+1
        i1=ncon

        facsig=1.0/sigma0(ipk)

        xmax=float(n1)
        dk=twopi/boxlen
        akmax=n1*dk

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            conc(k,i1)=cmplx(facsig,0.0)*gskrnc(k)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            conc(k,i1)=cmplx(facsig,0.0)*gskrnc(k)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              conc(k,i1)=cmplx(facsig,0.0)*gskrnc(k)
30          continue         
20        continue
10      continue


1001   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conk2(ncon,ipk)


      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              boxlen
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,
     2                  gamma,rstar
      complex           conc,gskrnc

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /conc/conc(ngdim0,nconmax)/gskrnc/gskrnc(ngdim0)
      common            /sigmap/sigma0(npeakm),sigma1(npeakm),
     2                          sigma2(npeakm),sigmav(npeakm),
     2                          sigmas(2,npeakm),     
     2                          gamma(npeakm),rstar(npeakm)    


        i1=ncon+1
        i2=ncon+2
        i3=ncon+3
        ncon=ncon+3

        facsig=1.0/sigma1(ipk)

        xmax=float(n1)
        dk=twopi/boxlen
        akmax=n1*dk

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            fact2=-facsig*ak2
            fact3=-facsig*ak3
            conc(k,i1)=cmplx(0.0,0.0)
            conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
            conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            fact1=-facsig*ak1
            fact2=-facsig*ak2
            fact3=-facsig*ak3
            conc(k,i1)=cmplx(0.0,fact1)*gskrnc(k)
            conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
            conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              fact1=-facsig*ak1
              fact2=-facsig*ak2
              fact3=-facsig*ak3
              conc(k,i1)=cmplx(0.0,fact1)*gskrnc(k)
              conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
              conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)
30          continue         
20        continue
10      continue


1002   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conk3(ncon,ipk)


      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              boxlen
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,gamma,rstar
      complex           conc,gskrnc

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /conc/conc(ngdim0,nconmax)/gskrnc/gskrnc(ngdim0)
      common            /sigmap/sigma0(npeakm),sigma1(npeakm),
     2                          sigma2(npeakm),sigmav(npeakm),
     2                          sigmas(2,npeakm),     
     2                          gamma(npeakm),rstar(npeakm)


        i1=ncon+1
        i2=ncon+2
        i3=ncon+3
        i4=ncon+4
        i5=ncon+5
        i6=ncon+6
        ncon=ncon+6

        facsig=1.0/sigma2(ipk)

        xmax=float(n1)
        dk=twopi/boxlen
        akmax=n1*dk

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            fact1=0.0
            fact2=-ak2*ak2*facsig
            fact3=-ak3*ak3*facsig
            fact4=0.0
            fact5=0.0
            fact6=-ak2*ak3*facsig
            conc(k,i1)=cmplx(fact1,0.0)
            conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
            conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
            conc(k,i4)=cmplx(fact4,0.0)
            conc(k,i5)=cmplx(fact5,0.0)
            conc(k,i6)=cmplx(fact6,0.0)*gskrnc(k)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            fact1=-ak1*ak1*facsig
            fact2=-ak2*ak2*facsig
            fact3=-ak3*ak3*facsig
            fact4=-ak1*ak2*facsig
            fact5=-ak1*ak3*facsig
            fact6=-ak2*ak3*facsig
            conc(k,i1)=cmplx(fact1,0.0)*gskrnc(k)
            conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
            conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
            conc(k,i4)=cmplx(fact4,0.0)*gskrnc(k)
            conc(k,i5)=cmplx(fact5,0.0)*gskrnc(k)
            conc(k,i6)=cmplx(fact6,0.0)*gskrnc(k)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              fact1=-ak1*ak1*facsig
              fact2=-ak2*ak2*facsig
              fact3=-ak3*ak3*facsig
              fact4=-ak1*ak2*facsig
              fact5=-ak1*ak3*facsig
              fact6=-ak2*ak3*facsig
              conc(k,i1)=cmplx(fact1,0.0)*gskrnc(k)
              conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
              conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
              conc(k,i4)=cmplx(fact4,0.0)*gskrnc(k)
              conc(k,i5)=cmplx(fact5,0.0)*gskrnc(k)
              conc(k,i6)=cmplx(fact6,0.0)*gskrnc(k)
30          continue         
20        continue
10      continue


1003   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conk4(ncon,ipk)


      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,gamma,rstar
      complex           conc,gskrnc

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common            /conc/conc(ngdim0,nconmax)
      common            /gskrnc/gskrnc(ngdim0)
      common            /sigmap/sigma0(npeakm),sigma1(npeakm),
     2                          sigma2(npeakm),sigmav(npeakm),
     2                          sigmas(2,npeakm),    
     2                          gamma(npeakm),rstar(npeakm)


        i1=ncon+1
        i2=ncon+2
        i3=ncon+3
        ncon=ncon+3

        facsig=1.0/sigmav(ipk)

        xmax=float(n1)
        dk=twopi/boxlen
        akmax=n1*dk

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            if (k.eq.1) then
              fact1=0.0
              fact2=0.0
              fact3=0.0
            else
              fact1=0.0
              fact2=-hf0*ak2*facsig/akk
              fact3=-hf0*ak3*facsig/akk
            endif
            conc(k,i1)=cmplx(0.0,fact1)*gskrnc(k)
            conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
            conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            fact1=-hf0*ak1*facsig/akk
            fact2=-hf0*ak2*facsig/akk
            fact3=-hf0*ak3*facsig/akk
            conc(k,i1)=cmplx(0.0,fact1)*gskrnc(k)
            conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
            conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              fact1=-hf0*ak1*facsig/akk
              fact2=-hf0*ak2*facsig/akk
              fact3=-hf0*ak3*facsig/akk
              conc(k,i1)=cmplx(0.0,fact1)*gskrnc(k)
              conc(k,i2)=cmplx(0.0,fact2)*gskrnc(k)
              conc(k,i3)=cmplx(0.0,fact3)*gskrnc(k)
30          continue         
20        continue
10      continue


1004   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine conk5(ncon,ipk)


      integer           dim
      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (dim=3,n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  sigma0,sigma1,sigma2,sigmav,sigmas,gamma,rstar
      double precision  facsig
      complex           conc,gskrnc

      common            /box/boxlen
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common            /conc/conc(ngdim0,nconmax)
      common            /gskrnc/gskrnc(ngdim0)
      common            /sigmap/sigma0(npeakm),sigma1(npeakm),
     2                          sigma2(npeakm),sigmav(npeakm),
     2                          sigmas(2,npeakm),    
     2                          gamma(npeakm),rstar(npeakm)
    

        i1=ncon+1
        i2=ncon+2
        i3=ncon+3
        i4=ncon+4
        i5=ncon+5
        ncon=ncon+5

        facsig=1.0/(hf0*sigma0(ipk))

        xmax=float(n1)
        dk=twopi/boxlen
        akmax=n1*dk

        do 10 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax

c Do k1=1 separately
            k23=(k2-1+(k3-1)*n1)
            akk=ak2*ak2+ak3*ak3
            k=1+n12*k23
            if (k.eq.1) then
              fact1=0.0
              fact2=0.0
              fact3=0.0
              fact4=0.0
              fact5=0.0
            else
              fact1=facsig*hf0/3.0
              fact2=-hf0*facsig*(ak2*ak2/akk-1.0/3.0)
              fact3=0.0
              fact4=0.0
              fact5=-hf0*facsig*ak2*ak3/akk
            endif
            conc(k,i1)=cmplx(fact1,0.0)*gskrnc(k)
            conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
            conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
            conc(k,i4)=cmplx(fact4,0.0)*gskrnc(k)
            conc(k,i5)=cmplx(fact5,0.0)*gskrnc(k)

c Do k1=n12+1, Nyquist frequency
            ak1=float(n12)*dk
            akk=ak1*ak1+ak2*ak2+ak3*ak3
            k=ngr2+1+k23
            fact1=-hf0*facsig*(ak1*ak1/akk-1.0/3.0)
            fact2=-hf0*facsig*(ak2*ak2/akk-1.0/3.0)
            fact3=-hf0*facsig*ak1*ak2/akk
            fact4=-hf0*facsig*ak1*ak3/akk
            fact5=-hf0*facsig*ak2*ak3/akk
            conc(k,i1)=cmplx(fact1,0.0)*gskrnc(k)
            conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
            conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
            conc(k,i4)=cmplx(fact4,0.0)*gskrnc(k)
            conc(k,i5)=cmplx(fact5,0.0)*gskrnc(k)

c Do k1=2,n12
            k23=(k2-1+(k3-1)*n1)*n12
            do 30 k1=2,n12
              ak1=(k1-1)*dk
              akk=ak1*ak1+ak2*ak2+ak3*ak3
              k=k1+k23
              fact1=-hf0*facsig*(ak1*ak1/akk-1.0/3.0)
              fact2=-hf0*facsig*(ak2*ak2/akk-1.0/3.0)
              fact3=-hf0*facsig*ak1*ak2/akk
              fact4=-hf0*facsig*ak1*ak3/akk
              fact5=-hf0*facsig*ak2*ak3/akk
              conc(k,i1)=cmplx(fact1,0.0)*gskrnc(k)
              conc(k,i2)=cmplx(fact2,0.0)*gskrnc(k)
              conc(k,i3)=cmplx(fact3,0.0)*gskrnc(k)
              conc(k,i4)=cmplx(fact4,0.0)*gskrnc(k)
              conc(k,i5)=cmplx(fact5,0.0)*gskrnc(k)
30          continue         
20        continue
10      continue


1005   return
      end


c     ********************************************************************
c     ********************************************************************


      subroutine concmp(g,grho,ncon,icmp)


c     ---------------------------------------------------------------------
c     Comparison of the values of the constraints for the field rhoc with 
c     the desired values. 
c     If (icmp.eq.1) the field rhoc is the unconstrained realization.
c     If (icmp.eq.2) the field rhoc is either the mean field or the 
c     constrained realization. 
c
c     The constraints conc are Fourier transformed so that the convolutions 
c     with the covariance matrix required to compute the mean field and the 
c     covariance of the constraints becomes multiplication by the power 
c     spectrum.
c     ---------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon,icmp
      complex           conc,rhoc
      complex*16        conch
      real              g(nconmax),grho(nconmax)

      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /rhoc/rhoc(ngdim0)
      common            /conc/conc(ngdim0,nconmax)


        if (icmp.eq.1) then
          write(*,*) 'Comparison sampled and desired constraints: '
          write(*,*)
          write(*,1002) 'Sampled','Desired'
          write(*,*)
          write(1,*) 'Comparison sampled and desired constraints: '
          write(1,*)
          write(1,1002) 'Sampled','Desired'
          write(1,*)
        endif
        if (icmp.eq.2) then
          write(*,*) 'Checking the values of the constraints for',
     2       ' constrained realization: '
          write(*,*)
          write(*,1002) 'Final','Desired'
          write(*,*)
          write(1,*) 'Checking the values of the constraints for',
     2       ' constrained realization: '
          write(1,*)
          write(1,1002) 'Final','Desired'
          write(1,*)
        endif

c       -----------------------------------------------------------
c       Compute the values of the constraints for the field rhoc.
c       realization (g0):
c       -----------------------------------------------------------
	do 50 i=1,ncon
	  sigma=0.0
c         -------------------------------------------------------------
c         Take twice the real part for harmonics with k1 > 0 to account 
c         for the conjugate harmonics with k1 < 0.
c         -------------------------------------------------------------
	  do 40 k=1,ngdim2
            conch=conjg(conc(k,i))
	    sigma=sigma+2.0*conch*rhoc(k)
40	  continue
c         --------------------------------------------------------
c         Correct for k1=1 and k1=n12+1, which are self-conjugate.
c         --------------------------------------------------------
	  do 46 k3=1,n1
	    do 43 k2=1,n1
	      k23=k2-1+(k3-1)*n1
	      k=1+n12*k23
              conch=conjg(conc(k,i))
	      sigma=sigma-conch*rhoc(k)
	      k=ngr2+1+k23
              conch=conjg(conc(k,i))
	      sigma=sigma-conch*rhoc(k)
43	    continue
46	  continue

	  grho(i)=sigma

	  write(*,1001) ' Constraint ',i,': ',real(sigma),g(i)
	  write(1,1001) ' Constraint ',i,': ',real(sigma),g(i)
50	continue

1001	format(1x,a,i3,a,f12.7,3x,f12.7)
1002    format(19x,a10,5x,a10)

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine covmat(qinv,ncon)


c     --------------------------------------------------------------------
c     Calculate the covariance matrix q of the constraints (\xi_ij in 
c     Hoffman-Ribak) and subsequently determining the inverse (\xi_ij^-1), 
c     qinv. During inversion the matrix q is destroyed.
c     --------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      double precision  twopi
      parameter         (twopi=6.283185307179586d0)
      parameter         (n0=128,npeakm=3,nconmax=18*npeakm)
      parameter         (n02=n0/2,n0n0=n0*n0,ngrid0=n0n0*n0)
      parameter         (ngr0=ngrid0/2,ngdim0=ngr0+n0n0,ngdimm=2*ngdim0)

      integer           ncon
      double precision  power
      double precision  q(nconmax,nconmax),qinv(nconmax,nconmax)
      complex           conc
      complex*16        conch1,conch2

      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim
      common            /power/power(ngdim0)
      common            /conc/conc(ngdim0,nconmax)


c       -----------------------------------------------------    
c       Calculate covariance matrix (\xi_ij in Hoffman-Ribak)
c       -----------------------------------------------------

	do 80 i2=1,ncon
	  do 70 i1=1,ncon
	    q(i1,i2)=0.0d0
	    do 60 k=1,ngdim2
              conch1=conjg(conc(k,i1))
              conch2=conc(k,i2)
              q(i1,i2)=q(i1,i2)+2.0*conch1*conch2*power(k)
60	    continue
c           --------------------------------
c           Correct q for k1=1 and k1=n12+1.
c           --------------------------------
	    do 66 k3=1,n1
	      do 63 k2=1,n1
	        k23=k2-1+(k3-1)*n1
	        k=1+n12*k23
                conch1=conjg(conc(k,i1))
                conch2=conc(k,i2)
                q(i1,i2)=q(i1,i2)-conch1*conch2*power(k)
	        k=ngr2+1+k23
                conch1=conjg(conc(k,i1))
                conch2=conc(k,i2)
                q(i1,i2)=q(i1,i2)-conch1*conch2*power(k)
63	      continue
66	    continue
70	  continue
80	continue

c       ------------------------------------------------------------------
c       Invert double precision matrix q to give qinv:  
c       Warning: q is destroyed! (this gives \xi_ij^{-1} in Hoffman-Ribak)
c       ------------------------------------------------------------------

        call matinv(q,ncon,nconmax,qinv)         


       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine conchs(qinv,g,g0,ncon)


c     ----------------------------------------------------------------------
c     Compute chi-squared for the constraints and printing this information.
c     ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter         (n0=128,npeakm=3,nconmax=18*npeakm)

      integer           ncon
      double precision  chisq,chisq0
      double precision  qinv(nconmax,nconmax)
      real              g(nconmax),g0(nconmax)


	chisq=0.0
	chisq0=0.0
	do 100 i2=1,ncon
	  do 90 i1=1,ncon
	    chisq=chisq+g(i1)*qinv(i1,i2)*g(i2)
	    chisq0=chisq0+g0(i1)*qinv(i1,i2)*g0(i2)
90	  continue
100	continue

	write(*,'(1x,a20,i3,a13)') 'Chi-squared for the ',ncon,
     2                             ' constraints:'
	write(*,1003) '   Sampled = ',chisq0
        write(*,1003) '   Desired = ',real(chisq)
	write(1,'(1x,a20,i3,a13)') 'Chi-squared for the ',ncon,
     2                             ' constraints:'
	write(1,1003) '   Sampled = ',chisq0
        write(1,1003) '   Desired = ',real(chisq)
1003    format(16x,a,f12.7)

	chisq0=chisq

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine spparam(rgauss,sigma0,sigma1,sigma2,gamma,rstar)


      implicit double precision (a-h,o-z)
      parameter (pi=3.141592653,twopi=2.*pi)

      external  pgssig

      double precision   rgauss,sigma0,sigma1,sigma2,gamma,rstar

      common /cosmos/h0,hubbl0,omega,omegahh,hf0
      common /spectr/rth,rgs,primn,amplit,bias,pindex,ispec


        pindex=2.0*0.0
        sigma0=pgssig(rgauss)
        pindex=2.0*1.0
        sigma1=pgssig(rgauss)
        pindex=2.0*2.0
        sigma2=pgssig(rgauss)

        gamma=sigma1**2/(sigma2*sigma0)
        rstar=sqrt(3.0)*sigma1/sigma2

       return
      end


C     ************************************************************************
C     ************************************************************************

      
      double precision function curvav(gamma,pnu)

      implicit double precision (a-h,o-z)

      external  theta


        wg=gamma*pnu
        curvav=wg+theta(gamma,wg)

1000   return
      end



C     ************************************************************************
C     ************************************************************************

      
      double precision function theta(gamma,wg)

      implicit double precision (a-h,o-z)
    

        gamma2=gamma*gamma
        gamma4=gamma2*gamma2

        wg2=wg/2.0

        help1=3.0*(1.0-gamma2)
        help2=1.216-0.9*gamma4
        help3=-(gamma/2.0)*wg2*wg2
        help3=help1+help2*exp(help3)

        help4=help1+0.45+wg2*wg2
        help4=sqrt(help4)+wg2

        theta=help3/help4


1000   return
      end


C     ************************************************************************
C     ************************************************************************


      subroutine xfunc(gamma,pnu)


c     -------------------------------------------------------------
c     Determining the distribution function of x and its cumulative
c     until the value of the cumulative function is equal to xperc. 
c     The corresponding value is xcfunc.  
c     -------------------------------------------------------------

      implicit double precision (a-h,o-z)

      external           pxpnuf

      integer            nxmax
      parameter          (nxmax=5000)

      double precision   gamma,pnu
      double precision   pxpnuc,pxpnu0,pxpnu1,xlow,xupp

      common /curvat/xlow,xdel,pxpnuc(nxmax),npx


        npx=nxmax
        xlow=0.0
        xupp=10.0
        xdel=(xupp-xlow)/float(npx-1)

        x=xlow
        pxpnu0=pxpnuf(x,gamma,pnu)
        pxpnuc(1)=0.0
        do 10 i=2,npx
          x=xlow+float(i-1)*xdel
          pxpnu1=pxpnuf(x,gamma,pnu)
          pxpnuc(i)=pxpnuc(i-1)+(pxpnu1+pxpnu0)*xdel/2.0
          pxpnu0=pxpnu1    
10      continue

       return
      end


C     ************************************************************************
C     ************************************************************************


      subroutine shfunc(xpeak)


      implicit double precision (a-h,o-z)

      external  pepx

      parameter (nemax=257,npmax=257,nepmax=nemax*npmax)

      integer             ne,np,nep,ihshp
      double precision    xpeak,e,p,emin,emax,pmin,pmax
      double precision    pshp(nemax,npmax),pshpcm,
     *                    hshp

      common /shape/pshpcm(nepmax),hshp(nepmax),ihshp(2,nepmax),nep
      common /shapep/emean,pmean,emax,pmax,eppred,pppred


        ne=101
        np=101

        emin=0.0
        emax=0.5
        pmin=-0.25
        pmax=0.25

        dele=(emax-emin)/float(ne-1)
        delp=(pmax-pmin)/float(np-1)

        do 10 i=1,ne
          e=emin+float(i-1)*dele
          do 20 j=1,np
            p=pmin+float(j-1)*delp
            pshp(i,j)=pepx(e,p,xpeak)
20        continue
10      continue

C ---------------------------------------------------------------------------

C       ------------------------------------------------
C       Determining the cumulative distribution function
C       ------------------------------------------------

        do 30 i=1,nepmax
          hshp(i)=0.0
          ihshp(1,i)=0
          ihshp(2,i)=0
30      continue

        iep=0
        do 40 i=1,ne
          do 50 j=1,np
            iep=iep+1
            hshp(iep)=pshp(i,j)
            ihshp(1,iep)=i
            ihshp(2,iep)=j
50        continue
40      continue

        nep=ne*np
        delaep=dele*delp
        call sort(hshp,ihshp,nep,store,istore)

        pshpcm(1)=hshp(1)*delaep
        do 60 i=2,nep
          pshpcm(i)=pshpcm(i-1)+hshp(i)*delaep
60      continue

C ---------------------------------------------------------------------------

C       ---------------------------------------
C       Determining some statistical parameters
C       ---------------------------------------

        is1=ihshp(1,1)
        is2=ihshp(2,1)
        emax=emin+float(is1-1)*dele
        pmax=pmin+float(is2-1)*delp

        emean=0.0
        pmean=0.0
        do 70 i=1,ne
          e=emin+float(i-1)*dele
          do 80 j=1,np
            p=pmin+float(j-1)*delp
            emean=emean+e*pshp(i,j)*delaep
            pmean=pmean+p*pshp(i,j)*delaep
80        continue
70      continue

        f1=6.0/(5.0*xpeak*xpeak)+1.0
        eppred=sqrt(f1)*xpeak*sqrt(5.0)
        eppred=1.0/eppred
        pppred=(f1**2)*5.0*xpeak*xpeak*xpeak*xpeak
        pppred=6.0/pppred

       return
      end


C     ************************************************************************
C     ************************************************************************


      subroutine vcfunc(sigmav)


c     ----------------------------------------------------------------
c     Determining the cumulative distribution function of the peculiar 
c     velocity of a peak. 
c     ----------------------------------------------------------------

      implicit double precision (a-h,o-z)

      external           pvfunc

      integer            nvmax
      parameter          (nvmax=5000)

      double precision   pv0,pv1,velpkc,vlow,vupp,sigmav

      common /pecvel/vlow,vdel,velpkc(nvmax),npv


        npv=nvmax
        vlow=0.0
        vupp=10.0*sigmav
        vdel=(vupp-vlow)/float(npv-1)

        v=vlow
        pv0=pvfunc(v,sigmav)
        velpkc(1)=0.0
        do 10 i=2,npv
          v=vlow+float(i-1)*vdel
          pv1=pvfunc(v,sigmav)
          velpkc(i)=velpkc(i-1)+(pv1+pv0)*vdel/2.0
          pv0=pv1
10      continue

       return
      end


C     ************************************************************************
C     ************************************************************************


      double precision function pxpnuf(x,gamma,pnu)


      implicit double precision (a-h,o-z)

      external fffx,gggn

      parameter (pi=3.141592653,twopi=2.*pi)


        xstar=gamma*pnu
  
        fc1=2.0*(1.0-gamma*gamma)

        f1=(x-xstar)*(x-xstar)/fc1
        f1=exp(-f1)
        f2=fffx(x)
        f3=sqrt(pi*fc1)
        f4=gggn(gamma,xstar)

        pxpnuf=f1*f2/f3/f4

1000   return
      end


C     ************************************************************************
C     ************************************************************************

      
      double precision function gggn(gamma,wg)

      implicit double precision (a-h,o-z)

      external  ag,bg,c1g,c2g,c3g

        
        help1=wg*wg*wg
        help2=-3.0*gamma*gamma*wg
        help3=(bg(gamma)*wg*wg+c1g(gamma))
        help4=ag(gamma)*wg*wg
        help4=exp(-help4)
        help1=help1+help2+help3*help4
        help5=c2g(gamma)
        help6=c3g(gamma)*wg
        help6=exp(-help6)
        help5=1.0+help5*help6

        gggn=help1/help5

1000   return
      end

C     ************************************************************************
C     ************************************************************************

      
      double precision function ag(gamma)

      implicit double precision (a-h,o-z)

        
        ag=9.0-5.0*gamma*gamma
        ag=5.0/(2.0*ag)

1000   return
      end

C     ************************************************************************
C     ************************************************************************

      
      double precision function bg(gamma)

      implicit double precision (a-h,o-z)

      parameter (pi=3.141592653)


        
        bg=9.0-5.0*gamma*gamma
        bg=bg**2.5
        bg=((10.0*pi)**0.5)*bg
        bg=432.0/bg

1000   return
      end

C     ************************************************************************
C     ************************************************************************

      
      double precision function c1g(gamma)

      implicit double precision (a-h,o-z)

        
        c1g=(1.0-gamma*gamma)**(5.72)
        c1g=1.84+1.13*c1g

1000   return
      end


C     ************************************************************************
C     ************************************************************************

      
      double precision function c2g(gamma)

      implicit double precision (a-h,o-z)



        c2g=6.51*gamma*gamma
        c2g=8.91+1.27*exp(c2g)        

1000   return
      end

C     ************************************************************************
C     ************************************************************************

      
      double precision function c3g(gamma)

      implicit double precision (a-h,o-z)

        
        c3g=1.05*gamma*gamma
        c3g=2.58*exp(c3g)

1000   return
      end


C     ************************************************************************
C     ************************************************************************


      double precision function fffx(x)


      implicit double precision (a-h,o-z)

      external   erf

      parameter (pi=3.141592653,twopi=2.*pi)


         x1=sqrt(2.5)*x
         x2=x1/2.0
         f1=(erf(x1)+erf(x2))/2.0
         f1=(x*x*x-3.0*x)*f1
     
         fc1=sqrt(2.0/(5.0*pi))
         f2=31.0*x*x/4.0+8.0/5.0
         f3=5.0*x*x/8.0
         f3=exp(-f3)
         f4=x*x/2.0-8.0/5.0
         f5=5.0*x*x/2.0
         f5=exp(-f5)
         f6=fc1*(f2*f3+f4*f5)

         fffx=f1+f6

1000    return
       end


C     ************************************************************************
C     ************************************************************************


      double precision function pepx(e,p,x)


      implicit double precision (a-h,o-z)

      external fffx,wep

      parameter (pi=3.141592653,twopi=2.0*pi)


        fc1=(5.0**2.5)*9.0/sqrt(twopi)

        f1=x**8.0
        f2=fffx(x)
        f3=2.5*x*x*(3.0*e*e+p*p)
        f3=exp(-f3)
        f4=wep(e,p)

        pepx=fc1*f1*f3*f4/f2

1000   return
      end


C     *************************************************************************
C     *************************************************************************


      double precision function wep(e,p)


      implicit double precision (a-h,o-z)

      external chiep


        f1=e
        f2=e*e-p*p
        f3=1.0-2.0*p
        f4=(1.0+p)*(1.0+p)-9.0*e*e
        f5=chiep(e,p)

        wep=f1*f2*f3*f4*f5

       return
      end


C     *************************************************************************
C     *************************************************************************


      double precision function chiep(e,p)

     
      implicit double precision (a-h,o-z)


        chiep=0.0

        if ((e.ge.0.0).and.(e.le.0.25)) then
          if ((p.ge.-e).and.(p.le.e)) then
            chiep=1.0
          endif
        endif

        if ((e.ge.0.25).and.(e.le.0.5)) then
          if ((p.ge.(3.0*e-1.0)).and.(p.le.e)) then
            chiep=1.0
          endif
        endif

       return
      end


C     *************************************************************************
C     *************************************************************************


      double precision function pvfunc(v,sigmav)      


      implicit double precision (a-h,o-z)

      double precision   twopi
      parameter          (twopi=6.283185307)      

      double precision   v,sigmav        
      double precision   fc1,fc2

    
        fc1=3.0*v*v/(2.0*sigmav*sigmav)
        fc1=v*v*exp(-fc1)
        fc2=(twopi*sigmav*sigmav/3.0)**1.5
        pvfunc=2.0*twopi*fc1/fc2

       return
      end


C     *************************************************************************
C     *************************************************************************


      double precision function pshpnu(ev,pv,pnu)


      implicit double precision (a-h,o-z)

      double precision  pi,twopi
      parameter         (pi=3.141592653d0,twopi=2.0*pi)


        v=pnu*ev
        w=pnu*pv

        fc0=twopi
        fc1=15.0*15.0*dsqrt(5.0d0)/(2.0*dsqrt(fc0))
        fc2=2.0*v*(v*v-w*w)
        fc3=2.5*w*w
        fc4=7.5*v*v
        fc5=exp(-fc3-fc4)
 
        pshpnu=fc1*fc2*fc3*fc4*fc5

1000   return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine shaprn(e,p,x)


      implicit double precision (a-h,o-z)

      external pepx,chiep

      parameter (nemax=257,npmax=257,nepmax=nemax*npmax)

      integer             nep,ihshp
      real                xr
      double precision    pshpcm,hshp

      common /shape/pshpcm(nepmax),hshp(nepmax),ihshp(2,nepmax),nep
      common /shapep/emean,pmean,emax,pmax,eppred,pppred
      


        prmax=1.0001*pepx(emax,pmax,x)

1110    call randa(xr)
        e=0.5*xr
        call randa(xr)
        p=-0.25+0.75*xr

        if ((e.ge.0.0).and.(e.le.0.25)) then
          if ((p.lt.-e).or.(p.gt.e)) goto 1110
        endif
        if ((e.ge.0.25).and.(e.le.0.5)) then
          if ((p.lt.(3.0*e-1.0)).and.(p.gt.e)) goto 1110
        endif

        call randa(xr)
        pran=prmax*xr
        prob=pepx(e,p,x)
        if (pran.gt.prob) goto 1110 

       return
      end


C     ************************************************************************
C     ************************************************************************


      subroutine axep(a32,a31,e,p)


      double precision   a31,a32,e,p

        a32=a32*a32
        a31=a31*a31
  
        e=(a31-1.0)/(2.0*(1.0+a31+a32))
        p=(a31-2.0*a32+1.0)/(2.0*(1.0+a31+a32))

       return
      end


C     ************************************************************************
C     ************************************************************************


      subroutine epax(e,p,a32,a31)


      double precision   a31,a32,e,p

        a32=(1.0-2.0*p)/(-3.0*e+p+1.0)
        a31=(3.0*e+p+1.0)/(-3.0*e+p+1.0)
        a32=sqrt(a32)
        a31=sqrt(a31)

       return
      end


C     *************************************************************************
C     *************************************************************************


      double precision function erf(x)

c     ---------------------------------------------------------
c     Returns the error function erf(x) with fractional error  
c     everywhere less than 1.2 x 10^-7. From Numerical Recipes, 
c     pg. 164. Based on Chebysev fitting.
c     ---------------------------------------------------------

      implicit double precision (a-h,o-z)

   
        z=abs(x)
        t=1./(1.+0.5*z)

        erfc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     2       t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     2       t*(1.48851587+t*(-.82215223+t*.17087277)))))))))      
        erf=1.-erfc

        if (x.lt.0.) erf=-erf

       return
      end


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE EULER(REULER,PHI,THETA,PSI)


      DOUBLE PRECISION   REULER(3,3),PHI,THETA,PSI
      DOUBLE PRECISION   COSPHI,SINPHI,COSTHT,SINTHT,COSPSI,SINPSI


        COSPHI=DCOSD(PHI)
        SINPHI=DSIND(PHI)
        COSTHT=DCOSD(THETA)
        SINTHT=DSIND(THETA)
        COSPSI=DCOSD(PSI)
        SINPSI=DSIND(PSI)

        REULER(1,1)=COSPSI*COSPHI-COSTHT*SINPHI*SINPSI
        REULER(2,1)=-SINPSI*COSPHI-COSTHT*SINPHI*COSPSI
        REULER(3,1)=SINTHT*SINPHI
        REULER(1,2)=COSPSI*SINPHI+COSTHT*COSPHI*SINPSI
        REULER(2,2)=-SINPSI*SINPHI+COSTHT*COSPHI*COSPSI
        REULER(3,2)=-SINTHT*COSPHI
        REULER(1,3)=SINPSI*SINTHT
        REULER(2,3)=COSPSI*SINTHT
        REULER(3,3)=COSTHT 

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE EULINV(REULIV,PHI,THETA,PSI)


      DOUBLE PRECISION   REULIV(3,3),PHI,THETA,PSI
      DOUBLE PRECISION   COSPHI,SINPHI,COSTHT,SINTHT,COSPSI,SINPSI


        COSPHI=DCOSD(PHI)
        SINPHI=DSIND(PHI)
        COSTHT=DCOSD(THETA)
        SINTHT=DSIND(THETA)
        COSPSI=DCOSD(PSI)
        SINPSI=DSIND(PSI)

        REULIV(1,1)=COSPSI*COSPHI-COSTHT*SINPHI*SINPSI
        REULIV(2,1)=COSPSI*SINPHI+COSTHT*COSPHI*SINPSI
        REULIV(3,1)=SINTHT*SINPSI
        REULIV(1,2)=-SINPSI*COSPHI-COSTHT*SINPHI*COSPSI
        REULIV(2,2)=-SINPSI*SINPHI+COSTHT*COSPHI*COSPSI
        REULIV(3,2)=SINTHT*COSPSI
        REULIV(1,3)=SINTHT*SINPHI
        REULIV(2,3)=-SINTHT*COSPHI
        REULIV(3,3)=COSTHT         

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE VTRANS(TRANS,PHI,THETA)


      DOUBLE PRECISION   TRANS(3,3),PHI,THETA
      DOUBLE PRECISION   COSPHI,SINPHI,COSTHT,SINTHT


        COSPHI=DCOSD(PHI)
        SINPHI=DSIND(PHI)
        COSTHT=DCOSD(THETA)
        SINTHT=DSIND(THETA)

        TRANS(1,1)=COSTHT*COSPHI
        TRANS(2,1)=-SINPHI
        TRANS(3,1)=SINTHT*COSPHI
        TRANS(1,2)=COSTHT*SINPHI
        TRANS(2,2)=COSPHI
        TRANS(3,2)=SINTHT*SINPHI
        TRANS(1,3)=-SINTHT
        TRANS(2,3)=0.0
        TRANS(3,3)=COSTHT

       RETURN
      END


c     *************************************************************************
c     *************************************************************************


      double precision function fpeebl(omega)

c  Computes the velocity function f(Omega)=d\log d1/d\log x.  Although this
c  is often approximated by Omega**0.6, the exact expression is used here.
c  See Peebles LSSU section 14.

      implicit double precision (a-h,o-z)


  	if (omega.eq.1.0) then
	  fpeebl=1.0
	  return
	else if (omega.le.0.0) then
	  write(*,*) 'Omega <= 0 in fpeebl!'
	  stop
  	end if

	x=1.0/omega-1.0
	if (x.gt.0.0) then
	  sqrx=sqrt(x)
	  sqrx1=sqrt(1.0+x)
	  sinh0=2.0*sqrx*sqrx1
	  eta=2.0*log(sqrx+sqrx1)
	  d1=0.75*sinh0*(sinh0-eta)/(x*x)-2.0
	  detadlx=2.0*x/sinh0
	  fpeebl=((1.0+2.0*x)/sinh0+2.0*x/(sinh0-eta))*detadlx-2.0
	  fpeebl=fpeebl*(1.0+2.0/d1)
	  return
	else 
          if (x.lt.0.0) then
	    sqrx=sqrt(-x)
	    sqrx1=sqrt(1.0+x)
	    sin0=2.0*sqrx*sqrx1
	    theta=2.0*atan2(sqrx,sqrx1)
	    d1=0.75*sin0*(theta-sin0)/(x*x)-2.0
	    dthedlx=-2.0*x/sin0
	    fpeebl=((1.0+2.0*x)/sin0-2.0*x/(theta-sin0))*dthedlx-2.0
	    fpeebl=fpeebl*(1.0+2.0/d1)
	    return
	  else
	    fpeebl=1.0
	    return
	  endif
        endif

      end


C     *************************************************************************
C     *************************************************************************


      subroutine spectrum


c     -------------------------
c     Setting up power spectrum
c     -------------------------

      implicit double precision (a-h,o-z)

      parameter         (n0=128,pi=3.141592653,twopi=2.*pi)
      parameter         (n0n0=n0*n0,ngr0=n0n0*n0/2,ngdim0=ngr0+n0n0)

      external          trans1,trans2,trans3,trans4,trans5,trans6
      real              boxlen
      double precision  h0,hubbl0,omega0,omegahh,hf0
      double precision  power

      common            /cosmos/h0,hubbl0,omega0,omegahh,hf0
      common            /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common            /box/boxlen
      common            /power/power(ngdim0)
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


        dk=twopi/boxlen

        write(*,*) 'Choice of fluctuation spectrum:'
        write(*,*)
        write(*,*) '   spectrum:                  ispec:'         
        write(*,*) '   power-law                    1'
        write(*,*) '   CDM, adiabatic, DEFW         2'
        write(*,*) '   CDM, adiabatic, BBKS         3'
        write(*,*) '   HDM, adiabatic, BBKS         4'
        write(*,*) '   WDM, adiabatic, BBKS         5'
        write(*,*) '   CDM, isocurvature, BBKS      6'
        write(*,*) 
        write(*,'('' ispec: '',$)')
        read(*,*) ispec
        write(*,*)
        write(1,*) 'Choice of fluctuation spectrum:'
        write(1,*)
        write(1,*) '   spectrum:                  ispec:'         
        write(1,*) '   power-law                    1'
        write(1,*) '   CDM, adiabatic, DEFW         2'
        write(1,*) '   CDM, adiabatic, BBKS         3'
        write(1,*) '   HDM, adiabatic, BBKS         4'
        write(1,*) '   WDM, adiabatic, BBKS         5'
        write(1,*) '   CDM, isocurvature, BBKS      6'
        write(1,*) 
        write(1,'('' ispec: '',$)')
        write(1,'(i3)') ispec
        write(1,*)

        if (ispec.ne.1) then
          write(*,'('' index N primordial spectrum [-3,2]: '',$)')
          write(1,'('' index N primordial spectrum [-3,2]: '',$)')
        else
          write(*,'('' Index n of power spectrum [-3,1]: '',$)')
          write(1,'('' Index n of power spectrum [-3,1]: '',$)')
        endif
        read(*,*) primn
        write(1,'(f6.1)') primn
        write(*,*)
        write(1,*)

        write(*,'('' Normalization of spectrum: '')')
        write(*,'('' Tophat window radius (h^{-1} Mpc):    '',$)')
        read(*,*) rth
        write(1,'('' Normalization of spectrum: '')')
        write(1,'('' Tophat window radius (h^{-1} Mpc):    '',$)')
        write(1,'(f10.3)') rth
        rth=rth/h0
        write(*,'('' Galaxy number sigma at Tophat radius: '',$)')
        read(*,*) sigma0
        write(*,'('' bias factor b:                        '',$)')
        read(*,*) bias
        write(1,'('' Galaxy number sigma at Tophat radius: '',$)')
        write(1,'(f10.3)') sigma0
        write(1,'('' bias factor b:                        '',$)')
        write(1,'(f10.3)') bias
        ak1=twopi/rth

        call normal(ispec,ak1,sigma0)
        
        if (ispec.eq.1) call pspect(trans1,dk)
        if (ispec.eq.2) call pspect(trans2,dk)
        if (ispec.eq.3) call pspect(trans3,dk)
        if (ispec.eq.4) call pspect(trans4,dk)
        if (ispec.eq.5) call pspect(trans5,dk)
        if (ispec.eq.6) call pspect(trans6,dk)

       return
      end


C     *************************************************************************
C     *************************************************************************


      subroutine pspect(trans,dk)


c     ---------------------------------------------------------------------
c     pspect computes the power in n1*n1*n1 frequency bands of 1-D width dk.
c     The power spectrum is assumed to depend only on the magnitude of the
c     wavevector.  N.B. power is dimensionless.
c     ---------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter (n0=128)
      parameter (n0n0=n0*n0,ngr0=n0n0*n0/2,ngdim0=ngr0+n0n0)

      external          primor,trans
      double precision  power

      common            /spectr/rth,rgs,primn,amplit,bias,pindex,ispec
      common            /power/power(ngdim0)
      common            /n1/n1,n12,n21,n1n1,ngrid,n11,ngr2,ngdim2,ngdim


	n12=n1/2
	d3k=dk*dk*dk
	akmax=float(n1)*dk

        do 30 k3=1,n1
          ak3=(k3-1)*dk
          if (k3.gt.n12) ak3=ak3-akmax
          ak33=ak3*ak3
          do 20 k2=1,n1
            ak2=(k2-1)*dk
            if (k2.gt.n12) ak2=ak2-akmax
            ak23=ak2*ak2+ak33
            k23=(k2-1+(k3-1)*n1)*n12
            do 10 k1=2,n12
c Do k1=1 and k1=n12+1 separately below.
              ak1=(k1-1)*dk
              k=k1+k23
              akk=ak1*ak1+ak23
              ak=sqrt(akk)
              power(k)=amplit*primor(ak)*(trans(ak)**2)*d3k
10          continue
            k23=k2-1+(k3-1)*n1
            k23c=n1+1-k2+(n1+1-k3)*n1
            if (k2.eq.1) k23c=k23c-n1
            if (k3.eq.1) k23c=k23c-n1n1
            if (k23.gt.k23c) go to 20
            do 15 i=1,2
c Do k1=1 and k1-n12+1
              if (i.eq.1) then
                k=1+n12*k23
                kc=1+n12*k23c
                ak1=0.0
              else
                k=ngr2+1+k23
                kc=ngr2+1+k23c
                ak1=-0.5*akmax
              endif
              akk=ak1*ak1+ak23
              ak=sqrt(akk)
              if (k.eq.1) then
                power(k)=0.0d0
              else
                power(k)=amplit*primor(ak)*(trans(ak)**2)*d3k
              endif
              power(kc)=power(k)
15          continue             
20        continue
30      continue

       return
      end


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE NORMAL(ISPEC,AK1,SIGMA0)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      EXTERNAL   MIDPNT,MIDINF,INTGT1,INTGT2,INTGT3,INTGT4,
     *           INTGT5,INTGT6

      PARAMETER  (PI=3.141592653,TWOPI=2.*PI)

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC



        IF (ISPEC.EQ.1) THEN
          CALL QROMO(INTGT1,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT1,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (ISPEC.EQ.2) THEN
          CALL QROMO(INTGT2,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT2,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (ISPEC.EQ.3) THEN
          CALL QROMO(INTGT3,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT3,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (ISPEC.EQ.4) THEN
          CALL QROMO(INTGT4,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT4,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (ISPEC.EQ.5) THEN
          CALL QROMO(INTGT5,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT5,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (ISPEC.EQ.6) THEN
          CALL QROMO(INTGT6,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGT6,AK1,1.0D30,S2,MIDINF)
        ENDIF

        AMPLIT=SIGMA0**2/((BIAS**2)*(S1+S2))

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT1(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS1


        INTGT1=PRIMOR(XK)*WINDTH(XK)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT2(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS2


        INTGT2=PRIMOR(XK)*WINDTH(XK)*(TRANS2(XK)**2)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT3(XK)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS3


        INTGT3=PRIMOR(XK)*WINDTH(XK)*(TRANS3(XK)**2)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT4(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS4


        INTGT4=PRIMOR(XK)*WINDTH(XK)*(TRANS4(XK)**2)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT5(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS5


        INTGT5=PRIMOR(XK)*WINDTH(XK)*(TRANS5(XK)**2)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGT6(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDTH,VOLUME,TRANS6


        INTGT6=PRIMOR(XK)*WINDTH(XK)*(TRANS6(XK)**2)*VOLUME(XK)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION PGSSIG(RGAUSS)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      EXTERNAL   INTGG1,INTGG2,INTGG3,INTGG4,
     *           INTGG5,INTGG6,MIDPNT,MIDINF

      PARAMETER  (PI=3.141592653,TWOPI=2.*PI)

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        RGS=RGAUSS
        AK1=TWOPI/RGS

        IF (KSPEC.EQ.1) THEN
          CALL QROMO(INTGG1,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG1,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (KSPEC.EQ.2) THEN
          CALL QROMO(INTGG2,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG2,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (KSPEC.EQ.3) THEN
          CALL QROMO(INTGG3,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG3,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (KSPEC.EQ.4) THEN
          CALL QROMO(INTGG4,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG4,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (KSPEC.EQ.5) THEN
          CALL QROMO(INTGG5,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG5,AK1,1.0D30,S2,MIDINF)
        ENDIF
        IF (KSPEC.EQ.6) THEN
          CALL QROMO(INTGG6,0.0D0,AK1,S1,MIDPNT)
          CALL QROMO(INTGG6,AK1,1.0D30,S2,MIDINF)
        ENDIF

        PGSSIG=SQRT(AMPLIT*(S1+S2))

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG1(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS1

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG1=PRIMOR(XK)*WINDGS(XK)*VOLUME(XK)
        INTGG1=INTGG1*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG2(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS2

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG2=PRIMOR(XK)*WINDGS(XK)*(TRANS2(XK)**2)*VOLUME(XK)
        INTGG2=INTGG2*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG3(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS3

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG3=PRIMOR(XK)*WINDGS(XK)*(TRANS3(XK)**2)*VOLUME(XK)
        INTGG3=INTGG3*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG4(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS4

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG4=PRIMOR(XK)*WINDGS(XK)*(TRANS4(XK)**2)*VOLUME(XK)
        INTGG4=INTGG4*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG5(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS5

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG5=PRIMOR(XK)*WINDGS(XK)*(TRANS5(XK)**2)*VOLUME(XK)
        INTGG5=INTGG5*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION INTGG6(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL PRIMOR,WINDGS,VOLUME,TRANS6

      COMMON    /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,KSPEC


        INTGG6=PRIMOR(XK)*WINDGS(XK)*(TRANS6(XK)**2)*VOLUME(XK)
        INTGG6=INTGG6*(XK**PINDEX)
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION WINDTH(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC

        XKW=XK*RTH
C       ______________
C       Tophat window:
C       --------------
        WINDTH=9.0D0*(SIN(XKW)-XKW*COS(XKW))**2/XKW**6

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION WINDGS(XK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC

        XKW=XK*RGS
C       ________________
C       Gaussian window:
C       ----------------
        WINDGS=EXP(-XKW**2)
        
       RETURN
      END
    

C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION PRIMOR(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC

        PRIMOR=XK**PRIMN
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS1(XK)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)


        TRANS1=1.0D0
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS2(XK)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC
      COMMON  /COSMOS/H0,HUBBL0,OMEGA0,OMEGAHH,HF0


C       ___________________________________________________
C       This is the Cold Dark Matter spectrum for adiabatic 
C       fluctuations as given by:
C          Davis, Efstathiou, Frenk and White, 
C          Astrophys. J. 292, 371 (1985).
C       ---------------------------------------------------
        Q=XK/(OMEGA0*H0*H0)
        T1=1.7D0*Q
        T2=9.0D0*Q*SQRT(Q)
        T3=1.0D0*Q*Q
        T4=1.0D0+T1+T2+T3
        TRANS2=1.0D0/T4
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS3(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC
      COMMON  /COSMOS/H0,HUBBL0,OMEGA0,OMEGAHH,HF0


C       ___________________________________________________
C       This is the Cold Dark Matter spectrum for adiabatic
C       fluctuations as given by:
C         Bardeen, Bond, Kaiser and Szalay,
C         Astrophys. J. 304, 15 (1986)
C       ---------------------------------------------------
        Q=XK/(OMEGA0*H0*H0)
        T1=2.34D0*Q
        T2=LOG(1.0D0+T1)/T1
        T3=3.89D0*Q
        T4=(16.1D0*Q)**2
        T5=(5.46D0*Q)**3
        T6=(6.71D0*Q)**4
        T7=(1.0D0+T3+T4+T5+T6)**(1./4.)
        TRANS3=T2/T7

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS4(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       
      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC
      COMMON  /COSMOS/H0,HUBBL0,OMEGA0,OMEGAHH,HF0


C       ________________________________________________
C       This is the spectrum for massive neutrinos with 
C       adiabatic fluctuations as given by:
C         Bardeen, Bond, Kaiser and Szalay,
C         Astrophys. J. 304, 15 (1986)
C       ------------------------------------------------
        Q=XK/omegahh
        RFNU=2.6D0/omegahh
        T1=0.16D0*RFNU*XK
        T2=0.5D0*(RFNU*XK)**2
        T3=EXP(-T1-T2)
        T4=1.6D0*Q
        T5=(4.0D0*Q)**(1.5)
        T6=(0.92D0*Q)**2
        T7=1.0D0+T4+T5+T6
        TRANS4=T3/T7

       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS5(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       
      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC
      COMMON  /COSMOS/H0,HUBBL0,OMEGA0,OMEGAHH,HF0


C       _____________________________________________________
C       This is the spectrum for Warm Dark Matter, adiabatic
C       fluctuations as given by:
C          Bardeen, Bond, Kaiser and Szalay,
C          Astrophys. J. 304, 15 (1986)
C       for a species with GXDEC effective number degrees of
C       freedom when the X particles decoupled; values in the 
C       range 60-300 are typical of minimal GUTs over the 
C       range of decoupling temperature T = 1-10**18 GeV.
C       -----------------------------------------------------

        GXDEC=100.0D0

        Q=XK/(OMEGA0*H0*H0)
        RFW=0.2D0*((100./GXDEC)**(4./3.))/(OMEGA0*H0*H0)
        T1=0.5D0*XK*RFW
        T2=0.5D0*((XK*RFW)**2)
        T3=EXP(-T1-T2)
        T4=1.7D0*Q
        T5=(4.3D0*Q)**(1.5)
        T6=Q*Q
        T7=1.0D0+T4+T5+T6
        TRANS5=T3/T7
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION TRANS6(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       
      COMMON  /SPECTR/RTH,RGS,PRIMN,AMPLIT,BIAS,PINDEX,ISPEC
      COMMON  /COSMOS/H0,HUBBL0,OMEGA0,OMEGAHH,HF0


C       _______________________________________________________
C       This is the spectrum for Cold Dark Matter, isocurvature
C       fluctuations as given by:
C          Bardeen, Bond, Kaiser and Szalay,
C          Astrophys. J. 304, 15 (1986)
C       -------------------------------------------------------

        Q=XK/(OMEGA0*H0*H0)
        T1=(5.6D0*Q)**2
        T2=(40.0D0*Q)**2
        T3=215.0D0*Q
        T4=(16.0D0*Q)**2
        T5=1.0D0+0.5D0*Q
        T6=T2/(1.0D0+T3+T4/T5)
        T7=(5.6D0*Q)**(8./5.)
        T8=(1.0D0+T6+T7)**(5./4.)
        TRANS6=T1/T8
        
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      DOUBLE PRECISION FUNCTION VOLUME(XK)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI=3.141592653)

        VOLUME=4.0D0*PI*(XK**2)
       RETURN
      END


c     *************************************************************************
c     *************************************************************************


	subroutine fft3r(a,n1)

c  fft3r performs a 3-D FFT from the spatial domain to the spectral
c  domain.   Array a is real of length n1*n1*(n1+2).
c  n1 must be divisible by 4 but it need not be a power of 2.
c  On input, the first n1*n1*n1 elements are filled with values of a(i,j,k)
c  in the spatial domain.  On output, the first n1*n1*n1 elements are filled
c  with values a(i,j,k), i=1,...,n1/2 of the complex Fourier transform.
c  The last 2*n1*n1 values are filled with complex transform values for
c  i=n/2+1.
c  N.B. The transform from spatial to spectral domain is defined with a
c  plus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n1*n1.
c

	parameter (twopi=-6.283185307179586d0,ndim=3,n0=512)
	implicit real (a-h,o-z)
	real a(*),wra(n0),wia(n0),work(2*n0)
	integer nn(ndim)
	double precision theta,wr,wi
	equivalence (wra(1),work(1)),(wia(1),work(n0+1))
c
	if (n1.gt.n0) then
	  write(*,*) 'Need to increase dimension n0 in fft3rinv!'
	  write(*,*) 'n0,n1=',n0,n1
	  write(1,*) 'Need to increase dimension n0 in fft3rinv!'
	  write(1,*) 'n0,n1=',n0,n1
	  stop
	end if
c
c  fourt requires contiguous storage, i.e., no strides allowed.
	n1d2=n1/2
	inc2y=n1/2
	inc3y=inc2y*n1
	n3eff=1+(n1-1)+2*((n1-1)*inc2y+(n1-1)*inc3y)
c
	nn(1)=n1d2
	nn(2)=n1
	nn(3)=n1
	call fourt(a,nn,ndim,-1,1,work)
c
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
	theta=twopi/n1
	  do 24 k1=1,n1
	  wra(k1)=cos(theta*(k1-1))
	  wia(k1)=sin(theta*(k1-1))
24	continue
c  Save k1=1,n14+1 for later.
	n14=n1d2/2
	  do 30 k3=1,n1
	  k3c=n1+2-k3
	  if (k3.eq.1) k3c=1
c  Swap k1 and k2 loops so that the inner loop is longer (for vectorization).
	    do 28 k1=2,n14
	    k1c=n1d2+2-k1
	    wr=wra(k1)
	    wi=wia(k1)
c*vdir prefer vector
	      do 26 k2=1,n1
	      k2c=n1+2-k2
	      if (k2.eq.1) k2c=1
	      ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	      indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	      ffer=0.5*(a(ind)+a(indc))
	      ffei=0.5*(a(ind+1)-a(indc+1))
	      ffor=0.5*(a(ind+1)+a(indc+1))
	      ffoi=-0.5*(a(ind)-a(indc))
	      tempr=wr*ffor-wi*ffoi
	      tempi=wi*ffor+wr*ffoi
	      a(ind)=ffer+tempr
	      a(ind+1)=ffei+tempi
	      a(indc)=ffer-tempr
	      a(indc+1)=-ffei+tempi
26	    continue
28	  continue
30	continue
	  do 34 k3=1,n1
	  k3c=n1+2-k3
	  if (k3.eq.1) k3c=1
	    do 32 k2=1,n1
	    k2c=n1+2-k2
	    if (k2.eq.1) k2c=1
c  k1=n14.
	    k1=n14+1
	    k1c=n1d2+2-k1
	    ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 31
	    ffer=0.5*(a(ind)+a(indc))
	    ffei=0.5*(a(ind+1)-a(indc+1))
	    ffor=0.5*(a(ind+1)+a(indc+1))
	    ffoi=-0.5*(a(ind)-a(indc))
	    wr=wra(k1)
	    wi=wia(k1)
	    tempr=wr*ffor-wi*ffoi
	    tempi=wi*ffor+wr*ffoi
	    a(ind)=ffer+tempr
	    a(ind+1)=ffei+tempi
	    a(indc)=ffer-tempr
	    a(indc+1)=-ffei+tempi
c  k1=1.
31	    k1=1
	    ind=1+2*((k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 32
	    k23=(k2-1)+(k3-1)*n1
	    k23c=(k2c-1)+(k3c-1)*n1
	    i3=n3eff+1+2*k23
	    i4=n3eff+1+2*k23c
	    ffer=0.5*(a(ind)+a(indc))
	    ffei=0.5*(a(ind+1)-a(indc+1))
	    ffor=0.5*(a(ind+1)+a(indc+1))
	    ffoi=-0.5*(a(ind)-a(indc))
	    a(ind)=ffer+ffor
	    a(ind+1)=ffei+ffoi
	    a(indc)=a(ind)
	    a(indc+1)=-a(ind+1)
	    a(i4)=ffer-ffor
	    a(i4+1)=-ffei+ffoi
	    a(i3)=a(i4)
	    a(i3+1)=-a(i4+1)
32	 continue
34	continue
	return
	end


c     *************************************************************************
c     *************************************************************************


	subroutine fft3rinv(a,n1)

c  fft3rinv performs a 3-D inverse FFT from the spectral domain to the
c  spatial domain.  Array a is real of length n1*n1*(n1+2).
c  n1 must be divisible by 4 but it need not be a power of 2.
c  On input, the first n1*n1*n1 elements are filled with values a(i,j,k),
c  i=1,...,n1/2 of the complex Fourier transform.  The last 2*n1*n1 values
c  are filled with complex transform values for i=n1/2+1.  On output, the
c  first n1*n1*n1 elements are filled with values of a(i,j,k) in the spatial
c  domain.
c  N.B. The transform from spectral to spatial domain is defined with a
c  minus sign in the complex exponential.  To change this, change the sign
c  of twopi.  After transforming and then performing the inverse transform,
c  a(i,j,k) must be divided by n1*n1*n1.
c

	parameter (twopi=6.283185307179586d0,ndim=3,n0=512)
	implicit real (a-h,o-z)
	real a(*),wra(n0),wia(n0),work(2*n0)
	integer nn(ndim)
	double precision theta,wr,wi
	equivalence (wra(1),work(1)),(wia(1),work(n0+1))
c
	if (n1.gt.n0) then
	  write(*,*) 'Need to increase dimension n0 in fft3r!'
	  write(*,*) 'n0,n1=',n0,n1
	  write(1,*) 'Need to increase dimension n0 in fft3r!'
	  write(1,*) 'n0,n1=',n0,n1
	  stop
	end if
c
c  fourt requires contiguous storage, i.e., no strides allowed.
	n1d2=n1/2
	inc2y=n1/2
	inc3y=inc2y*n1
	n3eff=1+(n1-1)+2*((n1-1)*inc2y+(n1-1)*inc3y)
c
c  Use symmetries to obtain the true transform.  This part is specially
c  written for 3-D.
	theta=twopi/n1
	  do 24 k1=1,n1
	  wra(k1)=cos(theta*(k1-1))
	  wia(k1)=sin(theta*(k1-1))
24	continue
c  Save k1=1,n14+1 for later.
	n14=n1d2/2
	  do 30 k3=1,n1
	  k3c=n1+2-k3
	  if (k3.eq.1) k3c=1
c  Swap k1 and k2 loops so that the inner loop is longer (for vectorization).
	    do 28 k1=2,n14
	    k1c=n1d2+2-k1
	    wr=wra(k1)
	    wi=wia(k1)
c*vdir prefer vector
	      do 26 k2=1,n1
	      k2c=n1+2-k2
	      if (k2.eq.1) k2c=1
	      ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	      indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
c  Multiply by 2 since this is a real FFT summing over only half
c  of the harmonics.
	      ffer=a(ind)+a(indc)
	      ffei=a(ind+1)-a(indc+1)
	      ffor=-(a(ind+1)+a(indc+1))
	      ffoi=a(ind)-a(indc)
	      tempr=wr*ffor-wi*ffoi
	      tempi=wi*ffor+wr*ffoi
	      a(ind)=ffer+tempr
	      a(ind+1)=ffei+tempi
	      a(indc)=ffer-tempr
	      a(indc+1)=-ffei+tempi
26	    continue
28	  continue
30	continue
	  do 34 k3=1,n1
	  k3c=n1+2-k3
	  if (k3.eq.1) k3c=1
	    do 32 k2=1,n1
	    k2c=n1+2-k2
	    if (k2.eq.1) k2c=1
c  k1=n14.
	    k1=n14+1
	    k1c=n1d2+2-k1
	    ind=1+2*((k1-1)+(k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k1c-1)+(k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 31
	    ffer=a(ind)+a(indc)
	    ffei=a(ind+1)-a(indc+1)
	    ffor=-(a(ind+1)+a(indc+1))
	    ffoi=a(ind)-a(indc)
	    wr=wra(k1)
	    wi=wia(k1)
	    tempr=wr*ffor-wi*ffoi
	    tempi=wi*ffor+wr*ffoi
	    a(ind)=ffer+tempr
	    a(ind+1)=ffei+tempi
	    a(indc)=ffer-tempr
	    a(indc+1)=-ffei+tempi
c  k1=1.
31	    k1=1
	    ind=1+2*((k2-1)*inc2y+(k3-1)*inc3y)
	    indc=1+2*((k2c-1)*inc2y+(k3c-1)*inc3y)
	    if (ind.gt.indc) go to 32
	    k23=(k2-1)+(k3-1)*n1
	    k23c=(k2c-1)+(k3c-1)*n1
	    i4=n3eff+1+2*k23c
	    ffer=a(ind)+a(i4)
	    ffei=a(ind+1)-a(i4+1)
	    ffor=-(a(ind+1)+a(i4+1))
	    ffoi=a(ind)-a(i4)
	    a(ind)=ffer+ffor
	    a(ind+1)=ffei+ffoi
	    a(indc)=ffer-ffor
	    a(indc+1)=-ffei+ffoi
32	 continue
34	continue
c
	nn(1)=n1d2
	nn(2)=n1
	nn(3)=n1
	call fourt(a,nn,ndim,+1,1,work)
c
	return
	end


c     *************************************************************************
c     *************************************************************************


	SUBROUTINE FOURT (DATA,NN,NDIM,ISIGN,IFORM,WORK)

C
C**********************************************************************
C*								      *
C*		THE COOLEY-TUKEY FAST FOURIER TRANSFORM		      *
C*								      *
C*   TRANSFORM(K1,K2,...) = SUM(DATA(J1,J2,...)*EXP(ISIGN*2*SQRT(-1)  *
C*   *((J1-1)*(K1-1)/NN(1)+(J2-1)*(K2-1)/NN(2)+...))),SUMMED FOR ALL  *
C*   J1, K1 BETWEEN 1 AND NN(1), J2, K2 BETWEEN 1 AND NN(2), ETC.     *
C*   THERE IS NO LIMIT TO THE NUMBER OF SUBSCRIPTS.  DATA IS A        *
C*   MULTIDIMENSIONAL COMPLEX ARRAY WHOSE REAL AND IMAGINARY	      *
C*   PARTS ARE ADJACENT IN STORAGE, SUCH AS FORTRAN IV PLACES THEM.   *
C*   IF ALL IMAGINARY PARTS ARE ZERO (DATA ARE DISGUISED REAL), SET   *
C*   IFORM TO ZERO TO CUT THE RUNNING TIME BY UP TO FORTY PERCENT.    *
C*   OTHERWISE, IFORM = +1. THE LENGTHS OF ALL DIMENSIONS ARE	      *
C*   STORED IN ARRAY NN, OF LENGTH NDIM. THEY MAY BE ANY POSITIVE     *
C*   INTEGERS, THO THE PROGRAM RUNS FASTER ON COMPOSITE INTEGERS, AND *
C*   ESPECIALLY FAST ON NUMBERS RICH IN FACTORS OF TWO. ISIGN IS +1   *
C*   OR -1.  IF A -1 TRANSFORM IS FOLLOWED BY A +1 ONE (OR A +1	      *
C*   BY A -1) THE ORIGINAL DATA REAPPEAR, MULTIPLIED BY NTOT (NN(1)   *
C*   *NN(2)*...).TRANSFORM VALUES ARE ALWAYS COMPLEX,AND ARE RETURNED *
C*   IN ARRAY DATA, REPLACING THE INPUT. IN ADDITION, IF ALL	      *
C*   DIMENSIONS ARE NOT POWERS OF TWO, ARRAY WORK MUST BE SUPPLIED,   *
C*   COMPLEX OF LENGTH EQUAL TO THE LARGEST NON 2**K DIMENSION.       *
C*   OTHERWISE, REPLACE WORK BY ZERO IN THE CALLING SEQUENCE.	      *
C*   NORMAL FORTRAN DATA ORDERING IS EXPECTED, FIRST SUBSCRIPT	      *
C*   VARYING FASTEST. ALL SUBSCRIPTS BEGIN AT ONE.		      *
C*								      *
C*   RUNNING TIME IS MUCH SHORTER THAN THE NAIVE NTOT**2, BEING       *
C*   GIVEN BY THE FOLLOWING FORMULA.  DECOMPOSE NTOT INTO	      *
C*   2**K2 * 3**K3 * 5**K5 * ....LET SUM2=2*K2,SUMF=3*K3+5*K5+...     *
C*   AND NF=K3+K5+...  THE TIME TAKEN BY A MULTIDIMENSIONAL	      *
C*   TRANSFORM ON THESE NTOT DATA IS T = T0 + NTOT*(T1+T2*SUM2+       *
C*   +T3*SUMF+T4*NF). ON THE CDC 3300 (FLOATING POINT ADD TIME OF     *
C*   SIX MICROSECONDS), T = 3000 + NTOT * (500+43*SUM2+68*SUMF+       *
C*   +320*NF) MICROSECONDS ON COMPLEX DATA. IN ADDITION, THE	      *
C*   ACCURACY IS GREATLY IMPROVED, AS THE RMS RELATIVE ERROR IS       *
C*   BOUNDED BY 3*2**(-B)*SUM(FACTOR(J)**1.5), WHERE B IS THE NUMBER  *
C*   OF BITS IN THE FLOATING POINT FRACTION AND FACTOR(J) ARE THE     *
C*   PRIME FACTORS OF NTOT.					      *
C*								      *
C*   THE DISCRETE FOURIER TRANSFORM PLACES THREE RESTRICTIONS UPON    *
C*   THE DATA:							      *
C*   1.  THE NUMBER OF INPUT DATA AND THE NUMBER OF TRANSFORM VALUES  *
C*	 MUST BE THE SAME.					      *
C*   2.  BOTH THE INPUT DATA AND THE TRANSFORM VALUES MUST REPRESENT  *
C*	 EQUISPACED POINTS IN THEIR RESPECTIVE DOMAINS OF TIME AND    *
C*	 FREQUENCY.CALLING THESE SPACINGS DELTAT AND DELTAF, IT MUST  *
C*	 BE TRUE THAT DELTAF=2*PI/(NN(I)*DELTAT).OF COURSE, DELTAT    *
C*	 NEED NOT BE THE SAME FOR EVERY DIMENSION.		      *
C*   3.  CONCEPTUALLY AT LEAST, THE INPUT DATA AND THE TRANSFORM      *
C*	 OUTPUT REPRESENT SINGLE CYCLES OF PERIODIC FUNCTIONS.	      *
C*								      *
C*	DC-COMPONENT IS MULTIPLIED BY NN, OTHERS BY NN/2	      *
C*								      *
C*   EXAMPLE 1. THREE-DIMENSIONAL FORWARD FFT OF A COMPLEX ARRAY      *
C*		DIMENSIONED 32 BY 25 BY 13 IN FORTRAN IV.	      *
C*		DIMENSION DATA(32,25,13),WORK(50),NN(3)		      *
C*		COMPLEX   DATA					      *
C*		DATA	  NN	/32,25,13/			      *
C*		DO 1 I = 1,32					      *
C*		DO 1 J = 1,25					      *
C*		DO 1 K = 1,13					      *
C*   1		DATA(I,J,K) = COMPLEX VALUE			      *
C*		CALL FOURT (DATA,NN,3,-1,1,WORK)		      *
C*								      *
C*   EXAMPLE 2. ONE-DIMENSIONAL FORWARD FFT OF A REAL ARRAY OF	      *
C*		LENGTH 64 IN FORTRAN II.			      *
C*		DIMENSION DATA(2,64)				      *
C*		DO 2 I = 1,64					      *
C*		DATA(1,I) = REAL PART				      *
C*   2		DATA(2,I) = 0.					      *
C*		CALL FOURT (DATA,64,1,-1,0,0)			      *
C*								      *
C**********************************************************************
C
	DIMENSION DATA(*),NN(*),IFACT(32),WORK(*)
	TWOPI = 6.2831853
	IF (NDIM-1) 920,1,1
    1	NTOT = 2
	DO 2 IDIM = 1,NDIM
	  IF (NN(IDIM)) 920,920,2
    2	NTOT = NTOT * NN(IDIM)
C
C		*** MAIN LOOP FOR EACH DIMENSION ***
C
	NP1 = 2
	DO 910 IDIM = 1,NDIM
	  N = NN(IDIM)
	  NP2 = NP1 * N
	  IF (N-1) 920,900,5
C
C		*** FACTOR N ***
C
    5	  M = N
	  NTWO = NP1
	  IF = 1
	  IDIV = 2
   10	  IQUOT = M / IDIV
	  IREM = M - IDIV * IQUOT
	  IF (IQUOT-IDIV) 50,11,11
   11	  IF (IREM) 20,12,20
   12	  NTWO = NTWO + NTWO
	  M = IQUOT
	  GO TO 10
   20	  IDIV = 3
   30	  IQUOT = M / IDIV
	  IREM = M - IDIV * IQUOT
	  IF (IQUOT-IDIV) 60,31,31
   31	  IF (IREM) 40,32,40
   32	  IFACT(IF) = IDIV
	  IF = IF + 1
	  M = IQUOT
	  GO TO 30
   40	  IDIV = IDIV + 2
	  GO TO 30
   50	  IF (IREM) 60,51,60
   51	  NTWO = NTWO + NTWO
	  GO TO 70
   60	  IFACT(IF) = M
C
C	SEPARATE FOUR CASES--
C	  1. COMPLEX TRANSFORM OR REAL TRANSFORM FOR THE 4TH,5TH,ETC.
C	     DIMENSIONS.
C	  2. REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION.  METHOD--
C	     TRANSFORM HALF THE DATA, SUPPLYING THE OTHER HALF BY
C	     CONJUGATE SYMMETRY.
C	  3. REAL TRANSFORM FOR THE 1ST DIMENSION, N ODD.  METHOD--
C	     TRANSFORM HALF THE DATA AT EACH STAGE, SUPPLYING THE
C	     OTHER HALF BY CONJUGATE SYMMETRY.
C	  4. REAL TRANSFORM FOR THE 1ST DIMENSION,N EVEN.  METHOD--
C	     TRANSFORM A COMPLEX ARRAY OF LENGTH N/2 WHOSE REAL PARTS
C	     ARE THE EVEN-NUMBERED REAL VALUES AND WHOSE IMAGINARY
C	     PARTS ARE THE ODD-NUMBERED REAL VALUES. SEPARATE AND
C	     SUPPLY THE SECOND HALF BY CONJUGATE SYMMETRY.
C
   70	  NON2 = NP1 * (NP2 / NTWO)
	  ICASE = 1
	  IF (IDIM-4) 71,90,90
   71	  IF (IFORM) 72,72,90
   72	  ICASE = 2
	  IF (IDIM-1) 73,73,90
   73	  ICASE = 3
	  IF (NTWO-NP1) 90,90,74
   74	  ICASE = 4
	  NTWO = NTWO / 2
	  N = N / 2
	  NP2 = NP2 / 2
	  NTOT = NTOT / 2
	  I = 3
	  DO 80 J = 2,NTOT
	    DATA(J) = DATA(I)
	    I = I + 2
   80	  CONTINUE
   90	  I1RNG = NP1
	  IF (ICASE-2) 100,95,100
   95	  I1RNG = NP0 * (1 + NPREV / 2)
C
C	*** SHUFFLE ON THE FACTORS OF TWO IN N. AS THE SHUFFLING
C	*** CAN BE DONE BY SIMPLE INTERCHANGE, NO WORK AREA IS NEEDED.
C
  100	  IF (NTWO-NP1) 600,600,110
  110	  NP2HF = NP2 / 2
	  J = 1
	  DO 150 I2 = 1,NP2,NON2
	    IF (J-I2) 120,130,130
  120	    I1MAX = I2 + NON2 - 2
	    DO 125 I1 = I2,I1MAX,2
	      DO 125 I3 = I1,NTOT,NP2
	        J3 = J + I3 - I2
	        TEMPR = DATA(I3)
		TEMPI = DATA(I3+1)
		DATA(I3) = DATA(J3)
		DATA(I3+1) = DATA(J3+1)
		DATA(J3) = TEMPR
		DATA(J3+1) = TEMPI
  125	    CONTINUE
  130	    M = NP2HF
  140	    IF (J-M) 150,150,145
  145	    J = J - M
	    M = M / 2
	    IF (M-NON2) 150,140,140
  150	  J = J + M
C
C	MAIN LOOP FOR FACTORS OF TWO. PERFORM FOURIER TRANSFORMS OF
C	LENGTH FOUR, WITH ONE OF LENGTH TWO IF NEEDED. THE
C	TWIDDLE FACTOR W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)). CHECK
C	FOR W=ISIGN*SQRT(-1) AND REPEAT FOR W=SQRT(-1)*CONJUGATE(W).
C
	  NON2T = NON2 + NON2
	  IPAR = NTWO / NP1
  310	  IF (IPAR-2) 350,330,320
  320	  IPAR = IPAR / 4
	  GO TO 310
  330	  CONTINUE
	  DO 340 I1 = 1,I1RNG,2
	    DO 340 J3 = I1,NON2,NP1
	      DO 340 K1 = J3,NTOT,NON2T
		K2 = K1 + NON2
		TEMPR = DATA(K2)
		TEMPI = DATA(K2+1)
		DATA(K2) = DATA(K1) - TEMPR
		DATA(K2+1) = DATA(K1+1) - TEMPI
		DATA(K1) = DATA(K1) + TEMPR
		DATA(K1+1) = DATA(K1+1) + TEMPI
  340	  CONTINUE
  350	  MMAX = NON2
  360	  IF (MMAX-NP2HF) 370,600,600
  370	  LMAX = MAX0 (NON2T,MMAX/2)
	  IF (MMAX-NON2) 405,405,380
  380	  THETA = -TWOPI * FLOAT(NON2)/FLOAT(4*MMAX)
	  IF (ISIGN) 400,390,390
  390	  THETA = -THETA
  400	  WR = COS (THETA)
	  WI = SIN (THETA)
	  WSTPR = -2. * WI * WI
	  WSTPI =  2. * WR * WI
  405	  CONTINUE
	  DO 570 L = NON2,LMAX,NON2T
	    M = L
	    IF (MMAX-NON2) 420,420,410
  410	    W2R = WR * WR - WI * WI
	    W2I = 2. * WR * WI
	    W3R = W2R * WR - W2I * WI
	    W3I = W2R * WI + W2I * WR
  420	   CONTINUE
	    DO 530 I1 = 1,I1RNG,2
	      DO 530 J3 = I1,NON2,NP1
		KMIN = J3 + IPAR * M
		IF (MMAX-NON2) 430,430,440
  430		KMIN = J3
  440		KDIF = IPAR * MMAX
  450		KSTEP = 4 * KDIF
		DO 520 K1 = KMIN,NTOT,KSTEP
		  K2 = K1 + KDIF
		  K3 = K2 + KDIF
		  K4 = K3 + KDIF
		  IF (MMAX-NON2) 460,460,480
  460		  U1R = DATA(K1) + DATA(K2)
		  U1I = DATA(K1+1) + DATA(K2+1)
		  U2R = DATA(K3) + DATA(K4)
		  U2I = DATA(K3+1) + DATA(K4+1)
		  U3R = DATA(K1) - DATA(K2)
		  U3I = DATA(K1+1) - DATA(K2+1)
		  IF (ISIGN) 470,475,475
  470		  U4R = DATA(K3+1) -DATA(K4+1)
		  U4I = DATA(K4) - DATA(K3)
		  GO TO 510
  475		  U4R = DATA(K4+1) - DATA(K3+1)
		  U4I = DATA(K3) - DATA(K4)
		  GO TO 510
  480		  T2R = W2R * DATA(K2) - W2I * DATA(K2+1)
		  T2I = W2R * DATA(K2+1) + W2I * DATA(K2)
		  T3R = WR * DATA(K3) - WI * DATA(K3+1)
		  T3I = WR * DATA(K3+1) + WI * DATA(K3)
		  T4R = W3R * DATA(K4) - W3I * DATA(K4+1)
		  T4I = W3R * DATA(K4+1) + W3I * DATA(K4)
		  U1R = DATA(K1) + T2R
		  U1I = DATA(K1+1) + T2I
		  U2R = T3R + T4R
		  U2I = T3I + T4I
		  U3R = DATA(K1) - T2R
		  U3I = DATA(K1+1) - T2I
		  IF (ISIGN) 490,500,500
  490		  U4R = T3I - T4I
		  U4I = T4R - T3R
		  GO TO 510
  500		  U4R = T4I - T3I
		  U4I = T3R - T4R
  510		  DATA(K1) = U1R + U2R
		  DATA(K1+1) = U1I + U2I
		  DATA(K2) = U3R + U4R
		  DATA(K2+1) = U3I + U4I
		  DATA(K3) = U1R - U2R
		  DATA(K3+1) = U1I - U2I
		  DATA(K4) = U3R - U4R
		  DATA(K4+1) = U3I - U4I
  520		CONTINUE
		KMIN = 4 * (KMIN-J3) + J3
		KDIF = KSTEP
		IF (KDIF-NP2) 450,530,530
  530	      CONTINUE
	      M = MMAX - M
	      IF (ISIGN) 540,550,550
  540	      TEMPR = WR
	      WR = -WI
	      WI = -TEMPR
	      GO TO 560
  550	      TEMPR = WR
	      WR = WI
	      WI = TEMPR
  560	      IF (M-LMAX) 565,565,410
  565	      TEMPR = WR
	      WR = WR * WSTPR - WI * WSTPI + WR
	      WI = WI * WSTPR + TEMPR * WSTPI + WI
  570	  CONTINUE
	  IPAR = 3 - IPAR
	  MMAX = MMAX + MMAX
	  GO TO 360
C
C	MAIN LOOP FOR FACTORS NOT EQUAL TO TWO. APPLY THE TWIDDLE
C	FACTOR W=EXP(ISIGN*2*PI*SQRT(-1)*(J2-1)*(J1-J2)/(NP2*IFP1)),
C	THEN PERFORM A FOURIER TRANSFORM OF LENGTH IFACT(IF),MAKING
C	USE OF CONJUGATE SYMMETRIES.
C
  600	  IF (NTWO-NP2) 605,700,700
  605	  IFP1 = NON2
	  IF = 1
	  NP1HF = NP1 / 2
  610	  IFP2 = IFP1 / IFACT(IF)
	  J1RNG = NP2
	  IF (ICASE-3) 612,611,612
  611	  J1RNG = (NP2+IFP1) / 2
	  J2STP = NP2 / IFACT(IF)
	  J1RG2 = (J2STP+IFP2) / 2
  612	  J2MIN = 1 + IFP2
	  IF (IFP1-NP2) 615,640,640
  615	  CONTINUE
	  DO 635 J2 = J2MIN,IFP1,IFP2
	    THETA = -TWOPI * FLOAT(J2-1) / FLOAT(NP2)
	    IF (ISIGN) 625,620,620
  620	    THETA = -THETA
  625	    SINTH = SIN (THETA/2.)
	    WSTPR = -2. * SINTH * SINTH
	    WSTPI = SIN (THETA)
	    WR = WSTPR + 1.
	    WI = WSTPI
	    J1MIN = J2 + IFP1
	    DO 635 J1 = J1MIN,J1RNG,IFP1
	      I1MAX = J1 + I1RNG - 2
	      DO 630 I1 = J1,I1MAX,2
		DO 630 I3 = I1,NTOT,NP2
		  J3MAX = I3 + IFP2 - NP1
		  DO 630 J3 = I3,J3MAX,NP1
		    TEMPR = DATA(J3)
		    DATA(J3) = DATA(J3) * WR - DATA(J3+1) * WI
		    DATA(J3+1) = TEMPR * WI + DATA(J3+1) * WR
  630	      CONTINUE
	      TEMPR = WR
	      WR = WR * WSTPR - WI * WSTPI + WR
	      WI = TEMPR * WSTPI + WI * WSTPR + WI
  635	  CONTINUE
  640	  THETA = -TWOPI / FLOAT(IFACT(IF))
	  IF (ISIGN) 650,645,645
  645	  THETA = -THETA
  650	  SINTH = SIN (THETA/2.)
	  WSTPR = -2. * SINTH * SINTH
	  WSTPI = SIN (THETA)
	  KSTEP = 2 * N / IFACT(IF)
	  KRANG = KSTEP * (IFACT(IF)/2) + 1
	  DO 698 I1 = 1,I1RNG,2
	    DO 698 I3 = I1,NTOT,NP2
	      DO 690 KMIN = 1,KRANG,KSTEP
		J1MAX = I3 + J1RNG - IFP1
		DO 680 J1 = I3,J1MAX,IFP1
		  J3MAX = J1 + IFP2 - NP1
		  DO 680 J3 = J1,J3MAX,NP1
		    J2MAX = J3 + IFP1 - IFP2
		    K = KMIN + (J3-J1+(J1-I3)/IFACT(IF)) / NP1HF
		    IF (KMIN-1) 655,655,665
  655		    SUMR = 0.
		    SUMI = 0.
		    DO 660 J2 = J3,J2MAX,IFP2
		      SUMR = SUMR + DATA(J2)
		      SUMI = SUMI + DATA(J2+1)
  660		    CONTINUE
		    WORK(K) = SUMR
		    WORK(K+1) = SUMI
		    GO TO 680
  665		    KCONJ = K + 2 * (N-KMIN+1)
		    J2 = J2MAX
		    SUMR = DATA(J2)
		    SUMI = DATA(J2+1)
		    OLDSR = 0.
		    OLDSI = 0.
		    J2 = J2 - IFP2
  670		    TEMPR = SUMR
		    TEMPI = SUMI
		    SUMR = TWOWR * SUMR - OLDSR + DATA(J2)
		    SUMI = TWOWR * SUMI - OLDSI + DATA(J2+1)
		    OLDSR = TEMPR
		    OLDSI = TEMPI
		    J2 = J2 - IFP2
		    IF (J2-J3) 675,675,670
  675		    TEMPR = WR * SUMR - OLDSR + DATA(J2)
		    TEMPI = WI * SUMI
		    WORK(K) = TEMPR - TEMPI
		    WORK(KCONJ) = TEMPR + TEMPI
		    TEMPR = WR * SUMI - OLDSI + DATA(J2+1)
		    TEMPI = WI * SUMR
		    WORK(K+1) = TEMPR + TEMPI
		    WORK(KCONJ+1) = TEMPR - TEMPI
  680		CONTINUE
		IF (KMIN-1) 685,685,686
  685		WR = WSTPR + 1.
		WI = WSTPI
		GO TO 690
  686		TEMPR = WR
		WR = WR * WSTPR - WI * WSTPI + WR
		WI = TEMPR * WSTPI + WI * WSTPR + WI
  690	      TWOWR = WR + WR
	      IF (ICASE-3) 692,691,692
  691	      IF (IFP1-NP2) 695,692,692
  692	      K = 1
	      I2MAX = I3 + NP2 - NP1
	      DO 693 I2 = I3,I2MAX,NP1
		DATA(I2) = WORK(K)
		DATA(I2+1) = WORK(K+1)
		K = K + 2
  693	      CONTINUE
	      GO TO 698
C
C	COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N ODD, BY
C	CONJUGATE SYMMETRIES AT EACH STAGE.
C
  695	      J3MAX = I3 + IFP2 - NP1
	      DO 697 J3 = I3,J3MAX,NP1
		J2MAX = J3 + NP2 - J2STP
		DO 697 J2 = J3,J2MAX,J2STP
		  J1MAX = J2 + J1RG2 - IFP2
		  J1CNJ = J3 + J2MAX + J2STP - J2
		  DO 697 J1 = J2,J1MAX,IFP2
		    K = 1 + J1 - I3
		    DATA(J1) = WORK(K)
		    DATA(J1+1) = WORK(K+1)
		    IF (J1-J2) 697,697,696
  696		    DATA(J1CNJ) = WORK(K)
		    DATA(J1CNJ+1) = -WORK(K+1)
  697	      J1CNJ = J1CNJ - IFP2
  698	    CONTINUE
	    IF = IF + 1
	    IFP1 = IFP2
	    IF (IFP1-NP1) 700,700,610
C
C	COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N EVEN, BY
C	CONJUGATE SYMMETRIES.
C
  700	  GO TO (900,800,900,701) , ICASE
  701	  NHALF = N
	  N = N + N
	  THETA = -TWOPI / FLOAT(N)
	  IF (ISIGN) 703,702,702
  702	  THETA = -THETA
  703	  SINTH = SIN (THETA/2.)
	  WSTPR = -2. * SINTH * SINTH
	  WSTPI = SIN (THETA)
	  WR = WSTPR + 1.
	  WI = WSTPI
	  IMIN = 3
	  JMIN = 2 * NHALF - 1
	  GO TO 725
  710	  J = JMIN
	  DO 720 I = IMIN,NTOT,NP2
	    SUMR = (DATA(I)+DATA(J)) / 2.
	    SUMI = (DATA(I+1)+DATA(J+1)) / 2.
	    DIFR = (DATA(I)-DATA(J)) / 2.
	    DIFI = (DATA(I+1)-DATA(J+1)) / 2.
	    TEMPR = WR * SUMI + WI * DIFR
	    TEMPI = WI * SUMI - WR * DIFR
	    DATA(I) = SUMR + TEMPR
	    DATA(I+1) = DIFI + TEMPI
	    DATA(J) = SUMR - TEMPR
	    DATA(J+1) = -DIFI + TEMPI
  720	  J = J + NP2
	  IMIN = IMIN + 2
	  JMIN = JMIN - 2
	  TEMPR = WR
	  WR = WR * WSTPR - WI * WSTPI + WR
	  WI = TEMPR * WSTPI + WI * WSTPR + WI
  725	  IF (IMIN-JMIN) 710,730,740
  730	  IF (ISIGN) 731,740,740
  731	  CONTINUE
	  DO 735 I = IMIN,NTOT,NP2
  735	  DATA(I+1) = -DATA(I+1)
  740	  NP2 = NP2 + NP2
	  NTOT = NTOT + NTOT
	  J = NTOT + 1
	  IMAX = NTOT / 2 + 1
  745	  IMIN = IMAX - 2 * NHALF
	  I = IMIN
	  GO TO 755
  750	  DATA(J) = DATA(I)
	  DATA(J+1) = -DATA(I+1)
  755	  I = I + 2
	  J = J - 2
	  IF (I-IMAX) 750,760,760
  760	  DATA(J) = DATA(IMIN) - DATA(IMIN+1)
	  DATA(J+1) = 0.
	  IF (I-J) 770,780,780
  765	  DATA(J) = DATA(I)
	  DATA(J+1) = DATA(I+1)
  770	  I = I - 2
	  J = J - 2
	  IF (I-IMIN) 775,775,765
  775	  DATA(J) = DATA(IMIN) + DATA(IMIN+1)
	  DATA(J+1) = 0.
	  IMAX = IMIN
	  GO TO 745
  780	  DATA(1) = DATA(1) + DATA(2)
	  DATA(2) = 0.
	  GO TO 900
C
C	COMPLETE A REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION
C	BY CONJUGATE SYMMETRIES.
C
  800	  IF (I1RNG-NP1) 805,900,900
  805	  CONTINUE
	  DO 860 I3 = 1,NTOT,NP2
	    I2MAX = I3 + NP2 - NP1
	    DO 860 I2 = I3,I2MAX,NP1
	      IMIN = I2 + I1RNG
	      IMAX = I2 + NP1 - 2
	      JMAX = 2 * I3 + NP1 - IMIN
	      IF (I2-I3) 820,820,810
  810	      JMAX = JMAX + NP2
  820	      IF (IDIM-2) 850,850,830
  830	      J = JMAX + NP0
	      DO 840 I = IMIN,IMAX,2
		DATA(I) = DATA(J)
		DATA(I+1) = -DATA(J+1)
  840	      J = J - 2
  850	      J = JMAX
	      DO 860 I = IMIN,IMAX,NP0
		DATA(I) = DATA(J)
		DATA(I+1) = -DATA(J+1)
  860	  J = J - NP0
C
C	END OF LOOP FOR EACH DIMENSION
C
  900	  NP0 = NP1
	  NP1 = NP2
  910	NPREV = N
  920	RETURN
	END


c     *************************************************************************
c     *************************************************************************


	subroutine matinv(a,n,np,y)

c  Invert the double precision nxn matrix a; physical dimensions are np x np.
c  Inverse is returned in y.  WARNING: a is destroyed!

	double precision a(np,np),y(np,np)
	parameter (nmax=1000)
	dimension indx(nmax)
c
	if (n.gt.nmax) then
	  write(*,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	  write(1,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	end if
	  do 12 j=1,n
	    do 11 i=1,n
	    y(i,j)=0.0d0
11	  continue
	  y(j,j)=1.0d0
12	continue
	call ludcmp(a,n,np,indx)
	  do 13 j=1,n
	  call lubksb(a,n,np,indx,y(1,j))
13	continue
	return
	end


c     *************************************************************************
c     *************************************************************************


	subroutine ludcmp(a,n,np,indx)

	parameter (tiny=0.0d-20,nmax=1000)
	double precision a(np,np),sum,temp
	double precision aamax,vv(nmax)
	dimension indx(n)

	if (n.gt.nmax) then
	  write(*,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	  write(1,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	end if
	  do 12 i=1,n
	  aamax=0.0d0
	  do 11 j=1,n
	  if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11        continue
	  if (aamax.eq.0.d0) pause 'singular matrix.'
	  vv(i)=1./aamax
12	continue
	  do 19 j=1,n
	  if (j.gt.1) then
	      do 14 i=1,j-1
	      sum=a(i,j)
	      if (i.gt.1) then
	          do 13 k=1,i-1
	          sum=sum-a(i,k)*a(k,j)
13	        continue
	        a(i,j)=sum
	      endif
14	    continue
	  endif
	  aamax=0.0d0
	    do 16 i=j,n
	    sum=a(i,j)
	    if (j.gt.1) then
	        do 15 k=1,j-1
	        sum=sum-a(i,k)*a(k,j)
15	      continue
	      a(i,j)=sum
	    endif
	    temp=vv(i)*abs(sum)
	    if (temp.ge.aamax) then
	      imax=i
	      aamax=temp
	    endif
16	  continue
	  if (j.ne.imax) then
	      do 17 k=1,n
	      temp=a(imax,k)
	      a(imax,k)=a(j,k)
	      a(j,k)=temp
17	    continue
	    vv(imax)=vv(j)
	  endif
	  indx(j)=imax
	  if (j.ne.n) then
c	    if (a(j,j).eq.0.d0) a(j,j)=tiny
	    temp=1.0d0/a(j,j)
	      do 18 i=j+1,n
	      a(i,j)=a(i,j)*temp
18	    continue
	  endif
19	continue
c	if (a(n,n).eq.0.0d0) a(n,n)=tiny
	return
	end


c     *************************************************************************
c     *************************************************************************


	subroutine lubksb(a,n,np,indx,b)

	double precision a(np,np),b(np),sum
	dimension indx(n)

	ii=0
	  do 12 i=1,n
	  ll=indx(i)
	  sum=b(ll)
	  b(ll)=b(i)
	  if (ii.ne.0) then
	      do 11 j=ii,i-1
	      sum=sum-a(i,j)*b(j)
11	    continue
	  else if (sum.ne.0.0d0) then
	    ii=i
	  endif
	  b(i)=sum
12	continue
	  do 14 i=n,1,-1
	  sum=b(i)
	  if (i.lt.n) then
	      do 13 j=i+1,n
	      sum=sum-a(i,j)*b(j)
13	    continue
	  endif
	  b(i)=sum/a(i,i)
14	continue
	return
	end


C     **************************************************************************
C     **************************************************************************


      SUBROUTINE JACOBI(A,N,D,V,NROT)

c     Computes all eigenvalues and eigenvectors of a real symmetric matrix A,
c     which is of size N by N, stored in a physical MAXDIM*MAXDIM array. On 
c     output, elements of A above the diagonal are destroyed. D returns the 
c     eigenvalues of A in its first N elements. V is a matrix with the same 
c     logical and physical dimensions as A whose columns contain, on output, 
c     the normalized eigenvectors of A. NROT returns the number of Jacobi 
c     rotations which were required. 

      parameter  (nmax=100,maxdim=3)

      double precision a(maxdim,maxdim),d(maxdim),v(maxdim,maxdim)
      double precision b(nmax),z(nmax)
      logical   lcond


c     ------------------------------
c     Initialize the identity matrix
c     ------------------------------
        do 12 ip=1,n
          do 11 iq=1,n
            v(ip,iq)=0.
11        continue
          v(ip,ip)=1.0
12      continue

c       ---------------------------------------
c       Initialize B and D to the diagonal of A
c       ---------------------------------------
        do 13 ip=1,n
          b(ip)=a(ip,ip)
          d(ip)=b(ip)
c         ------------------------------------------------------------------
c         This vector will accumulate terms of the form ta_pq as in equation
c         (11.1.14).
c         ------------------------------------------------------------------
          z(ip)=0.0
13      continue

        nrot=0

        do 24 i=1,50
          sm=0.0
c         -------------------------
c         Sum off-diagonal elements
c         -------------------------
          do 15 ip=1,n-1
            do 14 iq=ip+1,n
              sm=sm+abs(a(ip,iq))
14          continue
15        continue

          if (sm.eq.0.) return
c         -----------------------------------------------------------
c         The normal return, which relies on quadratic convergence to 
c         machine underflow.
c         -----------------------------------------------------------

          if (i.lt.4) then
            thresh=0.2*sm/n**2
c           ----------------------------
c           ...on the first three sweeps
c           ----------------------------
          else
            thresh=0.
c           --------------
c           ...thereafter.
c           --------------
          endif

          do 22 ip=1,n-1
            do 21 iq=ip+1,n
              g=100.*abs(a(ip,iq))
c             ----------------------------------------------------------------
c             After four sweeps, skip the rotation if the off-diagonal element
c             is small.
c             ----------------------------------------------------------------
              lcond=((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))))
              lcond=(lcond.and.(abs(d(iq))+g.eq.abs(d(iq))))
              if (lcond) then
                a(ip,iq)=0.0
              else if (abs(a(ip,iq)).gt.thresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g.eq.abs(h)) then
                  t=a(ip,iq)/h
c                 -------------
c                 t=1/(2 theta)
c                 -------------
                else
                  theta=0.5*h/a(ip,iq)
c                 -------------------
c                 Equation (11.1.10).
c                 -------------------
                  t=1./(abs(theta)+sqrt(1.+theta**2))
                  if (theta.lt.0.) t=-t
                endif
                c=1./sqrt(1+t**2)
                s=t*c
                tau=s/(1.+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.
  
c               -----------------------------
c               Case of rotations 1 <= j < p.
c               -----------------------------
                do 16 j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
16              continue

c               ------------------------
c               Case of rotations p<j<q.
c               ------------------------
                do 17 j=ip+1,iq-1
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
17              continue

c               -------------------------
c               Case of rotations q<j<=n.
c               -------------------------
                do 18 j=iq+1,n
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
18              continue

                do 19 j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
19              continue
                nrot=nrot+1
              endif
21          continue
22        continue

c         ------------------------------------------------
c         Update D with the sum ta_pq, and reinitialize Z.
c         ------------------------------------------------
          do 23 ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.0
23        continue
24      continue

        pause '50 iterations should never happen'

       return
      end     


C     **************************************************************************
C     **************************************************************************


      SUBROUTINE EIGSRT(D,V,N)


c     -------------------------------------------------------------------
c     Given the eigenvalues D and eigenvectors V as output from JACOBI or 
c     TQLI this routine sorts the eigenvalues into descending order, and 
c     rearranges the columns of V correspondingly. The method is straight 
c     insertion.
c     -------------------------------------------------------------------

      integer   maxdim
      parameter (maxdim=3)

      double precision  d(maxdim),v(maxdim,maxdim)


        do 13 i=1,n-1
          k=i
          p=d(i)
          do 11 j=i+1,n
            if (d(j).gt.p) then
              k=j
              p=d(j)
            endif
11        continue
          if (k.ne.i) then
            d(k)=d(i)
            d(i)=p
            do 12 j=1,n
              p=v(j,i)
              v(j,i)=v(j,k)
              v(j,k)=p
12          continue
          endif
13      continue

        if (v(3,3).lt.0.0) then
          do 15 m=1,3
            v(m,3)=-v(m,3)
15        continue
        endif

       return
      end


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE QROMO(FUNC,A,B,SS,CHOOSE)

C     ______________________________________________________________________
C     Romberg integration on an open interval. Returns as SS the integral
C     of the function FUNC from A to B, using any specified integrating 
C     subroutine CHOOSE and Romberg's method. Normally CHOOSE will be an 
C     open formula, not evaluating the function at the endpoints. It is 
C     assumed that CHOOSE triples the number of steps on each call, and that
C     its error series contains only even POWINTs of the number of steps. 
C     The routines MIDPNT,MIDINF,MIDSQL,MIDSQU, are possible choices for 
C     CHOOSE.
C     ----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      EXTERNAL  CHOOSE

      PARAMETER (EPS=1.0D-6,JMAX=15,JMAXP=JMAX+1,K=7,KM=K-1)
C     _________________________________________________________________
C     Here EPS is the fractional accuracy desired, as determined by the 
C     extrapolation error estimate; JMAX limits the total number of 
C     steps; K is the number of points used in the extrapolation. 
C     -----------------------------------------------------------------
      DIMENSION S(JMAXP),H(JMAXP)


        H(1)=1.0D0
        DO 11 J=1,JMAX
          CALL CHOOSE(FUNC,A,B,S(J),J)
          IF (J.GE.K) THEN
            CALL POLINT(H(J-KM),S(J-KM),K,0.0D0,SS,DSS)
            IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          ENDIF
          S(J+1)=S(J)
          H(J+1)=H(J)/9.0D0
C         _______________________________________________________________
C         This is where the assumption of step tripling and an even error
C         series is used.
C         ---------------------------------------------------------------
11      CONTINUE
        PAUSE 'Too many steps.'
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE MIDPNT(FUNC,A,B,S,N)


C     ____________________________________________________________________
C     This routine computes the N'th stage of refinement of an extended 
C     midpoint rule. FUNC is input as the name of the function to be 
C     integrated between limits A and B, also input. When called with 
C     N=1, the routine returns as S the crudest estimate of I[a,b] f(x)dx.
C     Subsequent calls with N=2,3,... (in that sequential order) will 
C     improve the accuracy of S by adding (2/3)*(3**(N-1)) additional 
C     interior points. S should not be modified between sequential calls.
C     --------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      EXTERNAL  FUNC
      

        IF (N.EQ.1) THEN
          S=(B-A)*FUNC(0.5D0*(A+B))
          IT=1
C         ________________________________________________
C         2*IT points will be added on the next refinement
C         ------------------------------------------------
        ELSE
          TNM=IT
          DEL=(B-A)/(3*TNM)
          DDEL=DEL+DEL
C         __________________________________________________________
C         The added points alternate in spacing between DEL and DDEL
C         ----------------------------------------------------------
          X=A+0.5D0*DEL
          SUM=0.0D0
          DO 11 J=1,IT
            SUM=SUM+FUNC(X)
            X=X+DDEL
            SUM=SUM+FUNC(X)
            X=X+DEL
11        CONTINUE
          S=(S+(B-A)*SUM/TNM)/3.0D0
C         ________________________________________________
C         The new sum is combined with the old integral to 
C         give a refined integral.
C         ------------------------------------------------
          IT=3*IT
        ENDIF
       RETURN
      END

C     *************************************************************************
C     *************************************************************************


      SUBROUTINE MIDINF(FUNK,AA,BB,S,N)


C     _____________________________________________________________________
C     This routine is an exact replacement for MIDPNT, ie. returns as S the 
C     Nth stage of refinement of the integral of FUNK from AA to BB, except 
C     that the function is evaluated at evenly spaced points in 1/x, rather 
C     than in x. This allows the upper limit BB to be as large and positive
C     as the computer allows, or the lower limit AA to be as large and 
C     negative, but not both. AA and BB must have the same sign.
C     ---------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      EXTERNAL   FUNK


        FUNC(X)=FUNK(1.0D0/X)/X**2
C       __________________________________________________________________
C       This is a statement function which effects the change of variable.
C       ------------------------------------------------------------------
        B=1.0D0/AA
        A=1.0D0/BB
C       __________________________________________________________________
C       These two statements change the limits of integration accordingly.
C       ------------------------------------------------------------------
C       _______________________________________________________________
C       From this point on, the routine is exactly identical to MIDPNT.
C       ---------------------------------------------------------------
        IF (N.EQ.1) THEN
          S=(B-A)*FUNC(0.5D0*(A+B))
          IT=1
        ELSE
          TNM=IT
          DEL=(B-A)/(3.0D0*TNM)
          DDEL=DEL+DEL
          X=A+0.5D0*DEL
          SUM=0.0D0
          DO 11 J=1,IT
            SUM=SUM+FUNC(X)
            X=X+DDEL
            SUM=SUM+FUNC(X)
            X=X+DEL
11        CONTINUE
          S=(S+(B-A)*SUM/TNM)/3.0D0
          IT=3*IT
        ENDIF     
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

C     ___________________________________________________________________
C     Given arrays XA and YA, each of length N, and given a value X, this 
C     routine returns a value Y, and an error estimate DY. If P(X) is the
C     polynomial of degree N-1 such that P(XA{i})=YA{i}), i=1,...,N, then
C     the returned value Y=P(X).
C     -------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (NMAX=10)
C     ________________________________________________________________
C     Change NMAX as desired to be the largest anticipated value of N.
C     ----------------------------------------------------------------

      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)


        NS=1
        DIF=ABS(X-XA(1))
        DO 11 I=1,N
C         _____________________________________________________
C         Here we find the index NS of the closest table entry.
C         -----------------------------------------------------
          DIFT=ABS(X-XA(I))
          IF (DIFT.LT.DIF) THEN
            NS=I
            DIF=DIFT
          ENDIF
C         __________________________________________
C         and initialize the tableau of C's and D's.
C         ------------------------------------------
          C(I)=YA(I)
          D(I)=YA(I)
11      CONTINUE
C       _______________________________________
C       This is the initial approximation to Y:
C       ---------------------------------------
        Y=YA(NS)
        NS=NS-1
        DO 13 M=1,N-1
C         ___________________________________________
C         For each column of the tableau we loop over 
C         the current C's and D's and update them.
C         -------------------------------------------
          DO 12 I=1,N-M
            HO=XA(I)-X
            HP=XA(I+M)-X
            W=C(I+1)-D(I)
            DEN=HO-HP
            IF (DEN.EQ.0.0D0) PAUSE
C           __________________________________________________________
C           This error can occur only if two input XA's are (to within 
C           roundoff) identical.
C           ----------------------------------------------------------
            DEN=W/DEN
C           _________________________________
C           Here the C's and D's are updated.
C           ---------------------------------
            D(I)=HP*DEN
            C(I)=HO*DEN
12        CONTINUE
C         ________________________________________________________________
C         After each column in the tableau is completed, we decide which 
C         correction, C or D, we want to add to our accumulating value of
C         Y, ie. which path to take through the tablau-forking up or down. 
C         We do this in such a way as to take the most "straight line" 
C         route through the tableau to its apex, updating NS accordingly 
C         to keep track of where we are. This route keeps the partial 
C         approximation centered (insofar as possible) on the target X.
C         The last DY added is thus the error indication.
C         ----------------------------------------------------------------
          IF (2*NS.LT.N-M) THEN
            DY=C(NS+1)
          ELSE
            DY=D(NS)
            NS=NS-1
          ENDIF
          Y=Y+DY
13      CONTINUE
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE LOCATE(XX,N,X,J)

C     _________________________________________________________________________
C     Given an array XX of length N, and given a value X, returns a value
C     J such that X is between XX(J) and XX(J+1). XX must be monotonic, either 
C     increasing or decreasing. J=0 or J=N is returned to indicate that X
C     is out of range.
C     Taken from Numerical Recipes, pg. 90; Press, Flannery, Teukolsky,
C     and Vetterling. 
C     -------------------------------------------------------------------------

      DOUBLE PRECISION  XX(N),X

        JL=0
        JU=N+1
C         _________________________________
C         Initialize lower and upper limits
C         ---------------------------------

10      IF (JU-JL.GT.1) THEN
          JM=(JU+JL)/2
C           __________________________________________
C           If we are not yet done compute a midpoint,
C           ------------------------------------------
          IF (XX(N).GT.XX(1)) THEN
            IF (X.GT.XX(JM)) THEN
C           ______________________________________________________
C           and replace either the lower limit or the upper limit,
C           as appropiate.
C           ------------------------------------------------------
              JL=JM
            ELSE
              JU=JM
            ENDIF
          ELSE
            IF (X.LT.XX(JM)) THEN
C           ______________________________________________________
C           and replace either the lower limit or the upper limit,
C           as appropiate.
C           ------------------------------------------------------
              JL=JM
            ELSE
              JU=JM
            ENDIF
          ENDIF
          GO TO 10
C           _______________________________________________
C           Repeat until the test condition 10 is satisfied
C           -----------------------------------------------
        ENDIF

C       ______________________________
C       Then set the output and return
C       ------------------------------
        J=JL
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE SORT(ARRAY,IARRAY,IROW,STORE,ISTORE)


      PARAMETER (NEMAX=257,NPMAX=257,NEPMAX=NEMAX*NPMAX)

      INTEGER           IROW
      INTEGER           IARRAY(2,NEPMAX),ISTORE(2)
      DOUBLE PRECISION  ARRAY(NEPMAX),STORE


        IF (IROW.LT.1000) THEN
          CALL SHELLS(ARRAY,IARRAY,IROW)
        ELSE
          CALL HEAPSR(ARRAY,IARRAY,IROW,STORE,ISTORE)
        ENDIF
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE SHELLS(ARRAY,IARRAY,IROW)


C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     +                                                                       +
C     + Sorts an (REAL) array ARRAY with IROW rows and ICOL columns into      +
C     + descending numerical order of column no. ISRT, by the Shell-Mezgar    +
C     + algorithm (diminishing increment sort).                               +
C     + ISRT, IROW and ICOL are input; ARRAY is replaced on output by its     +
C     + sorted rearrangement.                                                 +
C     +                                                                       +
C     + This program is based on the program presented on page 229 in:        +
C     +      Numerical Recipes                                                +
C     +      The Art of Scientific Computing                                  +
C     +      Press, Flannery, Teukolsky, Vetterling                           +
C     +      Cambridge University Press, 1986                                 +
C     +                                                                       +
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



      PARAMETER (NEMAX=257,NPMAX=257,NEPMAX=NEMAX*NPMAX)

      INTEGER            IROW,NSQR,NGAP,NCTR
      INTEGER            IARRAY(2,NEPMAX)
      DOUBLE PRECISION   ARRAY(NEPMAX)


        NSQR=2
        NGAP=IROW
1001    IF (NSQR.LE.IROW) THEN
C         ----------------------------
C         Loop over the partial sorts.
C         ----------------------------
          NGAP=NGAP/2
          NCTR=IROW-NGAP
          DO 10 J=1,NCTR
C           ---------------------------------
C           Outer loop of straight insertion.
C           ---------------------------------
            I=J
1002        CONTINUE
C             ---------------------------------
C             Inner loop of straight insertion.
C             ---------------------------------
            L=I+NGAP
            IF (ARRAY(L).GT.ARRAY(I)) THEN
              RHELP=ARRAY(I)
              ARRAY(I)=ARRAY(L)
              ARRAY(L)=RHELP
              IHELP=IARRAY(1,I)
              IARRAY(1,I)=IARRAY(1,L)
              IARRAY(1,L)=IHELP
              IHELP=IARRAY(2,I)
              IARRAY(2,I)=IARRAY(2,L)
              IARRAY(2,L)=IHELP
              I=I-NGAP
              IF (I.GE.1) GO TO 1002
            ENDIF
10        CONTINUE
          NSQR=NSQR*2
          GOTO 1001
        ENDIF
       RETURN
      END


C     *************************************************************************
C     *************************************************************************


      SUBROUTINE HEAPSR(ARRAY,IARRAY,IROW,STORE,ISTORE)


C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     +                                                                       +
C     + Sorts an (REAL) array ARRAY of IROW rows and ICOL columns into        +
C     + descending numerical order of column no. ISRT, using the Heapsort     +
C     + algorithm.                                                            +
C     + ISRT, IROW and ICOL are input; ARRAY is replaced on output by its     +
C     + sorted rearrangement.                                                 +
C     +                                                                       +
C     + This program is based on the program presented on page 229 in:        +
C     +      Numerical Recipes                                                +
C     +      The Art of Scientific Computing                                  +
C     +      Press, Flannery, Teukolsky, Vetterling                           +
C     +      Cambridge University Press, 1986                                 +
C     +                                                                       +
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      PARAMETER (NEMAX=257,NPMAX=257,NEPMAX=NEMAX*NPMAX)

      INTEGER           IROW,IHIRE,IPROM
      INTEGER           IARRAY(2,NEPMAX),ISTORE(2)
      DOUBLE PRECISION  ARRAY(NEPMAX),STORE


        IHIRE=IROW/2+1
        IPROM=IROW
C         --------------------------------------------------------------------
C         The index IPROM will be decremented from its initial value down to 1
C         during the "hiring" (heap creation) phase. Once it reaches 1, the
C         index IR will be decremented from its initial value down to 1 during
C         the "retirement-and-promotion" (heap selection) phase.
C         --------------------------------------------------------------------
1001    CONTINUE
          IF (IHIRE.GT.1) THEN
C           ----------------------
C           Still in hiring phase.
C           ----------------------
            IHIRE=IHIRE-1
            STORE=ARRAY(IHIRE)
            ISTORE(1)=IARRAY(1,IHIRE)
            ISTORE(2)=IARRAY(2,IHIRE)
          ELSE
C           ----------------------------------
C           In retirement-and-promotion phase.
C           ----------------------------------
            STORE=ARRAY(IPROM)
            ARRAY(IPROM)=ARRAY(1)
            ISTORE(1)=IARRAY(1,IPROM)
            IARRAY(1,IPROM)=IARRAY(1,1)
            ISTORE(2)=IARRAY(2,IPROM)
            IARRAY(2,IPROM)=IARRAY(2,1)
C           ------------------------------------------------------------------
C           Clear a space at end of array, retire the top of the heap into it.
C           ------------------------------------------------------------------
            IPROM=IPROM-1
C           ---------------------------------------
C           Decrease the size of the "corporation".
C           ---------------------------------------
            IF (IPROM.EQ.1) THEN
C             -----------------------------
C             Done with the last promotion.
C             -----------------------------
              ARRAY(1)=STORE
              IARRAY(1,1)=ISTORE(1)
              IARRAY(2,1)=ISTORE(2)
C             ----------------------------------
C             The least competent worker of all!
C             ----------------------------------
              RETURN
            ENDIF
          ENDIF
C         ------------------------------------------------------------------
C         Whether we are in the hiring phase or promotion phase, we here set
C         up to sift down element STORE to its proper level.
C         ------------------------------------------------------------------
          I=IHIRE
          J=IHIRE+IHIRE
1002      IF (J.LE.IPROM) THEN
C           ---------------------
C           "Do while J.LE.IPROM"
C           ---------------------
            IF (J.LT.IPROM) THEN
              IF (ARRAY(J).GT.ARRAY(J+1)) THEN
C               --------------------------------
C               Compare to the better underling.
C               --------------------------------
                J=J+1
              ENDIF
            ENDIF
            IF (STORE.GT.ARRAY(J)) THEN
C             -------------
C             Demote STORE.
C             -------------
              ARRAY(I)=ARRAY(J)
              IARRAY(1,I)=IARRAY(1,J)
              IARRAY(2,I)=IARRAY(2,J)
              I=J
              J=J+J
            ELSE
C             ----------------------------------------------------
C             This is STORE's level. Set J to terminate sift-down.
C             ----------------------------------------------------
              J=IPROM+1
            ENDIF
            GOTO 1002
          ENDIF
          ARRAY(I)=STORE
          IARRAY(1,I)=ISTORE(1)
          IARRAY(2,I)=ISTORE(2)
C         ------------------------
C         Put STORE into its slot.
C         ------------------------
        GOTO 1001
      END


c     *************************************************************************
c     *************************************************************************


      subroutine randini


c     -----------------------------------------------------------------------
c     randini is used to initialize the random number generators.  iseed is a
c     positive integer less than 1.0e9.  The basic random number generator is
c     taken from Press et al., Numerical Recipes, p. 199 and is based on
c     Knuth's suggestion for a portable random number generator.  mseed is
c     any large number less than m=1.0e9.
c     -----------------------------------------------------------------------

      parameter (m=1000000000,mseed=161803398,rm=1.0/m)

c     The dimension 55 is special and should not be modified; see Knuth.

      integer     marr(55)
      real        randar(500),xnorm

      common /randnos/ randar,irptr,nrand,marr,inext,inextp
      common /randnor/ xnorm,inorm
      save /randnos/,/randnor/


	nrand=157
        write(*,*) 'Initialization of random number generator: '           
        write(1,*) 'Initialization of random number generator: '        
        write(*,*)
        write(*,'('' Enter random number seed (9-digit integer): '',$)')
        read(*,*) isee
        write(1,*)
        write(1,'('' Enter random number seed (9-digit integer): '',$)')
        write(1,'(i10)') isee

	iseed=isee
	iseed=mod(iseed,m)
c  Initialize marr(55).
	mj=mseed-iseed
	mj=mod(mj,m)
	if (mj.lt.0) mj=mj+m
	marr(55)=mj
	mk=1
c  Now initialize the rest of the table, in a slightly random order,
c  with numbers that are not especially random.
	  do 10 i=1,54
	  ii=mod(21*i,55)
	  marr(ii)=mk
	  mk=mj-mk
	  if (mk.lt.0) mk=mk+m
	  mj=marr(ii)
10	continue
c  Randomize them by "warming up the generator."
	  do 30 k=1,4
	    do 20 i=1,55
	    marr(i)=marr(i)-marr(1+mod(i+30,55))
	    if (marr(i).lt.0) marr(i)=marr(i)+m
20	  continue
30	continue
	inext=0
	inextp=31
c  Exercise generator before storing in shuffling table.
	  do 40 i=1,nrand
	  inext=inext+1
	  if (inext.eq.56) inext=1
	  inextp=inextp+1
	  if (inextp.eq.56) inextp=1
	  mj=marr(inext)-marr(inextp)
	  if (mj.lt.0) mj=mj+m
	  marr(inext)=mj
40	continue
c  Now fill shuffling table.
	  do 50 i=1,nrand
	  inext=inext+1
	  if (inext.eq.56) inext=1
	  inextp=inextp+1
	  if (inextp.eq.56) inextp=1
	  mj=marr(inext)-marr(inextp)
	  if (mj.lt.0) mj=mj+m
	  marr(inext)=mj
	  randar(i)=mj*rm
50	continue
	inext=inext+1
	if (inext.eq.56) inext=1
	inextp=inextp+1
	if (inextp.eq.56) inextp=1
	mj=marr(inext)-marr(inextp)
	if (mj.lt.0) mj=mj+m
	marr(inext)=mj
	irptr=int(mj*rm*nrand)+1
	xnorm=0.0
	inorm=0

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine randa(x)


c     -------------------------------------------------------------
c     randa generates uniform random numbers in the interval [0,1).
c     It must be initialized with a call to randini.
c     -------------------------------------------------------------

      parameter (m=1000000000,mseed=161803398,rm=1.0/m)

c     ------------------------------------------------------------------
c     The dimension 55 is special and should not be modified; see Knuth.
c     ------------------------------------------------------------------

      integer     marr(55)
      real        randar(500),x

      common /randnos/ randar,irptr,nrand,marr,inext,inextp
      save /randnos/

c       -------------------------------------------
c       Extract random number from shuffling table.
c       -------------------------------------------
	x=randar(irptr)
	irptr=int(nrand*x)+1
c       -----------------------------
c       Generate a new random number.
c       -----------------------------
	inext=inext+1
	if (inext.eq.56) inext=1
	inextp=inextp+1
	if (inextp.eq.56) inextp=1
	mj=marr(inext)-marr(inextp)
	if (mj.lt.0) mj=mj+m
	marr(inext)=mj
	randar(irptr)=mj*rm

       return
      end


c     *************************************************************************
c     *************************************************************************


      subroutine jimcls


c     --------------------------
c     Closing remarks of program
c     --------------------------


        write(*,*)
        write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'    
        write(*,*)
        write(*,*)
        write(*,*)
        write(*,*) '                        This is the end'
        write(*,*) '                        My only friend, the end'
        write(*,*)
        write(*,*) '                               Jim Morrison, 1967'
        write(*,*)

        write(1,*)
        write(1,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     2             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'    
        write(1,*)
        write(1,*)
        write(1,*)
        write(1,*) '                        This is the end'
        write(1,*) '                        My only friend, the end'
        write(1,*)
        write(1,*) '                               Jim Morrison, 1967'

       return
      end






