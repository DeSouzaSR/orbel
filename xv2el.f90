
program el2xv
    ! Purpose: Convert cartesian coordinates and
    !   velocities to orbital elements
    !
    ! Usage:
    !   1) $ xv2el ! prompt ialpha, gm, a, e, inc, capom, omega, capm
    !   2) $ xv2el < input_file ! screen output
    !   3) $ xv2el < input_file > output_file ! output to file
    !
    ! Input:
    !   ialpha      -- conic section type: hyperbole (+1), parabola (0)
    !   and ellipse (-1).
    !   gm          -- factor solar mass. ex: 39.476926421373015 [au^3 a^-2]
    !   x, y, v     -- cartesian coordinats [ex. au]
    !   vx, vy, vz  -- velocities [ex. au a^-1]
    !
    ! Output:
    !   gm      -- factor solar mass [ex.: 39.476926421373015 au^3 a^-2]
    !   a       -- semi-major axis [ex.: au]
    !   e       -- eccentricity
    !   inc     -- inclination [deg]
    !   capom   -- longitude of the ascending node [deg]
    !   omega   -- argument of periapsis [deg]
    !   capm    -- mean anomaly [deg]
    !

    ! Define parameters and constants
    implicit none
    real, parameter :: PI = 4. * atan(1.) ! maximum precision
    real, parameter :: DEG2RAD = PI / 180.
    real, parameter :: RAD2DEG = 180. / PI

    ! Define variables
    character(len=20), parameter :: fmt = "(6f12.8)"
    integer :: ialpha
    real*8 :: gm
    real*8 :: a
    real*8 :: e
    real*8 :: inc
    real*8 :: capom
    real*8 :: omega
    real*8 :: capm
    real*8 :: x, y, z, vx, vy, vz

    ! Input data
    read(*,*) ialpha, gm, x, y, z, vx, vy, vz

    call orbel_xv2el(x, y, z, vx, vy, vz, gm, ialpha,a,e,inc,capom,omega,capm)

    ! Convert to radian
    inc = inc * RAD2DEG
    capom = capom * RAD2DEG
    omega = omega * RAD2DEG
    capm  = capm  * RAD2DEG

    write(*,fmt) a, e, inc, capom, omega, capm

end program el2xv


!----------------------------------------------------------------------
! Subroutines and functions
!----------------------------------------------------------------------

!***********************************************************************
!c                    ORBEL_XV2EL.F
!***********************************************************************
!!!*     PURPOSE:  Given the cartesian position and velocity of an orbit,
!!!*       compute the osculating orbital elements.
!*
!C       input:
!c            x,y,z    ==>  position of object (real scalars)
!c            vx,vy,vz ==>  velocity of object (real scalars)
!c            gmsum       ==> G*(M1+M2) (real scalar)
!c
!c       Output:
!c         ialpha   ==> conic section type ( see PURPOSE, integer scalar)
!C         a        ==> semi-major axis or pericentric distance if a parabola
!c                          (real scalar)
!c            e        ==> eccentricity (real scalar)
!C            inc      ==> inclination  (real scalar)
!C            capom    ==> longitude of ascending node (real scalar)
!C         omega    ==> argument of perihelion (real scalar)
!C         capm     ==> mean anomoly(real scalar)
!c
!!!*     ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
!!*     REMARKS:  If the inclination INC is less than TINY, we
!!*       arbitrarily choose the longitude of the ascending node LGNODE
!!*       to be 0.0 (so the ascending node is then along the X axis).  If
!!*       the  eccentricity E is less than SQRT(TINY), we arbitrarily
!!*       choose the argument of perihelion to be 0.
!!*     AUTHOR:  M. Duncan.
!!*     DATE WRITTEN:  May 8,1992.
!!*     REVISIONS: 12/8/2011
!***********************************************************************

subroutine orbel_xv2el(x,y,z,vx,vy,vz,gmsum, ialpha,a,e,inc,capom,omega,capm)

    include './swift.inc'

    ! !c...  Inputs Only:
    real*8 x,y,z,vx,vy,vz,gmsum

    !c...  Outputs
    integer ialpha
    real*8 a,e,inc,capom,omega,capm

    !c...  Internals:
    real*8 hx,hy,hz,h2,h,r,v2,v,vdotr,energy,fac,face,cape,capf,tmpf
    real*8 cw,sw,w,u

    !c----
    !c...  Executable code

    !* Compute the angular momentum H, and thereby the inclination INC.

    hx = y*vz - z*vy
    hy = z*vx - x*vz
    hz = x*vy - y*vx
    h2 = hx*hx + hy*hy +hz*hz
    h  = sqrt(h2)
    if(hz.gt.h) then                 ! Hal's fix
        hz = h
        hx = 0.0d0
        hy = 0.0d0
    endif
    inc = acos(hz/h)

    !* Compute longitude of ascending node CAPOM and the argument of
    !* latitude u.
    fac = sqrt(hx**2 + hy**2)/h

    if( (fac.lt. TINY ) .or. (inc.eq.0.0d0) ) then ! Hal's fix
        capom = 0.d0
        u = atan2(y,x)
        if(abs(inc - PI).lt. 10.d0*TINY) u = -u
    else
        capom = atan2(hx,-hy)
        u = atan2 ( z/sin(inc) , x*cos(capom) + y*sin(capom))
    endif

    if(capom .lt. 0.d0) capom = capom + 2.d0*PI
    if(u .lt. 0.d0) u = u + 2.d0*PI

    !*  Compute the radius R and velocity squared V2, and the dot
    !*  product RDOTV, the energy per unit mass ENERGY .

    r = sqrt(x*x + y*y + z*z)
    v2 = vx*vx + vy*vy + vz*vz
    v = sqrt(v2)
    vdotr = x*vx + y*vy + z*vz
    energy = 0.5d0*v2 - gmsum/r

    !*  Determine type of conic section and label it via IALPHA
    if(abs(energy*r/gmsum) .lt. sqrt(TINY)) then
    ialpha = 0
    else
        if(energy .lt. 0.d0) ialpha = -1
        if(energy .gt. 0.d0) ialpha = +1
    endif

!* Depending on the conic type, determine the remaining elements

!***
!c ELLIPSE :
      if(ialpha .eq. -1) then
        a = -0.5d0*gmsum/energy
        fac = 1.d0 - h2/(gmsum*a)

        if (fac .gt. TINY) then
           e = sqrt ( fac )
           face = (a-r)/(a*e)

    !c... Apr. 16/93 : watch for case where face is slightly outside unity
           if ( face .gt. 1.d0) then
              cape = 0.d0
           else
              if ( face .gt. -1.d0) then
                cape = acos( face )
              else
                cape = PI
              endif
           endif

           if ( vdotr .lt. 0.d0 ) cape = 2.d0*PI - cape
           cw = (cos( cape) -e)/(1.d0 - e*cos(cape))
           sw = sqrt(1.d0 - e*e)*sin(cape)/(1.d0 - e*cos(cape))
           w = atan2(sw,cw)
           if(w .lt. 0.d0) w = w + 2.d0*PI
        else
           e = 0.d0
           w = u
           cape = u
        endif

        capm = cape - e*sin (cape)
        omega = u - w
        if(omega .lt. 0.d0) omega = omega + 2.d0*PI
        omega = omega - int(omega/(2.d0*PI))*2.d0*PI
      endif
    !***
    !***
    !c HYPERBOLA
      if(ialpha .eq. +1) then

         a = +0.5d0*gmsum/energy
         fac = h2/(gmsum*a)

         if (fac .gt. TINY) then
           e = sqrt ( 1.d0 + fac )
           tmpf = (a+r)/(a*e)
           if(tmpf.lt.1.0d0) then
             tmpf = 1.0d0
           endif
           capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
           if ( vdotr .lt. 0.d0 ) capf = - capf
           cw = (e - cosh(capf))/(e*cosh(capf) - 1.d0 )
           sw = sqrt(e*e - 1.d0)*sinh(capf)/(e*cosh(capf) - 1.d0 )
           w = atan2(sw,cw)
           if(w .lt. 0.d0) w = w + 2.d0*PI
         else
    !c we only get here if a hyperbola is essentially a parabola
    !c so we calculate e and w accordingly to avoid singularities
           e = 1.d0
           tmpf = 0.5d0*h2/gmsum
           w = acos(2.d0*tmpf/r -1.d0)
           if ( vdotr .lt. 0.d0) w = 2.d0*PI - w
           tmpf = (a+r)/(a*e)
           capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
         endif

         capm = e * sinh(capf) - capf
         omega = u - w
         if(omega .lt. 0.d0) omega = omega + 2.d0*PI
         omega = omega - int(omega/(2.d0*PI))*2.d0*PI
      endif
    !***
    !***
    !c PARABOLA : ( NOTE - in this case we use "a" to mean pericentric distance)
      if(ialpha .eq. 0) then
         a =  0.5d0*h2/gmsum
         e = 1.d0
         w = acos(2.d0*a/r -1.d0)
         if ( vdotr .lt. 0.d0) w = 2.d0*PI - w
         tmpf = tan(0.5d0 * w)
         capm = tmpf* (1.d0 + tmpf*tmpf/3.d0)
         omega = u - w
         if(omega .lt. 0.d0) omega = omega + 2.d0*PI
         omega = omega - int(omega/(2.d0*PI))*2.d0*PI
      endif
    !***
    !***
      return
end    ! orbel_xv2el
!c------------------------------------------------------------------


!*****************************************************************************
!*                          ORBEL_EL2XV.F
!*****************************************************************************
!*     PURPOSE: To compute cartesian positions and velocities given
!*               central mass, ialpha ( = +1 for hyp., 0 for para. and
!*               -1 for ellipse), and orbital elements.
!C       input:
!c            gm       ==> G times central mass (real scalar)
!c         ialpha   ==> conic section type ( see PURPOSE, integer scalar)
!C         a        ==> semi-major axis or pericentric distance if a parabola
!c                          (real scalar)
!c            e        ==> eccentricity (real scalar)
!C            inc      ==> inclination  (real scalar)
!C            capom    ==> longitude of ascending node (real scalar)
!C         omega    ==> argument of perihelion (real scalar)
!C         capm     ==> mean anomoly(real scalar)
!!*
!c       Output:
!c            x,y,z    ==>  position of object (real scalars)
!c            vx,vy,vz ==>  velocity of object (real scalars)
!c
!!*     ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
!!*     REMARKS: All angles are in RADIANS
!!*
!!*     AUTHOR:  M. Duncan.
!!*     DATE WRITTEN:  May 11, 1992.
!!*     REVISIONS: May 26 - now use better Kepler solver for ellipses
!!*                 and hyperbolae called EHYBRID.F and FHYBRID.F
!***********************************************************************

subroutine orbel_el2xv(gm,ialpha,a,e,inc,capom,omega,capm,x,y,z,vx,vy,vz)

      include './swift.inc'

    !c...  Inputs Only:
      integer ialpha
      real*8 gm,a,e,inc,capom,omega,capm

    !c...  Outputs:
      real*8 x,y,z,vx,vy,vz

    !c...  Internals:
      real*8 cape,capf,zpara,em1
      real*8 sp,cp,so,co,si,ci
      real*8 d11,d12,d13,d21,d22,d23
      real*8 scap,ccap,shcap,chcap
      real*8 sqe,sqgma,xfac1,xfac2,ri,vfac1,vfac2
      real*8 orbel_ehybrid,orbel_fhybrid,orbel_zget

      !intent(in) gm,ialpha,a,e,inc,capom,omega,capm
      !intent(out) x,y,z,vx,vy,vz

    !c----
    !c...  Executable code

      if(e.lt.0.0) then
           write(*,*) ' ERROR in orbel_el2xv: e<0, setting e=0!!1'
           e = 0.0
      endif

    !c...    check for inconsistencies between ialpha and e
      em1 = e - 1.d0
      if( ((ialpha.eq.0) .and. (abs(em1).gt.TINY))  .or. &
        & ((ialpha.lt.0) .and. (e.gt.1.0d0))  .or. &
        & ((ialpha.gt.0) .and. (e.lt.1.0d0)) )  then
        write(*,*) 'ERROR in orbel_el2xv: ialpha and e inconsistent'
        write(*,*) '                       ialpha = ',ialpha
        write(*,*) '                            e = ',e
      endif

    !C Generate rotation matrices (on p. 42 of Fitzpatrick)
    !C
      call orbel_scget(omega,sp,cp)
      call orbel_scget(capom,so,co)
      call orbel_scget(inc,si,ci)
      d11 = cp*co - sp*so*ci
      d12 = cp*so + sp*co*ci
      d13 = sp*si
      d21 = -sp*co - cp*so*ci
      d22 = -sp*so + cp*co*ci
      d23 = cp*si

    !C--
    !C Get the other quantities depending on orbit type ( i.e. IALPHA)
    !C
      if (ialpha .eq. -1) then
      cape = orbel_ehybrid(e,capm)
      call orbel_scget(cape,scap,ccap)
      sqe = sqrt(1.d0 -e*e)
      sqgma = sqrt(gm*a)
      xfac1 = a*(ccap - e)
      xfac2 = a*sqe*scap
      ri = 1.d0/(a*(1.d0 - e*ccap))
      vfac1 = -ri * sqgma * scap
      vfac2 = ri * sqgma * sqe * ccap
      endif
    !c--
      if (ialpha .eq. +1) then
      capf = orbel_fhybrid(e,capm)
      call orbel_schget(capf,shcap,chcap)
      sqe = sqrt(e*e - 1.d0 )
      sqgma = sqrt(gm*a)
      xfac1 = a*(e - chcap)
      xfac2 = a*sqe*shcap
      ri = 1.d0/(a*(e*chcap - 1.d0))
      vfac1 = -ri * sqgma * shcap
      vfac2 = ri * sqgma * sqe * chcap
      endif
    !C--
      if (ialpha .eq. 0) then
      zpara = orbel_zget(capm)
      sqgma = sqrt(2.d0*gm*a)
      xfac1 = a*(1.d0 - zpara*zpara)
      xfac2 = 2.d0*a*zpara
      ri = 1.d0/(a*(1.d0 + zpara*zpara))
      vfac1 = -ri * sqgma * zpara
      vfac2 = ri * sqgma
      endif
    !C--
      x =  d11*xfac1 + d21*xfac2
      y =  d12*xfac1 + d22*xfac2
      z =  d13*xfac1 + d23*xfac2
      vx = d11*vfac1 + d21*vfac2
      vy = d12*vfac1 + d22*vfac2
      vz = d13*vfac1 + d23*vfac2

      return
end    ! orbel_el2xv

!c-----------------------------------------------------------------------



!***********************************************************************
!c                    ORBEL_EGET.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                           m ==> mean anomaly. (real scalar)
!!*             Returns:
!!*                  orbel_eget ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: Quartic convergence from Danby
!!*     REMARKS: For results very near roundoff, give it M between
!!*           0 and 2*pi. One can condition M before calling EGET
!!*           by calling my double precision function MOD2PI(M).
!!*           This is not done within the routine to speed it up
!!*           and because it works fine even for large M.
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 7, 1992.
!!*     REVISIONS: May 21, 1992.  Now have it go through EXACTLY two iterations
!!*                with the premise that it will only be called if
!!*             we have an ellipse with e between 0.15 and 0.8
!***********************************************************************

real*8 function orbel_eget(e,m)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,m

    !c...  Internals:
      real*8 x,sm,cm,sx,cx
      real*8 es,ec,f,fp,fpp,fppp,dx

    !c----
    !c...  Executable code

    !c Function to solve Kepler's eqn for E (here called
    !c x) for given e and M. returns value of x.
    !c MAY 21 : FOR e < 0.18 use ESOLMD for speed and sufficient accuracy
    !c MAY 21 : FOR e > 0.8 use EHIE - this one may not converge fast enough.

      call orbel_scget(m,sm,cm)

    !c  begin with a guess accurate to order ecc**3
      x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

    !c  Go through one iteration for improved estimate
      call orbel_scget(x,sx,cx)
      es = e*sx
      ec = e*cx
      f = x - es  - m
      fp = 1.d0 - ec
      fpp = es
      fppp = ec
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
      orbel_eget = x + dx

!    c Do another iteration.
!    c For m between 0 and 2*pi this seems to be enough to
!    c get near roundoff error for eccentricities between 0 and 0.8

      x = orbel_eget
      call orbel_scget(x,sx,cx)
      es = e*sx
      ec = e*cx
      f = x - es  - m
      fp = 1.d0 - ec
      fpp = es
      fppp = ec
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

      orbel_eget = x + dx

      return
end  ! orbel_eget

!c---------------------------------------------------------------------
!***********************************************************************
!c                    ORBEL_EHIE.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                           m ==> mean anomaly. (real scalar)
!!*             Returns:
!!*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: Use Danby's quartic for 3 iterations.
!!*                Eqn. is f(x) = x - e*sin(x+M). Note  that
!!*             E = x + M. First guess is very good for e near 1.
!!*             Need to first get M between 0. and PI and use
!!*         symmetry to return right answer if M between PI and 2PI
!!*     REMARKS: Modifies M so that both E and M are in range (0,TWOPI)
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 25,1992.
!!*     REVISIONS:
!***********************************************************************

real*8 function orbel_ehie(e,m)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,m

    !c...  Internals:
      integer iflag,nper,niter,NMAX
      real*8 dx,x,sa,ca,esa,eca,f,fp

      parameter (NMAX = 3)

    !c----
    !c...  Executable code

    !c In this section, bring M into the range (0,TWOPI) and if
    !c the result is greater than PI, solve for (TWOPI - M).
      iflag = 0
      nper = m/TWOPI
      m = m - nper*TWOPI
      if (m .lt. 0.d0) m = m + TWOPI

      if (m.gt.PI) then
       m = TWOPI - m
       iflag = 1
      endif

    !c Make a first guess that works well for e near 1.
      x = (6.d0*m)**(1.d0/3.d0) - m
      niter =0

    !c Iteration loop
      do niter =1,NMAX
        call orbel_scget(x + m,sa,ca)
        esa = e*sa
        eca = e*ca
        f = x - esa
        fp = 1.d0 -eca
        dx = -f/fp
        dx = -f/(fp + 0.5d0*dx*esa)
        dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx))
        x = x + dx
      enddo

      orbel_ehie = m + x

      if (iflag.eq.1) then
      orbel_ehie = TWOPI - orbel_ehie
      m = TWOPI - m
      endif

      return
end         !orbel_ehie


!c------------------------------------------------------------------

!***********************************************************************
!c                    ORBEL_EHYBRID.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                           m ==> mean anomaly. (real scalar)
!!*             Returns:
!!*              orbel_ehybrid ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: For e < 0.18 uses fast routine ESOLMD
!!*             For larger e but less than 0.8, uses EGET
!!*             For e > 0.8 uses EHIE
!!*     REMARKS: Only EHIE brings M and E into range (0,TWOPI)
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 25,1992.
!!*     REVISIONS: 2/26/93 hfl
!***********************************************************************

real*8 function orbel_ehybrid(e,m)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,m

    !c...  Internals:
      real*8 orbel_esolmd,orbel_eget,orbel_ehie

    !c----
    !c...  Executable code

      if(e .lt. 0.18d0) then
      orbel_ehybrid = orbel_esolmd(e,m)
      else
      if( e .le. 0.8d0) then
         orbel_ehybrid = orbel_eget(e,m)
      else
         orbel_ehybrid = orbel_ehie(e,m)
      endif
      endif

      return
end     ! orbel_ehybrid

!c--------------------------------------------------------------------

!***********************************************************************
!c                    ORBEL_ESOLMD.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn.   e is ecc.   m is mean anomaly.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                           m ==> mean anomaly. (real scalar)
!!*             Returns:
!!*                orbel_esolmd ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: Some sort of quartic convergence from Wisdom.
!!*     REMARKS: ONLY GOOD FOR SMALL ECCENTRICITY SINCE IT ONLY
!!*         ITERATES ONCE. (GOOD FOR PLANET CALCS.)
!!*            ALSO DOES NOT PUT M OR E BETWEEN 0. AND 2*PI
!!*     INCLUDES: needs SCGET.F
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 7, 1992.
!!*     REVISIONS: 2/26/93 hfl
!***********************************************************************

real*8 function orbel_esolmd(e,m)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,m

    !c...  Internals:
      real*8 x,sm,cm,sx,cx
      real*8 es,ec,f,fp,fpp,fppp,dx

    !c----
    !c...  Executable code

    !c...    Function to solve Kepler's eqn for E (here called
    !c...    x) for given e and M. returns value of x.

      call orbel_scget(m,sm,cm)
      x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

      call orbel_scget(x,sx,cx)
      es = e*sx
      ec = e*cx
      f = x - es  - m
      fp = 1.d0 - ec
      fpp = es
      fppp = ec
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

      orbel_esolmd = x + dx

      return   ! orbel_esolmd
end

!c--------------------------------------------------------------------
!***********************************************************************
!c                    ORBEL_FGET.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                        capn ==> hyperbola mean anomaly. (real scalar)
!!*             Returns:
!!*                  orbel_fget ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!!*           Cel. Mech. ".  Quartic convergence from Danby's book.
!!*     REMARKS:
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 11, 1992.
!!*     REVISIONS: 2/26/93 hfl
!***********************************************************************

real*8 function orbel_fget(e,capn)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,capn

    !c...  Internals:
      integer i,IMAX
      real*8 tmp,x,shx,chx
      real*8 esh,ech,f,fp,fpp,fppp,dx
      PARAMETER (IMAX = 10)

    !c----
    !c...  Executable code

    !c Function to solve "Kepler's eqn" for F (here called
    !c x) for given e and CAPN.

!c  begin with a guess proposed by Danby
      if( capn .lt. 0.d0) then
       tmp = -2.d0*capn/e + 1.8d0
       x = -log(tmp)
      else
       tmp = +2.d0*capn/e + 1.8d0
       x = log( tmp)
      endif

      orbel_fget = x

      do i = 1,IMAX
      call orbel_schget(x,shx,chx)
      esh = e*shx
      ech = e*chx
      f = esh - x - capn
    !c      write(6,*) 'i,x,f : ',i,x,f
      fp = ech - 1.d0
      fpp = esh
      fppp = ech
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
      orbel_fget = x + dx
    !c   If we have converged here there's no point in going on
      if(abs(dx) .le. TINY) RETURN
      x = orbel_fget
      enddo

      write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE'
      return
      end   ! orbel_fget


!c------------------------------------------------------------------
!***********************************************************************
!c                    ORBEL_FHYBRID.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                           n ==> hyperbola mean anomaly. (real scalar)
!!*             Returns:
!!*               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON
!!*             For larger N, uses FGET
!!*     REMARKS:
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 26,1992.
!!*     REVISIONS:
!!*     REVISIONS: 2/26/93 hfl
!***********************************************************************

real*8 function orbel_fhybrid(e,n)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,n

    !c...  Internals:
      real*8 abn
      real*8 orbel_flon,orbel_fget

    !c----
    !c...  Executable code

      abn = n
      if(n.lt.0.d0) abn = -abn

      if(abn .lt. 0.636d0*e -0.6d0) then
      orbel_fhybrid = orbel_flon(e,n)
      else
      orbel_fhybrid = orbel_fget(e,n)
      endif

      return
end  ! orbel_fhybrid

!c-------------------------------------------------------------------

!***********************************************************************
!c                    ORBEL_FLON.F
!***********************************************************************
!!*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
!*
!!*             Input:
!!*                           e ==> eccentricity anomaly. (real scalar)
!!*                        capn ==> hyperbola mean anomaly. (real scalar)
!!*             Returns:
!!*                  orbel_flon ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: Uses power series for N in terms of F and Newton,s method
!!*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 26, 1992.
!!*     REVISIONS:
!***********************************************************************

real*8 function orbel_flon(e,capn)

      include './swift.inc'

    !c...  Inputs Only:
      real*8 e,capn

    !c...  Internals:
      integer iflag,i,IMAX
      real*8 a,b,sq,biga,bigb
      real*8 x,x2
      real*8 f,fp,dx
      real*8 diff
      real*8 a0,a1,a3,a5,a7,a9,a11
      real*8 b1,b3,b5,b7,b9,b11
      PARAMETER (IMAX = 10)
      PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
      PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
      PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
      PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

    !c----
    !c...  Executable code


    !c Function to solve "Kepler's eqn" for F (here called
    !c x) for given e and CAPN. Only good for smallish CAPN

      iflag = 0
      if( capn .lt. 0.d0) then
       iflag = 1
       capn = -capn
      endif

      a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
      a0 = -6227020800.d0*capn/e
      b1 = a1

    !c  Set iflag nonzero if capn < 0., in which case solve for -capn
    !c  and change the sign of the final answer for F.
    !c  Begin with a reasonable guess based on solving the cubic for small F


      a = 6.d0*(e-1.d0)/e
      b = -6.d0*capn/e
      sq = sqrt(0.25*b*b +a*a*a/27.d0)
      biga = (-0.5*b + sq)**0.3333333333333333d0
      bigb = -(+0.5*b + sq)**0.3333333333333333d0
      x = biga + bigb
    !c    write(6,*) 'cubic = ',x**3 +a*x +b
      orbel_flon = x
    !c If capn is tiny (or zero) no need to go further than cubic even for
    !c e =1.
      if( capn .lt. TINY) go to 100

      do i = 1,IMAX
      x2 = x*x
      f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
      fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))
      dx = -f/fp
    !c      write(6,*) 'i,dx,x,f : '
    !c      write(6,432) i,dx,x,f
    432      format(1x,i3,3(2x,1p1e22.15))
      orbel_flon = x + dx
    !c   If we have converged here there's no point in going on
      if(abs(dx) .le. TINY) go to 100
      x = orbel_flon
      enddo

    !c Abnormal return here - we've gone thru the loop
    !c IMAX times without convergence
      if(iflag .eq. 1) then
       orbel_flon = -orbel_flon
       capn = -capn
      endif
      write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE'
      diff = e*sinh(orbel_flon) - orbel_flon - capn
      write(6,*) 'N, F, ecc*sinh(F) - F - N : '
      write(6,*) capn,orbel_flon,diff
      return

    !c  Normal return here, but check if capn was originally negative
    100   if(iflag .eq. 1) then
       orbel_flon = -orbel_flon
       capn = -capn
      endif

      return
end     ! orbel_flon

!    c------------------------------------------------------------------
!    ***********************************************************************
!    c                      ORBEL_SCGET.F
!    ***********************************************************************
!    !*     PURPOSE:  Given an angle, efficiently compute sin and cos.
!    *
!    !*        Input:
!    !*             angle ==> angle in radians (real scalar)
!    !*
!    !*        Output:
!    !*             sx    ==>  sin(angle)  (real scalar)
!    !*             cx    ==>  cos(angle)  (real scalar)
!    *
!    !*     ALGORITHM: Obvious from the code
!    !*     REMARKS: The HP 700 series won't return correct answers for sin
!    !*       and cos if the angle is bigger than 3e7. We first reduce it
!    !*       to the range [0,2pi) and use the sqrt rather than cos (it's faster)
!    !*       BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
!    !*     AUTHOR:  M. Duncan.
!    !*     DATE WRITTEN:  May 6, 1992.
!    !*     REVISIONS:
!    ***********************************************************************

subroutine orbel_scget(angle,sx,cx)

      include './swift.inc'

  !c...  Inputs Only:
      real*8 angle

  !c...  Output:
      real*8 sx,cx

  !c... Internals:
      integer nper
      real*8 x
      real*8 PI3BY2
      parameter(PI3BY2 = 1.5d0*PI)

  !c----
  !c...  Executable code

      nper = angle/TWOPI
      x = angle - nper*TWOPI
      if(x.lt.0.d0) then
           x = x + TWOPI
      endif
      sx = sin(x)
      cx= sqrt(1.d0 - sx*sx)
      if( (x .gt. PIBY2) .and. (x .lt.PI3BY2)) then
           cx = -cx
      endif

      return
end   ! orbel_scget
  !c-------------------------------------------------------------------
!***********************************************************************
!c                      ORBEL_SCHGET.F
!***********************************************************************
!*     PURPOSE:  Given an angle, efficiently compute sinh and cosh.
!*
!*        Input:
!*             angle ==> angle in radians (real scalar)
!*
!*        Output:
!*             shx    ==>  sinh(angle)  (real scalar)
!*             chx    ==>  cosh(angle)  (real scalar)
!*
!*     ALGORITHM: Obvious from the code
!*     REMARKS: Based on the routine SCGET for sine's and cosine's.
!*       We use the sqrt rather than cosh (it's faster)
!*       BE SURE THE ANGLE IS IN RADIANS AND IT CAN'T BE LARGER THAN 300
!*       OR OVERFLOWS WILL OCCUR!
!*     AUTHOR:  M. Duncan.
!*     DATE WRITTEN:  May 6, 1992.
!*     REVISIONS:
!***********************************************************************

subroutine orbel_schget(angle,shx,chx)

      include './swift.inc'

  !c...  Inputs Only:
      real*8 angle

  !c...  Output:
      real*8 shx,chx

  !c----
  !c...  Executable code

      shx = sinh(angle)
      chx= sqrt(1.d0 + shx*shx)

      return
end   ! orbel_schget
!  !c---------------------------------------------------------------------
!***********************************************************************
!c                    ORBEL_ZGET.F
!***********************************************************************
!!*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola
!!*          given Q (Fitz. notation.)
!*
!!*             Input:
!!*                           q ==>  parabola mean anomaly. (real scalar)
!!*             Returns:
!!*                  orbel_zget ==>  eccentric anomaly. (real scalar)
!*
!!*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!!*     REMARKS: For a parabola we can solve analytically.
!!*     AUTHOR: M. Duncan
!!*     DATE WRITTEN: May 11, 1992.
!!*     REVISIONS: May 27 - corrected it for negative Q and use power
!!*          series for small Q.
!***********************************************************************

real*8 function orbel_zget(q)

      include './swift.inc'

  !c...  Inputs Only:
      real*8 q

  !c...  Internals:
      integer iflag
      real*8 x,tmp

  !c----
  !c...  Executable code

      iflag = 0
      if(q.lt.0.d0) then
      iflag = 1
      q = -q
      endif

      if (q.lt.1.d-3) then
       orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
      else
       x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
       tmp = x**(1.d0/3.d0)
       orbel_zget = tmp - 1.d0/tmp
      endif

      if(iflag .eq.1) then
       orbel_zget = -orbel_zget
       q = -q
      endif

      return
end    ! orbel_zget
  !c----------------------------------------------------------------------

