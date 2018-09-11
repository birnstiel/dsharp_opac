MODULE fit_module
  IMPLICIT NONE
  !
  ! some constants
  !
  doubleprecision :: mu,m_p,k_b,pi,sig_h2,Grav       ! some mutual constants (see below)
  PARAMETER(pi       = 3.141593)                     ! PI
  PARAMETER(k_b      = 1.380658d-16)                 ! Boltzmann constant in erg/K
  PARAMETER(mu       = 2.3d0)                        ! mean molecular mass in proton masses
  PARAMETER(m_p      = 1.6726231d-24)                ! proton mass in g
  PARAMETER(sig_h2   = 2e-15)                        ! cross section of H2
  PARAMETER(Grav     = 6.67259e-8)                   ! gravitational constant in cm^3 g^-1 s^-2

CONTAINS


  ! __________________________________________________________________________________________________________________________________
  ! A new function that gives the velocities according to Ormel and Cuzzi (2007)
  !
  ! INPUT:   tau_1       =   stopping time 1
  !          tau_2       =   stopping time 2
  !          t0          =   large eddy turnover time
  !          v0          =   large eddy velocity
  !          ts          =   small eddy turnover time
  !          vs          =   small eddy velocity
  !          Reynolds    =   Reynolds number
  !
  ! RETURNS: v_rel_ormel =   relative velocity SQUARED
  !
  ! __________________________________________________________________________________________________________________________________
  doubleprecision FUNCTION v_rel_ormel(tau_1, tau_2,t0,v0,ts,vs,Reynolds)
    IMPLICIT NONE
    doubleprecision, INTENT(in)            :: tau_1, tau_2
    doubleprecision, INTENT(in)            :: t0,v0,ts,vs,Reynolds
    doubleprecision                        :: St1, St2, tau_mx, tau_mn, Vg2
    doubleprecision                        :: c0,c1,c2,c3, y_star, ya, eps
    doubleprecision                        :: hulp1, hulp2
    !
    ! sort tau's 1--> correspond to the max. now
    ! (a bit confusing perhaps)
    !
    IF (tau_1 .GE. tau_2) THEN
       tau_mx = tau_1
       tau_mn = tau_2
       St1 = tau_mx/t0
       St2 = tau_mn/t0
    ELSE
       tau_mx = tau_2
       tau_mn = tau_1
       St1 = tau_mx/t0
       St2 = tau_mn/t0
    ENDIF
    !
    ! note the square
    !
    Vg2 = 1.5 *v0**2.0
    !
    ! approximate solution for St*=y*St1; valid for St1 << 1.
    !
    ya = 1.6
    IF (tau_mx .LT. 0.2*ts) THEN
       !
       ! very small regime
       !
       v_rel_ormel = 1.5 *(vs/ts *(tau_mx - tau_mn))**2.0
    ELSEIF (tau_mx .LT. ts/ya) THEN
       v_rel_ormel = Vg2 *(St1-St2)/(St1+St2)*(St1**2.0/(St1+Reynolds**(-0.5)) - St2**2.0/(St2+Reynolds**(-0.5)))
    ELSEIF (tau_mx .LT. 5.0*ts) THEN
       !
       ! Eq. 17 of OC07. The second term with St_i**2.0 is negligible (assuming !Re>>1)
       ! hulp1 = Eq. 17
       ! hulp2 = Eq. 18
       !
       hulp1 = ( (St1-St2)/(St1+St2) * (St1**2.0/(St1+ya*St1) - St2**2.0/(St2+ya*St1)) )!note the -sign
       hulp2 = 2.0*(ya*St1-Reynolds**(-0.5)) + St1**2.0/(ya*St1+St1) - St1**2.0/(St1+Reynolds**(-0.5)) +&
            St2**2.0/(ya*St1+St2) - St2**2.0/(St2+Reynolds**(-0.5))
       v_rel_ormel = Vg2 *(hulp1 + hulp2)
    ELSEIF (tau_mx .LT. t0/5.0) THEN
       !
       ! stopping time ratio
       !
       eps=St2/St1
       !
       ! Full intermediate regime
       !
       v_rel_ormel = Vg2 *( St1*(2.0*ya - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+ya) + eps**3.0/(ya+eps) )) )
    ELSEIF (tau_mx .LT. t0) THEN
       !
       ! now y* lies between 1.6 (St1 << 1) and 1.0 (St1>=1).
       ! The fit below fits ystar to less than 1%
       !
       c3 =-0.29847604
       c2 = 0.32938936
       c1 =-0.63119577
       c0 = 1.6015125
       y_star = c0 + c1*St1 + c2*St1**2.0 + c3*St1**3.0
       !
       !stopping time ratio
       !
       eps=St2/St1
       !
       ! we can then employ the same formula as before
       !
       v_rel_ormel = Vg2 *( St1*(2.0*y_star - (1.0+eps) + 2.0/(1.0+eps) *(1.0/(1.0+y_star) + eps**3.0/(y_star+eps) )) )
    ELSE
       !
       ! heavy particle limit
       !
       v_rel_ormel = Vg2 *( 1.0/(1.0+St1) + 1.0/(1.0+St2) )
    ENDIF
  END FUNCTION v_rel_ormel
  ! ==================================================================================================================================

  ! __________________________________________________________________________________________________________________________________
  ! This function calculates the boundaries of the different regimes in the
  ! grain size distribution.
  !
  ! INPUT:
  !   nm          = # of size/mass points
  !   T           = temperature               [K]
  !   alpha       = turbulence parameter      []
  !   sigma_g     = gas surface density       [g/cm^2]
  !   rho_s       = dust grain volume density [g/cm^3]
  !   m_grid      = mass array                [g]
  !   a_grid      = grain size array          [cm]
  !   m_star      = stellar mass              [g]
  !   R           = radial distance to star   [cm]
  !   v_frag      = fragmentation velocity    [cm/s]
  !
  ! OUTPUT:
  ! position (in grain size) of the transisions:
  !   a_01    = first transition between brownian motion and turbulence
  !   a_12    = linear to intermediate regime of turbulende
  !   a_left  = left side of cratering peak
  !   a_peak  = position of cratering peak
  !   a_right = right side of cratering peak
  ! __________________________________________________________________________________________________________________________________
  SUBROUTINE get_boundaries(nm,T,alpha,sigma_g,rho_s,m_grid,a_grid,m_star,R,v_frag,a_01,a_12,a_left,a_peak,a_right,a_sett)
    IMPLICIT NONE
    INTEGER,         INTENT(in)  :: nm
    doubleprecision, INTENT(in)  :: T,alpha,sigma_g,rho_s,m_star,R,v_frag
    doubleprecision, INTENT(in)  :: m_grid(1:nm),a_grid(1:nm)
    doubleprecision, INTENT(out) :: a_01,a_12,a_left,a_peak,a_right,a_sett
    !
    ! used for polynomial solving
    !
    !    doubleprecision,external  :: polyn
    !    integer                   :: ierr
    !
    ! other
    !
    doubleprecision :: cs,Re,a1,ya,omega,ro,tn,ts,vn,vs,tau_1,tau_2,xL,xR,yL,yR
    doubleprecision :: dv_BM(1:nm,1:nm),dv_TM(1:nm,1:nm),dv_ii(1:nm),dv_i1(1:nm)
    INTEGER         :: i,j,i_peak,i_left,i_right
    !
    ! initialize output
    !
    a_01    = 0d0
    a_12    = 0d0
    a_left  = 0d0
    a_peak  = 0d0
    a_right = 0d0
    a_sett  = 0d0
    !
    ! sound speed, reynolds number and grainsize(1)
    !
    cs     = SQRT(k_b*T/mu/m_p)
    Re     = alpha*sig_h2*sigma_g/(2d0*mu*m_p)
    a1     = a_grid(1)
    !
    ! get the size where particles start to settle
    !
    a_sett = 2d0*alpha*sigma_g/(pi*rho_s)
    !
    ! get the transition from BM to turbulence
    !
    ! OLD VERSION
    !    A          = 0.1d0 * 32d0*k_b*T*sigma_g**2d0/(pi**4d0*rho_s**3d0*sqrt(Re)*alpha*cs**2d0)
    !    poly_coeff = (/1d0,-2d0*a1,a1**2d0,0d0,0d0,-A/)
    !    !
    !    ! solve polynomial
    !    !
    !    a_01 = zbrent(polyn,a_grid(1),a_grid(nm),1d-10,ierr)
    !    !
    !    ! if ZBRENT fails, then we give up
    !    !
    !    if (ierr/=0) then
    !        write(*, *) 'ERROR:   Failure by ZBRENT, exiting'
    !        stop 92378
    !    endif
    ! NEW APPROXIMATION
    a_01 = 0.7943d0 * ( 8d0*sigma_g/pi/rho_s*Re**(-0.25d0)*SQRT(mu*m_p/(3d0*pi*alpha))*(4d0/3d0*pi*rho_s)**(-0.5d0) )**(2d0/5d0)
    !
    ! get the turbulence bump
    !
    ya     = 1.6d0
    a_12   = 1d0/(ya*pi*rho_s)*SQRT(8d0*mu*m_p*sigma_g/(alpha*sig_h2))
    !
    ! get the cratering bump positions
    !
    !
    ! calculate relative velocities due to brownian motion
    !
    DO i=1,nm
       DO j=1,nm
          dv_BM(i,j) = SQRT(8d0*k_b*T*(m_grid(i)+m_grid(j))/pi/m_grid(i)/m_grid(j))
       ENDDO
    ENDDO
    !
    ! calculate relative velocities due to turbulence
    !
    DO i=1,nm
       DO j=1,nm
          omega = SQRT(Grav*m_star/R**3d0)
          ro    = sigma_g*omega/(SQRT(2*pi)*cs)
          tn    = 1d0/omega
          ts    = tn*Re**(-0.5d0)
          vn    = SQRT(alpha)*cs
          vs    = vn*Re**(-0.25d0)

          tau_1 = rho_s*a_grid(i)/sigma_g/omega*pi/2d0
          tau_2 = rho_s*a_grid(j)/sigma_g/omega*pi/2d0

          dv_TM(i,j) = v_rel_ormel(tau_1,tau_2,tn,vn,ts,vs,Re)
       ENDDO
    ENDDO
    !
    ! put the velocities together,
    ! note that the turbulent ones are already squared
    !
    DO i = 1,nm
       dv_ii(i) = SQRT(dv_BM(i,i)**2d0 + dv_TM(i,i))
       dv_i1(i) = SQRT(dv_BM(i,1)**2d0 + dv_TM(i,1))
    ENDDO
    !
    ! the position of the peak
    !
    IF (MAXVAL(dv_ii)<V_FRAG) THEN
       WRITE(*,*) 'WARNING:   particles will not fragment!'
       WRITE(*,*) '           Either increase turbulent velocities or decrease'
       WRITE(*,*) '           the fragmentation velocity'
       !stop 38332
    ENDIF
    i_peak = first(nm,(dv_ii-0.8d0*v_frag >=0d0))
    IF (i_peak<2) THEN
       i_peak = first(nm,(dv_ii-0.8d0*v_frag>=0d0).AND.(a_grid>1d-4))
       IF (i_peak<2) THEN
          WRITE (*,*) 'WARNING:   cannot find i_peak, putting everything in smallest bin'
          !stop 12783
          a_peak = a_grid(1)
       ELSE
          a_peak = a_grid(i_peak)
       ENDIF
    ELSE
       xL = a_grid(i_peak-1)
       xR = a_grid(i_peak)
       yL = dv_ii(i_peak-1)-0.8d0*v_frag
       yR = dv_ii(i_peak)  -0.8d0*v_frag
       a_peak = xL-yL*(xR-xL)/(yR-yL)
    ENDIF
    !
    ! the position of the left wing
    !
    i_left = first(nm,dv_i1-0.8d0*v_frag>=0)
    IF (i_left<2) THEN
       i_left    = first(nm,(dv_i1-0.8d0*v_frag>=0d0).AND.(a_grid>1d-4))
       IF (i_left<2) THEN
          WRITE(*,*) 'WARNING:   cannot find i_left, putting everything in smallest bin'
          !stop 23212
       ELSE
          a_left = a_grid(i_left)
       ENDIF
    ELSE
       xL = a_grid(i_left-1)
       xR = a_grid(i_left)
       yL = dv_i1(i_left-1)-0.8d0*v_frag
       yR = dv_i1(i_left)  -0.8d0*v_frag
       a_left = xL-yL*(xR-xL)/(yR-yL)
    ENDIF
    !
    ! the position of the right wing
    !
    i_right = first(nm,dv_ii-v_frag>=0d0)
    IF (i_right<2) THEN
       i_right = first(nm,(dv_ii-v_frag>=0d0).AND.(a_grid>1d-4))
       IF (i_right<2) THEN
          WRITE(*,*) 'WARNING:   cannot find i_right, putting everything in smallest bin'
          !stop 27431
          a_right = a_grid(1)
       ELSE
          a_right = a_grid(i_right)
       ENDIF
    ELSE
       xL = a_grid(i_right-1)
       xR = a_grid(i_right)
       yL = dv_ii(i_right-1)-v_frag
       yR = dv_ii(i_right)  -v_frag
       a_right = xL-yL*(xR-xL)/(yR-yL)
    ENDIF

  END SUBROUTINE get_boundaries
  ! ==================================================================================================================================

  ! __________________________________________________________________________________________________________________________________
  !                   FUNCTION: FIND LAST OF MASK
  !
  ! finds the last element where mask is true
  ! example:
  !
  ! A = (/0,1,2,3,4,5,6,7/)
  !
  ! last(size(A),A<2) = 2
  ! __________________________________________________________________________________________________________________________________
  INTEGER FUNCTION last(n,mask)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    LOGICAL,INTENT(in) :: mask(1:n)
    INTEGER :: i

    last = 0
    DO i=n,1,-1
       IF (mask(i)) THEN
          last=i
          RETURN
       ENDIF
    ENDDO

  END FUNCTION last
  ! ==================================================================================================================================

  ! __________________________________________________________________________________________________________________________________
  !                   FUNCTION: FIND FIRST OF MASK
  !
  ! finds the first element where mask is true
  ! example:
  !
  ! A = (/0,1,2,3,4,5,6,7/)
  !
  ! first(size(A),A>2) = 4
  ! __________________________________________________________________________________________________________________________________
  INTEGER FUNCTION first(n,mask)
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    LOGICAL,INTENT(in) :: mask(1:n)
    INTEGER :: i

    first = 0
    DO i=1,n
       IF (mask(i)) THEN
          first=i
          RETURN
       ENDIF
    ENDDO

  END FUNCTION first
  ! ==================================================================================================================================


  ! __________________________________________________________________________________________________________________________________
  ! This function smoothes the given 1D array of size n by a running average over
  ! 3 points. The averaging is performed n_smooth times.
  !
  ! USAGE:   smooth1D(n,n_smooth,arr_in)
  !
  ! INPUT:   n        = size of in/output array (1D)
  !          n_smooth = how much the array is to be smoothed
  !          arr_in   = array which is to be smoothed
  !
  ! OUTPUT:  the input array is smoothed
  ! __________________________________________________________________________________________________________________________________
  SUBROUTINE smooth1D(n,n_smooth,arr_in)
    IMPLICIT NONE
    INTEGER,         INTENT(in)    :: n,n_smooth
    doubleprecision, INTENT(inout) :: arr_in(1:n)
    doubleprecision                :: arr_dummy(1:n)
    INTEGER                        :: i_smooth,i
    !
    ! some checks
    !
    IF (n_smooth<1) THEN
       WRITE(*,*) 'ERROR:   n_smooth has to be larger than 0!'
       STOP 72311
    ENDIF
    IF (n_smooth>=0.5d0*n) THEN
       WRITE(*,*) 'WARNING: n_smooth is too large, this way you will average out'
       WRITE(*,*) '         any trends in your array!'
    ENDIF
    !
    ! average n_smooth times
    !
    arr_dummy = arr_in
    DO i_smooth = 1,n_smooth
       DO i = 1,n
          arr_in(i) = ( arr_dummy(MAX(1,i-1))+arr_dummy(i)+arr_dummy(MIN(i+1,n)) )/3d0
       ENDDO
       arr_dummy = arr_in
    ENDDO

  END SUBROUTINE smooth1D
  ! ==================================================================================================================================

  ! __________________________________________________________________________________________________________________________________
  ! Brent's Algorithm for root finding
  !
  ! USAGE:     zbrent(func,x1,x2,tol,ierr)
  !
  ! INPUT:     func   =    the double precision function to be worked with
  !            x1,x2  =    the interval around the root
  !            tol    =    the tolerance around the root
  !            ierr   =    0 means no error, 1 means error
  !
  ! OUTPUT:    zbrent =    the root of function "func" with a precision of "tol"
  !
  ! __________________________________________________________________________________________________________________________________
  doubleprecision FUNCTION zbrent(func,x1,x2,tol,ierr)
    IMPLICIT NONE
    EXTERNAL func
    doubleprecision             :: func
    doubleprecision, INTENT(in) :: tol,x1,x2
    INTEGER, INTENT(out)        :: ierr
    INTEGER                     :: iter
    doubleprecision             :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    INTEGER,PARAMETER         :: ITMAX=100
    doubleprecision,PARAMETER :: EPSS=3.d-8
    ierr = 0
    a    = x1
    b    = x2
    fa   = func(a)
    fb   = func(b)

    IF((fa.GT.0..AND.fb.GT.0.).OR.(fa.LT.0..AND.fb.LT.0.)) THEN
       WRITE(*,*) 'root must be bracketed for zbrent'
       ierr = 1
    ENDIF
    c  = b
    fc = fb
    DO iter=1,ITMAX
       IF((fb.GT.0..AND.fc.GT.0.).OR.(fb.LT.0..AND.fc.LT.0.))THEN
          c  = a
          fc = fa
          d  = b-a
          e  = d
       ENDIF
       IF(ABS(fc).LT.ABS(fb)) THEN
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       ENDIF
       tol1 = 2d0*EPSS*ABS(b)+0.5d0*tol
       xm   = 0.5d0*(c-b)
       IF(ABS(xm)<=tol1 .OR. fb<1e-12)THEN
          zbrent = b
          RETURN
       ENDIF
       IF (ABS(e)>=tol1 .AND. ABS(fa)>ABS(fb)) THEN
          s = fb/fa
          IF(ABS(a/c-1.0)<1e-8) THEN
             p = 2.d0*xm*s
             q = 1.d0-s
          ELSE
             q = fa/fc
             r = fb/fc
             p = s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
             q = (q-1.)*(r-1.)*(s-1.)
          ENDIF
          IF(p>0.) q = -q
          p = ABS(p)
          IF(2.d0 * p < MIN(3.*xm*q-ABS(tol1*q),ABS(e*q))) THEN
             e = d
             d = p/q
          ELSE
             d = xm
             e = d
          ENDIF
       ELSE
          d = xm
          e = d
       ENDIF
       a  = b
       fa = fb
       IF (ABS(d) > tol1) THEN
          b = b+d
       ELSE
          b = b + SIGN(tol1,xm)
       ENDIF
       fb = func(b)
    ENDDO
    WRITE(*,*) 'WARNING: zbrent exceeding maximum iterations'
    zbrent=b
    ierr = 2
    RETURN
  END FUNCTION zbrent
  ! ==================================================================================================================================


    ! __________________________________________________________________________________________________________________________________
    ! This subroutine provides a fit function to the grain size distribution of the
    ! given parameters BUT ONLY FOR A XI OF (about) 1.8!!
    !
    ! USAGE: fit_function18(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid,m_star,R,v_frag)
    !
    ! INPUT:    nm      = number of mass grid points   []
    !           xi      = fragmentation power law      []
    !           T       = mid-plane temperature        [K]
    !           alpha   = turbulence alpha             []
    !           sigma_g = gas surface density          [g cm^-2]
    !           sigma_d = dust surface density         [g cm^-2]
    !           rho_s   = solid density of the grains  [g cm^-3]
    !           m_grid  = mass grid                    [g]
    !           a_grid  = size grid                    [cm]
    !           m_star  = stellar mass                 [g]
    !           R       = mid-plane distance to star   [cm]
    !           v_frag  = fragmentation velocity       [cm s^-1]
    !
    ! OUTPUT:   fit     = fit-distribution in g cm^-2 in each mass bin
    !           a_01    = first regime boundary
    !           a_12    = second regime boundary
    !           a_l     = left wing of peak
    !           a_p     = center of peak
    !           a_r     = right wing of peak
    !           a_sett  = where settling starts to play a role
    ! __________________________________________________________________________________________________________________________________
    SUBROUTINE fit_function18_test(fit,a_01,a_12,a_l,a_p,a_r,a_sett,nm,xi,T,alpha,sigma_g,sigma_d,rho_s,m_grid,a_grid, &
        & m_star,R,v_frag)
      IMPLICIT NONE
      INTEGER,         INTENT(in)  :: nm
      doubleprecision, INTENT(in)  :: xi,T,alpha,sigma_g,sigma_d,rho_s,m_star,R,v_frag
      doubleprecision, INTENT(in)  :: m_grid(1:nm),a_grid(1:nm)
      doubleprecision, INTENT(out) :: fit(1:nm),a_01,a_12,a_l,a_p,a_r,a_sett
      INTEGER                      :: i_01,i_12,i_l,i_p,i_r,i_sett,i_inc,i,ir(1:4)
      doubleprecision              :: nu,L,N,sig
      doubleprecision              :: alpha_0,alpha_1,alpha_2,slope_0,slope_1,slope_2,dv,jumpfactor,inc_factor
      doubleprecision              :: alpha_0_s,alpha_1_s,alpha_2_s,slope_0_s,slope_1_s,slope_2_s
      doubleprecision              :: slopes_s(1:3),slopes(1:3)
      doubleprecision              :: bump(1:nm)
      !
      ! get the sizes of the boundaries and the according indices
      !
      CALL get_boundaries(nm,T,alpha,sigma_g,rho_s,m_grid,a_grid,m_star,R,v_frag,a_01,a_12,a_l,a_p,a_r,a_sett)
      i_01   = first(nm,a_grid>=a_01)
      i_12   = first(nm,a_grid>=a_12)
      i_l    = first(nm,a_grid>=a_l)
      i_p    = first(nm,a_grid>=a_p)
      i_r    = first(nm,a_grid>=a_r)
      i_sett = first(nm,a_grid>=a_sett)
      !
      ! some consistency checks
      !
      IF ((i_01==0).OR.(i_12==0).OR.(i_l==0).OR.(i_p==0).OR.(i_r==0)) THEN
         WRITE(*,*) 'ERROR:   indices are not correct!'
         WRITE(*,*) 'i_01   = ',i_01
         WRITE(*,*) 'i_12   = ',i_12
         WRITE(*,*) 'i_l    = ',i_l
         WRITE(*,*) 'i_p    = ',i_p
         WRITE(*,*) 'i_r    = ',i_r
         WRITE(*,*) 'i_sett = ',i_sett
         STOP 18072
      ENDIF
      IF (LOG10(MAXVAL((/a_l,a_p,a_r/)))-LOG10(a_12)<0.1d0) THEN
         WRITE(*,*) 'WARNING: maximum grain size is very small,'
         WRITE(*,*) '         fit results are not reliable!'
      ENDIF
      IF (i_p==1) THEN
         fit    = 0d0
         fit(1) = sigma_d
      ELSE
         !
         ! define the slopes without settling
         !
         nu = 1d0/6d0
         alpha_0 = (nu+xi+1d0)/2d0
         slope_0 = (6d0-3d0*alpha_0)

         nu = 1d0
         alpha_1 = (nu+xi+1d0)/2d0
         slope_1 = (6d0-3d0*alpha_1)

         nu = 5d0/6d0
         alpha_2 = (nu+xi+1d0)/2d0
         slope_2 = (6d0-3d0*alpha_2)
         !
         ! define the slopes with settling
         !
         nu        = 2d0/6d0
         alpha_0_s = (nu+xi+1d0)/2d0
         slope_0_s = 6d0-3d0*alpha_0_s

         nu        = 7d0/6d0
         alpha_1_s = (nu+xi+1d0)/2d0
         slope_1_s = 6d0-3d0*alpha_1_s

         nu        = 1d0
         alpha_2_s = (nu+xi+1d0)/2d0
         slope_2_s = 6d0-3d0*alpha_2_s
         !
         ! plug them together to an array each
         !
         slopes   = (/ slope_0,  slope_1,  slope_2   /)
         slopes_s = (/ slope_0_s,slope_1_s,slope_2_s /)
         !
         ! construct the fit function
         !
         ! EXPERIMENTAL: the jump factor is something between 1 and 2.5 it
         ! seems.
         !
         !dv = sqrt(3d0/2d0*alpha)*sqrt(3d0)*sqrt(k_b*T/mu/m_p)*sqrt(1d0/(sigma_g*2d0*1.6d0))*(8d0*mu*m_p*sigma_g/alpha/sig_h2)**0.25d0
         dv = SQRT(3d0/2d0*alpha)*SQRT(k_b*T/mu/m_p)*SQRT(1d0/(sigma_g*2d0*1.6d0))*(8d0*mu*m_p*sigma_g/alpha/sig_h2)**0.25d0
         !
         ! make the powerlaws with settling
         !
         ir = (/1,i_01,i_12,i_p/)
         fit    = 0d0
         fit(1) = 1
         DO i = 1,3
            IF (i==3) THEN
               ! old one:
               !jumpfactor = max(1,1.2+log(dv/30));
               ! this one works good for V_FRAG = 1 m/s
               !jumpfactor = min(2.5,max(1.1,1+dv/40));
               ! lets try to scale it with V_FRAG
               !jumpfactor = min(2.5d0,max(1.1d0,1d0 + dv/50d0 * 1d0/(v_frag/100d0)))
               ! just rewritten
               !jumpfactor = min(2.5d0,max(1.1d0,1d0 + sqrt(3d0)*2d0*dv/v_frag))
               ! now smoothed
               jumpfactor = (2.5d0**(-9d0)+(1.1d0**9d0+(1d0+2d0*SQRT(3d0)*dv/v_frag)**9d0)**(-1d0))**(-1d0/9d0)
            ELSE
               jumpfactor = 1
            ENDIF
            IF (i_sett<=ir(i)) THEN
               fit(ir(i):ir(i+1))  = fit(ir(i))/jumpfactor*(a_grid(ir(i):ir(i+1)) /a_grid(ir(i))) **slopes_s(i)
            ELSEIF ((ir(i)<=i_sett).AND.(i_sett<=ir(i+1))) THEN
               fit(ir(i):i_sett)   = fit(ir(i))/jumpfactor*(a_grid(ir(i):i_sett)  /a_grid(ir(i))) **slopes(i)
               fit(i_sett:ir(i+1)) = fit(i_sett)          *(a_grid(i_sett:ir(i+1))/a_grid(i_sett))**slopes_s(i)
            ELSEIF (i_sett>ir(i+1)) THEN
               fit(ir(i):ir(i+1))  = fit(ir(i))/jumpfactor*(a_grid(ir(i):ir(i+1)) /a_grid(ir(i))) **slopes(i)
            ELSE
               WRITE(*,*) 'ERROR:   something went wrong with the regimes!'
               STOP 19723
            ENDIF
         ENDDO
         !
         ! smooth the fit
         !
         !call smooth1D(nm,3,fit)
         !
         ! EXPERIMENTAL: SLOPE-INCREASE AT UPPER END OF DISTRIBUTION
         ! we will do this if a_12 and a_p are too close to each other
         !
         !
         IF (LOG10(a_p/a_12)<0.7) THEN
            i_inc = last(nm,a_grid<0.3d0*a_grid(i_p))
            inc_factor = 1d0
            IF (i_inc==0) THEN
               WRITE(*,*) 'WARNING: i_inc == 0 => i_p seems to be quite small, consider a smaller m_min'
               i_inc = 1
            ENDIF
            fit(i_inc:i_p) = fit(i_inc:i_p)+inc_factor*fit(i_inc:i_p)*(1-(a_grid(i_inc:i_p)-a_grid(i_p))/ &
            & (a_grid(i_inc)-a_grid(i_p)))
         ENDIF
         !
         ! produce a bump
         !
         ! GAUSSIAN:
         L    = fit(i_l)
         N    = 2d0*fit(i_l)
         sig  = MAX( ABS(a_grid(i_p-1)-a_grid(i_p)) , MIN(ABS(a_grid(i_r)-a_grid(i_p)),ABS(a_grid(i_l)-a_grid(i_p))))/SQRT(LOG(N/L))
         bump = N*EXP(-(a_grid-a_grid(i_p))**2d0/(sig**2d0))
         !
         ! add the bump
         !
         DO i = i_l,i_r
            fit(i) = MAX(fit(i),bump(i))
         ENDDO
         !
         ! normalize the thing
         !
         fit = sigma_d*fit/SUM(fit)
         IF (ANY(isnan(fit))) WRITE(*,*) 'NaN occured in renormalization'
      ENDIF
    END SUBROUTINE fit_function18_test
    ! ==================================================================================================================================


END MODULE fit_module
! ==================================================================================================================================
