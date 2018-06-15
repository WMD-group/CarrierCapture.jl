program nonrad
!===========================================================================
! A program to calculate the temperature-dependence
!           of the non-radiative recombination rate.
!           The electronic e-ph coupling matrix element
!           is assumed to be independent of temperature.
!           This is not a good assumption for capture of
!           free carriers. In this case the result should
!           be scaled according to a separate consideration.
!           1D approximation used
! Input needed:
!           DQ: difference in equilibrium geometrs, in amu^{1/2}*Angstrom
!           DE: energy difference between the states, in eV (ZPL)
!            V: electron-phonon coupling matrix element, in eV/(amu^{1/2}*Angstrom)
!       w1, w2: frequencies in the excited and ground states, in eV
!        sigma: smearing, optimal value ~0.8*w2
!        Omega: unit cell volume
! Output:
!         Rate: in s-1
! Caution should be exercised in choosing the smearing sigma
! Audrius Alkauskas, UCSB 2012 & 2013
! More info: PRB 90, 075202 (2014); PRL 109, 267401 (2012)
!===========================================================================
IMPLICIT NONE

! temperature range [T0 -- (T0 + dT*NT)]
REAL*16, PARAMETER               :: T0=0.0, dT=100.0
INTEGER, PARAMETER               :: NT=10
! NN1, NN2: number of vibrations in the excited and ground state to be included
! avoid numbers NN1+NN2>75, the factorial too large
! to be developed in the future
! INTEGER, PARAMETER               :: NN1=90, NN2=100
INTEGER, PARAMETER               :: NN1=60, NN2=1
INTEGER, PARAMETER               :: NN2_i=0, NN2_f=60
! energy difference in eV
! REAL*16, PARAMETER               :: DE=-0.352
REAL*16, PARAMETER               :: DE=-0.031
! mass-weighted total quadratic displacement in amu^{1/2}*Angstrom
REAL*16, PARAMETER               :: DQ=12.388333304692587
! frequencies in the excited and the ground state in eV
REAL*16, PARAMETER               :: w2=0.00353164678821, w1=0.0106953723964
! electron-phonon coupling in eV/(amu^{1/2}*Angstrom)
REAL*16, PARAMETER               :: V=0.0793538385464
! smearing of the delta-function in Fermi's golden rule in eV
! 0.75*w2 is always good
REAL*16, PARAMETER               :: sigma = 0.0075
! volume of the supercell
REAL*16, PARAMETER               :: Omega =1.28463E-21 !1E-21 !ZnO 1.136E-21 !GaN 1.10E-21 
! conversion factors, physical constants
REAL*16, PARAMETER               :: pi=3.14159265
REAL*16, PARAMETER               :: K2eV=8.61705E-5     ! k_B/e Kelvin to eV
REAL*16, PARAMETER               :: Factor=15.4669104   ! (e*amu)^{1/2}*1E-10/hbar 
REAL*16, PARAMETER               :: Factor2=6.351E12    ! hbar/1E-20/amu
REAL*16, PARAMETER               :: Factor3=1.52E15     !   e/hbar

! variables
INTEGER                          :: i, j, k
REAL*16                          :: TE, T, E, Ex
REAL*16                          :: R1, Rpart1, R2, Rpart2, R3, Rpart3, Rtot
REAL*16                          :: weight, Z, g
REAL*16                          :: sinfi, cosfi
REAL*16                          :: w, rho, S, Tot
! REAL*16, DIMENSION(0:NN1+1,0:NN2) :: Overlap
REAL*16, DIMENSION(0:NN1+1,NN2_i:NN2_f) :: Overlap
REAL*16                          :: control, control1, control2

!===========================================================================

! output files
! temperature dependance of non-radiative recombinations

! partial rate files could be written during debugging
! 2*pi*(V**2)*(<i|Q|f>**2)*delta(E), in s-1 
!OPEN (unit=13, file='./Temp-nonrad1.dat',  status='replace') 
! 2*pi*(V**2)*(<i|f>)**2*DQ**2*delta(E), in s-1
!OPEN (unit=14, file='./Temp-nonrad2.dat', status='replace')       
! 2*pi*(V**2)*(<i|Q+DQ|f>**2)*delta(E), in s-1
!OPEN (unit=15, file='./Temp-nonrad3.dat', status='replace')       

! total
OPEN (unit=20, file='./Temp-nonrad.dat', status='replace')       

! to calculate overlaps
w = w1*w2/(w1+w2)
rho = Factor*SQRT(w/2.0)*DQ
sinfi = SQRT(w2)/SQRT(w1+w2)
cosfi = SQRT(w1)/SQRT(w1+w2)

! overlap integrals
DO i=0, NN1+1
  DO j=NN2_i, NN2_f
    Overlap(i,j)=overl(i,j)
  END DO
END DO

! temperature cycle  
DO k=1,NT
  control = 0.0
  T = T0 + k*dT
  TE = T*K2eV
  ! statistical sum
  Z = 1.0/(1-exp(-w1/TE))
  R1 = 0.0
  R2 = 0.0
  R3 = 0.0
  Rpart1 = 0.0
  Rpart2 = 0.0
  Rpart3 = 0.0
  weight = 1/Z
  control2 = 0.0
  DO j=NN2_i, NN2_f
    g = exp(-(DE-j*w2)**2/(2.0*sigma**2))/(sigma*sqrt(2.0*pi))
    Rpart1 = Rpart1+weight*g*(Overlap(1,j))**2
    Rpart2 = Rpart2+weight*g*(Overlap(0,j))**2
    Rpart3 = Rpart3+weight*g*Overlap(1,j)*Overlap(0,j) ! to be checked
    control2 = control2 + (Overlap(0,j))**2
  END DO
  R1 = R1 + Rpart1
  R2 = R2 + Rpart2
  R3 = R3 + Rpart3
  control = control + weight*control2
  DO i=1,NN1
    control2 = 0.0
    Rpart1 = 0.0
    Rpart2 = 0.0
    Rpart3 = 0.0
    weight = exp(-i*w1/TE)/Z
    DO j=NN2_i, NN2_f
      g = exp(-(DE+i*w1-j*w2)**2/(2*sigma**2))/(sigma*sqrt(2*pi))
      Rpart1 = Rpart1+weight*g*(sqrt(REAL(i+1))*Overlap(i+1,j)+sqrt(REAL(i))*Overlap(i-1,j))**2
      Rpart2 = Rpart2+weight*g*(Overlap(i,j))**2
      Rpart3 = Rpart3+weight*g*Overlap(i,j)*(sqrt(REAL(i+1))*Overlap(i+1,j)+sqrt(REAL(i))*Overlap(i-1,j))
      control2 = control2 + (Overlap(i,j))**2
    END DO
    control = control + weight*control2
    R1 = R1 + Rpart1
    R2 = R2 + Rpart2
    R3 = R3 + Rpart3
  END DO
  R1 =    R1 * 2*pi * V**2 * Omega / (2*w1) * (Factor2) * 1.0
  R2 =    R2 * 2*pi * V**2 * Omega * DQ**2 * (Factor3) * 1.0
  R3 = 2* R3 * 2*pi * V**2 * Omega * DQ * (Factor3 * Factor2)**0.5 / (2*w1)**0.5 * 1.0
  Rtot = R1+R2+R3
  !WRITE(13, '(2e14.6)'), T,  R1*4.0!, R1*150.0/T**0.5 - "scaling factor"
  !WRITE(14, '(2e14.6)'), T,  R2*4.0!, R2*150.0/T**0.5 - "scaling factor"
  !WRITE(15, '(2e14.6)'), T,  R3*4.0!, R3*150.0/T**0.5 - "scaling factor"
  WRITE(20, '(4e14.6)'), T,  1/T, Rtot*1.0, Rtot*1.0 !150.0/T**0.5*4.0
  print *, control

END DO


CONTAINS

!======================================================================
FUNCTION overl(m, n)
! overlap between vibrational wavefunctions
! formulas from Zapol, CPL 93, 549 (1982)

INTEGER, INTENT(IN)  :: m, n
REAL*16              :: overl, Pr1, Pr2, Pr3
INTEGER              :: k, l
INTEGER              :: l1, l2, lx, k1, k2, kx
REAL*16              :: Ix, f, x

! Pr1 = (-1)**(n)*sqrt(2*cosfi*sinfi)*exp(-rho**2)  (-1)**n or (-1)**m??
  Pr1 = (-1)**(m)*sqrt(2*cosfi*sinfi)*exp(-rho**2)
  Ix = 0.0
  k1 = INT(m/2)
  k2 = MOD(m,2)
  l1 = INT(n/2)
  l2 = MOD(n,2)
  DO kx=0,k1
    DO lx=0,l1
      k=2*kx+k2
      l=2*lx+l2
      Pr2 = (fact(n)*fact(m))**(0.5)/(fact(k)*fact(l)*fact(k1-kx)*fact(l1-lx))*2.0**((k+l-m-n)/2)
      Pr3 = (sinfi**k)*(cosfi**l)
      CALL herm(k+l,rho,f)
      Ix = Ix + Pr1*Pr2*Pr3*f
    END DO
  END DO
  overl =  Ix
END FUNCTION overl
!======================================================================


!======================================================================
SUBROUTINE herm(n,x,f)
! Hermite polynomial of order n by recursion formula
! uses "physicist's" definition of the polynomials

INTEGER, INTENT(IN)                      :: n
REAL*16, INTENT(IN)                      :: x
REAL*16, INTENT(OUT)                     :: f
INTEGER                                  :: k
REAL*16                                  :: y0, y1, dy0, dy1
REAL*16                                  :: yn, dyn

  y0=1.0D0
  dy0=0.0D0
  y1=2.0D0*x
  dy1=2.0D0

IF (n==0) THEN
  yn=y0
ELSE IF (n==1) THEN
  yn=y1
ELSE 
  DO k=2,n
    yn = 2.0*x*y1-dy1
    dyn = 2.0*k*y1
    y1=yn
    dy1=dyn
  END DO
END IF

f = yn

END SUBROUTINE herm
!======================================================================


!======================================================================
FUNCTION fact(k)
! factorial of k

REAL*16              :: fact, f
INTEGER, INTENT(IN)  :: k
INTEGER              :: i
  f = 1.0
  DO i=1,k
    f = f*i
  END DO
  fact = f
END FUNCTION fact
!======================================================================

END PROGRAM nonrad



