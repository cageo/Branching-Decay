C     **********************
      FUNCTION DCMPLX(R1,I1)
c     **********************
C     This is a function to calculate a double precision complex number
c     DCMPLX = R1+i*I1
c     Coded by JP Keizer, 2002 05
c
      COMPLEX (KIND=8) :: DCMPLX
      REAL (KIND=8), INTENT(IN) :: R1,I1
c
      DCMPLX=CMPLX(R1,I1,KIND=8)
c
      END FUNCTION DCMPLX
c
C     ******************
      FUNCTION DREAL(I1)
c     ******************
C     This is a function that returns the real part of a
c     double precision complex number I1 = DREAL+i*DIMAG
C     Coded by CJ Neville, 2002 11
c
      REAL (KIND=8) :: DREAL
      COMPLEX (kind=8), INTENT(IN) :: I1
c
      DREAL=REAL(I1)
c
      END FUNCTION DREAL
c
C     ******************
      FUNCTION DIMAG(I1)
c     ******************
C     This is a function that returns the imaginary part of a
c     double precision complex number I1 = A+i*DIMAG
C     Coded by JP Keizer, 2002 05
c
      REAL (KIND=8) :: DIMAG
      COMPLEX (kind=8), INTENT(IN) :: I1
c
      DIMAG=AIMAG(I1)
c
      END FUNCTION DIMAG
c
C     *******************
      FUNCTION CDABS(Z1)
c     *******************
C     This is a function to calculate the absolute
c     value of a double precision complex number Z1 = A + iB
C     Coded by JP Keizer, 2002 11
c
      REAL (KIND=8) :: CDABS
      COMPLEX (KIND=8), INTENT(IN) :: Z1
      REAL (KIND=8) :: A   ! The real part of the complex number Z1
      REAL (KIND=8) :: B   ! The imaginary part of the complex number Z1
c
      A = REAL(Z1,KIND=8)
      B = AIMAG(Z1)

      CDABS=(A*A + B*B)**0.5
c
      END FUNCTION CDABS
c
c     ******************
      FUNCTION CDEXP(Z1)
c     ******************
C     This is a function to calculate the double precision exponential
c     of a complex number Z1=A+iB
c     Coded by JP Keizer, 2002 01
c
      COMPLEX (KIND=8) :: CDEXP
      COMPLEX (KIND=8), INTENT(IN) :: Z1
      REAL (KIND=8) :: A   ! The real part of the complex number Z
      REAL (KIND=8) :: B   ! The imaginary part of the complex number Z
c
      A = REAL(Z1,KIND=8)
      B = AIMAG(Z1)
      CDEXP = DEXP(A)*CMPLX(DCOS(B),DSIN(B),KIND=8)
c
      END FUNCTION CDEXP
c
c     ******************
      FUNCTION CDLOG(Z1)
c     ******************
C     This is a function to calculate the double precision
c     logarithm of a complex number Z1=A+iB
c     Coded by CJ Neville, 2002 11
c
      COMPLEX (KIND=8) :: CDLOG
      COMPLEX (KIND=8), INTENT(IN) :: Z1
      REAL (KIND=8) :: A   ! The real part of the complex number Z1
      REAL (KIND=8) :: B   ! The imaginary part of the complex number Z1
      REAL (KIND=8) :: R   ! The "ray" of the complex number Z1
      REAL (KIND=8) :: T   ! The "angle" of the complex number Z1
c
      A = REAL (Z1,KIND=8)
      B = AIMAG(Z1)
c
      R = dsqrt((a*a)+(b*b))
      T = datan(b/a)
c
      CDLOG = CMPLX(DLOG(R),T,KIND=8)
c
      END FUNCTION CDLOG
c
C     *******************
      FUNCTION CDSQRT(Z1)
c     *******************
C     This is a function to calculate the square root
c     of a double precision complex number Z1
C     Coded by JP Keizer, 2002 05
c
      COMPLEX (KIND=8) :: CDSQRT
      COMPLEX (KIND=8), INTENT(IN) :: Z1
c
      CDSQRT=Z1**0.5
c
      END FUNCTION CDSQRT
