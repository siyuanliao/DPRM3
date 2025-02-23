MODULE pchipmod

  ! subroutines and functions related to the calculation of a
  ! Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)

  IMPLICIT NONE
  SAVE

  INTEGER, PRIVATE       :: N
  REAL (KIND=8), PRIVATE :: h
  REAL (KIND=8), PRIVATE :: fprime_r, fprime_i

CONTAINS

  SUBROUTINE PCHIPinProgress( x, y, N1 )

    ! This is a rough draft; it hasn't been tested

    ! Monotone Piecewise Cubic Hermite Interpolating Polynomial of SSP data
    !  This is the variant of PCHIP described by Fritch-Butner in UCRL-87559
    !
    !  Also see: F. N. Fritsch and R. E. Carlson. "Monotone Piecewise Cubic
    !  Interpolation. Siam Journal on Numerical Analysis, 17:238-246, Apr 1980.

    ! x is a vector of nodal coordinates
    ! y is a vector with the value of the function at those coordinates
    ! N1 is the number of nodes

    INTEGER,           INTENT( IN ) :: N1
    INTEGER           :: ix, iSeg
    REAL     (KIND=8) :: h1, h2, xt, x( 5 )
    COMPLEX  (KIND=8) :: c0, c1, c2, c3, del1, del2, f1, f2, f1prime, f2prime, y( 5 ), Cubic

    N    = N1 - 1   ! number of segements is number of nodes - 1
    h    = ( x( N1 ) - x( 1 ) ) / N
    iSeg  = 1

    DO ix = 1, N1
       xt = x( 1 ) + ( ix - 1 ) * h
       IF ( ix == N1 ) xt = x( N1 )   ! Make sure no overshoot

       ! search through segments for active one
       DO WHILE ( xt > x( iSeg + 1 ) )
          iSeg = iSeg + 1
       END DO

       ! *** estimate the derivative at the left endpoint of this segment ***

       IF ( ix == 1 ) THEN        ! left endpoint (non-centered 3-=point formula, see UCRL-85104)
          CALL h_del( x, y, ix + 1, h1, h2, del1, del2 )
          f1prime = ( ( 2.0D0 * h1 + h2 ) * del1 - h1 * del2 ) / ( h1 + h2 );
          f1prime = fprime_left_end_Cmplx( del1, del2, f1prime )
       ELSE                       ! interior node (see UCRL-87559)
          CALL h_del( x, y, ix,     h1, h2, del1, del2 )
          f1prime = fprime_Cmplx( h1, h2, del1, del2 )
       END IF


       ! *** estimate the derivative at the right endpoint of this segment ***

       IF ( ix == N - 1 ) THEN   ! right endpoint (non-centered 3-point formula, see UCRL-85104)
          CALL h_del( x, y, ix,     h1, h2, del1, del2 )
          f2prime = ( -h2 * del1 + ( h1 + 2.0D0 * h2 ) * del2 ) / ( h1 + h2 );
          f2prime = fprime_right_end_Cmplx( del1, del2, f2prime )
       ELSE                       ! interior node (see UCRL-87559)
          CALL h_del( x, y, ix + 1, h1, h2, del1, del2 )
          f2prime = fprime_Cmplx( h1, h2, del1, del2 )
       END IF
     
       !!!!!y( ix ) = Cubic( xt - x( ix ), f1, f1prime, f2, f2prime, h )   ! here's the Hermite cubic

    END DO

    RETURN
  END SUBROUTINE PCHIPinProgress

  !**********************************************************************!

  SUBROUTINE h_del( x, y, ix, h1, h2, del1, del2 )
    INTEGER,          INTENT( IN  ) :: ix   ! index of the center point
    REAL    (KIND=8), INTENT( IN  ) :: x( * )
    COMPLEX (KIND=8), INTENT( IN  ) :: y( * )
    REAL    (KIND=8), INTENT( OUT ) :: h1, h2
    COMPLEX (KIND=8), INTENT( OUT ) :: del1, del2

    h1   =   x( ix     ) - x( ix - 1 )
    h2   =   x( ix + 1 ) - x( ix     )

    del1 = ( y( ix     ) - y( ix - 1 ) ) / h1
    del2 = ( y( ix + 1 ) - y( ix     ) ) / h2
  END SUBROUTINE h_del

  !**********************************************************************!

  FUNCTION fprime_Cmplx( h1, h2, del1, del2 )

    REAL    (KIND=8), INTENT( IN ) :: h1, h2
    COMPLEX (KIND=8), INTENT( IN ) :: del1, del2
    COMPLEX (KIND=8)               :: fprime_cmplx

    fprime_r = fprime( h1, h2, REAL(  del1 ), REAL(  del2 ) )
    fprime_i = fprime( h1, h2, AIMAG( del1 ), AIMAG( del2 ) )

    fprime_Cmplx = CMPLX( fprime_r, fprime_i, KIND=8 )

  END FUNCTION fprime_Cmplx

  !**********************************************************************!

  FUNCTION fprime_left_end_Cmplx( del1, del2, fprime )

    COMPLEX (KIND=8), INTENT( IN ) :: del1, del2, fprime
    COMPLEX (KIND=8)               :: fprime_left_end_Cmplx

    fprime_r = fprime_left_end( REAL(  del1 ), REAL(  del2 ), REAL(  fprime ) )
    fprime_i = fprime_left_end( AIMAG( del1 ), AIMAG( del2 ), AIMAG( fprime ) )

    fprime_left_end_Cmplx = CMPLX( fprime_r, fprime_i, KIND=8 )

  END FUNCTION fprime_left_end_Cmplx

  !**********************************************************************!

  FUNCTION fprime_right_end_Cmplx( del1, del2, fprime )

    COMPLEX (KIND=8), INTENT( IN ) :: del1, del2, fprime
    COMPLEX (KIND=8)               :: fprime_right_end_Cmplx

    fprime_r = fprime_right_end( REAL(  del1 ), REAL(  del2 ), REAL(  fprime ) )
    fprime_i = fprime_right_end( AIMAG( del1 ), AIMAG( del2 ), AIMAG( fprime ) )

    fprime_right_end_Cmplx = CMPLX( fprime_r, fprime_i, KIND=8 )

  END FUNCTION fprime_right_end_Cmplx

 !**********************************************************************!

  FUNCTION fprime( h1, h2, del1, del2 )

    REAL (KIND=8), INTENT( IN ) :: h1, h2, del1, del2
    REAL (KIND=8)               :: h, fprime

    h = ( h1 + 2.0D0 * h2 ) / ( 3.0D0 * ( h1 + h2 ) );

    IF ( del1 * del2 > 0.0 ) THEN
       ! slope of the secant lines have the same arithmetic sign
       fprime = ( del1 * del2 ) / ( ( 1.0D0 - h ) * del1 + h * del2 );
    ELSE
       ! set to zero if the slope of the secant lines change sign
       fprime = 0.0;
    END IF

  END FUNCTION fprime

  !**********************************************************************!

  FUNCTION fprime_left_end( del1, del2, fprime )

    REAL (KIND=8), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=8)               :: fprime_left_end

    fprime_left_end = fprime

    IF ( del1 * fprime <= 0.0D0 ) THEN
       ! set derivative to zero if the sign differs from sign of secant slope
       fprime_left_end = 0.0;
    ELSE IF ( ( del1 * del2 <= 0.0D0 ) .AND. ( ABS( fprime ) > ABS( 3.0D0 * del1 ) ) ) THEN
       ! adjust derivative value to enforce monotonicity
       fprime_left_end = 3.0D0 * del1;
    END IF

  END FUNCTION fprime_left_end

  !**********************************************************************!

  FUNCTION fprime_right_end( del1, del2, fprime )

    ! This is essentially the same as fprime_left_end( del2, del1, fprime )
    ! Written separately for clarity

    REAL (KIND=8), INTENT( IN ) :: del1, del2, fprime
    REAL (KIND=8)               :: fprime_right_end

    fprime_right_end = fprime

    IF ( del2 * fprime <= 0.0D0 ) THEN
       ! set derivative to zero if the sign differs from sign of secant slope
       fprime_right_end = 0.0;
    ELSE IF ( ( del1 * del2 <= 0.0D0 ) .AND. ( ABS( fprime ) > ABS( 3.0D0 * del2 ) ) ) THEN
       ! adjust derivative value to enforce monotonicity
       fprime_right_end = 3.0D0 * del2;
    END IF
    
  END FUNCTION fprime_right_end

  !**********************************************************************!

  FUNCTION Cubic( x, f1, f1prime, f2, f2prime, h )
    REAL    (KIND=8), INTENT( IN  ) :: x, h
    COMPLEX (KIND=8), INTENT( IN  ) :: f1, f1prime, f2, f2prime
    COMPLEX (KIND=8) :: Cubic, c0, c1, c2, c3

    !                                                           2      3
    ! compute coefficients of cubic polynomial: c0 + c1*x + c2*x + c3*x
    !
    ! x is a local coordinate in [ 0, h ]

    c0 = f1
    c1 = f1prime
    c2 = ( 3.0D0 * ( f2 - f1 ) - h * ( 2.0D0 * f1prime + f2prime ) ) / h**2
    c3 = ( h * ( f1prime + f2prime ) - 2.0D0 * ( f2 - f1 ) ) / h**3

    Cubic = c0 + ( c1 + ( c2 + c3 * x ) * x ) * x   ! here's the Hermite cubic

  END FUNCTION Cubic
END MODULE pchipmod
