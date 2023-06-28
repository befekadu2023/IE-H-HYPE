!************************************************************************
!*                                                                      *
!* Test examples for methods to solve first order ordinary systems of   *
!* differential equations.                                              *
!*                                                                      *
!* We supply the right hand sides, their alpha-numeric description and  *
!* the exact analytic solution, if known.                               *
!*                                                                      *
!* When running the main test program, the user can select among the    *
!* examples below.                                                      *
!*                                                                      *
!*                               F90 1.1 Release By J-P Moreau, Paris.  *
!* -------------------------------------------------------------------- *
!* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
!*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
!*                                                                      *
!* Release 1.1: added example #1 for test program tawp.                 *
!************************************************************************

! Note: here, only examples #0, #1 and #3 are implemented in Fortran

module gear_GlobVARs

							
IMPLICIT NONE
  SAVE
  
    REAL*8 	:: xxx   		!starting or end point of the time step......................
	REAL*8 	:: xend   	!desired end point (> xxx).....................................
	INTEGER :: bspn		!The identifier used for the set of ODEs to be solved........
	INTEGER :: nnn        !number of simultaneous ODEs.................................
	REAL*8	:: yyy(0:5) !initial value or solution ..................................
	REAL*8	::epsabs    !absolute error bound (to be provided by the user)...........
    REAL*8	::epsrel	!relative error bound .......................................
	REAL*8	:: h		!!starting or final step size................................
	INTEGER ::fmax       !maximal # of calls of  dgl() ..............................
    INTEGER :: aufrufe	!actual # of calls of  dgl() ................................
    INTEGER ::fehler                  !error code ...................................
	

end module gear_GlobVARs