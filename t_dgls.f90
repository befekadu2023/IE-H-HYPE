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

module t_dgls

USE gear_GlobVARs

USE MODVAR, ONLY : basin,soilthick,realzero, &
	                   numsubstances,   &
                       maxsoillayers,   &
                       genpar,soilpar,landpar, &            
                       classdata,		&
					   dayno,			&														!!Added by BTW JUly11_2020
					   currentdate																!!Added by BTW JUly12_2020
    USE HYPEVARIABLES, ONLY : basinlp, &
                              m_ttmp,  &
                              m_ttrig,m_tredA,m_tredB, &
							  m_perc1,m_perc2,    &
                              m_crate5,   &
                              m_onpercred, m_pppercred, &
							  soilrc

							
IMPLICIT NONE
  PRIVATE 
  
 PUBLIC dgl 
 
 CONTAINS
  
  Subroutine dgl(bspn,nnn,xxx,yyy,f,									&
	 &			ginfilt,barefrac,epot,epotfrac_1,						&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)
	 
	 

	
	
	!REAL*8, INTENT(IN) 	:: ginfilt
	INTEGER, INTENT(IN) :: bspn		!The identifier used for the set of ODEs to be solved...........
	INTEGER, INTENT(IN) :: nnn        !number of simultaneous ODEs....................................
	REAL*8, INTENT(IN) 	:: xxx   		!starting or end point of the time step.........................
	REAL*8, INTENT(IN)	:: yyy(0:nnn-1) !initial value or solution .....................................
    REAL*8, INTENT(OUT)	:: f(0:nnn-1)
	
	REAL, INTENT(IN)    :: ginfilt
	REAL, INTENT(IN)    :: barefrac
	REAL, INTENT(IN)    :: wp_1,wp_2,wp_3
	REAL, INTENT(IN)    :: fc_1,fc_2,fc_3
	REAL, INTENT(IN)    :: ep_1,ep_2,ep_3
	REAL, INTENT(IN)    :: epot
	REAL, INTENT(IN)    :: soilrc_1,soilrc_2,soilrc_3
	REAL, INTENT(IN)    :: epotfrac_1,epotfrac_2
    REAL, INTENT(IN)    :: basinlp_i
    REAL, INTENT(IN)    :: m_perc1,m_perc2
	REAL, INTENT(IN)    :: soiltemp_reduction_1,soiltemp_reduction_2
    
	INTEGER, INTENT(IN) :: isoil
	
	
	
	
	
	 
	 
	!real*8 xxx,yyy(0:nnn-1), f(0:nnn-1)
    !integer bspn, nnn
	
    Select Case(bspn)
      Case (0)
        f(0) = yyy(0) * yyy(1) + DCOS(xxx) - 0.5D0 * DSIN(2.d0 * xxx)
        f(1) = yyy(0)**2 + yyy(1)**2 - (1.d0 + DSIN(xxx))
      Case (1)
		f(0) = -yyy(0) + xxx / ((1.d0+xxx)*(1.d0+xxx))
	
		Case (2)												!!Added by BTW (to test FUSE550 state space equations)
		f(0) = 0- yyy(0)* 0.5*(0*0.5/(0.5*262.5)) - 						&
	 &		Min(yyy(0)*0.07208333*0.5/(0.5*262.5),0.07208333*0.5) - 		&
	 &		(yyy(0))**0.9*0.9*(1./262.5)**0.9						-		&
	 &		yyy(0)*(1.-0.5)*0.5*(1./((1.-0.5)*262.5))	
	 
        f(1) = (yyy(0)**0.9)*0.9*(1./262.5)**0.9 					-		&
	 &		Min(yyy(1)*0.07208333*0.5/(0.5*2750),0.07208333*0.5)	-		&
	 &		yyy(1)**5.5*500.0005*(1./2750)**5.5
      Case (3)
        f(0) = yyy(1)
        f(1) = yyy(2)
        f(2) = yyy(3)
        f(3) = yyy(4)
        f(4) = (45.d0 * yyy(2) * yyy(3) * yyy(4) - 40.d0 * yyy(3)**3) / (9.d0 * yyy(2)**2)
		
	  Case (4)											!!Added by BTW (to represent the HYPE state space formulation by BTW)
        f(0) = ginfilt - MIN(1.,barefrac)*epot*epotfrac_1*				&
	 &	soiltemp_reduction_1*Min(1.,(Max(0.,yyy(0)-wp_1)/				&
	 &	(basinlp_i * fc_1)))			-								&
	 &	Min(m_perc1*Max(0.,Min(1.,(yyy(0)-wp_1-fc_1)/(ep_1))),			&
     &	Max(0.,(wp_2+fc_2+ep_2) - yyy(1)))		-						&
	 &	Max(0.,yyy(0)-fc_1-wp_1)*soilrc_1
	 
		 f(1) = Min(m_perc1*Max(0.,Min(1.,(yyy(0)-wp_1-fc_1)/(ep_1))),					&
     &	Max(0.,(wp_2+fc_2+ep_2) - yyy(1)))		-				&
	 &	MIN(1.,barefrac)*(epot*epotfrac_2*soiltemp_reduction_2*			&
	 &	Min(1.,(Max(0.,yyy(1)-wp_2)/(basinlp_i * fc_2))))		-		&
	 &	Min(m_perc2*Max(0.,Min(1.,(yyy(1)-wp_2-fc_2)/(ep_2))),							&
	 &	Max(0.,(wp_3+fc_3+ep_3) - yyy(2)))					-			&
	 &	Max(0.,yyy(1)-fc_2-wp_2)*soilrc_2
	 
		f(2) =  Min(m_perc2*Max(0.,MIN(1.,(yyy(1)-fc_2-wp_2)/(ep_2))),					&
     &	Max(0.,(wp_3+fc_3+ep_3) - yyy(2)))					-			&
	 &	Max(0.,yyy(2)-fc_3-wp_3)*soilrc_3
	 
	 
	 Case (5)											!!Added by BTW (to represent the HYPE state space formulation by BTW)
		 f(0) = ginfilt - 												&
	 &	Min(m_perc1*Max(0.,Min(1.,(yyy(0)-wp_1-fc_1)/(ep_1))),			&
     &	Max(0.,(wp_2+fc_2+ep_2) - yyy(1)))		-						&
	 &	Max(0.,yyy(0)-fc_1-wp_1)*soilrc_1
	 
		 f(1) = Min(m_perc1*Max(0.,Min(1.,(yyy(0)-wp_1-fc_1)/(ep_1))),	&
     &	Max(0.,(wp_2+fc_2+ep_2) - yyy(1)))					-			&
	 &	Min(m_perc2*Max(0.,Min(1.,(yyy(1)-wp_2-fc_2)/(ep_2))),			&
	 &	Max(0.,(wp_3+fc_3+ep_3) - yyy(2)))			-			&
	 &	Max(0.,yyy(1)-fc_2-wp_2)*soilrc_2
	 
	     f(2) =  Min(m_perc2*Max(0.,Min(1.,(yyy(1)-wp_2-fc_2)/(ep_2))),	&
	 &	Max(0.,(wp_3+fc_3+ep_3) - yyy(2)))					-			&
	 &	Max(0.,yyy(2)-fc_3-wp_3)*soilrc_3
	 
		 
	 
	 

    End Select
  return
  End Subroutine dgl
  
  
  
  

  Subroutine settxt(bspn,nnn,dgltxt)
    integer bspn, nnn
	character*70 dgltxt(0:nnn-1)
    Select Case(bspn)
      Case (0)
        nnn=2
        dgltxt(0) = ' y1'' = y1 * y2 + cos(xxx) - 0.5 * sin(2.0*xxx)'
        dgltxt(1) = ' y2'' = y1 * y1 + y2 * y2 - (1 + sin(xxx))'
      Case (1)
        nnn=1
	dgltxt(0) = ' yyy'' = -yyy + xxx / ((1+xxx)*(1+xxx))'
      Case (3)
        nnn=5
        dgltxt(0) = ' y1'' = y2'
        dgltxt(1) = ' y2'' = y3'
        dgltxt(2) = ' y3'' = y4'
        dgltxt(3) = ' y4'' = y5'
        dgltxt(4) = ' y5'' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)'
    End Select
  return
  End Subroutine settxt
  
  
  !My example of the FSO550 state-space formulation
  
  ! 
! Hydrology_upper_layer <- function(t, state, parameters) {
!   with(as.list(c(state, parameters)), {
!     
!     dS1 <- (
!       Precip
!       -yyy*Phi_tens*(Precip*Ac_max/(Phi_tens*S1_max))
!       - min(yyy*PET*r_1/(Phi_tens*S1_max),PET*r_1)
!       - (yyy)^cx*k_u*(1/S1_max)^cx
!       - yyy*(1-Phi_tens)*ki*(1/((1-Phi_tens)*S1_max))
!       -q_ufof)
!     
!     list(c(dS1))
!   })
! }
! 
! Hydrology_lower_layer <- function(t, state, parameters) {
!   with(as.list(c(state, parameters)), {
!     
!     dS2 <- (
!       q_12 
!       - min(yyy*PET*r_2/(Phi_tens*S2_max),PET*r_2)
!       - yyy^nnn*k_s*(1/S2_max)^n
!       - q_sfof)
!     
!     list(c(dS2))
!   })
! }


!Other examples given in C

! ----------------------- DE system number 0 -----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl0(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(0) * yyy(1) + COS(xxx) - HALf * SIN(TWO * xxx);
!  f(1) = yyy(0) * yyy(0) + yyy(1) * yyy(1) - (ONE + SIN(xxx));
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt0(void)
!{
!  return
!    "y1' = y1 * y2 + cos(xxx) - 0.5 * sin(2.0*xxx)\n"
!    "y2' = y1 * y1 + y2 * y2 - (1 + sin(xxx))\n";
!}


! ----------------------- DE system number 1 -----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl1(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!  f(1) = (ONE + yyy(1) * yyy(1)) / yyy(0);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt1(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = (1 + y2 * y2) / y1\n";
!}


! ----------------------- DE system number 2 -----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl2(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!  f(1) = -fOUR * yyy(0) + EXP(xxx);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt2(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = -4 * y1 + exp(xxx)\n";
!}


!----------------------- DE system number 3 ------------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl3(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!  f(1) = yyy(2);
!  f(2) = yyy(3);
!  f(3) = yyy(4);
!  f(4) = ((REAL)45.0 * yyy(2) * yyy(3) * yyy(4) -
!          (REAL)40.0 * yyy(3) * yyy(3) * yyy(3)) / (NINE * yyy(2) * yyy(2));
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt3(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = y3\n"
!    "y3' = y4\n"
!    "y4' = y5\n"
!    "y5' = (45 * y3 * y4 * y5 - 40 * y4 * y4 * y4) / (9 * y3 * y3)\n";
!}



!----------------------- DE system number 4 ------------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl4(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = -yyy(1);
!  f(1) = yyy(0);
!  f(2) = yyy(4) - yyy(1) * yyy(2);
!  f(3) = yyy(0) * yyy(3) - yyy(5);
!  f(4) = -yyy(2) - yyy(1) * yyy(4);
!  f(5) = yyy(3) + yyy(0) * yyy(5);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt4(void)
!{
!  return
!    "y1' = -y2\n"
!    "y2' = y1\n"
!    "y3' = y5 - y2 * y3\n"
!    "y4' = y1 * y4 - y6\n"
!    "y5' = -y3 - y2 * y5\n"
!    "y6' = y4 + y1 * y6\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem  yyy(0) = (1,0,0,1,e,0)                   *
!*************************************************************************
!static void loesung4(REAL xxx, REAL *yyy)
!{
!  yyy(0) = COS(xxx);
!  yyy(1) = SIN(xxx);
!  yyy(2) = SIN(xxx) * EXP(COS(xxx));
!  yyy(3) = COS(xxx) * EXP(SIN(xxx));
!  yyy(4) = COS(xxx) * EXP(COS(xxx));
!  yyy(5) = SIN(xxx) * EXP(SIN(xxx));
!}


!--------------------------- DE system number 5 --------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl5(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = (REAL)-500.5 * yyy(0) + (REAL)499.5 * yyy(1);
!  f(1) =  (REAL)499.5 * yyy(0) - (REAL)500.5 * yyy(1);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt5(void)
!{
!  return
!    "y1' = -500.5 * y1 + 499.5 * y2\n"
!    "y2' =  499.5 * y1 - 500.5 * y2\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem  yyy(0) = (4,2)                           *
!*************************************************************************
!static void loesung5(REAL xxx, REAL *yyy)
!{
!  yyy(0) = THREE * EXP(-xxx) + EXP((REAL)-1000.0 * xxx);
!  yyy(1) = THREE * EXP(-xxx) - EXP((REAL)-1000.0 * xxx);
!}


!------------------------- DE system number 6 ----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl6(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = (REAL)15.0 * yyy(1);
!  f(1) = (REAL)-0.6 * yyy(0);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt6(void)
!{
!  return
!    "y1' = 15.0 * y2\n"
!    "y2' = -0.6 * y1\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem  yyy(0) = (0,1)                           *
!*************************************************************************
!static void loesung6(REAL xxx, REAL *yyy)
!{
!  yyy(0) = fIVE * SIN(THREE * xxx);
!  yyy(1) = COS(THREE * xxx);
!}



!------------------------- DE system number 7 ----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl7(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!  f(1) = yyy(2);
!  f(2) = yyy(3);
!  f(3) = fOUR * yyy(3) - SIX * yyy(2) + fOUR * yyy(1) - yyy(0);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt7(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = y3\n"
!    "y3' = y4\n"
!    "y4' = 4.0 * y4 - 6.0 * y3 + 4.0 * y2 - y1\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem  yyy(0) = (0,0,0,6)                       *
!*************************************************************************
!static void loesung7(REAL xxx, REAL *yyy)
!{
!  yyy(0) = xxx * xxx * xxx * EXP(xxx);
!  yyy(1) = xxx * xxx * (xxx + THREE) * EXP(xxx);
!  yyy(2) = xxx * (xxx * xxx + SIX * xxx + SIX) * EXP(xxx);
!  yyy(3) = (xxx * xxx * xxx + NINE * xxx * xxx + (REAL)18.0 * xxx + SIX) * EXP(xxx);
!}


!------------------------- DE system number 8 ----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl8(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = (REAL)-52.0 * yyy(0) + (REAL)50.0 * yyy(1);
!  f(1) =  (REAL)50.0 * yyy(0) - (REAL)52.0 * yyy(1);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt8(void)
!{
!  return
!    "y1' = -52 * y1 + 50 * y2\n"
!    "y2' =  50 * y1 - 52 * y2\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem                                         *
!*      yyy(0) = (1265588.55375228,-1265586.55375225)                      *
!*************************************************************************
!static void loesung8(REAL xxx, REAL *yyy)
!{
!  yyy(0) = (REAL)1.000000015       * EXP(        -TWO * xxx) +
!         (REAL)1265587.553752265 * EXP((REAL)-102.0 * xxx);
!  yyy(1) = (REAL)1.000000015 *       EXP(        -TWO * xxx) -
!         (REAL)1265587.553752265 * EXP((REAL)-102.0 * xxx);
!}


!------------------------- DE system number 9 ----------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl9(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!#if defined(BC3)
!  _fpreset();
!#endif
!  f(1) = -yyy(0) * COSH(xxx);
!}


!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt9(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = -y1 * cosh(xxx)\n";
!}


!------------------------- DE system number 10  --------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl10(REAL xxx, REAL *yyy, REAL *f)
!{
!  f(0) = yyy(1);
!  f(1) = -yyy(0) * yyy(0) * yyy(0);
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt10(void)
!{
!  return
!    "y1' = y2\n"
!    "y2' = -y1^3\n";
!}


!------------------------- DE system number 11 ---------------------------

!*************************************************************************
!* Compute the value f of the right hand side of a DE system at (xxx,yyy)    *
!*************************************************************************
!static void dgl11(REAL xxx, REAL *yyy, REAL *f)
!{
!#define dk  ((REAL)-10000.0)
!  f(0) = yyy(1);
!  f(1) = dk * yyy(0) + (dk - ONE) * yyy(1);
!#undef  dk
!}

!*************************************************************************
!* alpha-numeric description of above system                             *
!*************************************************************************
!static char *dgltxt11(void)
!{
!  return
!    "y1' =  y1                                    (dk = -10000)\n"
!    "y2' =  dk * y1 + (dk - 1) * y2\n";
!}

!*************************************************************************
!* Compute the value of the analytic solution yyy(xxx) for the above DE      *
!* for the initial value problem  (y0) = (1,1)                           *
!*************************************************************************
!static void loesung11(REAL xxx, REAL *yyy)
!{
!  yyy(0) = EXP(-xxx);
!  yyy(1) = -EXP(-xxx);
!}

!end of file t_dgls.f90


end module t_dgls