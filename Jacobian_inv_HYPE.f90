MODULE Jacobian_inv_HYPE

  USE Mat_inv_GlobVars
  
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

  CONTAINS
  
  
  function funcv(uxx,uxx_t0,dt,ginfilt,barefrac,epot,epotfrac_1,		&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)

    implicit none

    double precision, intent(in) :: uxx(dim_jacob)
	double precision, intent(in) :: uxx_t0(dim_jacob)
	double precision,intent(in) :: dt
    double precision, dimension(size(uxx)) :: funcv 
	
	
	
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
	
	
  
 
 
! Want Jacobian matric for the following functions:
	 !_____________________________________
	 
	 funcv(1) = uxx(1) - uxx_t0(1) -									&
	 &		dt*(ginfilt - MIN(1.,barefrac)*epot*epotfrac_1*				&
	 &	soiltemp_reduction_1*Min(1.,(Max(0.,uxx(1)-wp_1)/				&
	 &	(basinlp_i * fc_1)))			-								&
	 &	Min(m_perc1*Max(0.,Min(1.,(uxx(1))/(wp_1+fc_1+ep_1))),			&
	 &	Max(0.,uxx(1)-fc_1- wp_1),										&
     &	Max(0.,(wp_2+fc_2+ep_2) - uxx(2)))		-						&
	 &	Max(0.,uxx(1)-fc_1-wp_1)*soilrc_1)
	 
	 funcv(2) = uxx(2) - uxx_t0(2) -									&
	 &		dt*(Min(m_perc1*Max(0.,Min(1.,(uxx(1))/(wp_1+fc_1+ep_1))),	&
	 &	Max(0.,uxx(1)-fc_1- wp_1),										&
     &	Max(0.,(wp_2+fc_2+ep_2) - uxx(2)))		-				&
	 &	MIN(1.,barefrac)*(epot*epotfrac_2*soiltemp_reduction_2*			&
	 &	Min(1.,(Max(0.,uxx(2)-wp_2)/(basinlp_i * fc_2))))		-		&
	 &	Min(m_perc2*Max(0.,Min(1.,(uxx(2))/(wp_2+fc_2+ep_2))),			&
	 &	Max(0.,uxx(2)-fc_2-wp_2),										&
	 &	Max(0.,(wp_3+fc_3+ep_3) - uxx(3)))					-			&
	 &	Max(0.,uxx(2)-fc_2-wp_2)*soilrc_2)
	 
	 funcv(3) =  uxx(3) - uxx_t0(3) -									&
	 &		dt*(Min(m_perc2*Max(0.,MIN(1.,(uxx(2))/(fc_2+wp_2+ep_2))),	&
	 &			Max(0.,uxx(2)-fc_2-wp_2),								&
     &	Max(0.,(wp_3+fc_3+ep_3) - uxx(3)))					-			&
	 &	Max(0.,uxx(3)-fc_3-wp_3)*soilrc_3)
	 
	 !_____________________________________
	 
 
 end function funcv




subroutine Jacob_forward_Diff(uxx,uxx_t0,JacobMatrx,F_numerator,dt,		&
	 &	ginfilt,barefrac,epot,epotfrac_1,								&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)


   implicit none

   double precision, intent(inout) :: uxx(dim_jacob)
   double precision, intent(out) :: JacobMatrx(dim_jacob,dim_jacob)
   double precision, intent(out) :: F_numerator(dim_jacob)
   double precision, intent(in) :: uxx_t0(dim_jacob)
   double precision,intent(in) :: dt
  
   double precision, parameter :: EPS = 1.0e-04
   integer :: jjx,dim_jacobx
   double precision, dimension(size(uxx)) :: uxx_save, uxph1, uxph2, uh, uh_diff, ux2, ux1
   
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

     dim_jacobx = size(uxx)

     uxx_save = uxx
     uh = EPS*abs(uxx_save)

     where (uh <= EPS) uh = EPS

     uxph2 = uxx_save + uh              ! Trick to reduce finite precision error.
     uxph1 = uxx_save
     uh_diff = uxph2 - uxph1
	 
	 F_numerator = funcv(uxx,uxx_t0,dt,ginfilt,barefrac,epot,epotfrac_1,&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)
	 
     !
     ux2 = uxx
     ux1 = uxx
  do jjx = 1 , dim_jacobx
       ux2(jjx) = uxph2(jjx)
       ux1(jjx) = uxph1(jjx)
	   
	   
       ! Forward difference formula
	   
JacobMatrx(:,jjx)=((funcv(ux2,uxx_t0,dt,ginfilt,barefrac,epot,epotfrac_1,	&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)) - 		&
	 &(funcv(ux1,uxx_t0,dt,ginfilt,barefrac,epot,epotfrac_1,	&
	 &			soiltemp_reduction_1,wp_1,								&
	 &			basinlp_i,fc_1,m_perc1,isoil,ep_1,						&
	 &			wp_2,fc_2,ep_2,soilrc_1, soilrc_2, soilrc_3,       		&
	 &			epotfrac_2,												&
	 &			soiltemp_reduction_2,m_perc2,wp_3,fc_3,ep_3)))/uh_diff(jjx)
	 
       ux2(jjx) = uxx_save(jjx)
       ux1(jjx) = uxx_save(jjx)
	   
  end do

end subroutine Jacob_forward_Diff









  
  !subroutine Allocate_mat_inv_variables(uxx,JacobMatrx,JacobMatrxCopy,	&
  !   &inv_Jacob,row_permutat,matrx_multip,S_new,S_old,S_t0,F_numerator,	&
	! &	dim_jacob)
	 
 subroutine Allocate_mat_inv_variables

use Mat_inv_GlobVars
implicit none

 
  
INTEGER :: iix,jjx

  allocate(uxx(dim_jacob))
  allocate(JacobMatrx(dim_jacob,dim_jacob))
  allocate(JacobMatrxCopy(dim_jacob,dim_jacob))
  allocate(inv_Jacob(dim_jacob,dim_jacob))
  allocate(row_permutat(dim_jacob))
  allocate(matrx_multip(dim_jacob))
  allocate(S_new(dim_jacob))
  allocate(S_old(dim_jacob))
  allocate(S_t0(dim_jacob))
  allocate(F_numerator(dim_jacob))
 
  
  do iix=1, dim_jacob 

    do jjx=1, dim_jacob

	  inv_Jacob(iix,jjx) = 0.d0
    end do
    inv_Jacob(iix,iix) = 1.d0
  end do


end subroutine Allocate_mat_inv_variables


!subroutine deAllocate_mat_inv_variables(uxx,JacobMatrx,JacobMatrxCopy,	&
!     &inv_Jacob,row_permutat,matrx_multip,S_new,S_old,S_t0,F_numerator)

subroutine deAllocate_mat_inv_variables
	 
	use Mat_inv_GlobVars

	IF(ALLOCATED(uxx))    DEALLOCATE(uxx)
	IF(ALLOCATED(JacobMatrx))    DEALLOCATE(JacobMatrx)
	IF(ALLOCATED(JacobMatrxCopy))    DEALLOCATE(JacobMatrxCopy)
	IF(ALLOCATED(inv_Jacob))    DEALLOCATE(inv_Jacob)
	IF(ALLOCATED(row_permutat))    DEALLOCATE(row_permutat)
	IF(ALLOCATED(matrx_multip))    DEALLOCATE(matrx_multip)
	IF(ALLOCATED(S_new))    DEALLOCATE(S_new)
	IF(ALLOCATED(S_old))    DEALLOCATE(S_old)
	IF(ALLOCATED(S_t0))    DEALLOCATE(S_t0)
	IF(ALLOCATED(F_numerator))    DEALLOCATE(F_numerator)


end subroutine deAllocate_mat_inv_variables
  

  
  

!The following subroutines are adapted Base on LU modules http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/lu_f90.txt


!*******************************************************
!*    LU decomposition routines    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* Based on "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************

!  ***************************************************************
!  * Given an dim_jacob x dim_jacob matrix JacobMatrx, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. JacobMatrx is input and output, and dim_jacob   *
!  * are input. row_permutat is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return CODE_mtrx_mtrx is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(JacobMatrx,dim_jacob,row_permutat,CODE_mtrx)
 !use mat_inv_variables
 implicit none
  INTEGER,INTENT(IN):: dim_jacob
  double precision, INTENT(INOUT):: JacobMatrx(dim_jacob,dim_jacob)
  INTEGER,INTENT(OUT):: row_permutat(dim_jacob)
  INTEGER ,INTENT(OUT):: CODE_mtrx
  INTEGER :: NMAX,iix,jjx,kkx
  REAL*8 :: TINYx
 PARAMETER(NMAX=100,TINYx=1.5D-16)
 !REAL*8  AMAX,DUM, SUM, JacobMatrx(dim_jacob,dim_jacob),VV(NMAX)
 REAL*8  AMAX,DUM, SUM, VV(NMAX)
 INTEGER :: D,IMAX
 !INTEGER :: row_permutat(dim_jacob)

 D=1; CODE_mtrx=0

 DO iix=1,dim_jacob
   AMAX=0.d0
   DO jjx=1,dim_jacob
     IF (DABS(JacobMatrx(iix,jjx)).GT.AMAX) AMAX=DABS(JacobMatrx(iix,jjx))
   END DO ! jjx loop
   IF(AMAX.LT.TINYx) THEN
     CODE_mtrx = 1
     RETURN
   END IF
   VV(iix) = 1.d0 / AMAX
 END DO ! iix loop

 DO jjx=1,dim_jacob
   DO iix=1,jjx-1
     SUM = JacobMatrx(iix,jjx)
     DO kkx=1,iix-1
       SUM = SUM - JacobMatrx(iix,kkx)*JacobMatrx(kkx,jjx) 
     END DO ! kkx loop
     JacobMatrx(iix,jjx) = SUM
   END DO ! iix loop
   AMAX = 0.d0
   DO iix=jjx,dim_jacob
     SUM = JacobMatrx(iix,jjx)
     DO kkx=1,jjx-1
       SUM = SUM - JacobMatrx(iix,kkx)*JacobMatrx(kkx,jjx) 
     END DO ! kkx loop
     JacobMatrx(iix,jjx) = SUM
     DUM = VV(iix)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = iix
       AMAX = DUM
     END IF
   END DO ! iix loop  
   
   IF(jjx.NE.IMAX) THEN
     DO kkx=1,dim_jacob
       DUM = JacobMatrx(IMAX,kkx)
       JacobMatrx(IMAX,kkx) = JacobMatrx(jjx,kkx)
       JacobMatrx(jjx,kkx) = DUM
     END DO ! kkx loop
     D = -D
     VV(IMAX) = VV(jjx)
   END IF

   row_permutat(jjx) = IMAX
   IF(DABS(JacobMatrx(jjx,jjx)) < TINYx) JacobMatrx(jjx,jjx) = TINYx

   IF(jjx.NE.dim_jacob) THEN
     DUM = 1.d0 / JacobMatrx(jjx,jjx)
     DO iix=jjx+1,dim_jacob
       JacobMatrx(iix,jjx) = JacobMatrx(iix,jjx)*DUM
     END DO ! iix loop
   END IF 
 END DO ! jjx loop

 RETURN
 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of dim_jacob linear equations JacobMatrx . X = inv_Jacob.  Here JacobMatrx is     *
!  * input, not as the matrix JacobMatrx but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. row_permutat is input as the permuta-*
!  * tion vector returned by LUDCMP. inv_Jacob is input as the right-hand   *
!  * side vector inv_Jacob, and returns with the solution vector X. A, dim_jacob and*
!  * row_permutat are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(JacobMatrx,dim_jacob,row_permutat,out_vector)
 !use mat_inv_variables
 implicit none
 INTEGER, INTENT(IN):: dim_jacob
 double precision, INTENT(IN) :: JacobMatrx(dim_jacob,dim_jacob)
 INTEGER,INTENT(IN):: row_permutat(dim_jacob)
 double precision,INTENT(INOUT)::  out_vector(dim_jacob)
 
 INTEGER :: iix,jjx,II,LL
 REAL*8  SUM
 !REAL*8  SUM, out_vector(dim_jacob)
 !INTEGER row_permutat(dim_jacob)

 II = 0

 DO iix=1,dim_jacob
   LL = row_permutat(iix)
   SUM = out_vector(LL)
   out_vector(LL) = out_vector(iix)
   IF(II.NE.0) THEN
     DO jjx=II,iix-1
       SUM = SUM - JacobMatrx(iix,jjx)*out_vector(jjx)
     END DO ! jjx loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = iix
   END IF
   out_vector(iix) = SUM
 END DO ! iix loop
 

 DO iix=dim_jacob,1,-1
   SUM = out_vector(iix)
   IF(iix < dim_jacob) THEN
     DO jjx=iix+1,dim_jacob
       SUM = SUM - JacobMatrx(iix,jjx)*out_vector(jjx)
     END DO ! jjx loop
   END IF
   out_vector(iix) = SUM / JacobMatrx(iix,iix)
   
 END DO ! iix loop
 
 RETURN
 END subroutine LUBKSB
 
 
 
 
 SUBROUTINE MATMULT(matx1,matx2,matx_out,dim_1,dim_2)                                              
!*******************************************                                     
!*       MULTIPLY TWO REAL MATRICES        *
!* --------------------------------------- *                                     
!* INPUTS:    matx1  MATRIX dim_1*dim_1                *                                     
!*            matx2  MATRIX dim_1*dim_2                *                                     
!*            dim_1  INTEGER                   *                                     
!*            dim_2  INTEGER                   *                                     
!* --------------------------------------- *                                     
!* OUTPUT:    matx_out  MATRIX dim_1*dim_2, PRODUCT matx1*matx2   *                                     
!*                                         *                                     
!******************************************* 
  IMPLICIT NONE
  double precision, INTENT(IN):: matx1(dim_jacob,dim_jacob), matx2(dim_jacob,dim_jacob)
  INTEGER, INTENT(IN):: dim_1, dim_2
  double precision, INTENT(OUT):: matx_out(dim_jacob,dim_jacob)
  !REAL*8 matx1(dim_1,dim_1),matx2(dim_1,dim_2),matx_out(dim_1,dim_2),SUM  
  REAL*8 :: SUM  
  INTEGER :: iix,jjx,kkx
  DO iix=1,dim_1                                                                  
    DO jjx=1,dim_2                                                                
      SUM=0.                                                                
      DO kkx=1,dim_1                                                              
        SUM=SUM+matx1(iix,kkx)*matx2(kkx,jjx)                                               
      ENDDO                                                                 
      matx_out(iix,jjx)=SUM                                                            
    ENDDO                                                                   
  ENDDO                                                                     
  RETURN                                                                    
  END 
  
  
  
  !*****************************************
  ! Read a vector from unit u from index = 1
  !*****************************************
  Subroutine ReadVec1 (u, n, x) 
  INTEGER u
  REAL*8 x(1:n)
  read(u,*) (x(iix),iix=1,n)
  return
  end subroutine ReadVec1
  
  
  
  
  

 
 

END MODULE Jacobian_inv_HYPE

