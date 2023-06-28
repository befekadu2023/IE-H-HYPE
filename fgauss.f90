!*****************************************************************
!*         Gauss algorithm for solving linear equations          *
!* ------------------------------------------------------------- *
!* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges *
!*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].    *
!*                                                               *
!*                            F90 Release By J-P Moreau, Paris.  *
!*                                   (www.jpmoreau.fr)           *
!*****************************************************************

module fgauss

IMPLICIT NONE
  PRIVATE 
  
 PUBLIC gauss

CONTAINS

Subroutine ISWAP(a,b)
  integer a,b,tmp
  tmp=b; b=a; a=tmp
  return
End

Subroutine SWAP(a,b)
  real*8 a,b,tmp
  tmp=b; b=a; a=tmp
  return
End

! Gauss algorithm for solving linear equations .......................
Subroutine gauss(mode,n,mat,lumat,perm,b,x,signd,code)
  integer  mode                    ! Modus: 0, 1, 2, 3 ...............
  integer  n                       ! Dimension of matrix .............
  real*8   mat(0:n-1,0:n-1), &     ! Input matrix ....................
           lumat(0:n-1,0:n-1)      ! LU decomposition ................
  integer  perm(0:n-1,0:n-1)       ! row remutation vector ...........
  real*8   b(0:n-1),         &     ! right hand side .................
           x(0:n-1)                ! solution of the system ..........
  integer  signd,            &     ! sign of the permutation .........
           code                    ! return error code ...............
 !*====================================================================*
 !*                                                                    *
 !*  The procedure gauss solves a linear system :  mat * x = b.        *
 !*  Here mat is the nonsingular system matrix, b the right hand side  *
 !*  of the system and x the solution vector.                          *
 !*                                                                    *
 !*  gauss uses the Gauss algorithm and computes a triangular factori- *
 !*  zation of mat and scaled column pivot search.  (Crout method with *
 !*  row swaps).                                                       *
 !*                                                                    *
 !* ------------------------------------------------------------------ *
 !*                                                                    *
 !*   Application:                                                     *
 !*   ============                                                     *
 !*      Solve general linear system with a nonsingular coefficient    *
 !*      matrix.                                                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Control parameter:                                               *
 !*   ==================                                               *
 !*      mode     integer;                                             *
 !*               calling modus for gauss:                             *
 !*       = 0     Find factorization and solve linear system           *
 !*       = 1     Find factorization only.                             *
 !*       = 2     Solve linear system only; the factorization is       *
 !*               already available in lumat. This saves work when     *
 !*               solving a linear system repeatedly for several right *
 !*               hand sides and the same system matrix such as when   *
 !*               inverting the matrix.                                *
 !*       = 3     as under 2, additionally we improve the solution     *
 !*               via iterative refinement (not available here).       *
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      n        integer;  (n > 0)                                    *
 !*               Dimension of mat and lumat,                          *
 !*               size of the vector b, the right hand side, the       *
 !*               solution x and the permutation vector perm.          *
 !*      mat      matrix of the linear system.                         *
 !*      lumat    (for mode = 2, 3)                                    *
 !*               LU factors of mat                                    *
 !*               lumat can be stored in the space of mat.             *
 !*      perm     (for mode = 2, 3)                                    *
 !*               Permutation vector, of the row interchangfes in lumat*
 !*      b        (for mode = 0, 2, 3)                                 *
 !*               Right hand side of the system.                       *
 !*      signd    (for mode = 2, 3)                                    *
 !*               sign of the permutation in perm; the determinant of  *
 !*               mat can be computed as the product of the diagonal   *
 !*               entries of lumat times signd.                        *
 !*                                                                    *
 !*   Output parameters:                                               *
 !*   ==================                                               *
 !*      lumat    (for mode = 0, 1)                                    *
 !*               LU factorization of mat.                             *
 !*      perm     (for mode = 0, 1)                                    *
 !*               row ermutation vector                                *
 !*      x        (for mode = 0, 2, 3)                                 *
 !*               solution vector.                                     *
 !*      signd    (for mode = 0, 1)                                    *
 !*               sign of perm.                                        *
 !*                                                                    *
 !*   Return value (code):                                             *
 !*   ===================                                              *
 !*      =-1      Max. number (MAXITER) of iterative refinements       *
 !*               reached (MAXITER) while mod = 3                      *
 !*      = 0      all ok                                               *
 !*      = 1      n < 1 or other invalid input                         *
 !*      = 2      lack of memory                                       *
 !*      = 3      Matrix singular                                      *
 !*      = 4      Matrix numerically singular                          *
 !*      = 5      incorrect call                                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   subroutines used:                                                *
 !*   ================                                                 *
 !*                                                                    *
 !*      gaudec: determines  LU decomposition                          *
 !*      gausol: solves the linear system                              *
 !*                                                                    *
 !*====================================================================*
 integer rc

  if (n < 1) then
    code = 1
    return
  end if

  Select Case(mode)

    Case(0)    !Find factorization and solve system ...................
       call gaudec (n, mat, lumat, perm, signd,rc)
       if (rc.eq.0) then
         call gausol(n, lumat, perm, b, x,rc)
         code=rc
       else
         code=rc
         return
       end if

    Case(1)    !Find factorization only ...............................
       call gaudec(n, mat, lumat, perm, signd,rc)
       code=rc
       return  

    Case(2)    !Solve only ............................................
       call gausol(n, lumat, perm, b, x,rc)
       code=rc
       return

    Case(3)    !Solve and then use iterative refinement ...............
       print *,' fgauss: gausoli not implemented.'
	   code=5
       return

  End Select

  code = 5                                                 !Wrong call
  return
End


! Gauss decomposition ................................................
Subroutine gaudec(n,mat,lumat,perm,signd,rc)
  integer   n                      ! size of matrix ..................
  real*8    mat(0:n-1,0:n-1),  &   ! Input matrix ....................
            lumat(0:n-1,0:n-1)     ! matrix decomposition ............
  integer   perm(0:n-1),       &   ! row interchanges ................
            signd,             &   ! sign of perm ....................
            rc                     ! error code ......................
 !*====================================================================*
 !*                                                                    *
 !*  gaudec decomposes a nonsingular n x n matrix into a product of a  *
 !*  unit lower and an upper triangular matrix. Both triangular factors*
 !*  are stored in lumat (minus the unit diagonal, of course).         *
 !*                                                                    *
 !* ------------------------------------------------------------------ *
 !*                                                                    *
 !*   Input parameter:                                                 *
 !*   ================                                                 *
 !*      n        integer;  (n > 0)                                    *
 !*               Dimension of  mat and lumat,                         *
 !*               size of  b , x and perm.                             *
 !*      mat      original system matrix                               *
 !*                                                                    *
 !*   Output parameters:                                               *
 !*   ==================                                               *
 !*      lumat    LU factorization                                     *
 !*      perm     row permutation vector for lumat                     *
 !*      signd    sign of perm. The determinant of mat can be computed *
 !*               as the product of the diagonal entries of lumat times*
 !*               signd.                                               *
 !*                                                                    *
 !*   Return value (rc):                                               *
 !*   =================                                                *
 !*      = 0      all ok                                               *
 !*      = 1      n < 1 or invalid input                               *
 !*      = 2      lack of memory                                       *
 !*      = 3      Matrix is singular                                   *
 !*      = 4      Matrix numerically singular                          *
 !*                                                                    *
 !*====================================================================*
  real*8, parameter :: MACH_EPS = 1.2D-16
  integer m, j, i, j0
  real*8 piv, tmp, zmax
  real*8 d(0:n-1)              !scaling vector for pivoting

  if (n < 1) then
    rc=1                                !Invalid parameters
    return
  end if

  !copy mat to lumat
  do i=0, n-1
    do j=0, n-1
      lumat(i,j) = mat(i,j)
    end do
  end do

  do i = 0, n-1
    perm(i) = i                          ! Initialize perm
    zmax = 0.d0
    do j = 0, n-1                        !find row maxima
      tmp = DABS(lumat(i,j))
      if (tmp > zmax)  zmax = tmp
    end do

    if (zmax == 0.d0) then               !mat is singular
      rc = 3
      return
    end if
    d(i) = 1.d0 / zmax
  end do

  signd = 1                              !initialize sign of perm

  do i = 0, n-1
    piv = DABS(lumat(i,i)) * d(i)
    j0 = i                               !Search for pivot element
    do j = i+1, n-1
      tmp = DABS(lumat(j,i)) * d(j)
      if (piv < tmp) then
        piv = tmp                        !Mark pivot element and
        j0 = j                           !its location
      end if
    end do

    if (piv < MACH_EPS) then             !If piv is small, mat is
      signd = 0                          !nearly singular
      rc =4
      return
    end if

    if (j0.ne.i) then
      signd = - signd                    !update signd
      call ISWAP(perm(j0), perm(i))      !swap pivot entries
      call SWAP(d(j0), d(i))             !swap scaling vector
      do j=0, n-1                 
	    !swap j0-th and i-th rows of lumat
        call SWAP(lumat(j0,j), lumat(i,j)) 
      end do
    end if

    do j = i+1, n-1                      !Gauss elimination
      if (lumat(j,i).ne.0.D0) then
        lumat(j,i) = lumat(j,i) / lumat(i,i)
        tmp = lumat(j,i)
        do m = i+1, n-1
          lumat(j,m) = lumat(j,m) - tmp * lumat(i,m)
        end do
      end if
    end do
  end do !i loop

  rc=0   !all ok
  return
End


! Gauss solution .....................................................
Subroutine gausol(n,lumat,perm,b,x,rc)
    integer n                      ! size of matrix ..................
    real*8  lumat(0:n-1,0:n-1)     ! decomposed matrix (LU) ..........
    integer perm(0:n-1)            ! row permutation vector ..........
    real*8  b(0:n-1),        &     ! Right hand side .................
            x(0:n-1)               ! solution ........................
    integer rc                     ! error code ......................
 !====================================================================*
 !                                                                    *
 !  gausol  finds the solution x of the linear system  lumat * x = b  *
 !  for the product matrix lumat, that describes an LU decomposition, *
 !  as produced by gaudec.                                            *
 !                                                                    *
 !====================================================================*
 !                                                                    *
 !   Input parameters:                                                *
 !   ================                                                 *
 !      n        integer;  (n > 0)                                    *
 !               Dimension of lumat,                                  *
 !      lumat    LU factorization, as produced from gaudec            *
 !      perm     row permutation vector for lumat                     *
 !      b        right hand side of the system.                       *
 !                                                                    *
 !   Output parameter:                                                *
 !   ================                                                 *
 !      x        solution vector                                      *
 !                                                                    *
 !   Return value (rc):                                               *
 !   =================                                                *
 !      = 0      all ok                                               *
 !      = 1      n < 1 or other invalid input parameter               *
 !      = 3      improper LU decomposition ( zero diagonal entry)     *
 !                                                                    *
 !====================================================================*
  integer j,k
  real*8 sum, MACH_EPS

  if (n < 1) then
     rc=1                                     !Invalid input parameter
     return
  end if

  MACH_EPS = 1.2D-16

  do k=0, n-1                                                !update b
    x(k) = b(perm(k))
    do j=0, k-1
      x(k) = x(k) - lumat(k,j) * x(j)
    end do
  end do

  do k=n-1, 0, -1                                     !back substitute
    sum = 0.D0
    do j=k+1, n-1
      sum = sum + lumat(k,j) * x(j)
    end do

    if (DABS(lumat(k,k)) < MACH_EPS) then
      rc=3
      return
    end if
    x(k) = (x(k) - sum) / lumat(k,k)
  end do

  rc = 0  !all ok
  return

End


! Determinant  .......................................................
real*8 Function det(n, mat) 
       integer n                    !Dimension of the matrix .........
       real*8 mat(0:n-1,0:n-1)      !matrix ..........................
 !*====================================================================*
 !*                                                                    *
 !*  det computes the determinant of an n x n real matrix mat          *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameter:                                                 *
 !*   ================                                                 *
 !*      n        integer;  (n > 0)                                    *
 !*               Dimension of mat                                     *
 !*      mat      n x n matrix                                         *
 !*                                                                    *
 !*   Return value:                                                    *
 !*   =============                                                    *
 !*      REAL     Determinant of mat.                                  *
 !*               If the return value = 0, then the matrix is singular *
 !*               or the storage is insufficient                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   subroutine used:                                                 *
 !*   ================                                                 *
 !*                                                                    *
 !*          gaudec ():    LU decomposition of mat                     *
 !*                                                                    *
 !*====================================================================*
  integer i, rc, signd, perm(0:n-1)
  real*8 lu(0:n-1,0:n-1), tmpdet
  real*8 MACH_EPS, MAXROOT

  if (n < 1) then
    det = 0.D0
    return
  end if
 
  MACH_EPS = 1.2D-16
  MAXROOT  = 1.D16

  call gaudec(n, mat, lu, perm, signd, rc)                 !decompose

  if (rc.ne.0.or.signd.eq.0) then
    det=0.D0
    return
  end if

  tmpdet = 1.d0*signd

  do i=0, n-1
    if (DABS(tmpdet) < MACH_EPS) then
      det = 0.D0
      return
    else if (DABS(tmpdet) > MAXROOT.or.DABS(lu(i,i)) > MAXROOT) then
      det = MAXROOT
      return
    else
      tmpdet = tmpdet * lu(i,i)                          !compute det
    end if
  end do
  
  det = tmpdet

  return

End


! Gauss for multiple right hand sides ................................
Subroutine mgauss(n, k, mat, rmat, code)
    integer n,                &     !Dimension of system .............
            k                       !number of right hand sides ......
    real*8  mat(0:n-1,0:n-1), &     !original matrix .................
            rmat(0:n-1,0:n-1)       !Right hand sides/solutions ......
    integer code                    !Error code ......................
 !*====================================================================*
 !*                                                                    *
 !*  mgauss  finds the solution matrix x for the linear system         *
 !*  mat * x = rmat with an  n x n coefficient matrix mat and a        *
 !*  n x k matrix rmat of right hand sides. Here mat must be           *
 !*  nonsingular.                                                      *
 !*                                                                    *
 !* ------------------------------------------------------------------ *
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      n        integer;  (n > 0)                                    *
 !*               Dimension of mat.                                    *
 !*      k        integer k;  (k > 0)                                  *
 !*               number of right hand sides                           *
 !*      mat      n x n original system matrix                         *
 !*      rmat     matrix of right hand sides                           *
 !*                                                                    *
 !*   Output parameter:                                                *
 !*   ================                                                 *
 !*      rmat     solution matrix for the system.                      *
 !*               The input right hand sides are lost.                 *
 !*                                                                    *
 !*   Return value (code):                                             *
 !*   ===================                                              *
 !*      = 0      all ok                                               *
 !*      = 1      n < 1 or k < 1 or invalid input parameter            *
 !*      = 2      lack of memory (not used here)                       *
 !*      = 3      mat is numerically singular.                         *
 !*                                                                    *
 !* ------------------------------------------------------------------ *
 !*                                                                    *
 !*   subroutine used:                                                 *
 !*   ================                                                 *
 !*                                                                    *
 !*      gaudec:  LU decomposition of mat.                             *
 !*                                                                    *
 !*====================================================================*
  integer i, j, m, signd, rc, perm(0:n-1)
  real*8 lu(0:n-1,0:n-1), x(0:n-1), sum
  real*8 MACH_EPS

  if (n < 1.or.k < 1) then                      !Invalid parameter
    code = 1                   
    return
  end if

  MACH_EPS = 1.2D-16

  call gaudec (n, mat, lu, perm, signd, rc)     !compute factorization
                                                !in matrix lu         
  if (rc.ne.0.or.signd.eq.0) then               !if not possible      
    code = 3                                    !exit with code=3     
    return
  end if

  do m=0, k-1                          !Loop over the right hand sides
    do i=0, n-1                                      !Updating the b's
      x(i) = rmat(perm(i),m)
      do j=0, i-1
        x(i) = x(i) - lu(i,j) * x(j)
      end do
    end do

    do i=n-1, 0, -1                                 !back substitution
      sum = 0.D0
      do j=i+1, n-1
        sum = sum + lu(i,j) * x(j)
      end do

      if (DABS(lu(i,i)) < MACH_EPS) then     !invalid LU decomposition
        code = 2
        return
      end if
      x(i) = (x(i) - sum) / lu(i,i)
    end do

    do j=0, n-1                                           !Save result
      rmat(j,m) = x(j)
    end do
  end do

  code = 0
  return

End


End module fgauss

! ------------------------ END fgauss.f90 --------------------------