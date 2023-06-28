MODULE Mat_inv_GlobVars

double precision, dimension(:), allocatable ::  uxx   !real matrix (dim_jacob)
double precision, dimension(:,:), allocatable ::  JacobMatrx   !real matrix (dim_jacob x dim_jacob)
  double precision, dimension(:,:), allocatable ::  JacobMatrxCopy  !copy of matrix JacobMatrx
  double precision, dimension(:,:), allocatable ::  inv_Jacob   !real matrix (dim_jacob x dim_jacob)
  integer, allocatable ::  row_permutat(:)  !integer vector (dim_jacob)
  double precision, dimension(:), allocatable ::  matrx_multip   !real matrix (dim_jacob)
  double precision, dimension(:), allocatable ::  S_new, S_old,S_t0	!!real matrix (dim_jacob)
  double precision, dimension(:), allocatable :: F_numerator
  double precision :: dt
  integer :: dim_jacob,CODE_mtrx
  integer :: dim2_ode
  


  
  END MODULE Mat_inv_GlobVars