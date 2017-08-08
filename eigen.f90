!Eigen_solver is a F90 module which simplifies the eigen value sovling probelm for those
!first time users of lapack.
!author: Jiuzhou Tang   (tangjiuzhou@iccas.ac.cn)
!history:   2017-8-8(first eition) 

!Credits: this code is greatly inspired by the open source code(iQIST: An open
!source continuous-time quantum Monte Carlo impurity solver toolkit), with one
!of its author a close friend of mine.

!For detailed information of eigen solver, please refer to the lapack's user guide on line.
!this eigen module will follow matlab's eigen funtion style, i.e., [val,vec]=eigen(A), which means you can specify the output of eigen vectors as optional, which is necessary for many scientific applications.

!!!Very importmant notification:
! eigen values output are not ordered in  procedure eigen_g_r/eigen_g_c 
! eigen values output in Aescedent order in  procedure eigen_s_r/eigen_h_c

module eigen_solver
  integer,parameter :: dp = kind(1.0d0) !double precision 
  
  interface eigen
    module procedure eigen_g_r  ! for general real matrix, by default, using double real
    module procedure eigen_s_r  ! for symmetric real matrix, by default,double complex
    module procedure eigen_g_c  ! for general complex matrix
    module procedure eigen_h_c  ! for hermite complex matrix
  end interface
  
  contains
  
  subroutine eigen_g_r(A,val,vec,ierr)
    implicit none
    real(dp), intent(in)  :: A(:,:) ! input matrix A
    integer :: mat_dim ! matrix dimension, assuming square matrix
    real(dp),intent(out) :: val(:) ! output array of eigen value
    real(dp), allocatable, intent(out), optional :: vec(:,:) ! list of output arrays of eigen vector,
! stored in a matrix form. i.e., [[v1],[v2],...v[n]].
    integer, intent(out), optional :: ierr ! error info, 0 for sucess, otherwise
! error, labeled with different integer.
    integer :: dim(2) ! result of shape(A)
! lapack local variables
    integer :: ndim  !matrix dim used in lapack subroutine
    integer :: istat ! status flag
    integer :: info  ! info of lapack computing subroutine,if info=0, the eigenvalues will
!output in descending order
    integer :: lwork ! the length of the array work, lwork >= max(1,4*ndim)
    real(dp), allocatable :: work(:) ! workspace array
    real(dp), allocatable :: wr(:)
    real(dp), allocatable :: wi(:) ! auxiliary real(dp) matrix: real and imaginary parts of eigenvalues
    real(dp), allocatable :: vr(:,:)
    real(dp), allocatable :: vl(:,:) ! auxiliary real(dp) matrix: left and right eigenvectors
    real(dp), allocatable :: vec_tmp(:,:) ! a temporary place holder for input matrix and output eigen vec.

    if (present(ierr)) ierr = 0 ! set ierr = 0 if present    
    dim = shape(A)
      if (dim(1) .eq. dim(2)) then
        mat_dim = dim(1)
      else 
        if ( present(ierr))  ierr = 1
        return
      endif
! allocate memory
    if(present(vec))  then 
      allocate(vec(mat_dim,mat_dim), stat=ierr)
    else
      allocate(vec_tmp(mat_dim,mat_dim), stat=ierr)
    endif

    ndim = mat_dim
! initialize lwork
    lwork = 4 * ndim
    allocate(work(lwork),   stat=istat)
    allocate(wr(ndim),      stat=istat)
    allocate(wi(ndim),      stat=istat)
    allocate(vr(ndim,ndim), stat=istat)
    allocate(vl(ndim,ndim), stat=istat)
    if ( istat /= 0 ) then
      if ( present(ierr))  ierr = 2
      print*, 'can not allocate enough memory'
      return
    endif
   
! call the computational subroutine: dgeev
    if ( present(vec))  then
      vec = A
      call DGEEV('N', 'V', ndim, vec, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, info) 
    else
      vec_tmp = A
      call DGEEV('N', 'N', ndim, vec_tmp, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, info) 
    endif
      
! check the status
    if ( info /= 0 ) then
      if ( present(ierr))  ierr = 3
      print*, 'error in lapack subroutine dgeev'
      return
    endif ! 
! copy eigenvalues and eigenvectors
    val(1:ndim) = wr(1:ndim)
    if(present(vec)) then
      vec(1:ndim,1:ndim) = vr(1:ndim,1:ndim)
    else
      vec_tmp(1:ndim,1:ndim) = vr(1:ndim,1:ndim)
    endif
! dealloate memory for workspace array
    if ( allocated(work) ) deallocate(work)
    if ( allocated(wr  ) ) deallocate(wr  )
    if ( allocated(wi  ) ) deallocate(wi  )
    if ( allocated(vr  ) ) deallocate(vr  )
    if ( allocated(vl  ) ) deallocate(vl  )
    if(.not. present(vec)) then
      if ( allocated(vec_tmp  ) ) deallocate(vec_tmp  )
    endif  

  end subroutine eigen_g_r

  subroutine eigen_g_c(A,val,vec,ierr)
    implicit none
    complex(dp), intent(in)  :: A(:,:) ! input matrix A
    integer :: mat_dim ! matrix dimension, assuming square matrix
    complex(dp),intent(out) :: val(:) ! output array of eigen value
    complex(dp), allocatable, intent(out), optional :: vec(:,:) ! list of output arrays of eigen vector,
! stored in a matrix form. i.e., [[v1],[v2],...v[n]].
    integer, intent(out), optional :: ierr ! error info, 0 for sucess, otherwise
! error, labeled with different integer.
    integer :: dim(2) ! result of shape(A)
! lapack local variables
    integer :: ndim  !matrix dim used in lapack subroutine
    integer :: istat ! status flag
    integer :: info  ! info of lapack computing subroutine,if info=0, the eigenvalues will
!output in descending order
    integer :: lwork ! the length of the array work, lwork >= max(1,2*ndim)
    complex(dp), allocatable :: work(:) ! workspace array
    complex(dp), allocatable :: rwork(:)
    complex(dp), allocatable :: vr(:,:)
    complex(dp), allocatable :: vl(:,:) ! auxiliary  matrix: left and right eigenvectors
    complex(dp), allocatable :: vec_tmp(:,:) ! a temporary place holder for input matrix and output eigen vec.

    if (present(ierr)) ierr = 0 ! set ierr = 0 if present    
    dim = shape(A)
      if (dim(1) .eq. dim(2)) then
        mat_dim = dim(1)
      else 
        if ( present(ierr))  ierr = 1
        return
      endif
! allocate memory
    if(present(vec))  then 
      allocate(vec(mat_dim,mat_dim), stat=ierr)
    else
      allocate(vec_tmp(mat_dim,mat_dim), stat=ierr)
    endif

    ndim = mat_dim
! initialize lwork
    lwork = 2 * ndim
    allocate(work(lwork),   stat=istat)
    allocate(rwork(ndim),      stat=istat)
    allocate(vr(ndim,ndim), stat=istat)
    allocate(vl(ndim,ndim), stat=istat)
    if ( istat /= 0 ) then
      if ( present(ierr))  ierr = 2
      print*, 'can not allocate enough memory'
      return
    endif
   
! call the computational subroutine: zgeev
    if ( present(vec))  then
      vec = A
      call ZGEEV('N', 'V', ndim, vec, ndim,val, vl, ndim, vr, ndim, work, lwork,rwork, info) 
    else
      vec_tmp = A
      call ZGEEV('N', 'N', ndim, vec_tmp, ndim, val, vl, ndim, vr, ndim, work, lwork,rwork, info) 
    endif
 
! check the status
    if ( info /= 0 ) then
      if ( present(ierr))  ierr = 3
      print*, 'error in lapack subroutine zgeev'
      return
    endif ! 
    if(present(vec)) then
      vec = vr
    else
      vec_tmp = vr
    endif
! dealloate memory for workspace array
    if ( allocated(work) ) deallocate(work)
    if ( allocated(rwork  ) ) deallocate(rwork  )
    if ( allocated(vr  ) ) deallocate(vr  )
    if ( allocated(vl  ) ) deallocate(vl  )
    if(.not. present(vec)) then
      if ( allocated(vec_tmp  ) ) deallocate(vec_tmp  )
    endif  

  end subroutine eigen_g_c


  subroutine eigen_s_r(A,p,val,vec,ierr)
    implicit none
    real(dp), intent(in)  :: A(:,:) ! input matrix A
    character*1, intent(in) :: p ! using upper triangular or lower triangular of
!the symmetric matrix, 'U' for upper, 'L' for lower
    integer :: mat_dim ! matrix dimension, assuming square matrix
    real(dp),intent(out) :: val(:) ! output array of eigen value
    real(dp), allocatable, intent(out), optional :: vec(:,:) ! list of output arrays of eigen vector,
! stored in a matrix form. i.e., [[v1],[v2],...v[n]].
    integer, intent(out), optional :: ierr ! error info, 0 for sucess, otherwise
! error, labeled with different integer.
    integer :: dim(2) ! result of shape(A)
! lapack local variables
    integer :: ndim  !matrix dim used in lapack subroutine
    integer :: istat ! status flag
    integer :: info  ! info of lapack computing subroutine,if info=0, the eigenvalues will
!output in ascending order
    integer :: lwork ! the length of the array work, lwork >= max(1,4*ndim)
    real(dp), allocatable :: work(:) ! workspace array
    real(dp), allocatable :: vec_tmp(:,:) ! a temporary place holder for input matrix and output eigen vec.

    if (present(ierr)) ierr = 0 ! set ierr = 0 if present    
    dim = shape(A)
      if (dim(1) .eq. dim(2)) then
        mat_dim = dim(1)
      else 
        if ( present(ierr))  ierr = 1
        return
      endif
! allocate memory
    if(present(vec))  then 
      allocate(vec(mat_dim,mat_dim), stat=ierr)
    else
      allocate(vec_tmp(mat_dim,mat_dim), stat=ierr)
    endif

    ndim = mat_dim
! initialize lwork
    lwork = 3 * ndim -1
    allocate(work(lwork),   stat=istat)
    if ( istat /= 0 ) then
      if ( present(ierr))  ierr = 2
      print*, 'can not allocate enough memory'
      return
    endif
   
! call the computational subroutine: dgeev
    if ( present(vec))  then
      vec = A
      call DSYEV('V', p, ndim, vec, ndim, val, work, lwork, info)
    else
      vec_tmp = A
      call DSYEV('N', p, ndim, vec_tmp, ndim, val, work, lwork, info)
    endif

      
! check the status
    if ( info /= 0 ) then
      if ( present(ierr))  ierr = 3
      print*, 'error in lapack subroutine dsyev'
      return
    endif ! 

! dealloate memory for workspace array
    if ( allocated(work) ) deallocate(work)
    if(.not. present(vec)) then
      if ( allocated(vec_tmp  ) ) deallocate(vec_tmp  )
    endif  

  
  end subroutine eigen_s_r


  subroutine eigen_h_c(A,p,val,vec,ierr)
    implicit none
    complex(dp), intent(in)  :: A(:,:) ! input matrix A
    character*1, intent(in) :: p ! using upper triangular or lower triangular of
!the symmetric matrix, 'U' for upper, 'L' for lower
    integer :: mat_dim ! matrix dimension, assuming square matrix
    complex(dp),intent(out) :: val(:) ! output array of eigen value
    complex(dp), allocatable, intent(out), optional :: vec(:,:) ! list of output arrays of eigen vector,
! stored in a matrix form. i.e., [[v1],[v2],...v[n]].
    integer, intent(out), optional :: ierr ! error info, 0 for sucess, otherwise
! error, labeled with different integer.
    integer :: dim(2) ! result of shape(A)
! lapack local variables
    integer :: ndim  !matrix dim used in lapack subroutine
    integer :: istat ! status flag
    integer :: info  ! info of lapack computing subroutine,if info=0, the eigenvalues will
!output in ascending order
    integer :: lwork,lrwork ! the length of the array work
!! the length of the array work and rwork
! lwork >= max(1,2*ndim-1), lrwork >= max(1,3*ndim-2)
    complex(dp), allocatable :: work(:) ! workspace array
    complex(dp), allocatable :: rwork(:) ! workspace array
    complex(dp), allocatable :: vec_tmp(:,:) ! a temporary place holder for input matrix and output eigen vec.

    if (present(ierr)) ierr = 0 ! set ierr = 0 if present    
    dim = shape(A)
      if (dim(1) .eq. dim(2)) then
        mat_dim = dim(1)
      else 
        if ( present(ierr))  ierr = 1
        return
      endif
! allocate memory
    if(present(vec))  then 
      allocate(vec(mat_dim,mat_dim), stat=ierr)
    else
      allocate(vec_tmp(mat_dim,mat_dim), stat=ierr)
    endif

    ndim = mat_dim
! initialize lwork
    lwork = 2 * ndim - 1
    lrwork = 3 * ndim - 2
    allocate(work(lwork),   stat=istat)
    allocate(rwork(lrwork), stat=istat)
    if ( istat /= 0 ) then
      if ( present(ierr))  ierr = 2
      print*, 'can not allocate enough memory'
      return
    endif
   
! call the computational subroutine: dgeev
    if ( present(vec))  then
      vec = A
      call ZHEEV('V', p, ndim, vec, ndim, val, work, lwork, rwork,info)
    else
      vec_tmp = A
      call ZHEEV('N', p, ndim, vec_tmp, ndim, val, work, lwork, rwork,info)
    endif

      
! check the status
    if ( info /= 0 ) then
      if ( present(ierr))  ierr = 3
      print*, 'error in lapack subroutine zheev'
      return
    endif ! 

! dealloate memory for workspace array
    if ( allocated(work) ) deallocate(work)
    if ( allocated(rwork) ) deallocate(rwork)
    if(.not. present(vec)) then
      if ( allocated(vec_tmp  ) ) deallocate(vec_tmp  )
    endif  

  
  end subroutine eigen_h_c




end module eigen_solver












