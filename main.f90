module timer
real(8) :: t1,t2
contains

subroutine tic()
  implicit none
  call cpu_time(t1)
end subroutine tic

subroutine toc()
  implicit none
  call cpu_time(t2)
  print*,"Time Taken -->", real(t2-t1)
end subroutine toc
end module timer




program main
  use eigen_solver
  use timer
  implicit none
  integer,parameter :: ndim=1000
  real(dp) :: mat(ndim,ndim)
  real(dp) :: val(ndim)
  complex(dp) :: zmat(ndim,ndim)
  complex(dp) :: zval(ndim)
  !real(dp) :: vec(ndim,ndim)
  real(dp),allocatable :: vec(:,:)
  complex(dp),allocatable :: zvec(:,:)
  integer :: i,j,k,ierr
! generate a randm matrix

    do i = 1,ndim
      do j = 1,ndim
        mat(j,i) = i*j*1.0d0 + sin(j*i*1.0d0) ! an example of symetric matrix
        zmat(j,i) = cmplx(mat(j,i),0)
      enddo
    enddo

    allocate(vec(ndim,ndim))
    allocate(zvec(ndim,ndim))

   ! call eigen(mat,val,vec,ierr)
   ! if (ierr .eq. 0) then 
   !   print*, "val=:,",val
   ! else
   !   print*,"failed to solve eigen value,error info: ",ierr 
   ! endif
   !call eigen(zmat,zval,zvec,ierr)
   ! if (ierr .eq. 0) then 
   !   print*, "zval=:,",zval
   ! else
   !   print*,"failed to solve complex eigen value,error info: ",ierr 
   ! endif
   !call eigen(mat,'U',val,vec,ierr)
   ! if (ierr .eq. 0) then 
   !   print*, "sym val=:,",val
   ! else
   !   print*,"failed to solve sym real mat eigen value,error info: ",ierr 
   ! endif

   call tic()
   call eigen(mat,'U',val,vec,ierr)
   call toc()
    if (ierr .eq. 0) then 
      print*, " val=:,",val(1)
    else
      print*,"failed to solve hermite mat eigen value,error info: ",ierr 
    endif

end










