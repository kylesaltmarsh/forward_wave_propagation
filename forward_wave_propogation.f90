!--------------------------------------------------------------------------------------------------------
!  Kyle Saltmarsh   20741743
!--------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------------

program forwardmodel
  implicit none
  real(8) :: ds, srcval, dt, t, tf
  integer :: zdim, xdim, count
  real(8), dimension(:,:), allocatable :: v, uprev, unext, unow, usave, d2u
!initialising variables
ds = 10.d0
zdim = floor(3000.d0/ds)
xdim = floor(6000.d0/ds)
tf = 1.6d0
dt = 0.001
t = dt
count = 0
open(unit=2,file="shot.dat")
!allocating memory
allocate(v(xdim+2,zdim+2),unow(xdim+2,zdim+2),uprev(xdim+2,zdim+2),unext(xdim+2,zdim+2),d2u(xdim+2,zdim+2))
!create velocity model
call vmodel(xdim,zdim,ds,v)
open(unit=1, file="v.dat")
open(unit = 3, file="wavefield4.dat")
write(1,*)v
write(*,*)v(301,10)
write(*,*)v(301,202)
!intialise wavefield t=-0.1
uprev = 0.d0
unext = 0.d0
!wavefield for t=0 (only consists of source term)
unow = 0.d0
call source(0.d0,srcval)
unow(floor(xdim/2.d0),2) =  srcval
write(2,*)uprev(2:xdim+1,2)
write(2,*)unow(2:xdim+1,2)
!master loop over time
do while(t.le.tf)
   call source(t, srcval)
   call stencil(unow,d2u,ds,xdim,zdim)
   unext = (v**2)*(dt**2)*d2u + 2*unow - uprev
   !adding in source term
   unext(floor(xdim/2.d0),2) =unext(floor(xdim/2.d0),2)+srcval*(v(floor(xdim/2.d0),2)**2)*dt**2
   uprev = unow
   unow = unext
   t = t + dt
   write(2,*)unow(2:xdim+1,2)
   if(mod(count,50)==0) then
      write(3,*)unow
      write(*,*)count
   endif
   !write(3,*)unow
   count = count+1
end do

end program forwardmodel

!--------------------------------------------------------------------------------------------------------

!subroutine to generate velocity model
subroutine vmodel(xdim,zdim,ds,v)
  integer, intent(in) :: xdim, zdim
  real(8), intent(in) :: ds
  real(8), dimension(xdim+2,zdim+2), intent(out) :: v
  integer :: index
  v=0.d0
  index = ceiling(500.d0/ds)+1
  write(*,*)index
  v(2:xdim+1,2:index) = 1500.d0
  v(2:xdim+1,index:zdim+1) = 2000.d0
 !write(*,*)index
end subroutine vmodel

!--------------------------------------------------------------------------------------------------------

!subroutine to compute laplacian of wavefield
subroutine stencil(u,d2u,ds,xdim,zdim)
  integer, intent(in) :: xdim, zdim
  real(8), dimension(xdim+2,zdim+2), intent(in) :: u
  real(8), dimension(xdim+2,zdim+2), intent(out) :: d2u
  real(8), intent(in) :: ds
  integer :: i, j
d2u = 0.d0
  do i=2,xdim+1
     do j=2,zdim+1
        d2u(i,j) = -4*u(i,j) + u(i,j+1) + u(i+1,j) + u(i-1,j) + u(i,j-1)
     end do
  end do
  d2u = d2u/(ds**2)
end subroutine stencil

!--------------------------------------------------------------------------------------------------------

!subroutine to evaluate the source
subroutine source(t,f)
  real(8), intent(in) :: t
  real(8), intent(out) :: f
  real(8) :: A, fo, alpha, pi
  pi = 4.d0*datan(1.d0)
  fo = 20.d0
  A = 1.d0
  alpha = -30.7d0
  f = A*sin(2*pi*fo*t)*exp(alpha*t)
end subroutine source

!--------------------------------------------------------------------------------------------------------
! End of code
!--------------------------------------------------------------------------------------------------------
