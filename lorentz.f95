program lrnz
implicit none
integer:: i
integer,parameter::n=10000
real*8,dimension(1:3):: y0
real*8,dimension(n+1)::tout
real*8,dimension(1:3,n+1)::f
real*8::ti,tmax,h,r

write(*,*) 'Enter the value of r: '
read(*,*) r

!initial condition
do i=1,3
  y0(i)=1.0d0
end do

ti=0.0d0
tmax=25.0d0

!time step
h= tmax/n

!plot curve
call rkit(y0,ti,h,r,n,tout,f)
open(12, file="lorentz_curve.plt")
write(12,*)'variables ="t","x","y","z"'

do i=1,n+1
  write(12,*) tout(i),f(1,i),f(2,i),f(3,i)
  
end do
close(12)

end

!array of differential equations/functions

subroutine funcs(yin,t,r,f)
implicit none
real*8,dimension(1:3)::yin
real*8::t,r
real*8,dimension(1:3)::f

f(1) = 10.0d0 *(yin(2)- yin(1))
f(2) = (r*yin(1)) - yin(2) - (yin(1)*yin(3))
f(3) = (yin(1)*yin(2))-((8.0d0/3.0d0)*yin(3))

end

subroutine rk4(yin,ti,h,r,f)
implicit none
real*8,dimension(1:3)::yin
real*8::h,r,ti
real*8,dimension(1:3)::f
real*8,dimension(1:3)::k1,k2,k3,k4,yout

call funcs(yin,ti,r,yout)
k1 = h*yout

call funcs(yin+k1/2.0d0,ti+h/2.0d0,r,yout)
k2 = h*yout

call funcs(yin+k2/2.0d0,ti+h/2.0d0,r,yout)
k3 = h*yout

call funcs(yin+k3,ti+h,r,yout)
k4 = h*yout

f = yin + ((k1+2.0d0*(k2+k3)+k4)/6.0d0)

end

subroutine rkit(y0,ti,h,r,n,tout,f)
implicit none
integer::i,n
real*8,dimension(1:3)::y0
real*8,dimension(1:3,n+1)::f
real*8::h,ti,r
real*8,dimension(n+1)::tout

tout(1) = ti 
f(:,1) = y0(1)

do i = 2,n+1
tout(i) = tout(i-1) + h
call rk4(f(:,i-1),tout(i-1),h,r,f(:,i))
end do

end
