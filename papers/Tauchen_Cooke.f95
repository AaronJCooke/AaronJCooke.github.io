real function CDF(x)
real :: y,v,erf,sqrt2,pi,rn, dnml
data sqrt2/1.414213562373095/
data pi/3.141592653589793/
y = x/sqrt2
if (x<0) then
  y = -y
end if
v = 0.0
do n = 1,37
  rn = dfloat(n)
  v = v + exp(-rn*rn/25)/n*sin(2*n*y/5)
end do
!Use Strecok (1968) for erf approx
  v = v+y/5
  erf = 2*v/pi
  dnml = (1+erf)/2
if (x<0) then
  dnml = (1-erf)/2
  else if (x<-8.3d0) then 
  dnml = 0
  else if (x>8.3d0) then
  dnml = 1
end if
CDF = dnml
end function CDF

program Tauchen_Cooke
implicit none
integer :: i, s, n, k
real :: sigma_eps, m, sig_y, rho, y_state_0, y_state_n, w, CDF
real, allocatable :: y_state(:), P_state(:,:), sums(:)

s = 7
rho = .4
sigma_eps = .75! variance of normally distributed iid mean 0 epsilon in y(t+1)=rho*y(t) + epsilon

allocate(y_state(s),P_state(s,s))

!approximates an AR process with s state s by s dimensional Markov Chain

m = 3 !allowable standard deviation of the state variable from the mean
sig_y = sigma_eps * (1-(rho**2))**(-0.5) ! Standard deviation of the the state variable
y_state_0 = -1*sig_y*m ! The state with the smallest value
y_state_n = sig_y*m ! The state with the largest value
do i = 1,s! Equally spaced state space
  y_state(i) = y_state_0+(i-1)*(y_state_n-y_state_0)/(s-1)
end do
w = y_state(s) - y_state(s-1) ! The difference between consecutive state values

!print*,CDF(-3.0), CDF(0.0), CDF(3.0), 'should be approximately zero, 0.5, and 1'

do i = 1,s
  do k = 1,s
    if (k==1) then
      P_state(i,k)=CDF((y_state_0 - rho*y_state(i) + 0.5*w)/sigma_eps) ! see equation 3b of Tauchen 1986
    else if (k==s) then
      P_state(i,k)= 1.0 - CDF((y_state_n - rho*y_state(i) - 0.5* w)/sigma_eps) 
    else
      P_state(i,k)=CDF((y_state(k) - rho*y_state(i) + 0.5*w)/sigma_eps) - &
      CDF((y_state(k) - rho*y_state(i) - 0.5*w)/sigma_eps)
    end if
  end do
end do
open(11,file = 'tauchen.txt',status='replace') 
write(11,'(1F50.8)')P_state

y_state = exp(y_state)

open(11,file = 'ab.txt',status='replace') 
write(11,'(1F50.8)')y_state

allocate(sums(s))
	do i = 1,s
    sums(i) = sum(P_state(i,:))
    end do

print*,sums
end program Tauchen_Cooke