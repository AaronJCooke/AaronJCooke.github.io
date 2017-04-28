program valuefunctioniteration_cooke
implicit none
real, parameter :: b = .99, d = .06, a = .33, s = 1.8
integer :: grid_length, inc
real :: kmin, kinc, kmax, k, kp, c, delta, e
real, allocatable :: kgrid(:), value_1(:), value_2(:), value_opt(:)
real :: k_ss, c_ss, y_ss
integer :: i = 1, cont, iter, index_k, index_kp

inc = 100 !number of grid points in each dimension

!compute steady state
k_ss = (1/a*(1/b+d-1))**(1/(a-1))
c_ss = a*(k_ss)**(a-1)-d
y_ss = k_ss**a

print *, 'Steady State values'
print *, 'Output: ', y_ss, 'Capital: ', k_ss, 'Consumption: ', c_ss

kmax = k_ss*1.1
kmin = k_ss*.9
kinc = (kmax-kmin)/inc
e = kinc

grid_length = (kmax-kmin)/kinc+1
allocate(kgrid(grid_length), value_1(grid_length), value_2(grid_length), value_opt(grid_length))

do while (i<=grid_length) !assigning grid
  kgrid(i) = kmin + (i-1)*inc
  i = i+1
end do

value_1 = 0*kgrid !initial guess
cont = 1
iter = 1
do while (cont < 1000)
  iter = iter + 1
  do index_k = 1, grid_length !capital_t grid
    k = kgrid(index_k)
    do index_kp = 1, grid_length !capital_t+1 grid
      kp = kgrid(index_kp)
      c = k**a+(1-d)*k-kp
      if (c>0) then
        value_2(index_kp) = (c**(1-s))/(1-s)+b*value_1(index_kp)
      else
        value_2(index_kp) = -1000000000
      end if 
    enddo
    value_opt(index_k) = maxval(value_2)
  end do
  
 delta = sum(abs(value_opt-value_1))
 print*, 'Iteration =' ,iter-1, '  Absolute Difference', delta
 if(delta < e) then
  print*, 'Convergence reached'
  cont = cont+1000
 end if
 cont = cont+1
 value_1 = value_opt
end do
open(12,file='d.txt',status='replace')
write(12,'(1F50.8)') delta
open(12,file='vo.txt',status='replace')
write(12,'(1F50.8)') value_opt
endprogram valuefunctioniteration_cooke