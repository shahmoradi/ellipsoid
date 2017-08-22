include '../lib/mod_Misc'
include '../lib/mod_Parameters'
include '../lib/mod_SpecialFunctions'
include '../lib/mod_Matrix'
include '../lib/mod_Statistics'
include './mod_MultivariateUniformSampler'

program test

use Parameters, only: RK
use MultivariateUniformSampler, only: test_getMVUfromEllipsoids

implicit none

integer  :: ndMin, ndMax, npMin, npMax, ntest
integer  :: nd, np, fileUnit, stepsize
real(RK) :: start, finish

namelist /general/ ndMin, ndMax, npMin, npMax, ntest

open(newunit=fileUnit, file='main.in', status='old')
read(fileUnit,nml=general)

do nd = ndMin,ndMax,3
  stepsize = nd
  write(*,*) 'DIMENSION: ', nd
  do np = nd+1,nd*15,stepsize
    write(*,*) 'NP:        ', np
    call cpu_time(start)
    call test_getMVUfromEllipsoids(nd, np, ntest)
    call cpu_time(finish)
    write(*,*) 'Job Time:  ', finish - start, ' s'
  end do
end do

end program test
