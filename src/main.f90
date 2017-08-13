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
integer  :: nd, np, fileUnit
real(RK) :: start, finish

namelist /general/ ndMin, ndMax, npMin, npMax, ntest

open(newunit=fileUnit, file='main.in', status='old')
read(fileUnit,nml=general)

do nd = 2,2
  do np = 3,20
    !write(*,*) 'DIMENSION: ', nd
    !write(*,*) 'NP:        ', np
    call cpu_time(start)
    call test_getMVUfromEllipsoids(nd, np, ntest)
    call cpu_time(finish)
    write(*,*) 'Job Time:  ', finish - start, ' s'
  end do
end do

end program test
