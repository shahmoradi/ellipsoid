include '../lib/mod_Misc'
include '../lib/mod_Parameters'
include '../lib/mod_SpecialFunctions'
include '../lib/mod_Matrix'
include '../lib/mod_Statistics'
include './mod_MultivariateUniformSampler'

program test

use MultivariateUniformSampler, only: test_getMVUfromEllipsoids

implicit none
integer        :: nd, np, ntest
real           :: start, finish

ntest = 1000

do nd = 2,2
  do np = 3,20
    print *, 'DIMENSION: ', nd
    print *, 'NP:        ', np
    call cpu_time(start)
    call test_getMVUfromEllipsoids(nd, np, ntest)
    call cpu_time(finish)
    print *, 'Job Time:  ', finish - start, ' s'
  end do
end do

end program test
