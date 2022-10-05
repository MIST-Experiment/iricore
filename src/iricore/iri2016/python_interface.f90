subroutine iricore(jf, jmag, glat, glon, gsize, iyyyy, mmdd, dhour, heibeg, heiend, heistp, oarr, iri_res, datadir, ddsize)
  logical, intent(in) :: jf(50)
  logical, intent(in) :: jmag
  integer, intent(in) :: gsize
  real,  intent(in) :: glat(gsize), glon(gsize)
  integer, intent(in) :: iyyyy, mmdd
  real, intent(in) :: dhour, heibeg, heiend, heistp
  real, intent(inout) :: oarr(100)
  real, intent(inout) :: iri_res(20, 1000, gsize)
  integer, intent(in) :: ddsize
  character(ddsize), intent(in) :: datadir

  real :: outf(20,1000)
  character(256) :: datadir1
  common /folders/ datadir1
  integer :: i, k, j
  integer :: NITER

  datadir1 = datadir
  call read_ig_rz
  call readapf107
  NITER = SIZE(glat)

  do i = 1, NITER
      call IRI_SUB(jf, jmag, glat(i), glon(i), iyyyy, mmdd, dhour, heibeg, heiend, heistp, outf, oarr, datadir)
      do j = 1, 20
        do k = 1, 1000
          iri_res(j, k, i) = outf(j, k)
        end do
      end do
  end do
return
end
