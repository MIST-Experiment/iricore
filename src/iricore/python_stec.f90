subroutine stec(jf, jmag, lat, lon, heights, hsize, iyyyy, mmdd, dhour, oarr, iri_res, datadir, &
        ddsize)
    !                   ddsize, aap_, af107_, n_)
    logical, intent(in) :: jf(50)
    logical, intent(in) :: jmag
    integer, intent(in) :: hsize
    real, intent(in) :: lat(hsize), lon(hsize), heights(hsize)
    integer, intent(in) :: iyyyy, mmdd
    real, intent(in) :: dhour
    real, intent(inout) :: oarr(100)
    real, intent(inout) :: iri_res(20, 1000, hsize)
    integer, intent(in) :: ddsize
    character(ddsize), intent(in) :: datadir
    !  integer, intent(in) :: aap_(27000, 9)
    !  real, intent(in) :: af107_(27000, 3)
    !  integer, intent(in) :: n_

    integer :: aap(27000, 9)
    real :: af107(27000, 3)
    integer :: n
    common /apfa/aap, af107, n
    real :: outf(20, 1000)
    character(256) :: datadir1
    common /folders/ datadir1
    integer :: i, k, j
    integer :: NITER

  aap = aap_
  af107 = af107_
  n = n_
  datadir1 = datadir
  call read_ig_rz
    call readapf107
  NITER = SIZE(heights)


  do i = 1, NITER
      call IRI_SUB(jf, jmag, lat(i), lon(i), iyyyy, mmdd, dhour, heights(i), heights(i), 1., outf, oarr, datadir)
      do j = 1, 20
        iri_res(j, 1, i) = outf(j, 1)
      end do
  end do
return
end