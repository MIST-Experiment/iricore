module mod

real, allocatable, dimension(:, :, :) :: outf_out

contains
subroutine iri_core(JF,JMAG,GLAT,GLON,IYYYY,MMDD,DHOUR,HEIBEG,HEIEND,HEISTP,OARR,datadir)
  logical, intent(in) :: jf(50)
  logical, intent(in) :: jmag
  integer, intent(in) :: iyyyy, mmdd
  real :: OUTF(20,1000)
  real, dimension(:), intent(in) :: glat, glon
!         real, allocatable, dimension(:, :, :), intent(out) :: outf_out
  real, intent(inout) :: OARR(100)
  real, intent(in) :: dhour, heibeg, heiend, heistp
  character(*), intent(in) :: datadir
  character(256) :: datadir1
  common /folders/ datadir1
  integer :: i, k, j

  datadir1 = datadir
  call read_ig_rz
  call readapf107


  do i = 1, SIZE(glat), 1
      call IRI_SUB(JF,JMAG,GLAT(i),GLON(i),IYYYY,MMDD,DHOUR,HEIBEG,HEIEND,HEISTP, OUTF, OARR, datadir)
      do j = 1, 20
          do k = 1, 1000
              outf_out(j, k, i) = OUTF(j, k)
          end do
      end do
  end do
return
end

end module mod
