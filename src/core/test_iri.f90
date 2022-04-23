logical :: jf(50)
integer, parameter :: jmag = 0
integer :: iyyyy, mmdd
real :: glat, glon, dhour
real:: alt_km_range(3)

character(:), allocatable :: datadir

real :: oarr(100), outf(20,1000)
integer :: i

character(256) :: datadir1
common /folders/ datadir1

jf = .true.
jf(4:6) = .false.
jf(22:23) = .false.
jf(26) = .true.  ! jf(26) == jf(8) == .true. for foF2
jf(28:30) = .false.
jf(33:35) = .false.

glat = 65.
glon = -147.5
alt_km_range(1) = 500.
alt_km_range(2) = 500.
alt_km_range(3) = 1.
iyyyy = 2012
mmdd = 815
dhour = 6.
datadir = "/home/lap1dem/dev-python/iri_core/data/"
datadir1 = datadir

call read_ig_rz
call readapf107

do i = 1, 1000, 1
    call IRI_SUB(JF,JMAG,GLAT,GLON,IYYYY,MMDD,DHOUR,HEIBEG,HEIEND,HEISTP, OUTF, OARR, datadir)
end do


end program
