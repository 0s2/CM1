  MODULE sounding_module

  implicit none

  private
  public :: sharppy_reader

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sharppy_reader(nsnd,p_sfc,th_sfc,qv_sfc,zsnd,thsnd,qvsnd,usnd,vsnd)
      use input
      use constants

      implicit none

      integer, intent(out) :: nsnd
      real, intent(out) :: p_sfc, th_sfc, qv_sfc
      real, dimension(:), intent(out) :: zsnd,thsnd,qvsnd,usnd,vsnd
      real, dimension(:), allocatable :: pres,hght,temp,dwpt,wdir,wspd


      integer :: k,header,io
      real :: flag
      character(len=8)   :: eof  = '%END%'    !end of file
      character(len=8)   :: bad  = '-9999.00' !bad data flag
      integer, parameter :: nmax = 1000000    !should be enough!

      allocate( pres(nmax) )
      allocate( hght(nmax) )
      allocate( temp(nmax) )
      allocate( dwpt(nmax) )
      allocate( wdir(nmax) )
      allocate( wspd(nmax) )

      open(unit=40,file='input_sounding',status='old',err=520)

      ! read header of file:
      do header = 1,6
        read(40,*,err=521) 
      enddo

      ! now, read entire sounding until end of file is discovered
      nsnd=1
      do k=1,nmax
        read(40,*,end=445,err=521,iostat=io) pres(k+1),hght(k+1),temp(k+1),dwpt(k+1),wdir(k+1),wspd(k+1)
        if(io.ne.0)then
          exit
        else
          nsnd=nsnd+1
        endif
      enddo
445   continue
      if(dowr) write(outfile,*) '  Found ',(nsnd-1),'  levels'
      if(dowr) write(outfile,*)
      close(unit=40)

      ! check to see if there is any bad data in the sounding
      read(bad,*) flag
      do k=1,nsnd
        if((pres(k).eq.flag).or.(hght(k).eq.flag).or.(temp(k).eq.flag).or. &
           (dwpt(k).eq.flag).or.(wdir(k).eq.flag).or.(wspd(k).eq.flag))then
          goto 522
        endif
      enddo

      call calc_sndparms(pres, hght, temp, dwpt, wdir, wspd, &
                         zsnd, thsnd, qvsnd, usnd, vsnd, nsnd)

      p_sfc  = pres(2)
      th_sfc = thsnd(2)
      qv_sfc = qvsnd(2)

      deallocate( pres )
      deallocate( hght )
      deallocate( temp )
      deallocate( dwpt )
      deallocate( wdir )
      deallocate( wspd )

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(dowr) write(outfile,*) 'Leaving SHARPPY_READER'

      return

    !--------------------------------------------------------------

520   print *
      print *,'  520: error opening input_sounding file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

521   print *
      print *,'  521: error reading input_sounding file '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

522   print *
      print *,'  522: error reading input_sounding file '
      print *,'  incomplete or missing data levels '
      print *,'    ... stopping cm1 ... '
      print *
      call stopcm1

    !--------------------------------------------------------------

      return
      end subroutine sharppy_reader


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine calc_sndparms(pres,hght,temp,dwpt,wdir,wspd, &
	                           z,th,qv,u,v,nsnd)
      use constants
      implicit none

      integer, intent(in) :: nsnd

      real, dimension(nsnd), intent(inout)  :: pres,hght,temp,dwpt,wdir,wspd
      real, dimension(nsnd), intent(out) :: z,th,qv,u,v
      real, dimension(nsnd) :: qvpre,exner
      integer :: k
      real, parameter :: p0    = 1.0e3
      real, parameter :: kappa = 0.28562982892500527
      do k=1,nsnd
        qvpre(k) = 6.11 *(10**(7.5 * dwpt(k)/(237.7+dwpt(k))))
		wspd(k) = wspd(k)*0.514444
      enddo

      do k=1,nsnd
        z(k)  = hght(k)-hght(2)
        th(k) = (temp(k)+273.15)*((p0/pres(k))**kappa)
        qv(k) = (621.97*qvpre(k)/(pres(k)-qvpre(k)))/1000.0
        u(k)  = -wspd(k) * sin(wdir(k)*pi/180.0)
        v(k)  = -wspd(k) * cos(wdir(k)*pi/180.0)
      enddo

      z(1)  = 0.0
      th(1) = 0.0
      qv(1) = 0.0
      u(1)  = 0.0
      v(1)  = 0.0

      end subroutine calc_sndparms


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  END MODULE sounding_module