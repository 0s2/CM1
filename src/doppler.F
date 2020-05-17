  MODULE doppler_module

  implicit none

  private
  public :: doppler_velocity

  CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine doppler_velocity(xh,yh,zh,ua,va,wa,qdiag,velocity)
      use input
      use constants

      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: velocity
      real, intent(in), dimension(ibdq:iedq,jbdq:jedq,kbdq:kedq,nqdiag) :: qdiag
      real, dimension(ib:ie,jb:je,kb:ke) :: u,v,w

      integer :: i,j,k


      real :: delx,dely,refz,rhofact,s1,s2,vt
      real :: azim,dist,delz,eleva,range,dhdr,dsdr
      real :: uazmrad,vazmrad

      real, parameter :: zfrez  = 3000.      ! freezing level (m)
      real, parameter :: zice   = 8000.      ! level above which entirely ice
      real, parameter :: h0     = 7000.
      real, parameter :: denom  = 1.0/(zice-zfrez)
      real, parameter :: elvmin = -0.65
      real, parameter :: elvmax = 19.95

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        u(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
        v(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
        w(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
      enddo
      enddo
      enddo

      do j=1,nj
      do i=1,ni
        delx=xh(i)-xh(radarx)
        dely=yh(j)-yh(radary)
        dist=sqrt(delx*delx+dely*dely)
        do k=1,nk
          delz=zh(i,j,k)-zh(i,j,1)
          call beamelv(delz,dist,eleva,range)
          if(range.gt.rngmin.and.range.lt.rngmax)then
            if(eleva.gt.elvmin.and.eleva.lt.elvmax)then
              if(qdiag(i,j,k,qd_dbz).gt.thdbz)then
                refz=10.**(0.1*qdiag(i,j,k,qd_dbz))
                rhofact=exp(0.4*zh(i,j,k)/h0)
                if(zh(i,j,k).lt.zfrez)then
                  vt=2.6*(refz**0.107)*rhofact
                elseif(zh(i,j,k).lt.zice)then
                  s1=(zice-zh(i,j,k))*denom
                  s2=2.*(zh(i,j,k)-zfrez)*denom
                  vt=s1*2.6*(refz**0.107)*rhofact + s2
                else
                  vt=2.0
                endif

                call dhdrange(eleva,range,dhdr)
                dsdr=SQRT(AMAX1(0.,(1.-dhdr*dhdr)))
                uazmrad=delx/dist
                vazmrad=dely/dist
                velocity(i,j,k)=u(i,j,k)*uazmrad*dsdr                &
                               +v(i,j,k)*vazmrad*dsdr                &
                               +(w(i,j,k)-vt)*dhdr
              endif
            endif
          endif
        enddo
      enddo
      enddo

      return

      end subroutine doppler_velocity

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine beamelv(height,sfcrng,elvang,range)
      use constants

      implicit none
      real, intent(in) :: height
      real, intent(in) :: sfcrng
      real, intent(out) :: elvang
      real, intent(out) :: range

      double precision, parameter :: eradius  = 6371200.0d0
      double precision, parameter :: frthrde  = (4.0d0*eradius/3.0d0)
      double precision, parameter :: rad2deg  = 180.0d0/pi

      double precision :: elvrad,hgtdb,rngdb,drange

      if(sfcrng.gt.0.0d0)then

        hgtdb=frthrde+dble(height)
        rngdb=dble(sfcrng)/frthrde
      
        elvrad = atan((hgtdb*cos(rngdb) - frthrde)/(hgtdb * sin(rngdb)))
        drange = (hgtdb*sin(rngdb))/cos(elvrad)
        elvang=rad2deg*elvrad
        range=drange
      else
        elvang=90.0
        range=height
      endif
      return

      end subroutine beamelv

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dhdrange(elvang,range,dhdr)
      use constants

      implicit none
      real, intent(in) :: elvang
      real, intent(in) :: range
      real, intent(out) :: dhdr

      double precision, parameter :: eradius  = 6371200.0d0
      double precision, parameter :: frthrde  = (4.0d0*eradius/3.0d0)
      double precision, parameter :: eighthre = (8.0d0*eradius/3.0d0)
      double precision, parameter :: fthsq    = frthrde*frthrde
      double precision, parameter :: deg2rad  = pi/180.0d0

      double precision :: sinelv,dhdrdb,drange

      drange=dble(range)
      sinelv=sin(deg2rad*dble(elvang))
      dhdrdb = (drange+frthrde*sinelv)/                                     &
             sqrt(drange*drange + fthsq + eighthre*drange*sinelv)
      dhdr = dhdrdb

      return
      end subroutine dhdrange

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  END MODULE doppler_module