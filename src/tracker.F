      subroutine autotrans(nstep,nt,dt,mtime,zh,u0,ua,u3d,v0,va,v3d,wa,wten)
      use input
      implicit none

      integer, intent(in) :: nstep,nt
      real, intent(inout) :: dt
      double precision, intent(in) :: mtime
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa,wten

      integer :: i,j,k,kmin,kmax,ntlev
      real, dimension(ib:ie,jb:je,kb:ke) :: wa_new,wa_old
      real, dimension(ib:ie,jb:je,kb:ke) :: wa_fix
      real :: w1111,w1110,w2111,w0111,w1211,w1011
      real :: dwdt,dwdx,dwdy,bigu,bigv,speed
      real :: a,b,c,d,e,rdt,denom,velmax
      real :: uchange,vchange,oldumove,oldvmove
      real, save :: sumu,sumv
      integer, save :: ichnge

      data sumu /0.0/
      data sumv /0.0/
      data ichnge /0/

      uchange = 0.0
      vchange = 0.0
      rdt     = 1.0/dt
      velmax  = rdt*sqrt(dx*dy)
      a       = 0.0
      b       = 0.0
      c       = 0.0
      d       = 0.0
      e       = 0.0
      kmax    = 1 + chkdpth/dz

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        wa_fix(i,j,k) = 0.5*(wa(i,j,k)+wa(i,j,k+1))
        wa_old(i,j,k) = wa_new(i,j,k)
        wa_new(i,j,k) = 999.0
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        wa_new(i,j,k) = wa_fix(i,j,k)
      enddo
      enddo
      enddo

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        if(zh(i,j,k).ge.2000.0 .and. zh(i,j,k).le.5000.0)then
          w1111 = wa_new(i,j,k)
          w1110 = wa_old(i,j,k)
          w2111 = wa_new(i+1,j,k)
          w0111 = wa_new(i-1,j,k)
          w1211 = wa_new(i,j+1,k)
          w1011 = wa_new(i,j-1,k)

          dwdt = wten(i,j,k+1)

          dwdx = 0.25*(wa(i+1,j,k) + wa(i+1,j,k-1) - &
                       wa(i-1,j,k) + wa(i-1,j,k-1))*dx

          dwdy = 0.25*(wa(i,j+1,k) + wa(i,j+1,k-1) - &
                       wa(i,j-1,k) + wa(i,j-1,k-1))*dy

          a = a + dwdt*dwdx
          b = b + dwdx*dwdx
          c = c + dwdx*dwdy
          d = d + dwdt*dwdy
          e = e + dwdy*dwdy
        endif
      enddo
      enddo
      enddo


      denom = c*c - b*e

      if(denom.eq.0.0)then
        print *, 'Warning: denom = 0. bigu, bigv set to 0'
        bigu = 0.
        bigv = 0.
      else
        bigu = (a*e - c*d)/denom
        bigv = (b*d - c*a)/denom
        print *, 'bigu, bigv = ', bigu, bigv
      endif

      speed = sqrt(bigu**2 + bigv**2)

      if(speed.gt.velmax)then
        bigu = 0.
        bigv = 0.
        print *, 'Warning! speed= ', speed, ' exceeds velmax= ', velmax
        print *, 'Set bigu and bigv to zero.'
        print *, 'bigu, bigv = ', bigu, bigv
      endif

      sumu = sumu + bigu
      sumv = sumv + bigv

      ntlev = IFIX(twindow/dt)
      ichnge = MOD(nstep, ntlev)

      if(ichnge.eq.0)then

        uchange = sumu/ntlev
        vchange = sumv/ntlev

        oldumove = umove
        oldvmove = vmove

        print *, 'umove and vmove are being updated at time =', mtime

        print *, 'Old umove and vmove = ', oldumove, oldvmove
        print *, 'ntlev, uchange, vchange = ', ntlev, uchange, vchange

        umove = umove + uchange
        vmove = vmove + vchange

        print *, 'New umove and vmove will be  = ', umove, vmove

        print *, 'Updating U and V winds '
        do k=1,nk
          do j=jb,je
          do i=ib,ie+1
            u0(i,j,k) = u0(i,j,k) + (oldumove-umove)
            ua(i,j,k) = ua(i,j,k) + (oldumove-umove)
            u3d(i,j,k) = u3d(i,j,k) + (oldumove-umove)
          enddo
          enddo
          do j=jb,je+1
          do i=ib,ie
            v0(i,j,k) = v0(i,j,k) + (oldvmove-vmove)
            va(i,j,k) = va(i,j,k) + (oldvmove-vmove)
            v3d(i,j,k) = v3d(i,j,k) + (oldvmove-vmove)
          enddo
          enddo
        enddo

        sumu = 0.
        sumv = 0.

      endif
      return

      end subroutine autotrans