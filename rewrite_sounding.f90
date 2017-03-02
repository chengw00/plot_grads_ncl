	integer, parameter :: nz=56

	real, dimension(nz) :: pressure, height, tempc, u, v, td

        real, dimension(nz) :: windspeed, winddir
 
	open(unit=10, file='sounding.gdat', form='unformatted',  &
             access='direct', recl=4*nz, status='old')

        nrec=1
        read(10,rec=nrec) pressure
        nrec=nrec+1
        read(10,rec=nrec) height
        nrec=nrec+1
        read(10,rec=nrec) tempc
        nrec=nrec+1
        read(10,rec=nrec) u
        nrec=nrec+1
        read(10,rec=nrec) v
        nrec=nrec+1
        read(10,rec=nrec) td
        nrec=nrec+1

        ! ======= calculate wind speed =====
        windspeed=sqrt(u*u+v*v)

        do k=1,nz
         call calc_winddir(u(k),v(k),winddir(k),1,1,1)
        enddo

        call writeb3d(nz,1,1,pressure,'all.dat'//char(0))
        call writeb3d(nz,1,1,height,'all.dat'//char(0))
        call writeb3d(nz,1,1,tempc,'all.dat'//char(0))
        call writeb3d(nz,1,1,windspeed,'all.dat'//char(0))
        call writeb3d(nz,1,1,winddir,'all.dat'//char(0))
        call writeb3d(nz,1,1,td,'all.dat'//char(0))

        write(99,1000) 'zdef ', nz, 'levels', (pressure(k), k=1,nz)
1000    format(a4,1x,i3,1x,a6,1x,100(f6.1,1x))

	end

! ---------------------------
! subroutine to calculate meteorological wind direction in deg
! 0 deg = N
!
! input: u - zonal wind
!        v - meridional wind
!
! output: dir - meteorological wind direction (deg)

       subroutine calc_winddir(u,v,dir,nx,ny,nz)

       real u(nx,ny,nz), v(nx,ny,nz), dir(nx,ny,nz)

       do k=1,nz
        do j=1,ny
         do i=1,nx

          if (u(i,j,k).eq.0.) u(i,j,k)=1.e-11
          if (v(i,j,k).eq.0.) v(i,j,k)=1.e-11

          dir(i,j,k)=atan(v(i,j,k)/u(i,j,k))

          if ( (u(i,j,k).ge.0.).and.(v(i,j,k).ge.0.) ) then
           dir(i,j,k)=270.-57.29578*dir(i,j,k)
          elseif ( (u(i,j,k).ge.0.).and.(v(i,j,k).le.0.) ) then
           dir(i,j,k)=270.+57.29578*abs(dir(i,j,k))
          elseif ( (u(i,j,k).le.0.).and.(v(i,j,k).le.0.) ) then
           dir(i,j,k)=90.-57.29578*abs(dir(i,j,k))
          elseif ( (u(i,j,k).le.0.).and.(v(i,j,k).ge.0.) ) then
           dir(i,j,k)=90.+57.29578*abs(dir(i,j,k))
          endif

         enddo
        enddo
       enddo

       return
       end
