! -------- program to diagnose wrf output --------
!          calculates PV, vorticity and absolute vorticity
!          geopotential height, relative humidity
!
!          compile with cape.f write3d.o
!
!
!          2D fields output for web display (file: web.gdat)
!
!          Pressure-level fields are written into sequential binary
!          for subsequent re-writing into GRIB
!          (file: wrf.plev.fort.bin.dat)
!
!          Pressure-level fields are also written into direct
!          access binary for GrADS (file: all.dat)
!                          
      parameter (ptop = 5000)  ! model top in pa

      ! ===== d01 ========
      parameter (nx=279, ny=194, nz=36)  ! dimension of wrf original output in
                                         ! zeta_p (terrain-following) coordinate

      ! ===== d02 ========
      !parameter (nx=651, ny=483, nz=36)  ! dimension of wrf original output in
                                         ! zeta_p (terrain-following) coordinate

      ! ==== d03 =======
      !parameter (nx=513, ny=486, nz=36)  ! dimension of wrf original output in
                                         ! zeta_p (terrain-following) coordinate


      parameter (nzint=26)             ! number of interpolated levels
      parameter (nzint_select = 26)     ! number of vertical levels selected for output
      parameter (nzetap_select = 25)    ! number of zeta levels to output 

      parameter (dx=27000)
      !parameter (dx=9000)
      !parameter (dx=2700)             ! grid spacing in meters		

      ! --------- NCAR graphics parameters --------- !
      parameter (jprj = 3)             ! 1=Stereographic, 3=Lambert Conformal Conic
      parameter (plat = 10.)            ! these are the WRF grid projection parameters
      parameter (plon = 100.)
      parameter (rota = 40.) 
      ! ---------------------------------------------! 

      parameter (udef=-999999)        ! undefined value

!      character file2d*23, file3d*23

      ! ---- input variables ------------
      ! xlat: latitude    \
      ! xlon: longitude     2d
      ! xmpf: map factor  
      ! terr: terrain (m) /
      !
      ! press_hyd: basic state (hydrostatic) pressure (pa)  \
      ! press_p: perturbation pressure (Pa)
      ! qhydromet: cloud water mixing ratio (kg/kg)
      ! qice: ice water mixing ratio (kg/kg)
      ! qrain: rain water mixing ratio (kg/kg)
      ! qsnow: snow water mixing ratio (kg/kg)
      ! qv: water vapor mixing ratio (kg/kg)
      ! u: u-component of wind (m/s)                         3d in zeta_p      
      ! v: v-component of wind (m/s)                         coordinate
      ! w: vertical motion (m/s)                           
      ! z: height (m)   
      ! tempc: temperature (deg c)                          /                           
   
      real xlat(nx,ny), xlon(nx,ny), xmpf(nx,ny), terr(nx,ny)

      real press_hyd(nx,ny,nz), press_p(nx,ny,nz), press_total(nx,ny,nz)

      real qhydromet(nx,ny,nz), qhydromet_liq(nx,ny,nz), 
     *     qcloud(nx,ny,nz), qice(nx,ny,nz), qrain(nx,ny,nz), 
     *     qsnow(nx,ny,nz), qv(nx,ny,nz)

      real u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz),
     *     z(nx,ny,nz), tempc(nx,ny,nz)

      real zeta_full(nz+1), zeta_half(nz)

      ! ----- other input variables --------
      ! qv2:       mixing ratio at 2-m AGL (kg/kg) 
      ! rain_conv: subgrid-scale accumulated precip (mm)
      ! rain_grid: grid-scale accumulated precip (mm)
      ! pt2:       potential temperature at 2 m AGL (K)
      ! u10:       x-wind at 10 m AGL (m/s)
      ! v10:       y-wind at 10 m AGL (m/s) 
      ! slpwrf:    SLP from WRF (hPa)
      
      real qv2(nx,ny), rain_conv(nx,ny), 
     *     rain_grid(nx,ny), pt2(nx,ny),
     *     u10(nx,ny), v10(nx,ny), slpwrf(nx,ny),
     *     pblh(nx,ny), flux_lh(nx,ny), flux_sh(nx,ny),      
     *     filler(nx,ny)

      ! ------ intermediate variables -----------
      ! f: coriolis parameter (1/s)
      ! psfc: surface pressure (hpa)
      ! p_out: pressure levels to be interpolated to (hpa)
      
      real f(nx,ny), psfc(nx,ny)

      real p_out(nzint)

      real qv_zint(nx,ny,nzint)

      real exn(nx,ny,nz), theta(nx,ny,nz) 

      real p_in(nz), dumin(nz), vt1d(nz), 
     *     dumout(nzint), dum2d(nx,ny), dum3d(nx,ny,nz),
     *     rh(nx,ny,nz)

      real uvgrid(2)

      ! -------- output variables ----------
      ! utrop: u-component of wind on tropopause (m/s)  \
      ! vtrop: v-component of wind on tropopause (m/s)    2d
      ! ptrop: pressure on tropopause (hpa)             
      ! ppw:   precipitable water (mm)
      ! tempc2: temperature in deg C at 2-m AGL                  
      ! ept2:  equivalent potential temperature 
      !        at 2-m AGL (K)
      ! rh2:   relative humidity at 2-m AGL (%)         /
      ! 
      !
      !
      ! u_zint:       u-component of wind (m/s)  \
      ! v_zint:       v-component of wind (m/s)
      ! w_zint:       w-component of wind (m/s)
      ! dir_zint:     meteo wind direction (deg)
      ! dewptc_zint:  dew point temp (deg C)
      ! tempc_zint:   temperature (deg c)           3d in pressure levels 
      ! relhum_zint:  relative humidity (%) 
      ! absvort_zint: absolute vorticity (1/s)     
      ! vort_zint:    relative vorticity (1/s)      
      ! pv_zint:      potential vorticity (pvu)   /
      !
      ! 

      real u_zint(nx,ny,nzint), v_zint(nx,ny,nzint),
     *     w_zint(nx,ny,nzint),
     *     dir_zint(nx,ny,nzint), dewptc_zint(nx,ny,nzint),  
     *     tempc_zint(nx,ny,nzint),
     *     relhum_zint(nx,ny,nzint), z_zint(nx,ny,nzint),
     *     qhydromet_zint(nx,ny,nzint),
     *     qhydromet_liq_zint(nx,ny,nzint)

      real absvort_zint(nx,ny,nzint), relvort_zint(nx,ny,nzint), 
     *     pv_zint(nx,ny,nzint)

      real slpmm5(nx,ny), cape(nx,ny), cin(nx,ny), 
     *     utrop(nx,ny), vtrop(nx,ny), ptrop(nx,ny),
     *     ppw(nx,ny), tempc2(nx,ny), ept2(nx,ny), rh2(nx,ny),
     *     hydromet_path(nx,ny), hydromet_liq_path(nx,ny), 
     *     cloud_top_tempc(nx,ny) 

      ! ====== arrays for output at selected levels =======
      real u_zint_select(nx,ny,nzint_select),  
     *     v_zint_select(nx,ny,nzint_select),
     *     w_zint_select(nx,ny,nzint_select),
     *     z_zint_select(nx,ny,nzint_select),
     *     tempc_zint_select(nx,ny,nzint_select),
     *     relhum_zint_select(nx,ny,nzint_select),
     *     qhydromet_zint_select(nx,ny,nzint_select),
     *     qhydromet_liq_zint_select(nx,ny,nzint_select)

      real p_select(nzint_select)


      ! ================
      ! dummy variable for output
      real dummy_zetap(nx,ny,nzetap_select)

      ! ===== 2D average ===========
      real hydromet_path_cloudy_avg,  
     *     pblh_avg, 
     *     flux_lh_avg,
     *     flux_sh_avg,
     *     cfrac_avg                  

      ! ===========
      !data p_select /1000, 900, 850, 700, 500, 250/
      data p_select /1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
     *             600, 550, 500, 450, 400, 350, 300, 250, 200, 150,
     *             100,  70,  50,  30,  20,  10/

! ---- pressure levels ----

      data p_out /1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 
     *             600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 
     *             100,  70,  50,  30,  20,  10/ 

! ------ zeta coordinates (full levels) ------------

      data zeta_full /
     *    1.000000, 0.996200, 0.989737, 0.982460, 0.974381,
     *    0.965422, 0.955498, 0.944507, 0.932347, 0.918907,
     *    0.904075, 0.887721, 0.869715, 0.849928, 0.828211,
     *    0.804436, 0.778472, 0.750192, 0.719474, 0.686214,
     *    0.650339, 0.611803, 0.570656, 0.526958, 0.480854,
     *    0.432582, 0.382474, 0.330973, 0.278674, 0.226390,
     *    0.175086, 0.132183, 0.096211, 0.065616, 0.039773,
     *    0.018113, 0.000000/

! -----read wrf lat/long/map factor/terrain ---

       open(unit=30, 
     *   file='terr_mapfc_latlon_wrf_realtime.gdat',
     *       form='unformatted', access='direct',
     *       recl=4*nx*ny, status='old')


       nrec=1
       read(30,rec=nrec) terr
       nrec=nrec+1
       read(30,rec=nrec) xmpf
       nrec=nrec+1
       read(30,rec=nrec) xlat 
       nrec=nrec+1
       read(30,rec=nrec) xlon
       close(30)

c ---- read wrf output ----

       open(unit=40,
     *   file='wrfout_3d_realtime.gdat',
     *       form='unformatted', access='direct',
     *       recl=4*nx*ny*nz, status='old')

       nrec=1
       read(40,rec=nrec) u
       nrec=nrec+1
       read(40,rec=nrec) v
       nrec=nrec+1
       read(40,rec=nrec) w
       nrec=nrec+1
       read(40,rec=nrec) tempc
       nrec=nrec+1
       read(40,rec=nrec) press_total ! in hPa
       nrec=nrec+1
       read(40,rec=nrec) z
       nrec=nrec+1
       read(40,rec=nrec) qv
       nrec=nrec+1
       read(40,rec=nrec) qhydromet
       nrec=nrec+1
       read(40,rec=nrec) qhydromet_liq
       nrec=nrec+1
       read(40,rec=nrec) press_hyd  ! in Pa
       close(40)

       ! ======= define perturbation pressure in Pa =============!
       do k=1,nz
        do j=1,ny
         do i=1,nx
           press_p(i,j,k)=100*press_total(i,j,k)-press_hyd(i,j,k)
         enddo
        enddo
       enddo
       ! =================================================!

       open(unit=50,
     *   file='wrfout_2d_realtime.gdat',
     *       form='unformatted', access='direct',
     *       recl=4*nx*ny, status='old')

       nrec=1
       read(50,rec=nrec) qv2
       nrec=nrec+1
       read(50,rec=nrec) rain_conv
       nrec=nrec+1
       read(50,rec=nrec) rain_grid
       nrec=nrec+1
       read(50,rec=nrec) slpwrf
       nrec=nrec+1
       read(50,rec=nrec) pt2
       nrec=nrec+1
       read(50,rec=nrec) u10
       nrec=nrec+1
       read(50,rec=nrec) v10
       nrec=nrec+1
       read(50,rec=nrec) psfc
       nrec=nrec+1
       read(50,rec=nrec) pblh
       nrec=nrec+1
       read(50,rec=nrec) flux_lh
       nrec=nrec+1
       read(50,rec=nrec) flux_sh
       nrec=nrec+1

       close(50)

! ------- filler array ----------------------

        do j=1,ny
         do i=1,nx
          if (mod(i,2).eq.0) then
           filler(i,j)=-1.e-9
          else
           filler(i,j)=-5.e-9
          endif
         enddo
        enddo

c ------ calculate zeta at half levels -----

	do k=1,nz
         zeta_half(k)=0.5*(zeta_full(k)+zeta_full(k+1))
	enddo

c ------ calculate f (coriolis parameter) -----

	do i=1,nx
	 do j=1,ny
	  f(i,j)=2.0*2.0*(4.0*atan(1.0))/(24.*3600.)*
     *           sin(4.0*atan(1.0)*xlat(i,j)/180.)
         enddo
        enddo

! ======= find cloud top temperature ======

        call calc_cltop_tempc(nx,ny,nz,tempc,
     *                     qhydromet,cloud_top_tempc,udef)

! ------ calculate precipitable water ----------

	do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=press_p(i,j,k)+press_hyd(i,j,k)
            dumin(k)=qv(i,j,k)
           enddo         

           call precip_water(dumin,p_in,nz,ppw(i,j))

         enddo
        enddo

! ======= total hydrometeor path =======

        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=press_p(i,j,k)+press_hyd(i,j,k)
            dumin(k)=qhydromet(i,j,k)
           enddo

           call precip_water(dumin,p_in,nz,hydromet_path(i,j))

         enddo
        enddo

        ! ====== find average cloudy hydromet_path ======
        hydromet_path_cloudy_avg = 0
        cfrac_avg = 0
        ncount = 0

        do j=1,ny
         do i=1,nx
          if (hydromet_path(i,j).ge.10.e-3) then
           hydromet_path_cloudy_avg = hydromet_path_cloudy_avg 
     *                              + hydromet_path(i,j)
           cfrac_avg = cfrac_avg + 1
           ncount = ncount + 1
          endif
         enddo
        enddo 

        if (real(ncount).gt.0.5) then
         hydromet_path_cloudy_avg = hydromet_path_cloudy_avg/
     *                              real(ncount)
         cfrac_avg = cfrac_avg/real(nx*ny)
        endif

        ! ====== find average of SH, LH and PBLH ========
        pblh_avg = 0
        flux_lh_avg = 0
        flux_sh_avg = 0

        ncount = 0

        do j=1,ny
         do i=1,nx
          ncount = ncount + 1
          pblh_avg = pblh_avg + pblh(i,j)
          flux_lh_avg = flux_lh_avg + flux_lh(i,j)
          flux_sh_avg = flux_sh_avg + flux_sh(i,j)
         enddo
        enddo

        if (real(ncount).gt.0.5) then
         pblh_avg = pblh_avg/real(ncount)
         flux_lh_avg = flux_lh_avg/real(ncount)
         flux_sh_avg = flux_sh_avg/real(ncount)
        endif

        write(201,500) cfrac_avg,
     *                 hydromet_path_cloudy_avg,
     *                 pblh_avg,             
     *                 flux_lh_avg,
     *                 flux_sh_avg

500     format(5(e15.7,1x))
        ! ==============================================
        ! ====== liquid hydrometeor path =========

        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=press_p(i,j,k)+press_hyd(i,j,k)
            dumin(k)=qhydromet_liq(i,j,k)
           enddo

           call precip_water(dumin,p_in,nz,hydromet_liq_path(i,j))

         enddo
        enddo

c ------ calculate CAPE and CIN using subroutines from Greg Thompson ------

        do k=1,nz
	 do j=1,ny
          do i=1,nx
           theta(i,j,k)=(tempc(i,j,k)+273.16)* 
     *           (1.e5/(press_p(i,j,k)+press_hyd(i,j,k)))**(0.286)
            ! --- Exner function ----- !
           exn(i,j,k)=1004.*
     *        ((press_p(i,j,k)+press_hyd(i,j,k))/1.e5)**(0.286)
          enddo
         enddo
        enddo

        call stabcalc(nz,nx,ny,z,theta,exn,
     +       qv,cape,cin,udef)
       
c ------ calculate sea-level pressure from MM5 (GRAPH) package --------

	do j=1,ny
         do i=1,nx
        
          ! ------ need to flip array upside down to use SEAPRS_0 subroutine ----- 
          do k=1,nz
           kk=nz-k+1
           p_in(kk)=(press_p(i,j,k)+press_hyd(i,j,k))/100
           dumin(kk)=tempc(i,j,k)+273.16
          enddo   
          
          ! ------ estimate surface temperature using 6.5 deg/1000 m lapse rate ---- !
          t_sfc=tempc(i,j,1)
     *        +(z(i,j,1)-terr(i,j))*(6.5/1000.)+273.16

          !call SEAPRS_0(dumin,p_in,terr(i,j),psfc(i,j),
     *    !         t_sfc,1,1,nz,slpmm5(i,j))
          slpmm5(i,j)=undef
         enddo
        enddo

! ------ calculate RH ------------

        do j=1,ny
         do i=1,nx
          do k=1,nz
           p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
           call rhcalc(tempc(i,j,k),
     *                 qv(i,j,k),
     *                 p_in(k),
     *                 1, 1, 1,
     *                 rh(i,j,k))
          enddo
         enddo
        enddo

c ----- interpolate data --------
      
        ! ---- interpolate u ----------- ! 
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=u(i,j,k)    
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            u_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

        ! ---- interpolate v ----------- !
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=v(i,j,k)
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            v_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

        ! ---- interpolate w ----------- !
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=w(i,j,k)
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            w_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

        ! ---- interpolate qv ----------- !
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=qv(i,j,k)
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            qv_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

        ! ------ interpolate temperature ---- !

        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=tempc(i,j,k)+273.16
           enddo

           call interptemp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            tempc_zint(i,j,k)=dumout(k)-273.16
           enddo

         enddo
        enddo

        ! ------ interpolate height ---- !

        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            vt1d(k)=(tempc(i,j,k)+273.16)*(1+0.608*qv(i,j,k))
            dumin(k)=z(i,j,k)
           enddo

           zbot=terr(i,j)
           vt1000=(tempc_zint(i,j,1)+273.16)*
     *        (1+0.608*qv_zint(i,j,1))
           vt10=(tempc_zint(i,j,nzint)+273.16)*
     *        (1+0.608*qv_zint(i,j,nzint))

           call interpphgt(dumin,vt1d,p_in,zbot,nz,dumout,
     *           p_out,nzint,vt1000,vt10)

            do k=1,nzint
             z_zint(i,j,k)=dumout(k)
            enddo

         enddo
        enddo

        ! ---- interpolate qhydromet (total) ----------- !
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=qhydromet(i,j,k)
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            qhydromet_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

        ! ---- interpolate qhydromet (liquid) ----------- !
        do j=1,ny
         do i=1,nx

           do k=1,nz
            p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
            dumin(k)=qhydromet_liq(i,j,k)
           enddo

           call interpp(dumin,p_in,nz,dumout,p_out,nzint)

           do k=1,nzint
            qhydromet_liq_zint(i,j,k)=dumout(k)
           enddo

         enddo
        enddo

c ------ calculate relative humidity --------

        call rhcalc(tempc_zint,qv_zint,p_out,
     *      nx,ny,nzint,relhum_zint)

c ------ calculate pv -------------------------

	call pvcalc(u_zint,v_zint,tempc_zint,xmpf,
     *              f,p_out,dx,nx,ny,nzint,pv_zint)

c ------ calculate absolute vorticity (x 10^5) ----------

        print*, 'caclulate absolute vorticity in 10^5'
	call absvortcalc(u_zint,v_zint,xmpf,
     *          f,p_out,dx,nx,ny,nzint,absvort_zint)

        do k=1,nzint
         do j=1,ny
          do i=1,nx
           relvort_zint(i,j,k)=absvort_zint(i,j,k)-f(i,j)*1.e5
          enddo
         enddo
        enddo

c --------- find tropopause pressure, winds -----------------

        print*, 'find tropopause'
        call findtrop_pres(pv_zint,relhum_zint,
     *                     p_out,nx,ny,nzint,ptrop)

	do j=1,ny
         do i=1,nx

           if (ptrop(i,j).eq.udef) then

            utrop(i,j)=udef
            vtrop(i,j)=udef

           else
         
            do k=1,nz
             p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
             dumin(k)=u(i,j,k)
            enddo
            call interpp(dumin,p_in,nz,utrop(i,j),ptrop(i,j),1)

            do k=1,nz
             p_in(k)=(press_p(i,j,k)+press_hyd(i,j,k))/100
             dumin(k)=v(i,j,k)
            enddo
            call interpp(dumin,p_in,nz,vtrop(i,j),ptrop(i,j),1)

           endif
    
         enddo
        enddo

c -------- do some processing before output some 2D fields --------

        write(*,*) 'calculating 2-m, 10-m stuff'

        ! --- calculate RH at 2-m AGL ------- ! 

        do j=1,ny
         do i=1,nx
          dum_tempc=pt2(i,j)*((psfc(i,j)/1000.)**0.286)-273.16
          call rhcalc(dum_tempc,qv2(i,j),psfc(i,j),1,1,1,rh2(i,j))
         enddo
        enddo

        write(*,*) '2-m RH'
        ! -- calculate equivalent potential temperature at 2-m AGL -------- !

        do j=1,ny
         do i=1,nx
          dum_tempk=pt2(i,j)*((psfc(i,j)/1000.)**0.286)
          call THEQZH(dum_tempk,qv2(i,j),psfc(i,j),1,1,1,ept2(i,j),udef)
         enddo
        enddo

        write(*,*) '2-m EPT'

	! ------ calculate 2 m (AGL) temperature  --------

	do j=1,ny
         do i=1,nx
          tempc2(i,j)=pt2(i,j)*((psfc(i,j)/1000.)**0.286)
          tempc2(i,j)=tempc2(i,j)-273.16 
         enddo
        enddo	

        write(*,*) '2-m tempc'
        ! ------ calculating dewpoint ----------------
         do k=1,nzint
          do j=1,ny
           do i=1,nx
            press_temp=p_out(k)*100
            tempk_temp=tempc_zint(i,j,k)+273.16
            qv_min=amin1(qv_zint(i,j,k),rs(press_temp,tempk_temp))
            dewptc_zint(i,j,k)=td(press_temp,qv_min)-273.16
           enddo
          enddo
         enddo

        write(*,*) '2-m dewptc'

! ---------------------------------------------------------------------------
        ! ------ rotate winds to Earth-relative winds --------
        !call rotateuv(xlat,xlon,utrop,vtrop,nx,ny,1,udef)
        !call rotateuv(xlat,xlon,u10,v10,nx,ny,1,udef)
        !call rotateuv(xlat,xlon,u_zint,v_zint,nx,ny,nzint,udef)
        !call rotateuv(xlat,xlon,u,v,nx,ny,nz,udef)      ! rotate u,v before output
       
         pii = 3.14159265
         radians_per_degree = pii/180.
         tlat1 = plat
         tlat2 = rota
         clon  = plon 

         cone = 10**(cos(tlat1*radians_per_degree))  
     *         -10**(cos(tlat2*radians_per_degree))
         cone = cone/(10**(tan(45. -abs(tlat1/2.)*radians_per_degree)) -  
     *      10**(tan(45. -abs(tlat2/2.)*radians_per_degree))   )

        do j=1,ny
         do i=1,nx
          call dcomputeuvgrid(utrop(i,j),vtrop(i,j),uvgrid,      
     *       xlongca,xlongcb,xlon(i,j),xlat(i,j),                 
     *       clon,cone,radians_per_degree,1,1,1,1,1)
          utrop(i,j)=uvgrid(1)
          vtrop(i,j)=uvgrid(2)

          call dcomputeuvgrid(u10(i,j),v10(i,j),uvgrid,     
     *       xlongca,xlongcb,xlon(i,j),xlat(i,j),             
     *       clon,cone,radians_per_degree,1,1,1,1,1)
          u10(i,j)=uvgrid(1)
          v10(i,j)=uvgrid(2)

          do k=1,nzint
           call dcomputeuvgrid(u_zint(i,j,k),v_zint(i,j,k),uvgrid,
     *       xlongca,xlongcb,xlon(i,j),xlat(i,j),
     *       clon,cone,radians_per_degree,1,1,1,1,1)
           u_zint(i,j,k)=uvgrid(1)
           v_zint(i,j,k)=uvgrid(2)
          enddo


          do k=1,nz
           call dcomputeuvgrid(u(i,j,k),v(i,j,k),uvgrid,
     *       xlongca,xlongcb,xlon(i,j),xlat(i,j),
     *       clon,cone,radians_per_degree,1,1,1,1,1)
           u(i,j,k)=uvgrid(1)
           v(i,j,k)=uvgrid(2)
          enddo

         enddo
        enddo
                                                                                    
        ! ------ calculating meteorological wind direction --------
        !call winddir(u_zint,v_zint,dir_zint,nx,ny,nzint)

! ======== find levels to store output arrays ===========

       do j=1,ny
        do i=1,nx
         do k=1,nzint
          do kk=1,nzint_select
           if (int(p_out(k)).eq.int(p_select(kk))) then
            u_zint_select(i,j,kk)=u_zint(i,j,k)
            v_zint_select(i,j,kk)=v_zint(i,j,k)
            w_zint_select(i,j,kk)=w_zint(i,j,k)
            z_zint_select(i,j,kk)=z_zint(i,j,k)
            tempc_zint_select(i,j,kk)=tempc_zint(i,j,k)
            relhum_zint_select(i,j,kk)=relhum_zint(i,j,k)
            qhydromet_zint_select(i,j,kk)=qhydromet_zint(i,j,k)
            qhydromet_liq_zint_select(i,j,kk)=qhydromet_liq_zint(i,j,k)
            goto 1000
           endif
          enddo
1000      continue
         enddo
        enddo
       enddo

! -------- output in grads binary to file --------

        call writeb3d(1,ny,nx,terr,'all.dat'//char(0))               ! 1
        call writeb3d(1,ny,nx,ptrop,'all.dat'//char(0))              ! 2 
        call writeb3d(1,ny,nx,utrop,'all.dat'//char(0))              ! 3
        call writeb3d(1,ny,nx,vtrop,'all.dat'//char(0))              ! 4
        call writeb3d(1,ny,nx,slpmm5,'all.dat'//char(0))             ! 5
        call writeb3d(1,ny,nx,slpwrf,'all.dat'//char(0))             ! 6
        call writeb3d(1,ny,nx,psfc,'all.dat'//char(0))               ! 7
        call writeb3d(1,ny,nx,cape,'all.dat'//char(0))               ! 8
        call writeb3d(1,ny,nx,cin,'all.dat'//char(0))                ! 9
        call writeb3d(1,ny,nx,rain_conv,'all.dat'//char(0))          ! 10
        call writeb3d(1,ny,nx,rain_grid,'all.dat'//char(0))          ! 11
        call writeb3d(1,ny,nx,tempc2,'all.dat'//char(0))             ! 12
        call writeb3d(1,ny,nx,ept2,'all.dat'//char(0))               ! 13
        call writeb3d(1,ny,nx,rh2,'all.dat'//char(0))                ! 14
        call writeb3d(1,ny,nx,u10,'all.dat'//char(0))                ! 15
        call writeb3d(1,ny,nx,v10,'all.dat'//char(0))                ! 16
        call writeb3d(1,ny,nx,ppw,'all.dat'//char(0))                ! 17
        call writeb3d(1,ny,nx,hydromet_path,'all.dat'//char(0))      ! 18
        call writeb3d(1,ny,nx,hydromet_liq_path,'all.dat'//char(0))  ! 19
        call writeb3d(1,ny,nx,cloud_top_tempc,'all.dat'//char(0))    ! 20
        call writeb3d(1,ny,nx,pblh,'all.dat'//char(0))               ! 21
        call writeb3d(1,ny,nx,flux_lh,'all.dat'//char(0))            ! 22
        call writeb3d(1,ny,nx,flux_sh,'all.dat'//char(0))            ! 23


        call writeb3d(nzint,ny,nx,u_zint,'all.dat'//char(0))             ! 24
        call writeb3d(nzint,ny,nx,v_zint,'all.dat'//char(0))             ! 25
        call writeb3d(nzint,ny,nx,w_zint,'all.dat'//char(0))             ! 26
        call writeb3d(nzint,ny,nx,z_zint,'all.dat'//char(0))             ! 27
        call writeb3d(nzint,ny,nx,tempc_zint,'all.dat'//char(0))         ! 28
        call writeb3d(nzint,ny,nx,relhum_zint,'all.dat'//char(0))        ! 29
        call writeb3d(nzint,ny,nx,qhydromet_zint,'all.dat'//char(0))     ! 30
        call writeb3d(nzint,ny,nx,qhydromet_liq_zint,'all.dat'//char(0)) ! 31
        call writeb3d(nzint,ny,nx,pv_zint,'all.dat'//char(0))            ! 32
        call writeb3d(nzint,ny,nx,relvort_zint,'all.dat'//char(0))       ! 33
        call writeb3d(nzint,ny,nx,absvort_zint,'all.dat'//char(0))       ! 34

! ========== output the first five zeta_p levels ===========
        
        call fill_array(z,nx,ny,nz,dummy_zetap,nzetap_select)            ! 1
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *      'all_zetap.dat'//char(0))  

        call fill_array(press_total,nx,ny,nz,dummy_zetap,nzetap_select)  ! 2
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))  

        call fill_array(tempc,nx,ny,nz,dummy_zetap,nzetap_select)        ! 3
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(rh,nx,ny,nz,dummy_zetap,nzetap_select)           ! 4
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(u,nx,ny,nz,dummy_zetap,nzetap_select)            ! 5
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(v,nx,ny,nz,dummy_zetap,nzetap_select)            ! 6
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(w,nx,ny,nz,dummy_zetap,nzetap_select)            ! 7
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(qhydromet,nx,ny,nz,dummy_zetap,nzetap_select)    ! 8
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))

        call fill_array(qhydromet_liq,nx,ny,nz,dummy_zetap,
     *       nzetap_select)                       
        call writeb3d(nzetap_select,ny,nx,dummy_zetap,
     *       'all_zetap.dat'//char(0))                                   ! 9
! ========== output for storage ===========

        call writeb3d(1,ny,nx,terr,'all_2d.dat'//char(0))               ! 1
        call writeb3d(1,ny,nx,ptrop,'all_2d.dat'//char(0))              ! 2
        call writeb3d(1,ny,nx,utrop,'all_2d.dat'//char(0))              ! 3
        call writeb3d(1,ny,nx,vtrop,'all_2d.dat'//char(0))              ! 4
        call writeb3d(1,ny,nx,slpwrf,'all_2d.dat'//char(0))             ! 5
        call writeb3d(1,ny,nx,psfc,'all_2d.dat'//char(0))               ! 6
        call writeb3d(1,ny,nx,cape,'all_2d.dat'//char(0))               ! 7
        call writeb3d(1,ny,nx,cin,'all_2d.dat'//char(0))                ! 8
        call writeb3d(1,ny,nx,rain_conv,'all_2d.dat'//char(0))          ! 9
        call writeb3d(1,ny,nx,rain_grid,'all_2d.dat'//char(0))          ! 10
        call writeb3d(1,ny,nx,tempc2,'all_2d.dat'//char(0))             ! 11
        call writeb3d(1,ny,nx,ept2,'all_2d.dat'//char(0))               ! 12
        call writeb3d(1,ny,nx,rh2,'all_2d.dat'//char(0))                ! 13
        call writeb3d(1,ny,nx,u10,'all_2d.dat'//char(0))                ! 14
        call writeb3d(1,ny,nx,v10,'all_2d.dat'//char(0))                ! 15
        call writeb3d(1,ny,nx,ppw,'all_2d.dat'//char(0))                ! 16
        call writeb3d(1,ny,nx,hydromet_path,'all_2d.dat'//char(0))      ! 17
        call writeb3d(1,ny,nx,hydromet_liq_path,'all_2d.dat'//char(0))  ! 18
        call writeb3d(1,ny,nx,cloud_top_tempc,'all_2d.dat'//char(0))    ! 19
        call writeb3d(1,ny,nx,pblh,'all_2d.dat'//char(0))               ! 20
        call writeb3d(1,ny,nx,flux_lh,'all_2d.dat'//char(0))            ! 21
        call writeb3d(1,ny,nx,flux_sh,'all_2d.dat'//char(0))            ! 22

        call writeb3d(nzint_select,ny,nx,u_zint_select,                 ! 1
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,v_zint_select,                 ! 2
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,w_zint_select,                 ! 3
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,z_zint_select,                 ! 4
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,tempc_zint_select,             ! 5
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,relhum_zint_select,            ! 6
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,qhydromet_zint_select,         ! 7
     *       'all_3d.dat'//char(0))
        call writeb3d(nzint_select,ny,nx,qhydromet_liq_zint_select,     ! 8
     *       'all_3d.dat'//char(0))

      stop 
      end



c ---------- this subroutine interpolates a column of data in pressure ----
c            coordinate to specified pressure levels
c
c            a: column of data in the given pressure levels p(i)
c            ai: column of interpolated array to specified pressure levels ppi(i)
c
c            the pressure levels go from bottom-up.
c
        subroutine interpp(a,p,nz,ai,ppi,nzi)

        real a(nz), p(nz), ai(nzi), ppi(nzi)

        do 600 kk=1,nzi

         do 100 k=2,nz

          if (ppi(kk).gt.p(1)) then
           ai(kk)=a(1)
           go to 500
          elseif (ppi(kk).lt.p(nz)) then
           ai(kk)=a(nz)
           go to 500
          elseif (ppi(kk).eq.p(k)) then
           ai(kk)=a(k)
           go to 500
          elseif (ppi(kk).eq.p(k-1)) then
           ai(kk)=a(k-1)
           go to 500
          elseif ( (ppi(kk).lt.p(k-1)).and.(ppi(kk).gt.p(k)) ) then
           dp1=alog(ppi(kk)/p(k-1))
           dp2=alog(p(k)/ppi(kk))
           dp3=alog(p(k)/p(k-1))
           ai(kk)=( (a(k)-a(k-1))/dp3 )*dp1 + a(k-1)
           go to 500
          endif

100      continue

500      continue

600     continue

        return
        end

c ---------- this subroutine interpolates a column of temperature ------------------
c            to specified pressure levels 
c
c      input:
c        t: column of temperature (k) in the given pressure levels p(i)
c
c     output:
c       ai: column of interpolated virtual temperature to specified pressure levels ppi(i)
c
c            the pressure levels go from bottom-up.
c
        subroutine interptemp(t,p,nz,ai,ppi,nzi)

        real t(nz), p(nz), ai(nzi), ppi(nzi)

        r=287.0
        g=9.81
        gammadry=9.8/1000.0
        gammamoi=6.5/1000.0

        do 600 kk=1,nzi

          if ( ppi(kk).gt.p(1) ) then
           ai(kk)=t(1)*(ppi(kk)/p(1))**(r*gammamoi/g)
           go to 500

          elseif ( ppi(kk).lt.p(nz) ) then
           ai(kk)=t(nz)*(ppi(kk)/p(nz))**(r*gammamoi/g)
           go to 500

          endif


         do 100 k=2,nz

          if (ppi(kk).eq.p(k)) then
           ai(kk)=t(k)
           go to 500
          elseif (ppi(kk).eq.p(k-1)) then
           ai(kk)=t(k-1)
           go to 500
          elseif ( (ppi(kk).lt.p(k-1)).and.(ppi(kk).gt.p(k)) ) then
           dp1=alog(ppi(kk)/p(k-1))
           dp2=alog(p(k)/ppi(kk))
           dp3=alog(p(k)/p(k-1))
           ai(kk)=( (t(k)-t(k-1))/dp3 )*dp1 + t(k-1)
           go to 500
          endif

100      continue

500      continue

600     continue

        return
        end


c ---------- this subroutine interpolates a column of height in pressure ----
c            coordinate to specified pressure levels for height
c
c            a: column of height (m) in the given pressure levels p(i)
c           vt: column of virtual temp (k) in the given pressure levels p(i)
c            z: terrain
c       vt1000: 1000-hpa virtual temperature (k)
c         vt10:   10-hpa virtual temperature (k)
c
c            ai: column of interpolated array to specified pressure levels ppi(i)
c
c            the pressure levels go from bottom-up.
c
        subroutine interpphgt(a,vt,p,zbot,nz,ai,ppi,nzi,vt1000,vt10)

        real a(nz), vt(nz), p(nz), ai(nzi), ppi(nzi)
        real r, g

        r=287.0
        g=9.81

        do 600 kk=1,nzi

          if ( ppi(kk).gt.p(1) )  then
           tbar=0.5*(vt(1)+vt1000)
           ai(kk)=a(1) + (r*tbar/g)*alog(p(1)/ppi(kk))
           go to 500

          elseif ( ppi(kk).lt.p(nz) ) then
           tbar=0.5*(vt(nz)+vt10)
           ai(kk)=a(nz) + (r*tbar/g)*alog(p(nz)/ppi(kk))
           go to 500

          endif


         do 100 k=2,nz

          if (ppi(kk).eq.p(k)) then
           ai(kk)=a(k)
           go to 500
          elseif (ppi(kk).eq.p(k-1)) then
           ai(kk)=a(k-1)
           go to 500
          elseif ( (ppi(kk).lt.p(k-1)).and.(ppi(kk).gt.p(k)) ) then
           dp1=alog(ppi(kk)/p(k-1))
           dp2=alog(p(k)/ppi(kk))
           dp3=alog(p(k)/p(k-1))
           ai(kk)=( (a(k)-a(k-1))/dp3 )*dp1 + a(k-1)
c            dp1=ppi(kk)-p(k-1)
c            dp2=p(k)-ppi(kk)
c            dp3=p(k)-p(k-1)
c            ai(kk)=( (a(k)-a(k-1))/dp3 )*dp1 + a(k-1)

           go to 500
          endif

100      continue

500      continue

600     continue

        return
        end

! ---------- this subroutine calculates the relative humidity -------
!
!       input:
!         t: temperature (deg c)
!         q: water vapor mixing ratio (kg/kg)
!         p: pressure levels (hpa)
!
!      output:
!        rh: relative humidity (%)
!
        subroutine rhcalc(t,q,p,nx,ny,nz,rh)

        real t(nx,ny,nz), q(nx,ny,nz), rh(nx,ny,nz)
        real p(nz)

        do k=1,nz
         do j=1,ny
          do i=1,nx

           temp_k=t(i,j,k)+273.16              

           es=6.1078*
     *       exp(17.269*(temp_k-273.16)/(temp_k-35.86))      ! saturation vapor pressure
                                                             ! in hpa

           qs=0.622*es/(p(k)-es)                             ! saturation mixing ratio

           rh(i,j,k)=amax1(amin1(100.*q(i,j,k)/qs,100.),0.0) ! calculate rh in %

          enddo
         enddo
        enddo

        return
        end

! ---------------------------------------------------------------------------------
      ! ----- subroutine to calculate pv ---------------
      !
      ! input:
      !	    u: u-component of wind (m/s)
      !     v: v-component of wind (m/s)
      !     t: temperature (deg c)
      !   dmf: map factor
      !     p: pressure levels (hpa)
      !    ds: grid spacing (m)
      !
      ! output:
      !     pv: potential vorticity (pvu)

      subroutine pvcalc(u,v,t,dmf,f,p,ds,i1,j1,k1,pv)                      pvp.2
c                                                                          pvp.3
      real u(i1,j1,k1), v(i1,j1,k1), t(i1,j1,k1),                          pvp.4
     *          dmf(i1,j1), f(i1,j1),                                      pvp.5
     *          p(k1),                                                     pvp.6
     *          pv(i1,j1,k1)                                               pvp.7
c                                                                          pvp.13
c     ... everything is by slabs                                           pvp.14
c                                                                          pvp.15
      g=9.81                                                               pvp.16
      scale=-1.e6                                                          pvp.17
c                                                                          pvp.18
      do 5000 k=1,k1                                                       pvp.19
       do 5000 j=2,j1-1
 	do 5000 i=2,i1-1 	
c                                                                          pvp.20
c        ... compute vertical derivatives: du/dp, dv/dp, dtheta/dp         pvp.21
c                                                                          pvp.22
         if(k.eq.1) then                                                   pvp.23
            k0=k                                                           pvp.24
            k2=k+1                                                         pvp.25
         else if(k.eq.k1) then                                             pvp.26
            k0=k-1                                                         pvp.27
            k2=k                                                           pvp.28
         else                                                              pvp.29
            k0=k-1                                                         pvp.30
            k2=k+1                                                         pvp.31
         endif                                                             pvp.32

            dudp=(u(i,j,k0)-u(i,j,k2))                                     pvp.35
     *                  / (p(k0)-p(k2))                                    pvp.39
            dvdp=(v(i,j,k0)-v(i,j,k2))                                     pvp.40
     *                  / (p(k0)-p(k2))                                    pvp.44

            t1=(t(i,j,k0)+273.16)*(1000./p(k0))**0.286                     pvp.48
            t2=(t(i,j,k2)+273.16)*(1000./p(k2))**0.286                     pvp.49
            dtdp=(t1-t2)/(p(k0)-p(k2))                                     pvp.50
c                                                                          pvp.52
c        ... compute horizontal derivatives: dtheta/dx, dtheta/dy          pvp.53
c                                                                          pvp.54
         ds2r=1./(ds * 2.)                                                 pvp.55
            t1=(t(i+1,j  ,k)+273.16)*(1000./p(k))**0.286                   pvp.58
            t2=(t(i-1,j  ,k)+273.16)*(1000./p(k))**0.286                   pvp.59
            t3=(t(i  ,j+1,k)+273.16)*(1000./p(k))**0.286                   pvp.60
            t4=(t(i  ,j-1,k)+273.16)*(1000./p(k))**0.286                   pvp.61
            dtdx=(t1 - t2)*ds2r                                            pvp.62
            dtdy=(t3 - t4)*ds2r                                            pvp.63
c                                                                          pvp.67
c        ... compute slab absolute vorticity                               pvp.68
c                                                                          pvp.69
            u1=u(i,j-1,k)/dmf(i,j-1)
	    u2=u(i,j+1,k)/dmf(i,j+1)
	    v1=v(i-1,j,k)/dmf(i-1,j)
	    v2=v(i+1,j,k)/dmf(i+1,j)
            fx=f(i,j)                                                      pvp.80
            vor=dmf(i,j)**2 * ds2r*((v2-v1)-(u2-u1)) + fx                  pvp.81
c                                                                          pvp.84
c        ... group terms                                                   pvp.85
c                                                                          pvp.86
c --------- add map factors for horz derivative of theta ----
c --------- fixed the wrong sign on the horiz. components ------

c             pvhold(i,j,k)=g * scale /  100. *
              pv(i,j,k)=g * scale /100. * 
     *       (  vor * dtdp                                                 pvp.90
     *        - dvdp * dmf(i,j) * dtdx
     *        + dudp * dmf(i,j) * dtdy )
5000  continue                                                             pvp.93
c         
         call fillit(pv,i1,j1,k1,i1,j1,2,i1-1,2,j1-1)
c                                                                          pvp.100
      return                                                               pvp.101
      end


c ---------------------------------------------------------------------------------
      ! --------- subroutine to calculate absolute vorticity ------------
      !
      ! input:
      !     u: u-component of wind (m/s)
      !     v: v-component of wind (m/s)
      !     t: temperature (deg c)
      !   dmf: map factor
      !     p: pressure levels (hpa)
      !    ds: grid spacing (m)
      !
      ! output:
      !  vort: absolute vorticity (10^-5 s^-1)

      subroutine absvortcalc(u,v,dmf,f,p,ds,i1,j1,k1,vort)                 pvp.2
c                                                                          pvp.3
      real u(i1,j1,k1), v(i1,j1,k1),                                       pvp.4
     *          dmf(i1,j1), f(i1,j1),                                      pvp.5
     *          p(k1),                                                     pvp.6
     *          vort(i1,j1,k1)                                             pvp.7
c                                                                          pvp.13
c     ... everything is by slabs                                           pvp.14
c                                                                          pvp.15
c                                                                          pvp.18
      do 5000 k=1,k1                                                       pvp.19
       do 5000 j=2,j1-1
 	do 5000 i=2,i1-1 	
c                                                                          pvp.20
         ds2r=1./(ds * 2.)                                                 pvp.55
c                                                                          pvp.67
c        ... compute slab absolute vorticity                               pvp.68
c                                                                          pvp.69
            u1=u(i,j-1,k)/dmf(i,j-1)
	    u2=u(i,j+1,k)/dmf(i,j+1)
	    v1=v(i-1,j,k)/dmf(i-1,j)
	    v2=v(i+1,j,k)/dmf(i+1,j)
            fx=f(i,j)                                                      pvp.80
            vort(i,j,k)=
     *         dmf(i,j)**2 * ds2r*((v2-v1)-(u2-u1)) + fx                   pvp.81
c                                                                          pvp.84
c        ... group terms                                                   pvp.85
c                                                                          pvp.86

              vort(i,j,k)= vort(i,j,k) * 1.e5
5000  continue                                                             pvp.93
c         
         call fillit(vort,i1,j1,k1,i1,j1,2,i1-1,2,j1-1)
c                                                                          pvp.100
      return                                                               pvp.101
      end

c ----------------------------------------------------------------------------------
      subroutine fillit(f,ix,jx,kx,imx,jmx,ifirst,ilast,jfirst,jlast)      fillit.2
c                                                                          fillit.3
c     section  tools                                                       fillit.4
c     purpose  fill data out to imx,jmx from an interior domain            fillit.5
c                                                                          fillit.6
      dimension f(ix,jx,kx)                                                fillit.7
c                                                                          fillit.8
      do 1000 k=1,kx                                                       fillit.9
         do 300 j=jfirst,jlast                                             fillit.10
            do 100 i=1,ifirst-1                                            fillit.11
               f(i,j,k)=f(ifirst,j,k)                                      fillit.12
100         continue                                                       fillit.13
            do 200 i=ilast+1,imx                                           fillit.14
               f(i,j,k)=f(ilast,j,k)                                       fillit.15
200         continue                                                       fillit.16
300      continue                                                          fillit.17
         do 600 i=1,imx                                                    fillit.18
            do 400 j=1,jfirst-1                                            fillit.19
               f(i,j,k)=f(i,jfirst,k)                                      fillit.20
400         continue                                                       fillit.21
            do 500 j=jlast+1,jmx                                           fillit.22
               f(i,j,k)=f(i,jlast,k)                                       fillit.23
500         continue                                                       fillit.24
600      continue                                                          fillit.25
1000  continue                                                             fillit.26
      return                                                               fillit.27
      end

c ------------This subroutine rotates the wind in Lambert conformal projection --------------- 
c             (with only one standard parallel)
c             to true north and calculates the speed and meteorological wind
c             direction
c       Input: ugrd (Lambert grid u-component; in m/s)
c              vgrd (Lambert grid v-component; in m/s)
c
c       Output: ugrd
c               vgrd  (true lat/long wind components)


        SUBROUTINE ROTATEUV(XLAT,XLON,UGRD,VGRD,NX,NY,NZ,UDEF)
        REAL XLAT(NX,NY), XLON(NX,NY), UGRD(NX,NY,NZ), VGRD(NX,NY,NZ)

c ------- Parameters for the map projection ---------

      PARAMETER ( ROTCON_P   =  0.70710678    )
      PARAMETER ( LON_XX_P   = -114.0          )
C**  ROTCON_P          R  WIND ROTATION CONSTANT, = 1 FOR POLAR STEREO
C**                         AND SIN(LAT_TAN_P) FOR LAMBERT CONFORMAL
C**  LON_XX_P          R  MERIDIAN ALIGNED WITH CARTESIAN X-AXIS(DEG)
C**  LAT_TAN_P         R  LATITUDE AT LAMBERT CONFORMAL PROJECTION
C**                         IS TRUE (DEG)
      PARAMETER ( LAT_TAN_P  =  45.0          )

c ---------- rotate wind to true north ---

        DO 110 J=1,NY
         DO 110 I=1,NX
          ANGLE2 = ROTCON_P*(XLON(I,J)-LON_XX_P)*0.017453
          SINX2 = SIN(ANGLE2)
          COSX2 = COS(ANGLE2)

          DO 100 K=1,NZ
           IF ( (UGRD(I,J,K).EQ.UDEF).OR.(VGRD(I,J,K).EQ.UDEF) ) GOTO 90

           UT =  COSX2*UGRD(I,J,K)+SINX2*VGRD(I,J,K)
           VT = -SINX2*UGRD(I,J,K)+COSX2*VGRD(I,J,K)
           UGRD(I,J,K) = UT
           VGRD(I,J,K) = VT

90         CONTINUE

100       CONTINUE

110     CONTINUE

        RETURN
        END

c ---------------------------------------------------------------------------------
      SUBROUTINE SEAPRS_0(T,PP,TER,SFP,TS,IMX,JMX,KX,                            SEAPRS.25
     * SLP)                                                                      SEAPRS.26
C                                                                                SEAPRS.27
C     SECTION  DIAGNOSTIC                                                        SEAPRS.28
C     PURPOSE  COMPUTES SEA LEVEL PRESSURE FROM THE RULE                         SEAPRS.29
C              T1/T2=(P1/P2)**(GAMMA*R/G).                                       SEAPRS.30
C
C     *** LEVELS GO FROM TOP-DOWN ***
C                                                                                SEAPRS.31
C     INPUT       T        TEMPERATURE (Kelvin)                3D                SEAPRS.32
C                 TER      TERRAIN     (m)                     2D                SEAPRS.33
C                 SFP      SURFACE PRESSURE (hPa)              2D                SEAPRS.35
C                 IMX      DOT POINT DIMENSION N-S                               SEAPRS.37
C                 JMX      DOT POINT DIMENSION E-W                               SEAPRS.38
C                 KX       NUMBER OF VERTICAL LEVELS                             SEAPRS.39
C                                                                                SEAPRS.41
C     OUTPUT      SLP      SEA LEVEL PRESSURE (hPa)            2D                SEAPRS.42
C                                                                                SEAPRS.43
      DIMENSION T(IMX,JMX,KX), PP(IMX,JMX,KX),                                   SEAPRS.44
     *          PS(IMX,JMX)  ,SFP(IMX,JMX) ,                                     SEAPRS.45
     *          TER(IMX,JMX)                                                     SEAPRS.46
      DIMENSION PL(IMX,JMX),T0(IMX,JMX),TS(IMX,JMX),                             SEAPRS.47
     *          XKLEV(IMX,JMX)                                                   SEAPRS.48
      DIMENSION SLP(IMX,JMX)                                                     SEAPRS.49
      PARAMETER (R=287.04,G=9.8,GAMMA=6.5E-3)                                    SEAPRS.50
      PARAMETER (TC=273.16+17.5) ! T CRITICAL IN PSFC/PSLV                       SEAPRS.51
      PARAMETER (PCONST=100.)                                                    SEAPRS.52
C                                                                                SEAPRS.53
      LOGICAL L1,L2,L3,L4                                                        SEAPRS.54
C                                                                                SEAPRS.55
C                                                                                SEAPRS.57
C                                                                                SEAPRS.58
C     ... SEA LEVEL PRESSURE                                                     SEAPRS.59
C                                                                                SEAPRS.60
      XTERM=GAMMA*R/G                                                            SEAPRS.61
C                                                                                SEAPRS.62
C     ... COMPUTE PRESSURE AT PCONST MB ABOVE SURFACE (PL)                       SEAPRS.63
C                                                                                SEAPRS.64
      KUPTO=KX/2                                                                 SEAPRS.65
99    CONTINUE                                                                   SEAPRS.66
      DO 100 J=1,JMX                                                             SEAPRS.67
      DO 100 I=1,IMX                                                             SEAPRS.68
         PL(I,J)=SFP(I,J)-PCONST                                                 SEAPRS.69
         XKLEV(I,J)=0.                                                           SEAPRS.70
100   CONTINUE                                                                   SEAPRS.71
C                                                                                SEAPRS.72
C     ... FIND 2 LEVELS ON SIGMA SURFACES SURROUNDING PL AT EACH I,J             SEAPRS.73
C                                                                                SEAPRS.74
      DO 150 J=1,JMX                                                             SEAPRS.75
      DO 150 I=1,IMX                                                             SEAPRS.76
         DO 125 K=KX-1,KUPTO,-1                                                  SEAPRS.77
            XK=FLOAT(K)                                                          SEAPRS.78
            XKHOLD=XKLEV(I,J)                                                    SEAPRS.79
            XKLEV(I,J)=CVMGT(XK,XKHOLD,                                          SEAPRS.80
     *         (((PP(I,J,K)).LT.PL(I,J)) .AND.                                   SEAPRS.81
     *          ((PP(I,J,K+1)).GE.PL(I,J))))                                     SEAPRS.82
125      CONTINUE                                                                SEAPRS.83
         IF(XKLEV(I,J).LT.1.) THEN                                               SEAPRS.84
            PRINT *,'ERROR FINDING PRESSURE LEVEL ',PCONST,' MB ',               SEAPRS.85
     *              'ABOVE THE SURFACE'                                          SEAPRS.86
            PRINT *,'LAST K LEVEL =',KUPTO                                       SEAPRS.87
            IF(KUPTO.NE.1) THEN                                                  SEAPRS.88
               PRINT *,'TRYING AGAIN WITH KUPTO=1'                               SEAPRS.89
               KUPTO=1                                                           SEAPRS.90
               GOTO 99                                                           SEAPRS.91
            ELSE                                                                 SEAPRS.92
               PRINT *,'I,J=',I,J                                                SEAPRS.93
               PRINT *,'PL=',PL(I,J)                                             SEAPRS.94
               PRINT *,'PSFC=',SFP(I,J)                                          SEAPRS.95
               STOP                                                              SEAPRS.96
            END IF                                                               SEAPRS.97
         END IF                                                                  SEAPRS.98
150   CONTINUE                                                                   SEAPRS.99
C                                                                                SEAPRS.100
C     ... GET TEMPERATURE AT PL (TL), EXTRAPOLATE T AT SURFACE (TS)              SEAPRS.101
C         AND T AT SEA LEVEL (T0) WITH 6.5 K/KM LAPSE RATE                       SEAPRS.102
C                                                                                SEAPRS.103
      DO 200 J=1,JMX                                                             SEAPRS.104
      DO 200 I=1,IMX                                                             SEAPRS.105
         KLO=NINT(XKLEV(I,J))+1                                                  SEAPRS.106
         KHI=NINT(XKLEV(I,J))                                                    SEAPRS.107
         PLO=PP(I,J,KLO)                                                         SEAPRS.108
         PHI=PP(I,J,KHI)                                                         SEAPRS.109
         TLO=T(I,J,KLO)                                                          SEAPRS.110
         THI=T(I,J,KHI)                                                          SEAPRS.111
         TL=THI-(THI-TLO)*ALOG(PL(I,J)/PHI)/ALOG(PLO/PHI)                        SEAPRS.112
         TS(I,J)=TL*(SFP(I,J)/PL(I,J))**XTERM                                    SEAPRS.113
         TBAR=(TS(I,J)+TL)*0.5                                                   SEAPRS.114
         HL=TER(I,J)-R/G*ALOG(PL(I,J)/SFP(I,J))*TBAR                             SEAPRS.115
         T0(I,J)=TL+GAMMA*HL                                                     SEAPRS.116
200   CONTINUE                                                                   SEAPRS.117
C                                                                                SEAPRS.118
C     ... CORRECT SEA LEVEL TEMPERATURE IF TOO HOT                               SEAPRS.119
C                                                                                SEAPRS.120
      DO 400 J=1,JMX                                                             SEAPRS.121
      DO 400 I=1,IMX                                                             SEAPRS.122
         L1=T0(I,J).LT.TC                                                        SEAPRS.123
         L2=TS(I,J).LE.TC                                                        SEAPRS.124
         L3=.NOT.L1                                                              SEAPRS.125
         T0HOLD=T0(I,J)                                                          SEAPRS.126
         T0(I,J)=CVMGT(T0HOLD,                                                   SEAPRS.127
     *      CVMGT(TC,TC-0.005*(TS(I,J)-TC)**2,L2.AND.L3),                        SEAPRS.128
     *      L1.AND.L2)                                                           SEAPRS.129
400   CONTINUE                                                                   SEAPRS.130
C                                                                                SEAPRS.131
C     ... COMPUTE SEA LEVEL PRESSURE                                             SEAPRS.132
C                                                                                SEAPRS.133
      DO 600 J=1,JMX                                                             SEAPRS.134
      DO 600 I=1,IMX                                                             SEAPRS.135
         SLP(I,J)=SFP(I,J)*EXP(2.*G*TER(I,J)/(R*(TS(I,J)+T0(I,J))))              SEAPRS.136
600   CONTINUE                                                                   SEAPRS.137
      RETURN                                                                     SEAPRS.138
      END                                                                        SEAPRS.139

C -----------------------------------------------------------------------------------------

      FUNCTION CVMGT(X1,X2,L1)                                                   CVMGT.1
      LOGICAL L1                                                                 CVMGT.2
      IF (L1) THEN                                                               CVMGT.3
         CVMGT=X1                                                                CVMGT.4
      ELSE                                                                       CVMGT.5
         CVMGT=X2                                                                CVMGT.6
      ENDIF                                                                      CVMGT.7
      RETURN                                                                     CVMGT.8
      END                                                                        CVMGT.9

c -------- subroutine to find the dynamic tropopause pressure level --------
c
c input: pv: potential vorticity (PVU)
c        rh: relative humidity (%)
c        plev: pressure levels (hPa)
c    
c output:
c        ptrop: dynamic tropopause pressure (hPa)
 
	subroutine findtrop_pres(pv,rh,plev,nx,ny,nz,ptrop)

	parameter (pcrit=2)  ! value of PV representing tropopause in PVU
	parameter (udef=-999999) ! udefined value

	real pv(nx,ny,nz), rh(nx,ny,nz)
	real plev(nz)
	real ptrop(nx,ny)

	do j=1,ny
         do i=1,nx

          do k=nz,2,-1    ! search from top-down

           if ( (pv(i,j,k).ge.pcrit).and.
     *          (pv(i,j,k-1).le.pcrit) ) then

             dp1=(pcrit-pv(i,j,k-1))
             dp2=alog(plev(k)/plev(k-1))
             dp3=pv(i,j,k)-pv(i,j,k-1)

            ptrop(i,j)=plev(k-1)*exp(dp1*dp2/dp3)
            go to 100
           endif 
          enddo

          ptrop(i,j)=undef

100       continue

          if (ptrop(i,j).gt.600.) ptrop(i,j)=udef
         
         enddo
        enddo

        return
        end

! ---------- subroutine to calculate thetae ------------
!            algorithm from Bolton
!
! input:
!      T: temperature (K)
!      Q: water vapor mixing ratio (kg/kg)
!    PRS: pressure (hPa)
!
! output:
!      thetae: equivalent potential temperature (K)
!
        SUBROUTINE THEQZH(T,Q,PRS,IX,JX,KX,THETAE,UDEF)

        DIMENSION T(IX,JX,KX),Q(IX,JX,KX),THETAE(IX,JX,KX),
     *    PRS(IX,JX,KX)

      DO 40 K=1,KX
       DO 40 J=1,JX
        DO 40 I=1,IX

         IF (T(I,J,K).LT.4.) THEN
          T(I,J,K)=UDEF
          GO TO 30
         ENDIF

         IF(T(I,J,K).GT.273.16) THEN
          ES=6.11*EXP(19.84659-5418.12/T(I,J,K))
         ELSE
          ES=6.11*EXP(22.514-6.15E3/T(I,J,K))
         ENDIF

         QS=.622*ES/(PRS(I,J,K)-ES)
         Q(I,J,K)=AMAX1(1.E-10,Q(I,J,K))
         RHH=Q(I,J,K)/QS
         RHH=AMIN1(1.,RHH)

          IF (RHH.LE.0.) THEN
           PRINT*, 'RHH=',RHH
           PRINT*, I, J, K
           STOP
          ENDIF

         THETA=T(I,J,K)*(1000./PRS(I,J,K))**0.286
         FACT=1./(1./(T(I,J,K)-55.)-ALOG(RHH)/2840.)+55.
         THETAE(I,J,K)=THETA*EXP(2675.*Q(I,J,K)/FACT)

30     CONTINUE
40     CONTINUE

        RETURN
        END

! ---------- preciptatable water in mm -----------------
!
! input:
!      q: water vapor mixing ratio (kg/kg)
!      p: pressure (Pa)
!
! output:
!    ppw: precipitable water (mm)
!     
	subroutine precip_water(q,p,nz,ppw)
        real q(nz), p(nz)
        
        denw=1000
        g=9.8

	ppw=0
	do k=nz,2,-1
         ppw=ppw+0.5*(q(k)+q(k-1))*(p(k-1)-p(k))
	enddo

        ppw=1000*ppw/(denw*g)   ! convert to mm

	return
	end

! ------------------------------

	subroutine fill2d(filler,q,nx,ny,nz,klev)

	real filler(nx,ny), q(nx,ny,nz)

	do j=1,ny
	 do i=1,nx
          filler(i,j)=q(i,j,klev)
         enddo
        enddo

        return
	end

! ---------------------------
! subroutine to calculate meteorological wind direction in deg
! 0 deg = N
!
! input: u - zonal wind
!        v - meridional wind
!
! output: dir - meteorological wind direction (deg)

       subroutine winddir(u,v,dir,nx,ny,nz)

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

! ----- function to calculate dewpoint --------
! input: p - pressure (Pa)
!        rs - saturation water vapor mixing ratio (kg/kg)
!
! output: td - dewpoint (K)
       
      function td(p,rs)

      implicit none
      real rr,rs,es,esln,p,td

      rr=rs+1e-8
      es=p*rr/(.622+rr)
      esln=log(es)
      td=(35.86*esln-4947.2325)/(esln-23.6837)
      return
      end

! ------ function to calculate saturation water vapor mixing ratio -------
!
! input: p - pressure (Pa)
!        T - temperature (K)
!
! output: rs: saturation water vapor mixing ratio (kg/kg)
!
      FUNCTION RS(P,T)
      ES=610.78*EXP(17.269*(T-273.16)/(T-35.86))
      RS=.622*ES/(P-ES)
      RETURN
      END

! ---------- interface for interpolation -------------- !

	subroutine interp_inf(a,nx,ny,nz,bx,by,a_1d)

	real a(nx,ny,nz), a_1d(nz)
        real bx, by

	real dum2d(nx,ny)

         do k=1,nz
          do j=1,ny
           do i=1,nx
            dum2d(i,j)=a(i,j,k)
           enddo
          enddo
          call intr(dum2d,nx,ny,bx,by,a_1d(k),.false.)  ! bi-quadratic doesn't work too well
         enddo

         return
         end


c ---------------------------------------------------------------------------------
c note: bi-quadratic inerpolation sometimes does not work well
c
         SUBROUTINE INTR(P,NX,NY,BX,BY,BB,QUAD)                            INTRR.2
C                                                                          INTRR.3
C      PERFORMS BI-QUADRATIC INTERPOLATION WHERE POSSIBLE,                 INTRR.4
C      LINEAR INTERPOLATION IN OUTSIDE GRID INTERVAL, AND                  INTRR.5
C      LINEAR EXTRAPOLATION OUTSIDE GRID.                                  INTRR.6
C                                                                          INTRR.7
C            JULY, 1968            GLAHN, HOLLENBAUGH                      INTRR.8
C            AUGUST, 1977            ADAPTED BY W. LOTTES                  INTRR.9
C                                                                          INTRR.10
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   INTRR.11
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   INTRR.12
C                                                                          INTRR.13
C      ARGUMENTS                                                           INTRR.14
C            P    =  GRID FIELD TO INTERPOLATE FROM                        INTRR.15
C            NX   =  X DIMENSION OF GRID (JMAX IN RAWINS)                  INTRR.16
C            NY   =  Y DIMENSION OF GRID (IMAX IN RAWINS)                  INTRR.17
C            BX   =  X-COORDINATE, FROM LEFT                               INTRR.18
C            BY   =  Y-COORDINATE, FROM BOTTOM                             INTRR.19
C            BB   =  INTERPOLATED (OR EXTRAPOLATED) VALUE RETURNED         INTRR.20
C                  TO CALLING PROGRAM                                      INTRR.21
C            QUAD =  LOGICAL FLAG TO INDICATE LINEAR OR                    INTRR.22
C                  BI-QUADRATIC INTERPOLATION                              INTRR.23
C                                                                          INTRR.24
C                                                                          INTRR.25
      LOGICAL QUAD                                                         INTRR.26
      REAL P(NX,NY),B(4)                                                   INTRR.27
C                                                                          INTRR.28
104   NBX=INT(BX)                                                          INTRR.29
      NBY=INT(BY)                                                          INTRR.30
      IF (NBX-1) 114,120,111                                               INTRR.31
111   IF (NBX-(NX-1)) 112,120,115                                          INTRR.32
112   IF (NBY-1) 121,130,113                                               INTRR.33
113   IF (NBY-(NY-1)) 140,130,123                                          INTRR.34
114   NBX=1                                                                INTRR.35
      GO TO 120                                                            INTRR.36
115   NBX=NX-1                                                             INTRR.37
120   IF (NBY-1) 121,130,122                                               INTRR.38
121   NBY=1                                                                INTRR.39
      GO TO 130                                                            INTRR.40
122   IF (NBY-NY) 130,123,123                                              INTRR.41
123   NBY=(NY-1)                                                           INTRR.42
C                                                                          INTRR.43
C---STATEMENT 130 STARTS BI-LINEAR INTERPOLATION-EXTRAPOLATION.            INTRR.44
C                                                                          INTRR.45
130   NBXP1=NBX+1                                                          INTRR.46
      NBYP1=NBY+1                                                          INTRR.47
      DX=BX-FLOAT(NBX)                                                     INTRR.48
      DY=BY-FLOAT(NBY)                                                     INTRR.49
      BB=P(NBX,NBY)+(P(NBXP1,NBY)-P(NBX,NBY))*DX+(P(NBX,NBYP1)-            INTRR.50
     1 P(NBX,NBY))*DY+(P(NBX,NBY)+P(NBXP1,NBYP1)-P(NBX,NBYP1)-             INTRR.51
     2 P(NBXP1,NBY))*DX*DY                                                 INTRR.52
      RETURN                                                               INTRR.53
C                                                                          INTRR.54
C---STATEMENT 140 STARTS BI-QUADRATIC INTERPOLATION.                       INTRR.55
C                                                                          INTRR.56
140   IF(.NOT.QUAD) GO TO 130            ! DO PURE LINEAR ?                INTRR.57
      DX=BX-FLOAT(NBX)                                                     INTRR.58
      DY=BY-FLOAT(NBY)                                                     INTRR.59
      NBYP2=NBY+2                                                          INTRR.60
      NBYP1=NBY+1                                                          INTRR.61
      NBYM1=NBY-1                                                          INTRR.62
      FCT=(DY**2-DY)/4.                                                    INTRR.63
      FET=(DX**2-DX)/4.                                                    INTRR.64
      DO 145 J=1,4                                                         INTRR.65
      N=NBX-2+J                                                            INTRR.66
      B(J)=P(N,NBY)+(P(N,NBYP1)-P(N,NBY))*DY+(P(N,NBYM1)+P(N,NBYP2)-       INTRR.67
     1 P(N,NBY)-P(N,NBYP1))*FCT                                            INTRR.68
145      CONTINUE                                                          INTRR.69
      BB=B(2)+(B(3)-B(2))*DX+(B(1)+B(4)-B(2)-B(3))+FET                     INTRR.70
      RETURN                                                               INTRR.71
      END

! ===============
      subroutine calc_cltop_tempc(n1,n2,n3,c,        
     *                     b,a,undef)

       ! c is the temperature in Celcius
       ! b is the total condensate mixing ratio in kg/kg
       ! a is the cloud top temperature

      integer :: n1, n2, n3

      real :: undef

      real, dimension(n1,n2,n3) :: c, b
      real, dimension(n1,n2) :: a

      a = undef

       do 100 j=1,n2
        do 70 i=1,n1
         do 50 k=n3-1,2,-1
          if (b(i,j,k).ge.0.1e-3) then
           a(i,j)=c(i,j,k)
           go to 70
          endif
50       continue
!           a(i,j,1)=0.5*(c(i,j,1)+c(i,j,2))
70      continue
100     continue

        return
        end

! ===============
	subroutine fill_array(a_in,nx,ny,nz,a_out,nz_out)

        integer nx,ny,nz
        integer nz_out

        real a_in(nx,ny,nz)
        real a_out(nx,ny,nz_out)

        do k=1,nz_out
         do j=1,ny
          do i=1,nx
           a_out(i,j,k)=a_in(i,j,k)
          enddo
         enddo
        enddo

        return
        end
