c ------- for calculating CAPE, CIN, etc ------------------------------------------
c         subroutine written by Greg Thompson
c
c         ARRAYS GO FROM BOTTOM-UP
c
c input:
c       zpts: number of points in z direction
c       xpts: number of points in x direction
c       ypts: number of points in y direction
c         zt: height of the grid points (m)
c      theta: potential temperature (K)
c      totpi: Exner function (J/(K kg))
c      gmixr: water vapor mixing ratio (kg/kg)
c      TYPE: 1: CAPE (J/kg)  
c            2: CIN (J/kg)
c            3: height of LCL (m)
c            4: pressure of LCL (hPa)
c            5: height of LFC (m)
c            6: pressure of LFC (hPa)
c      
c      udef: undefined value
c
c output:
c      svar (depends on TYPE)

c      subroutine stabcalc(zpts,xpts,ypts,zt,theta,totpi,
c     +       gmixr,svar,TYPE,udef)

      subroutine stabcalc(zpts,xpts,ypts,zt,theta,totpi,
     +       gmixr,cape_out,cin_out,udef)

      integer zpts,xpts,ypts,i,j,k,TYPE
      real theta(xpts,ypts,zpts),totpi(xpts,ypts,zpts),
     +     rtgt(xpts,ypts),gmixr(xpts,ypts,zpts),svar(xpts,ypts),
     +     gtmpc(zpts),gpres(zpts),ptemp(zpts),grmix(zpts),
     +     gdwpt(zpts),ghgt(zpts),vapprs,zt(xpts,ypts,zpts)

      ! --- added by WC: 2003-03-20 ------- !
      real cape_out(xpts,ypts) ,cin_out(xpts,ypts)

      parameter (cp=1004.)

      do i=1,xpts
       do j=1,ypts
        do k=1,zpts
         gpres(k)=1000.*(totpi(i,j,k)/1004.0)**(1004.0/287.0)
         gtmpc(k)=theta(i,j,k)*totpi(i,j,k)/cp-273.16
         if(gmixr(i,j,k).gt.0.)then
          vapprs=gpres(k)*gmixr(i,j,k)/(0.622+gmixr(i,j,k))
          gdwpt(k)=(240.97*alog(vapprs/6.1365))/
     +               (17.502-alog(vapprs/6.1365))
         else
           gdwpt(k)=-999.0
         endif
         ghgt(k)=zt(i,j,k)
         ptemp(k)=theta(i,j,k)
         grmix(k)=gmixr(i,j,k)
         enddo
         call sndganls(gpres,gtmpc,gdwpt,ghgt,zpts,ptemp,grmix,
     +                 gcape,gcin,hlcl,plcl,hlfc,plfc,udef)

c        IF ( TYPE .EQ. 1 ) THEN
c         svar(i,j)=gcape
c        ELSEIF ( TYPE .EQ. 2 ) THEN
c         svar(i,j)=gcin
c        ELSEIF (type.eq.3) then
c         svar(i,j)=hlcl
c        elseif (type.eq.4) then
c         svar(i,j)=plcl
c        elseif (type.eq.5) then
c         svar(i,j)=hlfc
c        elseif (type.eq.6) then
c         svar(i,j)=plfc
c        ENDIF

        cape_out(i,j)=gcape
        cin_out(i,j)=gcin

       enddo
      enddo
      return
      end
c ------------------------------------------------------
      SUBROUTINE SNDGANLS(PRES,TEMP,DEWP,HT,NLEVELS,THETA1,RMIX,
     +                    CAPE,CIN,CLCLH,CLCLP,CLFCH,CLFCP,UDEF)
c
C---This subroutine computes some sounding analysis parmeters like
C---Lifted, K, Totals, and Sweat indicies; Convective temp,
C---info about LCL,CCL, LFC and EL; Richardson number, vertical
C---velocity profile, and precipitable water.  Tc,CCL,Ri and vv
C---may not be working yet.    - Greg Thompson 4/21/93
C
      PARAMETER(MAXLEV=200,XKELVIN=273.16)
      parameter(R=287.,
     +          CP=1004.,
     +          P00=100000.,
     +          G=9.80,
     +          XPI=3.1415927,
     +          PI180=3.1415927/180.,
     +          SPCON=111120.,
     +          ALVL=2.35E6,
     +          ERAD=6367000.,
     +          RCP=R/CP,
     +          P00I=1./P00,
     +          CPR=CP/R,
     +          GRAV=9.8,
     +          RGAS=287.,
     +          XEPS=0.62197,
     +          CTV=1.-XEPS)
C
      REAL PRES(NLEVELS),TEMP(NLEVELS),DEWP(NLEVELS),
     & HT(NLEVELS),
     & THETA1(NLEVELS),THETA(MAXLEV),
     & RMIX(NLEVELS),
     & TXP(MAXLEV),PTP(MAXLEV),TVP(MAXLEV),
     & PRES2(MAXLEV),HT2(MAXLEV)
C
C
C
C---Compute surface stuff, also open a file to dump all the good stuff
C
C
C---If input levels do not provide better than 50mb increments then
C---interpolate data to provide 50mb increments by calling pres interpolation
C
c     CALL PINTERP(PRES,TEMP,DEWP,WDIR,WSPD,HT,NLEVELS,NLEVS)
      nlevs=nlevels
C
C---Compute all derivable quantitites from the input variables
C---Theta, mixing ratio, vapor pres, sat vapor pres, virtual temp,
C--- u-winds, v-winds, density, wet bulb theta, wet bulb temp, prec water.
C
      TWMIN=999.
      DO 10 K=1,NLEVS
         THETA(K) = THETA1(K) - XKELVIN
 10   CONTINUE
C
      HT15 = HT(1)+1000.
      SDH = 0.
      DO 60 K=2,NLEVS
         IF(HT(K).GE.HT15) GOTO 22
 60   CONTINUE
      RETURN
 22   CONTINUE
      IHT15 = K
      SDELTP=0.
      SDELTT=0.
      SDELTD=0.
      DO 70 K=2,IHT15
         DH = HT(K) - HT(K-1)
         SDH = SDH+DH
         SDELTP = SDELTP + PRES(K-1)-PRES(K)
         SDELTT = SDELTT + (0.5*DH*(TEMP(K-1)+TEMP(K)))
         SDELTD = SDELTD + (0.5*DH*(DEWP(K-1)+DEWP(K)))
 70   CONTINUE
      SPPR        = PRES(1)-SDELTP*0.5
      SPTT        = SDELTT/SDH
      SPTD        = SDELTD/SDH
      SPRMIX      = WMRJB1(SPPR,SPTD)
C
      TLCL=TCONJB1(SPTT,SPTD)
      PLCL=PCONJB1(SPPR,SPTT,TLCL)
      RLCL=WMRJB1(PLCL,TLCL)
      WBPTLCL     = POWTJB1(SPTT,SPPR,SPTD)
      WBTLCL      = SATLFTJB1(WBPTLCL,PLCL)
      DO 80 K=1,NLEVS
         IF(PRES(K).LE.PLCL) GOTO 23
 80   CONTINUE
      RETURN
 23   CONTINUE
      ILCL        = K
      HLCL        = HT(ILCL) - (RGAS*(TLCL+XKELVIN)/GRAV)
     &                           *ALOG(PLCL/PRES(ILCL))

      CLCLH=HLCL
      CLCLP=PLCL
C
      DO 100 K=1,NLEVS
         IF(PRES(K).LE.PLCL) THEN
            TXP(K) = SATLFTJB1(WBPTLCL,PRES(K))
         ELSE
            TXP(K) = SATLFTJB1(WBPTLCL,PRES(K))
         ENDIF
         SEP = ESWJB1(TXP(K))
         WP  = WMRJB1(PRES(K),TXP(K))
         PTP(K)    = (TXP(K)+XKELVIN)*(1000./PRES(K))**RCP - XKELVIN
 100  CONTINUE
C
C
C---Find EL (equilibrium level), LFC (level of free convection) and
C---CCL (convective condensation level).
C
C---First Pressure at EL
C
      DB = 0.
      SUMB = 0.
      SUMBEL = 0.
      SUMBLFC = 0.
      SUMBC = 0.
      DO 250 K=NLEVS,1,-1
         IF( (PTP(K)-THETA(K)).GT.0.0) GOTO 251
 250  CONTINUE
 251  CONTINUE
      KEL=K
      IF (K.EQ.0) GOTO 300
      DIFFMIN=ABS(THETA(KEL+1)-PTP(KEL+1))
      KEL2=KEL+1
C      PRINT*, 'FINDING MORE EXACT EL'
      DO 260 K=NLEVS+1,NLEVS+25
         PINC2 = (PRES(KEL)-PRES(KEL+1))/25.
         PRES2(K) = PRES(KEL+1)+(K-NLEVS)*PINC2
         TXP(K) = SATLFTJB1(WBPTLCL,PRES2(K))
         PTP(K) = (TXP(K)+XKELVIN)*(1000./PRES2(K))**RCP - XKELVIN
         THETA(K) = THETA(KEL+1) + ( (THETA(KEL)-THETA(KEL+1))*
     &                 (PRES2(K)-PRES(KEL+1)) )/(PRES(KEL)-PRES(KEL+1))
         IF(K.EQ.NLEVS+1) THEN
            THBAR = (THETA(K)+THETA(KEL+1))/2.
            TBAR= (THBAR+XKELVIN)*(PRES2(K)/1000.)**RCP
            HT2(K) = HT(KEL+1)-(RGAS*TBAR/G)*LOG(PRES2(K)/PRES(KEL+1))
         ELSE
            THBAR = (THETA(K)+THETA(K-1))/2.
            TBAR= (THBAR+XKELVIN)*(PRES2(K)/1000.)**RCP
            HT2(K) = HT2(K-1)-(RGAS*TBAR/G)*LOG(PRES2(K)/PRES2(K-1))
         ENDIF
C         PRINT*, 'PRES,THETA,PTP,HT ARE ',PRES(K),
C     &             THETA(K),PTP(K),HT(K)
         IF(ABS(THETA(K)-PTP(K)).LT.DIFFMIN) THEN
            DIFFMIN=ABS(THETA(K)-PTP(K))
            KEL2=K
         ENDIF
 260  CONTINUE
      PEL        = PRES2(KEL2)
      HTEL       = HT2(KEL2)
C
C---Sum the positive buoyancy between PRES(KEL) and PEL
C
      DO 265 K=NLEVS+1,NLEVS+25
         IF(PRES2(K).GT.PEL .AND. PRES2(K).LE.PRES(KEL)) THEN
            DB = (PTP(K)-THETA(K))*(HT2(K-1)-HT2(K))/(PTP(K)+XKELVIN)
            ! ------ added by WC: 2003-03-21 ---  !
            IF (DB.GT.0.) SUMBEL = SUMBEL + DB
!            SUMBEL = SUMBEL + DB
         ENDIF
 265  CONTINUE
      DB = 0.
C
C---Now find pressure at LFC
C
      DO 270 K=KEL,1,-1
         IF( (PTP(K)-THETA(K)).LE.0.0) GOTO 271
 270  CONTINUE
 271  CONTINUE
      KLFC=K
      DIFFMIN=(PTP(KLFC+1)-THETA(KLFC+1))
      KLFC2=KLFC+1
C      PRINT*, 'FINDING MORE EXACT LFC'
      DO 280 K=NLEVS+1,NLEVS+25
         PINC2 = (PRES(KLFC)-PRES(KLFC+1))/25.
         PRES2(K) = PRES(KLFC+1)+(K-NLEVS)*PINC2
         TXP(K) = SATLFTJB1(WBPTLCL,PRES2(K))
         PTP(K) = (TXP(K)+XKELVIN)*(1000./PRES2(K))**RCP - XKELVIN
         THETA(K) = THETA(KLFC+1) + ( (THETA(KLFC)-THETA(KLFC+1))*
     &             (PRES2(K)-PRES(KLFC+1)) )/(PRES(KLFC)-PRES(KLFC+1))
         IF(K.EQ.NLEVS+1) THEN
            THBAR = (THETA(K)+THETA(KLFC+1))/2.
            TBAR= (THBAR+XKELVIN)*(PRES2(K)/1000.)**RCP
            HT2(K)=HT(KLFC+1)-(RGAS*TBAR/G)*LOG(PRES2(K)/PRES(KLFC+1))
         ELSE
            THBAR = (THETA(K)+THETA(K-1))/2.
            TBAR= (THBAR+XKELVIN)*(PRES2(K)/1000.)**RCP
            HT2(K) = HT2(K-1)-(RGAS*TBAR/G)*LOG(PRES2(K)/PRES2(K-1))
         ENDIF
         IF(ABS(THETA(K)-PTP(K)).LT.DIFFMIN) THEN
            DIFFMIN=ABS(THETA(K)-PTP(K))
            KLFC2=K
         ENDIF
 280  CONTINUE
      PLFC       = PRES2(KLFC2)
      HTLFC      = HT2(KLFC2)

      CLFCP=PLFC
      CLFCH=HTLFC
C
C---Sum the positive buoyancy between PLFC and PRES(KLFC+1)
C
      DO 284 K=NLEVS+1,NLEVS+25
         IF(PRES2(K).LE.PLFC .AND. PRES2(K).GT.PRES(KLFC+1)) THEN
            DB = (PTP(K)-THETA(K)) * (HT2(K-1)-HT2(K))/(PTP(K)+XKELVIN)
         ! ----- added by WC: 2003-03-21 ------- !
         IF (DB.GT.0.) SUMBLFC = SUMBLFC + DB
!            SUMBLFC = SUMBLFC + DB
         ENDIF
 284  CONTINUE
      DB = 0.
C
C---Sum the positive buoyancy between PRES(KLFC+1) and PRES(KEL)
C
      DO 285 K=KLFC+1,KEL
         APTP = (PTP(K)+PTP(K-1))/2.
         ATHETA = (THETA(K)+THETA(K-1))/2.
         DB = (APTP-ATHETA) * (HT(K)-HT(K-1))/(PTP(K)+XKELVIN)
         ! ----- added by WC: 2003-03-21 ------- !
         IF (DB.GT.0.) SUMBC = SUMBC + DB
!         SUMBC = SUMBC + DB
 285  CONTINUE
C
C---Now total all three of these positive areas to obtain cape
C
      CAPE       = (SUMBEL+SUMBLFC+SUMBC)*G

C
C*** COMPUTE THE NEGATIVE AREA FROM LCL TO LFC
C
      SUMN = 0.0
      IF(ILCL .GE. KLFC+1)THEN
         SUMN = 0.0
         GOTO 297
      ELSE
       DO 295 K=ILCL,KLFC+1
          APTP = (PTP(K)+PTP(K-1))/2.
          ATHETA = (THETA(K)+THETA(K-1))/2.
          DN = (APTP-ATHETA)*(HT(K)-HT(K-1))/(PTP(K)+XKELVIN)
          ! ----- added by WC: 2003-03-21 ---- !
          IF (DN.LT.0.) SUMN = SUMN + DN
!          SUMN = SUMN + DN
 295  CONTINUE
      ENDIF
 297  CONTINUE

       IF ( SUMN .GT. 0.0 ) THEN
           SUMN = 0.0
       ENDIF

          CIN = G*SUMN

      RETURN

 300  CONTINUE
      CAPE=0.
      CIN=UDEF
      CLFCH=UDEF
      CLFCP=UDEF 

      RETURN
      END

C-------------------------------------------------------------
C
      FUNCTION ESWJB1(T)
C
      PARAMETER(ESO = 6.1078)
      POL = 0.99999683      + T*(-0.90826951E-02 +
     &   T*(0.78736169E-04  + T*(-0.61117958E-06 +
     &   T*(0.43884187E-08  + T*(-0.29883885E-10 +
     &   T*(0.21874425E-12  + T*(-0.17892321E-14 +
     &   T*(0.11112018E-16  + T*(-0.30994571E-19)))))))))
      ESWJB1 = ESO / POL**8
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION WMRJB1(P, T)
C 
      PARAMETER (XEPS = 0.62197)
      X = 0.02 * (T - 12.5 + 7500. / P)
      WFW = 1. + 4.5E-06 * P + 1.4E-03 * X * X
      FWESW = WFW * ESWJB1(T)
      R = XEPS * FWESW / (P - FWESW)
      WMRJB1 = 1000. * R
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION TCONJB1(T, D)
      S = T - D
      DLT = S * (1.2185 + 1.278E-03 * T +
     &      S * (-2.19E-03 + 1.173E-05 * S - 5.2E-06 * T))
      TCONJB1 = T - DLT
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION PCONJB1(P, T, TC)
      PARAMETER (AKAPI = 3.5037)
      TK = T + 273.16
      TCK = TC + 273.16
      PCONJB1 = P * (TCK / TK) ** AKAPI
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION WOBFJB1(T)
      X = T - 20.
      IF (X .GT. 0.) GO TO 10
      POL = 1.                     + X * (-8.8416605E-03
     &     + X * ( 1.4714143E-04   + X * (-9.6719890E-07
     &     + X * (-3.2607217E-08   + X * (-3.8598073E-10)))))
      WOBFJB1 = 15.130 / POL ** 4
      RETURN
 10   CONTINUE
      POL = 1.                     + X * ( 3.6182989E-03
     &     + X * (-1.3603273E-05   + X * ( 4.9618922E-07
     &     + X * (-6.1059365E-09   + X * ( 3.9401551E-11
     &     + X * (-1.2588129E-13   + X * ( 1.6688280E-16)))))))
      WOBFJB1 = 29.930 / POL ** 4 + 0.96 * X - 14.8
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION POWTJB1(T, P, TD)
      PARAMETER (CTA = 273.16, AKAP = 0.28541)
      PT = (T + CTA) * (1000. / P) ** AKAP - CTA
      TC = TCONJB1(T, TD)
      POWTJB1 = PT - WOBFJB1(PT) + WOBFJB1(TC)
      RETURN
      END
C-------------------------------------------------------------
      FUNCTION SATLFTJB1(THW, P)
      PARAMETER (CTA = 273.16, AKAP = 0.28541)
      IF (P .NE. 1000.) GO TO 5
      SATLFTJB1 = THW
      RETURN
 5    CONTINUE
      PWRP = (P / 1000.) ** AKAP
      TONE = (THW + CTA) * PWRP - CTA
      EONE = WOBFJB1(TONE) - WOBFJB1(THW)
      RATE = 1.
      GO TO 15
 10   CONTINUE
      RATE = (TTWO - TONE) / (ETWO - EONE)
      TONE = TTWO
      EONE = ETWO
 15   CONTINUE
      TTWO = TONE - EONE * RATE
      PT = (TTWO + CTA) / PWRP - CTA
      ETWO = PT + WOBFJB1(TTWO) - WOBFJB1(PT) - THW
      DLT = ETWO * RATE
      IF (ABS(DLT) .GT. 0.1) GO TO 10
      SATLFTJB1 = TTWO - DLT
      RETURN
      END
C -----------------------------------------------------------------
      FUNCTION  ESATJB1(T)
C  ESAT(MILLIBARS),T(KELVIN)
      TC = T-273.16
      ESATJB1 = 6.1078*EXP((17.2693882*TC)/(TC+237.3))
      RETURN
      END
C -----------------------------------------------------------------
      FUNCTION  WJB1(T,P)
C  W(GRAMS WATER VAPOR/KILOGRAM DRY AIR ), P(MILLIBAR )
      IF(T.GE.999.)GO TO 10
      X = ESATJB1(T)
      WJB1 = 621.97*X/(P-X)
      RETURN
 10   WJB1=0.0
      RETURN
      END
