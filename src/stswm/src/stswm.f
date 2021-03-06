      PROGRAM STSWM
C
C THIS IS THE MAIN PROGRAM FOR A GLOBAL SPECTRAL SHALLOW WATER MODEL.
C THE NUMERICAL METHODS FOR THIS CODE ARE DOCUMENTED IN THE NCAR
C TECHNICAL NOTE TN-343+STR "DESCRIPTION OF A GLOBAL SHALLOW WATER
C MODEL BASED ON THE SPECTRAL TRANSFORM METHOD" BY              
C JAMES J. HACK AND RUEDIGER JAKOB.
C EVERY EFFORT WAS MADE TO USE STANDARD FORTRAN 77 TO ACHIEVE
C PORTABILITY. THE ONLY EXTENSIONS USED ARE INCLUDE FILES
C AND NAMELIST-SPECIFICATION OF MODEL PARAMETERS.
C SUBROUTINES FROM THE FOLLOWING LIBRARIES HAVE BEEN USED:
C NCAR GRAPHICS PACKAGE 3.01
C NAG NUMERICAL INTEGRATION ROUTINE D01AHE
C ECMWF FAST FOURIER TRANSFORM FFT991
C MACHINE INDEPENDENT FILE FORMAT NETCDF
C
C THE MAIN PROGRAM CALLS INPUT TO READ NAMELIST MODEL PARAMETERS
C AND SET UP ARRAYS FOR SPECTRAL TRANSFORMATIONS, AND THEN CALLS
C INIT TO SET UP THE TEST CASE PARAMETERS. ROUTINES ERRANL AND
C NRGTCS ARE CALLED ONCE BEFORE THE MAIN TIMESTEPPING LOOP FOR
C ERROR NORMALIZATION, ONCE AFTER THE MAIN TIMESTEPPING FOR 
C PLOTTING THE TIME DEPENDENCE OF ENERGETICS DATA AND ERRORS,
C AND PERIODICALLY DURING THE TIMESTEPPING.
C THE PROGNOSTIC FIELDS ARE INITIALIZED USING ROUTINE ANLYTC,
C WHICH PROVIDES THE ANALYTIC SOLUTION. EACH CALL TO STEP 
C ADVANCES THE COMPUTED FIELDS BY A TIMESTEP DT.
C COMPUTED FIELDS ARE PLOTTED BY ROUTINE PLOTS, OR 
C WRITTEN TO A NETCDF FILE USING ROUTINE OUTPTP.
C
C CALLED BY:
C CALLS: ANLYTC, ERRANL, GACWK, GCLKS, GCLWK, GDAWK, GLAT, 
C        GOPKS, GOPWK, INIT, INPUT, NRGTCS, 
C        OUTPTP, PLOTS, SHTRNS, SPCANL, STEP, ZD
C
C REVISIONS:
C 7-10-92 CHANGE TO CCM CODING CONVENTIONS (R. JAKOB)
C
C---- Model Parameters -------------------------------------------------
C
      INCLUDE 'params.i'
C
C---- Common Blocks ----------------------------------------------------
C
C CONSTANTS & TIMESTEPS
      INCLUDE 'consts.i'
C TRANSFORM ARRAYS
      INCLUDE 'trnsfm.i'
C TIME DEPENDENT FIELDS
      INCLUDE 'tdvars.i'
C WORKSPACE
      INCLUDE 'wrkspc.i'
C MOUNTAIN DATA
      INCLUDE 'finit.i'
C
C------ Local Parameters -----------------------------------------------
C
C     FLAG TO ENABLE ORTHOGONALITY/NORMALITY CHECK OF
C     BASISFUNCTIONS
C
      LOGICAL LORTHO
      PARAMETER (LORTHO = .FALSE.)
C
C     FLAG FOR POSTPROCESSING OF SPECTRAL COEFFICIENTS IN NETCDF FILE
C     (INSTEAD OF COMPUTING FIELDS, THEY ARE READ IN FROM FILE FNIN)
C
      LOGICAL LDIFF
      PARAMETER (LDIFF = .FALSE.)
C
C------ Local Variables ------------------------------------------------
C
C     OUTPUT AND ANALYSIS FLAGS
C
      LOGICAL OUT, GPHS, ENERGY, CALLEA 
C    
C     ENERGETICS, ERROR ANALYSIS AND OUTPUT INDEX
C
      INTEGER ECTR,L2CTR,OCTR
C
C     LATITUDE
C
      REAL RLAT
C
C     GRIDPOINT WEIGHT
C
      REAL GWTS
C
C     MEAN GEOPOTENTIAL AND MAX VELOCITY FOR PARAMETER CHECK
C
      REAL MAXV,COUR,VLEN
C
C     TEMPORARIES
C
      INTEGER I,J,ISZ
C     
C     YGYU add 
      INTEGER ip, k2
      CHARACTER(len=20) :: filename
C     scalar not array
      REAL    DIV_s,ZETA_s,PHI_s,U_s,V_s,TOPO_s
      INTEGER nx, ny
      REAL    dx, dy, ymin, ymax, h_A2, u_A2, v_A2
      REAL    RADIUS, DIST, H_mount, MOUNTA
      REAL PI
      PARAMETER (PI=3.141592653589793)

C
C----- External Functions ----------------------------------------------
C
C GRID LATITUDE AND WEIGHTS
C
      EXTERNAL GLAT, WEIGHT
      REAL GLAT, WEIGHT
C
C----- Statement Function ----------------------------------------------
C
C     Time Interval Check
C
      REAL XTAU
      LOGICAL EVENT
      EVENT(XTAU) = MOD(NSTEP,MAX0(1,
     $              IFIX(XTAU * 3600.0/DT + 1000.0*EPS))) .EQ. 0
C                                                                               
C----- Executable Statements -------------------------------------------        
C
C     COPYRIGHT MESSAGE
C
      WRITE(6,10)
   10 FORMAT(/,' SPECTRAL TRANSFORM SHALLOW WATER MODEL, Version 2.0',
     $      /,' Copyright (C) 1992',
     $      /,' University Corporation for Atmospheric Research',
     $      /,' All Rights Reserved',/)
C                                                                               
C     TIMESTEP COUNTER
C
      NSTEP = 0                                                                 
C
C     MODEL TIME
C
      TAU   = 0.0                                                               
C
C     INDEX INTO CIRCULAR BUFFER FOR TIME DEPENDENT VARIABLES
C
      LN = 1
C
C     MEAN GEOPOTENTIAL AND MAX. VELOCITY FOR PARAMETER CHECK
C
      PHIBAR = 0.0
      MAXV = 0.0
C
C     INPUT ROUTINE (SET UP CONSTANTS, ETC.)                               
C     READ MODEL PARAMETERS AS NAMELIST,
C     SET UP ARRAYS FOR SPECTRAL TRANSFORM PROCEDURE 
C                                                                               
      CALL INPUT                                                                
C
C     OPEN NCAR GRAPHICS
C
      IF (LGPHS) THEN
C        CALL OPNGKS
         CALL GOPKS(0,ISZ)
         CALL GOPWK(1,2,1)
         CALL GACWK(1)
      ENDIF
C
C     ORTHOGONALITY CHECK
C
      IF (LORTHO) THEN
         CALL ORTHO
      ENDIF
C
C     INITIALIZATION ROUTINE FOR TEST CASES 
C                                                                               
      CALL INIT
C
C     INITIALIZE FIELD VARIABLES AT TIME = 0.0
C
      IF (LDIFF) THEN 
         CALL ANLYTC(LVLS,LN,0.0,7,WS1,WS2,WS3,DIV,ZETA)
      ELSE 
         CALL ANLYTC(LVLS,LN,0.0,ICOND,WS1,WS2,WS3,DIV,ZETA)
      ENDIF
C
C     SPECTRAL TRANSFORM ALGORITHM USES MODIFIED VARIABLES:
C     PHI  := G*(H-MOUNT)-PHIBAR    (USE FLUID DEPTH)
C     UCOS := U*COS(RLAT)           (REDEFINITION BECAUSE OF
C     VCOS := V*COS(RLAT)            MULTIVALUED U,V AT POLE)
C
      DO 1810 J=1,NLAT
         RLAT = GLAT(J)
C
C        LATITUDE = RLAT = GLAT(J)
C
C        WEIGHT INDEPENDENT OF LONGITUDE
C
         GWTS = WEIGHT(1,J)
C
         DO 1805 I=1,NLON
            UCOS(I,J)   = WS1(I,J)*COS(RLAT)
            VCOS(I,J)   = WS2(I,J)*COS(RLAT)
            VLEN  = SQRT(WS1(I,J)**2+WS2(I,J)**2)
            IF (VLEN .GT. MAXV) THEN
               MAXV = VLEN
            ENDIF
            PHI(I,J,LN) = GRAV*WS3(I,J)
            IF (FTOPO) THEN
               PHI(I,J,LN) = PHI(I,J,LN) - GRAV*MOUNT(I,J)
            ENDIF
            PHIBAR = PHIBAR + PHI(I,J,LN)*GWTS
C           PHIBAR = PHIBAR + PHI(I,J,LN)/REAL(NLON*NLAT)
 1805    CONTINUE
 1810 CONTINUE
C
C     SUBTRACT GLOBAL MEAN GEOPOTENTIAL (FOR SEMI-IMPLICIT TIMESTEPPING)
C     PHI' := PHI - PHIBAR
C
      DO 1830 J = 1,NLAT
         DO 1820 I = 1,NLON
            PHI(I,J,LN) = PHI(I,J,LN) - PHIBAR
C
C           SEPARATE COPY FOR INITIAL SPECTRAL COEFFICIENTS
C
            WS5(I,J)    = PHI(I,J,LN)
 1820    CONTINUE
 1830 CONTINUE
C
C     COMPUTE COURANT NUMBER
C
C     FOR EXPLICIT TIMESTEPPING, FASTEST WAVES ARE DETERMINED
C     BY GRAVITY WAVES OF SPEED SQRT(PHIBAR), FOR SEMI-
C     IMPLICIT TIMESTEPPING ROSSBY WAVES OF SPEED MAXV
C
C     THIS FORMULA IS ONLY CORRECT FOR SPECTRAL TRANSFORM CODES
C     IT SHOULD BE LESS THAN 1 AT ALL TIMES
C
      IF (SITS) THEN
         COUR = MAXV*DT*KK/A
      ELSE
         COUR = SQRT(PHIBAR)*DT*KK/A
      ENDIF
      WRITE (6,129) MAXV,COUR
  129 FORMAT (' MAX. WIND      = ',1PE16.9,/,
     $        ' COURANT NUMBER = ',0PF6.4)
C
C     WARNING FOR TOO LONG TIMESTEP
C
      IF (COUR .GE. 1.0) THEN
         WRITE(0,130) COUR
  130    FORMAT(/,' STSWM: WARNING FROM MAIN PROGRAM:',
     $          /,' TIMESTEP TOO LONG FOR EXPERIMENT',
     $          /,' COURANT NUMBER = ',0PF8.4)
      ENDIF
C
C     WRITE MEAN GEOPOTENTIAL
C
      WRITE(6,135) PHIBAR
  135 FORMAT(/,' GLOBAL MEAN STEADY GEOPOTENTIAL = ',1PE16.9)
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     PLOT INITIAL FIELDS
C    
      CALL PLOTS(0,WS3,WS4,WS1,WS2,ZETA,DIV,LVLS,LN)
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''        
C     GET INITIAL CONSERVATION ESTIMATE
C     
      ECTR  = 1
      IF (EGYFRQ .LE. TAUE) THEN
         CALL NRGTCS(WS3,WS1,WS2,ZETA,DIV,MOUNT,LN,ECTR,.FALSE.) 
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C'    GET INITIAL ERROR ESTIMATE                                       '
C'                                                                     '
      L2CTR = 1
      IF (ERRFRQ .LE. TAUE) THEN
         CALL ERRANL(WS3,WS1,WS2,DIV,ZETA,LN,.FALSE.,L2CTR,.FALSE.)
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C
C     COMPUTE INITIAL SPECTRAL COEFFICIENTS
C     DESTROY WS1,WS2,WS5 !!!
C
      CALL ZD(WS1,WS2,DIVSC,ZETASC)
      CALL SHTRNS(1,1,-1,WS5,PHISC)
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     SAVE INITIAL SPECTRAL COEFFICIENTS AS REFERENCE SOLUTION
C
      IF (TAUO .LE. TAUE) THEN
         OCTR = 1
         CALL OUTPTP (ZETASC,DIVSC,PHISC,PHIBAR,LDIAG,.TRUE.,OCTR)
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     INITIAL SPECTRAL ANALYSIS
C
      IF (SPCFRQ .LE. TAUE) THEN
         CALL SPCANL(PHISC,ZETASC,DIVSC)
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C                                                                               
C     MAIN COMPUTATIONAL CONTROL                                                
C                                                                               
  500 NSTEP  = NSTEP + 1                                                        
      TAU    = NSTEP*DT/3600.0                                             
      IF (LDIFF) THEN
C
C        READ DATA FROM FILE
C
         CALL ANLYTC(LVLS,LN,NSTEP*DT,7,WS1,WS2,WS3,DIV,ZETA)
C
      ELSE
C
C        COMPUTE NEXT TIMESTEP
C
         CALL STEP  
C
C        TRANSFORM BACK TO STANDARD VARIABLES
C
         DO 100 J=1,NLAT
            RLAT = GLAT(J)
C
C           LATITUDE = RLAT = GLAT(J)
C
            DO 90 I=1,NLON
               WS1(I,J) = UCOS(I,J)/COS(RLAT)
               WS2(I,J) = VCOS(I,J)/COS(RLAT)
               WS3(I,J) = (PHI(I,J,LN) + PHIBAR)/GRAV 
               IF (FTOPO) THEN
                  WS3(I,J) = WS3(I,J) + MOUNT(I,J)
               ENDIF
   90       CONTINUE
  100    CONTINUE
      ENDIF
C
C     PLOT FIELD VALUES ?
C
      GPHS = EVENT(GPHFRQ)
C
C''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     STANDARD ENERGETICS INFORMATION                                           
C                                                                               
      ENERGY = EVENT(EGYFRQ)                                                    
      IF (ENERGY) THEN 
         ECTR = ECTR + 1                                                        
         CALL NRGTCS(WS3,WS1,WS2,ZETA,DIV,MOUNT,LN,ECTR,.FALSE.) 
      ENDIF
C
C''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     SPECTRAL ANALYSIS
C
      IF (EVENT(SPCFRQ)) THEN
         CALL SPCANL(PHISC,ZETASC,DIVSC)
      ENDIF                                                                     
C                                                                               
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''        
C'    ERROR ANALYSIS (COMPARE WITH ANALYTIC SOLUTION)
C'                                                                             

      CALLEA = EVENT(ERRFRQ)                                                    
      IF (CALLEA) THEN 
         L2CTR = L2CTR + 1                                                      
         CALL ERRANL(WS3,WS1,WS2,DIV,ZETA,LN,GPHS,L2CTR,.FALSE.)  
      ENDIF                                                                     
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     STANDARD GRAPHICAL OUTPUT
C
      IF (GPHS) THEN
         CALL PLOTS(0,WS3,WS4,WS1,WS2,ZETA,DIV,LVLS,LN)
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''        
C     WRITE SPECTRAL COEFFICIENTS FOR REFERENCE SOLUTION
C
      OUT = EVENT(TAUO)                                                         
      IF (OUT) THEN
         OCTR = OCTR + 1
         CALL OUTPTP (ZETASC,DIVSC,PHISC,PHIBAR,LDIAG,.FALSE.,OCTR)
      ENDIF
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     SIMULATION TIME OVER ?
C
      IF (TAU + EPS .GE. TAUE) GOTO 1000
      GO TO 500                                           
C
C     POSTPROCESSING CODE
C
 1000 CONTINUE
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
C     PLOT ENERGETICS OVER TIME  
C     (ONLY IF RELEVANT DATA)
C                                                                               
      IF (EGYFRQ .LE. TAUE) THEN
         CALL NRGTCS(WS3,WS1,WS2,ZETA,DIV,MOUNT,LN,ECTR,.TRUE.) 
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''        
C     PLOT ERRORS OVER TIME
C     (ONLY IF RELEVANT DATA)
C                                                                      '        
      IF (ERRFRQ .LE. TAUE) THEN
         CALL ERRANL(WS3,WS1,WS2,DIV,ZETA,LN,.FALSE.,L2CTR,.TRUE.)
      ENDIF
C
C'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''        
C                                                                               
C     CLOSE NCAR GRAPHICS
C
      IF (LGPHS) THEN
C        CALL CLSGKS
         CALL GDAWK(1)
         CALL GCLWK(1)
         CALL GCLKS
      ENDIF

C
C     NORMAL TERMINATION SHALLOW WATER MODEL
C
C-----------------------------------------------------------------------
C     s1.1  ouput  h/U/V on icos lon-lat grid
C
C      glvl=5
C      nip = 10 * 4**glvl + 2
C      niE = 30 * 4**glvl

      write(filename, '(a,i2.2,a)') 'sph_var_G',glvl,'.dat' 
      open (11, file=trim(filename), status='unknown', 
     $      form='unformatted')
      read (11) icos_grid(1:2,1:nip)
      do ip=1,nip
         RLON= icos_grid(1,ip)
         RLAT= icos_grid(2,ip)
C                  MMH,NNH,KKH --> MM,NN,KK  H= high resolution
         CALL EVAL(DIVSC,ZETASC,PHISC,MM,NN,KK,A,ALPHA,
     $        OMEGA,RLON,RLAT,DIV_s,ZETA_s,PHI_s,U_s,V_s)
         h_A(ip)= (PHI_s + PHIBAR)/GRAV
         u_A(ip)= U_s
         v_A(ip)= V_s
C
C-- YGYU: donot need H_s  surface
C         IF (FTOPO) THEN
C            h_A(ip)= h_A(ip) + TOPO_s
C         ENDIF
      enddo
      !
      !- test
      ! 
      do ip=1,nip, 4**glvl/4
         write(6,'(2x,a,i7,5f10.2)') 
     $           'ip,icos_grid(1:2,ip),h_A(ip),u_A(ip),v_A(ip)',
     $            ip,icos_grid(1:2,ip),h_A(ip),u_A(ip),v_A(ip)
      enddo

      close(11)
      open(unit=11,file='stswm_icos_huv.dat',status='unknown', 
     $     form="unformatted", action='write')
      write(11) icos_grid(1:2,1:nip), h_A(1:nip),
     $     u_A(1:nip), v_A(1:nip)
      close(11)
      write(6,*) 'Finish up writing stswm_icos_huv.dat'


C-----------------------------------------------------------------------
C     s1.2  ouput  h/U/V on regular lon-lat grid
C
C      glvl=5
C      nip = 10 * 4**glvl + 2
C      niE = 30 * 4**glvl

      MOUNTA = 2000.0
      RADIUS = PI/9.0
      open(unit=11,file='stswm_reg_latlon.dat',status='unknown', 
     $     form="unformatted", action='write')
C      
      write(6, *) ' t1'
      nx= 10*2**(glvl-1)
      ny= nx
      dx=2.d0*pi/dble(nx-1)
      ymin=-pi/2.d0*0.9999d0    ! avoid pole problem   vs. 0.995
      ymax=-ymin
      dy=(ymax-ymin)/dble(ny-1)
      k2=0
C     total points      
      write(11) nx, ny   
      do j=1,ny
         do i=1,nx
            k2=k2+1
C            g_ll(k2,1)= ymin + dy*dble(j-1) ! lat
C            g_ll(k2,2)= dx*dble(i-1) ! lon
C             g_ll(k2,2)= -pi/2.d0 + dy*dble(j-1)
C             write(6,141) 'k2, g_ll(k2,1),g_ll(k2,2)=', k, g_ll(k2,1),g_ll(k2,2)
            RLON= dx*dble(i-1)  ! lon
            RLAT= ymin + dy*dble(j-1) ! lat
C                  MMH,NNH,KKH --> MM,NN,KK  H= high resolution
            CALL EVAL(DIVSC,ZETASC,PHISC,MM,NN,KK,A,ALPHA,
     $           OMEGA,RLON,RLAT,DIV_s,ZETA_s,PHI_s,U_s,V_s)
            h_A2= (PHI_s + PHIBAR)/GRAV
            u_A2= U_s
            v_A2= V_s
            IF (FTOPO) THEN
C
C            h is  h + hb_a : total height
               DIST = SQRT((RLON -  1.5*PI)**2 + (RLAT - PI/6.0)**2)
               IF (DIST .LT. RADIUS) THEN
                  H_mount = MOUNTA*(1.0 - DIST/RADIUS)
               ELSE
                  H_mount = 0.0
               ENDIF
               h_A2= h_A2 + H_mount
            ENDIF
C     
C           reuse  icos_grid  as  lat_lon grid
            write(11) RLON, RLAT, h_A2, u_A2, v_A2
            if ( mod(k2, nx*ny/10) .EQ. 1 ) then
               write(6,'(2x,a,i7,8f10.2)') 
     $              'k2,lon,lat,h,u,v', k2, RLON, RLAT, h_A2, u_A2, v_A2
            endif
         enddo
      enddo
      close(11)
      write(6,*) 'Finish up writing stswm_reg_latlon.dat'

      STOP 'stswm finishes'
C-----------------------------------------------------------------------
      END PROGRAM
