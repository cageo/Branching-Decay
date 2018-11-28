C     *******************************************************************
C     ***                                                             ***
C     ***       B R A N C H I N G   D E C A Y  - MAX. 3 LEVELS        ***
C     ***                                                             ***
C     ***               XIAOMIN WANG, FEBRUARY 10, 2016               ***
C     *******************************************************************
C
C     DESCRIPTION
C     ===========
C     LAPLACE-TRANSFORMED ANALYTICAL SOLUTIONS FOR 1-D TRANSPORT
C     WITH BRANCHING FIRST-ORDER DECAY REACTIONS.
C
C     REAL-TIME SOLUTIONS ARE GENERATED USING THE DE HOOG ET AL. (1982)
C     INVERSION ROUTINE.
c
c     --> This version assumes TYPE I and III INFLOW BOUNDARY CONDITION,
C         with general time-varying concentrations
C
C
C     DEFINTION OF INPUT PARAMETERS
C     =============================
C     Q      : DARCY FLUX
C     POR    : POROSITY (SATURATED WATER CONTENT)
C     NLEVELS: TOTAL NUMBER OF BRANCHING LEVELS (MAX. = 3)
C     D      : DISPERSION COEFFICIENT
C     KD     : SORPTION COEFFICIENTS
C     LAMDA1 : FIRST-ORDER DECAY COEFFICIENTS FOR DISSOLVED PHASE
C     LAMDA2 : FIRST-ORDER DECAY COEFFICIENTS FOR SORBED PHASE
C     LEVEL  : THE LEVEL ON THE BRANCHING CHAIN(STARTING FROM ZERO)
C              example:
C                         (2)  => cis-DCE  =  (3)            
C                  (1)        =              =>
C              PCE ====> TCE ====> 1,1-DCE ====> VC
C                             =              =>  
C                         (2)  =>trans-DCE =  (3)              
C
C     NUMP   : TOTAL NUMBER OF PARENTS
C     ID     : SPECIES NUMBER OF PARENTS
C     Y      : YIELD COEFFICIENTS
C     ETA    : FRACTION OF PARENT THAT TRANSORMS INTO DAUGHTER SPEICES
C     PATH   : DECAY PATH OF SPECIES, PATH1,PATH2,PATH3 FOR EACH LEVEL
C     DELTA  : PARAMETER IN THE GENERAL BOUNDARY CONDITION
C
c     IBC    : INFLOW BOUNDARY CONDITION  
C              = 1 : TYPE I   (DIRICHLET, DELTA = 0)
C              = 3 : TYPE III (CAUCHY, DELTA = 1)
C
C     NP()    : NUMBER OF POINTS DESCRIBING INFLOW CONCENTRATION HISTORY (NSPECIES)
C     TI()    : INFLOW CONCENTRATION HISTORY POINT (TIME) (NSPECIES)
C     CI()    : INFLOW CONCENTRATION HISTORY POINT (CONCENTRATION) (NSPECIES)
C     MAXPT   : MAXIMUM NUMBER OF POINT DESCRIBING INFLOW CONCENTRATION HISTORY
C     MX_S    : MAXIMUM NUMBER OF SPECIES DURING THE DECAY PROCESS
C
C     DECLARATION OF VARIABLES
C     ========================
      IMPLICIT NONE
C 
      INTEGER NSPECIES,ISPEC,N,NLEVELS
      INTEGER NT,I
      INTEGER NX,J
      INTEGER K,M
C
      DOUBLE PRECISION Q,POR,RHOB
      DOUBLE PRECISION XMIN,XMAX,DX,XX
      DOUBLE PRECISION TT
C
      DOUBLE PRECISION DELTA
      INTEGER IBC
c
C     VARIABLE DIMENSIONS
C     -------------------
      INTEGER MX_T,MX_X,MX_S
      PARAMETER(MX_T=20,MX_X=2000,MX_S=10)
      INTEGER MAXPT
      PARAMETER(MAXPT=2000)
C
      DOUBLE PRECISION T(MX_T)
      DOUBLE PRECISION X(MX_X)
C
      DOUBLE PRECISION TI(MX_S,MAXPT),CI(MX_S,MAXPT),TS(MX_S,MAXPT),
     +                 DELC(MX_S,MAXPT)
      INTEGER NP(MX_S)    
C
      DOUBLE PRECISION D(MX_S)
      DOUBLE PRECISION LAMDA1(MX_S),LAMDA2(MX_S)
      DOUBLE PRECISION KD(MX_S)
      INTEGER          LEVEL(MX_S),NUMP(MX_S),ID(MX_S,MX_S)
      DOUBLE PRECISION Y(MX_S,MX_S),ETA(MX_S,MX_S)
      DOUBLE PRECISION FT(MX_S)
C
      CHARACTER*80 FILENAME
C
C     PARAMETERS FOR LAPLACE-TRANSFORMED SOLUTIONS
C     ============================================
      COMMON/INVDATA/Q,POR,RHOB,D,KD,LAMDA1,LAMDA2,LEVEL,NUMP,ID,Y,
     +               ETA,XX,ISPEC,NSPECIES
      COMMON/BCS/DELTA,NP,DELC,TS,TT
C
C     OPEN INPUT AND OUTPUT FILES
C     ===========================
      PRINT *,' FILENAME FOR INPUT DATA FILE : '
      READ(*,100) FILENAME
      OPEN(UNIT=55,FILE=FILENAME,STATUS='OLD')
      PRINT*,'(WARNING: PARAMETERS D, kd, LAMDA1 AND LAMDA2 CANNOT BE
     + COMPLETELY IDENTICAL' 
      PRINT*,'          BETWEEN TWO SPECIES ON THE SAME DECAY PATH)'
C
      PRINT *,' FILENAME FOR OUTPUT LISTING  : '
      READ(*,100) FILENAME
      OPEN(UNIT=66,FILE=FILENAME,STATUS='UNKNOWN')
c
      PRINT *,' FILENAME FOR PLOTTING OUTPUT : '
      READ(*,100) FILENAME
      OPEN(UNIT=67,FILE=FILENAME,STATUS='UNKNOWN')
C
C     READ PROBLEM DATA
C     =================
C     TRANSPORT PARAMETERS
C     --------------------
      READ(55,*) Q
      READ(55,*) POR
      READ(55,*) RHOB
      READ(55,*) NSPECIES
      READ(55,*) NLEVELS
      READ(55,*) (D(N),N=1,NSPECIES)
      READ(55,*) (KD(N),N=1,NSPECIES)
      READ(55,*) (LAMDA1(N),N=1,NSPECIES)
      READ(55,*) (LAMDA2(N),N=1,NSPECIES)
      READ(55,*) (LEVEL(N),N=1,NSPECIES)
      READ(55,*) (NUMP(N),N=1,NSPECIES)
      
      DO N=1,NSPECIES
        IF(NUMP(N) .EQ. 0) THEN
          READ(55,*) ID(N,1)
        ELSE
          READ(55,*) (ID(N,M),M=1,NUMP(N)) 
        END IF
      END DO
      DO N=1,NSPECIES
        IF(NUMP(N) .EQ. 0) THEN
          READ(55,*) Y(N,1)
        ELSE
          READ(55,*) (Y(N,M),M=1,NUMP(N))
        END IF
      END DO
      DO N=1,NSPECIES
        IF(NUMP(N) .EQ. 0) THEN
          READ(55,*) ETA(N,1)
        ELSE
          READ(55,*) (ETA(N,M),M=1,NUMP(N))
        END IF
      END DO
C
      READ(55,*) IBC
C
C     INFLOW BOUNDARY CONDITION DATA (NSPECIES)
C     ------------------------------------------
      DO N=1,NSPECIES
        READ(55,*) NP(N)
        do k=1,NP(N)
           READ(55,*) TI(N,K),CI(N,K)
        end do
      END DO
C
C     DISCRETIZATION PARAMETERS
C     -------------------------
      READ(55,*) XMIN,XMAX,DX
C
C     NUMBER OF TIMES
C     ---------------
      READ(55,*) NT
      READ(55,*) (T(I), I=1,NT)
C
C     ECHO INPUT DATA
C     ===============
      WRITE(66,200)
      WRITE(66,201)
      WRITE(66,210) Q
      WRITE(66,220) POR
      WRITE(66,221) RHOB
      WRITE(66,222) NLEVELS
      
      DO N=1,NSPECIES
         WRITE(66,230) N,D(N),KD(N),LAMDA1(N),LAMDA2(N),LEVEL(N),
     +                 NUMP(N)
        IF(NUMP(N) .EQ. 0) THEN
          WRITE(66,233) 0
        ELSE
          DO M=1,NUMP(N)
            WRITE(66,233) ID(N,M)
          END DO
        END IF
      END DO
      WRITE(66,231)
      DO N=2,NSPECIES
         IF(NUMP(N) .GT. 0) THEN
           DO M=1,NUMP(N)
             WRITE(66,232) ID(N,M),N,Y(N,M),ETA(N,M)  
           END DO
         END IF
      END DO
C      
C     SET INFLOW BOUNDARY CONDITION
C     =============================
      IF(IBC.EQ.3) THEN
         DELTA = 1.D0
      ELSE
         DELTA = 0.D0
      END IF
C
c     calculate inflow concentration history (NSPECIES)
c     ==================================================
      DO N=1,NSPECIES
        ts(N,1)   = ti(N,1)
        delc(N,1) = ci(N,1)
        if(np(N).gt.1) then
          do k=2,np(N)
             ts(N,k)   = (ti(N,k)+ti(N,k-1))/2.d0
             delc(N,k) = ci(N,k)-ci(N,k-1)
          END DO
        end if
      END DO
c
c
c     echo concentration history (NSPECIES)
c     --------------------------------------
      DO N=1,NSPECIES
        WRITE(66,235) IBC
        write(66,1234) N
        do k=1,np(N)
           if (k.lt.np(N)) then
              write(66,1250) ts(N,k),ts(N,k+1),ci(N,k)
           else
              write(66,1260) ts(N,k),ci(N,k)
           end if
        END DO
      END DO
c
C
C     COMPUTE X-COORDINATES
C     =====================
      NX = IDINT((XMAX-XMIN)/DX+0.5D0)+1
C
      DO J=1,NX
         X(J)=XMIN + DBLE(J-1)*DX
      END DO
C
C     LOOP OVER THE NUMBER OF TIMES
C     =============================
      DO I=1,NT
         TT=T(I)
         WRITE(66,300) TT
         WRITE(66,310)
C
C        LOOP OVER THE PROFILE
C        =====================
         DO J=1,NX
            XX=X(J)    
C
C           COMPUTE THE REAL TIME ANALYTICAL SOLUTIONS
C           ==========================================
            do n=1,nspecies
               ispec = n
               CALL HOOGD(TT,FT(n))
            end do
c
            write(66,320) xx,(ft(n),n=1,nspecies)
            write(67,320) xx,(ft(n),n=1,nspecies)
c
         END DO
      END DO

C     FORMAT STATEMENTS
C     =================
100   FORMAT(A80)
200   FORMAT(5X,'BRANCHING CHAIN',/,5X,7('=')//)
201   FORMAT(5X,'INPUT DATA',/,5X,10('-'))
210   FORMAT(5X,'DARCY FLUX                    Q  : ',1PE15.8)
220   FORMAT(5X,'POROSITY                    POR  : ',1PE15.8)
221   FORMAT(5X,'BULK DENSITY               RHOB  : ',1PE15.8)
222   FORMAT(5X,'TOTAL NUMBER OF LEVELS  NLEVELS  :  ',I1)
230   FORMAT(/5X,'SPECIES                          :  ',i1,/,
     +       5X,'DISPERSION COEFFICIENT,       D  : ',1PE15.8,/,
     +       5X,'SORPTION COEFFICIENT,         KD : ',1PE15.8,/,
     +       5X,'LAMDA1 (DISSOLVED PHASE)         : ',1PE15.8,/,
     +       5X,'LAMDA2 (SORBED PHASE)            : ',1PE15.8,/,
     +       5X,'LEVEL                            :  ',I1,/,
     +       5X,'NUMBER OF PARENTS           NUMP :  ',I1)
233   FORMAT(5X,'SPECIES NO. OF PARENTS        ID :  ',I1)
231   FORMAT(/5X,'YIELD COEFFICIENTS (Y) AND FRACTION (ETA)')
232   FORMAT(5X,'SPECIES',1X,I2,' - SPECIES',1X,I2,8X,' : ',1PE15.8,
     +       8X,1PE15.8)
235   FORMAT(/5X,'INFLOW BOUNDARY CONDITION DATA',/,
     +       5X,'------------------------------',/,
     +       5X,'INFLOW BOUNDARY CONDITION TYPE  ',I2)
300   FORMAT(/5X,'CONCENTRATION PROFILE AT TIME   : ',1PE15.8,/)
310   FORMAT(10X,'X',16X,'C(X,T)',/,5X,60('-'))
320   FORMAT(4X,1PE15.8,10(2X,1PE15.8))
1234  FORMAT(//5X,
     +       'CONSTRUCTED INFLOW CONCENTRATION HISTOGRAM (species',i3,
     +       ')',/,15X,'TIME INTERVAL',11X,'CONCENTRATION',/,5X,48('-'))
1250  FORMAT(5X,1PE13.6,' - ',1PE13.6,6X,1PE13.6)
1260  FORMAT(5X,1PE13.6,' --> INFINITY ',7X,1PE13.6)
C
C     END THE MAIN PROGRAM
C     ====================
      CLOSE(55)
      CLOSE(66)
      close(67)
C
      STOP
      END
C
C     ***************************
      COMPLEX*16 FUNCTION FBAR(P)
C     ***************************

C     DECLARATION OF VARIABLES
C     ========================
      IMPLICIT NONE
C
      INTEGER N,MX_S,ISPEC,K,NSPECIES,MX_ID,M,M1,M2,M3,TEMP1,TEMP2,
     +        N1,N2,N3 
C
      DOUBLE PRECISION Q,POR,RHOB,XX
C
      PARAMETER(MX_S=10)
      DOUBLE PRECISION D(MX_S),KD(MX_S),
     +                 LAMDA1(MX_S),LAMDA2(MX_S)
      INTEGER          LEVEL(MX_S),NUMP(MX_S),ID(MX_S,MX_S)
      INTEGER          PATH1(MX_S,MX_S),PATH2(MX_S,MX_S),
     +                 PATH3(MX_S,MX_S) 
      DOUBLE PRECISION Y(MX_S,MX_S),ETA(MX_S,MX_S)
      DOUBLE PRECISION R(MX_S),MU(MX_S)
C
      DOUBLE PRECISION DELTA
C
      INTEGER MAXPT
      PARAMETER(MAXPT=2000)
      DOUBLE PRECISION TT,TS(MX_S,MAXPT),DELC(MX_S,MAXPT)
      INTEGER NP(MX_S)
C
      COMPLEX*16 P,COMP1,COMP2,COMP3,COMP4
      COMPLEX*16 SUM(MX_S),TERM(MX_S),FS(MX_S)
      COMPLEX*16 B(MX_S),A(MX_S)
      COMPLEX*16 A2,A3,B1,B2,B3
      COMPLEX*16 KSTAR1,KSTAR2,KSTAR3,KSTAR4,
     +           KSTAR5,KSTAR6,KSTAR7,KSTAR8,
     +           KSTAR9,KSTAR10
c
      complex*16 cdexp,cdsqrt,dcmplx
C
C     PARAMETERS FOR LAPLACE-TRANSFORMED SOLUTIONS
C     ============================================
      COMMON/INVDATA/Q,POR,RHOB,D,KD,LAMDA1,LAMDA2,LEVEL,NUMP,ID,Y,
     +               ETA,XX,ISPEC,NSPECIES
      COMMON/BCS/DELTA,NP,DELC,TS,TT
C
c     compute transformed inflow concentration history (NSPECIES)
c     ============================================================
      DO N=1,NSPECIES
        sum(N) = dcmplx(0.d0,0.d0)
        do k=1,np(N)
           if(ts(N,k).lt.TT) then
              if(ts(N,k).le.0.d0) then
                 term(N) = delc(N,k)/p
              else
                 term(N) = (delc(N,k)/p)*cdexp(-p*ts(N,k))
              end if
           else
              term(N) = dcmplx(0.d0,0.d0)
           end if
           sum(N) = sum(N)+term(N)
        END DO
        fs(N) = sum(N)
      END DO
c
C      
C     CALCULATE DEGRADATION PARAMETERS
C     ================================
      DO N=1,NSPECIES
         R(N)  = 1.D0+(RHOB/POR)*KD(N)
         MU(N) = LAMDA1(N) + (RHOB/POR)*KD(N)*LAMDA2(N)
      END DO
C
C     CALCULATE SOLUTIONS FOR HOMOGENEOUS PDE
C     ========================================
      DO N=1,NSPECIES
        B(N) = (Q/(2.D0*POR*D(N))) - CDSQRT((Q/(POR*D(N)))**2
     +    +(4.D0/(POR*D(N)))*(R(N)*POR*P + POR*MU(N)))/2.d0
c
        A(N) = (Q/(2.D0*POR*D(N))) + CDSQRT((Q/(POR*D(N)))**2
     +    +(4.D0/(POR*D(N)))*(R(N)*POR*P + POR*MU(N)))/2.d0
      END DO
C
C     CALCULATE THE DECAY PATH (PARENTS, GRANDPARENTS...)
C     ===================================================
      DO N=1,NSPECIES
        IF (NUMP(N) .GT. 0) THEN
          DO M1=1,NUMP(N)
            PATH1(N,M1) = ID(N,M1)
            TEMP1 = ID(N,M1)
            IF (NUMP(TEMP1) .GT. 0) THEN
              DO M2=1,NUMP(TEMP1)
                PATH2(N,M2) = ID(TEMP1,M2)
                TEMP2 = ID(TEMP1,M2)
                IF (NUMP(TEMP2) .GT. 0) THEN
                  DO M3=1,NUMP(TEMP2)
                    PATH3(N,M3) = ID(TEMP2,M3)
                  END DO
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
C
C
C     LAPLACE DOMAIN SOLUTIONS
C     ========================
C     DECAY LEVEL 0
C     ---------
      IF(LEVEL(ISPEC) .EQ. 0) THEN
         KSTAR1 = 1.d0/(-POR*DELTA*D(ISPEC)*B(ISPEC) + Q)
     +               *(Q*fs(ISPEC))    
         FBAR  = KSTAR1*CDEXP(B(ISPEC)*XX)
      END IF
C
C     
C     DECAY LEVEL 1
C     ---------
      IF(LEVEL(ISPEC) .EQ. 1) THEN
         COMP1 = dcmplx(0.d0,0.d0)
         COMP2 = dcmplx(0.d0,0.d0)
         DO M1=1,NUMP(ISPEC)
           N1 = PATH1(ISPEC,M1)
           B1 = B(PATH1(ISPEC,M1))
           IF (B1 .eq. B(ISPEC)) then
	     write(*,*) '!WARNING: DENOMINATOR EQUALS TO 0!'
           END IF
        
           KSTAR1 = 1.d0/(-POR*DELTA*D(N1)*B1 + Q)*(Q*fs(N1))
           KSTAR2 = 1.d0/(-POR*DELTA*D(ISPEC)*B(ISPEC)+Q)*(Q*fs(ISPEC)
     +              -POR*DELTA*MU(N1)*Y(ISPEC,M1)*ETA(ISPEC,M1)*KSTAR1
     +              *B1/((B1-B(ISPEC))*(B1-A(ISPEC)))
     +              +Q*MU(N1)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +              *KSTAR1/((B1-B(ISPEC))*(B1-A(ISPEC))))
           KSTAR3 = -MU(N1)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)*KSTAR1
     +              /((B1-B(ISPEC))*(B1-A(ISPEC)))
           COMP1=COMP1+KSTAR2
           COMP2=COMP2+KSTAR3*CDEXP(B1*XX)
         END DO
         COMP1 = COMP1-(NUMP(ISPEC)-1)*Q*fs(ISPEC)/(-POR*DELTA*D(ISPEC)
     +           *B(ISPEC)+Q) 
         FBAR = COMP1*CDEXP(B(ISPEC)*XX)+COMP2
      END IF
C
C
C     DECAY LEVEL 2
C     ---------
      IF(LEVEL(ISPEC) .EQ. 2) THEN
         COMP1 = dcmplx(0.d0,0.d0)
         COMP2 = dcmplx(0.d0,0.d0)
         COMP3 = dcmplx(0.d0,0.d0)
         DO M1=1,NUMP(ISPEC)
           N2 = PATH1(ISPEC,M1)
           B2 = B(PATH1(ISPEC,M1))
           A2 = A(PATH1(ISPEC,M1))           
           DO M2=1,NUMP(N2)
             N1 = PATH2(ISPEC,M2)
             B1 = B(PATH2(ISPEC,M2))
             IF (B1 .eq. B(ISPEC) .OR. B2 .eq. B(ISPEC)) then
	       write(*,*) '!WARNING: DENOMINATOR EQUALS TO 0!'
             END IF
       
             KSTAR1 = 1.d0/(-POR*DELTA*D(N1)*B1+Q)*(Q*fs(N1))
             KSTAR2 = 1.d0/(-POR*DELTA*D(N2)*B2+Q)*(Q*fs(N2)-POR*DELTA
     +               *MU(N1)*Y(N2,M2)*ETA(N2,M2)*KSTAR1*B1/((B1-B2)
     +               *(B1-A2))+Q*MU(N1)*Y(N2,M2)*ETA(N2,M2)/D(N2)         
     +               *KSTAR1/((B1-B2)*(B1-A2))) 
             KSTAR3 = -MU(N1)*Y(N2,M2)*ETA(N2,M2)/D(N2)
     +                *KSTAR1/((B1-B2)*(B1-A2))
             KSTAR4 = 1.d0/(-POR*DELTA*D(ISPEC)*B(ISPEC)+Q)*(Q*fs(ISPEC)
     +                -POR*DELTA*MU(N2)*Y(ISPEC,M1)*ETA(ISPEC,M1)
     +                *(KSTAR2*B2/((B2-B(ISPEC))*(B2
     +                -A(ISPEC)))+KSTAR3*B1/((B1-B(ISPEC))
     +                *(B1-A(ISPEC))))+Q*MU(N2)*Y(ISPEC,M1)
     +                *ETA(ISPEC,M1)/D(ISPEC)*(KSTAR2/((B2-B(ISPEC))
     +                *(B2-A(ISPEC)))+KSTAR3/((B1-B(ISPEC))
     +                *(B1-A(ISPEC)))))
	     KSTAR5 = -MU(N2)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +                *KSTAR2/((B2-B(ISPEC))*(B2-A(ISPEC)))
             KSTAR6 = -MU(N2)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +                *KSTAR3/((B1-B(ISPEC))*(B1-A(ISPEC)))           
             COMP1=COMP1+KSTAR4
             COMP2=COMP2+KSTAR5*CDEXP(B2*XX) 
             COMP3=COMP3+KSTAR6*CDEXP(B1*XX)  
           END DO
         END DO
         COMP1=COMP1-(NUMP(ISPEC)-1)*Q*fs(ISPEC)/(-POR*DELTA*D(ISPEC)
     +           *B(ISPEC)+Q) 
         FBAR = COMP1*CDEXP(B(ISPEC)*XX)+COMP2+COMP3
      END IF
C
C
C     DECAY LEVEL 3
C     ---------
      IF(LEVEL(ISPEC) .EQ. 3) THEN
        COMP1 = dcmplx(0.d0,0.d0)
        COMP2 = dcmplx(0.d0,0.d0)
        COMP3 = dcmplx(0.d0,0.d0)
        COMP4 = dcmplx(0.d0,0.d0)
        DO M1=1,NUMP(ISPEC)
          N3 = PATH1(ISPEC,M1)
          B3 = B(PATH1(ISPEC,M1))
          A3 = A(PATH1(ISPEC,M1))        
          DO M2=1,NUMP(N3)
            N2 = PATH2(ISPEC,M2)
            B2 = B(PATH2(ISPEC,M2))
            A2 = A(PATH2(ISPEC,M2))
            DO M3=1,NUMP(N2)
              N1 = PATH3(ISPEC,M3)
              B1 = B(PATH3(ISPEC,M3))
              IF (B1.eq.B(ISPEC) .OR. B2.eq.B(ISPEC) .OR. 
     +            B3.eq.B(ISPEC)) then     
	      	write(*,*) '!WARNING: DENOMINATOR EQUALS TO 0!'
              END IF
              KSTAR1 = 1.d0/(-POR*DELTA*D(N1)*B1+Q)*(Q*fs(N1))
              KSTAR2 = 1.d0/(-POR*DELTA*D(N2)*B2+Q)*(Q*fs(N2)
     +                 -POR*DELTA*MU(N1)*Y(N2,M3)*ETA(N2,M3)
     +                 *KSTAR1*B1/((B1-B2)*(B1-A2))
     +                 +Q*MU(N1)*Y(N2,M3)*ETA(N2,M3)/D(N2)
     +                 *KSTAR1/((B1-B2)*(B1-A2)))
              KSTAR3 = -MU(N1)*Y(N2,M3)*ETA(N2,M3)/D(N2)
     +                 *KSTAR1/((B1-B2)*(B1-A2))
              KSTAR4 = 1.d0/(-POR*DELTA*D(N3)*B(N3)+Q)*(Q*fs(N3)
     +                 -POR*DELTA*MU(N2)*Y(N3,M2)*ETA(N3,M2)*(KSTAR2*B2
     +                 /((B2-B3)*(B2-A3))+KSTAR3*B1/((B1-B3)*(B1-A3)))              
     +                 +Q*MU(N2)*Y(N3,M2)*ETA(N3,M2)/D(N3)*(KSTAR2
     +                 /((B2-B3)*(B2-A3))+KSTAR3/((B1-B3)*(B1-A3))))               
	      KSTAR5 = -MU(N2)*Y(N3,M2)*ETA(N3,M2)/D(N3)
     +                 *KSTAR2/((B2-B3)*(B2-A3))
              KSTAR6 = -MU(N2)*Y(N3,M2)*ETA(N3,M2)/D(N3)
     +                 *KSTAR3/((B1-B3)*(B1-A3)) 
              KSTAR7 =1.d0/(-POR*DELTA*D(ISPEC)*B(ISPEC)+Q)*(Q*fs(ISPEC)
     +                 -POR*DELTA*MU(N3)*Y(ISPEC,M1)*ETA(ISPEC,M1)
     +                 *(KSTAR4*B3/((B3-B(ISPEC))*(B3-A(ISPEC)))+KSTAR5
     +                 *B2/((B2-B(ISPEC))*(B2-A(ISPEC)))+KSTAR6*B1
     +                 /((B1-B(ISPEC))*(B1-A(ISPEC))))+Q*MU(N3)
     +                 *Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)*(KSTAR4
     +                 /((B3-B(ISPEC))*(B3-A(ISPEC)))+KSTAR5
     +                 /((B2-B(ISPEC))*(B2-A(ISPEC)))+KSTAR6               
     +                 /((B1-B(ISPEC))*(B1-A(ISPEC)))))
              KSTAR8 =-MU(N3)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +                *KSTAR4/((B3-B(ISPEC))*(B3-A(ISPEC)))
              KSTAR9 =-MU(N3)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +                *KSTAR5/((B2-B(ISPEC))*(B2-A(ISPEC)))
              KSTAR10=-MU(N3)*Y(ISPEC,M1)*ETA(ISPEC,M1)/D(ISPEC)
     +                *KSTAR6/((B1-B(ISPEC))*(B1-A(ISPEC)))          
              COMP1=COMP1+KSTAR7
              COMP2=COMP2+KSTAR8*CDEXP(B3*XX)
              COMP3=COMP3+KSTAR9*CDEXP(B2*XX) 
              COMP4=COMP4+KSTAR10*CDEXP(B1*XX)  
            END DO
          END DO
        END DO
        COMP1=COMP1-(NUMP(ISPEC)-1)*Q*fs(ISPEC)/(-POR*DELTA*D(ISPEC)
     +        *B(ISPEC)+Q)        
        FBAR = COMP1*CDEXP(B(ISPEC)*XX)+COMP2+COMP3+COMP4
      END IF
c
      RETURN
      END
C
C     *******************
      INCLUDE 'HOOGD.FOR'
      include 'functions.for'
