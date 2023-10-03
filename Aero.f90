!PROGRAM VLM
!VORTEX LATTICE METHOD FOR UNCHEMBERED,SWEEP AND TAPERED WINGS
    DIMENSION GAMA(100),GAM(100)
     COMMON DX,DY,AR,N1,ALPHA,CT,B
     COMMON /COF/ A(100,101),NS
    ! OPEN(3,FILE='X1',STATUS='NEW')
    ! OPEN(4,FILE='Y1',STATUS='NEW')
    ! OPEN(5,FILE='Z1',STATUS='NEW')
    ! OPEN(6,FILE='A',STATUS='NEW')
    ! OPEN(8,FILE='X2',STATUS='NEW')
    ! OPEN(9,FILE='X3',STATUS='NEW')
    ! OPEN(7,FILE='X4',STATUS='NEW')
    ! OPEN(10,FILE='X5',STATUS='NEW')
    ! OPEN(11,FILE='X6',STATUS='NEW')
    ! OPEN(12,FILE='X7',STATUS='NEW')
    ! OPEN(13,FILE='X8',STATUS='NEW')
    ! OPEN(14,FILE='X9',STATUS='NEW')

!WE HAVE TO SUPPLY ASPECT RATIO,SWEPT ANGLE,TAPPER RATIO AND ANGLE OF ATTACK IN DEGREES(ALPHA)

        PI=3.141592
        WRITE(*,*)'GIVE THE VALUE OF AR'
        READ(*,*)AR
        WRITE(*,*)'GIVE THE VALUE OF N1,N2'
        READ(*,*)NX,NY
        WRITE(*,*)'GIVE THE VALUE OF ALPHA,SWEEP'
        READ(*,*)ALPHA,SWEEP
        WRITE(*,*)'GIVE THE VALUE OF CT'
        READ(*,*)CT
        B=AR
        DX=1.0/FLOAT(N1)
 36     DY=AR/(2.0*N2+0.5)
        NS=N1*N2
        WRITE(*,*)'DX= ', DX,'  DY= ',DY
        Y1=(0.5*B-0.25*DY)*TAN(SWEEP*PI/180.0)
        COSALF=COS(ALPHA*PI/180.0)
        SINALF=SIN(ALPHA*PI/180.0)
        WRITE(*,46)
 46     FORMAT (//,T2,'*',T3,35('-'),T38,'*')
        WRITE(*,45)
 45     FORMAT(T2,'LATTICE LOCATION',T24,'*',T30,'GAMA',T38,'*'/T2,'*',T3,35('-'),T38,'*')
        WRITE(8,47)
 47     FORMAT(//,T60,'*',T61,35('-'),T96,'*','*'//T60,'*',T61,'LATTICE LOCATION','*'//T82,'*',T88, &
        'GAMA',T96,'*'//T60,'*',T61,35('-'),T96,'*')

!SET CO-EFFICIENTS OF EQUATIONS FOR VORTEX STRENGTHS

    do 50  I = 1, N2
        do 50 J = 1, N1
            M1=(I-1)*N1+J
            A(M1,NS+1)=SINALF
            do 50 k = 1,N2
            do 50 L = 1,N1
            M2=(K-1)*N1+l
            call DW(I,J,K,L,A1,PI)
            A(M1,M2)=A1
    50 CONTINUE;

!SOLVE FOR VORTEX STRENGTH

    NNN=1
    call GE(NNN)
    do 60 I = 1,N2
        do 60 J = 1,N1
            M1=(I-1)*N1+J
            GAMA(M1)=A(M1,NS+1)
            GAM(M1)=ABS(GAMA (M1))
        WRITE(*,9)I,J,GAMA(M1)
        WRITE(3,8)I,J,GAMA(M1)
    9 FORMAT(//,T8,'*',T10,I2,T12,',',T13,I3,T25,'*', T26, F8.4,T36,'*')
    8 FORMAT(//,T63,'*', T66,I2,T68,',',T69,I3,T83,'*',T84,F8.4,T94,'*')
    60 CONTINUE
    WRITE(*,61)
    61 FORMAT(//,T2,'*',T3,38('-'),T40,'*')
    WRITE(9,62)
    62 FORMAT(//,T60,'*',T61,38('-'),T98,'*')
    WRITE(*,82)
    WRITE(*,81)
    82 FORMAT(//,T2,'*',T3,77('-'),T80,'*')
    81 FORMAT(T7,'Y',T15,'CL',T25,'CD',T35,'CM',T47,'CLT',T64,'CDT',T75,'CMT'/T2,'*',T3,77('-'),T80,'*')
    WRITE(7,83)
    83 FORMAT(//,T40,'*',T41,77('-'),T118,'*'T45,'Y',T53,'CL',&
    T63,'CD',T73,'CM',T85,'CLT',T102,'CDT',T113,'CMT'/T40,'*',T41,77('-'),T118,'*')

! COMPUTE FORCE AND MOMENT CO-EFFICIENTS
    CMT=0.0
    CDT=0.0
    CLT=0.0
    do 70 I=1,N2
        CX = 0.0
        CZ = 0.0
        CM = 0.0
        do 80 J=1,N1
            M1 = (I-1)*N1+J
            WX=0.0
            do 90 K=1,N2
                do 90 L=1,N1
                    M2=(K-1)*N1+l
                    call DW1(K,L,I,J,DELW,PI)
                    WX=WX+DELW*GAMA(M2)
                90 continue
                CX=CX+GAMA(M1)*(WX-SINALF)*2
                CZ=CZ+GAMA(M1)*COSALF*2
                CM=CM-GAMA(M1)*DX*(j-0.75)*COSALF*2
        80 continue
        CL=CZ*COSALF-CX*SINALF
        CD=CZ*SINALF+CX*COSALF
        CLT=CLT+CL*DY*2.0/AR
        CDT=CDT+CD*DY*2.0/AR
        CMT=CMT+CM*DY*2.0/AR
        XCPY=-CM/CL
        Y=(I-0.5)*DY/(0.5*AR)
        WRITE(*,100)Y,CL,CD,CM,CLT,CDT,CMT
        WRITE(4,107)Y,CL,CD,CM,CLT,CDT,CMT
    70 continue
    WRITE(10,103)
    103 FORMAT(//,T40,'*',T41,77('-'),T118,'*')
    WRITE(*,102)
    102 FORMAT(//,T2,'*',T3,77('-'),T80,'*')
    XCP=-CMT/CLT
    CO2=CDT/CLT**2
    WRITE(*,85)
    WRITE(*,91)
    85 FORMAT(//,T2,'*',T3,56('-'),T59,'*')
    91 FORMAT(T2,'LIFT COEFFICIENT CLT',T30,'INDUCED DRAG COFT. CDTI')
    WRITE(*,84)
    84 FORMAT(T2,'*',T3,56('-'),T59,'*')
    WRITE(11,98)
    98 FORMAT(//,T40,'*',T41,56('-'),T97,'*',T40,'LIFT COEFFICENT CLT',T68,'INDUCED DRAG COFT. CDTI',T40,'*',T41,56('-'),T97,'*')
    WRITE(13,988)
    988 FORMAT(//,T40,'*',T41,56('-'),T97,'*')
    WRITE(*,110)CLT,CDT
    WRITE(*,121)CLT,CDT
    WRITE(*,115)
    115 FORMAT(//,T2,'*',T3,44('-'),T46,'*')
    WRITE(*,999)
    WRITE(*,112)
    999 FORMAT(////,T2,'*',T3,36('-'),T39,'*')
    112 FORMAT(T5,'CO2',T20,'CMT',T30,'XCP'/T2,'*',36('-'),T39,'*')
    WRITE(12,120)
    120 FORMAT(////,T40,'*',T41,36('-'),T77,'*',T43,'CO2',T58,'CMT',T68,'XCP'/T40,'*',36('-'),T77,'*')
    WRITE(14,420)
    420 FORMAT(////,T40,'*',T41,36('-'),T77,'*')
    WRITE(*,111)CDOCL2,CMT,XCP
    WRITE(6,122)CDOCL2,CMT,XCP
    WRITE(*,117)
    117 FORMAT(//,T2,'*',T3,34('-'),T37,'*')
    100 FORMAT(//,T2,'*',T4,F6.3,T11,'*',T13,F6.3,T20,'*',T22,F6.3,T30,'*',&
    T32,F6.3,T42,'*',T44,F6.3,T59,'*',T61,T61,F6.3,T67,'*',T71,'*',T73,F6.3,T80,'*')
    107 FORMAT(//T40,'*',T42,F6.3,T49,'*',T51,F6.3,T58,'*',T60,F6.3,T68,'*',&
    T70,F6.3,T80,'*',T82,F6.3,T97,'*',T99,F6.3,T107,'*',T109,'*',T111,F6.3,T118,'*')
    110 FORMAT(//,T2,'*',T3,F6.3,T11,'*',T36,'*',T38,F6.3,T44,'*')
    121 FORMAT(//,T40,'*',T41,F6.3,T49,'*',T74,'*'T76,F6.3,T82,'*')
    111 FORMAT(//,T2,'*',T3,F6.3,T16,'*',T18,F6.3,T26,'*',T28,F6.3,T36,'*')
    122 FORMAT(//,T40,'*',T41,F6.3,T54,'*',T56,F6.3,T64,'*',T66,F6.3)
    ! CLOSE(3)
    ! CLOSE(4)
    ! CLOSE(5)
    ! CLOSE(6)
    ! CLOSE(7)
    ! CLOSE(8)
    ! CLOSE(9)
    ! CLOSE(10)
    ! CLOSE(11)
    STOP
    END

 SUBROUTINE DW(I,J,K,L,W,PI)
!COMPUTE DOWN WASH ON PANEL CENTERED AT(L-0.5)DX(DX1),(K-0.5)DY
!DUE TO VORTICES AT PANELS CENTRED AT (J-0.5)DX, +(I-0.5)DY
    COMMON DX,DY,AR,N1,SWEEP,CT,B
    CT1=CT-(DX*0.25)
    CT2=CT1+DY*0.25*TAN(SWEEP*PI/180.0)
    DCT2=CT2/FLOAT(N1)
    Y1=(0.5*B-0.25*DY)*TAN(SWEEP*PI/180.0)
    SHI = ATAN((Y1+(J-1)*DCT2-(J-1)*DX)/(0.5*B-0.25*DY))
    SII = ATAN((Y1+(2*J-1)*DCT2*0.5-0.5*(2*J-1)*DX)/(0.5*B-0.25*DY))
    SHIR=SHI
    SIIR=SII
    XL=DX*(J-0.75)+(I-1)*TAN(SHIR)*DY
    XR=DX*(J-0.75)+I*TAN(SHIR)*DY
    YL=DY*(I-1)
    YR=I*DY
    XP=DX*(L-0.75)+TAN(SHIR)*(2*I-1)*DY*0.5
    XP=DX*(L-0.25)+TAN(SIIR)*(2*I-1)*DY*0.5
    YP=DY*(k-0.5)

!CALCULATING THE DOWN WASH
        WW1= WH(XP,XL)
        WW2=WH(XP,XR)
        WW3=WH(XR,XL)
        WW4=WV(YP,YL)
        WW5=WV(YP,YR)
        WW6=WV(YR,YL)
        WW7=WV(YP,-YL)
        WW8=WV(YP,-YR)
        WW9=WH(-YR,-YL)
        W15=WW1/WW4
        W16=WW2/WW5
        W17=WW3/WW6
        W18=WW1/WW7
        W19=WW2/WW8
        W20=WW3/WW9
        W9=WH(XP,XL)**2+WV(YP,-YL)**2
        WO9=SQRT(W9)
        W10=WH(XP,XR)**2+WV(YP,-YR)**2
        WO10=SQRT(W10)
        W1=WH(XP,XL)**2+WV(YP,YL)**2
        W11=SQRT(W1)
        W2=WH(XP,XR)**2+WV(YP,YR)**2
        W22=SQRT(W2)
        W7=WV(YP,YR)
        W5=WV(YP,YL)
        IF(W15.EQ.W16)GO TO 130
        IF(W15.EQ.W17)GO TO 130
        IF(W16.EQ.W17)GO TO 130
        IF(W18.EQ.W19)GO TO 130
        IF(W18.EQ.W20)GO TO 130
        IF(W19.EQ.W20)GO TO 130
        W3=WH(XP,XL)*WV(YR,YL)-WH(XR,XL)*WV(YP,YL)
        W4=WV(YR,YL)
        WO1=(1.0+W4*W11/W3)/W5
        W6=WH(XP,XR)*WV(YR,YL)-WH(XR,XL)*WV(YP,YR)
        WO2=(1.+W4*W22/W6)/W7
        W8=WV(-YR,-YL)
        W111=WH(XP,XL)*WV(-YR,-YL)-WH(XR,XL)*WV(YP,-YL)
        W12=WV(YP,-YL)
        WO3=(1.0+W8*WO9/W111)/W12
        W13=WH(XP,XR)*WV(-YR,-YL)-WH(XR,XL)*WV(YP,-YR)
        W14=WV(YP,-YR)
        WO4=(1.0+W8*WO10/W13)/W14
        WO=WO1-WO2+WO3-WO4
        130 WO=(1.0+WH(XP,XL)/W11)/W5-(1.0+WH(XP,XR)/W22)/W7+(1.0+WH(XP,XL)/WO9)/WV(YP,-YL)-(1.0+WH(XP,XR)/WO10)/WV(YP,-YR)
        W=WO*0.25/PI
        RETURN
        END

        SUBROUTINE DW1(I,J,K,L,W,PI)
!COMPUTE DOWNWASH ON PANEL CENTERED AT(L-0.5)DX(DX1),(K-0.5)DY
!DUE TO VORTRICES AT PANELS CENTERED AT(J-0.5)DX,+-(I-0.5)DY

        COMMON DX,DY,AR,N1,SWEEP,CT,B
        CT1=CT-DX*0.25
        CT2=CT1+DY*0.25*TAN(SWEEP*PI/180.0)
        DCT2=CT2/FLOAT(N1)
        Y1=(0.5*B-0.25*DY)*TAN(SWEEP*PI/180.0)
        SHI=ATAN(Y1+(J-1)*DCT2-(J-1)*DX)/(0.5*B-0.25*DY)
        SII=ATAN((Y1+(2*J-1)*DCT2*0.5-0.5*(2*J-1)*DX)/(0.5*B-0.25*DY))
        SHIR=SHI
        SIIR=SII
        XL=DX*(J-0.75)+(I-1)*TAN(SHIR)*DY
        XR=DX*(J-0.75)+I*TAN(SHIR)*DY
        YL=DY*(I-1)
        YR=I*DY
        XP=DX*(L-0.25)+TAN(SIIR)*(2*I-1)*DY*0.5
            XP=DX*(L-0.75)+TAN(SIIR)*(2*I-1)*DY*0.5
            YP=DY*(K-0.5)

!CALCULATING DOWN WASH
        WW1=WH(XP,XL)
        WW2=WH(XP,XR)
        WW3=WH(YP,YL)
        WW4=WH(XP,YL)
        WW5=WH(YP,YR)
        WW6=WH(YR,YL)
        WW7=WH(YP,-YL)
        WW8=WH(YP,-YR)
        WW9=WH(-YR,-YL)
        W15=WW1/WW4
        W16=WW2/WW5
        W17=WW3/WW6
        W18=WW1/WW7
        W19=WW2/WW8
        W20=WW3/WW9
        W9=WH(XP,XL)**2+WV(YP,-YL)**2
        WO9=SQRT(W9)
        W10=WH(XP,XR)**2+WV(YP,-YR)**2
        WO10=SQRT(W10)
        W1=WH(XP,XL)**2+WV(YP,YL)**2
        W11=SQRT(W1)
        W2=WH(XP,XR)**2+WV(YP,YR)**2
        W22=SQRT(W2)
        W7=WV(YP,YR)
        W5=WV(YP,YL)
        IF(W15.EQ.W16) GO TO 130
        IF(W15.EQ.W17) GO TO 130
        IF(W16.EQ.W17) GO TO 130
        IF(W18.EQ.W19) GO TO 130
        IF(W18.EQ.W20) GO TO 130
        IF(W19.EQ.W20) GO TO 130
        W3=WH(XP,XL)*WV(YR,YL)-WH(XR,XL)*WV(YP,YL)
        W4=WV(YR,YL)
        WO1=(1.0+W4*W11/W3)/W5
        W6=WH(XP,XR)*WV(YR,YL)-WH(XR,XL)*WV(YP,YR)
        WO2=(1.+W4*W22/W6)/W7
        W8=WV(-YR,-YL)
        W111=WH(XP,XL)*WV(-YR,-YL)-WH(XR,XL)*WV(YP,-YL)
        W12=WV(YP,-YL)
        WO3=(1.0+W8*WO9/W111)/W12
        W13=WH(XP,XR)*WV(-YR,-YL)-WH(XR,XL)*WV(YP,-YR)
        W14=WV(YP,-YR)
        WO4=(1.0+W8*WO10/W13)/W14
        WO=WO1-WO2+WO3-WO4
 130    WO=(1.0+WH(XP,XL)/W11)/W5-(1.0+WH(XP,XR)/W22)/W7+(1.0+WH(XP,XL)/WO9)/WV(YP,-YL)-(1.0+WH(XP,XR)/WO10)/WV(YP,-YR)
        W=WO*0.25/PI
        RETURN
        END
!GAUSS ELIMINATION METHOD TO SOLVE ALGEBRIC EQUATION
        SUBROUTINE GE(NRHS)
!(A)=COEFFICIENT MATRIX
!NS=NUMBER OF EQUATIONS
!NRHS=NUMBER OF RIGHT HAND SIDES
!RIGHT HAND SIDES AND SOLUTIONS STORED IN COLUMNS NEQUANS+1
!THRU NSEQANS+NRHS OF (A)

        COMMON/COF/A(100,101),NS
        NP=NS+1
        NTOT=NS+NRHS
!GAUSS REDUCTION
        DO 151 I=2,NS
!SEARCH FOR LARGEST ENTRY IN(I-1)TH COLUMN ON OR BELOW MAIN DIAGONAL
        IM=I-1
        IMAX=IM
        IMAX=IM
        BIGG=ABS(A(IM,IM))
        DO 111 J=I,NS
        IF(BIGG.GE.ABS(A(J,IM))) GO TO 111
        IMAX= J
        BIGG=ABS(A(J,IM))
 111    CONTINUE

!SWITCH(I-1)TH AND IMAXTH EQUATIONS
         IF(IMAX.NE.IM)GO TO 141
         DO 131 J = IM ,NTOT
         TEMP=A(IM,J)
         A(IM,J)=A(IMAX,J)
         A(IMX,J)=TEMP
 131     CONTINUE
!ELIMINATE (I-1)TH UNKNOWN FROM ITH THRU (NEQNS)TH EQUATIONS
 141     DO 151 J=I,NS
         R=A(J,IM)/A(IM,IM)
            DO 151 K=I,NTOT
 151        A(J,K)=A(J,K)-R*A(IM,K)
!BACK SUBSTITUTION
         DO 221 K=NP,NTOT
         A(NS,K)=A(NS,K)/A(NS,NS)
         DO 211 L=2,NS
         I=NS+1-L
         IP=I+1
         DO 201 J=IP,NS
 201     A(I,K)=A(I,K)-A(I,J)*A(J,K)
 211     A(I,K)=A(I,K)/A(I,I)
 221     CONTINUE
         RETURN
         END

         FUNCTION WH(X01,X02)
         WH=X01-X02
         RETURN
         END

         FUNCTION WV(Y01,Y02)
         WV=Y01-Y02
         RETURN
         END