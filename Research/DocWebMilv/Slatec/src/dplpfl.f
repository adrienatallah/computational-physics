*DECK DPLPFL
      SUBROUTINE DPLPFL (MRELAS, NVARS, IENTER, ILEAVE, IBASIS, IND,
     +   IBB, THETA, DIRNRM, RPRNRM, CSC, WW, BL, BU, ERP, RPRIM,
     +   PRIMAL, FINITE, ZEROLV)
C***BEGIN PROLOGUE  DPLPFL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DSPLP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SPLPFL-S, DPLPFL-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
C     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
C
C     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
C     /REAL (12 BLANKS)/DOUBLE PRECISION/.
C
C     THIS SUBPROGRAM IS PART OF THE DSPLP( ) PACKAGE.
C     IT IMPLEMENTS THE PROCEDURE (CHOOSE VARIABLE TO LEAVE BASIS).
C     REVISED 811130-1045
C     REVISED YYMMDD-HHMM
C
C***SEE ALSO  DSPLP
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   811215  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890605  Removed unreferenced labels.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DPLPFL
      INTEGER IBASIS(*),IND(*),IBB(*)
      DOUBLE PRECISION CSC(*),WW(*),BL(*),BU(*),ERP(*),RPRIM(*),
     * PRIMAL(*),BOUND,DIRNRM,RATIO,RPRNRM,THETA,ZERO
      LOGICAL FINITE,ZEROLV
C***FIRST EXECUTABLE STATEMENT  DPLPFL
      ZERO=0.D0
C
C     SEE IF THE ENTERING VARIABLE IS RESTRICTING THE STEP LENGTH
C     BECAUSE OF AN UPPER BOUND.
      FINITE=.FALSE.
      J=IBASIS(IENTER)
      IF (.NOT.(IND(J).EQ.3)) GO TO 20002
      THETA=BU(J)-BL(J)
      IF(J.LE.NVARS)THETA=THETA/CSC(J)
      FINITE=.TRUE.
      ILEAVE=IENTER
C
C     NOW USE THE BASIC VARIABLES TO POSSIBLY RESTRICT THE STEP
C     LENGTH EVEN FURTHER.
20002 I=1
      N20005=MRELAS
      GO TO 20006
20005 I=I+1
20006 IF ((N20005-I).LT.0) GO TO 20007
      J=IBASIS(I)
C
C     IF THIS IS A FREE VARIABLE, DO NOT USE IT TO
C     RESTRICT THE STEP LENGTH.
      IF (.NOT.(IND(J).EQ.4)) GO TO 20009
      GO TO 20005
C
C     IF DIRECTION COMPONENT IS ABOUT ZERO, IGNORE IT FOR COMPUTING
C     THE STEP LENGTH.
20009 IF (.NOT.(ABS(WW(I)).LE.DIRNRM*ERP(I))) GO TO 20012
      GO TO 20005
20012 IF (.NOT.(WW(I).GT.ZERO)) GO TO 20015
C
C     IF RPRIM(I) IS ESSENTIALLY ZERO, SET RATIO TO ZERO AND EXIT LOOP.
      IF (.NOT.(ABS(RPRIM(I)).LE.RPRNRM*ERP(I))) GO TO 20018
      THETA=ZERO
      ILEAVE=I
      FINITE=.TRUE.
      GO TO 20008
C
C     THE VALUE OF RPRIM(I) WILL DECREASE ONLY TO ITS LOWER BOUND OR
C     ONLY TO ITS UPPER BOUND.  IF IT DECREASES TO ITS
C     UPPER BOUND, THEN RPRIM(I) HAS ALREADY BEEN TRANSLATED
C     TO ITS UPPER BOUND AND NOTHING NEEDS TO BE DONE TO IBB(J).
20018 IF (.NOT.(RPRIM(I).GT.ZERO)) GO TO 10001
      RATIO=RPRIM(I)/WW(I)
      IF (.NOT.(.NOT.FINITE)) GO TO 20021
      ILEAVE=I
      THETA=RATIO
      FINITE=.TRUE.
      GO TO 20022
20021 IF (.NOT.(RATIO.LT.THETA)) GO TO 10002
      ILEAVE=I
      THETA=RATIO
10002 CONTINUE
20022 CONTINUE
      GO TO 20019
C
C     THE VALUE RPRIM(I).LT.ZERO WILL NOT RESTRICT THE STEP.
10001 CONTINUE
C
C     THE DIRECTION COMPONENT IS NEGATIVE, THEREFORE THE VARIABLE WILL
C     INCREASE.
20019 GO TO 20016
C
C     IF THE VARIABLE IS LESS THAN ITS LOWER BOUND, IT CAN
C     INCREASE ONLY TO ITS LOWER BOUND.
20015 IF (.NOT.(PRIMAL(I+NVARS).LT.ZERO)) GO TO 20024
      RATIO=RPRIM(I)/WW(I)
      IF (RATIO.LT.ZERO) RATIO=ZERO
      IF (.NOT.(.NOT.FINITE)) GO TO 20027
      ILEAVE=I
      THETA=RATIO
      FINITE=.TRUE.
      GO TO 20028
20027 IF (.NOT.(RATIO.LT.THETA)) GO TO 10003
      ILEAVE=I
      THETA=RATIO
10003 CONTINUE
20028 CONTINUE
C
C     IF THE BASIC VARIABLE IS FEASIBLE AND IS NOT AT ITS UPPER BOUND,
C     THEN IT CAN INCREASE TO ITS UPPER BOUND.
      GO TO 20025
20024 IF (.NOT.(IND(J).EQ.3 .AND. PRIMAL(I+NVARS).EQ.ZERO)) GO TO 10004
      BOUND=BU(J)-BL(J)
      IF(J.LE.NVARS) BOUND=BOUND/CSC(J)
      RATIO=(BOUND-RPRIM(I))/(-WW(I))
      IF (.NOT.(.NOT.FINITE)) GO TO 20030
      ILEAVE=-I
      THETA=RATIO
      FINITE=.TRUE.
      GO TO 20031
20030 IF (.NOT.(RATIO.LT.THETA)) GO TO 10005
      ILEAVE=-I
      THETA=RATIO
10005 CONTINUE
20031 CONTINUE
      CONTINUE
10004 CONTINUE
20025 CONTINUE
20016 GO TO 20005
20007 CONTINUE
C
C     IF STEP LENGTH IS FINITE, SEE IF STEP LENGTH IS ABOUT ZERO.
20008 IF (.NOT.(FINITE)) GO TO 20033
      ZEROLV=.TRUE.
      I=1
      N20036=MRELAS
      GO TO 20037
20036 I=I+1
20037 IF ((N20036-I).LT.0) GO TO 20038
      ZEROLV=ZEROLV.AND. ABS(THETA*WW(I)).LE.ERP(I)*RPRNRM
      IF (.NOT.(.NOT. ZEROLV)) GO TO 20040
      GO TO 20039
20040 GO TO 20036
20038 CONTINUE
20039 CONTINUE
20033 CONTINUE
      RETURN
      END
