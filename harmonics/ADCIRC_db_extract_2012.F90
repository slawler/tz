!-------------------------------------------------------------------------------
!
!                         ADCIRC_db_extract_2012.F90
! 
!                               modified from  
!
!                             ADCIRC_db_extract.F
!                            last mod: Jesse Feyen
!                                Version 2.0
!
!                 added modified KDTREE adcirc routines from
!
!                              GridCompare.F90
!                              Zachary Cobell
!                            Version 10.0_alpha4
!
!                                    by
!
!                             Christine Szpilka
!                          University of Oklahoma
!                            cmszpilka@ou.edu
!                         Last Modified: 04/05/2013
!
! Program to extract tidal information from ADCIRC 2DDI tidal databases
! for use as boundary conditions in another model run.
!
! This version has been updated to use the KDTREE2 set of algorithms for
! faster search times. It contains the KDTREE2 modules and can be compiled 
! as a stand-alone code with no dependancies.
!
! All hard-coded array dimensions and filenames have been updated to user
! input and dynamic allocation for ease of use and flexibility. 
! Original comments and history included below as appropriate.
!
! The ADCIRC_db_extract GLOBAL module, which contains all of the subroutines
! used to find the element, calculate interpolation factors, and compute
! harmonic values at the extraction points, begins on line 1449 after the
! KDTREE2 modules.
!
! The ADCIRC_db_extract_2012 main code, which contains the input and output
! and calls subroutines from GLOBAL, begins on line 1757.
!
!-------------------------------------------------------------------------------
! Requires four files as input:
!
! 1) infile     input file containing list of lon/lat locations for extraction
! 2) fort.14    ADCIRC mesh used to create the tidal database
! 3) fort.53    Global elevation harmonic output for database (opt: if elevs = TRUE)
! 4) fort.54    Global velocity harmonic output for database (opt: if vels = TRUE)
!
! Input file descriptions and format
!
!          UNIT 12 : infile : file listing the X and Y coordinates of all
!                  locations for extracted harmonics in lon/lat format
!                  
!		   NOUT = number of points in file
!                  XOUT(i),YOUT(i)  i = 1, NOUT
!
!          UNIT 14 : ADCIRC finite element grid file used for database
!
!			 AGRID = Alphanumeric grid description (<= 24 characters)
!			 NE, NP = number of elements and nodes in mesh (respectively)
!			 JKI,X(JKI),Y(JKI),DP(JKI)   JKI = 1, NP
!				node number, X and Y coordinates, bathymetry
!			 JKI,NHY,NM(JKI,1),NM(JKI,2),NM(JKI,3)  JKI=1,NE  
!				element number, type and connectivity (specified with a
!				counterclockwise orientation)
!
!          UNIT 53 : Global harmonic constituents for elevations at all NP nodes
!                    from fort.14 file.
!
!          		 NHARFR = number of constituents in file
!          		 For J=1,NHARF
!             	   FREQ, NFACT, EQARG, NAME  (only save NAME, rest read in as junk)
!			 end J loop
!          		 NP = number of nodes in grid (should match fort.14)
!                  J=1,NP
!			   J
!             	   K=1,NHARF
!                      HAMP(J,K),HPHASE(J,K)
!			   end K loop
!			 end J loop
!
!          UNIT 54 : Global harmonic constituents for velocities at all NP nodes
!                    from fort.14 file.
!
!          		 NHARFR = number of constituents in file
!          		 For J=1,NHARF
!             	   FREQ, NFACT, EQARG, NAME  (only save NAME, rest read in as junk)
!			 end J loop
!          		 NP = number of nodes in grid (should match fort.14)
!                  J=1,NP
!			   J
!             	   K=1,NHARF
!                      UAMP(J,K),UPHASE(J,K),VAMP(J,K),VPHASE(J,K)
!			   end K loop
!			 end J loop
!
!-------------------------------------------------------------------------------
!
! Creates two to three output files, depending upon input selections:
!
! 1) tides.dia    contains diagnostic messages from extraction program
! 2) elev_hc.out  contains extracted elevation harmonics for NOUT points (opt: if elevs = TRUE)
! 3) vel_hc.out   contains extracted velocity harmonics for NOUT points (opt: if vels = TRUE)
!
! Output file descriptions and format
!
!          UNIT 40 : tides.dia : output location and found element for each NOUT point and any
!			 error messages
!
!          UNIT 50 : elev_hc.out : elevation harmonics at extraction points
!
!          UNIT 60 : vel_hc.out : velocity harmonics at extraction points
!-------------------------------------------------------------------------------


!----------------------------------------------------------------------------------!
!                                                                                  !
!                           START OF KDTREE ROUTINES                               !
!                                                                                  !
!----------------------------------------------------------------------------------!
        MODULE KDTREE2_PRECISION_MODULE
          INTEGER, PARAMETER :: SP = KIND(0.0)
          INTEGER, PARAMETER :: DP = KIND(0.0D0)

          PRIVATE :: SP, DP
          !INTEGER, PARAMETER :: KDKIND = SP
          INTEGER, PARAMETER :: KDKIND = DP
          PUBLIC :: KDKIND

          END MODULE KDTREE2_PRECISION_MODULE

          MODULE KDTREE2_PRIORITY_QUEUE_MODULE
          USE KDTREE2_PRECISION_MODULE

          TYPE KDTREE2_RESULT
             REAL(KDKIND)    :: DIS!=0.0
             INTEGER :: IDX!=-1   INITIALIZERS CAUSE SOME BUGS IN COMPILERS.
          END TYPE KDTREE2_RESULT

          TYPE PQ
             INTEGER :: HEAP_SIZE = 0
             TYPE(KDTREE2_RESULT), POINTER :: ELEMS(:)
          END TYPE PQ

          PUBLIC :: KDTREE2_RESULT

          PUBLIC :: PQ
          PUBLIC :: PQ_CREATE
          PUBLIC :: PQ_DELETE, PQ_INSERT
          PUBLIC :: PQ_EXTRACT_MAX, PQ_MAX, PQ_REPLACE_MAX, PQ_MAXPRI
          PRIVATE

          CONTAINS


          FUNCTION PQ_CREATE(RESULTS_IN) RESULT(RES)

          TYPE(KDTREE2_RESULT), TARGET:: RESULTS_IN(:)
          TYPE(PQ) :: RES

          INTEGER :: NALLOC

          NALLOC = SIZE(RESULTS_IN,1)
          IF (NALLOC .LT. 1) THEN
             WRITE (*,*) 'PQ_CREATE: ERROR, ARRAYS MUST BE ALLOCATED.'
          END IF
          RES%ELEMS => RESULTS_IN
          RES%HEAP_SIZE = 0
          RETURN
          END FUNCTION PQ_CREATE

          SUBROUTINE HEAPIFY(A,I_IN)
          TYPE(PQ),POINTER   :: A
          INTEGER, INTENT(IN) :: I_IN
          !
          INTEGER :: I, L, R, LARGEST

          REAL(KDKIND)    :: PRI_I, PRI_L, PRI_R, PRI_LARGEST

          TYPE(KDTREE2_RESULT) :: TEMP

          I = I_IN

          BIGLOOP:  DO
             L = 2*I ! LEFT(I)
             R = L+1 ! RIGHT(I)

             IF (L .GT. A%HEAP_SIZE) THEN
                EXIT
             ELSE
                PRI_I = A%ELEMS(I)%DIS
                PRI_L = A%ELEMS(L)%DIS
                IF (PRI_L .GT. PRI_I) THEN
                   LARGEST = L
                   PRI_LARGEST = PRI_L
                ELSE
                   LARGEST = I
                   PRI_LARGEST = PRI_I
                ENDIF
                IF (R .LE. A%HEAP_SIZE) THEN
                   PRI_R = A%ELEMS(R)%DIS
                   IF (PRI_R .GT. PRI_LARGEST) THEN
                      LARGEST = R
                   ENDIF
                ENDIF
             ENDIF

             IF (LARGEST .NE. I) THEN
                TEMP = A%ELEMS(I)
                A%ELEMS(I) = A%ELEMS(LARGEST)
                A%ELEMS(LARGEST) = TEMP
                I = LARGEST
                CYCLE BIGLOOP
             ELSE
                RETURN
             END IF
          ENDDO BIGLOOP
          RETURN
          END SUBROUTINE HEAPIFY

          SUBROUTINE PQ_MAX(A,E)

          TYPE(PQ),POINTER :: A
          TYPE(KDTREE2_RESULT),INTENT(OUT)  :: E

          IF (A%HEAP_SIZE .GT. 0) THEN
             E = A%ELEMS(1)
          ELSE
             WRITE (*,*) 'PQ_MAX: ERROR, HEAP_SIZE < 1'
             STOP
          ENDIF
          RETURN
          END SUBROUTINE PQ_MAX

          REAL(KDKIND) FUNCTION PQ_MAXPRI(A)
          TYPE(PQ), POINTER :: A

          IF (A%HEAP_SIZE .GT. 0) THEN
             PQ_MAXPRI = A%ELEMS(1)%DIS
          ELSE
             WRITE (*,*) 'PQ_MAX_PRI: ERROR, HEAPSIZE < 1'
             STOP
          ENDIF
          RETURN
          END FUNCTION PQ_MAXPRI

          SUBROUTINE PQ_EXTRACT_MAX(A,E)
          TYPE(PQ),POINTER :: A
          TYPE(KDTREE2_RESULT), INTENT(OUT) :: E

          IF (A%HEAP_SIZE .GE. 1) THEN
             E = A%ELEMS(1)
             A%ELEMS(1) = A%ELEMS(A%HEAP_SIZE)
             A%HEAP_SIZE = A%HEAP_SIZE-1
             CALL HEAPIFY(A,1)
             RETURN
          ELSE
             WRITE (*,*) 'PQ_EXTRACT_MAX: ERROR,',&
                        ' ATTEMPTED TO POP NON-POSITIVE PQ'
             STOP
          END IF

          END SUBROUTINE PQ_EXTRACT_MAX


          REAL(KDKIND) FUNCTION PQ_INSERT(A,DIS,IDX)

          TYPE(PQ),POINTER  :: A
          REAL(KDKIND), INTENT(IN) :: DIS
          INTEGER, INTENT(IN) :: IDX
          INTEGER :: I, ISPARENT
          REAL(KDKIND)    :: PARENTDIS
          A%HEAP_SIZE = A%HEAP_SIZE + 1
          I = A%HEAP_SIZE

          DO WHILE (I .GT. 1)
             ISPARENT = INT(I/2)
             PARENTDIS = A%ELEMS(ISPARENT)%DIS
             IF (DIS .GT. PARENTDIS) THEN
                A%ELEMS(I)%DIS = PARENTDIS
                A%ELEMS(I)%IDX = A%ELEMS(ISPARENT)%IDX
                I = ISPARENT
             ELSE
                EXIT
             ENDIF
          END DO

          A%ELEMS(I)%DIS = DIS
          A%ELEMS(I)%IDX = IDX

          PQ_INSERT = A%ELEMS(1)%DIS
          RETURN

          END FUNCTION PQ_INSERT

          SUBROUTINE PQ_ADJUST_HEAP(A,I)
          TYPE(PQ),POINTER  :: A
          INTEGER, INTENT(IN) :: I
          REAL(KDKIND) :: PRICHILD
          INTEGER :: PARENT, CHILD, N

          TYPE(KDTREE2_RESULT) :: E

          E = A%ELEMS(I)

          PARENT = I
          CHILD = 2*I
          N = A%HEAP_SIZE

          DO WHILE (CHILD .LE. N)
             IF (CHILD .LT. N) THEN
                IF (A%ELEMS(CHILD)%DIS .LT. A%ELEMS(CHILD+1)%DIS) THEN
                   CHILD = CHILD+1
                ENDIF
             ENDIF
             PRICHILD = A%ELEMS(CHILD)%DIS
             IF (E%DIS .GE. PRICHILD) THEN
                EXIT
             ELSE
                A%ELEMS(PARENT) = A%ELEMS(CHILD)
                PARENT = CHILD
                CHILD = 2*PARENT
             END IF
          END DO
          A%ELEMS(PARENT) = E
          RETURN
          END SUBROUTINE PQ_ADJUST_HEAP


          REAL(KDKIND) FUNCTION PQ_REPLACE_MAX(A,DIS,IDX)
          TYPE(PQ),POINTER         :: A
          REAL(KDKIND), INTENT(IN) :: DIS
          INTEGER, INTENT(IN) :: IDX
          INTEGER :: PARENT, CHILD, N
          REAL(KDKIND)    :: PRICHILD, PRICHILDP1

          TYPE(KDTREE2_RESULT) :: ETMP

          IF (.TRUE.) THEN
             N=A%HEAP_SIZE
             IF (N .GE. 1) THEN
                PARENT =1
                CHILD=2

                LOOP: DO WHILE (CHILD .LE. N)
                   PRICHILD = A%ELEMS(CHILD)%DIS
                   IF (CHILD .LT. N) THEN
                      PRICHILDP1 = A%ELEMS(CHILD+1)%DIS
                      IF (PRICHILD .LT. PRICHILDP1) THEN
                        CHILD = CHILD+1
                        PRICHILD = PRICHILDP1
                      ENDIF
                   ENDIF

                   IF (DIS .GE. PRICHILD) THEN
                      EXIT LOOP
                   ELSE
                      A%ELEMS(PARENT) = A%ELEMS(CHILD)
                      PARENT = CHILD
                      CHILD = 2*PARENT
                   END IF
                END DO LOOP
                A%ELEMS(PARENT)%DIS = DIS
                A%ELEMS(PARENT)%IDX = IDX
                PQ_REPLACE_MAX = A%ELEMS(1)%DIS
             ELSE
                A%ELEMS(1)%DIS = DIS
                A%ELEMS(1)%IDX = IDX
                PQ_REPLACE_MAX = DIS
             ENDIF
          ELSE
             CALL PQ_EXTRACT_MAX(A,ETMP)
             ETMP%DIS = DIS
             ETMP%IDX = IDX
             PQ_REPLACE_MAX = PQ_INSERT(A,DIS,IDX)
          ENDIF
          RETURN
          END FUNCTION PQ_REPLACE_MAX

          SUBROUTINE PQ_DELETE(A,I)
          TYPE(PQ),POINTER :: A
          INTEGER           :: I

          IF ((I .LT. 1) .OR. (I .GT. A%HEAP_SIZE)) THEN
             WRITE (*,*) 'PQ_DELETE: ERROR, ATTEMPT TO REMOVE',&
                        ' OUT OF BOUNDS ELEMENT.'
             STOP
          ENDIF

          A%ELEMS(I) = A%ELEMS(A%HEAP_SIZE)
          A%HEAP_SIZE = A%HEAP_SIZE - 1

          CALL HEAPIFY(A,I)

          END SUBROUTINE PQ_DELETE

          END MODULE KDTREE2_PRIORITY_QUEUE_MODULE


          MODULE KDTREE2_MODULE

          USE KDTREE2_PRECISION_MODULE
          USE KDTREE2_PRIORITY_QUEUE_MODULE

          PUBLIC :: KDKIND
          PUBLIC :: KDTREE2, KDTREE2_RESULT, TREE_NODE
          PUBLIC :: KDTREE2_CREATE, KDTREE2_DESTROY
          PUBLIC :: KDTREE2_N_NEAREST,KDTREE2_N_NEAREST_AROUND_POINT
          PUBLIC :: KDTREE2_R_NEAREST, KDTREE2_R_NEAREST_AROUND_POINT
          PUBLIC :: KDTREE2_SORT_RESULTS
          PUBLIC :: KDTREE2_R_COUNT, KDTREE2_R_COUNT_AROUND_POINT
          PUBLIC :: KDTREE2_N_NEAREST_BRUTE_FORCE
          PUBLIC :: KDTREE2_R_NEAREST_BRUTE_FORCE

          INTEGER, PARAMETER :: BUCKET_SIZE = 12

          TYPE INTERVAL
          REAL(KDKIND) :: LOWER,UPPER
          END TYPE INTERVAL

          TYPE :: TREE_NODE
             PRIVATE
             INTEGER :: CUT_DIM
             REAL(KDKIND) :: CUT_VAL
             REAL(KDKIND) :: CUT_VAL_LEFT, CUT_VAL_RIGHT
             INTEGER :: L, U
             TYPE (TREE_NODE), POINTER :: LEFT, RIGHT
             TYPE(INTERVAL), POINTER :: BOX(:) => NULL()
          END TYPE TREE_NODE

          TYPE :: KDTREE2
             INTEGER :: DIMEN=0, N=0
             REAL(KDKIND), POINTER :: THE_DATA(:,:) => NULL()
             INTEGER, POINTER :: IND(:) => NULL()
             LOGICAL       :: SORT = .FALSE.
             LOGICAL       :: REARRANGE = .FALSE.
             REAL(KDKIND), POINTER :: REARRANGED_DATA(:,:) => NULL()
             TYPE (TREE_NODE), POINTER :: ROOT => NULL()
          END TYPE KDTREE2


          TYPE :: TREE_SEARCH_RECORD

             PRIVATE

             INTEGER           :: DIMEN
             INTEGER           :: NN, NFOUND
             REAL(KDKIND)      :: BALLSIZE
             INTEGER           :: CENTERIDX=999, CORRELTIME=9999
             INTEGER           :: NALLOC
             LOGICAL           :: REARRANGE
             LOGICAL           :: OVERFLOW
             REAL(KDKIND), POINTER :: QV(:)
             TYPE(KDTREE2_RESULT), POINTER :: RESULTS(:)
             TYPE(PQ) :: PQ
             REAL(KDKIND), POINTER :: DATA(:,:)
             INTEGER, POINTER      :: IND(:)
          END TYPE TREE_SEARCH_RECORD

          PRIVATE

          TYPE(TREE_SEARCH_RECORD), SAVE, TARGET :: SR

          CONTAINS

          FUNCTION KDTREE2_CREATE(INPUT_DATA,IDIM2,SORT,REARRANGE)RESULT(MR)
          TYPE (KDTREE2), POINTER :: MR
          INTEGER, INTENT(IN), OPTIONAL      :: IDIM2
          LOGICAL, INTENT(IN), OPTIONAL      :: SORT
          LOGICAL, INTENT(IN), OPTIONAL      :: REARRANGE
          REAL(KDKIND), TARGET :: INPUT_DATA(:,:)
          INTEGER :: I
          ALLOCATE (MR)
          MR%THE_DATA => INPUT_DATA

          IF (PRESENT(IDIM2)) THEN
             MR%DIMEN = IDIM2
          ELSE
             MR%DIMEN = SIZE(INPUT_DATA,1)
          END IF
          MR%N = SIZE(INPUT_DATA,2)

          IF (MR%DIMEN > MR%N) THEN
             WRITE (*,*) 'KD_TREE_TRANS: LIKELY USER ERROR.'
             WRITE (*,*) 'KD_TREE_TRANS: YOU PASSED IN MATRIX WITH D=',&
                         MR%DIMEN
             WRITE (*,*) 'KD_TREE_TRANS: AND N=',MR%N

             WRITE (*,*) 'KD_TREE_TRANS: NOTE, THAT NEW FORMAT IS',&
                        ' DATA(1:D,1:N)'
             WRITE (*,*) 'KD_TREE_TRANS: WITH USUALLY N >> D.  ',&
                        'IF N =APPROX= D, THEN A K-D TREE'
             WRITE (*,*) 'KD_TREE_TRANS: IS NOT AN APPROPRIATE DATA',&
                        ' STRUCTURE.'
             STOP
          END IF

          CALL BUILD_TREE(MR)

          IF (PRESENT(SORT)) THEN
             MR%SORT = SORT
          ELSE
             MR%SORT = .FALSE.
          ENDIF

          IF (PRESENT(REARRANGE)) THEN
             MR%REARRANGE = REARRANGE
          ELSE
             MR%REARRANGE = .TRUE.
          ENDIF

          IF (MR%REARRANGE) THEN
             ALLOCATE(MR%REARRANGED_DATA(MR%DIMEN,MR%N))
             DO I=1,MR%N
             MR%REARRANGED_DATA(:,I) = MR%THE_DATA(:,&
                   MR%IND(I))
             ENDDO
          ELSE
             NULLIFY(MR%REARRANGED_DATA)
          ENDIF

          END FUNCTION KDTREE2_CREATE

          SUBROUTINE BUILD_TREE(TP)
             TYPE (KDTREE2), POINTER :: TP
             ! ..
             INTEGER :: J
             TYPE(TREE_NODE), POINTER :: DUMMY => NULL()
             ! ..
             ALLOCATE (TP%IND(TP%N))
             FORALL (J=1:TP%N)
             TP%IND(J) = J
             END FORALL
             TP%ROOT => BUILD_TREE_FOR_RANGE(TP,1,TP%N, DUMMY)
          END SUBROUTINE BUILD_TREE

          RECURSIVE FUNCTION BUILD_TREE_FOR_RANGE(TP,L,U,PARENT)&
           RESULT (RES)
             TYPE (TREE_NODE), POINTER :: RES
             TYPE (KDTREE2), POINTER :: TP
             TYPE (TREE_NODE),POINTER           :: PARENT
             INTEGER, INTENT (IN) :: L, U
             INTEGER :: I, C, M, DIMEN
             LOGICAL :: RECOMPUTE
             REAL(KDKIND)    :: AVERAGE
             DIMEN = TP%DIMEN
             ALLOCATE (RES)
             ALLOCATE(RES%BOX(DIMEN))

             IF ( U < L ) THEN
               NULLIFY(RES)
               RETURN
             END IF

             IF ((U-L)<=BUCKET_SIZE) THEN
             DO I=1,DIMEN
                CALL SPREAD_IN_COORDINATE(TP,I,L,U,RES%BOX(I))
             END DO
             RES%CUT_DIM = 0
             RES%CUT_VAL = 0.0
             RES%L = L
             RES%U = U
             RES%LEFT =>NULL()
             RES%RIGHT => NULL()
             ELSE
             DO I=1,DIMEN
                RECOMPUTE=.TRUE.
                IF (ASSOCIATED(PARENT)) THEN
                  IF (I .NE. PARENT%CUT_DIM) THEN
                     RECOMPUTE=.FALSE.
                  END IF
                ENDIF
                IF (RECOMPUTE) THEN
                  CALL SPREAD_IN_COORDINATE(TP,I,L,U,RES%BOX(I))
                ELSE
                  RES%BOX(I) = PARENT%BOX(I)
                ENDIF
             END DO

             C = MAXLOC(RES%BOX(1:DIMEN)%UPPER-RES%BOX(1:DIMEN)%LOWER,1)

             IF (.FALSE.) THEN
                M = (L+U)/2
                CALL SELECT_ON_COORDINATE(TP%THE_DATA,TP%IND,C,M,L,U)
             ELSE
                IF (.TRUE.) THEN
                AVERAGE = SUM(TP%THE_DATA(C,TP%IND(L:U))) / &
                 REAL(U-L+1,KDKIND)
                ELSE
                AVERAGE = (RES%BOX(C)%UPPER + RES%BOX(C)%LOWER)/2.0
                ENDIF

                RES%CUT_VAL = AVERAGE
                M = SELECT_ON_COORDINATE_VALUE(&
                                        TP%THE_DATA,TP%IND,C,AVERAGE,L,U)
             ENDIF

             RES%CUT_DIM = C
             RES%L = L
             RES%U = U

             RES%LEFT => BUILD_TREE_FOR_RANGE(TP,L,M,RES)
             RES%RIGHT => BUILD_TREE_FOR_RANGE(TP,M+1,U,RES)

             IF (ASSOCIATED(RES%RIGHT) .EQV. .FALSE.) THEN
                RES%BOX = RES%LEFT%BOX
                RES%CUT_VAL_LEFT = RES%LEFT%BOX(C)%UPPER
                RES%CUT_VAL = RES%CUT_VAL_LEFT
             ELSEIF (ASSOCIATED(RES%LEFT) .EQV. .FALSE.) THEN
                RES%BOX = RES%RIGHT%BOX
                RES%CUT_VAL_RIGHT = RES%RIGHT%BOX(C)%LOWER
                RES%CUT_VAL = RES%CUT_VAL_RIGHT
             ELSE
                RES%CUT_VAL_RIGHT = RES%RIGHT%BOX(C)%LOWER
                RES%CUT_VAL_LEFT = RES%LEFT%BOX(C)%UPPER
                RES%CUT_VAL = (RES%CUT_VAL_LEFT + RES%CUT_VAL_RIGHT)/2

                RES%BOX%UPPER = MAX(RES%LEFT%BOX%UPPER,RES%RIGHT%BOX%UPPER)
                RES%BOX%LOWER = MIN(RES%LEFT%BOX%LOWER,RES%RIGHT%BOX%LOWER)
             ENDIF
             END IF
          END FUNCTION BUILD_TREE_FOR_RANGE

          INTEGER FUNCTION SELECT_ON_COORDINATE_VALUE(V,IND,C,ALPHA,LI,UI)RESULT(RES)
             INTEGER, INTENT (IN) :: C, LI, UI
             REAL(KDKIND), INTENT(IN) :: ALPHA
             REAL(KDKIND) :: V(1:,1:)
             INTEGER :: IND(1:)
             INTEGER :: TMP
             INTEGER :: LB, RB

             LB = LI; RB = UI

             DO WHILE (LB < RB)
                IF ( V(C,IND(LB)) <= ALPHA ) THEN
                   LB = LB+1
                ELSE
                   TMP = IND(LB); IND(LB) = IND(RB); IND(RB) = TMP
                   RB = RB-1
                ENDIF
             END DO

             IF (V(C,IND(LB)) <= ALPHA) THEN
                RES = LB
             ELSE
                RES = LB-1
             ENDIF

          END FUNCTION SELECT_ON_COORDINATE_VALUE

          SUBROUTINE SELECT_ON_COORDINATE(V,IND,C,K,LI,UI)
             INTEGER, INTENT (IN) :: C, K, LI, UI
             INTEGER :: I, L, M, S, T, U
             REAL(KDKIND) :: V(:,:)
             INTEGER :: IND(:)
             L = LI
             U = UI
             DO WHILE (L<U)
             T = IND(L)
             M = L
             DO I = L + 1, U
                IF (V(C,IND(I))<V(C,T)) THEN
                M = M + 1
                S = IND(M)
                IND(M) = IND(I)
                IND(I) = S
                END IF
             END DO
             S = IND(L)
             IND(L) = IND(M)
             IND(M) = S
             IF (M<=K) L = M + 1
             IF (M>=K) U = M - 1
             END DO
          END SUBROUTINE SELECT_ON_COORDINATE

          SUBROUTINE SPREAD_IN_COORDINATE(TP,C,L,U,INTERV)
             TYPE (KDTREE2), POINTER :: TP
             TYPE(INTERVAL), INTENT(OUT) :: INTERV
             INTEGER, INTENT (IN) :: C, L, U
             REAL(KDKIND) :: LAST, LMAX, LMIN, T, SMIN,SMAX
             INTEGER :: I, ULOCAL
             REAL(KDKIND), POINTER :: V(:,:)
             INTEGER, POINTER :: IND(:)
             ! ..
             V => TP%THE_DATA(1:,1:)
             IND => TP%IND(1:)
             SMIN = V(C,IND(L))
             SMAX = SMIN

             ULOCAL = U

             DO I = L + 2, ULOCAL, 2
             LMIN = V(C,IND(I-1))
             LMAX = V(C,IND(I))
             IF (LMIN>LMAX) THEN
                T = LMIN
                LMIN = LMAX
                LMAX = T
             END IF
             IF (SMIN>LMIN) SMIN = LMIN
             IF (SMAX<LMAX) SMAX = LMAX
             END DO
             IF (I==ULOCAL+1) THEN
             LAST = V(C,IND(ULOCAL))
             IF (SMIN>LAST) SMIN = LAST
             IF (SMAX<LAST) SMAX = LAST
             END IF

             INTERV%LOWER = SMIN
             INTERV%UPPER = SMAX

          END SUBROUTINE SPREAD_IN_COORDINATE


          SUBROUTINE KDTREE2_DESTROY(TP)
          TYPE (KDTREE2), POINTER :: TP
          ! ..
          CALL DESTROY_NODE(TP%ROOT)

          DEALLOCATE (TP%IND)
          NULLIFY (TP%IND)

          IF (TP%REARRANGE) THEN
             DEALLOCATE(TP%REARRANGED_DATA)
             NULLIFY(TP%REARRANGED_DATA)
          ENDIF

          DEALLOCATE(TP)
          RETURN

          CONTAINS
          RECURSIVE SUBROUTINE DESTROY_NODE(NP)
             TYPE (TREE_NODE), POINTER :: NP
             INTRINSIC ASSOCIATED
             IF (ASSOCIATED(NP%LEFT)) THEN
             CALL DESTROY_NODE(NP%LEFT)
             NULLIFY (NP%LEFT)
             END IF
             IF (ASSOCIATED(NP%RIGHT)) THEN
             CALL DESTROY_NODE(NP%RIGHT)
             NULLIFY (NP%RIGHT)
             END IF
             IF (ASSOCIATED(NP%BOX)) DEALLOCATE(NP%BOX)
             DEALLOCATE(NP)
             RETURN

          END SUBROUTINE DESTROY_NODE

          END SUBROUTINE KDTREE2_DESTROY

          SUBROUTINE KDTREE2_N_NEAREST(TP,QV,NN,RESULTS)
          TYPE (KDTREE2), POINTER      :: TP
          REAL(KDKIND), TARGET, INTENT (IN)    :: QV(:)
          INTEGER, INTENT (IN)         :: NN
          TYPE(KDTREE2_RESULT), TARGET :: RESULTS(:)


          SR%BALLSIZE = HUGE(1.0)
          SR%QV => QV
          SR%NN = NN
          SR%NFOUND = 0
          SR%CENTERIDX = -1
          SR%CORRELTIME = 0
          SR%OVERFLOW = .FALSE.

          SR%RESULTS => RESULTS

          SR%NALLOC = NN

          SR%IND => TP%IND
          SR%REARRANGE = TP%REARRANGE
          IF (TP%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF
          SR%DIMEN = TP%DIMEN

          CALL VALIDATE_QUERY_STORAGE(NN)
          SR%PQ = PQ_CREATE(RESULTS)

          CALL SEARCH(TP%ROOT)

          IF (TP%SORT) THEN
             CALL KDTREE2_SORT_RESULTS(NN, RESULTS)
          ENDIF
          RETURN
          END SUBROUTINE KDTREE2_N_NEAREST

          SUBROUTINE KDTREE2_N_NEAREST_AROUND_POINT(TP,IDXIN,CORRELTIME,&
                                                   NN,RESULTS)
          TYPE (KDTREE2), POINTER        :: TP
          INTEGER, INTENT (IN)           :: IDXIN, CORRELTIME, NN
          TYPE(KDTREE2_RESULT), TARGET   :: RESULTS(:)

          ALLOCATE (SR%QV(TP%DIMEN))
          SR%QV = TP%THE_DATA(:,IDXIN)
          SR%BALLSIZE = HUGE(1.0)
          SR%CENTERIDX = IDXIN
          SR%CORRELTIME = CORRELTIME

          SR%NN = NN
          SR%NFOUND = 0

          SR%DIMEN = TP%DIMEN
          SR%NALLOC = NN

          SR%RESULTS => RESULTS

          SR%IND => TP%IND
          SR%REARRANGE = TP%REARRANGE

          IF (SR%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF

          CALL VALIDATE_QUERY_STORAGE(NN)
          SR%PQ = PQ_CREATE(RESULTS)

          CALL SEARCH(TP%ROOT)

          IF (TP%SORT) THEN
             CALL KDTREE2_SORT_RESULTS(NN, RESULTS)
          ENDIF
          DEALLOCATE (SR%QV)
          RETURN
          END SUBROUTINE KDTREE2_N_NEAREST_AROUND_POINT

          SUBROUTINE KDTREE2_R_NEAREST(TP,QV,R2,NFOUND,NALLOC,RESULTS)
          TYPE (KDTREE2), POINTER      :: TP
          REAL(KDKIND), TARGET, INTENT (IN)    :: QV(:)
          REAL(KDKIND), INTENT(IN)             :: R2
          INTEGER, INTENT(OUT)         :: NFOUND
          INTEGER, INTENT (IN)         :: NALLOC
          TYPE(KDTREE2_RESULT), TARGET :: RESULTS(:)

          SR%QV => QV
          SR%BALLSIZE = R2
          SR%NN = 0
          SR%NFOUND = 0
          SR%CENTERIDX = -1
          SR%CORRELTIME = 0

          SR%RESULTS => RESULTS

          CALL VALIDATE_QUERY_STORAGE(NALLOC)
          SR%NALLOC = NALLOC
          SR%OVERFLOW = .FALSE.
          SR%IND => TP%IND
          SR%REARRANGE= TP%REARRANGE

          IF (TP%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF
          SR%DIMEN = TP%DIMEN

          CALL SEARCH(TP%ROOT)
          NFOUND = SR%NFOUND
          IF (TP%SORT) THEN
             CALL KDTREE2_SORT_RESULTS(NFOUND, RESULTS)
          ENDIF

          IF (SR%OVERFLOW) THEN
             WRITE (*,*) 'KD_TREE_TRANS: WARNING! RETURN FROM',&
                        ' KDTREE2_R_NEAREST FOUND MORE NEIGHBORS'
             WRITE (*,*) 'KD_TREE_TRANS: THAN STORAGE WAS PROVIDED FOR.',&
                        '  ANSWER IS NOT SMALLEST BALL'
             WRITE (*,*) 'KD_TREE_TRANS: WITH THAT NUMBER OF NEIGHBORS!',&
                        '  I.E. IT IS WRONG.'
          ENDIF

          RETURN
          END SUBROUTINE KDTREE2_R_NEAREST

          SUBROUTINE KDTREE2_R_NEAREST_AROUND_POINT(TP,IDXIN,CORRELTIME,R2,&
               NFOUND,NALLOC,RESULTS)

          TYPE (KDTREE2), POINTER      :: TP
          INTEGER, INTENT (IN)         :: IDXIN, CORRELTIME, NALLOC
          REAL(KDKIND), INTENT(IN)             :: R2
          INTEGER, INTENT(OUT)         :: NFOUND
          TYPE(KDTREE2_RESULT), TARGET :: RESULTS(:)

          INTRINSIC HUGE
          ! ..
          ALLOCATE (SR%QV(TP%DIMEN))
          SR%QV = TP%THE_DATA(:,IDXIN)
          SR%BALLSIZE = R2
          SR%NN = 0
          SR%NFOUND = 0
          SR%CENTERIDX = IDXIN
          SR%CORRELTIME = CORRELTIME

          SR%RESULTS => RESULTS

          SR%NALLOC = NALLOC
          SR%OVERFLOW = .FALSE.

          CALL VALIDATE_QUERY_STORAGE(NALLOC)

          SR%IND => TP%IND
          SR%REARRANGE = TP%REARRANGE

          IF (TP%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF
          SR%REARRANGE = TP%REARRANGE
          SR%DIMEN = TP%DIMEN

          CALL SEARCH(TP%ROOT)
          NFOUND = SR%NFOUND
          IF (TP%SORT) THEN
             CALL KDTREE2_SORT_RESULTS(NFOUND,RESULTS)
          ENDIF

          IF (SR%OVERFLOW) THEN
             WRITE (*,*) 'KD_TREE_TRANS: WARNING! RETURN FROM',&
                        ' KDTREE2_R_NEAREST FOUND MORE NEIGHBORS'
             WRITE (*,*) 'KD_TREE_TRANS: THAN STORAGE WAS PROVIDED FOR.',&
                        '  ANSWER IS NOT SMALLEST BALL'
             WRITE (*,*) 'KD_TREE_TRANS: WITH THAT NUMBER OF NEIGHBORS!',&
                        '  I.E. IT IS WRONG.'
          ENDIF

          DEALLOCATE (SR%QV)
          RETURN
          END SUBROUTINE KDTREE2_R_NEAREST_AROUND_POINT

          FUNCTION KDTREE2_R_COUNT(TP,QV,R2) RESULT(NFOUND)
          TYPE (KDTREE2), POINTER   :: TP
          REAL(KDKIND), TARGET, INTENT (IN) :: QV(:)
          REAL(KDKIND), INTENT(IN)          :: R2
          INTEGER                   :: NFOUND
          INTRINSIC HUGE
          ! ..
          SR%QV => QV
          SR%BALLSIZE = R2

          SR%NN = 0
          SR%NFOUND = 0
          SR%CENTERIDX = -1
          SR%CORRELTIME = 0

          NULLIFY(SR%RESULTS)

          SR%NALLOC = 0

          SR%IND => TP%IND
          SR%REARRANGE = TP%REARRANGE
          IF (TP%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF
          SR%DIMEN = TP%DIMEN

          SR%OVERFLOW = .FALSE.

          CALL SEARCH(TP%ROOT)

          NFOUND = SR%NFOUND

          RETURN
          END FUNCTION KDTREE2_R_COUNT

          FUNCTION KDTREE2_R_COUNT_AROUND_POINT(TP,IDXIN,CORRELTIME,R2)RESULT(NFOUND)

          TYPE (KDTREE2), POINTER :: TP
          INTEGER, INTENT (IN)    :: CORRELTIME, IDXIN
          REAL(KDKIND), INTENT(IN)        :: R2
          INTEGER                 :: NFOUND

          INTRINSIC HUGE
          ! ..
          ALLOCATE (SR%QV(TP%DIMEN))
          SR%QV = TP%THE_DATA(:,IDXIN)
          SR%BALLSIZE = R2

          SR%NN = 0
          SR%NFOUND = 0
          SR%CENTERIDX = IDXIN
          SR%CORRELTIME = CORRELTIME
          NULLIFY(SR%RESULTS)

          SR%NALLOC = 0

          SR%IND => TP%IND
          SR%REARRANGE = TP%REARRANGE

          IF (SR%REARRANGE) THEN
             SR%DATA => TP%REARRANGED_DATA
          ELSE
             SR%DATA => TP%THE_DATA
          ENDIF
          SR%DIMEN = TP%DIMEN
          SR%OVERFLOW = .FALSE.

          CALL SEARCH(TP%ROOT)

          NFOUND = SR%NFOUND

          RETURN
          END FUNCTION KDTREE2_R_COUNT_AROUND_POINT


          SUBROUTINE VALIDATE_QUERY_STORAGE(N)

          INTEGER, INTENT(IN) :: N

          IF (SIZE(SR%RESULTS,1) .LT. N) THEN
             WRITE (*,*) 'KD_TREE_TRANS:  YOU DID NOT PROVIDE ENOUGH',&
                        ' STORAGE FOR RESULTS(1:N)'
             STOP
             RETURN
          ENDIF

          RETURN
          END SUBROUTINE VALIDATE_QUERY_STORAGE

          FUNCTION SQUARE_DISTANCE(D, IV,QV) RESULT (RES)

          REAL(KDKIND) :: RES
          INTEGER :: D
          REAL(KDKIND) :: IV(:),QV(:)

          RES = SUM( (IV(1:D)-QV(1:D))**2 )
          END FUNCTION SQUARE_DISTANCE

          RECURSIVE SUBROUTINE SEARCH(NODE)

          TYPE (TREE_NODE), POINTER          :: NODE
          ! ..
          TYPE(TREE_NODE),POINTER            :: NCLOSER, NFARTHER
          !
          INTEGER                            :: CUT_DIM, I
          ! ..
          REAL(KDKIND)                               :: QVAL, DIS
          REAL(KDKIND)                               :: BALLSIZE
          REAL(KDKIND), POINTER           :: QV(:)
          TYPE(INTERVAL), POINTER :: BOX(:)

          IF ((ASSOCIATED(NODE%LEFT) .AND. ASSOCIATED(NODE%RIGHT))&
              .EQV. .FALSE.) THEN
             IF (SR%NN .EQ. 0) THEN
             CALL PROCESS_TERMINAL_NODE_FIXEDBALL(NODE)
             ELSE
             CALL PROCESS_TERMINAL_NODE(NODE)
             ENDIF
          ELSE
             QV => SR%QV(1:)
             CUT_DIM = NODE%CUT_DIM
             QVAL = QV(CUT_DIM)

             IF (QVAL < NODE%CUT_VAL) THEN
                 NCLOSER => NODE%LEFT
                 NFARTHER => NODE%RIGHT
                 DIS = (NODE%CUT_VAL_RIGHT - QVAL)**2
             ELSE
                 NCLOSER => NODE%RIGHT
                 NFARTHER => NODE%LEFT
                 DIS = (NODE%CUT_VAL_LEFT - QVAL)**2
             ENDIF

             IF (ASSOCIATED(NCLOSER)) CALL SEARCH(NCLOSER)

             IF (ASSOCIATED(NFARTHER)) THEN
             BALLSIZE = SR%BALLSIZE
             IF (DIS <= BALLSIZE) THEN
                BOX => NODE%BOX(1:)
                DO I=1,SR%DIMEN
                IF (I .NE. CUT_DIM) THEN
                   DIS = DIS + DIS2_FROM_BND(QV(I),BOX(I)%LOWER,&
                        BOX(I)%UPPER)
                   IF (DIS > BALLSIZE) THEN
                        RETURN
                   ENDIF
                ENDIF
                END DO
                CALL SEARCH(NFARTHER)
             ENDIF
             ENDIF
          END IF
          END SUBROUTINE SEARCH


          REAL(KDKIND) FUNCTION DIS2_FROM_BND(X,AMIN,AMAX) RESULT (RES)
          REAL(KDKIND), INTENT(IN) :: X, AMIN,AMAX

          IF (X > AMAX) THEN
             RES = (X-AMAX)**2;
             RETURN
          ELSE
             IF (X < AMIN) THEN
             RES = (AMIN-X)**2;
             RETURN
             ELSE
             RES = 0.0
             RETURN
             ENDIF
          ENDIF
          RETURN
          END FUNCTION DIS2_FROM_BND

          LOGICAL FUNCTION BOX_IN_SEARCH_RANGE(NODE, SR) RESULT(RES)
          TYPE (TREE_NODE), POINTER :: NODE
          TYPE (TREE_SEARCH_RECORD), POINTER :: SR

          INTEGER :: DIMEN, I
          REAL(KDKIND)    :: DIS, BALLSIZE
          REAL(KDKIND)    :: L, U

          DIMEN = SR%DIMEN
          BALLSIZE = SR%BALLSIZE
          DIS = 0.0
          RES = .TRUE.
          DO I=1,DIMEN
             L = NODE%BOX(I)%LOWER
             U = NODE%BOX(I)%UPPER
             DIS = DIS + (DIS2_FROM_BND(SR%QV(I),L,U))
             IF (DIS > BALLSIZE) THEN
             RES = .FALSE.
             RETURN
             ENDIF
          END DO
          RES = .TRUE.
          RETURN
          END FUNCTION BOX_IN_SEARCH_RANGE


          SUBROUTINE PROCESS_TERMINAL_NODE(NODE)

          TYPE (TREE_NODE), POINTER          :: NODE

          REAL(KDKIND), POINTER          :: QV(:)
          INTEGER, POINTER       :: IND(:)
          REAL(KDKIND), POINTER  :: DATA(:,:)

          INTEGER            :: DIMEN, I, INDEXOFI, K, CENTERIDX, CORRELTIME
          REAL(KDKIND)                   :: BALLSIZE, SD, NEWPRI
          LOGICAL                :: REARRANGE
          TYPE(PQ), POINTER      :: PQP

          QV => SR%QV(1:)
          PQP => SR%PQ
          DIMEN = SR%DIMEN
          BALLSIZE = SR%BALLSIZE
          REARRANGE = SR%REARRANGE
          IND => SR%IND(1:)
          DATA => SR%DATA(1:,1:)
          CENTERIDX = SR%CENTERIDX
          CORRELTIME = SR%CORRELTIME


          MAINLOOP: DO I = NODE%L, NODE%U
             IF (REARRANGE) THEN
             SD = 0.0
             DO K = 1,DIMEN
                SD = SD + (DATA(K,I) - QV(K))**2
                IF (SD>BALLSIZE) CYCLE MAINLOOP
             END DO
             INDEXOFI = IND(I)
             ELSE
             INDEXOFI = IND(I)
             SD = 0.0
             DO K = 1,DIMEN
                SD = SD + (DATA(K,INDEXOFI) - QV(K))**2
                IF (SD>BALLSIZE) CYCLE MAINLOOP
             END DO
             ENDIF

             IF (CENTERIDX > 0) THEN
             IF (ABS(INDEXOFI-CENTERIDX) < CORRELTIME) CYCLE MAINLOOP
             ENDIF

             IF (SR%NFOUND .LT. SR%NN) THEN

             SR%NFOUND = SR%NFOUND +1
             NEWPRI = PQ_INSERT(PQP,SD,INDEXOFI)
             IF (SR%NFOUND .EQ. SR%NN) BALLSIZE = NEWPRI

             ELSE
                BALLSIZE = PQ_REPLACE_MAX(PQP,SD,INDEXOFI)
             ENDIF
          END DO MAINLOOP

          SR%BALLSIZE = BALLSIZE

          END SUBROUTINE PROCESS_TERMINAL_NODE

          SUBROUTINE PROCESS_TERMINAL_NODE_FIXEDBALL(NODE)

          TYPE (TREE_NODE), POINTER          :: NODE

          REAL(KDKIND), POINTER          :: QV(:)
          INTEGER, POINTER       :: IND(:)
          REAL(KDKIND), POINTER          :: DATA(:,:)
          INTEGER                :: NFOUND
          INTEGER                :: DIMEN, I, INDEXOFI, K
          INTEGER                :: CENTERIDX, CORRELTIME, NN
          REAL(KDKIND)                   :: BALLSIZE, SD
          LOGICAL                :: REARRANGE

          QV => SR%QV(1:)
          DIMEN = SR%DIMEN
          BALLSIZE = SR%BALLSIZE
          REARRANGE = SR%REARRANGE
          IND => SR%IND(1:)
          DATA => SR%DATA(1:,1:)
          CENTERIDX = SR%CENTERIDX
          CORRELTIME = SR%CORRELTIME
          NN = SR%NN
          NFOUND = SR%NFOUND


          MAINLOOP: DO I = NODE%L, NODE%U

             IF (REARRANGE) THEN
             SD = 0.0
             DO K = 1,DIMEN
                SD = SD + (DATA(K,I) - QV(K))**2
                IF (SD>BALLSIZE) CYCLE MAINLOOP
             END DO
             INDEXOFI = IND(I)
             ELSE
             INDEXOFI = IND(I)
             SD = 0.0
             DO K = 1,DIMEN
                SD = SD + (DATA(K,INDEXOFI) - QV(K))**2
                IF (SD>BALLSIZE) CYCLE MAINLOOP
             END DO
             ENDIF

             IF (CENTERIDX > 0) THEN
             IF (ABS(INDEXOFI-CENTERIDX)<CORRELTIME) CYCLE MAINLOOP
             ENDIF

             NFOUND = NFOUND+1
             IF (NFOUND .GT. SR%NALLOC) THEN
             SR%OVERFLOW = .TRUE.
             ELSE
             SR%RESULTS(NFOUND)%DIS = SD
             SR%RESULTS(NFOUND)%IDX = INDEXOFI
             ENDIF
          END DO MAINLOOP

          SR%NFOUND = NFOUND
          END SUBROUTINE PROCESS_TERMINAL_NODE_FIXEDBALL

          SUBROUTINE KDTREE2_N_NEAREST_BRUTE_FORCE(TP,QV,NN,RESULTS)

          TYPE (KDTREE2), POINTER :: TP
          REAL(KDKIND), INTENT (IN)       :: QV(:)
          INTEGER, INTENT (IN)    :: NN
          TYPE(KDTREE2_RESULT)    :: RESULTS(:)

          INTEGER :: I, J, K
          REAL(KDKIND), ALLOCATABLE :: ALL_DISTANCES(:)
          ! ..
          ALLOCATE (ALL_DISTANCES(TP%N))
          DO I = 1, TP%N
             ALL_DISTANCES(I) = SQUARE_DISTANCE(TP%DIMEN,QV,&
                TP%THE_DATA(:,I))
          END DO

          DO I = 1, NN
             RESULTS(I)%DIS =  HUGE(1.0)
             RESULTS(I)%IDX = -1
          END DO
          DO I = 1, TP%N
             IF (ALL_DISTANCES(I)<RESULTS(NN)%DIS) THEN

             DO J = 1, NN
                IF (ALL_DISTANCES(I)<RESULTS(J)%DIS) EXIT
             END DO

             DO K = NN - 1, J, -1
                RESULTS(K+1) = RESULTS(K)
             END DO
             RESULTS(J)%DIS = ALL_DISTANCES(I)
             RESULTS(J)%IDX = I
             END IF
          END DO
          DEALLOCATE (ALL_DISTANCES)
          END SUBROUTINE KDTREE2_N_NEAREST_BRUTE_FORCE


          SUBROUTINE KDTREE2_R_NEAREST_BRUTE_FORCE(TP,QV,R2,NFOUND,RESULTS)

          TYPE (KDTREE2), POINTER :: TP
          REAL(KDKIND), INTENT (IN)       :: QV(:)
          REAL(KDKIND), INTENT (IN)       :: R2
          INTEGER, INTENT(OUT)    :: NFOUND
          TYPE(KDTREE2_RESULT)    :: RESULTS(:)

          INTEGER :: I, NALLOC
          REAL(KDKIND), ALLOCATABLE :: ALL_DISTANCES(:)
          ! ..
          ALLOCATE (ALL_DISTANCES(TP%N))
          DO I = 1, TP%N
             ALL_DISTANCES(I) = SQUARE_DISTANCE(TP%DIMEN,QV,&
                TP%THE_DATA(:,I))
          END DO

          NFOUND = 0
          NALLOC = SIZE(RESULTS,1)

          DO I = 1, TP%N
             IF (ALL_DISTANCES(I)< R2) THEN

             IF (NFOUND .LT. NALLOC) THEN
                NFOUND = NFOUND+1
                RESULTS(NFOUND)%DIS = ALL_DISTANCES(I)
                RESULTS(NFOUND)%IDX = I
             ENDIF
             END IF
          ENDDO
          DEALLOCATE (ALL_DISTANCES)

          CALL KDTREE2_SORT_RESULTS(NFOUND,RESULTS)


          END SUBROUTINE KDTREE2_R_NEAREST_BRUTE_FORCE

          SUBROUTINE KDTREE2_SORT_RESULTS(NFOUND,RESULTS)

          INTEGER, INTENT(IN)          :: NFOUND
          TYPE(KDTREE2_RESULT), TARGET :: RESULTS(:)

          IF (NFOUND .GT. 1) CALL HEAPSORT_STRUCT(RESULTS,NFOUND)

          RETURN
          END SUBROUTINE KDTREE2_SORT_RESULTS

          SUBROUTINE HEAPSORT(A,IND,N)

          INTEGER,INTENT(IN)          :: N
          REAL(KDKIND), INTENT(INOUT)         :: A(:)
          INTEGER, INTENT(INOUT)      :: IND(:)

          !
          !
          REAL(KDKIND)        :: VALUE
          INTEGER     :: IVALUE

          INTEGER     :: I,J
          INTEGER     :: ILEFT,IRIGHT

          ILEFT=N/2+1
          IRIGHT=N

          IF(N.EQ.1) RETURN

          DO
             IF(ILEFT > 1)THEN
             ILEFT=ILEFT-1
             VALUE=A(ILEFT); IVALUE=IND(ILEFT)
             ELSE
             VALUE=A(IRIGHT); IVALUE=IND(IRIGHT)
             A(IRIGHT)=A(1); IND(IRIGHT)=IND(1)
             IRIGHT=IRIGHT-1
             IF (IRIGHT == 1) THEN
                A(1)=VALUE;IND(1)=IVALUE
                RETURN
             ENDIF
             ENDIF
             I=ILEFT
             J=2*ILEFT
             DO WHILE (J <= IRIGHT)
             IF(J < IRIGHT) THEN
                IF(A(J) < A(J+1)) J=J+1
             ENDIF
             IF(VALUE < A(J)) THEN
                A(I)=A(J); IND(I)=IND(J)
                I=J
                J=J+J
             ELSE
                J=IRIGHT+1
             ENDIF
             END DO
             A(I)=VALUE; IND(I)=IVALUE
          END DO
          END SUBROUTINE HEAPSORT

          SUBROUTINE HEAPSORT_STRUCT(A,N)

          INTEGER,INTENT(IN)                 :: N
          TYPE(KDTREE2_RESULT),INTENT(INOUT) :: A(:)

          !
          !
          TYPE(KDTREE2_RESULT) :: VALUE

          INTEGER     :: I,J
          INTEGER     :: ILEFT,IRIGHT

          ILEFT=N/2+1
          IRIGHT=N

          IF(N.EQ.1) RETURN

          DO
             IF(ILEFT > 1)THEN
             ILEFT=ILEFT-1
             VALUE=A(ILEFT)
             ELSE
             VALUE=A(IRIGHT)
             A(IRIGHT)=A(1)
             IRIGHT=IRIGHT-1
             IF (IRIGHT == 1) THEN
                A(1) = VALUE
                RETURN
             ENDIF
             ENDIF
             I=ILEFT
             J=2*ILEFT
             DO WHILE (J <= IRIGHT)
                 IF(J < IRIGHT) THEN
                    IF(A(J)%DIS < A(J+1)%DIS) J=J+1
                 ENDIF
                 IF(VALUE%DIS < A(J)%DIS) THEN
                    A(I)=A(J);
                    I=J
                    J=J+J
                 ELSE
                    J=IRIGHT+1
                 ENDIF
             END DO
             A(I)=VALUE
          END DO
          END SUBROUTINE HEAPSORT_STRUCT

        END MODULE KDTREE2_MODULE

!----------------------------------------------------------------------------------!
!                                                                                  !
!                          END OF KDTREE ROUTINES                                  !
!                                                                                  !
!----------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------!
!                                                                                  !
!                  START ADCIRC_db_extract GLOBAL ROUTINES                         !
!                                                                                  !
!----------------------------------------------------------------------------------!

        MODULE GLOBAL

            USE KDTREE2_MODULE

            !...SOME GLOBAL VARIABLES
            REAL(8),ALLOCATABLE      :: WEIGHTS(:,:)
            INTEGER                  :: NHARM,NE,NP,NOUT
            REAL(8), PARAMETER       :: PI=3.141592653589793

            CONTAINS 

            SUBROUTINE MakeElementKDTREE(MyNodes,MyElements,MyTree)

            !Modified from Zach's version in GridCompareX
            !Removed Grid Type parts

                IMPLICIT NONE

                INTEGER               :: I
                INTEGER               :: N1,N2,N3
                INTEGER,INTENT(IN)    :: MyElements(:,:)
                REAL(8)               :: X1,X2,X3
                REAL(8)               :: Y1,Y2,Y3
                REAL(8),ALLOCATABLE   :: XYtree(:,:)
                REAL(8),INTENT(IN)    :: MyNodes(:,:)
                TYPE(KDTREE2),POINTER,INTENT(OUT) :: MyTree

                !...Compute Element centers
                ALLOCATE(XYtree(1:2,1:NE))
                DO I = 1,NE
                    N1 = MyElements(I,1)
                    N2 = MyElements(I,2)
                    N3 = MyElements(I,3)
                    X1 = MyNodes(N1,1)
                    X2 = MyNodes(N2,1)
                    X3 = MyNodes(N3,1)
                    Y1 = MyNodes(N1,2)
                    Y2 = MyNodes(N2,2)
                    Y3 = MyNodes(N3,2)
                    XYtree(1,I) = ( X1 + X2 + X3 ) / 3D0
                    XYtree(2,I) = ( Y1 + Y2 + Y3 ) / 3D0
                ENDDO

                MyTree => KDTREE2_CREATE(xytree,SORT=.TRUE.,REARRANGE=.TRUE.)

                RETURN

            END SUBROUTINE


            SUBROUTINE FindElement(X,Y,MyNodes,MyElements,tree,E,W1,W2,W3,ElementFound,USearchDepth)

            !Modified from Zach's version in GridCompareX
            !Removed Grid Type parts

                IMPLICIT NONE

                !...OPTIONAL PARAMETERS
                INTEGER,INTENT(OUT),OPTIONAL     :: E             !...Return element that was found
                INTEGER,INTENT(IN),OPTIONAL      :: USearchDepth  !...Set the maximum search depth
                REAL(8),INTENT(OUT),OPTIONAL     :: W1,W2,W3      !...Return the weights of the 3 nodes
                LOGICAL,INTENT(OUT),OPTIONAL     :: ElementFound  !...Return if the x,y resides in an element

                !...REQUIRED PARAMETERS
                REAL(8),INTENT(IN)               :: X
                REAL(8),INTENT(IN)               :: Y
                REAL(8),INTENT(IN)               :: MyNodes(:,:)
                INTEGER,INTENT(IN)               :: MyElements(:,:)
                TYPE(KDTREE2),POINTER,INTENT(IN) :: TREE

                INTEGER                          :: SearchDepth = 20 !...Maximum number of near elements to check
                INTEGER                          :: MyE

                REAL(8)                          :: X1,X2,X3
                REAL(8)                          :: Y1,Y2,Y3
                REAL(8)                          :: S1,S2,S3
                REAL(8)                          :: TA
                INTEGER                          :: N1,N2,N3
                INTEGER                          :: I
                TYPE(KDTREE2_RESULT),ALLOCATABLE :: KDRESULTS(:)
                TYPE(KDTREE2_RESULT),ALLOCATABLE :: KDRESULTS2(:)
                LOGICAL                          :: SmallSearch
                LOGICAL                          :: Found
                LOGICAL                          :: BADIN
                LOGICAL                          :: EXTENDEDINFOWEIGHT


                !...User specified search depth, otherwise 20
                IF(PRESENT(USEARCHDEPTH))THEN
                    SearchDepth=USEARCHDEPTH
                ENDIF

                !...Sanity Check on depth of search
                IF(Tree%N.LT.SearchDepth)THEN
                    SearchDepth = Tree%N
                ENDIF

                ALLOCATE(KDRESULTS(1:SearchDepth))
                CALL KDTREE2_N_NEAREST(TP=TREE,QV=(/X,Y/),NN=SearchDepth,&
                        RESULTS=KDRESULTS)

                !...SANITY CHECK ON INPUTS
                BADIN=.FALSE.
                IF(PRESENT(W1).AND.&
                    ((.NOT.PRESENT(W2)).OR.&
                     (.NOT.PRESENT(W3)).OR.&
                     (.NOT.PRESENT(E))))THEN
                        BADIN=.TRUE.
                ELSEIF(PRESENT(W2).AND.&
                    ((.NOT.PRESENT(W1)).OR.&
                     (.NOT.PRESENT(W3)).OR.&
                     (.NOT.PRESENT(E))))THEN
                        BADIN=.TRUE.
                ELSEIF(PRESENT(W3).AND.&
                    ((.NOT.PRESENT(W2)).OR.&
                     (.NOT.PRESENT(W1)).OR.&
                     (.NOT.PRESENT(E))))THEN
                        BADIN=.TRUE.
                ENDIF

                IF(PRESENT(W1))THEN
                    EXTENDEDINFOWEIGHT=.TRUE.
                ELSE
                    EXTENDEDINFOWEIGHT=.FALSE.
                ENDIF

                IF(BADIN)THEN
                    WRITE(*,'(A)') "ERROR: Please check input parameters to 'FindElement' subroutine."
                    STOP
                ENDIF

                Found = .FALSE.
                FindEL: DO I = 1,SearchDepth
                    MyE = KDRESULTS(I)%IDX
                    IF(PRESENT(E))E=MyE
                    N1 = MyElements(MyE,1)
                    N2 = MyElements(MyE,2)
                    N3 = MyElements(MyE,3)
                    X1 = MyNodes(N1,1)
                    X2 = MyNodes(N2,1)
                    X3 = MyNodes(N3,1)
                    Y1 = MyNodes(N1,2)
                    Y2 = MyNodes(N2,2)
                    Y3 = MyNodes(N3,2)
                    S1 = ABS((X2*Y3-X3*Y2)-(X*Y3-X3*Y)+(X*Y2-X2*Y))
                    S2 = ABS((X*Y3-X3*Y)-(X1*Y3-X3*Y1)+(X1*Y-X*Y1))
                    S3 = ABS((X2*Y-X*Y2)-(X1*Y-X*Y1)+(X1*Y2-X2*Y1))
                    TA = ABS((X2*Y3-X3*Y2)-(X1*Y3-X3*Y1)+(X1*Y2-X2*Y1))
                    IF((S1+S2+S3).LE.1.001D0*TA)THEN
                        IF(EXTENDEDINFOWEIGHT)THEN
                            W1 = ((X-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y))/TA
                            W2 = ((X-X1)*(Y3-Y1)-(Y-Y1)*(X3-X1))/TA
                            W3 = ((Y-Y1)*(X2-X1)-(X-X1)*(Y2-Y1))/TA
                        ENDIF
                        Found = .TRUE.
                        EXIT FindEL
                    ENDIF
                ENDDO FindEL

                IF(.NOT.Found)THEN
                    IF(PRESENT(E))E  = KDRESULTS(1)%IDX
                    IF(EXTENDEDINFOWEIGHT)THEN
                        W1 = 1D0/3D0
                        W2 = 1D0/3D0
                        W3 = 1D0/3D0
                    ENDIF
                ENDIF

                IF(PRESENT(ELEMENTFOUND))ELEMENTFOUND=FOUND

                RETURN

            END SUBROUTINE


            SUBROUTINE CreateTranslationTableEWeight(MyNodes,MyElements,MyPoints,W)

            !Modified from Zach's version in GridCompareX
            !Removed Grid Type parts
                
                IMPLICIT NONE

                TYPE(KDTREE2),POINTER    :: TREE
                REAL(8),INTENT(IN)       :: MyNodes(:,:)
                REAL(8),INTENT(IN)       :: MyPoints(:,:)
                REAL(8),INTENT(INOUT)    :: W(:,:)
                REAL(8)                  :: X,Y
                REAL(8)                  :: W1,W2,W3
                INTEGER, INTENT(IN)      :: MyElements(:,:)
                INTEGER                  :: I,E
                LOGICAL                  :: ElementFound

                ALLOCATE(WEIGHTS(NP,1:4))

                !...In this sub, find the 10 nearest element centroids
                !   Then figure out which point the node lies in and calculate
                !   interpolation factors.
                WRITE(*,'(A,$)') "Building KDTREE.."
                CALL MakeElementKDTREE(MyNodes,MyElements,Tree)
                WRITE(*,'(A)') "done!"

                WRITE(*,'(A,$)') "Creating translation table..."
                DO I = 1,NOUT
                    X  = MyPoints(I,1)
                    Y  = MyPoints(I,2)
                    CALL FindElement(X,Y,MyNodes,MyElements,tree,E,W1,W2,W3,ElementFound,10)
                    WRITE(40,1999) 'Results for node: ',I,' located at ',X,'  ',Y
                    IF(.NOT.ElementFound) THEN
                       WRITE(40,2000) E,MyElements(E,1),MyElements(E,2),MyElements(E,3)
                    END IF
                    IF(ElementFound) THEN
                       WRITE(40,2100) E,MyElements(E,1),MyElements(E,2),MyElements(E,3)
                    END IF
                    Weights(I,1) = DBLE(E)
                    Weights(I,2) = W1
                    Weights(I,3) = W2
                    Weights(I,4) = W3
                    WRITE(40,2200) W1,W2,W3
                ENDDO    
                WRITE(*,'(A)') "done!"

 1999   FORMAT (A,I6,A,F8.4,A,F8.4)
 2000   FORMAT(/,' WARNING -  SPECIFIED LOCATION DOES NOT LIE', &
                 ' WITHIN ANY ELEMENT IN THE DOMAIN.',/,' CHECK THE', &
                 ' LONGITUDE AND LATITUDE FOR THIS LOCATION',/, &
                 ' PROGRAM WILL USE NEAREST ELEMENT ',I8,/, &
            5X,' WHICH IS MADE UP OF NODES ',3I8)
 2100   FORMAT(/,' SPECIFIED LOCATION WAS FOUND IN ELEMENT ',I8,/, &
            5X,' WHICH IS MADE UP OF NODES ',3I8)
 2200   FORMAT(5X,' INTERPOLATION WEIGHTS ARE  ',F6.4,2x,F6.4,2x,F6.4,/)

                RETURN

            END SUBROUTINE


            SUBROUTINE ComputeHarmonics(MyAmp,MyPhase,MyElements,NewAmp,NewPhase,W,K)
                
                IMPLICIT NONE

                INTRINSIC                :: SQRT,COS,SIN,ACOS

                REAL(8),INTENT(IN)       :: MyAmp(:,:)
                REAL(8),INTENT(IN)       :: MyPhase(:,:)
                REAL(8),INTENT(INOUT)    :: NewAmp(:)
                REAL(8),INTENT(INOUT)    :: NewPhase(:)
                REAL(8),INTENT(IN)       :: W(:,:)
                INTEGER,INTENT(IN)       :: K
                INTEGER,INTENT(IN)       :: MyElements(:,:)
                INTEGER                  :: J

                REAL(8)                  :: C1R,C2R,C3R,CTR
                REAL(8)                  :: C1I,C2I,C3I,CTI
                REAL(8)                  :: W1,W2,W3
                INTEGER                  :: MyE,N1,N2,N3

                DO J=1,NOUT
                   MyE=W(J,1)
                   W1=W(J,2)
                   W2=W(J,3)
                   W3=W(J,4)
                   N1=MyElements(MyE,1)
                   N2=MyElements(MyE,2)
                   N3=MyElements(MyE,3)

                   C1R=MyAmp(N1,K)*COS(MyPhase(N1,K))
                   C1I=MyAmp(N1,K)*SIN(MyPhase(N1,K))
                   C2R=MyAmp(N2,K)*COS(MyPhase(N2,K))
                   C2I=MyAmp(N2,K)*SIN(MyPhase(N2,K))
                   C3R=MyAmp(N3,K)*COS(MyPhase(N3,K))
                   C3I=MyAmp(N3,K)*SIN(MyPhase(N3,K))
                   CTR=C1R*W1+C2R*W2+C3R*W3
                   CTI=C1I*W1+C2I*W2+C3I*W3

                   NewAmp(J)=SQRT(CTR*CTR+CTI*CTI)
                   IF(NewAmp(J).EQ.0.) THEN
                      NewPhase(J)=0.
                   ELSE
                      NewPhase(J)=180.*ACOS(CTR/NewAmp(J))/PI
                      IF(CTI.LT.0.) NewPhase(J)=360.-NewPhase(J)
                   ENDIF
                END DO

                RETURN

            END SUBROUTINE


        END MODULE

!----------------------------------------------------------------------------------!
!                                                                                  !
!                 END OF ADCIRC_db_extract GLOBAL ROUTINES                         !
!                                                                                  !
!----------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------!
!                                                                                  !
!                             START MAIN PROGRAM                                   !
!                                                                                  !
!----------------------------------------------------------------------------------!

        PROGRAM adcirc_db_extract_2012

            USE GLOBAL

            IMPLICIT NONE

            INTEGER                  :: MyFunction
            INTEGER                  :: NODE,I,J,K,JKI,JUNKI
            INTEGER, ALLOCATABLE     :: NM(:,:)                  !dimension NE,3

            REAL(8)                  :: JUNKR
            REAL(8),ALLOCATABLE      :: XY(:,:)                  !dimension NP,2
            REAL(8),ALLOCATABLE      :: XYOUT(:,:)               !dimension NOUT,2
            REAL(8),ALLOCATABLE      :: EAMP(:,:),EPHA(:,:)      !dimension NP,NHARM
            REAL(8),ALLOCATABLE      :: UAMP(:,:),UPHA(:,:)      !dimension NP,NHARM
            REAL(8),ALLOCATABLE      :: VAMP(:,:),VPHA(:,:)      !dimension NP,NHARM
            REAL(8),ALLOCATABLE      :: ETAMP(:),ETPHA(:)        !dimension NOUT
            REAL(8),ALLOCATABLE      :: UTAMP(:),UTPHA(:)        !dimension NOUT
            REAL(8),ALLOCATABLE      :: VTAMP(:),VTPHA(:)        !dimension NOUT

            CHARACTER(10)           :: JUNKC
            CHARACTER(50)           :: db_grid
            CHARACTER(50)           :: extract_points
            CHARACTER(50)           :: elev_harm
            CHARACTER(50)           :: vel_harm
            CHARACTER(10),ALLOCATABLE :: NAMES(:)

            LOGICAL                  :: exists


            WRITE(*,'(A)') '-------------------------------------------'
            WRITE(*,'(A)') "ADCIRC_db_extract_2012 (pre-release)2013/04"
            WRITE(*,'(A)') "    -C. Szpilka"
            WRITE(*,'(A)') ""
            WRITE(*,'(A)') "This program extracts harmonic constants"
            WRITE(*,'(A)') "from any ADCIRC tidal database."
            WRITE(*,'(A)') ""
            WRITE(*,'(A)') "It uses the KDTREE2 algorithm to search for"
            WRITE(*,'(A)') "the element for each extraction point."
            WRITE(*,'(A)') "-------------------------------------------"
            WRITE(*,'(A)') ""
            WRITE(*,'(A)') "What do you want to do?"
            WRITE(*,'(A)') "(1) Extract elevation harmonics only"
            WRITE(*,'(A)') "(2) Extract velocity harmonics only"
            WRITE(*,'(A)') "(3) Extract both elevation and velocity"
            WRITE(*,'(A,$)') "==> "
            READ(*,*) MyFunction

            WRITE(*,'(A)') "Enter file containing extraction points: "
            WRITE(*,'(A)') "  Note: first line lists number of points"
            WRITE(*,'(A)') "  followed by NOUT lines with their",&
                           " locations"
            WRITE(*,'(A)') "  in east-positive long and lat"
            WRITE(*,'(A,$)') "==> "
            READ(*,'(A)') extract_points
            INQUIRE(FILE=TRIM(extract_points),EXIST=exists)
            IF(.NOT.exists)THEN
               WRITE(*,'(3A)') "Specified file ",TRIM(extract_points),&
                               " does not exist."
               STOP
            ENDIF

            WRITE(*,'(A)') "Enter name of grid file used to create",&
                             " tidal database: "
            WRITE(*,'(A,$)') "==> "
            READ(*,'(A)') db_grid
            INQUIRE(FILE=TRIM(db_grid),EXIST=exists)
            IF(.NOT.exists)THEN
               WRITE(*,'(3A)') "Specified file ",TRIM(db_grid),&
                               " does not exist."
               STOP
            ENDIF

            SELECT CASE(MyFunction)
            CASE(1)
                WRITE(*,'(A)') "Enter name of file containing ",&
                                 "global elevation harmonics (fort.53):"
                WRITE(*,'(A,$)') "==> "
                READ(*,'(A)') elev_harm
                INQUIRE(FILE=TRIM(elev_harm),EXIST=exists)
                IF(.NOT.exists)THEN
                   WRITE(*,'(3A)') "Specified file ",TRIM(elev_harm),&
                                   " does not exist."
                   STOP
                ENDIF
                OPEN(40,FILE="tides.dia",ACTION="WRITE")
                OPEN(50,FILE="elev_hc.out",ACTION="WRITE")
            CASE(2)
                WRITE(*,'(A)') "Enter name of file containing ",&
                                 "global velocity harmonics (fort.54):"
                WRITE(*,'(A,$)') "==> "
                READ(*,'(A)') vel_harm
                INQUIRE(FILE=TRIM(vel_harm),EXIST=exists)
                IF(.NOT.exists)THEN
                   WRITE(*,'(3A)') "Specified file ",TRIM(vel_harm),&
                                   " does not exist."
                   STOP
                ENDIF
                OPEN(40,FILE="tides.dia",ACTION="WRITE")
                OPEN(60,FILE="vel_hc.out",ACTION="WRITE")
            CASE(3)
                WRITE(*,'(A)') "Enter name of file containing ",&
                                 "global elevation harmonics (fort.53):"
                WRITE(*,'(A,$)') "==> "
                READ(*,'(A)') elev_harm
                WRITE(*,'(A)') "Enter name of file containing ",&
                                 "global velocity harmonics (fort.54):"
                WRITE(*,'(A,$)') "==> "
                READ(*,'(A)') vel_harm
                INQUIRE(FILE=TRIM(elev_harm),EXIST=exists)
                IF(.NOT.exists)THEN
                   WRITE(*,'(3A)') "Specified file ",TRIM(elev_harm),&
                                   " does not exist."
                   STOP
                ENDIF
                INQUIRE(FILE=TRIM(vel_harm),EXIST=exists)
                IF(.NOT.exists)THEN
                   WRITE(*,'(3A)') "Specified file ",TRIM(vel_harm),&
                                   " does not exist."
                   STOP
                ENDIF
                OPEN(40,FILE="tides.dia",ACTION="WRITE")
                OPEN(50,FILE="elev_hc.out",ACTION="WRITE")
                OPEN(60,FILE="vel_hc.out",ACTION="WRITE") 
            CASE DEFAULT
                WRITE(*,'(A, $)') "Error: you entered ", MyFunction
                WRITE(*,'(A,$)') "You must enter 1 2 or 3. Program ",&
                                 "will terminate." 
                STOP
            END SELECT
            WRITE(*,*) ' '

            !---------------------------------------------------------------
            ! Open extraction point file, allocate arrays and read in points
            !---------------------------------------------------------------
            WRITE(*,'(A,$)') 'Reading in extraction point locations ...'
            OPEN(12,FILE=TRIM(extract_points),ACTION="READ")
            READ(12,*) NOUT
            WRITE(40,*) 'There are ', NOUT, ' points to process.'
            WRITE(40,*) ' '
            ALLOCATE(XYOUT(1:NOUT,1:2))
            DO I=1,NOUT
               READ(12,*) XYOUT(I,1),XYOUT(I,2)
            END DO
            CLOSE(12,STATUS="KEEP")
            WRITE(*,'(A,/)') 'Done!' 

            !----------------------------------------------------------
            ! Open database grid file, allocate arrays and read in data
            !----------------------------------------------------------
            WRITE(*,'(A,$)') 'Reading in database grid...'
            OPEN(14,FILE=TRIM(db_grid),ACTION="READ")
            READ(14,*) JUNKC
            READ(14,*) NE,NP
            ALLOCATE(XY(1:NP,1:2))
            ALLOCATE(NM(1:NE,1:3)) 
            DO I=1,NP
               READ(14,*) JKI,XY(JKI,1),XY(JKI,2),JUNKR
            END DO
            DO J=1,NE
               READ(14,*) JKI,JUNKI,NM(JKI,1),NM(JKI,2),NM(JKI,3)
            END DO
            CLOSE(14,STATUS="KEEP")
            WRITE(*,'(A,/)') 'Done!'

            !---------------------------------------------------------------
            ! Open harmonic database files, allocate arrays and read in data
            !---------------------------------------------------------------
            NHARM=0
            IF(MyFunction .EQ. 1 .OR. MyFunction .EQ. 3) THEN
               WRITE(*,'(A,$)') 'Reading in elevation harmonics...'
               OPEN(53,FILE=TRIM(elev_harm),ACTION="READ")
               READ(53,*) NHARM
               ALLOCATE(EAMP(1:NP,1:NHARM))
               ALLOCATE(EPHA(1:NP,1:NHARM))
               ALLOCATE(ETAMP(1:NOUT))
               ALLOCATE(ETPHA(1:NOUT))
               ALLOCATE(NAMES(1:NHARM))
               DO I=1,NHARM
                  READ(53,*) JUNKR,JUNKR,JUNKR,NAMES(I)
               END DO
               READ(53,*) JUNKI
               IF(JUNKI .NE. NP) THEN
                  WRITE(*,*) 'Number of nodes does not match nodes in',&
                             ' database grid:'
                  WRITE(*,'(I8,A,A)') JunkI, ' nodes in ',&
                                      TRIM(elev_harm)
                  WRITE(*,'(I8,A,A)') NP, ' nodes in ', TRIM(db_grid)
                  WRITE(*,'(A)') 'Program will terminate.'
                  STOP
               END IF
               DO K=1,NP
                  READ(53,*) NODE
                  DO J=1,NHARM
                     READ(53,*)EAMP(NODE,J),EPHA(NODE,J)
                     EPHA(NODE,J)=PI*EPHA(NODE,J)/180.
                  END DO
               END DO
               CLOSE(53,STATUS="KEEP")
               WRITE(*,'(A,/)') 'Done!'
            END IF

            IF(MyFunction .EQ. 2 .OR. MyFunction .EQ. 3) THEN
               WRITE(*,'(A,$)') 'Reading in velocity harmonics...'
               OPEN(54,FILE=TRIM(vel_harm),ACTION="READ")
               IF(NHARM .EQ. 0) THEN
                  READ(54,*) NHARM
               ELSE
                  READ(54,*) JUNKI
               END IF
               ALLOCATE(UAMP(1:NP,1:NHARM))
               ALLOCATE(UPHA(1:NP,1:NHARM))
               ALLOCATE(VAMP(1:NP,1:NHARM))
               ALLOCATE(VPHA(1:NP,1:NHARM))
               ALLOCATE(UTAMP(1:NOUT))
               ALLOCATE(UTPHA(1:NOUT))
               ALLOCATE(VTAMP(1:NOUT))
               ALLOCATE(VTPHA(1:NOUT))
               IF(.NOT.ALLOCATED(NAMES)) THEN
                  ALLOCATE(NAMES(1:NHARM))
                  DO I=1,NHARM
                     READ(54,*) JUNKR,JUNKR,JUNKR,NAMES(I)
                  END DO
               ELSE
                  DO I=1,NHARM
                     READ(54,*) JUNKR,JUNKR,JUNKR,JUNKC
                  END DO
               END IF
               READ(54,*) JUNKI
               IF(JUNKI .NE. NP) THEN
                  WRITE(*,*) 'Number of nodes does not match nodes in',&
                             ' database grid:'
                  WRITE(*,'(I8,A,A)') JunkI, ' nodes in ',&
                                      TRIM(vel_harm)
                  WRITE(*,'(I8,A,A)') NP, ' nodes in ', TRIM(db_grid)
                  WRITE(*,'(A)') 'Program will terminate.'
                  STOP
               END IF
               DO K=1,NP
                  READ(54,*) NODE
                  DO J=1,NHARM
                     READ(54,*)UAMP(NODE,J),UPHA(NODE,J),VAMP(NODE,J),&
                               VPHA(NODE,J)
                     UPHA(NODE,J)=PI*UPHA(NODE,J)/180.
                     VPHA(NODE,J)=PI*VPHA(NODE,J)/180.
                  END DO
               END DO
               CLOSE(54,STATUS="KEEP")
               WRITE(*,'(A,/)') 'Done!'
            END IF

            !-------------------------------------------------------------
            !                   For all NOUT points
            ! Find element and calculate interpolation weights
            ! Compute harmonic constituents and output to UNITS 50 and 60
            !-------------------------------------------------------------

            CALL CreateTranslationTableEWeight(XY,NM,XYOUT,Weights)

            IF(MyFunction .EQ. 1 .OR. MyFunction .EQ. 3) THEN
               DO K=1,NHARM
                  CALL ComputeHarmonics(EAMP,EPHA,NM,ETAMP,ETPHA,&
                       Weights,K)
                  WRITE(50,3998) NAMES(K)
                  DO I=1,NOUT
                     WRITE(50,4002) ETAMP(I),ETPHA(I)
                  END DO
               END DO
            END IF
            IF(MyFunction .EQ. 2 .OR. MyFunction .EQ. 3) THEN
               DO K=1,NHARM
                  CALL ComputeHarmonics(UAMP,UPHA,NM,UTAMP,UTPHA,&
                       Weights,K)
                  CALL ComputeHarmonics(VAMP,VPHA,NM,VTAMP,VTPHA,&
                       Weights,K)
                  WRITE(60,3998) NAMES(K)
                  DO I=1,NOUT
                     WRITE(60,8002) UTAMP(I),UTPHA(I),VTAMP(I),VTPHA(I)
                  END DO
               END DO
            END IF

            !-------------------------------------------------------
            ! Clean up: dealocate arrays, close files, destroy trees?
            !-------------------------------------------------------

            IF(ALLOCATED(NM)) DEALLOCATE(NM)
            IF(ALLOCATED(XY)) DEALLOCATE(XY)
            IF(ALLOCATED(XYOUT)) DEALLOCATE(XYOUT)
            IF(ALLOCATED(EAMP)) DEALLOCATE(EAMP)
            IF(ALLOCATED(EPHA)) DEALLOCATE(EPHA)
            IF(ALLOCATED(UAMP)) DEALLOCATE(UAMP)
            IF(ALLOCATED(UPHA)) DEALLOCATE(UPHA)
            IF(ALLOCATED(VAMP)) DEALLOCATE(VAMP)
            IF(ALLOCATED(VPHA)) DEALLOCATE(VPHA)
            IF(ALLOCATED(ETAMP)) DEALLOCATE(ETAMP)
            IF(ALLOCATED(ETPHA)) DEALLOCATE(ETPHA)
            IF(ALLOCATED(UTAMP)) DEALLOCATE(UTAMP)
            IF(ALLOCATED(UTPHA)) DEALLOCATE(UTPHA)
            IF(ALLOCATED(VTAMP)) DEALLOCATE(VTAMP)
            IF(ALLOCATED(VTPHA)) DEALLOCATE(VTPHA)
            IF(ALLOCATED(NAMES)) DEALLOCATE(NAMES)
            IF(ALLOCATED(Weights)) DEALLOCATE(Weights)

            WRITE(*,4100)
            CLOSE(40,STATUS="KEEP")
            CLOSE(50,STATUS="KEEP")
            CLOSE(60,STATUS="KEEP")
         
            STOP

!---------------------------
! VARIOUS FORMAT STATEMENTS
!---------------------------
 3998 FORMAT(1X,A9)
 4002 FORMAT(2(X,F11.6))
 4100 FORMAT('****',/,' RESULTS HAVE BEEN STORED IN OUTPUT FILES',/, &
                      ' elev_hc.out AND vel_hc.out AS APPROPRIATE',//, &
                      ' DIAGNOSTIC INFORMATION WRITTEN TO tides.dia')
 8002 FORMAT(2(2X,f12.6,2X,F8.3))


      END PROGRAM adcirc_db_extract_2012
