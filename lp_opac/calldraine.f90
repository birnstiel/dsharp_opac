program call_draine
implicit none

INTEGER,PARAMETER :: MXNANG=1000

! Arguments:

INTEGER :: NANG=3
REAL :: X=110433.26495898841
COMPLEX :: REFREL = (1.437, 0.401)
REAL :: GSCA,QBACK,QEXT,QSCA,QABS
COMPLEX :: S1(2*MXNANG-1),S2(2*MXNANG-1)



call BHMIE(X,REFREL,NANG,S1,S2,QEXT,QABS,QSCA,QBACK,GSCA)

write(*,*) ' Qext = ',QEXT
write(*,*) ' QABS = ',QABS
write(*,*) ' QSCA = ',QSCA
write(*,*) ' GSCA = ',GSCA

end program call_draine
