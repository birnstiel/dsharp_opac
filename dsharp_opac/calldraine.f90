program call_draine
implicit none

! Arguments:

INTEGER, PARAMETER :: NANG=3
REAL :: X=110433.26495898841
COMPLEX :: REFREL = (1.437, 0.401)
REAL :: GSCA,QBACK,QEXT,QSCA,QABS
COMPLEX :: S1(2*NANG-1),S2(2*NANG-1)

call BHMIE_FORTRAN(X,REFREL,NANG,S1,S2,QEXT,QABS,QSCA,QBACK,GSCA)

write(*,*) ' Qext = ',QEXT
write(*,*) ' QABS = ',QABS
write(*,*) ' QSCA = ',QSCA
write(*,*) ' GSCA = ',GSCA

end program call_draine
