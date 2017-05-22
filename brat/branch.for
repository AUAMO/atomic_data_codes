C FILE: branch_sub.f
      SUBROUTINE csec(istart, ifinish, imeta, ntemps, CSEC_DATA)
      implicit none

      INTEGER, allocatable :: LAT1(:),ISAT1(:),POINTER(:),indx(:)
      REAL, allocatable  :: EMESH1(:),ENAT(:), OMEGAR(:),
     &                        MATA(:,:),br(:),ctemp(:)
      INTEGER :: i,j,k,l,n, ierr, NAST1, MXE1, NOMT1, NZED1, NELC1
      INTEGER :: istart, ifinish, imeta, ntemps
      REAL, dimension(4, ntemps) :: CSEC_DATA
Cf2py intent(in) istart, ifinish, imeta, ntemps
Cf2py intent(out) CSEC_DATA

      CSEC_DATA = 0
      open(10,file='omega',form='formatted',status='unknown')
      open(20,file='bratin',form='formatted',status='unknown')
      read(20,*)
      read(20,*)

      allocate(indx(istart:ifinish))
      allocate(br(istart:ifinish))
      allocate(ctemp(istart:ifinish))

      DO I=istart,ifinish
        read(20,*)indx(i),br(i)
      ENDDO

C This sections reads the first two lines of every
C OMEGA to find out exactly how many points are present
C and to allocate accordingly

      READ(10,*)NZED1,NELC1
      READ(10,*)NAST1,MXE1,NOMT1

      allocate(LAT1(NAST1))         !! stores the L,S or 2J+1 values
      allocate(ISAT1(NAST1))
      allocate(EMESH1(MXE1))        !! stores the energy mesh of size mxe1
      allocate(POINTER(MXE1))
      allocate(ENAT(nast1))    !! stores energy values for the nast1 levels
      allocate(OMEGAR((nast1*(nast1-1))/2))
      allocate(MATA((NAST1*(NAST1-1)) !! stores A-values for all the
     &         /2,mxe1),stat=ierr)     !! transitions
      if(ierr.ne.0)stop ' allocation fails for MATA'

      READ(10,*)(ISAT1(I),LAT1(I),I=1,NAST1)
      READ(10,*)(ENAT(I),I=1,NAST1)

C Circulating round the energies for omega file
380   FORMAT(1pe14.8,6(1pe11.3)/(14x,6(e11.3)))
      K=0
      DO L=1,MXE1
        K=K+1
        READ(10,380)EMESH1(K),(MATA(N,K),N=1,NOMT1)
      ENDDO

C Create an array of energies in eV
      DO i=1,ntemps
        emesh1(i) = emesh1(i) - enat(imeta)
        CSEC_DATA(1,i) = emesh1(i) * 13.6058
      ENDDO
      close(10)

C Calculate the apparent cross section
      ctemp=0.0
      DO i=1,ntemps
       DO j=istart,ifinish
         ctemp(j) = MATA(j,i) * br(j) * 87.974 /
     &  (abs(float(LAT1(imeta))) *  emesh1(i))
       ENDDO
       CSEC_DATA(2,i) = sum(ctemp)
      ENDDO

C Calculate the direct cross section
      ctemp=0.0
      DO i=1,ntemps
       DO j=istart,istart
         ctemp(j) = MATA(j,i) * br(j) * 87.974 /
     &  (abs(float(LAT1(imeta))) *  emesh1(i))
       ENDDO
       CSEC_DATA(3,i) = sum(ctemp)
      ENDDO

C Calculate the cascade cross section
      ctemp=0.0
      DO i=1,ntemps
       DO j=istart+1,ifinish
         ctemp(j) = MATA(j,i) * br(j) * 87.974 /
     &  (abs(float(LAT1(imeta))) *  emesh1(i))
       ENDDO
       CSEC_DATA(4,i) = sum(ctemp)
      ENDDO

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c     SUBROUTINE cmax(istart, ifinish, imeta, ntemps, CSEC_MAX)
c      implicit none
c
c      INTEGER, allocatable :: LAT1(:),ISAT1(:),POINTER(:),indx(:)
c      REAL, allocatable  :: EMESH1(:),ENAT(:), OMEGAR(:),
c     &                        MATA(:,:),br(:)
c      INTEGER :: i,j,k,l,n, ierr, NAST1, MXE1, NOMT1, NZED1, NELC1
c      INTEGER :: istart, ifinish, imeta, ntemps
c      REAL, dimension(istart:ifinish) :: CSEC_MAX
cCf2py intent(in) istart, ifinish, imeta, ntemps
cCf2py intent(out) CSEC_MAX
c
c      CSEC_MAX = 0
c      open(10,file='omega',form='formatted',status='unknown')
c      open(20,file='bratin',form='formatted',status='unknown')
c      read(20,*)
c      read(20,*)
c
c      allocate(indx(istart:ifinish))
c      allocate(br(istart:ifinish))
c
c      DO I=istart,ifinish
c        read(20,*)indx(i),br(i)
c      ENDDO
c
cC This sections reads the first two lines of every
cC OMEGA to find out exactly how many points are present
cC and to allocate accordingly
c
c      READ(10,*)NZED1,NELC1
c      READ(10,*)NAST1,MXE1,NOMT1
c
c      allocate(LAT1(NAST1))         !! stores the L,S or 2J+1 values
c      allocate(ISAT1(NAST1))
c      allocate(EMESH1(MXE1))        !! stores the energy mesh of size mxe1
c      allocate(POINTER(MXE1))
c      allocate(ENAT(nast1))    !! stores energy values for the nast1 levels
c      allocate(OMEGAR((nast1*(nast1-1))/2))
c      allocate(MATA((NAST1*(NAST1-1)) !! stores A-values for all the
c     &         /2,mxe1),stat=ierr)     !! transitions
c      if(ierr.ne.0)stop ' allocation fails for MATA'
c
c      READ(10,*)(ISAT1(I),LAT1(I),I=1,NAST1)
c      READ(10,*)(ENAT(I),I=1,NAST1)
c
cC Circulating round the energies for omega file
c381   FORMAT(1pe14.8,6(1pe11.3)/(14x,6(e11.3)))
c      K=0
c      DO L=1,MXE1
c        K=K+1
c        READ(10,381)EMESH1(K),(MATA(N,K),N=1,NOMT1)
c      ENDDO
c
cC Create an aray of energies in eV
c      DO i=1,ntemps
c        emesh1(i) = emesh1(i) - enat(imeta)
c        CSEC_DATA(1,i) = emesh1(i) * 13.6058
c      ENDDO
c      close(10)
c
cC Calculate the apparent cross section
c       DO j=istart,ifinish
c         CSEC_MAX(j) = MATA(j,:) * br(j) * 87.974 /
c     &  (abs(float(LAT1(imeta))) *  emesh1(:))
c       ENDDO
c
c      END
c
cC END FILE: branch_sub.f
