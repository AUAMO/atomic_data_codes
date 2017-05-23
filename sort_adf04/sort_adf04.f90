!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     program: sort_adfo4
!     author:  ivan
!
!     why????? We have two adf04 files. file1 is a 'good' adf04 file
!     that was generated using electric dipole (E1) transitions. At a
!     later time it was decided that we would like to include electric 
!     quadrupole (E2), magnetic dipole (M1), and magnetic 
!     quadrupole (M2) transitions. file2 is a 'dummy' adf04 file created
!     by AUTOSTRUCTURE that contains an A-value calculated using 
!     E1/E2/M1/M2transitions. file2 is out of order and contains a much 
!     greater number of configurations than file1. We are not interested 
!     in these extra transitions
!
!     This program reads the portion of each file that 
!     contains the total A-value into character arrays AA(E1 rates) and 
!     BB(E1,E2,M1,M2 rates). It then compares each element of AA to 
!     the elements of BB. If an atomic transition from AA is found 
!     in BB, the rate (in AA) is replaced by its counterpart from BB. 
!     The resulting combined array is then written to the file 'output.dat'
!
!     variables....
!     rdline, c, d      :  used to store transient string values
!     file1, file2      :  the filenames
!     AA, BB            :  character arrays used to store the A-values
!     a1, a2, b1, b2    :  the beginning and ending row indices for AA, BB
!     i, j              :  loop counters       
!
!     The namelist 'input.dat' contains values for file1, file2, a1, a2, b1, b2

      program sort_adf04
!     Declare
      implicit none
      character (len=16) :: rdline, c, d
      Character (len=16) :: file1, file2
      character (len=16), allocatable, dimension(:) :: AA, BB
      integer :: a1, a2, b1, b2, i, j
!     Call the namelist      
      namelist/input/ file1, file2, a1, a2, b1, b2
      open(4, file='input.dat', form='formatted', status='unknown')
      read(4,input)
      close(4)
!     Allocate the arrays
      allocate(AA(a2))
      allocate(BB(b2))
!     Open the in/out files
      open(5, file=file1, form ='formatted')
      open(6, file=file2, form ='formatted')
      open(7, file='output.dat', form='formatted')
!     Read in file1, place into array AA
      do i=1,a2
       read(5,10) rdline
       AA(i) = rdline
      enddo
!     Truncate AA to only include region of interest
      AA = AA(a1:a2)     
!     Uncomment the following section to write AA to file 'fort.30'
!      do i=1,(a2-a1+1)
!        write(30,10) AA(i)
!      end do
!     Read in file 2, place into array BB
      do i=1,b2
       read(6,10) rdline
       BB(i) = rdline
      enddo
!     Truncte BB to only include the region of interest
      BB = BB(b1:b2)
!     Uncomment the following section to write AA to file 'fort.40'
!      do i=1,(b2-b1+1)
!       write(40,11) BB(i)
!      enddo
!     For each element of AA, search BB for a matching transition. If 
!     the transition exists, overwrite AA(i) with BB(j)
      do i = 1,(a2-a1+1)
       c = AA(i)
       c = c(1:8)
        do j = 1,(b2-b1+1)
         d = BB(j)
         d = d(1:8)
          if(c.eq.d) then
            AA(i) = BB(j)
          endif
        enddo
      enddo
!     Write the result to 'output.dat'
      do i=1,(a2-a1+1)
        write(7,10) AA(i)
      end do
!     We're outta here..........
      close(5)
      close(6)
      close(7)

10    format(a16)
  
      end program

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

