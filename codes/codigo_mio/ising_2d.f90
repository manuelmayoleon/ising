program main

    implicit none 
    REAL(kind = 8), EXTERNAL :: dran_u
    integer, external :: i_dran
    integer, dimension(:), allocatable :: s, n1, n2, n3, n4
    real(kind = 8), dimension(-4:4) :: h
    real(kind = 8) :: T  
    integer:: M, M0, mc , contador 
    integer :: L, N , i ,iseed, j ,ij, ib, im 
    real(kind = 8) :: c, rm, rm2, rm1 , rm0 , tau ,error , chi ,uaux, U, U2, uaux2,cv, rm4,cbinder,rm2aux
    real(kind = 8) :: U4, konstaux ,errorchi, errorE,errorcv,errorE2,errorm4,errorrm2
    real(kind = 16) :: errorbinder
    real(kind = 8) :: varm2, varm , Uspin ,varU, varU2
    character(len=10)::long,temp,newn


    L = 164  !Discretization
    N = L*L

    ! M = 1e4 !total measures 
    ! M0 = 1000 !inital thermalization steps
    ! mc = 1e5 !monte carlo steps per measure 


    M = 8192 !total measures 
    M0 = 1000 !inital thermalization steps
    mc = 1 !monte carlo steps per measure 



    tau = 0
    uaux = 0
    uaux2 = 0
    U2 = 0
    U = 0
    cv = 0 
    cbinder = 0
    rm4 = 0 
    iseed = 1210
    rm2aux = 0
    U4 = 0
    konstaux = 0 
    errorE = 0
    errorcv = 0
    errorE2 = 0
    errorm4 = 0
    errorbinder = 0 
    errorrm2 = 0 
    Uspin = 0 


    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ISING_2D_SIMULATION'
    write ( *, '(a)' ) '  FORTRAN90 version'
    write ( *, '(a)' ) '  Monte Carlo simulation of a 2D Ising model.'
  !
  !  Get input from commandline or user.



    call dran_ini(iseed)

    allocate(s(N), n1(N), n2(N), n3(N), n4(N))

    call neighbors(n1,n2,n3,n4,L )



 s = initial_state(N,dran_u)

!!!!!!!! PLOT THE INITIAL STATE
 WRITE(long,'(i0)')   L

! PLOT FILE create .txt to plot with gnuplot
 
 call plot_file ( L, s, 'Initial Configuration', 'ising_2d_initial_L' // trim(adjustl(long)) // '.txt', &
 'ising_2d_initial_L' // trim(adjustl(long))  // '.png' )

    T = 0.49 !initial T 
    do while (T >= 0.49)
        do j = -4, 4, 2 
            h(j) = min(1.0, exp(-2d0*j/T))
        end do 
        !!! GLAUBERT acceptance probability
        ! do j=-4,4,2
        !   h(j)=1.0d/(1.0+exp(2*j/T))
        ! end do
 

        do ij = 1, M0*N !Initial thermalization steps.
            i = i_dran(N)


          
            ib = s(i)* (s(n1(i)) + s(n2(i)) + s(n3(i)) + s(n4(i)))
            
            if (dran_u() < h(ib)) then
                s(i) = -s(i)
            endif

        end do 
        

        !Thermalization has finished

        ! !Initialization of averages
        c = 0.d00 !avg correlation
        rm = 0.d00 !avg magnetization
        rm2 = 0.d00 !avg magnetization squared 
        cv = 0.d00
        rm4 = 0.d00
        U = 0.d00
        U2 = 0.d00 
        U4 = 0.d00
        errorm4 = 0.d00
        

   contador = 0 

        rm1 = real(abs(sum(s)))/N

        ! open(66, file = 'fort.66', status="unknown", action="write") 
        do im = 1, M !M is total measurements
            contador = contador+1
!alternative is to run sequentially through the sites to update instead of selecting them sequentially            
            ! do ij=1,mc
            !   do i=1,N
            !   ib=s(i)*(s(n1(i))+s(n2(i))+s(n3(i))+s(n4(i)))
            !   if (ran_u() < h(ib)) s(i)=-s(i)
            !   enddo
            !   enddo

            do ij = 1, mc*N  !here we evolve mc times in order to avoid correlations
                i = i_dran(N)
                ib = s(i)*(s(n1(i)) + s(n2(i)) + s(n3(i)) + s(n4(i)))
                
                if (dran_u() < h(ib)) then
                    s(i) = -s(i)
                endif
            end do 
        

            !Measurements  
            uaux= 0 
            do i=1,N 
                uaux = uaux  - s(i)* (s(n1(i)) + s(n2(i)) + s(n3(i)) + s(n4(i)))
            enddo 

        
!!! TAKE MEASURES EACH 20 STEPS
        
            if (contador ==20) then 
                contador = 0
                U = U + uaux 
                U2 = U2 + uaux**2
                U4 = U4 + uaux **4
                rm0 = real(abs(sum(s)))/N !current abs(magnetization)
                rm = rm + rm0 
                rm2 = rm2 + rm0*rm0 !sum of magnetization

                rm4 = rm4 + rm0**4 !sum of magnetization^2
                c = c +rm0*rm1
                rm1 = rm0 

            endif
        end do 

        !final avg 

        rm = rm/M  !mean magnetization 
        rm2= rm2/M  !mean magnetization squared		
        rm4 = rm4 /M !mean magnetization^4
        U = U/M !mean Total Energy
        U2 = U2/M !mean Energy^2
        U4 = U4/M !mean Energy^4

        cbinder =  1-rm4/(3*rm2**2) !binder ratio

        varm2 = rm4 - (rm2)**2  !variance of magnetization ^2
        varm = rm2 - rm*rm  !variance of mean magnetization 
        varU = U2 - U*U  ! variance of U 
        varU2 = U4 - (U2)**2 !Variance of U^2

		 
        c = (c/M - rm*rm)/varm  !correlation time
        chi  = varm * N/T !magnetic susceptibility 
        cv = varU/(N*T**2 ) !specific heat 
		
        Uspin = U/ N !Energy per spin 
		
		
        print*, 'c',  c

        if (c /= 1.0) then
            tau = c/(1.0 - c)
        endif 
        konstaux = sqrt((2*tau + 1)/M)

        error = sqrt(varm*(2*tau + 1)/M)
        errorrm2  =sqrt(varm2)*konstaux

        errorE  = sqrt(varU)*konstaux  !error de U 
        errorE2 = sqrt(varU2)*konstaux
        



        errorchi = N/T * ( (errorrm2)**2 + (2*rm)**2 *error**2)**0.5
       
	
        errorcv = (1/(N*T**2d0))*((errorE2)**2d0 + (2d0*U)**2d0 * (errorE)**2d0)**0.5d0

		
        ! write(66,*) T,rm,error,chi,errorchi,Uspin,errorE/N,cv,errorcv,cbinder
        T = T - 0.05 
        print*, 'T', T 
        
          !  Write the final state to a gnuplot file.
          !
        if (T<=0.49) then
              WRITE(temp,'(F10.2)')   T 

            newn= trim(long)// "_"// trim(adjustl(temp)) 
          call plot_file ( L, s, 'Final Configuration', 'ising_2d_final_'// trim(adjustl(newn))  // '.txt', &
          'ising_2d_final_' // trim(adjustl(newn)) // '.png' )
        end if

    
    end do 
    ! close(88)    
    ! close(66)

contains 

function initial_state(N,dran_u) result(s) 
    implicit none  
    REAL(kind = 8), EXTERNAL :: dran_u 
    integer, dimension(N)  :: s
    integer ::i,N

    do i = 1, N 
        if (dran_u() <0.5) then 
            s(i) = +1
        else 
            s(i) = -1
        endif
    end do
    
END FUNCTION 

end program 

!!!!!!!!!!!!!!!!!!!! NEIGBOURS SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Compute the neighbors on the right, top, left, and bottom sites of i

subroutine neighbors(n1, n2, n3, n4, L)
    implicit none 
    integer, intent(in):: L 
    integer, dimension(L*L)  :: n1,n2,n3,n4
    integer :: ix, iy, i, ix1, ix3,iy2,iy4

    do ix =1, L 
        do iy = 1, L
            i = (iy - 1)*L + ix 
            ix1 = ix + 1 

            if (ix1 == (L+1)) then 
                ix1 = 1
            endif
            n1(i)= (iy - 1)*L + ix1

            iy2 = iy + 1 

            if (iy2 == (L+1)) then 
                iy2 = 1
            endif
            n2(i)= (iy2 - 1)*L + ix

            ix3 = ix - 1 

            if (ix3 == (0)) then 
                ix3 = L
            endif
            n3(i)= (iy - 1)*L + ix3

            iy4 = iy - 1 

            if (iy4 == (0)) then 
                iy4 = L
            endif
            n4(i)= (iy4 - 1)*L + ix
        end do 
    end do 
end subroutine



subroutine plot_file (L, s, title, plot_filename, png_filename )

!*****************************************************************************
!
!! PLOT_FILE writes the current configuration to a GNUPLOT plot file.
!
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) C1(M,N), the current state of the system.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
!    Input, character ( len = * ) PLOT_FILENAME, a name for the GNUPLOT
!    command file to be created.
!
!    Input, character ( len = * ) PNG_FILENAME, the name of the PNG graphics
!    file to be created.
!
  implicit none
  integer ( kind = 4 ) L
!   integer ( kind = 4 ) s(N)
  integer, dimension(L*L) :: s
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j,k
  character ( len = * ) plot_filename
  integer ( kind = 4 ) plot_unit
  character ( len = * ) png_filename
  real ( kind = 8 ) ratio
  character ( len = * ) title
  integer ( kind = 4 ) x1
  integer ( kind = 4 ) x2
  integer ( kind = 4 ) y1
  integer ( kind = 4 ) y2

  call get_unit ( plot_unit )
 
  open ( unit = plot_unit, file = plot_filename, status = 'replace' )

  ratio = real ( L, kind = 8 ) / real ( L, kind = 8 )

  write ( plot_unit, '(a)' ) 'set term png'
  write ( plot_unit, '(a)' ) 'set output "' // trim ( png_filename ) // '"'
  write ( plot_unit, '(a,i4,a)' ) 'set xrange [ 0 : ', L , ' ]'
  write ( plot_unit, '(a,i4,a)' ) 'set yrange [ 0 : ', L , ' ]'
  write ( plot_unit, '(a)' ) 'set nokey'
  write ( plot_unit, '(a)' ) 'set title "' // trim ( title ) // '"'
  write ( plot_unit, '(a)' ) 'unset tics'

  write ( plot_unit, '(a,g14.6)' ) 'set size ratio ', ratio
  k=1
  do j = 1, L
    y1 = j - 1
    y2 = j
    do i = 1, L
      x1 = L - i
      x2 = L + 1 - i
      k = (j - 1)*L + i 
      if ( s(k) < 0 ) then
        write ( plot_unit, '(a,i3,a,i3,a,i3,a,i3,a)' ) 'set object rectangle from ', &
          x1, ',', y1, ' to ', x2, ',', y2, ' fc rgb "black"'
      else
        write ( plot_unit, '(a,i3,a,i3,a,i3,a,i3,a)' ) 'set object rectangle from ', &
          x1, ',', y1, ' to ', x2, ',', y2, ' fc rgb "white"'
      end if
    end do
  end do

  write ( plot_unit, '(a)' ) 'plot 1'
  write ( plot_unit, '(a)' ) 'quit'

  close ( unit = plot_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created the gnuplot graphics file "' // trim ( plot_filename ) // '".'

  return
end

subroutine get_unit ( iunit )
  
    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is a value between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is a value between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
      implicit none
    
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ios
      integer ( kind = 4 ) iunit
      logical lopen
    
      iunit = 0
    
      do i = 1, 99
    
        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
    
          inquire ( unit = i, opened = lopen, iostat = ios )
    
          if ( ios == 0 ) then
            if ( .not. lopen ) then
              iunit = i
              return
            end if
          end if
    
        end if
    
      end do
    
      return
    end