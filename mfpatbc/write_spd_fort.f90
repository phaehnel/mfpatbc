!! ==============================================================================
!! This is only to test the subroutines independent of Python

!program main
!    implicit none

!    integer :: kper = 3
!    integer :: perspd(3) = (/3,4,6/)
!    integer :: lenper = size(perspd, 1)
!    real :: spd(2,6)
!    integer :: nrow = size(spd, 1)
!    integer :: ncol = size(spd, 2)
!    logical :: verbose = .TRUE. 
!    integer :: mxactb = 10
!    integer :: ighbcb = 0
!    
!    spd(1,1) = 2
!    spd(1,2) = 10
!    spd(1,3) = 89
!    spd(1,4) = 1.23
!    spd(1,5) = 864
!    spd(1,6) = 1025
!    
!    spd(2,1) = 0
!    spd(2,2) = 4
!    spd(2,3) = 182
!    spd(2,4) = 0.34
!    spd(2,5) = 8640
!    spd(2,6) = 1025
!    
!    call write_stress_period_data(kper, perspd, lenper, spd, nrow, ncol, verbose)

!    return

!end

!! ==============================================================================


! file: write_stress_period_data.f90

! ==============================================================================
! Write stress period data to GHB or DRN files
subroutine write_spd_fort(fname, spd, nrow, ncol, itmp, kper, verbose, do_write)

    ! initialize array in subroutine
    character(len=*), intent(in) :: fname
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real, intent(in) :: spd(nrow, ncol)
    integer, intent(in) :: itmp
    integer, intent(in) :: kper
    logical, intent(in) :: verbose
    logical, intent(in) :: do_write
    
    open(42, file = fname, status = 'old', position = 'append', action = 'write')

    ! Line 5
    write(42,1) itmp, 0, kper
  1 format(2I9, ' # stress period ', I9)

    ! Loop through model cells
    if (do_write) then
        do idx = 1, nrow
            
            ! Add 1 to get 1-based indices as required by MODFLOW
    !        spd(idx, 1) = spd(idx, 1) + 1
    !        spd(idx, 2) = spd(idx, 2) + 1
    !        spd(idx, 3) = spd(idx, 3) + 1
            
            ! Line 6: Write cell specific data to file
            write(42,3) int(spd(idx, 1:3)), spd(idx, 4:)
          3 format(3I9, 100E15.7)
          
        end do
    end if
    
    close(42)

end

! end file write_stress_period_data.f90
