! file: write_spd_fort.f90

! Write stress period data to GHB or DRN files
subroutine write_spd_fort(fname, spd, nrow, ncol, itmp, kper, do_write)

    ! initialize array in subroutine
    character(len=*), intent(in) :: fname
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real, intent(in) :: spd(nrow, ncol)
    integer, intent(in) :: itmp
    integer, intent(in) :: kper
    logical, intent(in) :: do_write
    
    open(42, file = fname, status = 'old', position = 'append', action = 'write')

    ! Line 5
    write(42,1) itmp, 0, kper
  1 format(2I9, ' # stress period ', I9)

    ! Loop through model cells
    if (do_write) then
        do idx = 1, nrow
            
            ! Line 6: Write cell specific data to file
            write(42,3) int(spd(idx, 1:3)), spd(idx, 4:)
          3 format(3I9, 100E15.7)
          
        end do
    end if
    
    close(42)

end

! end file write_spd_fort.f90
