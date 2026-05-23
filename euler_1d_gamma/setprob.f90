subroutine setprob

    implicit none

    character*12 fname
    integer :: iunit
    real(kind=8) :: rhol, vl, pl, gammal
    real(kind=8) :: rhor, vr, pr, gammar
    common /cparam/ rhol, vl, pl, gammal, rhor, vr, pr, gammar

    iunit = 7
    fname = 'setprob.data'
    call opendatafile(iunit, fname)

    read(iunit,*) rhol
    read(iunit,*) rhor
    read(iunit,*) vl
    read(iunit,*) vr
    read(iunit,*) pl
    read(iunit,*) pr
    read(iunit,*) gammal
    read(iunit,*) gammar

    close(iunit)

end subroutine setprob
