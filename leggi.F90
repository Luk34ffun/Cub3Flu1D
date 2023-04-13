subroutine leggi()
use mod_fluid
implicit none

open(unit=1,file=opfile,form='unformatted')
read(1) ux,uy,uz
close(1)

OPEN(UNIT=1,FILE=opfileth,FORM='unformatted')
READ(1) th
CLOSE(1)

WRITE(*,'(1x,A)') ' The I.C. were read from the previous simulation DATA ',opfile(1:32),'and',opfileth
return
endsubroutine leggi
