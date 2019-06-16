program fcc_mov
use MD
implicit none
include 'mpif.h'
character(len=9)::      FILE1,FILE2
integer::               T,i,iter,ierror
real*8,allocatable::    xc(:),yc(:),zc(:),xp(:),yp(:),zp(:),v(:)
real::                  t1,t2

!!!!BUG-FREE: USE INTEL COMPILER!!!!
CALL MPI_INIT(ierror)

        CALL IMPORT_CONSTANT
        CALL EQ_FCCFILE       !!!INCLUDE THIS LINE IF YOU DONT HAVE THE INPUT FCC FILE!!!  
        allocate(xc(atom),yc(atom),zc(atom),xp(atom),yp(atom),zp(atom),v(atom))

        T=20       !K
        iter=1000  !fs

        CALL INITIALIZE(T,xc,yc,zc,xp,yp,zp)
        write(FILE1,'("Traj_",i2,"K")') T
        write(FILE2,'("Engy_",i2,"K")') T
        open(unit=1,file=FILE1,status='unknown')
        open(unit=2,file=FILE2,status='unknown')

t1 = MPI_WTIME()
        do i=1,iter
                CALL VERLET(xc,yc,zc,xp,yp,zp,v,T)
                CALL OUTPUT_DATA(xc,yc,zc,v,i)
        enddo
t2 = MPI_WTIME()
        
        write(*,*) "JOB FINISHED WITH", "1", "PROCESSORS."
        write(*,*) "IT TAKES",t2-t1,"SECOND."
        close(1,status='keep')
        close(2,status='keep')

stop
end program
