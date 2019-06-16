program fcc_mov
use MD
implicit none
include 'mpif.h'
character(len=9)::      FILE1,FILE2
integer::               T,i,j,k,iter,START,FINISH,size,rank,ierror,status(MPI_STATUS_SIZE)
real*8,allocatable::    xc(:),yc(:),zc(:),xp(:),yp(:),zp(:),v(:),xcc(:),ycc(:),zcc(:),xpp(:),ypp(:),zpp(:),vc(:)
real*8,allocatable::    fx(:),fy(:),fz(:)
real::                  t1,t2

!!!!BUG-FREE: USE INTEL COMPILER!!!!
CALL MPI_INIT(ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

        T=20       !K
        iter=1000  !fs

        CALL IMPORT_CONSTANT
        allocate(xc(atom),yc(atom),zc(atom),xp(atom),yp(atom),zp(atom),v(atom))
        allocate(fx(atom),fy(atom),fz(atom),neighbor_list(atom))

if (rank.EQ.0) then
        CALL EQ_FCCFILE       !!!INCLUDE THIS LINE IF YOU DONT HAVE THE INPUT FCC FILE!!!  
        CALL INITIALIZE(T,xc,yc,zc,xp,yp,zp)
        write(FILE1,'("Traj_",i2,"K")') T
        write(FILE2,'("Engy_",i2,"K")') T
        open(unit=1,file=FILE1,status='unknown')
        open(unit=2,file=FILE2,status='unknown')
endif
CALL MPI_BCAST(atom,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(xc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(yc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(zc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(xp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(yp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(zp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BCAST(ae,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

        allocate(xcc(atom/size),ycc(atom/size),zcc(atom/size), &
                 xpp(atom/size),ypp(atom/size),zpp(atom/size), &
                 vc(atom/size))

START=1+(rank)*(atom/size)
FINISH=(rank+1)*(atom/size)

        CALL CREATE_NEIGHBOR_LIST(START,FINISH,atom,xc,yc,zc,cutoff,neighbor_list)
t1= MPI_WTIME()
        do j=1,iter
                 
                if (MOD(iter,10).EQ.0) then
                        do k=START,FINISH
                                deallocate(neighbor_list(k)%list)
                        enddo
                        CALL CREATE_NEIGHBOR_LIST(START,FINISH,atom,xc,yc,zc,cutoff,neighbor_list) 
                endif

                CALL VERLET(xc,yc,zc,xp,yp,zp,fx,fy,fz,v,T,START,FINISH,neighbor_list)
                do i=START,FINISH
                        xcc(i-(atom/size)*rank)=xc(i)
                        ycc(i-(atom/size)*rank)=yc(i)
                        zcc(i-(atom/size)*rank)=zc(i)
                        xpp(i-(atom/size)*rank)=xp(i)
                        ypp(i-(atom/size)*rank)=yp(i)
                        zpp(i-(atom/size)*rank)=zp(i)
                        vc(i-(atom/size)* rank)=v(i)
                enddo

                CALL MPI_GATHER(xcc, atom/size, MPI_DOUBLE_PRECISION, xc,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_GATHER(ycc, atom/size, MPI_DOUBLE_PRECISION, yc,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_GATHER(zcc, atom/size, MPI_DOUBLE_PRECISION, zc,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_GATHER(xpp, atom/size, MPI_DOUBLE_PRECISION, xp,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)                
                CALL MPI_GATHER(ypp, atom/size, MPI_DOUBLE_PRECISION, yp,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_GATHER(zpp, atom/size, MPI_DOUBLE_PRECISION, zp,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_GATHER(vc , atom/size, MPI_DOUBLE_PRECISION, v ,atom/size,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

                CALL MPI_BCAST(xc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(yc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(zc, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(xp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(yp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(zp, atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BCAST(v , atom, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
                CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

                if ((rank.EQ.0)) then
                        write(1,*) atom
                        write(1,*) "  "
                        do i=1,atom
                                write(1,*) "Ar",xc(j),yc(j),zc(j)
                        enddo
                        CALL OUTPUT_DATA(xc,yc,zc,v,j)
                endif

        enddo
t2 = MPI_WTIME()
        if (rank.EQ.0) then
                close(1,status='keep')
                close(2,status='keep')
                write(*,*) "JOB FINISH WITH",size,"PROCESSORS."
                write(*,*) "ITERATION TOOK",t2-t1,"SECOND."
        endif

CALL MPI_FINALIZE(ierror)
stop
end program
