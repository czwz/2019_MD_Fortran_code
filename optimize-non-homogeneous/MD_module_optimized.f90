module MD
implicit none
integer,public    ::line,atom
real*8,public     ::sigma,epsilon,cutoff,ae,mass,timestep,kb,eV

type neighbor_list_type
    integer, allocatable:: list(:)
    integer::              length
end type
type(neighbor_list_type), allocatable:: neighbor_list(:)

contains

subroutine IMPORT_CONSTANT
sigma=3.40
epsilon=0.01
kb=1.3806485279
eV=1.60217662
line=8
atom=(line**3)*4
ae=5.315
mass=39.948                                              
timestep=0.001                                           
cutoff=2.5*ae
end subroutine 

subroutine fcc_generator(line,a)
real*8,intent(in):: a
integer,intent(in):: line
integer:: i,j,k
        open(unit=2,file='fcc',status='unknown')
        do i=0,(line-1)
                do j=0,(line-1)
                        do k=0,(line-1)
                                write(2,*) (0+i)*a,(0+j)*a,(0+k)*a
                                write(2,*) (0.5+i)*a,(0.5+j)*a,(0+k)*a
                                write(2,*) (0+i)*a,(0.5+j)*a,(0.5+k)*a
                                write(2,*) (0.5+i)*a,(0+j)*a,(0.5+k)*a
                        enddo 
                enddo
        enddo
        close(2,status='keep')
end subroutine

subroutine new_fcc_generator(line,a)
real*8,intent(in):: a
integer,intent(inout):: line
integer:: i,j,k,half
        half = (line-MOD(line,2))/2
        open(unit=2,file='fcc',status='unknown')
        do i=0,line-1
                do j=0,line-1
                        do k=0,line-1

                                if ((i.GT.half-1))then
                                write(2,*) (0+i)*3*a + 40*a,   (0+j)*3*a + 40*a,      (0+k)*3*a + 40*a
                                write(2,*) (0.5+i)*3*a + 40*a, (0.5+j)*3*a +40*a,     (0+k)*3*a + 40*a
                                write(2,*) (0+i)*3*a + 40*a,   (0.5+j)*3*a +40*a,     (0.5+k)*3*a + 40*a
                                write(2,*) (0.5+i)*3*a + 40*a, (0+j)*3*a +40*a,       (0.5+k)*3*a + 40*a
                                else
                                write(2,*) (0+i)*3*a + 20*a,   (0+j)*3*a + 20*a,      (0+k)*3*a + 20*a
                                write(2,*) (0.5+i)*3*a + 20*a, (0.5+j)*3*a + 20*a,    (0+k)*3*a + 20*a
                                write(2,*) (0+i)*3*a + 20*a,   (0.5+j)*3*a + 20*a,    (0.5+k)*3*a + 20*a
                                write(2,*) (0.5+i)*3*a + 20*a, (0+j)*3*a + 20*a,      (0.5+k)*3*a + 20*a
                                endif

                        enddo
                enddo
        enddo

        line = 20+3*line+20+20

        close(2,status='keep')
end subroutine


subroutine EQ_FCCFILE
real*8,allocatable::   xx(:),yy(:),zz(:)
real*8,allocatable::   rn(:),Un(:)
real*8::               r,a,space,U
integer::              i,j,k,step
space=(1/1000.)
        allocate(xx((3**3)*4),yy((3**3)*4),zz((3**3)*4))
        open(unit=3,file='Un-r',status='unknown')
        step=5000
        do k=1,step
                a=3.+space*k
                cutoff=2.5*a

                CALL fcc_generator(3,a)
                open(unit=2,file='fcc',status='unknown')
                do i=1,(3**3)*4
                        read(2,*) xx(i),yy(i),zz(i)
                enddo
                close(2,status='keep')
                U=0
                do i=1,(3**3)*4
                        do j=(i+1),(3**3)*4
                                r=(pbc((xx(i)-xx(j)),a)**2+pbc((yy(i)-yy(j)),a)**2+pbc((zz(i)-zz(j)),a)**2)**0.5
                                if (r.LE.cutoff) then
                                        U=U+epsilon*((sigma/r)**12-(sigma/r)**6)
                                endif
                        enddo
                enddo
            write(3,*) a,U
        enddo
        close(3,status='keep')
        deallocate(xx,yy,zz)
        
        open(unit=3,file='Un-r',status='unknown')
        allocate(rn(step),Un(step))
        do i=1,step
                read(3,*) rn(i),Un(i)
        enddo
        U=minval(Un)
        do i=1,step
                if (Un(i).EQ.U) then
                        ae=rn(i)
                        CALL new_fcc_generator(line,ae)
                endif
        enddo
        deallocate(rn,Un)
        close(3,status='keep')
end subroutine

subroutine INITIALIZE(T,xc,yc,zc,xp,yp,zp)
real*8,intent(inout):: xc(:),yc(:),zc(:),xp(:),yp(:),zp(:)
real*8::               A,B,C,v1,v2,v3,v_factor,Ek
integer::              i,T                
        OPEN(unit=4,file='fcc',status='unknown')
do i=1,atom
        read(4,*) xc(i),yc(i),zc(i)                      !update current position rc
enddo
        CALL random_seed()
do i=1,atom
        CALL random_number(A)
        CALL random_number(B)
        CALL random_number(C)
        v1=2*A-1
        v2=2*B-1
        v3=2*C-1
        Ek=0.5*mass*(v1**2+v2**2+v3**2)
        v_factor=((3*kb*real(T)*6.02*1000)/(2*Ek))**(0.5)
        xp(i)=-v1*v_factor*(0.01)*timestep+xc(i)
        yp(i)=-v2*v_factor*(0.01)*timestep+yc(i)
        zp(i)=-v3*v_factor*(0.01)*timestep+zc(i)
enddo
        CLOSE(4,status='keep')
end subroutine

subroutine VERLET(xc,yc,zc,xp,yp,zp,fx,fy,fz,v,T,START,FINISH,neighbor_list)
real*8,intent(inout)::  xc(:),yc(:),zc(:),xp(:),yp(:),zp(:),fx(:),fy(:),fz(:),v(:)
real*8,allocatable::    xf(:),yf(:),zf(:)
real*8::                r,factor,vx,vy,vz,v_factor
integer::               i,j,k,l
integer,intent(in)::    T,START,FINISH
type(neighbor_list_type),intent(in):: neighbor_list(:)
factor=eV*6.02*1000                                                         !anstromg/ps2
        allocate(xf(atom),yf(atom),zf(atom))
        CALL force_calculator(neighbor_list,xc,yc,zc,fx,fy,fz,START,FINISH)
        do j=START,FINISH
                xf(j)=2*xc(j)-xp(j)+factor*(fx(j)/mass)*(timestep)**(2)     !update future position rf
                yf(j)=2*yc(j)-yp(j)+factor*(fy(j)/mass)*(timestep)**(2)
                zf(j)=2*zc(j)-zp(j)+factor*(fz(j)/mass)*(timestep)**(2)
                vx=(pbc(xf(j)-xc(j),ae))/(timestep)                         !anstromg/ps
                vy=(pbc(yf(j)-yc(j),ae))/(timestep)
                vz=(pbc(zf(j)-zc(j),ae))/(timestep)
                v(j)=(vx**2+vy**2+vz**2)**0.5
                CALL MIRROR(xf(j),yf(j),zf(j))
                xp(j)=xc(j)
                xc(j)=xf(j)
                yp(j)=yc(j)
                yc(j)=yf(j)
                zp(j)=zc(j)
                zc(j)=zf(j)
        enddo
        deallocate(xf,yf,zf)
end subroutine

subroutine OUTPUT_DATA(x,y,z,v,step)
real*8,intent(in)::    x(:),y(:),z(:),v(:)
real*8::               r,Ek,Untot
integer,intent(in)::   step
integer::              i,j
        Ek=0
        Untot=0
        do i=1,atom
                Ek=Ek+0.5*mass*(v(i)**2)*(10**(-3.)/(6.02*eV))    !eV
                do j=(i+1),atom
                        r=(pbc((x(i)-x(j)),ae)**2+pbc((y(i)-y(j)),ae)**2+pbc((z(i)-z(j)),ae)**2)**0.5
                        if (r.LE.cutoff) then
                                Untot=Untot+epsilon*((sigma/r)**12-(sigma/r)**6)
                        endif
                enddo
        enddo
        100 format(I5," ",f20.4," ",f10.4," ",f20.4," ",f10.4)
        write(2,100) step,Untot,Ek,Untot+Ek,(Ek)*(eV/(kb*atom))*10000*(2/3.)  !eV, eV, eV, K
end subroutine

subroutine MIRROR(xx,yy,zz)
real*8,intent(inout):: xx,yy,zz
        if (xx.GT.ae*line) then
                xx=xx-ae*line
        elseif (xx.LT.0) then
                xx=xx+ae*line
        endif

        if (yy.GT.ae*line) then
                yy=yy-ae*line
        elseif (yy.LT.0) then
                yy=yy+ae*line
        endif

        if (zz.GT.ae*line) then
                zz=zz-ae*line
        elseif (zz.LT.0) then
                zz=zz+ae*line
        endif
end subroutine

subroutine force_calculator(neighbor_list,x,y,z,fx,fy,fz,START,FINISH)
real*8,intent(inout):: x(:),y(:),z(:),fx(:),fy(:),fz(:)
real*8:: r,rx,ry,rz
integer:: i,j,k
integer,intent(in):: START,FINISH
type(neighbor_list_type),intent(in):: neighbor_list(:)
do i=START,FINISH
        fx(i)=0
        fy(i)=0
        fz(i)=0
        do k=1,neighbor_list(i)%length
                j=neighbor_list(i)%list(k)
                        rx=pbc(x(j)-x(i),ae)
                        ry=pbc(y(j)-y(i),ae)
                        rz=pbc(z(j)-z(i),ae)
                        r=(rx**2+ry**2+rz**2)**0.5
                                if (r.LE.cutoff) then
                                        fx(i)=fx(i)-(12*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(rx/r)
                                        fy(i)=fy(i)-(12*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(ry/r)
                                        fz(i)=fz(i)-(12*epsilon/sigma)*((sigma/r)**13-0.5*(sigma/r)**7)*(rz/r)
                                endif
        enddo
enddo
end subroutine

function pbc(r,a)
real*8,intent(in):: r,a
real*8:: pbc,l
l=((line)*a)
        if ((r.GT.0).AND.(r**2.GT.(0.5*l)**2)) then
                pbc=(r-l)
        elseif ((r.LT.0).AND.(r**2.GT.(0.5*l)**2)) then 
                pbc=(r+l)
        elseif ((r**2).LE.(0.5*l)**2) then
                pbc=(r)
        endif
return
end function

subroutine NEIGHBOR(Ntotal,ID,x,y,z,radius,neighbor_list)
real*8,intent(in)::                   x(:),y(:),z(:),radius
real*8::                                 rx, ry, rz, r
type(neighbor_list_type),intent(inout):: neighbor_list(:)
integer,intent(in)::                  ID, Ntotal
integer::                             N_neighbors,i

    N_neighbors=0
    neighbor_list(ID)%length=0
    allocate(neighbor_list(ID)%list(Ntotal))

    do i=1,Ntotal
        if (i.NE.ID) then
            rx=pbc(x(i)-x(ID),ae)
            ry=pbc(y(i)-y(ID),ae)
            rz=pbc(z(i)-z(ID),ae)
            r=(rx**2.+ry**2.+rz**2)**0.5
            if (r.LT.radius) then
                N_neighbors = N_neighbors + 1
                neighbor_list(ID)%list(N_neighbors)=i
            endif
        endif
    enddo
    neighbor_list(ID)%length=N_neighbors

end subroutine NEIGHBOR

subroutine CREATE_NEIGHBOR_LIST(IDstart,IDend,Ntotal,x,y,z,radius,neighbor_list)
real*8,intent(in)::                   x(:),y(:),z(:),radius
type(neighbor_list_type),intent(inout):: neighbor_list(:)
integer,intent(in)::                  IDstart,IDend,Ntotal
integer::                                i
    do i=IDstart,IDend
        CALL NEIGHBOR(Ntotal,i,x,y,z,radius,neighbor_list)
    enddo
end subroutine CREATE_NEIGHBOR_LIST

end module
