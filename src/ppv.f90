subroutine makePPV(n,pos,vel,mass,rho,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,ppvcube)

  implicit none

  integer,intent(in) :: n
  double precision,intent(in) :: pos(3,n),vel(3,n),mass(n),rho(n),xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in) :: nx,ny,nz
  double precision,intent(out) :: ppvcube(nx,ny,nz)
  integer :: ix,iy,iz
  double precision :: dx,dy,dz,rholim
  double precision :: posx(n),posy(n),velz(n)
  integer :: i

  dx = (xmax-xmin)/nx
  dy = (ymax-ymin)/ny
  dz = (zmax-zmin)/nz

  ppvcube = 0.

  rholim = 1e4 * 1.4*1.67e-24/1.99e-18 ! convert cm-3 to code units

  posx = pos(1,:)
  posy = pos(2,:)
  velz = -vel(3,:)

  do i=1,n
     if (rho(i) .lt. rholim) cycle
     if ((posx(i) .le. xmin) .or. (posx(i) .ge. xmax) .or. (posy(i) .le. ymin) .or. (posy(i) .ge. ymax) &
          .or. (velz(i) .le. zmin) .or. (velz(i) .ge. zmax)) cycle
     ix = int((posx(i)-xmin)/dx) + 1
     iy = int((posy(i)-ymin)/dy) + 1
     iz = int((velz(i)-zmin)/dz) + 1
     ppvcube(ix,iy,iz) = ppvcube(ix,iy,iz) + mass(i)
  end do

end subroutine makePPV

subroutine writePPV(filename,ppvcube,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz)

  use cellData

  implicit none

  character(len=*),intent(in) :: filename
  double precision,intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in) :: nx,ny,nz
  double precision,intent(in) :: ppvcube(nx,ny,nz)
  double precision :: dx,dy,dz,vcen
  character(len=1) :: newline
  integer :: i,j,k

  newline = ' '

  dx = (xmax-xmin)/nx * udist
  dy = (ymax-ymin)/ny * udist
  dz = (zmax-zmin)/nz ! in code units for now

  open(unit=16,file=trim(filename),status='replace')

  i = 1
  write(16,'(I1)') i
  write(16,'(2(I3,2X))') nx,ny
  write(16,'(I3)') nz
  write(16,'(2(E10.3,2X))') dx,dy

  do i=1,nz
     vcen = zmin + (real(i)-0.5) * dz
     vcen = vcen * udist / utime / 1e5 ! convert to km s-1
     write(16,'(E10.3)') vcen
  end do

  write(16,*) newline

  do k=1,nz
     do j=1,ny
        do i=1,nx
           write(16,'(E10.3)') ppvcube(i,j,k)
        end do
     end do
     write(16,*) newline
  end do

  close(unit=16)

end subroutine writePPV
