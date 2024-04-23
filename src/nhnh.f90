subroutine make_nhnh(n, nimages, points, density, cell_size, &
                                    x0, x1,  y0, y1,  z0, z1, nx, ny, nz, image)
   !
   ! load the tree module
   use kdtree2_module
   use kdtree2_precision_module
 
   implicit none

   !
   ! Input variables
   integer :: n
   integer :: nimages
   real(kdkind), intent(in), dimension(3, n) ::  points
   real(kdkind), intent(in), dimension(n) :: density
   real(kdkind), intent(in), dimension(n) :: cell_size
   integer, intent(in) :: nx, ny
   double precision, intent(in) :: x0, x1,  y0, y1,  z0, z1

   !
   ! Output variables
   double precision, dimension(nx, ny, nimages) :: image
   double precision, dimension(n) :: coldens

   !
   ! Internal variables
   real(kdkind), dimension(3) :: ray_point
   type(kdtree2), pointer :: tree  ! this is how you declare a tree in your main program
   type(kdtree2_result) :: results(1:3)
   integer :: i, j, nz, ic
   integer :: ip
   integer :: p,xx,yy
   integer :: ix,iy,iz
   double precision :: dx, dy, dz
   double precision, parameter :: d_cell_size = 0.2 ! could make this an input parameter
   integer, dimension(nx, ny) :: nsteps
   character*1 creturn 

   creturn = achar(13)

   !
   ! set pixel size  
 
   dx = (x1 - x0)/real(nx)
   dy = (y1 - y0)/real(ny)

   ! Define x-y plane of image and dimension to integrate through
   ix = 1
   iy = 2
   iz = 3

   !
   ! build the tree
   write (*,*) 'Building Tree...'
   tree => kdtree2_create(points,sort=.false.,rearrange=.false.) 
   write(*, *) 'Done'

   !
   ! set up arrays
   write(*,*) 'Cleaning image arrays..' 
   image = 0
   nsteps = 0
   write(*,*) 'Done'

   !
   ! walk the tree
   write (*,*) 'walking tree to get image'
   do j = 1, ny
      ray_point(iy) = y0 + (j-1)*dy + 0.5*dy
      write(*, 101, advance='no') creturn, int(real(j) / real(ny) * 100)
101   format(a, 'Percent complete: ', i3, ' %')
      do i = 1, nx
         ray_point(ix) = x0 + (i-1)*dx + 0.5*dx 
         !
         ! walk along the ray with adaptive steps 
         nz = 0
         ray_point(iz) = z0
         do while (ray_point(iz) .lt. z1)
            !
            ! query tree
            call kdtree2_n_nearest(tree, ray_point, 1, results) 
            ip = results(1)%idx
            !
            ! set the step size
            dz = cell_size(ip)*d_cell_size
            !
            ! column density (also used for weighting)
            image(i, j, 1) = image(i, j, 1) + dz*density(ip)
            !
            ray_point(iz) = ray_point(iz) + dz 
            ! update couter
            nz = nz + 1
         end do
         nsteps(i, j) = nz
      end do
   end do
   print *, 'Done'

   !
   ! this releases memory for the tree 
   call kdtree2_destroy(tree)

   open(unit=42,file='nhNh.dat',status='replace')

  do p=1,n
     if ((points(ix,p) .le. x0) .or. (points(ix,p) .ge. x1) .or. (points(iy,p) .le. y0) .or. (points(iy,p) .ge. y1) &
          & .or. (points(iz,p) .le. z0) .or. (points(iz,p) .ge. z1)) cycle
     xx = int((points(ix,p)-x0)/dx) + 1
     yy = int((points(iy,p)-y0)/dx) + 1
     coldens(p) = image(xx,yy,1)
     write(42,'(ES10.3,2X,ES10.3)') density(p),coldens(p)
  end do

  close(unit=42)

end subroutine make_nhnh