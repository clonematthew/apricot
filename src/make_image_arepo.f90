!
! 18/01/2023 PCC
!
! This is a subroutine for making column density images from arepo data 
! It's a simple, brute force method, that steps along 'rays' in z to make x, y plane image.
! At each position along the ray it finds the nearest cell point using a tree, and uses the properties (in this case, rho) to build the sum.
! The steps are adaptive, using a fraction of the local cell size (hence why it's included). 
! XXX: Note that the error in the image isn't controlled: it will be improved by reducing the step size.
!
! I/O layout 
! IN: 
!       n		- number of points in data set
! 	     points 		- (3, n) array
! 	     density 	- (n) array
!	     cell_size	- (n) array
!       x0		- min x coord of image pixels
!       x1		- max y coord of image pixels
!       y0 		- etc
!  	  y1		-
!       z0		-
!	     z1		-
!       nx 		- num of pixel rows in image
!       ny 		- num of pixel columns in image
! OUT:
!	image		- (nx, ny) array holding the final image
!
subroutine makeImageArepo(n, points, density, cell_size, x0, x1,  y0, y1,  z0, z1, nx, ny, nz, image, xAxis, yAxis, zAxis)
   
   ! Load the tree module
   use kdtree2_module
   use kdtree2_precision_module
 
   implicit none

   ! Input variables
   integer :: n
   real(kdkind), intent(in), dimension(3,n) ::  points
   real(kdkind), intent(in), dimension(n) :: density
   real(kdkind), intent(in), dimension(n) :: cell_size
   integer, intent(in) :: nx, ny
   double precision, intent(in) :: x0, x1,  y0, y1,  z0, z1

   ! Output variables
   double precision, dimension(nx, ny) :: image

   ! Internal variables
   real(kdkind), dimension(3) :: ray_point
   type(kdtree2), pointer :: tree  ! this is how you declare a tree in your main program
   type(kdtree2_result) :: results(1:3)
   integer :: i, j, nz
   integer :: ip
   double precision :: dx, dy, dz
   double precision, parameter :: d_cell_size = 0.2 ! could make this an input parameter
   integer, dimension(nx, ny) :: nsteps
   character*1 creturn 

   integer :: xAxis, yAxis, zAxis
   double precision :: x0i, x1i, y0i, y1i, z0i, z1i

   creturn = achar(13)

   ! Select pixel size using axes we've chosen
   select case(xAxis)
      case(1)
         dx = (x1 - x0)/real(nx)
         x0i = x0
         x1i = x1 
      case(2)
         dx = (y1 - y0)/real(nx)
         x0i = y0
         x1i = y1
      case(3)
         dx = (z1 - z0)/real(nx)
         x0i = z0
         x1i = z1
   end select

   select case(yAxis)
      case(1)
         dy = (x1 - x0)/real(ny)
         y0i = x0
         y1i = x1 
      case(2)
         dy = (y1 - y0)/real(ny)
         y0i = y0
         y1i = y1
      case(3)
         dy = (z0 - z1)/real(ny)
         y0i = z1
         y1i = z0
   end select

   select case(zAxis)
      case(1)
         z0i = x0
         z1i = x1
      case(2)
         z0i = y0
         z1i = y1
      case(3)
         z0i = z0
         z1i = z1
   end select

   ! Build the tree
   write (*,*) "apricot: Building KDTRee."
   tree => kdtree2_create(points,sort=.false.,rearrange=.false.) 
   write(*, *) "apricot: TreeBuild Completed."

   ! Set up arrays
   image = 0
   nsteps = 0

   ! Walk the tree
   write (*,*) "apricot: Walking tree to get image."
   do j = 1, ny
      ray_point(yAxis) = y0i + (j-1)*dy + 0.5*dy
      write(*, 101, advance='no') creturn, int(real(j) / real(ny) * 100)
101   format(a, "apricot: Percent complete: ", i3, " %")
      do i = 1, nx
         ray_point(xAxis) = x0i + (i-1)*dx + 0.5*dx 

         ! Walk along the ray with adaptive steps 
         nz = 0
         ray_point(zAxis) = z0i
         do while (ray_point(zAxis) .lt. z1i)
            
            ! Query tree
            call kdtree2_n_nearest(tree, ray_point, 1, results) 
            ip = results(1)%idx
            
            ! Set the step size
            dz = cell_size(ip)*d_cell_size
            
            !
            image(i, j) = image(i, j) + dz*density(ip) 
            !
            ! Move to new point
            ray_point(zAxis) = ray_point(zAxis) + dz 

            ! Update couter
            nz = nz + 1
         end do
         nsteps(i, j) = nz
      end do
   end do

   !print *, 'image', image(1, 1), image(nx, ny), image(204, 145)
   !print *, 'Max nz', maxval(nsteps)
   !print *, 'Min nz', minval(nsteps)
   !print *, 'mean nz', sum(nsteps) / real(nx) / real(ny)

   ! This releases memory for the tree 
   call kdtree2_destroy(tree)

end subroutine makeImageArepo
