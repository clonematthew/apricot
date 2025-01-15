!
! 06/03/2023 PCC
!
! This is a subroutine for making column density images of the sgchem arepo data.
! It's a simple, brute force method, that steps along 'rays' in z to make x, y plane image.
! At each position along the ray it finds the nearest cell point using a tree, and uses the properties (in this case, rho) to build the sum.
! The steps are adaptive, using a fraction of the local cell size (hence why it's included). 
! XXX: Note that the error in the image isn't controlled: it will be improved by reducing the step size.
!
! I/O layout 
! IN: 
!       n		- number of points in data set
!       nimages         - the number of images in the 'image' cube
! 	     points 		- (3, n) array
! 	     density 	- (n) array
!       chem            - (9, n) array
!       tdust           - (n) array
!       tgas            - (n) array
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
!	image		- (nx, ny, nimages) array holding the final image cube
!
subroutine makeSGChemImagesArepo(n, nimages, points, density, chem, tdust, tgas, cell_size, &
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
   real(kdkind), intent(in), dimension(9, n) :: chem
   real(kdkind), intent(in), dimension(n) :: tdust
   real(kdkind), intent(in), dimension(n) :: tgas
   real(kdkind), intent(in), dimension(n) :: cell_size
   integer, intent(in) :: nx, ny
   double precision, intent(in) :: x0, x1,  y0, y1,  z0, z1

   !
   ! Output variables
   double precision, dimension(nx, ny, nimages) :: image

   !
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

   creturn = achar(13)

   !
   ! set pixel size  
 
   dx = (x1 - x0)/real(nx)
   dy = (y1 - y0)/real(ny)

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
      ray_point(2) = y0 + (j-1)*dy + 0.5*dy
      write(*, 101, advance='no') creturn, int(real(j) / real(ny) * 100)
101   format(a, 'Percent complete: ', i3, ' %')
      do i = 1, nx
         ray_point(1) = x0 + (i-1)*dx + 0.5*dx 
         !
         ! walk along the ray with adaptive steps 
         nz = 0
         ray_point(3) = z0
         do while (ray_point(3) .lt. z1)
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
            ! column-weighted gas temp 
            image(i, j, 2) = image(i, j, 2) + tgas(ip)*dz*density(ip) 
            !
            ! column-weighted dust temp
            image(i, j, 3) = image(i, j, 3) + tdust(ip)*dz*density(ip)
            !
            ! column-weighted H2 *abundance*
            image(i, j, 4) = image(i, j, 4) + chem(1, ip)*dz*density(ip)
            !
            ! column-weighted H+ *abundance*
            image(i, j, 5) = image(i, j, 5) + chem(2, ip)*dz*density(ip)
            !
            ! column-weighted CO *abundance*
            image(i, j, 6) = image(i, j, 6) + chem(5, ip)*dz*density(ip)            
            !
            ! column-weighted C+ *abundance*
            image(i, j, 7) = image(i, j, 7) + chem(3, ip)*dz*density(ip)            
            !
            ! move to new point
            ray_point(3) = ray_point(3) + dz 
            ! update couter
            nz = nz + 1
         end do
         nsteps(i, j) = nz
      end do
   end do
   print *, 'Done'

   !
   ! divide through by column image to complete the column weighting
   image(:, :, 2) = image(:, :, 2) / image(:, :, 1)   
   image(:, :, 3) = image(:, :, 3) / image(:, :, 1)
   image(:, :, 4) = image(:, :, 4) / image(:, :, 1)
   image(:, :, 5) = image(:, :, 5) / image(:, :, 1)
   image(:, :, 6) = image(:, :, 6) / image(:, :, 1)
   image(:, :, 7) = image(:, :, 7) / image(:, :, 1)

   !
   ! few diagostics....
   print *, 'Max nz', maxval(nsteps)
   print *, 'Min nz', minval(nsteps)
   print *, 'mean nz', sum(nsteps) / real(nx) / real(ny)

   !
   ! this releases memory for the tree 
   call kdtree2_destroy(tree)

end subroutine makeSGChemImagesArepo
