! read UCLCHEM stuff
subroutine read_uclchem_abuns(snapshot_stem, snapshot_number)
  use cell_data
  use kdtree2_module
  use kdtree2_precision_module

  implicit none

  character(*) :: snapshot_stem
  character(len=:), allocatable :: arepo_file_name
  character (len = 3) :: num
  integer :: snapshot_number
  integer :: npost,nchem
  double precision :: xyzcen(3),rmax,r2max,rcell(3),r2cell
  double precision,allocatable :: uclpos(:,:),uclabun(:,:)
  real(kdkind),allocatable :: points(:,:),offset(:,:)
  real(kdkind) :: cellpos(3)
  type(kdtree2), pointer :: tree
  type(kdtree2_result) :: results(1:3)
  integer :: i,j,ip

  write(num, 100) snapshot_number
100 format(I3.3)
  arepo_file_name = trim(snapshot_stem) // '_' // num // '_abundance.dat'

  open(20,file=arepo_file_name,status='old')

  read(20,*) npost,nchem,(xyzcen(i),i=1,3),rmax

  r2max = rmax*rmax

  write(*,*) npost,' chemical tracers, ',nchem,' species'

  arepo%nspecies = nchem

  allocate(uclpos(3,npost))
  allocate(uclabun(nchem,npost))
  allocate(arepo%uclchem(arepo%nspecies,arepo%ngas))

  write(*,*) 'Reading UCLCHEM data'

  do i=1,npost
     read(20,*) (uclpos(j,i),j=1,3),(uclabun(j,i),j=1,nchem)
  end do

  close(unit=20)

  write(*,*) 'Done'

  allocate(points(3,npost),offset(3,npost))

  call random_number(offset)

  offset = minval(arepo%cellsize)*offset*0.1

  points = uclpos + offset

  write(*,*) 'Building tree for chemical tracers'

  tree => kdtree2_create(points,sort=.false.,rearrange=.false.)

  write(*,*) 'Done'

  write(*,*) 'Assigning abundances to gas cells'

  arepo%uclchem = 0.

  do i=1,arepo%ngas
     rcell = arepo%pos(:,i) - xyzcen
     r2cell = sum(rcell*rcell)
     if (r2cell .gt. r2max) cycle ! don't assign abundances for cells outside max tracer radius
     cellpos = arepo%pos(:,i)
     call kdtree2_n_nearest(tree,cellpos,1,results)
     ip = results(1)%idx
     arepo%uclchem(:,i) = uclabun(:,ip)
  end do

  write(*,*) 'Gas cells all assigned'

  deallocate(uclpos,uclabun,points,offset)

  call kdtree2_destroy(tree)
  
end subroutine read_uclchem_abuns

subroutine make_uclchem_images_arepo(n, nimages, nchem, points, density, chem, cell_size, &
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
   integer :: nchem
   real(kdkind), intent(in), dimension(3, n) ::  points
   real(kdkind), intent(in), dimension(n) :: density
   real(kdkind), intent(in), dimension(nchem, n) :: chem
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
   integer :: i, j, nz, ic
   integer :: ip
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
   ix = 3
   iy = 2
   iz = 1

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
            ! abundances for UCLCHEM species
            do ic=1,nchem
               image(i, j, 1+ic) = image(i, j, 1+ic) + chem(ic, ip)*dz*density(ip)
            end do
            ray_point(iz) = ray_point(iz) + dz 
            ! update couter
            nz = nz + 1
         end do
         nsteps(i, j) = nz
      end do
   end do
   print *, 'Done'

   !
   ! few diagostics....
   print *, 'Max nz', maxval(nsteps)
   print *, 'Min nz', minval(nsteps)
   print *, 'mean nz', sum(nsteps) / real(nx) / real(ny)

   !
   ! this releases memory for the tree 
   call kdtree2_destroy(tree)


end subroutine make_uclchem_images_arepo
