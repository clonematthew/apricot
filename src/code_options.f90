module code_options

   double precision :: xmin, xmax
   double precision :: ymin, ymax
   double precision :: zmin, zmax
   integer :: npixx, npixy, npixz
   integer :: image_mode, include_sinks, fileType
   character(15) :: image_quantity    
   character(200) :: snap_stem
   integer, allocatable :: snap_number(:)
   integer :: xAxis, yAxis, zAxis

contains
   subroutine read_options
      ! Open image options file
      open(10, file='image_arepo_options.dat',status='old')
      
      ! Read the image mode
      read(10, *) image_mode

      ! Read the type of image we're making
      read(10, *) image_quantity      

      ! Do we want to keep the sinks?
      read(10, *) include_sinks

      ! Read the snapshot stem
      read(10, *) snap_stem

      ! Read in the image filenumber, or the stop and start filenumbers
      if (image_mode .eq. 1 .or. image_mode .eq. 2) then
         allocate(snap_number(1:1))
         read(10, *) snap_number(1)
      else if (image_mode .eq. 3) then
         allocate(snap_number(2:1))
         read(10, *) snap_number(1), snap_number(2)
      end if

      ! Read file type
      read(10, *) fileType

      ! Read limits
      read(10, *) xmin, xmax
      read(10, *) ymin, ymax
      read(10, *) zmin, zmax

      ! Read the pixels of each dimension
      read(10, *) npixx, npixy, npixz

      ! Read the orientation of the axes
      read(10, *) xAxis, yAxis, zAxis

      ! Close the options file
      close(10)
   end subroutine read_options
end module
