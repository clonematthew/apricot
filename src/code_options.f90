module code_options

   double precision :: xmin, xmax, xminEnd, xmaxEnd
   double precision :: ymin, ymax, yminEnd, ymaxEnd
   double precision :: zmin, zmax, zminEnd, zmaxEnd
   double precision :: image_size
   integer :: npixx, npixy, npixz
   integer :: image_mode, include_sinks, fileType
   character(15) :: image_quantity    
   character(200) :: snap_stem
   integer, allocatable :: snap_number(:)
   integer :: xAxis, yAxis, zAxis
   integer :: sinkID, tracerID

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
      else if (image_mode .gt. 2) then
         allocate(snap_number(2:1))
         read(10, *) snap_number(1), snap_number(2)
      end if

      ! Read file type
      read(10, *) fileType

      ! Read limits, or start and end limits for a movie
      if (image_mode .eq. 1 .or. image_mode .eq. 2)  then
         read(10, *) xmin, xmax
         read(10, *) ymin, ymax
         read(10, *) zmin, zmax
      else if (image_mode .eq. 3) then
         read(10, *) xmin, xmax, xminEnd, xmaxEnd
         read(10, *) ymin, ymax, yminEnd, ymaxEnd
         read(10, *) zmin, zmax, zminEnd, zMaxEnd
      else if (image_mode .eq. 4 .or. image_mode .eq. 5) then
         read(10, *) image_size
      end if

      ! Read the pixels of each dimension
      read(10, *) npixx, npixy, npixz

      ! Read the orientation of the axes
      read(10, *) xAxis, yAxis, zAxis

      ! Read the ID of the sink we want to follow
      if (image_mode .eq. 4) then
         read(10, *) sinkID
      end if

      ! Read the ID of the tracer particle we want to follow
      if (image_mode .eq. 5) then
         read(10, *) tracerID
      end if

      ! Close the options file
      close(10)
   end subroutine read_options
end module
