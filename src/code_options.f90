module codeOptions

   double precision :: xMin, xMax, xMinEnd, xMaxEnd
   double precision :: yMin, yMax, yMinEnd, yMaxEnd
   double precision :: zMin, zMax, zMinEnd, zMaxEnd
   double precision, allocatable :: imageSize(:)
   integer :: nPixX, nPixY, nPixZ
   integer :: imageMode, include_sinks, fileType
   character(15) :: imageQuantity, movieMode   
   character(200) :: snapStem
   integer, allocatable :: snapNumber(:)
   integer :: xAxis, yAxis, zAxis, followID

contains
   subroutine readOptions
      ! Open image options file
      open(10, file='image_arepo_options.dat',status='old')
      
      ! Read the image mode
      read(10, *) imageMode

      ! Read the mode of the image/movie
      if (imageMode .eq. 1) then
         read(10, *) imageQuantity
      else if (imageMode .eq. 2) then
         read(10, *) movieMode
      end if  

      ! Do we want to keep the sinks?
      read(10, *) include_sinks

      ! Read the snapshot stem
      read(10, *) snapStem

      ! Read in the image filenumber, or the stop and start filenumbers
      if (imageMode .eq. 1) then
         allocate(snapNumber(1:1))
         read(10, *) snapNumber(1)
      else if (imageMode .eq. 2) then
         allocate(snapNumber(2:1))
         read(10, *) snapNumber(1), snapNumber(2)
      end if

      ! Read file type
      read(10, *) fileType

      ! Read limits, or start and end limits for a movie
      if (imageMode .eq. 1)  then
         read(10, *) xMin, xMax
         read(10, *) yMin, yMax
         read(10, *) zMin, zMax
      else if (imageMode .eq. 2) then
         ! Linearly move between two points across the movie (or static if same)
         if (movieMode .eq. "linear") then
            read(10, *) xMin, xMax, xMinEnd, xMaxEnd
            read(10, *) yMin, yMax, yMinEnd, yMaxEnd
            read(10, *) zMin, zMax, zMinEnd, zMaxEnd
         ! Follow a sink or tracer particle and show a given radius around it
         else if (movieMode .eq. "followSink" .or. movieMovie .eq. "followTracer") then
            allocate(imageSize(2:1))
            read(10, *) imageSize(1), imageSize(2)
         end if
      end if

      ! Read the pixels of each dimension
      read(10, *) nPixX, nPixY, nPixZ

      ! Read the orientation of the axes
      read(10, *) xAxis, yAxis, zAxis

      ! Read the ID of the sink or tracer we want to follow
      if (imageMode .eq. 2) then
         if (movieMode .eq. "followSink" .or. movieMode .eq. "followTracer") then
            read(10, *) followID
         end if
      end if 

      ! Close the options file
      close(10)
   end subroutine readOptions
end module
