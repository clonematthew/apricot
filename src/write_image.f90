!
! 21/01/2023 PCC
!
! This is a subroutine for writing the images to file. 
!
! I/O layout 
! IN: 
!       image_file_name - file to write out to 
!       x0		- min x coord of image pixels
!       x1		- max y coord of image pixels
!       y0 		- etc
!  	  y1		-
!       z0		-
!	     z1		-
!       nx 		- num of pixel rows in image
!       ny 		- num of pixel columns in image
!	     image	- (nx, ny) array holding the image
!
subroutine writeImage(image_file_name, x0, x1,  y0, y1,  z0, z1, nx, ny, nz, image)
   ! Access the arepo data to get time, sinks, units, etc.
   use cellData
   use codeOptions
 
   implicit none

   ! Input variables
   double precision, intent(in) :: x0, x1,  y0, y1,  z0, z1
   integer, intent(in) :: nx, ny, nz
   double precision, intent(in), dimension(nx, ny) :: image

   ! Internal variables
   integer :: i, j
   character(*) :: image_file_name

   write (*,*) "apricot: Writing image to file ", image_file_name

   ! Open file
   open(30, file=image_file_name, access='stream', form='unformatted',status='replace')
   
   ! Write header   
   write(30) arepo%time
   write(30) nx, ny, nz
   write(30) x0, x1,  y0, y1,  z0, z1
   write(30) umass, udist, utime 
   
   ! Write the image
   write(30) image
   
   ! Write out sink info if present
   write(30) arepo%nsink
   if ( arepo%nsink .gt. 0 .and. include_sinks.eq.1 ) then
      write(30) arepo%sinkx
      write(30) arepo%sinky
      write(30) arepo%sinkz
   end if 
   
   ! Close the file
   close(30)
   
end subroutine writeImage
