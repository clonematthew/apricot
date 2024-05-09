!
! 21/01/2023 -- PCC
!
! Code to make images from Arepo snapshots. 
! Rather than building the mesh, the code just uses a tree to find nearest neighbours. 
! So can be fast but inaccurate, or slower but more accurate, depending on the purpose at hand.
! The code produces image files that need to be read into python / IDL for plotting
! The files contain all the data needed to make a full image (including range, units, npix, and sink locations, etc)
! 
program image_arepo

   use code_options
   use cell_data 

   implicit none

   double precision, allocatable :: image(:, :)
   double precision, allocatable :: image_cube(:, :, :)
   character(len=:), allocatable :: image_file_name
   integer :: nimages, ic, sn
   character(3) :: num
   character(2) :: cnum

   ! Say hello and read options
   write(*,*) "apricot: Reading code options from file"
   call read_options
   print *, 'apricot: Done.'

   ! Do stuff depending on options
   select case(image_mode)
      ! Image mode = 1 is just a single image (integrated, such as column, or slice)
      case(1)
         call fileReader(fileType, snap_stem, snap_number(1))

         ! now call different image makers, depending on what we want in the image
         if (trim(image_quantity) .eq. 'column') then
            !
            ! allocate image
            print *, 'apricot: Allocating image array of size', npixx, 'by ', npixy
            allocate ( image(1:npixx, 1:npixy) )
            print *, 'apricot: Done.'
            !
            ! make the actual image here
            call make_image_arepo(arepo%ngas, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellsize, & 
                                    xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image,  &
                                    xAxis, yAxis, zAxis)
            !
            ! write the image.  Format is e.g. 'column_MYSNAPSHOT_034'
            write(num, 100) snap_number(1)
            image_file_name = trim(image_quantity) // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            ! deallocate the image array
            deallocate(image)
         else if (trim(image_quantity) .eq. 'sgchem') then
            !
            ! allocate image *cube*
            nimages = 7
            print *, 'apricot: Allocating image cube of size', npixx, ' by', npixy, 'by', nimages
            allocate ( image_cube(1:npixx, 1:npixy, 1:nimages) )
            allocate ( image(1:npixx, 1:npixy) )
            print *, 'apricot: Done.'
            !
            ! make all the image in this subroutine...
            call make_sgchem_images_arepo(arepo%ngas, nimages, arepo%pos(:, 1:arepo%ngas), arepo%rho, &
                                          arepo%chem, arepo%tdust, arepo%temp, arepo%cellsize, &
                                    xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image_cube)
            !
            ! write the images...
            write(num, 100) snap_number(1)
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 1)
            image_file_name = trim(image_quantity) // '_' // 'column' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 2)
            image_file_name = trim(image_quantity) // '_' // 'tgas' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 3)
            image_file_name = trim(image_quantity) // '_' // 'tdust' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 4)
            image_file_name = trim(image_quantity) // '_' // 'xH2' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 5)
            image_file_name = trim(image_quantity) // '_' // 'xHP' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 6)
            image_file_name = trim(image_quantity) // '_' // 'xCO' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 7)
            image_file_name = trim(image_quantity) // '_' // 'xCP' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            !
            ! deallocate the image cube
            deallocate(image_cube)
            deallocate(image)
         else if (trim(image_quantity) .eq. 'uclchem') then
            !
            ! allocate image *cube* - all chemical species plus column density
            nimages = arepo%nspecies + 1
            print *, 'apricot: Allocating image cube of size', npixx, ' by', npixy, 'by', nimages
            allocate ( image_cube(1:npixx, 1:npixy, 1:nimages) )
            allocate ( image(1:npixx, 1:npixy) )
            print *, 'apricot: Done.'
            !
            ! make all the image in this subroutine...
            call make_uclchem_images_arepo(arepo%ngas, nimages, arepo%nspecies, arepo%pos(:, 1:arepo%ngas), arepo%rho, &
                                          arepo%uclchem, arepo%cellsize, &
                                    xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image_cube)
            
            ! Write the images...
            write(num, 100) snap_number(1)
            image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 1)
            image_file_name = trim(image_quantity) // '_' // 'column' // '_' // trim(snap_stem) // '_'// num // '.dat'
            call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            
            ! Write UCLCHEM column densities
            do ic=1,arepo%nspecies
               write(cnum,'(I2.2)') ic
               image(1:npixx, 1:npixy) = image_cube(1:npixx, 1:npixy, 1+ic)
               image_file_name = trim(image_quantity) // '_abun' // cnum // '_' // trim(snap_stem) // '_'// num // '.dat'
               call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
            end do
            
            ! Deallocate the image cube
            deallocate(image_cube)
            deallocate(image)           
         else
            print *, "apricot: Image quantity", trim(image_quantity),' is not supported. Feel free to add code!'
            stop
         end if
      ! case 2: PPV cube of gas mass, in z direction by default
      case(2)
         call fileReader(fileType, snap_stem, snap_number(1))

         if (trim(image_quantity) .eq. 'ppv') then
            ! read snapshot data
            call read_arepo_snap_type2(snap_stem,snap_number(1))
            ! make PPV cube
            allocate(image_cube(npixx,npixy,npixz))
            call makeppv(arepo%ngas,arepo%pos(:,1:arepo%ngas),arepo%vel(:,1:arepo%ngas),arepo%mass(1:arepo%ngas), &
                  arepo%rho(1:arepo%ngas), xmin,xmax,ymin,ymax,zmin,zmax,npixx,npixy,npixz,image_cube)
            ! write PPV cube in RADMC format
            write(num, 100) snap_number(1)
            image_file_name = 'image_ppv.out'
            call writeppv(image_file_name,image_cube,xmin,xmax,ymin,ymax,zmin,zmax,npixx,npixy,npixz)
            deallocate(image_cube)
         else if (trim(image_quantity) .eq. 'nhNh') then
            ! read snapshot data
            call read_arepo_snap_type2(snap_stem,snap_number(1))
            nimages = 1
            print *, 'apricot: Allocating image cube of size', npixx, ' by', npixy, 'by', nimages
            allocate ( image_cube(1:npixx, 1:npixy, 1:nimages) )
            allocate ( image(1:npixx, 1:npixy) )
            print *, 'apricot: Done.'

            ! make column density map, assign to cells
            call make_nhnh(arepo%ngas, nimages, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellsize, &
                                    xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image_cube)
            ! deallocate the image cube
            deallocate(image_cube)
            deallocate(image)
         else
            print *, "apricot: Image quantity", trim(image_quantity),' is not supported. Feel free to add code!'
            stop
         end if
      ! Case 3: Make a series of images for a movie
      case(3)
         if (trim(image_quantity) .eq. "column") then
            ! Loop through all the snapshots we want to make an image of
            do sn = snap_number(1), snap_number(2)
               call fileReader(fileType, snap_stem, sn)

               ! Allocate image
               print *, 'apricot: Allocating image array of size', npixx, 'by ', npixy
               allocate ( image(1:npixx, 1:npixy) )
               print *, 'apricot: Done.'
               
               ! Make the actual image here
               call make_image_arepo(arepo%ngas, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellsize, & 
                                       xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image,  &
                                       xAxis, yAxis, zAxis)
               ! Write the image.  Format is e.g. 'column_MYSNAPSHOT_034'
               write(num, 100) sn
               image_file_name = trim(image_quantity) // '_' // trim(snap_stem) // '_'// num // '.dat'
               call write_image(image_file_name, xmin, xmax,  ymin, ymax,  zmin, zmax, npixx, npixy, npixz, image)
               
               ! Deallocate the image array
               deallocate(image)
               
               ! Deallocate the arepo params
               deallocate(arepo%pos)
               deallocate(arepo%vel)
               deallocate(arepo%chem)
               deallocate(arepo%mass)
               deallocate(arepo%u)
               deallocate(arepo%rho)
               deallocate(arepo%temp)
               deallocate(arepo%tdust)
               deallocate(arepo%cellsize)
               deallocate(arepo%ids)
            end do
         end if
      case default
         write(*,*) 'apricot: No rule for image_mode = ', image_mode, '. We better stop!'
         stop
   end select

   ! bye bye...
   write(*, *) "apricot: Done."

   ! any format specifiers
   100       format(I3.3)
end program image_arepo

! Function to read the file based on its type 
subroutine fileReader(fileType, snap_stem, snap_number)
   use cell_data
   implicit none 

   integer :: fileType, snap_number
   character(200) :: snap_stem

   ! Read in the file based on filetype
   select case(fileType)
      ! Type 2 snapshot, with UCLCHEM
      case(1)
         call read_arepo_snap_type2(snap_stem, snap_number)
         call read_uclchem_abuns(snap_stem, snap_number) 
      ! Type 2 snapshot only
      case(2)
         call read_arepo_snap_type2(snap_stem, snap_number)
      ! Type 3, hdf5, snapshot only
      case(3)
         call read_arepo_snap_type3(snap_stem, snap_number)
   end select  
end subroutine fileReader