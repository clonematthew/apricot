!
! 15/01/2025 -- MC
! Code to make images and movies from AREPO snapshots
! 

program imageArepo
    ! Load modules for reading code options and cell data
    use codeOptions
    use cellData

    implicit none

    ! Define variables to use later
    double precision, allocatable :: image(:, :)
    double precision, allocatable :: imageCube(:, :, :)
    character(len=:), allocatable :: imageFilename
    double precision :: dxFrame, dyFrame, dzFrame, dImageSize
    integer :: nImages, ic, sn, si, followIndex
    character(3) :: num
    character(2) :: cnum

    ! Array of chemical species names for SGChem image mode
    character(len=3) :: speciesNames(6)
    speciesNames = ["tGs", "tDu", "xH2", "xHP", "xCO", "xCP"]

    ! Initialise and read options
    write(*,*) "apricot: Reading code options from file."
    call readOptions
    write(*,*) "apricot: Code options read successfully."

    ! Begin the decision tree
    select case(imageMode)
        ! Image mode 1, single image
        case(1)
            ! Load the snapshot
            call fileReader(fileType, snapStem, snapNumber(1), imageMode)

            ! Allocate image array (common to all)
            write(*,*) "apricot: Allocating image array, size ", nPixX, "by", nPixY
            allocate(image(1:nPixX, 1:nPixY))
            allocate(imageCube(1:nPixX, 1:nPixY, 1))

            ! Convert the snap number into a string for writing
            write(num, 100) snapNumber(1)

            ! Generate different types of images based on the requested type
            if (trim(imageQuantity) .eq. "column") then
                ! Call the image maker
                call makeImageArepo(arepo%ngas, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellSize, & 
                                    xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image,  &
                                    xAxis, yAxis, zAxis)

                ! Write out the image
                imageFilename = trim(imageQuantity) // "_" // "column" // "_" // trim(snapStem) // "_" // num // '.dat'
                call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)

            else if (trim(imageQuantity) .eq. "sgchem") then
                ! Allocate image cube
                nImages = 7
                allocate(imageCube(1:nPixX, 1:nPixY, 1:nImages))

                ! Make all the SGChem images
                call makeSGChemImagesArepo(arepo%ngas, nImages, arepo%pos(:, 1:arepo%ngas), arepo%rho, &
                                           arepo%chem, arepo%tdust, arepo%temp, arepo%cellSize, xMin, &
                                           xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, imageCube)

                ! Write out the column density image
                image(1:nPixX, 1:nPixY) = imageCube(1:nPixX, 1:nPixY, 1)
                imageFilename = trim(imageQuantity) // "_" // "column" // "_" // trim(snapStem) // "_" // num // '.dat'
                call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)

                ! Write the abundance images to file
                do ic=1, nImages-1
                    image(1:nPixX, 1:nPixY) = imageCube(1:nPixX, 1:nPixY, 1+ic)
                    imageFilename = trim(imageQuantity) // "_" // speciesNames(ic) // "_" // trim(snapStem) // "_" // num // ".dat"
                    call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)
                end do

            else if (trim(imageQuantity) .eq. "uclchem") then
                ! Allocate image cube
                nImages = arepo%nSpecies + 1
                allocate(imageCube(1:nPixX, 1:nPixY, 1:nImages))

                ! Make all the UCLChem images
                call makeUCLChemImagesArepo(arepo%ngas, nImages, arepo%nSpecies, arepo%pos(:, 1:arepo%ngas), &
                                            arepo%rho, arepo%uclchem, arepo%cellSize, xMin, xMax,  yMin, &
                                            yMax, zMin, zMax, nPixX, nPixY, nPixZ, imageCube)

                ! Write out the column density image
                image(1:nPixX, 1:nPixY) = imageCube(1:nPixX, 1:nPixY, 1)
                imageFilename = trim(imageQuantity) // "_" // "column" // "_" // trim(snapStem) // "_" // num // '.dat'
                call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)

                ! Write the abundance images to file
                do ic=1, arepo%nSpecies
                    image(1:nPixX, 1:nPixY) = imageCube(1:nPixX, 1:nPixY, 1+ic)
                    imageFilename = trim(imageQuantity) // "_" // "_abun" // cnum // "_" // trim(snapStem) // "_" // num // ".dat"
                    call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)
                end do
            
            else if (trim(imageQuantity) .eq. "ppv") then
                ! Allocate image cube
                allocate(imageCube(1:nPixX, 1:nPixY, 1:nPixZ))

                ! Make the PPV cube
                call makePPV(arepo%ngas,arepo%pos(:,1:arepo%ngas),arepo%vel(:,1:arepo%ngas), &
                             arepo%mass(1:arepo%ngas), arepo%rho(1:arepo%ngas), xMin, xMax, &
                             yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, imageCube)

                ! Write out the PPV cube in RADMC format
                imageFilename = "image_ppv.out"
                call writePPV(imageFilename, imageCube, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ)

            else if (trim(imageQuantity) .eq. "nhNh") then
                ! Allocate image cube
                nImages = 1
                allocate(imageCube(1:nPixX, 1:nPixY, 1:nImages))

                ! Make Nhnh Image
                call makeNHNH(arepo%ngas, nImages, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellSize, &
                              xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, imageCube)

            else
                write(*,*) "apricot: Image quantity", trim(imageQuantity), "is not supported."
                stop
            end if

            ! Deallocate arrays
            deallocate(imageCube)
            deallocate(image)

        ! Image mode 2, series of images (i.e a movie)
        case(2)
            ! Work out the dx per frame if using a linear follow
            if (movieMode .eq. "linear") then
                dxFrame = (xMaxEnd - xMax) / (snapNumber(2) - snapNumber(1))
                dyFrame = (yMaxEnd - yMax) / (snapNumber(2) - snapNumber(1))
                dzFrame = (zMaxEnd - zMax) / (snapNumber(2) - snapNumber(1))
            ! Work out how much the radius changes if using a sink/tracer follow
            else if (movieMode .eq. "followSink" .or. movieMode .eq. "followTracer") then
                dImageSize = (imageSize(2) - imageSize(1)) / (snapNumber(2) - snapNumber(1))
            end if

            ! Loop through the snapshot range
            do sn = snapNumber(1), snapNumber(2)
                ! Load this frame's snapshot data
                call fileReader(fileType, snapStem, sn, movieMode)

                ! Find the index of the tracer or sink particle we're following
                if (movieMode .eq. "followSink") then
                    do si = 1, arepo%nsink
                        if (arepo%sinkids(si) .eq. followID) then
                            followIndex = si
                        end if
                    end do
                else if (movieMode .eq. "followTracer") then
                    do si = 1, arepo%nzoomtracer
                        if (arepo%zoomids(si) .eq. followID) then
                            followIndex = si
                        end if
                    end do
                end if   

                ! Convert the snap number into a string for writing
                write(num, 100) sn

                ! Allocate image array
                write(*,*) "apricot: Allocating image array, size ", nPixX, "by", nPixY
                allocate(image(1:nPixX, 1:nPixY))

                ! Update the limits of the image
                if (movieMode .eq. "linear") then
                    xMin = xMin + dxFrame
                    xMax = xMax + dxFrame
                    yMin = yMin + dyFrame
                    yMax = yMax + dyFrame
                    zMin = zMin + dzFrame
                    zMax = zMax + dzFrame
                else if (movieMode .eq. "followSink") then 
                    ! Re-centre on the sink
                    xMin = arepo%sinkpos(1, followIndex) - imageSize(1)
                    xMax = arepo%sinkpos(1, followIndex) + imageSize(1)
                    yMin = arepo%sinkpos(2, followIndex) - imageSize(1)
                    yMax = arepo%sinkpos(2, followIndex) + imageSize(1)
                    zMin = arepo%sinkpos(3, followIndex) - imageSize(1)
                    zMax = arepo%sinkpos(3, followIndex) + imageSize(1)
                    imageSize(1) = imageSize(1) + dImageSize

                    deallocate(arepo%sinkids)
                    deallocate(arepo%sinkpos)
                else if (movieMode .eq. "followTracer") then
                    ! Re-centre on the zoom tracer
                    xMin = arepo%zoompos(1, followIndex) - imageSize(1)
                    xMax = arepo%zoompos(1, followIndex) + imageSize(1)
                    yMin = arepo%zoompos(2, followIndex) - imageSize(1)
                    yMax = arepo%zoompos(2, followIndex) + imageSize(1)
                    zMin = arepo%zoompos(3, followIndex) - imageSize(1)
                    zMax = arepo%zoompos(3, followIndex) + imageSize(1)
                    imageSize(1) = imageSize(1) + dImageSize

                    deallocate(arepo%zoomids)
                    deallocate(arepo%zoompos)
                end if
                write (*, *) "apricot: Image dimensions: ", xMin, xMax, yMin, yMax, zMin, zMax

                ! Generate the image
                call makeImageArepo(arepo%ngas, arepo%pos(:, 1:arepo%ngas), arepo%rho, arepo%cellSize, & 
                                    xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image,  &
                                    xAxis, yAxis, zAxis)

                ! Write out the image
                imageFilename = "column" // "_" // trim(snapStem) // "_" // num // '.dat'
                call writeImage(imageFilename, xMin, xMax, yMin, yMax, zMin, zMax, nPixX, nPixY, nPixZ, image)

                ! Deallocate arrays ready to load next frame
                deallocate(image)
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
        case default
            write(*,*) "apricot: Unknown image mode, try again!"
    end select

    ! Finish up
    write(*,*) "apricot: Complete, bye!"
    100         format(I3.3)
end program imageArepo

! Function to read the file based on its type 
subroutine fileReader(fileType, snapStem, snapNumber, movieMode)
   use cellData
   implicit none 

   integer :: fileType, snapNumber
   character(200) :: snapStem, movieMode

   ! Read in the file based on filetype
   select case(fileType)
      ! Type 2 snapshot, with UCLCHEM
      case(1)
         call read_arepo_snap_type2(snapStem, snapNumber)
         call read_uclchem_abuns(snapStem, snapNumber) 
      ! Type 2 snapshot only
      case(2)
         call read_arepo_snap_type2(snapStem, snapNumber)
      ! Type 3, hdf5, snapshot only
      case(3)
         call read_arepo_snap_type3(snapStem, snapNumber, movieMode)
   end select  
end subroutine fileReader