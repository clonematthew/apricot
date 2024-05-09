module cell_data

    implicit none
  
    type :: arepo_data_type
        integer, dimension(6) :: npart, nall, massarr
        integer :: ngas, nsink, ntracer, ntotal, nreal
        double precision :: time
        double precision, allocatable, dimension(:,:) :: pos, vel, chem
        double precision, allocatable, dimension(:) :: mass, u, rho
        double precision, allocatable, dimension(:) :: temp, tdust, cellsize
        double precision, allocatable, dimension(:) :: sinkx, sinky, sinkz
        double precision, allocatable, dimension(:) :: sinkvx, sinkvy, sinkvz
        double precision, allocatable, dimension(:) :: sinkmass
        integer, allocatable, dimension(:) :: ids, sinkids
        ! UCLCHEM stuff
        integer :: nspecies
        double precision, allocatable, dimension(:,:) :: uclchem
    end type
 
    ! Declare variables we need eveywhere
    type(arepo_data_type) :: arepo

    character(200) :: input_snapshot
    double precision, parameter ::  udist = 1.0d17
    double precision, parameter ::  umass = 1.991d33
    double precision, parameter ::  utime = 2.743689806d12
    double precision, parameter ::  ABHE = 0.1 
    double precision, parameter ::  uerg = 2.64481d+42
    double precision, parameter ::  udens = 1.991d-18
    double precision, parameter ::  k_B =1.3806581e-16
    double precision, parameter ::  mp = 1.6726575e-24

    contains

    subroutine read_arepo_snap_type3(snapshot_stem, snapshot_number)
        use hdf5
        implicit none

        ! Assignment of filename variables
        character(*) :: snapshot_stem
        integer :: snapshot_number
        character(len=:), allocatable :: arepo_file_name
        character(len=3) :: num

        ! Variables for reading the file
        integer :: errorID
        integer(HID_T) :: fileID, groupID, attrID, datasetID, spaceID
        integer(HSIZE_T), dimension(1) :: dims
        integer(HSIZE_T), dimension(2) :: data_dims, max_dims

        ! Assign temp values for wokring out temperature
        double precision, allocatable, dimension(:) :: yn, yntot, energy
        
        ! Create the actual name of the file
        write(num, 100) snapshot_number
    100     format(I3.3)
        arepo_file_name = trim(snapshot_stem) // '_' // num // ".hdf5"
        print *, "apricot: Filename: ", arepo_file_name

        ! Initialise the library
        call h5open_f(errorID)

        ! Open the file
        call h5fopen_f(arepo_file_name, H5F_ACC_RDONLY_F, fileID, errorID)

        ! Open the header 
        print *, "apricot: Reading Header"
        call h5gopen_f(fileID, "Header", groupID, errorID)

        ! Read the data from the header
        dims = 6
        call h5aopen_f(groupID, "NumPart_Total", attrID, errorID)
        call h5aread_f(attrID, H5T_NATIVE_INTEGER, arepo%npart, dims, errorID)

        call h5aopen_f(groupID, "MassTable", attrID, errorID)
        call h5aread_f(attrID, H5T_NATIVE_INTEGER, arepo%massarr, dims, errorID)

        dims = 1
        call h5aopen_f(groupID, "Time", attrID, errorID)
        call h5aread_f(attrID, H5T_NATIVE_DOUBLE, arepo%time, dims, errorID)

        ! Close the header
        call h5gclose_f(groupID, errorID)

        ! Assign variables
        arepo%ntotal = sum(arepo%npart(1:6))
        arepo%ngas = arepo%npart(1)
        arepo%ntracer = arepo%npart(4)
        arepo%nsink = arepo%npart(6)
        arepo%nreal = arepo%ntotal - arepo%ntracer
        print *, "apricot: Header Read, NPart: ", arepo%ntotal 

        ! Open the particle data
        call h5gopen_f(fileID, "PartType0", groupID, errorID)

        ! Read the positions
        print *, "apricot: Reading Positions"
        call h5dopen_f(groupID, "Coordinates", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        call h5sget_simple_extent_dims_f(spaceID, data_dims, max_dims, errorID)
        allocate(arepo%pos(data_dims(1), data_dims(2)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%pos, data_dims, errorID)

        ! Read the velocities 
        print *, "apricot: Reading Velocities"
        call h5dopen_f(groupID, "Velocities", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        allocate(arepo%vel(data_dims(1), data_dims(2)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%vel, data_dims, errorID)

        ! Read the masses
        print *, "apricot: Reading Masses"
        call h5dopen_f(groupID, "Masses", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        call h5sget_simple_extent_dims_f(spaceID, data_dims, max_dims, errorID)
        allocate(arepo%mass(data_dims(1)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%mass, data_dims, errorID)

        ! Read the internal energy
        print *, "apricot: Reading Internal Energy"
        call h5dopen_f(groupID, "InternalEnergy", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        allocate(arepo%u(data_dims(1)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%u, data_dims, errorID)

        ! Reading density 
        print *, "apricot: Reading Densities"
        call h5dopen_f(groupID, "Density", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        allocate(arepo%rho(data_dims(1)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%rho, data_dims, errorID)

        ! Reading dust temperatures
        print *, "apricot: Reading Dust Temperatures"
        call h5dopen_f(groupID, "DustTemperature", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        allocate(arepo%tdust(data_dims(1)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%tdust, data_dims, errorID)

        ! Reading IDs
        print *, "apricot: Reading IDs"
        call h5dopen_f(groupID, "ParticleIDs", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        allocate(arepo%ids(data_dims(1)))
        call h5dread_f(datasetID, H5T_NATIVE_INTEGER, arepo%ids, data_dims, errorID)

        ! Reading chemistry
        print *, "apricot: Reading Chemistry"
        call h5dopen_f(groupID, "ChemicalAbundances", datasetID, errorID)
        call h5dget_space_f(datasetID, spaceID, errorID)
        call h5sget_simple_extent_dims_f(spaceID, data_dims, max_dims, errorID)
        allocate(arepo%chem(data_dims(1), data_dims(2)))
        call h5dread_f(datasetID, H5T_NATIVE_DOUBLE, arepo%chem, data_dims, errorID)

        ! Close the file 
        call h5gclose_f(groupID, errorID)
        call h5fclose_f(fileID, errorID)

        ! Calculate the cellsize
        allocate(arepo%cellsize(1:arepo%ngas))
        arepo%cellsize = (arepo%mass(1:arepo%ngas) / arepo%rho)**(1./3.)

        ! Calculate the gas temperature
        if (allocated(arepo%chem)) then
            allocate(arepo%temp(1:arepo%ngas))
            allocate(yn(1:arepo%ngas))
            allocate(energy(1:arepo%ngas))
            allocate(yntot(1:arepo%ngas))

            yn = arepo%rho * udens / ((1.0 + 4.0 * ABHE) * mp)
            energy = arepo%u * arepo%rho * uerg / udist**3
            yntot = (1.0 + ABHE - arepo%chem(1, :) + arepo%chem(2, :)) * yn
            arepo%temp = 2.0 * energy / (3.0 * yntot * k_B)

            deallocate(yn)
            deallocate(yntot)
            deallocate(energy)
        endif
    end subroutine read_arepo_snap_type3

    subroutine read_arepo_snap_type2(snapshot_stem, snapshot_number)

        implicit none

        character(*) :: snapshot_stem
        character(len=:), allocatable :: arepo_file_name
        character (len = 4) :: data_tag
        character (len = 3) :: num
        integer :: IOstatus
        integer :: snapshot_number
        integer, dimension(64) :: unused
        double precision, allocatable, dimension(:) :: dummy_1d
        double precision, allocatable, dimension(:,:) :: dummy_2d   
        integer, allocatable, dimension(:) :: dummy_int
        double precision, allocatable, dimension(:) :: yn, yntot, energy
        integer :: n_not_sink

        ! Open file
        write(num, 100) snapshot_number
    100    format(I3.3)
        arepo_file_name = trim(snapshot_stem) // '_' // num
        write(*,*) 'apricot: Opening snaphot file ', arepo_file_name
        open(20, file=arepo_file_name, form='unformatted',status='old')

        ! Loop over records and store what we find
        read(20, IOSTAT=IOstatus) data_tag
        do while (IOstatus==0)
            !
            ! branch depending on what's found
            select case (data_tag)
            case ('HEAD')
                write(*,*) 'apricot: Reading header'

                ! Read header
                read(20) arepo%npart, arepo%massarr, arepo%time, unused(1:4), arepo%nall, unused(1:24)
                arepo%ntotal = sum( arepo%npart(1:6) ) 
                arepo%ngas = arepo%npart(1)
                arepo%nsink = arepo%npart(6)
                arepo%ntracer = arepo%npart(4)
                arepo%nreal = arepo%ntotal - arepo%ntracer
                write(*,*) 'npart:', arepo%npart
                write(*,*) 'massarr:', arepo%massarr
                write(*,*) 'Time:', arepo%time
                write(*,*) 'ngas ', arepo%ngas, ' nsink', arepo%nsink
            case ('POS ')
                write(*,*) 'apricot: Reading positions'

                ! Read variables
                allocate( arepo%pos(1:3, 1:arepo%nreal) )
                read(20) arepo%pos
            case ('VEL ')
                write(*,*) 'apricot: Reading velocities'
                allocate( arepo%vel(1:3, 1:arepo%nreal) )
                read(20) arepo%vel
            case ('ID ')
                write(*,*) 'apricot: Reading IDs'
                allocate( arepo%ids(1:arepo%nreal) )
                read(20) arepo%ids
            case ('MASS')
                write(*,*) 'apricot: Reading masses'
                allocate( arepo%mass(1:arepo%nreal) )
                read(20) arepo%mass
            case ('U ')
                write(*,*) 'apricot: Reading u'
                allocate( arepo%u(1:arepo%ngas) )
                read(20) arepo%u
            case ('RHO ')
                write(*,*) 'apricot: Reading densities'
                allocate( arepo%rho(1:arepo%ngas) )
                read(20) arepo%rho
            case ('DUST')
                write(*,*) 'apricot: Reading dust temperatures'
                allocate( arepo%tdust(1:arepo%ngas) )
                read(20) arepo%tdust
            case ('CHEM')
                write(*,*) 'apricot: Reading chemistry'
                allocate( arepo%chem(1:9, 1:arepo%ngas) )
                read(20) arepo%chem
            case default
                write(*,*) 'apricot: No rule for record tag ', data_tag, ' . Skipping!'              
                read(20)
            end select
            
            ! Read next record tag...
            read(20, IOSTAT=IOstatus) data_tag
        end do
        write(*,*) 'apricot: Reached end of file'
        write(*,*) 'apricot: Min / Max x:', minval( arepo%pos(1, 1:arepo%nreal) ), maxval( arepo%pos(3, 1:arepo%nreal) )
        write(*,*) 'apricot: Min / Max density:',  minval( arepo%rho(1:arepo%ngas) ), maxval( arepo%rho(1:arepo%ngas) )
        
        ! Clean out sinks if preset, and put them in their own arrays
        if (arepo%nsink .gt. 0) then
            print *, "apricot: Sinks present: removing from arrays"
            
            ! Allocate the sink arrays and copy the data in  
            allocate( arepo%sinkx(1:arepo%nsink) )
            allocate( arepo%sinky(1:arepo%nsink) )
            allocate( arepo%sinkz(1:arepo%nsink) )
            allocate( arepo%sinkvx(1:arepo%nsink) )
            allocate( arepo%sinkvy(1:arepo%nsink) )
            allocate( arepo%sinkvz(1:arepo%nsink) )
            allocate( arepo%sinkmass(1:arepo%nsink) )
            allocate( arepo%sinkids(1:arepo%nsink) )
            n_not_sink = sum( arepo%npart(1:5) ) - arepo%ntracer
            print *, "apricot: total non-sink items in arrays", n_not_sink
            print *, "apricot: Moving sinks to their own arrays..."
            arepo%sinkx(1:arepo%nsink) = arepo%pos(1, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinky(1:arepo%nsink) = arepo%pos(2, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkz(1:arepo%nsink) = arepo%pos(3, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkvx(1:arepo%nsink) = arepo%vel(1, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkvy(1:arepo%nsink) = arepo%vel(2, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkvz(1:arepo%nsink) = arepo%vel(3, n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkmass(1:arepo%nsink) = arepo%mass(n_not_sink+1:n_not_sink+arepo%nsink)
            arepo%sinkids(1:arepo%nsink) = arepo%ids(n_not_sink+1:n_not_sink+arepo%nsink)
            print *, 'apricot: Done.'
            
            ! Resize the pos / vel / mass / ids arrays to only include the gas 
            allocate( dummy_1d(1:arepo%ngas) )
            allocate( dummy_2d(1:3, 1:arepo%ngas) )

            ! Posisitions
            print *, "apricot: Pos of cell 1", arepo%pos(1, 1), arepo%pos(2, 1), arepo%pos(3, 1)
            print *, "apricot: Pos of cell 2", arepo%pos(1, 2), arepo%pos(2, 2), arepo%pos(3, 2)
    end if 

        ! Make any other arrays that you might need for plotting
        allocate( arepo%cellsize(1:arepo%ngas) )
        arepo%cellsize = (arepo%mass(1:arepo%ngas) / arepo%rho)**(1./3.) 
        if ( allocated(arepo%chem) ) then
            print *, "apricot: Chemistry present, so proving a gas temperature"
            allocate( arepo%temp(1:arepo%ngas) )
            allocate( yn(1:arepo%ngas) )
            allocate( energy(1:arepo%ngas) )
            allocate( yntot(1:arepo%ngas) )
            yn = arepo%rho * udens / ((1.0 + 4.0 * ABHE) * mp)
            energy = arepo%u * arepo%rho * uerg / udist**3
            yntot = (1.0 + ABHE - arepo%chem(1, :) + arepo%chem(2, :)) * yn
            arepo%temp = 2.0 * energy / (3.0 * yntot * k_B)
            deallocate(yn)
            deallocate(yntot)
            deallocate(energy)
        end if
    
        ! Close file
        print *, 'apricot: density, cellsize', arepo%rho(5), arepo%cellsize(5)
        print *, 'apricot: Closing the file reader. Should have everything we need'
        close(20)
    end subroutine read_arepo_snap_type2

end module cell_data
