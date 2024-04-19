#### Apricot 
##### Column Density Images of AREPO snapshots

Creates column density images using AREPO snapshots.
Current supports snapshot type 2 and 3. Type 3 requires the HDF5 library to be installed.

#### Code Options

Once made the code reads image_arepo_options.dat, which follows the template given in the repo.

Image mode - Imaging mode, currently 1 (single image) and 2 (ppv) are supported.
Image type - Type of image, currently "column" for column density is supported.
Flag to keep sinks - 1 to keep them, 0 to not. Unused.

Snapshot stem - Stem of the snapshot filename
Snapshot number - Number of the snapshot
Snapshot type - Snapshot type (2 or 3)

x, y and z limits - Max and Min values 
Number of Pixels - X, Y and Z pixel numbers

Axes orientation - Where to put what axis, first number is x (so for x on x put 1, for y on x put 2, etc), second is y, third is z.

#### Dependencies

Uses the KDTree algorithim created by Matthew Kennel, Institute for Nonlinear Science (2004).
