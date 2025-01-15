### Apricot Image Generator

##### Column Density Images of AREPO Snapshots

This module creates column density images using snapshots produced by the AREPO code. It "walks" the path along rays that travel perpendicular to the desired image plane (i.e z for an x-y image). At each position along the ray the nearest cell is found and its contribution added to the cumulative sum. The step sizes along the ray are adaptive and respond to the size of the cell the ray is intersecting. The image maker is in principle usable for any type of 3D simulation data, but only AREPO snapshot types 2 (binary) and 3 (HDF5) can currently be read by it. 

### Code Options and Functionality

apricot has two "master" functions: image generation and movie generation. Image genertion is further split by the type of image, and movie generation is split by the movement of each frame's limits. 

| Option | Options File Line | Values | Notes |
| ------ | ----------------- | ------ | ---- |
| Image Mode | Line 1| 1 or 2 | 1 is for image generation, 2 is for movie frame generation |
| Image Quantity | Line 2 |column, sgchem, uclchem, ppv, nhnh | The type of image to be generated |
| Movie Type | Line 2 | linear, followSink, followTracer | Method by which to move the frame limits, either linearly moving between a start/end point, or following a sink/tracer particle inside the snapshot |
| Include Sinks | Line 3 | 0 | Depreciated |
| Snapshot Stem | Line 4 | Any | The stem of whatever snapshot file you are using, i.e snapshots called "snapshot_001.hdf5" would have "snapshot" as the stem |
| Snapshot Number | Line 5 | Either 1 or 2 numbers | 1 Number for images, 2 numbers for movies (start/end frame) |
| File Type | Line 6 | 2 or 3 | Snapshot filetype, 2 is binary, 3 is HDF5 |
| Limits | Line 7-9 | minimum, maximum | The minimum and maximum limits of the image or frame in each dimension (x, y, z). Movies using "linear" add end values for the limits as well: xMin, xMax, xMinEnd, xMaxEnd |
| Radius | Line 7 | start, end | For movies following a tracer or a sink, the radius around it to include in the frame. Set start = end for a static radius value |
| Image Pixels | Line 10 or 8 | xPix, yPix, zPix | The size of the final image in pixels |
| Image Axes | Line 11 or 9 | 1, 2, 3 | Where to put each axis in the final image. Axes are labeled as x (1), y (2), z (3). An z-y image would be 3, 2, 1 |
| Follow ID | Line 12 | ID | If using a follow movie mode, the ID of the sink/tracer particle to follow across the movie |
