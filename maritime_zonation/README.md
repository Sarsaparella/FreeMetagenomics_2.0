# How to distinguish maritime zones?

To sort samples from a specific maritime zone, the coordinates of samples should be overlapped with reference areas of these zones.

<br />
 
For that, a reference file is required. In our case, shapefiles from [Marine Regions](https://www.marineregions.org) under the CC BY 4.0 licence are used.
> Due to copyright reasons, please, obtain these files on your own :shipit:

<br />
 
The pipeline has the following major steps:
- uploading reference shapefiles (as a path to a folder with zips);
- uploading a table with samples (it should contain columns Latitude and Longitude);
- overlapping them to see how many areas samples fall into;
- visualisation of results;
- saving a sample table with a new column reflecting maritime zones.

<br />

🔥 See a demonstration in the notebook 🔥

> Please, consider this code only as a recommendation since it fits specific input files

![geo_coord_sample_coverage](https://user-images.githubusercontent.com/107003133/236578540-0e26bd51-766d-48d8-94c4-4c0351158dcd.png)
