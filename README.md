# KaggleHeart

## Top Level

> runGroup.m

Main function, runs through all subjects in a group folder, outputs dia,sys LV volume in vol.txt
<br><br>

> test_*.m

Tests a post-processing component on a list of cines
<br><br>

> loadCine.m

Loads a single cine, used by all scripts that load cine.
<br><br>

> checkall.m

Loops through DICOMs of all SAX folders of all subjects in a group to check/print some property of each cine.
<br><br>

## Core components used in final script:

> cineClusterBldMyo.m

k-means clustering of blood/myocardial intensities<br><br>
> findLVarea_RegGro.m

Region Growth to find LV area<br><br>
> findLVctr_cluster.m

Finds LV center of a single cine from patches blood pixels (determined by k-means) based on their roundness, etc.<br><br>
> maskRegRndmap.m

Finds roundness of each region in a binary image
<br><br>

## Naming convention

Name prefix describes the kind of data the function operates on:

mask*: 2D binary images

img*: 2D, sometimes 3D

cine*: 3D stack of images in time
