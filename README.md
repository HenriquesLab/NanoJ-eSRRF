# NanoJ-eSRRF

## Adaptive image reconstruction for high-fidelity, fast and easy-to-use 3D live-cell super-resolution microscopy

eSRRF (enhanced Super-Resolution Radial Fluctuation) is an extension of the SRRF method developped by the Henriques lab (**[SRRF paper](http://www.nature.com/articles/ncomms12471)**). 

eSRRF aims at improving the fidelity of SRRF images with respect to the underlying true structure. Below is shown a representative dataset obtained from high-density emitters for which the underlying structure was obtained via DNA-PAINT (SMLM). 


<img src="https://github.com/HenriquesLab/NanoJ-eSRRF/blob/master/wiki_files/eSRRF_showcase.png" width="500"/>

The (e)SRRF approach is based on
* **a spatial analysis** of the high-density emitter data using Radial symmetry transform
* **a temporal analysis** of the obtained temporal stack using similar appraoch to SOFI.

<img src="https://github.com/HenriquesLab/NanoJ-eSRRF/blob/master/wiki_files/eSRRF_method.png" width="500"/>

## Features of eSRRF

Some of the new features available in eSRRF are as follows:
* improved fidelity of reconstructions
* adaptive reconstruction schemes allowing to explore the compromises between **fidelity** and **resolution**. This is enabled by an integration of our SQUIRREL approach (**[SQUIRREL paper](https://doi.org/10.1038/nmeth.4605)**).
* estimation of the maximum number of frames to use for eSRRF analysis from a dataset, based on SSIM (or optical flow magnitude) calculation
* rolling analysis
* integrated drift correction
* better memory management
* full OpenCL integration
* direct saving to disk for large dataset analysis
* **eSRRF 3D** recontruction!

## Getting eSRRF plugin on Fiji

The latest stable version of eSRRF can be directrly obtained from our Fiji update site: **https://sites.imagej.net/NanoJ-LiveSRRF/**

Information about update sites can be found [here](https://imagej.net/update-sites/).


## Tools coming with eSRRF plugin

eSRRF comes packed with useful Tools plugins to perform a range of things, such as (but not limited to):
* fluorescence fluctuation simulator
* rescale individual slices within a stack and convert it to RGB (useful to visualise the parameter sweep output)
* Save all current open images as Tiff files
* perform linear rescaling on stack
* etc.
* 

## People involved

Many people are involved in developping and testing this method, here are some of the key players:
* Romain F. Laine ([@LaineBioImaging](https://twitter.com/LaineBioImaging))
* Ricardo Henriques ([@HenriquesLab](https://twitter.com/HenriquesLab))
* Guillaume Jacquemet ([@guijacquemet](https://twitter.com/guijacquemet))
* Christophe Leterrier ([@christlet](https://twitter.com/christlet))
* Siân Culley ([@SuperResoluSian](https://twitter.com/SuperResoluSian))
* Bassam Hajj ([@Bassam_A_HAJJ](https://twitter.com/Bassam_A_HAJJ))
* ...
