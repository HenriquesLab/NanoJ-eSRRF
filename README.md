# NanoJ-eSRRF

## Adaptive image reconstruction for high-fidelity, fast and easy-to-use 3D live-cell super-resolution microscopy

eSRRF (enhanced Super-Resolution Radial Fluctuations) is an extension of the SRRF method developed by the Henriques lab, described in **[Gustafsson _et al._ (2016)](http://www.nature.com/articles/ncomms12471)**. For more details check out our **[preprint](https://doi.org/10.1101/2022.04.07.487490)** on bioRXiv. 

eSRRF aims at improving the fidelity of SRRF images with respect to the underlying true structure. Below is shown a representative dataset obtained from high-density emitters for which the underlying structure was obtained via DNA-PAINT (SMLM). 


<img src="https://github.com/HenriquesLab/NanoJ-eSRRF/blob/master/wiki_files/eSRRF_promo.gif" width="500"/>

The (e)SRRF approach is based on
* **A spatial analysis** of the high-density emitter data using Radial symmetry transform;
* **A temporal analysis** of the obtained temporal stack using similar appraoch to SOFI.

<img src="https://github.com/HenriquesLab/NanoJ-eSRRF/blob/master/wiki_files/eSRRF_method.png" width="500"/>

## Features of eSRRF

Some of the new features available in eSRRF include:
* Improved fidelity of reconstructions;
* Adaptive reconstruction schemes allowing to explore the compromise between **fidelity** and **resolution**. This is enabled by an integration of our SQUIRREL approach, described in **[Culley _et al._ (2018)](https://doi.org/10.1038/nmeth.4605)**;
* Estimation of the maximum number of frames to use for eSRRF analysis from a dataset, based on SSIM (or optical flow magnitude) calculation;
* Rolling analysis;
* Integrated drift correction;
* Better memory management;
* Full OpenCL integration, enabling GPU acceleration;
* Direct saving to disk for large dataset analysis;
* **eSRRF 3D** reconstruction!

## Getting the eSRRF plugin on Fiji

The latest stable version of eSRRF can be directly obtained from our Fiji update site: **https://sites.imagej.net/NanoJ-LiveSRRF/**

Information about update sites can be found [here](https://imagej.net/update-sites/).

Video guide: Installation

[<img alt="Video guide: Installation" width="400px" src="https://github.com/HenriquesLab/NanoJ-eSRRF/blob/master/wiki_files/eSRRF-Guide%20installation_screeshot.png" />](https://www.youtube.com/watch?v=3Oa2ADEa-qY&t=4s)


## Test datasets
We have published test datasets including eSRRF parameter suggestions on [Zenodo](https://doi.org/10.5281/zenodo.6466472). Download and get started right away!

## Tools included in the eSRRF plugin

eSRRF comes packed with useful Tools plugins to perform a range of things, such as (but not limited to):
* Fluorescence fluctuation simulator;
* Rescale individual slices within a stack and convert it to RGB (useful to visualise the parameter sweep output);
* Save all current open images as Tiff files;
* Perform linear rescaling on stack;

## People involved

Many people are involved in developing and testing this method, here are some of the key players:
* Romain F. Laine ([@LaineBioImaging](https://twitter.com/LaineBioImaging))
* Ricardo Henriques ([@HenriquesLab](https://twitter.com/HenriquesLab))
* Guillaume Jacquemet ([@guijacquemet](https://twitter.com/guijacquemet))
* Christophe Leterrier ([@christlet](https://twitter.com/christlet))
* Siân Culley ([@SuperResoluSian](https://twitter.com/SuperResoluSian))
* Bassam Hajj ([@Bassam_A_HAJJ](https://twitter.com/Bassam_A_HAJJ))
* Hannah S. Heil ([@Hannah_SuperRes](https://twitter.com/hannah_superres))
* Simao Coelho ([@simaopc](https://twitter.com/simaopc))
* Jonathon Nixon-Abell ([@AbellJonny](https://twitter.com/AbellJonny))
* Angélique Jimenez 
* Tommaso Galgani
* Aki Stubb ([@akistub](https://twitter.com/akistub))
* Gautier Follain ([@Follain_Ga](https://twitter.com/Follain_Ga))
