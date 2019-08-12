# dsa_antenna_alignment

This repo contains code to solve for the astrometric pointing location of a wide-angle image from a digital camera of the night sky. The particular aim is to use such images from a camera mounted on a radio dish to align an axis of motion.

### Requirements

* python 3.x with usual astro packages (astropy etc)
* astrometry.net software. Install as described [here](http://astrometry.net/doc/readme.html) with the wide-angle index files.
* netpbm tools for use of astrometry.net

### Usage

Edit and run *run.py* in the scripts directory, following the comments. For example, if in the demo directory, run >python ../scripts/run.py . Note that only png and jpeg input images have been tested, and so e.g., bitmaps will have to be converted before running. 



