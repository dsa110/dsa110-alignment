from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description="Script that inputs two Polaris images (up and down orientation) and outputs pointing error from the north/south axis for the DSA dish on which the camera was mounted. Make sure to calculate current NCP coordinates and edit file for accurate results!")
parser.add_argument('--upfile',type=str, help="File path for image taken with camera in right-side up orientation")
parser.add_argument('--downfile',type=str, help="File path for image taken with camera in upside down orientation")

args = parser.parse_args()

# path to image files
upimg = args.upfile
downimg = args.downfile

# name of output files for solve-field command
up_files = os.path.splitext(os.path.basename(upimg))[0]
down_files = os.path.splitext(os.path.basename(downimg))[0]

# check if __temp directory exists. If no, create it 
# this directory will store files and then will be deleted at the end of the script
outdir = '__temp'

if not os.path.isdir(outdir):
    os.system('mkdir %s' %(outdir,))

# run solve_field
print('\nSolving fields...\n')


cmd_solve_up = 'solve-field -D %s -o %s --no-plots -O -d 30 --ra 0 --dec 90 --radius 20 %s' %(outdir, up_files, upimg)
print(cmd_solve_up+'\n')
os.system(cmd_solve_up)

cmd_solve_down = 'solve-field -D %s -o %s --no-plots -O -d 30 --ra 0 --dec 90 --radius 20 %s' %(outdir, down_files, downimg)
print(cmd_solve_down+'\n')
os.system(cmd_solve_down)

upwcs = os.path.join(outdir, up_files+'.wcs')
downwcs = os.path.join(outdir, down_files+'.wcs')

# get pixel scale for image
print('Getting pixel scale for images...')
wcsinfo = os.path.join(outdir, up_files+'_wcsinfo.txt')
cmd_wcsinfo = 'wcsinfo %s >> %s' %(upwcs, wcsinfo)

print(cmd_wcsinfo+'\n')
os.system(cmd_wcsinfo)

with open(wcsinfo) as f:
	data = f.readlines()
	pix_ = data[16]
	pixscale = float(pix_.split()[1]) 

# get pix coords for true north
print('Solving for locations of north celestial pole...')

up_northpixfile = os.path.join(outdir, up_files+'_northpix.txt')
down_northpixfile = os.path.join(outdir, down_files+'_northpix.txt')

# IMPORTANT NOTE: Want to calculate updated NCP to do calculations! #
cmd_northpix_up = 'wcs-rd2xy -w %s -r 359.8700747 -d 89.88773118 >> %s' %(upwcs, up_northpixfile)
cmd_northpix_down = 'wcs-rd2xy -w %s -r 359.8700747 -d 89.88773118  >> %s' %(downwcs, down_northpixfile)

print(cmd_northpix_up+'\n')
os.system(cmd_northpix_up)

print(cmd_northpix_down+'\n')
os.system(cmd_northpix_down)

with open(up_northpixfile) as f_up:
	data = f_up.readlines()
	nx_up = float(data[0].split()[-2][1:-1])

with open(down_northpixfile) as f_down:
	data = f_down.readlines()
	nx_down = float(data[0].split()[-2][1:-1])


# get difference between the pixels in x-coordinates and divide by 2 to get distance from symmetry axis
nx_mid = 0.5*(nx_up+nx_down)
offset_pix = 0.5*np.abs(nx_up-nx_down)

offset_arcmin = offset_pix*pixscale/60

# if the symmetry axis is to the left of north celestial pole in the 'up' image, it means the telescope is pointed west of north
if (nx_mid-nx_up) < 0:
	orient = 'west of north'
else:
	orient = 'east of north'

print('Pixel scale of images = %04.2f arcsec/pix' %(pixscale,))
print('Pointing Offset: %04.2f arcmin %s\n' %(offset_arcmin, orient))

os.system('rm -r %s' %(outdir,))
