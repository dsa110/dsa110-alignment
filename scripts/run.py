import sys, numpy as np, matplotlib.pyplot as plt
import funcs

# defaults
file_name = 'polaris.jpg' # name of input file to process
output_name = 'p1' # output dir and base name for files created by astrometry.net
ref_ra = 0.0 # deg
ref_dec = 90.0 # deg
ref_rad = 20.0 # for use in searching around ref point. Set to None if unwanted
tim = '2019-08-11T08:00:00.0' # time image was taken - UTC
plot = True # make diagnostic plot

# solve astrometry
f = funcs.solve_field(fin=file_name,oname=output_name,refRA=ref_ra,refDEC=ref_dec,refRAD=ref_rad)

# return position offset
funcs.calc_field(fin=f,pra=ref_ra,pdec=ref_dec,tim=tim,plot=plot)

