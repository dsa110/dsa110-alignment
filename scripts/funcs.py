import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
import subprocess, sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
import matplotlib.pyplot as plt
import os

 
# function to read bright stars catalog and provide all sources within a certain separation of a position.
# cat: path to input catalog
# cra, cdec (deg): ra and dec of position (can be None to return all sources)
# sep (deg): separation from positon (can be None to return all sources)
def read_bsc(cat='/Users/vikram/Projects/DSA-110/CODE/dsa_antenna_alignment/data/bsc5.dat',cra=0.0,cdec=45.0,sep=10.0):

    output = []
    rr = []
    dd = []

    c1 = None
    if cra is not None:
        if cdec is not None:
            if sep is not None:
                c1 = SkyCoord(cra*u.degree,cdec*u.degree)
    
    f = open(cat,'r')
    for cnt, line in enumerate(f):

        try:
            hr = line[0:4].replace(" ","")
            name = line[4:14].replace(" ","")
            HD = line[25:31].replace(" ","")
            SAO = line[31:37].replace(" ","")
            FK5 = line[37:41].replace(" ","")
            ra = (int(line[75:77])+int(line[77:79])/60.+float(line[79:83])/3600.)*15.
            dec = int(line[83:86])+int(line[86:88])/60.+float(line[88:90])/3600.
            vmag = float(line[102:107])
            sptype = line[127:147].replace(" ","")
            pmra = float(line[148:154])
            pmdec = float(line[154:160])            

            d = {"hr":hr,"name":name,"HD":HD,"SAO":SAO,"FK5":FK5,"ra":ra,"dec":dec,"vmag":vmag,"sptype":sptype,"pmra":pmra,"pmdec":pmdec}
            output.append(d)
            rr.append(ra)
            dd.append(dec)                
                    
        except:
            None


    cat = SkyCoord(ra = rr*u.degree,dec=dd*u.degree)

    if c1 is not None:        
        seps = cat.separation(c1)
        wrs = np.array(np.where(seps.degree<sep)[0]).astype('int')
        cat = cat[wrs]

        oo = []
        for ct in wrs:
            oo.append(output[ct])

        output = oo
    
    return output,cat


    


# function to run solve-field and return fits file
# fin: input file name (tested png and jpeg)
# oname: basename and dir name for files created by solve-field
def solve_field(fin='polaris.png',oname='p1',refRA=None,refDEC=None,refRAD=None):

    # run command
    cmd = 'solve-field -D p1 -o '+oname+' --no-plots -O -d 30 '+fin
    if refRA is not None:
        if refDEC is not None:
            if refRAD is not None:
                cmd = 'solve-field -D p1 -o '+oname+' --no-plots -O -d 30 --ra '+str(refRA)+' --dec '+str(refDEC)+' --radius '+str(refRAD)+' '+fin
    
    os.system(cmd)

    # try to return newly created fits file
    
    try:
        fitsfile = fits.open(oname+'/'+oname+'.new')
        return fitsfile
    
    except:
        print('Could not open new fits file -- exiting')
        return
        

# function to return position offset
# fin: input fits file
# pra, pdec (deg): ra and dec of position
# tim (UTC -- ISO 8601 format string): time to calculate field rotation angle at
# plot: show diagnostic plot
def calc_field(fin=None,pra=0.0,pdec=90.0,tim='2019-08-11T08:00:00.0',plot=True):

    # extract stuff from fits file
    try:
        h = fin[0].header
        d = fin[0].data
        xl = h['NAXIS1']
        yl = h['NAXIS2']
        try:
            del h['COMMENT']
        except:
            print()
        try:
            del h['HISTORY']
        except:
            print()
        h['NAXIS'] = 2
        wcs = WCS(h)

    except:
        print('Cannot read fits file:',fin)
        return 

    # get alt-az frame
    time = Time(tim)
    ovro = EarthLocation(lat=37.2314*u.deg,lon=-118.2941*u.deg,height=1222*u.m)
    aaz = AltAz(obstime=time,location=ovro)

    # get positions of field center and reference
    ra_0,dec_0 = wcs.all_pix2world(xl/2,yl/2,0)
    p0 = SkyCoord(ra=ra_0*u.degree,dec=dec_0*u.degree)
    pref = SkyCoord(ra=pra*u.degree,dec=pdec*u.degree)

    # get differences in alt and az
    aaz0 = p0.transform_to(aaz)
    aazref = pref.transform_to(aaz)
    alt_diff = aaz0.alt - aazref.alt
    az_diff = aaz0.az - aazref.az

    print('Elevation offset (center - reference, deg):',alt_diff.deg)
    print('Azimuth offset (center - reference, deg):',az_diff.deg)
    
    # make plot if needed
    if plot is True:

        # get ra, dec, and image scale
        sep = (np.sqrt(h['CD1_1']**2.+h['CD1_2']**2.))*np.max(xl)/2.
        
        # extract stars from BSC
        output,cat = read_bsc(cra=ra_0,cdec=dec_0,sep=sep)
        print('Got stars from BSC:',len(output))
        print()

        # get min and max of image
        m1 = np.median(d)-.05*np.std(d)
        m2 = np.median(d)+7.5*np.std(d)

        plt.figure(figsize=(8.0,8.0))
        ax = plt.subplot(projection=wcs)
        plt.imshow(d,vmin=m1,vmax=m2,origin='lower',cmap='afmhot')
        plt.grid(color='white', ls='solid')
        plt.xlabel('RA (deg J2000)')
        plt.ylabel('DEC (deg J2000)')

        # overlay stars
        ax.scatter(cat.ra.degree,cat.dec.degree, transform=ax.get_transform('fk5'), s=50,edgecolor='yellow', facecolor='none')
        for o in output:
            s = o['name']+' HD'+o['HD']+'_V'+str(o['vmag'])
            ax.text(o['ra'],o['dec'],s,transform=ax.get_transform('fk5'),fontsize=7,color='cyan')

        # plot positions of field centre and reference
        r1 = np.asarray([p0.ra.degree,pref.ra.degree])
        d1 = np.asarray([p0.dec.degree,pref.dec.degree])
        ax.scatter(r1,d1, s=300,marker='+',transform=ax.get_transform('fk5'), color='white')
        ax.text(p0.ra.degree,p0.dec.degree,'CENTER',transform=ax.get_transform('fk5'),fontsize=12,color='white')
        ax.text(pref.ra.degree,pref.dec.degree,'REFERENCE',transform=ax.get_transform('fk5'),fontsize=12,color='white')

        # plot line towards zenith
        newcoord = SkyCoord(alt=aaz0.alt+2.*u.deg,az=aaz0.az,obstime=time,frame='altaz',location=ovro)
        newicrs = newcoord.transform_to('icrs')
        r1 = np.asarray([p0.ra.degree,newicrs.ra.degree])
        d1 = np.asarray([p0.dec.degree,newicrs.dec.degree])
        ax.plot(r1,d1,transform=ax.get_transform('fk5'))
        ax.text(r1[1],d1[1],'Zenith',transform=ax.get_transform('fk5'),fontsize=9,color='white')    
            
        plt.show()
    
        
       
