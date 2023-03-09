import os, sys
import numpy as np
import matplotlib.pyplot as plt
import aplpy
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack
import astropy.units as u

from mylib.utils.cutout import Legacy_cutout_downloader, PS1_cutout_downloader

class finder(object):
    def __init__(self, ra, dec, image):
        self.ra = ra
        self.dec = dec
        self.image = image

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=2)
        self.jname = 'J' + coordstr.replace(' ', '')
        self.jname_short = self.jname[:5] + self.jname[10:15]

    def plot(self, output, offset=None, cmap='gist_yarg', color_scale='cyan'):
        fig = plt.figure(figsize=(9, 10))
        fig.suptitle(self.jname, fontsize=24)
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.9, wspace=0, hspace=0)

        si = aplpy.FITSFigure(self.image, figure=fig, aspect='equal')
        si.show_colorscale(interpolation='nearest', aspect='equal', cmap=cmap,\
                           stretch='linear', pmin=5.0,pmax=95.,smooth=1)
        #si.recenter(ra, dec, radius=size/2.0 / 3600.)
        si.show_circles(self.ra, self.dec, radius=3.0 / 3600., color='lime',linewidth=1)
        si.add_label(self.ra, self.dec + 7.0 / 3600, self.jname_short, color='lime',weight='bold')

        # the cos(dec) factor;
        # dRA * dec_factor = distance
        # distance / dec_factor = dRA

        dec_factor = np.cos((self.dec/180.0)*np.pi)

        delta = 30.0 / dec_factor
        si.show_arrows(self.ra-(90.0 / dec_factor)/3600,\
                        self.dec+90.0/3600, delta/3600,0,color='cyan')
        si.add_label(self.ra-(70.0 / np.cos((self.dec/180.0)*np.pi))/3600, self.dec+85.0/3600,\
                        'E', color='cyan', weight='bold',size='large')

        si.show_arrows(self.ra-(90.0 / dec_factor)/3600,\
                        self.dec+90.0/3600,0,30.0/3600, color='cyan')
        si.add_label(self.ra-(95.0 / dec_factor)/3600, self.dec+110.0/3600,\
                        'N', color='cyan', weight='bold',size='large')
        if offset:
            offra = float(offset['ra'])
            offdec = float(offset['dec'])
            offcoord = SkyCoord(ra=offra, dec=offdec, unit='deg')

            si.show_circles(offra, offdec, radius=3.0 / 3600., color='yellow',linewidth=1)

            # delta_ra>0: TARGET is at the EAST of the offset
            delta_ra = (self.ra - offra) * dec_factor * 3600
            # delta_dec>0: TARGET is at the NORTH of the offset
            delta_dec = (self.dec - offdec) * 3600

            si.add_label(self.ra + (100.0 / dec_factor) / 3600, self.dec + 90.0 / 3600, \
                         'From OFFSET to Target:', color='yellow', weight='bold', size='large')

            if delta_ra<0:
                si.add_label(self.ra + (90.0 / dec_factor)/3600, self.dec+80.0/3600, \
                             '%.2f arcsec West'%abs(delta_ra), color='yellow', weight='bold',size='large')
            else:
                si.add_label(self.ra + (90.0 / dec_factor)/3600, self.dec+80.0/3600, \
                             '%.2f arcsec East'%abs(delta_ra), color='yellow', weight='bold',size='large')
            if delta_dec<0:
                si.add_label(self.ra + (90.0 / dec_factor)/3600, self.dec+70.0/3600, \
                             '%.2f arcsec South'%abs(delta_dec), color='yellow', weight='bold',size='large')
            else:
                si.add_label(self.ra + (90.0 / dec_factor)/3600, self.dec+70.0/3600, \
                             '%.2f arcsec North'%abs(delta_dec), color='yellow', weight='bold',size='large')

            si.tick_labels.hide()
            si.ticks.hide()
            delta = 30.0 / dec_factor
            si.show_arrows(self.ra-(90.0 / dec_factor)/3600,\
                           self.dec+90.0/3600,delta/3600,0,color='cyan')
            si.add_label(self.ra-(70.0 / dec_factor)/3600, self.dec+85.0/3600,\
                         'E', color='cyan', weight='bold',size='large')

            si.show_arrows(self.ra-(90.0 / dec_factor)/3600,\
                           self.dec+90.0/3600,0,30.0/3600,color='cyan')
            si.add_label(self.ra-(95.0 /dec_factor)/3600, self.dec+110.0/3600,\
                         'N', color='cyan', weight='bold',size='large')

            si.show_lines([np.array([[self.ra-40.0/3600,self.ra-40.0/3600-delta/3600],\
                                     [self.dec-60.0/3600,self.dec-60.0/3600]])],color=color_scale)
            si.add_label(self.ra-55.0/3600, self.dec-55.0/3600,'30 arcsec',color=color_scale,weight='bold')

        plt.savefig(output)
        plt.close('all')

class legacyfinder(finder):
    def __init__(self, ra, dec, direct, band='z'):
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=2)
        jname = 'J' + coordstr.replace(' ', '')

        image = direct + '/' + jname + '.%s.fits'%band

        super().__init__(ra, dec, image)

        self.band = band
        self.direct = direct

    def download_cutout(self, size=360, format='fits'):
        Legacy_cutout_downloader.download_images(self.ra, self.dec, filters=self.band, outname_root='J', outdir=self.direct,\
                                                size=size, format=format)

class PS1finder(finder):
    def __init__(self, ra, dec, direct, band='z'):
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=2)
        jname = 'J' + coordstr.replace(' ', '')

        image = direct + '/' + jname + '.%s.fits'%band

        super().__init__(ra, dec, image)

        self.band = band
        self.direct = direct

    def download_cutout(self):
        PS1_cutout_downloader.download_images(self.ra, self.dec, filters=self.band, outname_root='J', outdir=self.direct,\
                                            size=360)
