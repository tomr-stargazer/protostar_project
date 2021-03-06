"""
We downloaded the Spitzer (CASSIS) and Herschel PACS spectra for B1-a.

Let's access it.

"""

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.table
from astroquery.simbad import Simbad

from access_irac_mips_data import E09_ids, young15_table

cassis_data_path = os.path.expanduser("~/Documents/Data/cassis/")
pacs_data_path = os.path.expanduser("~/Documents/Data/joelgreen_pacs/CDF_archive/B1-a/pacs/data")
spire_data_path = os.path.expanduser("~/Documents/Data/joelgreen_pacs/CDF_archive/B1-a/spire/data")
bolocam_data_path = os.path.expanduser("~/Documents/Data/c2d_catalog")

# load the cassis data

cassis_table_lores = astropy.table.Table.read(os.path.join(cassis_data_path, "cassis_tbl_spcfw_15918080t.tbl"), format='ascii.ipac')
cassis_table_hires = astropy.table.Table.read(os.path.join(cassis_data_path, "cassis_tbl_opt_15918080.tbl"), format='ascii.ipac') 
pacs_table = astropy.table.Table.read(
    os.path.join(pacs_data_path, "B1-a_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt"), format='ascii.basic')
spire_table = astropy.table.Table.read(
    os.path.join(spire_data_path, "B1-a_spire_corrected.txt"), format='ascii.basic')
bolocam_table = astropy.table.Table.read(os.path.join(bolocam_data_path, "bolocam_enoch06.fits.fit"))

# extract the IRAC/MIPS data
young_table_subset = young15_table[young15_table['E09'] == E09_ids['B1-a']]
irac_mips_bands = ['3.6', '4.5', '5.8', '8.0', '24', '70']
irac_mips_bands_float = np.array([float(x) for x in irac_mips_bands])
irac_mips_photometry_mJy = np.array([young_table_subset['F{0}'.format(x)][0] for x in irac_mips_bands])

fig = plt.figure()

long_cassis_hires = cassis_table_hires[cassis_table_hires['wavelength'] > 14]

try:
    simbad_query = Simbad.query_object('BARN 1 IRS')
except Exception as e:
    raise e

plt.plot(np.log10(cassis_table_lores['wavelength']), np.log10(cassis_table_lores['flux']), 'k.', label='Spitzer IRS (Cassis), low-res')
plt.plot(np.log10(long_cassis_hires['wavelength']), np.log10(long_cassis_hires['flux']), 'g.', label='Spitzer IRS (Cassis), hi-res')
plt.plot(np.log10(pacs_table['Wavelength(um)']), np.log10(pacs_table['Flux_Density(Jy)']), 'b-', label='Herschel PACS (Green+16)')
plt.plot(np.log10(spire_table['Wavelength(u']), np.log10(spire_table['Flux_Density']), 'r-', label='Herschel SPIRE (Green+16)')
plt.plot(np.log10(irac_mips_bands_float), np.log10(irac_mips_photometry_mJy/1000), 'rd', label='Spitzer IRAC & MIPS (Young+15)')

# a hack:
scuba_850um = 850 # microns
scuba_flux = 9.92 # janskys
plt.plot(np.log10(scuba_850um), np.log10(scuba_flux), 'ms', label='SCUBA 850um (Kirk+06,07)')


plt.xlabel(r"log($\lambda$ / $\mu m$)", fontsize=18)
plt.ylabel("log(Flux / Jy)", fontsize=18)

plt.title("SED of B1-a", fontsize=18, family='serif')
plt.legend(loc='lower right')

plt.show()

