"""
We downloaded the Spitzer (CASSIS) and Herschel PACS spectra for this source.

Let's access it.

"""

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.table

from access_irac_mips_data import E09_ids, young15_table

cassis_data_path = os.path.expanduser("~/Documents/Data/cassis/")
pacs_data_path = os.path.expanduser("~/Documents/Data/joelgreen_pacs/CDF_archive/IRAS03245/pacs/data")
spire_data_path = os.path.expanduser("~/Documents/Data/joelgreen_pacs/CDF_archive/IRAS03245/spire/data")
bolocam_data_path = os.path.expanduser("~/Documents/Data/c2d_catalog")

# load the cassis data

cassis_table_lores = astropy.table.Table.read(os.path.join(cassis_data_path, "cassis_tbl_spcf_6368000t.tbl"), format='ascii.ipac')
cassis_table_hires = astropy.table.Table.read(os.path.join(cassis_data_path, "cassis_tbl_opt_6368000.tbl"), format='ascii.ipac') 
pacs_table = astropy.table.Table.read(
    os.path.join(pacs_data_path, "IRAS03245_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt"), format='ascii.basic')
spire_table = astropy.table.Table.read(
    os.path.join(spire_data_path, "IRAS03245_spire_corrected.txt"), format='ascii.basic')
bolocam_table = astropy.table.Table.read(os.path.join(bolocam_data_path, "bolocam_enoch06.fits.fit"))

# extract the IRAC/MIPS data
young_table_subset = young15_table[young15_table['E09'] == E09_ids['L1455-IRS 1']]
irac_mips_bands = ['3.6', '4.5', '5.8', '8.0', '24', '70']
irac_mips_bands_float = np.array([float(x) for x in irac_mips_bands])
irac_mips_photometry_mJy = np.array([young_table_subset['F{0}'.format(x)][0] for x in irac_mips_bands])

fig = plt.figure()

long_cassis_hires = cassis_table_hires[cassis_table_hires['wavelength'] > 14]

plt.plot(np.log10(cassis_table_lores['wavelength']), np.log10(cassis_table_lores['flux']), 'k.', label='Spitzer IRS (Cassis), low-res')
plt.plot(np.log10(long_cassis_hires['wavelength']), np.log10(long_cassis_hires['flux']), 'g.', label='Spitzer IRS (Cassis), hi-res')
plt.plot(np.log10(pacs_table['Wavelength(um)']), np.log10(pacs_table['Flux_Density(Jy)']), 'b-', label='Herschel PACS (Green+16)')
plt.plot(np.log10(spire_table['Wavelength(u']), np.log10(spire_table['Flux_Density']), 'r-', label='Herschel SPIRE (Green+16)')
plt.plot(np.log10(irac_mips_bands_float), np.log10(irac_mips_photometry_mJy/1000), 'rd', label='Spitzer IRAC & MIPS (Young+15)')

# a hack:
scuba_850um = 850 # microns
scuba_flux = 1.59 # janskys
plt.plot(np.log10(scuba_850um), np.log10(scuba_flux), 'ms', label='SCUBA 850um (Kirk+06,07)')

# Bolocam! 1.1 mm
bolo_lambda = 1100 # microns
bolo_flux = bolocam_table[21]['F40'] # janskys
plt.plot(np.log10(bolo_lambda), np.log10(bolo_flux), 'co', label='Bolocam 1.1mm (Enoch+06)')


plt.xlabel(r"log($\lambda$ / $\mu m$)", fontsize=18)
plt.ylabel("log(Flux / Jy)", fontsize=18)

plt.title("SED of L1455-IRS 1", fontsize=18, family='serif')
plt.legend(loc='lower right')

plt.show()

