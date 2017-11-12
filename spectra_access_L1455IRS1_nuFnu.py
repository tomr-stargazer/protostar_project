"""
We downloaded the Spitzer (CASSIS) and Herschel PACS spectra for L1455-IRS 1.

Let's access it.

Notes from Green+2013:
Comparison to previous photometry in the PACS wavelength
range is problematic. The IRAS beam was much larger, and
MIPS data often suffered from saturation on these sources. The
typical MIPS 70um flux density is almost always low compared
with the PACS spectra extracted from the full array by ~20%,
except in the brighter half of our sample, where the MIPS value
is even lower-a factor of two-compared with the Herschel
data.

"""

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.table
from astroquery.simbad import Simbad

from access_irac_mips_data import E09_ids, young15_table

import astropy.units as u

def wavelength_to_freq(wavelength_array, lambda_units=u.um):
    lambda_array = u.Quantity(wavelength_array, unit=lambda_units)
    freq_array = lambda_array.to(u.Hz, equivalencies=u.spectral())
    return freq_array.value

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

try:
    simbad_query = Simbad.query_object('LDN 1455 IRS 1')
except Exception as e:
    raise e

def plot_nuFnu_vs_wavelength(wavelength_array, jansky_array, fmt=None, **kwargs):

    jansky_array = u.Quantity(jansky_array, unit=u.Jy)
    nuFnu_array = jansky_array * wavelength_to_freq(wavelength_array)*u.Hz

    return plt.plot(wavelength_array, nuFnu_array.to(u.erg/u.s/u.cm**2).value, fmt, **kwargs)


plot_nuFnu_vs_wavelength((cassis_table_lores['wavelength']), (cassis_table_lores['flux']), 'k.', ms=2, label='Spitzer IRS (Cassis), low-res')
plot_nuFnu_vs_wavelength((long_cassis_hires['wavelength']), (long_cassis_hires['flux']), 'g.', ms=2, label='Spitzer IRS (Cassis), hi-res')
plot_nuFnu_vs_wavelength((pacs_table['Wavelength(um)']), (pacs_table['Flux_Density(Jy)']), 'b-', lw=1, label='Herschel PACS (Green+16)')
plot_nuFnu_vs_wavelength((spire_table['Wavelength(u']), (spire_table['Flux_Density']), 'r-', lw=1, label='Herschel SPIRE (Green+16)')
plot_nuFnu_vs_wavelength((irac_mips_bands_float), (irac_mips_photometry_mJy/1000), 'rd', label='Spitzer IRAC & MIPS (Young+15)')

# a hack:
# Let's replace this with real table access
sharc_350um = 350 # microns
sharc_flux = 13.0 # janskys
plot_nuFnu_vs_wavelength((sharc_350um), (sharc_flux), 'y*', ms=10, label='SHARC 350um (Wu+07)')

# a hack:
# Let's replace this with real table access
scuba_850um = 850 # microns
scuba_flux = 1.59 # janskys
plot_nuFnu_vs_wavelength((scuba_850um), (scuba_flux), 'ms', label='SCUBA 850um (Kirk+06,07)')

# Bolocam! 1.1 mm
bolo_lambda = 1100 # microns
bolo_flux = bolocam_table[21]['F40'] # janskys
plot_nuFnu_vs_wavelength((bolo_lambda), (bolo_flux), 'co', label='Bolocam 1.1mm (Enoch+06)')

plt.loglog()
plt.xlabel(r"$\lambda$ ($\mu m$)", fontsize=15, family='serif')
plt.ylabel(r"$\nu F_\nu$ (erg$\cdot$s$^{-1}\cdot$cm$^{-2}$)", fontsize=15, family='serif')

plt.title("SED of L1455-IRS 1", fontsize=18, family='serif')
plt.legend(loc='lower right')

plt.savefig("L1455IRS1_SED_nuFnu.png", bbox_inches='tight')

plt.show()



def output_flux_for_transphere(wavelength, flux, label='', **kwargs):
    wavelength = u.Quantity(wavelength, u.um)
    flux = u.Quantity(flux, u.Jy)
    return_string = "{0:.3f} {1:.3f} {2}".format(wavelength.to(u.um).value, flux.to(u.Jy).value, label)

    return return_string


def iterable_output_flux_for_transphere(wavelength_list, flux_list, label='', **kwargs):
    lambda_array = u.Quantity(wavelength_list, unit=u.um)
    jansky_array = u.Quantity(flux_list, unit=u.Jy)
    return_list = [output_flux_for_transphere(x, y, label=label) for x, y in zip(lambda_array, jansky_array)]

    return return_list

print(iterable_output_flux_for_transphere((cassis_table_lores['wavelength']), (cassis_table_lores['flux']), label='Spitzer IRS (Cassis), low-res'))
print(iterable_output_flux_for_transphere((long_cassis_hires['wavelength']), (long_cassis_hires['flux']),label='Spitzer IRS (Cassis), hi-res'))
print(iterable_output_flux_for_transphere((pacs_table['Wavelength(um)']), (pacs_table['Flux_Density(Jy)']),label='Herschel PACS (Green+16)'))
print(iterable_output_flux_for_transphere((spire_table['Wavelength(u']), (spire_table['Flux_Density']), label='Herschel SPIRE (Green+16)'))
print(iterable_output_flux_for_transphere((irac_mips_bands_float), (irac_mips_photometry_mJy/1000), label='Spitzer IRAC & MIPS (Young+15)'))
print(output_flux_for_transphere((sharc_350um), (sharc_flux), label='SHARC 350um (Wu+07)'))
print(output_flux_for_transphere((scuba_850um), (scuba_flux), label='SCUBA 850um (Kirk+06,07)'))
print(output_flux_for_transphere((bolo_lambda), (bolo_flux), label='Bolocam 1.1mm (Enoch+06)'))
