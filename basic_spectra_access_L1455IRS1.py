"""
We downloaded the Spitzer (CASSIS) and Herschel PACS spectra for this source.

Let's access it.

"""

from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.table

cassis_data_path = os.path.expanduser("~/Documents/Data/cassis/")
pacs_data_path = os.path.expanduser("~/Documents/Data/joelgreen_pacs/CDF_archive/IRAS03245/pacs/data")

# load the cassis data

cassis_table = astropy.table.Table.read(os.path.join(cassis_data_path, "cassis_tbl_spcf_6368000t.tbl"), format='ascii.ipac')
pacs_table = astropy.table.Table.read(
    os.path.join(pacs_data_path, "IRAS03245_centralSpaxel_PointSourceCorrected_CorrectedYES_trim.txt"), format='ascii.basic')

fig = plt.figure()

plt.plot(np.log10(cassis_table['wavelength']), np.log10(cassis_table['flux']), 'k.')
plt.plot(np.log10(pacs_table['Wavelength(um)']), np.log10(pacs_table['Flux_Density(Jy)']), 'b.')

plt.xlabel(r"log($\lambda$ / $\mu m$)", fontsize=18)
plt.ylabel("log(Flux / Jy)", fontsize=18)
plt.show()

