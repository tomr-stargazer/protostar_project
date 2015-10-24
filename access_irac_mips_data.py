"""
Displays the IRAC and MIPS data for our protostars.

Data is from Young et al. (2015)
"THE SPITZER c2d SURVEY OF LARGE, NEARBY, INTERSTELLAR CLOUDS. XII. THE PERSEUS YSO
POPULATION AS OBSERVED WITH IRAC AND MIPS"
AJ 2015..150..40
http://adsabs.harvard.edu/abs/2015AJ....150...40Y

"""

from __future__ import division

import os
import numpy as np
import astropy.table

path_to_data = os.path.expanduser("~/Documents/Data/c2d_catalog/")

young15_table = astropy.table.Table.read(
    os.path.join(path_to_data, "table3.dat"), 
    readme=os.path.join(path_to_data, "ReadMe"), 
    format='ascii.cds')

E09_ids = {}

# not directly identified in Simbad - I had to use the 'SMM J033327+31078' name & compare to Young+15
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40639592&Name=NAME%20BARN%201%20IRS&submit=submit
E09_ids['B1-a'] = 296
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%403879193&Name=IRAS%2003235%2b3004&submit=submit
E09_ids['IRAS 03235'] = 127
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40673515&Name=NAME%20BARN%205%20IRS%201&submit=submit
E09_ids['B5 IRS1'] = 505
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%403879195&Name=NAME%20LDN%201455%20smm%201&submit=submit
E09_ids['L1455-IRS 4'] = 132
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40639757&Name=NAME%20LDN%201455%20IRS%201&submit=submit
E09_ids['L1455-IRS 1'] = 130

for key, value in E09_ids.items():

    table_subset = young15_table[young15_table['E09']==value]

    print key, table_subset


