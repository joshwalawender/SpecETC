# SpecETC




## Sky data downloaded from Gemini
https://www.gemini.edu/observing/telescopes-and-sites/sites#SkyBackground

```
sky_file = Path(f'skybg_50_10.dat').expanduser()
with open(sky_file, 'r') as f:
    sky_contents = f.readlines()
sky_table = Table.read(sky_contents[13:], format='ascii')
# synphot.units.PHOTLAM is phot/s/cm^2/A
# table is phot/s/nm/arcsec^2/m^2
A = sky_table['nm']*10
sky_table.add_column(A, name='A')
sky_table.remove_column('nm')
# synphot.units.PHOTLAM is phot/s/cm^2/A
# table is phot/s/nm/arcsec^2/m^2
photlam = sky_table['phot/s/nm/arcsec^2/m^2']/1e5 # phot/s/A/arcsec^2/cm^2
sky_table.add_column(photlam, name='photlam/acrsec^2')
sky_table.remove_column('phot/s/nm/arcsec^2/m^2')
sky_file2 = Path(f'skybg_50_10_photlam.dat').expanduser()
sky_table.write(sky_file2, overwrite=True, format='ascii')
```
