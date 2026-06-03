from pathlib import Path
from astropy.table import Table
import synphot


path_to_data = Path(__file__).parent.parent / 'data' / 'pickles'


def get_list_of_stellar_models(verbose=False):
    pickles_readme_file = path_to_data / 'AA_README'
    with open(pickles_readme_file) as f:
        readme = f.readlines()
    pickles_readme_table = Table.read(readme[113:193], format='ascii', names=('file', 'SPTYPE', 'Teff'))
    types = [str(spt) for spt in set(pickles_readme_table['SPTYPE'])]

    if verbose:
        for ty in ['O', 'B', 'A', 'F', 'G', 'K', 'M']:
            subtypes = [st for st in types if st[0] == ty]
            print(ty, subtypes)

    return pickles_readme_table


def get_stellar_spectrum(sptype='F8V', mag=10):
    pickles_readme_table = get_list_of_stellar_models()
    files = pickles_readme_table[pickles_readme_table['SPTYPE'] == sptype]
    file = str(files[0]['file'])+'.fits'
    if file.index('uk') > 0:
        file = path_to_data / 'dat_uvk' / file
    else:
        file = path_to_data / 'dat_uvi' / file
    star = synphot.SourceSpectrum.from_file(str(file))*10**(-mag/2.5)
    return star
