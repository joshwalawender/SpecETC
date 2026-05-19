## Get passband for GaiaG filter

# p = Path('data/GaiaEDR3_passbands_zeropoints_version2/')
# pb = QTable.read(p / 'passband.dat', format='ascii.no_header',
#                  names=('wav', 'GPb', 'e_GPb', 'BPb', 'e_BPb', 'RPb', 'e_RPb'))
# t = Table(pb, masked=True)
# t['GPb'].mask = t['GPb'] > 98
# t['e_GPb'].mask = t['GPb'].mask
# t['BPb'].mask = t['BPb'] > 98
# t['e_BPb'].mask = t['BPb'].mask
# t['RPb'].mask = t['RPb'] > 98
# t['e_RPb'].mask = t['RPb'].mask
# t['wav'].unit = u.nm
# t.replace_column('wav', t['wav'].to(u.angstrom))

# t['wav', 'GPb'][~t['GPb'].mask].write(p / 'GaiaG.txt', format='ascii.csv', overwrite=True)
# bp = SpectralElement.from_file(str(p / 'GaiaG.txt'))
