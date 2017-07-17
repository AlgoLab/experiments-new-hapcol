usage = '''

usage: python sitesinfo.py [-s S] file.sites

given file file.sites with sites (and per-site coverage), like the
file.wif.info_/sites_ file returned by 'wiftools -i file.wif', add
various other info to this based on the files provided via optional
parameters

  -s S        bed file containing per-site switch error

'''

import sys

#
# load sites with per-site coverage and bridge coverage
def load_sites(lines) :

    site_cov = {}
    site_brcov = {}
    for line in lines :
        if line.startswith('#') :
            continue

        a,b,c,d = line.split()
        site_cov[int(a)] = int(b)
        site_brcov[int(a)] = int(d)

    return site_cov, site_brcov

#
# load per-site switch error from .bed file
def load_swerrs(lines) :

    swerrs = set([])
    for line in lines :
        a,b,c,d = line.split()

        swerrs.add(int(b))

    if not swerrs : # bed file is empty, exit quietly
        print('bed file is empty, exiting ..', file = sys.stderr)
        exit(0)

    return swerrs

#
# PARSER
#----------------------------------------------------------------------

site_cov = None
swerrs = None

a = sys.argv[1:]
assert len(a), usage
i = 0
while i < len(a) : # bash-style argparse

    if a[i].startswith('-') :
        if a[i] == '-s' :
            i += 1
            swerrs = load_swerrs(open(a[i],'r'))
        else :
            assert False, usage
    else :
        site_cov, site_brcov = load_sites(open(a[i],'r'))

    i += 1 # shift once by default

#
# MAIN
#----------------------------------------------------------------------

print('#site', 'coverage', 'bridgecov', 'swerr', sep = '\t')
for site in sorted(site_cov) :

    swerr = 1 if site in swerrs else 0
    print(site, site_cov[site], site_brcov[site], swerr, sep = '\t')
