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
# load sites and per-site coverage
def load_sites(lines) :

    site_cov = {}
    for line in lines :
        if line.startswith('#') :
            continue

        a,b,c = line.split()
        site_cov[int(a)] = int(b)

    return site_cov

#
# load per-site switch error from .bed file
def load_swerrs(lines) :

    swerrs = set([])
    for line in lines :
        a,b,c,d = line.split()

        swerrs.add(int(b))

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
        site_cov = load_sites(open(a[i],'r'))

    i += 1 # shift once by default

#
# MAIN
#----------------------------------------------------------------------

print('#site', 'coverage', 'swerr', sep = '\t')
for site in sorted(site_cov) :

    swerr = 1 if site in swerrs else 0
    print(site, site_cov[site], swerr, sep = '\t')
