#! /usr/local/bin/python3
# itcfile.py - convert ITC run data files to Matteoli format

import sys
import argparse


def readitc(infile):
    """Read data from an ITC run"""

    itc = []
    with open(infile, 'r') as f:
        title = f.readline()
        tok = f.readline().strip().split()
        if tok[0].startswith('amp') and len(tok) >= 3:
                n1c, n2c = [ float(tok[i]) for i in [1, 2] ]   # mol
        else:
            print('error: line 2 in ITC file {0} is not\n'
                'ampoule: n1c n2c   # amounts in ampoule (mol)'.format(infile))
            sys.exit(0)
        tok = f.readline().strip().split()
        if tok[0].startswith('dis') and len(tok) >= 5:
                x1d, x2d, rhod, vd = [ float(tok[i]) for i in [1, 2, 3, 4] ]
                n1d = x1d * rhod * vd * 1.0e-6    # mol
                n2d = x2d * rhod * vd * 1.0e-6    # mol
        else:
            print('error: line 3 in ITC file {0} is not\n'
                'dispenser: x1 x2 rho v   # composition rho/(mol/L) '
                'v/(muL)'.format(infile))
            sys.exit(0)
            
        for line in f.readlines():
            if line.startswith('#'):
                continue
            tok = line.strip().split()
            q = float(tok[0])             # mJ
            
            pt = {}
            pt['n1c'] = n1c               # mol
            pt['n2c'] = n2c               # mol
            pt['n1d'] = n1d * 1000.0      # mmol
            pt['n2d'] = n2d * 1000.0      # mmol
            pt['q'] = q * 0.001           # J
            itc.append(pt)
            
            n1c += n1d
            n2c += n2d
    return itc


def main():
    parser = argparse.ArgumentParser(description = \
        'Read data files from ITC runs and write in the format of\n'
        'E. Matteoli, L. Lepori, Fluid Phase Eq. 174 (2000) 115-131\n\n'
        'Format for input files:\n\n'
        'title\n'
        'n1c/mol n2c/mol            # amounts in ampoule\n'
        'x1d x2d rho/(mol/L) v/muL  # mole fractions, density, volume of '
        'addition\n'
        'Q/mJ\n'
        '...',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infiles', nargs='+', help = 'infile [...], '
        'files with ITC data')
    args = parser.parse_args()

    itc = []
    for infile in args.infiles:
        itci = readitc(infile)
        itc += itci

    print('# n1c/mol     n2c/mol      n1d/mmol     n2d/mmol     Q/J')
    for pt in itc:
        print('{0:12.5e} {1:12.5e} {2:12.5e} {3:12.5e} {4:12.5e}'.format(
            pt['n1c'], pt['n2c'], pt['n1d'], pt['n2d'], pt['q']))
    

if __name__ == '__main__':
    main()
