#!/usr/bin/env python3
# hmix.py - fit Redlich-Kister polynomial to titration calorimetry data

import sys
import argparse
import numpy as np
from lmfit import minimize, Parameters, Parameter, fit_report

# --------------------------------------


def readitc(infile):
    """Read data from ITC experiment
       n1c and n2c are amounts of components in ampoule
       n1d and n2d are amounts of components in dispenser
       file format:
       n1c/mol n2c/mol n1d/mmol n2d/mmol Q/J
       ...
    """

    itc = []
    for line in open(infile, 'r'):
        if line.startswith('#'):
            continue
        tok = line.strip().split()
        n1c, n2c, n1d, n2d, q = [ float(tok[i]) for i in range(5) ]
        pt = {}
        pt['n1c'] = n1c
        pt['n2c'] = n2c
        pt['n1d'] = n1d * 1.0e-3
        pt['n2d'] = n2d * 1.0e-3
        pt['q'] = q
        itc.append(pt)
    return itc


# --------------------------------------


def herk2(x2, a):
    """HE Redlich-Kister in mole fraction of component 2"""

    n = len(a)
    he = 0.0
    for i in range(n):
        he += a[i] * (1.0 - 2.0*x2)**i
    he *= (1.0 - x2) * x2
    return he

def herk1(x1, a):
    """HE Redlich-Kister in mole fraction of component 1"""

    n = len(a)
    he = 0.0
    for i in range(n):
        he += a[i] * (2.0*x1 - 1.0)**i
    he *= x1 * (1.0 - x1)
    return he


def h2rk(x2, a):
    """Partial molar enthalpy of component 2 from RK"""
    
    n = len(a)
    h2_a = h2_b = 0.0
    for i in range(n):
        if i > 0:
            h2_a += -2.0 * i * a[i] * (1.0 - 2.0*x2)**(i-1)
        h2_b += a[i] * (1.0 - 2.0*x2)**i
    h2 = (x2 - 1.0)**2 * (x2 * h2_a + h2_b)
    return h2

def h1rk(x1, a):
    """Partial molar enthalpy of component 2 from RK"""

    n = len(a)
    h1_a = h1_b = 0.0
    for i in range(n):
        if i > 0:
            h1_a += 2.0 * i * a[i] * (2.0*x1 - 1.0)**(i-1)
        h1_b += a[i] * (2.0*x1 - 1.0)**i
    h1 = (1.0 - x1)**2 * (x1 * h1_a + h1_b)
    return h1


def qcalc(pt, a):
    """Calculate heat effect due to addition"""

    n1c = pt['n1c']
    n2c = pt['n2c']
    n1d = pt['n1d']
    n2d = pt['n2d']
    
    x2ci = n2c / (n1c + n2c)
    x1ci = 1.0 - x2ci
    x2cf = (n2c + n2d) / (n1c + n1d + n2c + n2d)
    x1cf = 1.0 - x2cf
    x2d = n2d / (n1d + n2d)
    x1d = 1.0 - x2d

    q = n1d * (h1rk(x1cf, a) - h1rk(x1d, a)) + \
        n1c * (h1rk(x1cf, a) - h1rk(x1ci, a)) + \
        n2d * (h2rk(x2cf, a) - h2rk(x2d, a)) + \
        n2c * (h2rk(x2cf, a) - h2rk(x2ci, a))
    return q

# --------------------------------------

def qfit(params, ict, q):
    """function to minimize"""

    a = 5 * [0.0 ]
    a[0] = params['a0'].value
    a[1] = params['a1'].value
    a[2] = params['a2'].value
    a[3] = params['a3'].value
    a[4] = params['a4'].value

    npts = len(q)
    qc = np.zeros(npts)
    for i in range(npts):
        qc[i] = qcalc(ict[i], a)
    return qc - q

# --------------------------------------
        
def main():
    parser = argparse.ArgumentParser(description = \
        'Fit Redlich-Kister polynomial to titration calorimetry experiment\n'
        'and calculate partial molar enthaplies and enthalpy of mixing.\n'
        'Input file with same format as data tables in\n'
        'E. Matteoli, L. Lepori, Fluid Phase Eq. 174 (2000) 115-131:\n'
        'n1c/mol n2c/mol n1d/mmol n2s/mmol Q/J\n...',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--degree', type = int, default = 4,
                        help = 'degree of RK polynomial (default 4)')
    parser.add_argument('--plot', action = 'store_true',
                        help = 'plot results using matplotlib')
    parser.add_argument('itcfile', help = 'file with ITC data')
    args = parser.parse_args()

    if args.degree < 0 or args.degree > 4:
        print('choose RK degree up to 4')
        sys.exit(0)
    plotlib = True
    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except:
            plotlib = False
            print('matplotlib not imported, skipping plots')

    fname = args.itcfile.split('.')[0]

    itc = readitc(args.itcfile)
    npts = len(itc)
    x2 = np.zeros(npts)
    q = np.zeros(npts)
    for i in range(npts):
        pt = itc[i]
        x2[i] = pt['n2c'] / (pt['n1c'] + pt['n2c'])
        q[i] = pt['q']

    vary = 5 * [ True ]
    for i in range(args.degree + 1, 5):
        vary[i] = False
        
    params = Parameters()
    params.add('a0', value = 0.0, vary = vary[0])
    params.add('a1', value = 0.0, vary = vary[1])
    params.add('a2', value = 0.0, vary = vary[2])
    params.add('a3', value = 0.0, vary = vary[3])
    params.add('a4', value = 0.0, vary = vary[4])

    result = minimize(qfit, params, args = (itc, q))

    print(fit_report(result))

    a = 5 * [ 0.0 ]
    a[0] = result.params['a0'].value
    a[1] = result.params['a1'].value
    a[2] = result.params['a2'].value
    a[3] = result.params['a3'].value
    a[4] = result.params['a4'].value

    qc = np.zeros(npts)
    for i in range(npts):
        qc[i] = qcalc(itc[i], a)

    qfile = fname + '_q.out'
    with open(qfile, 'w') as f:
        f.write('#     x2        Qexp/J       Qcalc/J\n')
        for i in range(npts):
            f.write('{0:5d} {1:8.6f} {2:12.5e} {3:12.5e}\n'.format(i, x2[i],
                  q[i], qc[i]))
    print('calculated heats:\n  {0}'.format(qfile))

    if args.plot and plotlib:
        plt.axis([0, npts, -0.1, 0.1])
        plt.plot([0, npts], [0.0, 0.0], 'k--')
        plt.plot(qc - q, 'r+')
        plt.xlabel('data')
        plt.ylabel('Q (J)')
        plt.show()

    x1exp = []
    h1exp = []
    x2exp = []
    h2exp = []
    for i in range(npts):
        pt = itc[i]
        if pt['n1d'] > 0.0 and pt['n2d'] == 0.0:
            x1exp.append(x2[i])
            h1exp.append(pt['q'] / pt['n1d'])
        if pt['n2d'] > 0.0 and pt['n1d'] == 0.0:
            x2exp.append(x2[i])
            h2exp.append(pt['q'] / pt['n2d'])
    h1file = fname + '_h1.out'
    h2file = fname + '_h2.out'
    with open(h1file, 'w') as f:
        f.write('# x2      h1exp/(J/mol)\n')
        for i in range(len(h1exp)):
            f.write('{0:8.6f} {1:12.5e}\n'.format(x1exp[i], h1exp[i]))
    with open(h2file, 'w') as f:
        f.write('# x2      h2exp/(J/mol)\n')
        for i in range(len(h2exp)):
            f.write('{0:8.6f} {1:12.5e}\n'.format(x2exp[i], h2exp[i]))
    print('experimental partial molar h:\n  {0}\n  {1}'.format(h1file,
                                                               h2file))

    # fitted parameters from Matteoli
    # a = [608.72, 3954.6, -950.93, 3618.5, -1120.9]

    nx = 101
    x2rk = np.linspace(0.0, 1.0, nx)
    h1fit = np.zeros(nx)
    h2fit = np.zeros(nx)
    he = np.zeros(nx)
    for i in range(nx):
        h1fit[i] = h1rk(1.0 - x2rk[i], a)
        h2fit[i] = h2rk(x2rk[i], a)
        he[i] = herk2(x2rk[i], a)
    hrkfile = fname + '_hrk.out'
    with open(hrkfile, 'w') as f:
        f.write('# x2      h1/(J/mol)   h2/(J/mol)   he/(J/mol)\n')
        for i in range(101):
            f.write('{0:8.6f} {1:12.5e} {2:12.5e} {3:12.5e}\n'.format(
                    x2rk[i], h1fit[i], h2fit[i], he[i]))
    print('calculated partial molar and excess H:\n  {0}'.format(hrkfile))

    if args.plot and plotlib:
        plt.axis([0, 1.0, -11000.0, 7000.0])
        plt.plot([0.0, 1.0], [0.0, 0.0], 'k--')
        plt.plot(x1exp, h1exp, 'ro')
        plt.plot(x2exp, h2exp, 'bo')
        plt.plot(x2rk, h1fit, 'r')
        plt.plot(x2rk, h2fit, 'b')
        plt.plot(x2rk, he, 'k')
        plt.xlabel('x2')
        plt.ylabel('H / (kJ/mol)')
        plt.show()


if __name__ == '__main__':
    main()
