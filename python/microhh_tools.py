import netCDF4 as nc4
import numpy   as np
import struct  as st
import matplotlib.pylab as pl

import glob
import fileinput
import re

class _Empty: pass

def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if ('.' in value or 'e' in value):
            return float(value)
        else:
            return int(value)
    except:
        return value.rstrip()


def _convert_value(value):
    """ Helper function: convert namelist value or list """
    if ',' in value:
        value = value.split(',')
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)


def _find_namelist_file():
    """ Helper function: automatically find the .ini file in the current directory """
    namelist_file = glob.glob('*.ini')
    if len(namelist_file) == 0:
        raise RuntimeError('Can\'t find any .ini files in the current directory!')
    if len(namelist_file) > 1:
        raise RuntimeError('There are multiple .ini files: {}'.format(namelist_file))
    else:
        return namelist_file[0]


class Read_namelist:
    """ Reads a MicroHH .ini file to memory """
    def __init__(self, namelist_file=None):
        if (namelist_file is None):
            namelist_file = _find_namelist_file()

        curr_group = None
        with open(namelist_file) as f:
            for line in f:
                if len(line.strip()) > 0:
                    if line.strip()[0] == '[' and line.strip()[-1] == ']':
                        curr_group_name = line.strip()[1:-1]
                        curr_group = _Empty()
                        setattr(self, curr_group_name, curr_group)
                    elif curr_group is not None:
                        setattr(curr_group, line.split('=')[0], _convert_value(line.split('=')[-1]))


def replace_namelist_var(variable, new_value, namelist_file=None):
    """ Replace a variable value in an existing namelist """
    if namelist_file is None:
        namelist_file = _find_namelist_file()

    with open(namelist_file, "r") as source:
        lines = source.readlines()
    with open(namelist_file, "w") as source:
        for line in lines:
            source.write(re.sub(r'({}).*'.format(variable), r'\1={}'.format(new_value), line))


class Read_statistics:
    """ Read all statistics to memory """
    def __init__(self, stat_file):
        f = nc4.Dataset(stat_file, 'r')

        self.data       = {}
        self.units      = {}
        self.names      = {}
        self.dimensions = {}

        for var in f.variables:
            self.data[var]       = f.variables[var].__array__()
            self.units[var]      = f.variables[var].units
            self.names[var]      = f.variables[var].long_name
            self.dimensions[var] = f.variables[var].dimensions

        f.close()

    def __getitem__(self, name):
        """ Retrieve data arrays using notation s['th'] """
        if name in self.data:
            return self.data[name]
        else:
            raise RuntimeError('Can\'t find variable {} in statistics file'.format(name))

    def __getattr__(self, name):
        """ Retrieve data arrays using notation s.th """
        if name in self.data:
            return self.data[name]


def get_cross_indices(variable, mode):
    """ Get the cross-section indices """
    files = glob.glob('{}.{}.*.*'.format(variable, mode))
    time  = files[0].split('.')[-1]
    files = glob.glob('{}.{}.*.{}'.format(variable, mode, time))
    indices = [int(f.split('.')[-2]) for f in files]
    indices.sort()
    return indices

# Thermodynamics, as in MicroHH
T0   = 273.15
Rd   = 287.04
Rv   = 461.5
cp   = 1005.
ep   = Rd/Rv

c00  = +6.1121000000E+02;
c10  = +4.4393067270E+01;
c20  = +1.4279398448E+00;
c30  = +2.6415206946E-02;
c40  = +3.0291749160E-04;
c50  = +2.1159987257E-06;
c60  = +7.5015702516E-09;
c70  = -1.5604873363E-12;
c80  = -9.9726710231E-14;
c90  = -4.8165754883E-17;
c100 = +1.3839187032E-18;

def esat(T):
    x = np.maximum(T-T0, -75.)
    return c00+x*(c10+x*(c20+x*(c30+x*(c40+x*(c50+x*(c60+x*(c70+x*(c80+x*(c90+x*c100)))))))))

def qsat(p, T):
    return ep*esat(T)/(p-(1-ep)*esat(T))

def exner(p, p0=1e5):
    return pow((p/p0),(Rd/cp))

if __name__ == "__main__":
    """ Examples """

    # Read namelist file:
    nl = Read_namelist()
    print('Namelist values: itot={}, endtime={}'.format(nl.grid.itot, nl.time.endtime))

    # Read all statistics
    s = Read_statistics('patch.default.0000000.nc')
    print('Statistics heights: {}, ....'.format(s.z[:2]))

