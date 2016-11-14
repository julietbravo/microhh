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


if __name__ == "__main__":
    """ Examples """

    # Read namelist file:
    nl = Read_namelist()
    print('Namelist values: itot={}, endtime={}'.format(nl.grid.itot, nl.time.endtime))

    # Read all statistics
    s = Read_statistics('patch.default.0000000.nc')
    print('Statistics heights: {}, ....'.format(s.z[:2]))
