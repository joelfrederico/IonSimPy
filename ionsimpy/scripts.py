import argparse as _argparse
from .visualization import readfield as _readfield
from matplotlib import pyplot as _plt


def _readfield_script():
    
    parser = _argparse.ArgumentParser()
    parser.add_argument('filename', nargs='?', default='output.h5')
    
    args = parser.parse_args()
    
    _readfield(args.filename)
