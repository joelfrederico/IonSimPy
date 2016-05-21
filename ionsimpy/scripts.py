import argparse
from .visualization import readfield as _readfield


def _readfield_script():
    
    parser = argparse.ArgumentParser()
    # parser.add_argument(['-f', 'filename'], dest='filename', default='output.h5')
    parser.add_argument('-f', dest='filename', default='output.h5')
    
    args = parser.parse_args()
    
    _readfield(args.filename)
