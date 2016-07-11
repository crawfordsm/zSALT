
"""
GALEXTRACT

extract galaxy spectra from image or images and combine them

"""

from agnextract import salt_extract

if __name__=='__main__':
    import argparse 
    parser = argparse.ArgumentParser(description='Extract spectra from SALT 2D rectified image')
    parser.add_argument('objfile', help='SALT 2D rectified image')
    parser.add_argument('yc', help='Central row of object', type=int)
    parser.add_argument('dy', help='Half width of object', type=int)
    parser.add_argument('--spst', dest='cal_file', default='',
                   help='SPST calibration file')
    parser.add_argument('--f', dest='format', default='ascii', choices=['ascii','lcogt'],
                   help='Format for output file')
    parser.add_argument('--e', dest='ext', default=1, type=int,
                   help='Extension to extract')
    args = parser.parse_args()

   
    salt_extract(args.objfile, args.yc, args.dy, specformat=args.format, ext=args.ext, calfile=args.cal_file, convert=True, normalize=False)
