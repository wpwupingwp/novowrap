#!/usr/bin/python3

from sys import argv
from novowrap.assembly import assembly_main
from novowrap.utils import get_all_third_party
import logging


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('novowrap')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger('novowrap')


def main():
    if argv[-1] == 'init':
        get_all_third_party()
        return
    elif argv[-1] == '-h':
        assembly_main()
        return
    elif len(argv) > 1:
        assembly_main()
        return
    try:
        from novowrap.ui import ui_main
        ui_main()
    except Exception:
        assembly_main()


if __name__ == '__main__':
    main()
