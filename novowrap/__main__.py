#!/usr/bin/python3

from sys import argv
from novowrap.assembly import assembly_main
from novowrap.ui import ui_main
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
    elif len(argv) > 1 and argv[-1] != '-h':
        get_all_third_party()
        argv.pop()
        assembly_main()
        return
    try:
        ui_main()
    except Exception:
        get_all_third_party()
        assembly_main()


if __name__ == '__main__':
    main()
