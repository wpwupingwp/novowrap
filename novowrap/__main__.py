#!/usr/bin/python3

from sys import argv
from novowrap.assembly import assembly_main
from novowrap.ui import ui_main


def main():
    if argv[-1] == '-h':
        assembly_main()
    try:
        ui_main()
    except Exception:
        assembly_main()


if __name__ == '__main__':
    main()
