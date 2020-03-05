def main():
    from novowrap.assembly import assembly_main
    from novowrap.ui import ui_main
    try:
        ui_main()
    except Exception:
        assembly_main()


if __name__ == '__main__':
    main()
