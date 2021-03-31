#!/usr/bin/python3

from logging import handlers
from pathlib import Path
import logging
import queue
import threading
try:
    from tkinter import messagebox, filedialog, scrolledtext
    import tkinter as tk
    from tkinter import ttk
except ImportError:
    print('Cannot find the module tkinter, please follow the README to install it.')
    print('Or try to use portable version https://github.com/wpwupingwp/novowrap/releases')
    print('Or, use conda instead (see README for details)')
    raise

from novowrap.utils import accessible
from novowrap.merge import merge_main
from novowrap.assembly import assembly_main
from novowrap.validate import validate_main


# define logger
# ui needn't coloredlogs
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('novowrap')


def scroll_text(window):
    """
    ScrolledText that shows logs.
    """
    def poll():
        while True:
            try:
                msg = log_queue.get(block=False)
                level = msg.levelname
                msg = formatter.format(msg) + '\n'
                scroll.insert('end', msg, level)
                scroll.yview('end')
            except queue.Empty:
                break
        # to avoid orphan poll()
        if log.hasHandlers():
            scroll.after(10, poll)
        else:
            return

    # clean old handlers
    for i in log.handlers:
        log.removeHandler(i)
    log_queue = queue.Queue()
    formatter = logging.Formatter(fmt=FMT, datefmt=DATEFMT)
    # do not add formatter to queuehandler, or msg will be formatted twice
    queue_handler = handlers.QueueHandler(log_queue)
    # give poll() time to quit
    root.after(100, log.addHandler(queue_handler))
    scroll = scrolledtext.ScrolledText(window)
    scroll.tag_config('INFO', foreground='black')
    scroll.tag_config('WARNING', foreground='orange')
    scroll.tag_config('ERROR', foreground='red')
    scroll.tag_config('CRITICAL', foreground='red')
    scroll.tag_config('EXCEPTION', foreground='red')
    scroll.pack(fill='both')
    scroll.after(0, poll)


def wlabel(window, text, row, column=0, width=25, padx=0, pady=0, sticky='ew',
           **kargs):
    """
    Generate and pack labels.
    """
    label = ttk.Label(window, text=text, width=width, anchor=tk.CENTER,
                      **kargs)
    label.grid(row=row, column=column, padx=padx, pady=pady, sticky=sticky)
    return label


def fentry(window, row, column, default='', padx=0, pady=0):
    """
    Generate and pack entrys.
    Fill with default string.
    """
    entry = ttk.Entry(window)
    entry.insert(0, default)
    entry.grid(row=row, column=column, padx=padx, pady=pady)
    return entry


def info(message):
    """
    For shorter.
    """
    messagebox.showinfo(message=message)


def after_close(frame):
    """
    Deiconify root before destroy.
    """
    def func():
        root.deiconify()
        frame.destroy()
    return func


def open_file(title, entry, single=True, entry2=None):
    """
    Set title, fill entry 1, empty entry 2.
    """
    def func():
        if single:
            a = filedialog.askopenfilename(title=title)
        else:
            a = filedialog.askopenfilenames(title=title)
        entry.delete(0, 'end')
        entry.insert(0, a)
        if entry2 is not None:
            entry2.delete(0, 'end')
    return func


def open_folder(title, entry):
    def func():
        a = filedialog.askdirectory(title=title)
        entry.delete(0, 'end')
        entry.insert(0, a)
    return func


def check_output(output, out_entry):
    """
    Make sure given path does not exist and writable.
    """
    output = Path(output).absolute()
    if output.exists():
        info('Output folder exists. Please use another folder.')
        out_entry.configure(style='hi.TButton')
        return False
    else:
        ok = accessible(output, 'folder')
        if not ok:
            info(f'You do not have permission to write in given output '
                 'folder. Please try another path.')
            out_entry.configure(style='hi.TButton')
            return False
    return True


def thread_wrap(function, arg_str, window):
    """
    Wrap for callback.
    The validate and merge share same return structure.
    Args:
        function(callable): function to call
        arg_str(str): string for fuction's argparse
        window(Toplevel): window to hide
    """
    try:
        result = function(arg_str)
    except Exception:
        log.exception('Abort.')
        info('Abort.')
        root.deiconify()
        return
    if result[0]:
        info(f'Done. See {result[1]} for details.')
    else:
        info(f'Fail. See {result[1]} for details.')
    window.withdraw()
    root.deiconify()
    return


def assembly_ui():
    """
    UI of assembly.
    """
    def show_adv():
        wroot.geometry(big_size)
        advance.grid(advance_info)

    def submit_assembly():
        # prepare arg_str
        arg_str = ''
        arg_input = input_entry.get()
        arg_input_list = list_entry.get()
        if arg_input and arg_input_list:
            info('Please only use one of "Input file" and "Input list"!')
            return
        elif arg_input == '' and arg_input_list == '':
            info('Input is required!')
            return
        # use dialog window to choose input, usually needn't check whether
        # exists
        elif arg_input:
            arg_str += f'-input {arg_input}'
        else:
            arg_str += f'-list {arg_input_list}'
        arg_ref = ref_entry.get()
        # use underscore in taxon name, which NCBI could handle
        arg_taxon = taxon_entry.get()
        if arg_ref and arg_taxon:
            info('Please only use one of "Genbank file" and "Taxonomy"!')
            return
        elif arg_ref == '' and len(arg_taxon) == 0:
            info('Reference is required!')
            return
        elif arg_ref != '' and len(arg_taxon) == 0:
            arg_str += f' -ref {arg_ref}'
        elif arg_ref == '' and len(arg_taxon) != 0:
            arg_str += f' -taxon {arg_taxon}'
        arg_out = out_entry.get()
        if arg_out == '"Current folder"':
            out_path = Path('.').absolute()
        else:
            out_path = Path(arg_out).absolute()
        if out_path.exists():
            if out_path.is_dir():
                out_path = out_path / 'Output'
            else:
                info('Invalid output folder name.')
                return
        if not check_output(out_path, out_entry):
            return
        arg_str += f' -out {out_path}'
        # advanced options
        arg_split = split_entry.get()
        if arg_split:
            arg_str += f' -split {arg_split}'
        arg_insert = insert_entry.get()
        if arg_insert:
            arg_insert += f' -insert_size {arg_insert}'
        arg_str += f' -platform {platform.get()}'
        arg_kmer = int(kmer_entry.get())
        if arg_kmer > 39 or arg_kmer < 23 or (arg_kmer % 2 != 1):
            info('K-mer should be an odd number in (23, 39)!')
            kmer_entry.configure(style='hi.TButton')
            return
        arg_str += f' -kmer {arg_kmer}'
        arg_size = size_entry.get()
        if '-' not in arg_size:
            info('Genome size should be "min-max" format!')
            size_entry.configure(style='hi.TButton')
            return
        min_size, max_size = arg_size.split('-')
        arg_str += f' -min {min_size} -max {max_size}'
        arg_seed = seed_entry.get().strip(',')
        if arg_seed:
            arg_str += f' -seed {arg_seed}'
        arg_seed_file = seed_file_entry.get()
        if arg_seed_file:
            arg_str += f' -seed_file {arg_seed_file}'
        # handle checkbutton
        if mt_mode.get():
            simple_validation.set(True)
            size_ok = messagebox.askyesno(
                title='Warning',
                message=(f'Given genome length range is {min_size}-{max_size}, '
                        f'ok for mitochondria?'))
            if not size_ok:
                return
            arg_str += ' -mt'
        else:
            if simple_validation.get():
                arg_str += ' -simple_validate'
        if simple_validation.get():
            info('The program will skip detecting quadripartite structure and '
                 'adjustment.')
        arg_s = float(s_entry.get())
        arg_l = float(l_entry.get())
        if max(arg_s, arg_l) > 1 or min(arg_s, arg_l) <= 0:
            info('Bad value!')
            s_entry.configure(style='hi.TButton')
            l_entry.configure(style='hi.TButton')
            return
        else:
            arg_str += f' -len_diff {arg_l} -perc_identity {arg_s}'
        # info(arg_str)
        # call validate
        wroot.withdraw()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(assembly_main, arg_str, run),
                             daemon=True)
        r.start()

    root.iconify()
    wroot = tk.Toplevel(root)
    wroot.geometry(size)
    wroot.title('Assembly')
    w = ttk.Frame(wroot)
    w.place(relx=0.5, rely=0.5, anchor='center')
    # use variable for easily edit
    row = 0
    inputs = ttk.LabelFrame(w, text='Input')
    inputs.grid(row=row, columnspan=4)
    row += 1
    wlabel(inputs, 'Input', row=row, column=1)
    input_entry = fentry(inputs, row=row, column=2)
    input_button = ttk.Button(inputs, text='Open', command=open_file(
        'Input file', input_entry, single=False))
    input_button.grid(row=row, column=3)
    row += 1
    # ttk label do not support 'fg'
    label1 = tk.Label(inputs, text='OR', fg='red')
    label1.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(inputs, 'Input list', row=row, column=1)
    list_entry = fentry(inputs, row=row, column=2)
    list_button = ttk.Button(inputs, text='Open',
                             command=open_file('List file', list_entry,
                                               entry2=input_entry))
    list_button.grid(row=row, column=3)
    row += 1
    ref = ttk.LabelFrame(w, text='Reference')
    ref.grid(row=row, columnspan=4)
    row += 1
    wlabel(ref, 'Taxonomy', row=row, column=1)
    taxon_entry = fentry(ref, row=row, column=2, default='Nicotiana tabacum')
    row += 1
    label2 = tk.Label(ref, text='OR', fg='red')
    label2.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(ref, 'Genbank file', row=row, column=1)
    ref_entry = fentry(ref, row=row, column=2)
    r_button = ttk.Button(ref, text='Open',
                          command=open_file('Reference file', ref_entry,
                                            entry2=taxon_entry))
    r_button.grid(row=row, column=3)
    row += 1
    wlabel(w, 'Output', row=row)
    out_entry = fentry(w, row=row, column=1, default=str(Path('.').absolute()))
    o_button = ttk.Button(w, text='Open',
                          command=open_folder('Output folder', out_entry))
    o_button.grid(row=row, column=2)
    row += 1

    advance = ttk.LabelFrame(w, text='Advanced settings')
    advance.grid(row=row, columnspan=3)
    advance_info = advance.grid_info()
    row += 1
    adv_input = ttk.LabelFrame(advance, text='Input')
    adv_input.grid(row=row)
    row += 1
    simple_validation = tk.BooleanVar()
    simple_validation.set(False)
    check1 = ttk.Checkbutton(adv_input, text='Simple validation',
                             variable=simple_validation, onvalue=True,
                             offvalue=False)
    check1.grid(row=row, column=1, sticky='w')
    row += 1
    mt_mode = tk.BooleanVar()
    mt_mode.set(False)
    check2 = ttk.Checkbutton(adv_input, text='Mitochondria genome',
                             variable=mt_mode, onvalue=True, offvalue=False)
    check2.grid(row=row, column=1, sticky='w')
    row += 1
    wlabel(adv_input, 'Split reads', row=row, column=0)
    split_entry = fentry(adv_input, row=row, column=1)
    row += 1
    wlabel(adv_input, 'Insert size', row=row, column=0)
    insert_entry = fentry(adv_input, row=row, column=1)
    row += 1
    platform = tk.StringVar()
    platform.set('illumina')
    wlabel(adv_input, 'Sequence platform', row=row, column=0)
    radio1 = ttk.Radiobutton(adv_input, text='illumina', variable=platform,
                             value='illumina')
    radio1.grid(row=row, column=1, sticky='w')
    radio2 = ttk.Radiobutton(adv_input, text='ion torrent', variable=platform,
                             value='ion')
    radio2.grid(row=row, column=2, sticky='w')
    row += 1
    adv_assembly = ttk.LabelFrame(advance, text='Assembly')
    adv_assembly.grid(row=row)
    row += 1
    wlabel(adv_assembly, 'K-mer (23-39, odd)', row=row, column=0)
    kmer_entry = fentry(adv_assembly, row=row, column=1, default='39')
    row += 1
    wlabel(adv_assembly, 'Genome Size range (bp)', row=row, column=0)
    size_entry = fentry(adv_assembly, row=row, column=1,
                        default='100000-200000')
    row += 1
    wlabel(adv_assembly, 'Seed genes', row=row, column=0)
    seed_entry = fentry(adv_assembly, row=row, column=1,
                        default='rbcL,psaB,psaC,rrn23')
    row += 1
    wlabel(adv_assembly, 'Seed file', row=row, column=0)
    seed_file_entry = fentry(adv_assembly, row=row, column=1)
    seed_button = ttk.Button(adv_assembly, text='Open',
                             command=open_file('Seed file', seed_file_entry,
                                               single=True))
    seed_button.grid(row=row, column=2)

    adv_validate = ttk.LabelFrame(advance, text='Validate')
    adv_validate.grid(row=row)
    row += 1
    wlabel(adv_validate, 'Sequence similarity (0-1)', row=row, column=0)
    s_entry = fentry(adv_validate, row=row, column=1, default='0.7')
    row += 1
    wlabel(adv_validate, 'Length difference (0-1)', row=row, column=0)
    l_entry = fentry(adv_validate, row=row, column=1, default='0.2')

    row += 1
    show_more = ttk.Button(w, text='More options', command=show_adv)
    show_more.grid(row=row, column=0, pady=10)
    ok = ttk.Button(w, text='Enter', command=submit_assembly)
    ok.grid(row=row, column=1, columnspan=2, sticky='EW', padx=40)
    advance.grid_forget()
    wroot.protocol('WM_DELETE_WINDOW', after_close(wroot))
    return


def merge_ui():
    """
    UI of merge.
    """
    def submit_merge():
        arg_str = ''
        inputs = input_entry.get()
        if inputs == '':
            info('Input is required!')
            return
        arg_str += f'-input {inputs}'
        arg_out = out_entry.get()
        if arg_out == '"Current folder"':
            out_path = Path('.').absolute()
        else:
            out_path = Path(arg_out).absolute()
        if out_path.exists():
            if out_path.is_dir():
                out_path = out_path / 'Output'
            else:
                info('Invalid output folder name.')
                return
        if not check_output(out_path, out_entry):
            return
        arg_str += f' -out {out_path}'
        print(arg_str)
        wroot.withdraw()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(merge_main, arg_str, wroot),
                             daemon=True)
        r.start()

    root.iconify()
    wroot = tk.Toplevel(root)
    wroot.geometry(small_size)
    wroot.title('Merge')
    frame = ttk.Frame(wroot)
    frame.place(relx=0.5, rely=0.5, anchor='center')
    row = 0
    wlabel(frame, 'Input', row=row, column=1)
    input_entry = fentry(frame, row=row, column=2)
    input_button = ttk.Button(frame, text='Open',
                              command=open_file('Input file', input_entry,
                                                single=False))
    input_button.grid(row=row, column=3)
    row += 1
    wlabel(frame, 'Output', row=row, column=1, pady=10)
    out_entry = fentry(frame, row=row, column=2, default='"Current folder"')
    o_button = ttk.Button(frame, text='Open',
                          command=open_folder('Output folder', out_entry))
    o_button.grid(row=row, column=3)
    row += 1
    ok = ttk.Button(frame, text='Enter', command=submit_merge)
    ok.grid(row=row, column=1, columnspan=3, sticky='EW', pady=10)
    wroot.protocol('WM_DELETE_WINDOW', after_close(wroot))
    return


def validate_ui():
    """
    UI of validate.
    """
    def submit_validate():
        # prepare arg_str
        arg_str = ''
        inputs = i_entry.get()
        if inputs == '':
            info('Input is required!')
            return
        else:
            arg_str += f'-input {inputs}'
        arg_ref = r_entry.get()
        # use underscore in taxon name, which NCBI could handle
        arg_taxon = t_entry.get().replace(' ', '_')
        if arg_ref and arg_taxon:
            info('Please only use one of "Reference file" and "Taxonomy"!')
            return
        elif arg_ref != '' and arg_taxon == '':
            arg_str += f' -ref {arg_ref}'
        elif arg_ref == '' and arg_taxon != '':
            arg_str += f' -taxon {arg_taxon}'
        arg_out = o_entry.get()
        if arg_out == '"Current folder"':
            out_path = Path('.').absolute()
        else:
            out_path = Path(arg_out).absolute()
        if out_path.exists():
            if out_path.is_dir():
                out_path = out_path / 'Output'
            else:
                info('Invalid output folder name.')
                return
        if not check_output(out_path, o_entry):
            return
        if mt_mode.get():
            simple_validation.set(True)
            arg_str += ' -mt'
        else:
            if simple_validation.get():
                arg_str += ' -simple_validate'
        if simple_validation.get():
            info('The program will skip detecting quadripartite structure and '
                 'adjustment.')
        # validate need existing out if called by others
        out_path.mkdir()
        out_tmp = out_path / 'Temp'
        out_tmp.mkdir()
        arg_str += f' -out {out_path}'
        arg_s = float(s_entry.get())
        arg_l = float(l_entry.get())
        if max(arg_s, arg_l) > 1 or min(arg_s, arg_l) <= 0:
            info('Bad value!')
            s_entry.configure(style='hi.TButton')
            l_entry.configure(style='hi.TButton')
            return
        else:
            arg_str += f' -len_diff {arg_l} -perc_identity {arg_s}'
        # call validate
        wroot.withdraw()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = ttk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(validate_main, arg_str, run),
                             daemon=True)
        r.start()

    root.iconify()
    wroot = tk.Toplevel(root)
    wroot.geometry(size)
    wroot.title('Validate')
    # on top
    # wroot.wm_transient(root)
    w = ttk.Frame(wroot)
    w.place(relx=0.5, rely=0.5, anchor='center')
    # use variable for easily edit
    row = 0
    wlabel(w, 'Input', row=row, padx=15, pady=10)
    i_entry = fentry(w, row=row, column=1)
    i_button = ttk.Button(w, text='Open',
                          command=open_file('Input file', i_entry))
    i_button.grid(row=row, column=2)
    ref = ttk.LabelFrame(w, text='Reference')
    row += 1
    ref.grid(row=row, column=0, columnspan=5, padx=5, pady=8)
    row += 1
    wlabel(ref, 'Taxonomy', row=row, column=1)
    t_entry = fentry(ref, row=row, column=2, default='Nicotiana tabacum')
    row += 1
    label1 = tk.Label(ref, text='OR', fg='red')
    label1.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(ref, 'Genbank/FASTA file', row=row, column=1, padx=13)
    r_entry = fentry(ref, row=row, column=2)
    r_button = ttk.Button(ref, text='Open',
                          command=open_file('Reference file (gb/fasta format)',
                                            r_entry, entry2=t_entry))
    r_button.grid(row=row, column=3)
    row += 1
    wlabel(w, 'Output', row=row, padx=10, pady=8)
    o_entry = fentry(w, row=row, column=1, default=str(Path('.').absolute()))
    o_button = ttk.Button(w, text='Open', command=open_folder('Output folder',
                                                              o_entry))
    o_button.grid(row=row, column=2)
    row += 1
    options = ttk.LabelFrame(w, text='Options')
    options.grid(row=row, padx=5, columnspan=5)
    row += 1
    simple_validation = tk.BooleanVar()
    simple_validation.set(False)
    check1 = ttk.Checkbutton(options, text='Simple validation',
                             variable=simple_validation, onvalue=True,
                             offvalue=False)
    check1.grid(row=row, column=1, sticky='w')
    row += 1
    mt_mode = tk.BooleanVar()
    mt_mode.set(False)
    check2 = ttk.Checkbutton(options, text='Mitochondria genome',
                             variable=mt_mode, onvalue=True, offvalue=False)
    check2.grid(row=row, column=1, sticky='w')
    row += 1
    wlabel(options, 'Sequence similarity (0-1)', row=row, padx=10)
    s_entry = fentry(options, row=row, column=1, default='0.7')
    row += 1
    wlabel(options, 'Length difference (0-1)', row=row, padx=10)
    l_entry = fentry(options, row=row, column=1, default='0.2')
    row += 1
    ok = ttk.Button(w, text='Enter', command=submit_validate)
    ok.grid(row=row, column=0, columnspan=3, sticky='EW', padx=50, pady=10)
    wroot.protocol('WM_DELETE_WINDOW', after_close(wroot))
    return


def ui_main():
    # init window
    global root
    root = tk.Tk()
    root.attributes('-topmost', 'true')
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    s = min(w, h) // 2
    global size
    size = f'{s}x{int(s*0.618)}+{w//3}+{h//3}'
    global small_size
    small_size = f'{s}x{int(s*0.618/2)}+{w//3}+{h//3}'
    global big_size
    big_size = f'{s}x{int(s*0.618*2)}'
    style = ttk.Style()
    # don't know why 'hi.TEntry' can't work
    style.configure('hi.TButton', background='red')
    root.geometry(small_size)
    root.title('novowrap')
    # 1366x768
    if h < 800:
        root.tk.call('tk', 'scaling', 0.9)
    # 1440x900
    elif h < 1000:
        root.tk.call('tk', 'scaling', 0.95)
    # 2k
    elif h > 1400:
        root.tk.call('tk', 'scaling', 2.0)
    # cross-platform font?
    # default_font = font.nametofont('TkDefaultFont')
    # default_font_cfg = default_font.configure()
    # default_font.config(size=12)
    root_frame = ttk.Frame(root)
    # in center
    root_frame.place(relx=0.5, rely=0.5, anchor='center')
    # high dpi
    # root.tk.call('tk', 'scaling', 2.0)
    # btn->button
    assembly_btn = ttk.Button(root_frame, text='Assembly sequences',
                              command=assembly_ui)
    assembly_btn.grid(row=0, pady=10)
    validate_btn = ttk.Button(root_frame, text='Validate assemblies',
                              command=validate_ui)
    validate_btn.grid(row=1, pady=10)
    merge_button = ttk.Button(root_frame, text='Merge contigs',
                              command=merge_ui)
    merge_button.grid(row=2, pady=10)
    root.mainloop()


if __name__ == '__main__':
    ui_main()
