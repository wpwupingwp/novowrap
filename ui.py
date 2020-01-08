#!/usr/bin/python3

import logging
from logging import handlers
import threading
import queue
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog, scrolledtext

from pathlib import Path

from merge import merge_main
from validate import validate_main
from novowrap import assembly_main


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
                msg = formatter.format(msg) + '\n'
                scroll.insert('end', msg)
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
    scroll.pack(fill='both')
    scroll.after(0, poll)


def wlabel(window, text, row, column=0, width=25, padx=0, pady=0, sticky='EW',
           **kargs):
    """
    Generate and pack labels.
    """
    label = tk.Label(window, text=text, width=width, **kargs)
    label.grid(row=row, column=column, padx=padx, pady=pady, sticky=sticky)
    return label


def fentry(window, row, column, default='', padx=0, pady=0):
    """
    Generate and pack entrys.
    Fill with default string.
    """
    entry = tk.Entry(window)
    entry.insert(0, default)
    entry.grid(row=row, column=column, padx=padx, pady=pady)
    return entry


def info(message):
    """
    For shorter.
    """
    messagebox.showinfo(message=message)


def open_file(title, entry, single=True, entry2=None):
    """
    Set title, fill entry 1, empty entry 2.
    """
    def func():
        if single:
            a = filedialog.askopenfilename(title=title)
        else:
            a = filedialog.askopenfilenames(title=title)
        entry.insert(0, a)
        if entry2 is not None:
            entry2.delete(0, 'end')
    return func


def open_folder(title, entry):
    def func():
        a = filedialog.askdirectory(title=title)
        entry.insert(0, a)
    return func


def thread_wrap(function, arg_str, window):
    """
    Wrap for callback.
    The validate and merge share same return structure.
    Args:
        function(callable): function to call
        arg_str(str): string for fuction's argparse
        window(Toplevel): window to hide
    """
    result = function(arg_str)
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
        elif arg_input:
            arg_str += f'-i {arg_input}'
        else:
            arg_str += f'-l {arg_input_list}'
        arg_ref = ref_entry.get()
        # use underscore in taxon name, which NCBI could handle
        arg_taxon = taxon_entry.get().replace(' ', '_')
        if arg_ref and arg_taxon:
            info('Please only use one of "Refence" and "Taxonomy"!')
            return
        elif arg_ref != '' and arg_taxon == '':
            arg_str += f' -ref {arg_ref}'
        elif arg_ref == '' and arg_taxon != '':
            arg_str += f' -taxon {arg_taxon}'
        arg_out = out_entry.get()
        if arg_out == '"Current folder"':
            arg_out = str(Path('.').absolute())
        else:
            arg_out = str(Path(arg_out).absolute())
        arg_str += f' -out {arg_out}'
        # advanced options
        arg_split = split_entry.get()
        if arg_split:
            arg_str += f' -split {arg_split}'
        arg_insert = insert_entry.get()
        if arg_insert:
            arg_insert += f' -insert_size {arg_insert}'
        arg_str += f' -p {platform.get()}'
        arg_kmer = int(kmer_entry.get())
        if arg_kmer > 39 or arg_kmer < 23 or (arg_kmer % 2 != 1):
            info('K-mer should be an odd number in (23, 39)!')
            return
        arg_str += f' -kmer {arg_kmer}'
        arg_size = size_entry.get()
        if '-' not in arg_size:
            info('Genome size should be "min-max" format!')
            return
        min_size, max_size = arg_size.split('-')
        arg_str += f' -min {min_size} -max {max_size}'
        arg_seed = seed_entry.get().strip(',')
        if arg_seed:
            arg_str += f' -seed {arg_seed}'
        arg_seed_file = seed_file_entry.get()
        if arg_seed_file:
            arg_str += f' -seed_file {arg_seed_file}'
        arg_s = float(s_entry.get())
        arg_l = float(l_entry.get())
        if max(arg_s, arg_l) > 1 or min(arg_s, arg_l) <= 0:
            info('Bad value!')
            s_entry.configure(bg='red')
            l_entry.configure(bg='red')
            return
        else:
            arg_str += f' -len_diff {arg_l} -perc_identity {arg_s}'
        # call validate
        wroot.withdraw()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = tk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        # to be continued
        r = threading.Thread(target=thread_wrap,
                             args=(assembly_main, arg_str, run))
        r.start()

    root.iconify()
    wroot = tk.Toplevel(root)
    wroot.geometry(size)
    wroot.title('Assembly')
    w = tk.Frame(wroot)
    w.grid(row=0, sticky='WENS', padx=30)
    # use variable for easily edit
    row = 0
    inputs = tk.LabelFrame(w, text='Input')
    inputs.grid(row=row, columnspan=4)
    row += 1
    wlabel(inputs, 'Input', row=row, column=1)
    input_entry = fentry(inputs, row=row, column=2)
    input_button = tk.Button(inputs, text='Open', command=open_file(
        'Input file', input_entry, single=False))
    input_button.grid(row=row, column=3)
    row += 1
    label1 = tk.Label(inputs, text='OR', fg='red')
    label1.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(inputs, 'Input list', row=row, column=1)
    list_entry = fentry(inputs, row=row, column=2)
    list_button = tk.Button(inputs, text='Open',
                            command=open_file('List file', list_entry,
                                              input_entry))
    list_button.grid(row=row, column=3)
    row += 1
    ref = tk.LabelFrame(w, text='Reference')
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
    r_button = tk.Button(ref, text='Open',
                         command=open_file('Reference file', ref_entry,
                                           taxon_entry))
    r_button.grid(row=row, column=3)
    row += 1
    wlabel(w, 'Output', row=row, pady=8)
    out_entry = fentry(w, row=row, column=1, default='"Current folder"')
    o_button = tk.Button(w, text='Open', command=open_folder('Output folder',
                                                             out_entry))
    o_button.grid(row=row, column=2)
    row += 1

    advance = tk.LabelFrame(w, text='Advanced settings')
    advance.grid(row=row, columnspan=3)
    advance_info = advance.grid_info()
    row += 1
    adv_input = tk.LabelFrame(advance, text='Input')
    adv_input.grid(row=row)
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
    radio1 = tk.Radiobutton(adv_input, text='illumina', variable=platform,
                            value='illumina')
    radio1.grid(row=row, column=1)
    radio2 = tk.Radiobutton(adv_input, text='ion torrent', variable=platform,
                            value='ion')
    radio2.grid(row=row, column=2)
    row += 1
    adv_assembly = tk.LabelFrame(advance, text='Assembly')
    adv_assembly.grid(row=row, padx=5)
    row += 1
    wlabel(adv_assembly, 'K-mer (23-39, odd)', row=row, column=0)
    kmer_entry = fentry(adv_assembly, row=row, column=1, default=39)
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
    seed_button = tk.Button(adv_assembly, text='Open',
                            command=open_file('Seed file', seed_file_entry,
                                              single=True))
    seed_button.grid(row=row, column=2)

    adv_validate = tk.LabelFrame(advance, text='Validate')
    adv_validate.grid(row=row)
    wlabel(adv_validate, 'Sequence similarity (0-1)', row=row, column=0)
    s_entry = fentry(adv_validate, row=row, column=1, default='0.7')
    row += 1
    wlabel(adv_validate, 'Length difference (0-1)', row=row)
    l_entry = fentry(adv_validate, row=row, column=1, default='0.2')

    row += 1
    show_more = tk.Button(w, text='More options', command=show_adv)
    show_more.grid(row=row, column=0, pady=10)
    ok = tk.Button(w, text='Enter', command=submit_assembly)
    ok.grid(row=row, column=1, columnspan=2, sticky='EW', padx=40)
    advance.grid_forget()
    return


def merge_ui():
    """
    UI of merge.
    """
    root.iconify()
    inputs = filedialog.askopenfilenames()
    wroot = tk.Toplevel(root)
    wroot.geometry(small_size)
    wroot.title('Merge')
    frame = tk.Frame(wroot)
    frame.pack(fill='both')
    scroll_text(frame)
    arg_str = ' '.join(inputs)
    r = threading.Thread(target=thread_wrap, args=(merge_main, arg_str, wroot))
    r.start()
    return


def validate_ui():
    """
    UI of validate.
    """
    def submit_validate():
        # prepare arg_str
        arg_str = ''
        arg_input = i_entry.get()
        if arg_input == '':
            info('Input is required!')
            return
        else:
            arg_str += arg_input
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
            arg_out = str(Path('.').absolute())
        arg_str += f' -out {arg_out}'
        arg_s = float(s_entry.get())
        arg_l = float(l_entry.get())
        if max(arg_s, arg_l) > 1 or min(arg_s, arg_l) <= 0:
            info('Bad value!')
            s_entry.configure(bg='red')
            l_entry.configure(bg='red')
            return
        else:
            arg_str += f' -len_diff {arg_l} -perc_identity {arg_s}'
        # call validate
        wroot.withdraw()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = tk.Frame(run)
        frame.pack(fill='both')
        scroll_text(frame)
        r = threading.Thread(target=thread_wrap,
                             args=(validate_main, arg_str, run))
        r.start()

    root.iconify()
    wroot = tk.Toplevel(root)
    wroot.geometry(size)
    wroot.title('Validate')
    # on top
    # wroot.wm_transient(root)
    w = tk.Frame(wroot)
    w.grid(row=0, sticky='WENS', padx=30)
    # use variable for easily edit
    row = 0
    wlabel(w, 'Input', row=row, padx=15, pady=10)
    i_entry = fentry(w, row=row, column=1)
    i_button = tk.Button(w, text='Open', command=open_file('Input file',
                                                           i_entry))
    i_button.grid(row=row, column=2)
    ref = tk.LabelFrame(w, text='Reference')
    row += 1
    ref.grid(row=row, column=0, columnspan=5, padx=5, pady=8)
    row += 1
    wlabel(ref, 'Taxonomy', row=row, column=1)
    t_entry = fentry(ref, row=row, column=2, default='Nicotiana tabacum')
    row += 1
    label1 = tk.Label(ref, text='OR', fg='red')
    label1.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(ref, 'File', row=row, column=1, padx=13)
    r_entry = fentry(ref, row=row, column=2)
    r_button = tk.Button(ref, text='Open',
                         command=open_file('Reference file (gb/fasta format)',
                                           r_entry, t_entry))
    r_button.grid(row=row, column=3)
    row += 1
    wlabel(w, 'Output', row=row, padx=10, pady=8)
    o_entry = fentry(w, row=row, column=1, default='"Current folder"')
    o_button = tk.Button(w, text='Open', command=open_folder('Output folder',
                                                             o_entry))
    o_button.grid(row=row, column=2)
    row += 1
    options = tk.LabelFrame(w, text='Options')
    options.grid(row=row, padx=5, columnspan=5)
    row += 1
    wlabel(options, 'Sequence similarity (0-1)', row=row, padx=10)
    s_entry = fentry(options, row=row, column=1, default='0.7')
    row += 1
    wlabel(options, 'Length difference (0-1)', row=row, padx=10)
    l_entry = fentry(options, row=row, column=1, default='0.2')
    row += 1
    ok = tk.Button(w, text='Enter', command=submit_validate)
    ok.grid(row=row, column=0, columnspan=3, sticky='EW', padx=50, pady=10)
    return


def o():
    pass
    # d = simpledialog.SimpleDialog(
    #     dialog, title='dialog test', text='this is a dialog',
    #     buttons=['yes', 'no', 'wo bu zhidao'], cancel=3, default=1)
    # value = d.go()
    # print(value)


# init window
root = tk.Tk()
root.attributes('-topmost', 'true')
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
s = min(w, h) // 2
size = f'{s}x{int(s*0.618)}'
small_size = f'{s}x{int(s*0.618/2)}'
big_size = f'{s}x{int(s*0.618*2)}'
root.geometry(small_size)
root.title('novowrap')
assembly = tk.LabelFrame(root, text='')
assembly.pack(side='left', padx=20)
a_button1 = tk.Button(assembly, text='Assembly sequences', command=assembly_ui)
a_button1.pack()
merge = tk.LabelFrame(root, text='')
merge.pack(side='left', padx=10, pady=50)
m_button1 = tk.Button(merge, text='Merge contigs', command=merge_ui)
m_button1.pack()
validate = tk.LabelFrame(root, text='')
validate.pack(side='left', padx=20)
v_button1 = tk.Button(validate, text='Validate assembly', command=validate_ui)
v_button1.pack()
# lf = tk.LabelFrame(validate, text='Unique')
# lf.pack(padx=20, pady=20)
# for i in 'first,longest,none'.split(','):
#    r = ttk.Radiobutton(lf, text=i, variable=tk.IntVar, width=10)
#    r.pack()
va = tk.IntVar(root)
vb = tk.StringVar(root)
va.set(100)
vc = tk.DoubleVar(root, 3.14159)
# print([i.get() for i in (va, vb, vc)])
root.mainloop()
