#!/usr/bin/python3

import logging
import threading
import queue
import tkinter as tk
import time
from tkinter import messagebox, simpledialog, filedialog, scrolledtext

from pathlib import Path

from merge import merge_main
from validate import validate_main


class QueueHandler(logging.Handler):
    """
    For queue.
    """
    def __init__(self, log_queue_):
        super().__init__()
        self.log_queue = log_queue_

    def emit(self, msg):
        self.log_queue.put(msg)


def stext(window):
    """
    ScrolledText that shows logs.
    """
    scroll = scrolledtext.ScrolledText(window)
    scroll.pack(fill='both')

    def poll():
        while True:
            try:
                msg = log_queue.get(block=False)
                msg = formatter.format(msg) + '\n'
                scroll.insert('end', msg)
                scroll.yview('end')
            except queue.Empty:
                break
        scroll.after(100, poll)

    scroll.after(100, poll)


def wlabel(window, text, row, column=0, width=25, padx=0, pady=0, sticky='EW',
           **kargs):
    """
    Generate and pack labels.
    """
    label = tk.Label(window, text=text, width=width, **kargs)
    label.grid(row=row, column=column, padx=padx, pady=pady, sticky=sticky)
    return label


def fentry(window, row, column, default=''):
    """
    Generate and pack entrys.
    Fill with default string.
    """
    entry = tk.Entry(window)
    entry.insert(0, default)
    entry.grid(row=row, column=column)
    return entry


def info(message):
    """
    For shorter.
    """
    messagebox.showinfo(message=message)


def validate_wrap(arg_str, window):
    """
    Wrap for callback.
    Args:
        arg_str(str): strings for validate()
        window(Toplevel): window to destroy
    """
    validated, output_info = validate_main(arg_str)
    info(f'Done. See {output_info} for details.')
    window.destroy()
    return


def open_single(title, entry, entry2=None):
    """
    Set title, fill entry 1, empty entry 2.
    Only open one file.
    """
    def func():
        a = filedialog.askopenfilename(title=title)
        entry.insert(0, a)
        if entry2 is not None:
            entry2.delete(0, 'end')
    return func


def open_folder(title, entry):
    def func():
        a = filedialog.askdirectory(title=title)
        entry.insert(0, a)
    return func


def assembly_single():
    w = tk.Toplevel(root)
    w.geometry(size)
    w_label = tk.Label(w, text='Input:')
    w_label.pack()
    w_entry = tk.Entry(w)
    w_entry.pack(side='right')
    w_button = tk.Button(w, text='Open', command=open_single)
    w_button.pack(side='right')


def assembly_batch():
    w = tk.Toplevel(root)
    w_label = tk.Label(w, text='Input:')
    w_label.pack()
    w_entry = tk.Entry(w)
    w_entry.pack(side='right')
    w_button = tk.Button(w, text='Open', command=open_single)
    w_button.pack(side='right')


def merge_single():
    w = tk.Toplevel(root)
    w_label = tk.Label(w, text='Input:')
    w_label.pack()
    w_entry = tk.Entry(w)
    w_entry.pack(side='right')
    w_button = tk.Button(w, text='Open', command=open_single)
    w_button.pack(side='right')


def validate_ui():
    """
    UI of validate.
    """
    def submit_validate():
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

        wroot.destroy()
        run = tk.Toplevel(root)
        run.geometry(size)
        run.title('Running...')
        run.wm_transient()
        frame = tk.Frame(run)
        frame.pack(fill='both')
        stext(frame)
        r = threading.Thread(target=validate_wrap, args=(arg_str, run))
        r.start()
        run.mainloop()

    wroot = tk.Toplevel(root)
    wroot.geometry(size)
    wroot.title('Validate')
    # on top
    wroot.wm_transient(root)
    w = tk.Frame(wroot)
    w.grid(row=0, sticky='WENS', padx=30)
    # use variable for easily edit
    row = 0
    wlabel(w, 'Input', row=row, padx=15, pady=10)
    i_entry = fentry(w, row=row, column=1)
    i_button = tk.Button(w, text='Open', command=open_single('Input file',
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
    wlabel(ref, 'Genbank file', row=row, column=1, padx=13)
    r_entry = fentry(ref, row=row, column=2)
    r_button = tk.Button(ref, text='Open',
                         command=open_single('Reference file', r_entry,
                                             t_entry))
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


# define logger
formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
                              datefmt='%H:%M:%S')
log = logging.getLogger('novowrap')
log_queue = queue.Queue()
queue_handler = QueueHandler(log_queue)
queue_handler.setFormatter(formatter)
log.addHandler(queue_handler)

root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
s = min(w, h) // 2
size = f'{s}x{int(s*0.618)}'
root.geometry(size)
root.title('novowrap')
assembly = tk.LabelFrame(root, text='Assembly')
assembly.pack(side='left', padx=50)
a_button1 = tk.Button(assembly, text='Single Mode', command=assembly_single)
a_button1.pack()
a_button2 = tk.Button(assembly, text='Batch Mode', command=assembly_batch)
a_button2.pack()
merge = tk.LabelFrame(root, text='Merge')
merge.pack(side='left', padx=10, pady=50)
m_button1 = tk.Button(merge, text='Merge contigs', command=merge_single)
m_button1.pack()
validate = tk.LabelFrame(root, text='Validate')
validate.pack(side='left', padx=50)
v_button1 = tk.Button(validate, text='Validate', command=validate_ui)
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
