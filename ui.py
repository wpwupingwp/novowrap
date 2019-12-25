#!/usr/bin/python3

import logging
import queue
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog, scrolledtext

from pathlib import Path

from merge import merge_main
from validate import validate_main

log = logging.getLogger('novowrap')


class QueueHandler(logging.Handler):
    def __init__(self, log_queue):
        super().__init__()
        self.log_queue = log_queue

    def emit(self, msg):
        self.log_queue.put(msg)


def stext(root): 
    s = scrolledtext.ScrolledText(root)
    s.pack(fill='both')
    log_queue = queue.Queue()
    log = logging.getLogger('novowrap')
    log.addHandler(QueueHandler(log_queue))
    while log_queue.not_empty:
        msg = log_queue.get(block=False)
        s.insert('end', msg+'\n')
        s.yview('end')



def wlabel(root, text, row, column=0, width=25, padx=0, pady=0, sticky='EW',
           **kargs):
    label = tk.Label(root, text=text, width=width, **kargs)
    label.grid(row=row, column=column, padx=padx, pady=pady, sticky=sticky)
    return label


def entry(root, row, column, default=''):
    entry = tk.Entry(root)
    entry.insert(0, default)
    entry.grid(row=row, column=column)
    return entry


def info(message):
    messagebox.showinfo(message=message)


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
    def open_single():
        a = filedialog.askopenfilenames(title='Sequence files')
        w_entry.insert(0, a)

    w = tk.Toplevel(root)
    w.geometry(size)
    w_label = tk.Label(w, text='Input:')
    w_label.pack()
    w_entry = tk.Entry(w)
    w_entry.pack(side='right')
    w_button = tk.Button(w, text='Open', command=open_single)
    w_button.pack(side='right')


def assembly_batch():
    def open_single():
        a = filedialog.askopenfilenames(title='Table files')
        w_entry.insert(0, a)

    w = tk.Toplevel(root)
    w_label = tk.Label(w, text='Input:')
    w_label.pack()
    w_entry = tk.Entry(w)
    w_entry.pack(side='right')
    w_button = tk.Button(w, text='Open', command=open_single)
    w_button.pack(side='right')


def merge_single():
    def open_single():
        a = filedialog.askopenfilenames(title='fasta files')
        w_entry.insert(0, a)

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
        run.geometry('600x400')
        run.title('Running...')
        run.wm_transient()
        frame = tk.Frame(run)
        frame.pack(fill='both')
        stext(frame)
        r = validate_main(arg_str)
        info('Done.')
        run.destroy()
        # ok = messagebox.askokcancel(message=arg_str)

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
    i_entry = entry(w, row=row, column=1)
    i_button = tk.Button(w, text='Open', command=open_single('Input file',
                                                             i_entry))
    i_button.grid(row=row, column=2)
    ref = tk.LabelFrame(w, text='Reference')
    row += 1
    ref.grid(row=row, column=0, columnspan=5, padx=5, pady=8)
    row += 1
    wlabel(ref, 'Taxonomy', row=row, column=1)
    t_entry = entry(ref, row=row, column=2, default='Nicotiana tabacum')
    row += 1
    label1 = tk.Label(ref, text='OR', fg='red')
    label1.grid(row=row, column=0, columnspan=2)
    row += 1
    wlabel(ref, 'Genbank file', row=row, column=1, padx=13)
    r_entry = entry(ref, row=row, column=2)
    r_button = tk.Button(ref, text='Open',
                         command=open_single('Reference file', r_entry,
                                             t_entry))
    r_button.grid(row=row, column=3)
    row += 1
    wlabel(w, 'Output', row=row, padx=10, pady=8)
    o_entry = entry(w, row=row, column=1, default='"Current folder"')
    o_button = tk.Button(w, text='Open', command=open_folder('Output folder',
                                                             o_entry))
    o_button.grid(row=row, column=2)
    row += 1
    options = tk.LabelFrame(w, text='Options')
    options.grid(row=row, padx=5, columnspan=5)
    row += 1
    wlabel(options, 'Sequence similarity (0-1)', row=row, padx=10)
    s_entry = entry(options, row=row, column=1, default='0.7')
    row += 1
    wlabel(options, 'Length difference (0-1)', row=row, padx=10)
    l_entry = entry(options, row=row, column=1, default='0.2')
    row += 1
    ok = tk.Button(w, text='Enter', command=submit_validate)
    ok.grid(row=row, column=0, columnspan=3, sticky='EW', padx=50, pady=10)


def o():
    pass
    # d = simpledialog.SimpleDialog(
    #     dialog, title='dialog test', text='this is a dialog',
    #     buttons=['yes', 'no', 'wo bu zhidao'], cancel=3, default=1)
    #value = d.go()
    #print(value)


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
#lf = tk.LabelFrame(validate, text='Unique')
#lf.pack(padx=20, pady=20)
# for i in 'first,longest,none'.split(','):
#    r = ttk.Radiobutton(lf, text=i, variable=tk.IntVar, width=10)
#    r.pack()
va = tk.IntVar(root)
vb = tk.StringVar(root)
vc = tk.DoubleVar(root)
va.set(100)
vc = tk.DoubleVar(root, 3.14159)
#print([i.get() for i in (va, vb, vc)])
root.mainloop()
