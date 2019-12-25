#!/usr/bin/python3

import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog

from pathlib import Path


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


def open_single(title, entry, entry2=None):
    """
    Set title, fill entry 1, empty entry 2.
    """
    def func():
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


def submit_validate():
    ok = messagebox.askokcancel('Run?')
    if ok:
        wroot.destroy()


def validate_single():
    global wroot
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


def ffo():
    a = filedialog.askdirectory(title='gb files')
    return a


def c():
    raise SystemExit


def o():
    pass
    # d = simpledialog.SimpleDialog(
    #     dialog, title='dialog test', text='this is a dialog',
    #     buttons=['yes', 'no', 'wo bu zhidao'], cancel=3, default=1)
    #value = d.go()
    #print(value)


def ow():
    w = tk.Toplevel(root)
    w.geometry('500x309')


def pp(event):
    pass
    #print('ent', ent.get())


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
v_button1 = tk.Button(validate, text='Validate', command=validate_single)
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
