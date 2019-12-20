#!/usr/bin/python3

import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog
from tkinter import ttk


def quit(event):
    messagebox.showinfo(message=str(event))


def fff():
    a = filedialog.askopenfilenames(
        title='gb files', filetypes=[('Genbank file', '*.gb'),  ('All',
                                                                 '*.*')])
    ent.insert(0, a)
    return a


def ffo():
    a = filedialog.askdirectory( title='gb files') 
    return a


def c():
    raise SystemExit


def o():
    d = simpledialog.SimpleDialog(
        dialog, title='dialog test', text='this is a dialog',
        buttons = ['yes', 'no', 'wo bu zhidao'], cancel=3, default=1)
    value = d.go()
    print(value)


def ow():
    w = tk.Toplevel(root)
    w.geometry('500x309')

def pp(event):
    print('ent', ent.get())

root = tk.Tk()
w, h = root.winfo_screenwidth(), root.winfo_screenheight()
s = min(w, h) // 2
size = f'{s}x{int(s*0.618)}'
root.geometry(size)
root.title('Window')
n = ttk.Notebook(root)
f1 = ttk.Frame(n)
f2 = ttk.Frame(n)
n.add(f1, text='Advance', state='normal')
n.add(f2, text='Other')
n.pack()
f1.pack()
f2.pack()
assembly = tk.LabelFrame(root, text='assembly')
assembly.pack(side='left', padx=10, pady=10)
label = tk.Label(assembly, text='Test')
label.bind('<Button-1>', quit)
label.pack()
button = tk.Button(assembly, text='quit', command=c)
button.pack(anchor=tk.E, padx=30, pady=20, side=tk.LEFT)
validate = tk.LabelFrame(root, text='validate')
validate.pack()
lf = tk.LabelFrame(validate, text='Unique')
lf.pack(padx=20, pady=20)
for i in 'first,longest,none'.split(','):
    r = ttk.Radiobutton(lf, text=i, variable=tk.IntVar, width=10)
    r.pack()
dialog = tk.Frame(root)
dialog.pack()
fileframe = tk.Frame(dialog)
fileframe.pack()
ent = tk.Entry(fileframe)
ent.bind('<Button-1>', pp)
ent.pack(fill='both')
b2 = ttk.Button(fileframe, text='open files', command=fff)
b2.pack(side='right')
b3 = ttk.Button(dialog, text='output folder', command=ffo)
b3.pack(side='right')
bw = tk.Button(root, text='open new window', command=ow)
bw.pack(side='bottom')
va = tk.IntVar(root)
vb = tk.StringVar(root)
vc = tk.DoubleVar(root)
va.set(100)
vc = tk.DoubleVar(root, 3.14159)
print([i.get() for i in (va, vb, vc)])
root.mainloop()
