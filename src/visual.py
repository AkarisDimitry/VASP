import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

import random, time, sys,  pickle

if sys.version_info[0] == 3:    # for Python3
    from tkinter import *   ## notice lowercase 't' in tkinter here
else:   # for Python2
    from Tkinter import *   ## notice capitalized T in Tkinter
import tkinter as tk
from tkinter import ttk
import tkinter.filedialog

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import matplotlib.animation as animation
from matplotlib.figure import Figure

from PIL import ImageTk, Image
'''
import PROCAR as pc
import DOSCAR as dc
import POSCAR as poscar
from ase.visualize.plot import plot_atoms
from ase.io import read, write
'''
class SeaofBTCapp(tk.Tk):
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs) 

        #tk.Tk.iconbitmap(self, default="clienticon.ico")
        tk.Tk.wm_title(self, "JML - Chemometrics")
        #self.geometry('500x500')
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.data = None
        self.frames = {}

        for F in (StartPage, ):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage,)

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()


class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)

        self.controller = controller 
        self.controller.protocol('WM_DELETE_WINDOW', self.exit)
        
        # **** **** **** **** **** **** **** LABELS organization **** **** **** **** **** **** **** # 
        # -- info(TOP) || plot (TOP)-- #
        self.label_top = tk.Label(self, )
        self.label_top.pack(side=TOP, pady=2,padx=2)

        # -- info(BOT) || plot (BOT)-- #
        self.label_bot = tk.Label(self, )
        self.label_bot.pack(side=BOTTOM, pady=2,padx=2)

        self.button_load_procar = tk.Button(self.label_bot, text="LOAD PROCAR",
                            command=self.load_procar, bg='#AAAAAA', font=("Verdana", 10))
        self.button_load_procar.pack(fill=X, )

        self.button_load_doscar = tk.Button(self.label_bot, text="LOAD DOSCAR",
                            command=self.load_doscar, bg='#AAAAAA', font=("Verdana", 10))
        self.button_load_doscar.pack(fill=X, )

    def load_doscar(self, ):
        DC = dc.DOSCAR()
        self.load_file(DC)
        self.add_plot_doscar(obj=[DC])
        return DC

    def load_procar(self, ):
        PC = pc.PROCAR()
        self.load_file(PC)
        self.add_plot(obj=PC)
        return PC

    def load_file(self, obj):
        input_file_name = tk.filedialog.askopenfilename(filetypes=[("All Files", "*"), ("Text Documents", "*.txt")])
        if input_file_name:
            file_name = input_file_name
            obj.load(file_name)


    def add_plot(self, event=None, obj=None):
        def act_plot(self,):
            fig.clear()
            obj.plot(figure=fig)
            canvas.draw()

        def command_bar():
            if orbitals[-2].get() == 1:
                for n in orbitals: n.set(0)

            if atoms[-1].get() == 1:
                for n in atoms: n.set(1)
                atoms[-2].set(0)
                atoms[-1].set(0)

            obj.plot_spins = [] 
            for i, n in enumerate(spins):
                if n.get(): obj.plot_spins.append(i)

            obj.plot_orbitals = []
            for i, n in enumerate(orbitals):
                if n.get(): obj.plot_orbitals.append(i)

            obj.plot_ions = []
            for i, n in enumerate(atoms):
                if n.get(): obj.plot_ions.append(i)

            act_plot(self,)

        try:
            # ---- General FONT ---- #
            font00 = ("Verdana", 10)

            # ****** Label to store all plots ****** #
            label_general = tk.Label(self.label_top )
            label_general.pack(pady=1,padx=1, side=LEFT)
            label_general.option_add('*tearOff', False)

            label_plot_top = tk.Label(label_general )
            label_plot_top.pack(pady=1,padx=1, side=TOP)

            label_bottom = tk.Label(label_general )
            label_bottom.pack(pady=1,padx=1, fill=X, side=BOTTOM)

            label_info_rigth = tk.Label(label_general )
            label_info_rigth.pack(pady=1,padx=1, side=RIGHT)

            label_info_left = tk.Label(label_general )
            label_info_left.pack(pady=1,padx=1, side=LEFT)
            #,width=700,bg="white",text="test",borderwidth=0

            # ****** remove buttom ****** #
            self.button_erase = tk.Button(label_bottom, text="Remove",
                                   command=label_general.destroy, bg='#AAAAAA', font=("Verdana", 10))
            self.button_erase.pack(fill=X, side=BOTTOM  ) 


            # ****** Select orbitals to plot ****** #
            # - combo - #
            orbitals_names = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 'all', 'None']
            orbitals = []

            mb=  Menubutton ( label_bottom, text="Select orbitals", relief=RAISED )
            mb.pack(fill=X, )
            mb.menu  =  Menu ( mb, tearoff = 0 )
            mb.option_add('*tearOff', False) 
            mb["menu"]  =  mb.menu

            for i, n in enumerate(orbitals_names):
                temp_var = IntVar()
                mb.menu.add_checkbutton ( label=n, variable=temp_var, command=command_bar)
                orbitals.append(temp_var)

            # ****** Select ATOMS to plot ****** #
            # - combo - #
            atoms_names = list(range(obj.n_oins))+['all', 'None']
            atoms = []

            mb=  Menubutton ( label_bottom, text="Select atoms", relief=RAISED )
            mb.pack(fill=X, )
            mb.menu  =  Menu ( mb, tearoff = 0 )
            mb.option_add('*tearOff', False) 
            mb["menu"]  =  mb.menu

            for i, n in enumerate(atoms_names):
                temp_var = IntVar()
                mb.menu.add_checkbutton ( label=str(n), variable=temp_var, command=command_bar)
                if (i+1)%11 == 0:  mb.menu.entryconfigure(i, columnbreak=1)
                atoms.append(temp_var)

            # ****** Select SPIN to plot ****** #
            # - combo - #
            spin_names = ['UP', 'DOWN']
            spins = []

            mb=  Menubutton ( label_bottom, text="Select spin", relief=RAISED )
            mb.pack(fill=X, )
            mb.menu  =  Menu ( mb, tearoff = 0 )
            mb.option_add('*tearOff', False) 
            mb["menu"]  =  mb.menu

            for i, n in enumerate(spin_names):
                temp_var = IntVar()
                mb.menu.add_checkbutton ( label=n, variable=temp_var, command=command_bar)
                spins.append(temp_var)


            # ****** PLT GRAPH ****** # in label_general
            fig = Figure(figsize=(5, 6), dpi=100) # figure size and resolution # 

            fig_plot = fig.add_subplot(111)
            fig.axes[0].format_coord = lambda x, y: ""

            canvas = FigureCanvasTkAgg(fig, label_plot_top)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.BOTTOM, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, label_plot_top,)
            toolbar.update()
            canvas._tkcanvas.pack(fill=X, side=tk.TOP, expand=True )
            act_plot(self,)
            # *********************************************** #
          
        except: print('ERRO :: code 000 :: main.add_plot() :: can NOT plot ')



    def add_plot_doscar(self, event=None, obj=None):
        def add_plot():
            obj.append( load_doscar() )
            act_plot()

        def load_doscar():
            DC = dc.DOSCAR()
            load_file(DC)
            return DC

        def load_file(obj):
            input_file_name = tk.filedialog.askopenfilename(filetypes=[("All Files", "*"), ("Text Documents", "*.txt")])
            if input_file_name:
                file_name = input_file_name
                obj.load(file_name)

        def act_plot():
            fig.clear()

            if plot_sum.get() == 1:  plot_sum_logic = True
            else:                    plot_sum_logic = False

            if plot_legend.get() == 1:  plot_legend_logic = True
            else:                       plot_legend_logic = False

            for i, n in enumerate(obj):
                n.plot(figure=fig, color=i*2, sum_plot=plot_sum_logic, legend_plot=plot_legend_logic)
            canvas.draw()

        def command_bar():

            if orbitals[-1].get() == 1:
                for n in orbitals: n.set(0)

            if orbitals[-2].get() == 1:
                for n in orbitals: n.set(1)
                orbitals[-2].set(0);  orbitals[-1].set(0)

            if atoms[-1].get() == 1:
                for n in atoms: n.set(0)

            if atoms[-1].get() == 1:
                for n in atoms: n.set(1)
                atoms[-2].set(0);    atoms[-1].set(0)

            for o in obj: 
                o.plot_orbitals = []
                for i, n in enumerate(orbitals):
                    if n.get(): o.plot_orbitals.append(i+1)

            for o in obj: 
                o.plot_ions = []
                for i, n in enumerate(atoms):
                    if n.get(): o.plot_ions.append(i)

            act_plot()

        def DOSCAR_export():
            input_file_name = tk.filedialog.asksaveasfilename(defaultextension=".dat",
                                                   filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")])
        
            if input_file_name:
                file_name = input_file_name
                for i, n in enumerate(obj):       n.export(file_name='{}_{}'.format(file_name, i) )

        def select_group_s():
            # orbital S #  # UP/DOWN #
            orbitals[0].set(1)
            orbitals[1].set(1)
            command_bar()

        def select_group_p():
            # orbital px # # UP/DOWN #
            orbitals[2].set(1)
            orbitals[3].set(1)
            # orbital py # # UP/DOWN #
            orbitals[4].set(1)
            orbitals[5].set(1)
            # orbital pz # # UP/DOWN #
            orbitals[6].set(1)
            orbitals[7].set(1)
            command_bar()

        def select_group_d():
            # orbital dxy # # UP/DOWN #
            orbitals[8].set(1)
            orbitals[9].set(1)
            # orbital dyz # # UP/DOWN #
            orbitals[10].set(1)
            orbitals[11].set(1)
            # orbital dz2 # # UP/DOWN #
            orbitals[12].set(1)
            orbitals[13].set(1)
            # orbital dxz # # UP/DOWN #
            orbitals[14].set(1)
            orbitals[15].set(1)
            # orbital dx2 # # UP/DOWN #
            orbitals[16].set(1)
            orbitals[17].set(1)
            command_bar()
        #try:
        if 1 ==1 :
            # ---- General FONT ---- #
            font00 = ("Verdana", 10)

            # ****** Label to store all plots ****** #
            label_general = tk.Label(self.label_top )
            label_general.pack(pady=1,padx=1, side=LEFT)
            label_general.option_add('*tearOff', False)

            label_plot_top = tk.Label(label_general )
            label_plot_top.pack(pady=1,padx=1, side=TOP)

            label_bottom = tk.Label(label_general )
            label_bottom.pack(pady=1,padx=1, fill=X, side=BOTTOM)

            label_info_rigth = tk.Label(label_general )
            label_info_rigth.pack(pady=1,padx=1, side=RIGHT)

            label_info_left = tk.Label(label_general )
            label_info_left.pack(pady=1,padx=1, side=LEFT)
            #,width=700,bg="white",text="test",borderwidth=0

            # ****** remove buttom ****** #
            self.button_erase = tk.Button(label_bottom, text="Remove",
                                   command=label_general.destroy, bg='#AAAAAA', font=("Verdana", 10))
            self.button_erase.pack(fill=X, side=BOTTOM  ) 

            # ****** ADD plot buttom ****** #
            self.button_addplot = tk.Button(label_bottom, text="Add DOSCAR",
                                   command=add_plot, bg='#AAAAAA', font=("Verdana", 10))
            self.button_addplot.pack(fill=X, side=BOTTOM  ) 



            # ****** Select orbitals to plot ****** #
            label_orbitals_atoms = tk.Label(label_bottom )
            label_orbitals_atoms.pack(fill=X, pady=1,padx=1, )

            # - combo - #
            orbitals_names = [  's(up)',    's(down)',      'px(up)',   'px(down)',     'py(up)',   'py(down)',     'pz(up)',       'pz(down)',     'dxy(up)', 'dxy(down)', 
                                'dyz(up)',  'dyz(down)',    'dz2(up)',  'dz2(down)',    'dxz(up)',  'dxz(down)',    'dx2(down)',    'dx2(down)',   
                                'all', 'None']
            orbitals = []
            mb=  Menubutton ( label_orbitals_atoms, text="Select orbitals", relief=RAISED, font=("Verdana", 15), padx=80 )
            mb.pack(side=RIGHT, fill=X)
            mb.menu  =  Menu ( mb, tearoff = 0 )
            mb.option_add('*tearOff', False) 
            mb["menu"]  =  mb.menu

            for i, n in enumerate(orbitals_names):
                temp_var = IntVar()
                mb.menu.add_checkbutton ( label=n, variable=temp_var, command=command_bar)
                orbitals.append(temp_var)

            # ****** Select ATOMS to plot ****** #
            # - combo - #
            atoms_names = list(range(obj[0].n_oins))+['all', 'None']
            atoms = []

            mb=  Menubutton ( label_orbitals_atoms, text="Select atoms", relief=RAISED, font=("Verdana", 15), padx=80 )
            mb.pack(side=RIGHT, fill=X)
            mb.menu  =  Menu ( mb, tearoff = 0 )
            mb.option_add('*tearOff', False) 
            mb["menu"]  =  mb.menu

            for i, n in enumerate(atoms_names):
                temp_var = IntVar()
                mb.menu.add_checkbutton ( label=str(n), variable=temp_var, command=command_bar)
                if (i+1)%11 == 0:  mb.menu.entryconfigure(i, columnbreak=1)
                atoms.append(temp_var)

            # ****** Select orbitals GROUPs to plot ****** #
            label_info_bottom = tk.Label(label_bottom )
            label_info_bottom.pack(fill=X, pady=1,padx=1, )

            self.button_erase = tk.Button(label_info_bottom, text="d(up) + d(down)",
                                   command=select_group_d, bg='#AAAAAA', font=("Verdana", 10), padx=45)
            self.button_erase.pack(fill=X, side=RIGHT  ) 

            self.button_erase = tk.Button(label_info_bottom, text="p(up) + p(down)",
                                   command=select_group_p, bg='#AAAAAA', font=("Verdana", 10), padx=45)
            self.button_erase.pack(fill=X, side=RIGHT  ) 

            self.button_erase = tk.Button(label_info_bottom, text="s(up) + s(down)",
                                   command=select_group_s, bg='#AAAAAA', font=("Verdana", 10), padx=45)
            self.button_erase.pack(fill=X, side=LEFT  )

            # ****** plot preSets ****** #
            label_plot_set = tk.Label(label_bottom )
            label_plot_set.pack(fill=X, pady=1,padx=1, )

            plot_sum = IntVar()
            self.button_sumplot = tk.Checkbutton(label_plot_set, text="Plot SUM", variable=plot_sum,
                                   command=act_plot, bg='#AAAAAA', font=("Verdana", 10))
            self.button_sumplot.pack(fill=X, side=LEFT  ) 

            plot_legend = IntVar()
            self.button_refplot = tk.Checkbutton(label_plot_set, text="Show legend", variable=plot_legend,
                                   command=act_plot, bg='#AAAAAA', font=("Verdana", 10))
            self.button_refplot.pack(fill=X, side=LEFT  ) 

            # ****** SAVE plot data ****** #
            self.button_export = tk.Button(label_bottom, text="Export data",
                                   command=DOSCAR_export, bg='#AAAAAA', font=("Verdana", 10))
            self.button_export.pack(fill=X, side=BOTTOM  ) 


            # ****** PLT GRAPH ****** # in label_general
            fig = Figure(figsize=(5, 6), dpi=100) # figure size and resolution # 

            fig_plot = fig.add_subplot(111)
            fig.axes[0].format_coord = lambda x, y: ""

            canvas = FigureCanvasTkAgg(fig, label_plot_top)
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.BOTTOM, expand=True)

            toolbar = NavigationToolbar2Tk(canvas, label_plot_top,)
            toolbar.update()
            canvas._tkcanvas.pack(fill=X, side=tk.TOP, expand=True )
            act_plot()
            # *********************************************** #
          
        #except: print('ERRO :: code 000 :: main.add_plot() :: can NOT plot ')


    def exit(self, event=None):
        if tkinter.messagebox.askokcancel("Quit?", "Do you want to QUIT for sure?\n Make sure you've saved your current work."):
            self.controller.destroy()


POSCAR = poscar.POSCAR()
DC = dc.DOSCAR()

app = SeaofBTCapp()
app.mainloop()









