from datetime import datetime
from collections import OrderedDict
import warnings
import logging

import numpy as np
import scipy as scipy
import scipy.stats as stats
import tkinter as tk
import tkinter.messagebox
from textwrap import wrap
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('error', category=RuntimeWarning)
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=scipy.linalg.LinAlgWarning)

today = datetime.today().isoformat().split("T")[0]

class GUI(tk.Frame):
    # variables
    start_run = 0
    quit_run = 0

    # functions
    def __init__(self, parent):
        self.root = tk.Tk()
        self.parent = parent
        self.root.title('MOLUSC')
        self.create_widgets()

    def start(self):
        self.root.mainloop()

    def create_widgets(self):

        self.root.configure(bg='#ECECEC')

        # Create main panels
        Analysis_Options = tk.Frame(self.root, relief=tk.SUNKEN, borderwidth=4, bg='#ECECEC')
        Output = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
        Stellar_Info = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
        Advanced_Options = tk.Frame(self.root, width=200, height=200, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
        Buttons = tk.Frame(self.root, bg='#ECECEC', height=50)
        Display = tk.Frame(self.root, borderwidth=2, relief=tk.SUNKEN, bg='#ECECEC')

        # Place main panels
        Analysis_Options.grid(column=0, row=0, sticky='nsew')
        Stellar_Info.grid(column=0, row=1, sticky='ew')
        Output.grid(column=0, row=2,  sticky='ew')
        Advanced_Options.grid(column=1, row=0, columnspan=2, rowspan=2, sticky='nsew')
        Buttons.grid(column=2, row=2, sticky='se')
        Display.grid(column=3, row=0, rowspan=4, sticky='nsew')

        # Analysis Options
        analysis_label = tk.Label(Analysis_Options, text="Analysis Options", bg='#ECECEC')
        analysis_label.grid(column=0, row=0, columnspan=3)

        vcmd1 = (self.root.register(self.validate_resolution), '%P')

        #  toggle row
        checkbox_frame = tk.Frame(Analysis_Options, borderwidth=0, bg='#ECECEC')
        self.__ao_check = tk.BooleanVar()
        self.__rv_check = tk.BooleanVar()
        self.__ruwe_check = tk.BooleanVar()
        self.__gaia_check = tk.BooleanVar()

        ao_box = tk.Checkbutton(checkbox_frame, text="HRI", variable=self.__ao_check, command=self.activate_ao, bg='#ECECEC')
        rv_box = tk.Checkbutton(checkbox_frame, text="RV", variable=self.__rv_check, command=self.activate_rv, bg='#ECECEC')
        ruwe_box = tk.Checkbutton(checkbox_frame, text='RUWE', variable=self.__ruwe_check, bg='#ECECEC')
        gaia_box = tk.Checkbutton(checkbox_frame, text='Gaia', variable=self.__gaia_check, bg='#ECECEC')

        # sub frame for AO, necessary from multiple AO files
        self.ao_frame = tk.Frame(Analysis_Options, borderwidth=0)

        # AO variables and widgets go into the separate AO frame so that I can add multiple rows of them if needed
        self.ao_rows = 1
        self.ao_file_boxes = []
        self.ao_filter_menus = []
        # label column
        rv_label = tk.Label(Analysis_Options, text='RV File:', bg='#ECECEC')
        ao_label = tk.Label(self.ao_frame, text='Contrast File:', bg='#ECECEC')
        # entry column
        self.__ao_file = tk.StringVar()
        self.__rv_file = tk.StringVar()
        self.__ao_file_box = tk.Entry(self.ao_frame, textvariable=self.__ao_file, state=tk.DISABLED, width=25, disabledbackground='#ECECEC', highlightthickness=0, bd=1)
        self.ao_file_boxes.append(self.__ao_file_box) # adds the first file box to the list
        self.__rv_file_box = tk.Entry(Analysis_Options, textvariable=self.__rv_file, state=tk.DISABLED, width=25, disabledbackground='#ECECEC', highlightthickness=0, bd=1)
        # AO filter menu
        self.__filter_str = tk.StringVar()
        self.__filter_menu = tk.OptionMenu(self.ao_frame, self.__filter_str, 'Filter', 'J', 'H', 'K', 'G','Bp','Rp', 'R', 'I', 'L', 'LL', 'M')
        self.__filter_str.set('Filter')
        self.__filter_menu.config(state=tk.DISABLED, highlightthickness=0)
        self.ao_filter_menus.append(self.__filter_str)
        # AO add button
        self.__add_button = tk.Button(self.ao_frame, text='+', command=self.add_ao, state=tk.DISABLED, width=1)
        # RV resolution Box
        resolution_label = tk.Label(Analysis_Options, text=' Resolution:', bg='#ECECEC')
        self.__resolution = 0.
        self.__resolution_box = tk.Entry(Analysis_Options, validate='focusout', vcmd=vcmd1, state=tk.DISABLED, width=6, disabledbackground='#ECECEC', highlightthickness=0, bd=1)

        # Grid
        # checkbox row
        checkbox_frame.grid(column=0, row=1, columnspan=4, pady=10)
        ao_box.grid(column=0, row=0)
        rv_box.grid(column=1, row=0)
        ruwe_box.grid(column=3, row=0)
        gaia_box.grid(column=4, row=0)
        #  ao
        self.ao_frame.grid(column=0, columnspan=4, row=2)
        ao_label.grid(column=0, row=1, sticky='e')
        self.__ao_file_box.grid(column=1, row=1)
        self.__filter_menu.grid(column=2, row=1, padx=2)
        self.__add_button.grid(column=3, row=1)
        #  rv row
        rv_label.grid(column=0, row=3, sticky='w')
        self.__rv_file_box.grid(column=1, row=3, pady=15)
        resolution_label.grid(column=2, row=3, pady=15)
        self.__resolution_box.grid(column=3, row=3, pady=15)

        # Stellar Info/Basic Options
        stellar_label = tk.Label(Stellar_Info, text='Star Info', bg='#ECECEC')
        stellar_label.grid(column=0, row=0, columnspan=4)

        #  label columns
        number_label = tk.Label(Stellar_Info, text='Generated Companions*:', bg='#ECECEC')
        ra_label = tk.Label(Stellar_Info, text='RA (hms)*:', bg='#ECECEC')
        dec_label = tk.Label(Stellar_Info, text='DEC (dms)*:', bg='#ECECEC')
        mass_label = tk.Label(Stellar_Info, text='Star Mass (M_sun)*:', bg='#ECECEC')
        age_label = tk.Label(Stellar_Info, text='Star Age (Gyr):', bg='#ECECEC')
        added_jitter_label = tk.Label(Stellar_Info, text='Added Jitter (m/s):', bg='#ECECEC')
        rv_floor_label = tk.Label(Stellar_Info, text='RV Floor (m/s):', bg='#ECECEC')

        #  entry columns
        self.ra_str = ''
        self.dec_str = ''
        self.mass = 0
        self.age = 5
        self.num_generated = 0
        self.added_jitter = 20
        self.rv_floor = 20
        vcmd1 = (self.root.register(self.validate_ra), '%P')
        vcmd2 = (self.root.register(self.validate_dec), '%P')
        vcmd3 = (self.root.register(self.validate_num), '%P')
        vcmd4 = (self.root.register(self.validate_mass), '%P')
        vcmd5 = (self.root.register(self.validate_age), '%P')
        vcmd6 = (self.root.register(self.validate_jitter), '%P')
        vcmd7 = (self.root.register(self.validate_floor), '%P')

        self.number_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd3, highlightthickness=0, bd=1)
        self.ra_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd1, highlightthickness=0, bd=1)
        self.dec_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd2, highlightthickness=0, bd=1)
        self.mass_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd4, highlightthickness=0, bd=1)
        self.age_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd5, highlightthickness=0, bd=1)
        self.age_box.insert(-1, '5')
        self.age_box.config(fg='gray')
        self.added_jitter_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd6, highlightthickness=0, bd=1)
        self.added_jitter_box.insert(-1, '20')
        self.added_jitter_box.config(fg='gray')
        self.rv_floor_box = tk.Entry(Stellar_Info, width=15, validate='focusout', vcmd=vcmd7, highlightthickness=0, bd=1)
        self.rv_floor_box.insert(-1, '20')
        self.rv_floor_box.config(fg='gray')


        #  help button column
        age_help = tk.Button(Stellar_Info, text='?', command=self.help_age, height=1, width=2, highlightthickness=0, bd=3)
        ra_help = tk.Button(Stellar_Info, text='?', command=self.help_coord, height=1, width=2, highlightthickness=0, bd=3)
        dec_help = tk.Button(Stellar_Info, text='?', command=self.help_coord, height=1, width=2, highlightthickness=0, bd=3)
        added_jitter_help = tk.Button(Stellar_Info, text='?', command=self.help_added_jitter, height=1, width=2, highlightthickness=0, bd=3)
        rv_floor_help = tk.Button(Stellar_Info, text='?', command=self.help_floor, height=1, width=2, highlightthickness=0, bd=3)

        # Grid
        Stellar_Info.grid_columnconfigure(0, weight=1)
        Stellar_Info.grid_columnconfigure(2, weight=1)
        Stellar_Info.grid_rowconfigure(0, weight=1)
        Stellar_Info.grid_rowconfigure(7, weight=1)
        # labels
        number_label.grid(column=0, row=1, sticky='e')
        ra_label.grid(column=0, row=2, sticky='e')
        dec_label.grid(column=0, row=3, sticky='e')
        mass_label.grid(column=0, row=4, sticky='e')
        age_label.grid(column=0, row=5, sticky='e')
        added_jitter_label.grid(column=0, row=6, sticky='e')
        rv_floor_label.grid(column=0, row=7, sticky='e')
        # boxes
        self.number_box.grid(column=1, row=1)
        self.ra_box.grid(column=1, row=2)
        self.dec_box.grid(column=1, row=3)
        self.mass_box.grid(column=1, row=4)
        self.age_box.grid(column=1, row=5)
        self.added_jitter_box.grid(column=1, row=6)
        self.rv_floor_box.grid(column=1, row=7)
        # buttons
        ra_help.grid(column=2, row=2, sticky='w')
        dec_help.grid(column=2, row=3, sticky='w')
        age_help.grid(column=2, row=5, sticky='w')
        added_jitter_help.grid(column=2, row=6, sticky='w')
        rv_floor_help.grid(column=2, row=7, sticky='w')

        # Output
        self.__prefix = tk.StringVar()
        self.__extra = tk.IntVar()
        self.__all_out = tk.IntVar()
        output_label = tk.Label(Output, text='Output', bg='#ECECEC')
        prefix_label = tk.Label(Output, text='File Prefix*:', bg='#ECECEC')
        self.prefix_box = tk.Entry(Output, textvariable=self.__prefix, width=20, highlightthickness=0, bd=1)
        extra_box = tk.Checkbutton(Output, text='Extra Output', variable=self.__extra, pady=9, bg='#ECECEC')
        all_out_box = tk.Checkbutton(Output, text='Write Out All', variable=self.__all_out, bg='#ECECEC')

        Output.grid_columnconfigure(0, weight=1)
        Output.grid_columnconfigure(2, weight=1)
        Output.grid_rowconfigure(0, weight=1)
        output_label.grid(column=0, row=0, columnspan=3)
        prefix_label.grid(column=0, row=1, sticky='e')
        self.prefix_box.grid(column=1, row=1, columnspan=2, sticky='w')
        extra_box.grid(column=0, row=3, sticky='e')
        all_out_box.grid(column=2, row=3)

        # Advanced Options
        self.__P_fixed = None
        self.__P_min = None
        self.__P_max = None
        self.__i_fixed = None
        self.__i_min = None
        self.__i_max = None
        self.__e_fixed = None
        self.__e_min = None
        self.__e_max = None
        self.__arg_peri_fixed = None
        self.__arg_peri_min = None
        self.__arg_peri_max = None
        self.__m_fixed = None
        self.__m_min = None
        self.__m_max = None
        self.__a_fixed = None
        self.__a_min = None
        self.__a_max = None
        self.__phi_fixed = None
        self.__phi_min = None
        self.__phi_max = None
        self.__pd_mu = 5.03
        self.__pd_sig = 2.28
        self.__q_exp = 0.0
        self.__gaia_limit = tk.IntVar()
        self.__gaia_limit.set(18)

        advanced_label = tk.Label(Advanced_Options, text='Advanced Options', bg='#ECECEC')

        # Frames
        grid_frame = tk.Frame(Advanced_Options, borderwidth=0, bg='#ECECEC')
        dist_frame = tk.Frame(Advanced_Options, borderwidth=0, bg='#ECECEC')

        #  Grid Labels
        fixed_label = tk.Label(grid_frame, text='fixed', bg='#ECECEC', width=7)
        min_label = tk.Label(grid_frame, text='min', bg='#ECECEC', width=7)
        max_label = tk.Label(grid_frame, text='max', bg='#ECECEC', width=7)
        period_label = tk.Label(grid_frame, text='Period (days)', bg='#ECECEC')
        inc_label = tk.Label(grid_frame, text='cos(i)', bg='#ECECEC')
        ecc_label = tk.Label(grid_frame, text='Eccentricity', bg='#ECECEC')
        arg_peri_label = tk.Label(grid_frame, text='Arg Periapsis (radians)', bg='#ECECEC')
        mass_ratio_label = tk.Label(grid_frame, text='Mass Ratio (M2/M1)', bg='#ECECEC')
        a_label = tk.Label(grid_frame, text='Semi-Major Axis (AU)', bg='#ECECEC')
        phase_label = tk.Label(grid_frame, text='Phase (radians)', bg='#ECECEC')
        # Split Label
        split_label = tk.Label(Advanced_Options, text='- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -', justify='center', bg='#ECECEC')
        # Distribution Labels
        period_dist_label = tk.Label(dist_frame, text='Period Distribution', justify='center', bg='#ECECEC')
        mu_label = tk.Label(dist_frame, text=u'\u03BC:', bg='#ECECEC')
        sig_label = tk.Label(dist_frame, text=u'\u03C3:', bg='#ECECEC')
        mass_exp_label = tk.Label(dist_frame, text='Mass Ratio Distribution', justify='center', bg='#ECECEC')
        exp_label = tk.Label(dist_frame, text=u'\u03B3:', bg='#ECECEC')
        #  Validation Command registration
        vcmd_p = (self.root.register(self.validate_P), '%P', '%W')
        vcmd_i = (self.root.register(self.validate_i), '%P', '%W')
        vcmd_e = (self.root.register(self.validate_e), '%P', '%W')
        vcmd_w = (self.root.register(self.validate_arg_peri), '%P', '%W')
        vcmd_m = (self.root.register(self.validate_m), '%P', '%W')
        vcmd_a = (self.root.register(self.validate_a), '%P', '%W')
        vcmd_phi = (self.root.register(self.validate_phi), '%P', '%W')
        vcmd_pd_mu = (self.root.register(self.validate_pd_mu), '%P')
        vcmd_pd_sig = (self.root.register(self.validate_pd_sig), '%P')
        vcmd_mass_exp = (self.root.register(self.validate_mass_exp), '%P')
        #  Grid Entries
        entry_width = 7
        self.P1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
        self.P2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
        self.P3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_p, highlightthickness=0, bd=1)
        self.i1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
        self.i2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
        self.i3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_i, highlightthickness=0, bd=1)
        self.e1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
        self.e2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
        self.e3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_e, highlightthickness=0, bd=1)
        self.w1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
        self.w2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
        self.w3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_w, highlightthickness=0, bd=1)
        self.m1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
        self.m2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
        self.m3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_m, highlightthickness=0, bd=1)
        self.a1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
        self.a2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
        self.a3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_a, highlightthickness=0, bd=1)
        self.phi1_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
        self.phi2_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
        self.phi3_box = tk.Entry(grid_frame, width=entry_width, validate='focusout', vcmd=vcmd_phi, highlightthickness=0, bd=1)
        # Distribution Entries
        self.PD_mu_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_pd_mu, highlightthickness=0, bd=1)
        self.PD_mu_box.insert(-1, '5.03')
        self.PD_mu_box.config(fg='gray')
        self.PD_sig_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_pd_sig, highlightthickness=0, bd=1)
        self.PD_sig_box.insert(-1, '2.28')
        self.PD_sig_box.config(fg='gray')
        self.mass_exp_box = tk.Entry(dist_frame, width=entry_width, validate='focusout', vcmd=vcmd_mass_exp, highlightthickness=0, bd=1)
        self.mass_exp_box.insert(-1, '0.0')
        self.mass_exp_box.config(fg='gray')
        # Gaia Limit
        gaia_label = tk.Label(dist_frame, text='Gaia Completeness Limit:', bg='#ECECEC')
        gaia_18 = tk.Radiobutton(dist_frame, text="18th", variable=self.__gaia_limit, value=18)
        gaia_20 = tk.Radiobutton(dist_frame, text="20th", variable=self.__gaia_limit, value=20)

        #  Help buttons
        help_button = tk.Button(Advanced_Options, text='?', command=self.help_advanced, bg='#ECECEC', width=3, highlightthickness=0, relief=tk.RAISED, bd=3)
        period_help = tk.Button(dist_frame, text='?', command=self.help_period, height=1, width=2, highlightthickness=0, bd=3)
        mass_help = tk.Button(dist_frame, text='?', command=self.help_mass, height=1, width=2, highlightthickness=0, bd=3)

        #  Full Grid
        advanced_label.grid(column=0, columnspan=3, row=0)
        help_button.grid(column=2, row=0)
        grid_frame.grid(column=0, row=1, columnspan=3)
        split_label.grid(column=0, row=2, columnspan=3)
        dist_frame.grid(column=0, row=3, columnspan=3)
        #  limit grid
        fixed_label.grid(column=2, row=1)
        min_label.grid(column=3, row=1)
        max_label.grid(column=4, row=1)
        period_label.grid(column=0, row=2)
        inc_label.grid(column=0, row=3)
        ecc_label.grid(column=0, row=4)
        arg_peri_label.grid(column=0, row=5)
        mass_ratio_label.grid(column=0, row=6)
        a_label.grid(column=0, row=7)
        phase_label.grid(column=0, row=8)
        self.P1_box.grid(column=2, row=2)
        self.P2_box.grid(column=3, row=2)
        self.P3_box.grid(column=4, row=2)
        self.i1_box.grid(column=2, row=3)
        self.i2_box.grid(column=3, row=3)
        self.i3_box.grid(column=4, row=3)
        self.e1_box.grid(column=2, row=4)
        self.e2_box.grid(column=3, row=4)
        self.e3_box.grid(column=4, row=4)
        self.w1_box.grid(column=2, row=5)
        self.w2_box.grid(column=3, row=5)
        self.w3_box.grid(column=4, row=5)
        self.m1_box.grid(column=2, row=6)
        self.m2_box.grid(column=3, row=6)
        self.m3_box.grid(column=4, row=6)
        self.a1_box.grid(column=2, row=7)
        self.a2_box.grid(column=3, row=7)
        self.a3_box.grid(column=4, row=7)
        self.phi1_box.grid(column=2, row=8)
        self.phi2_box.grid(column=3, row=8)
        self.phi3_box.grid(column=4, row=8)
        # distribution
        period_dist_label.grid(column=0, row=0)
        mass_exp_label.grid(column=0, row=1)
        mu_label.grid(column=1, row=0)
        self.PD_mu_box.grid(column=2, row=0)
        sig_label.grid(column=3, row=0)
        self.PD_sig_box.grid(column=4, row=0)
        period_help.grid(column=5, row=0)
        exp_label.grid(column=1,row=1)
        self.mass_exp_box.grid(column=2, row=1)
        mass_help.grid(column=5, row=1)
        # gaia
        gaia_label.grid(column=0, row=2, columnspan=2)
        gaia_18.grid(column=2, row=2)
        gaia_20.grid(column=4, row=2)

        # Buttons
        self.run_button = tk.Button(Buttons, text='Run', command=self.run_code, width=6, height=2, highlightthickness=0, bd=3, bg='#ECECEC')
        close_button = tk.Button(Buttons, text='Close', command=self.close, width=8, height=2, highlightthickness=0, bd=3, bg='#ECECEC')

        self.run_button.grid(column=0, row=0)
        close_button.grid(column=1, row=0, padx=4)

        # Display
        self.display_box_text = tk.StringVar()
        self.status_box_text = tk.StringVar()
        self.status_box_text.set('Status: Waiting for Input')

        self.display_box = tk.Text(Display, height=29, width=40, background='white')
        self.status_box = tk.Label(Display, textvariable=self.status_box_text, width=33, borderwidth=4, relief=tk.SUNKEN, bg='#ECECEC')
        scrollbar = tk.Scrollbar(Display, orient=tk.VERTICAL, command=self.display_box.yview, width=15)

        self. display_box.configure(yscrollcommand=scrollbar.set)

        self.status_box.grid(column=0, row=0, columnspan=2)
        self.display_box.grid(column=0, row=1)
        scrollbar.grid(column=1, row=1)

        # end create_widgets

    def activate_ao(self):
        if self.__ao_check.get() == 1:  # whenever checked
            self.__ao_file_box.config(state=tk.NORMAL)
            self.__filter_menu.config(state=tk.NORMAL)
            self.__add_button.config(state=tk.NORMAL)
        elif self.__ao_check.get() == 0:  # whenever unchecked
            self.__ao_file_box.config(state=tk.DISABLED)
            self.__filter_menu.config(state=tk.DISABLED)
            self.__add_button.config(state=tk.DISABLED)

    def activate_rv(self):
        if self.__rv_check.get() == 1:  # whenever checked
            self.__rv_file_box.config(state=tk.NORMAL)
            self.__resolution_box.config(state=tk.NORMAL)
        elif self.__rv_check.get() == 0:  # whenever unchecked
            self.__rv_file_box.config(state=tk.DISABLED)
            self.__resolution_box.config(state=tk.DISABLED)

    def add_ao(self):
        # Adds another row with entry & filter for mutliple AO files
        self.ao_rows += 1
        ao_file = tk.StringVar()
        filter_str = tk.StringVar()
        # Create Widgits
        label = tk.Label(self.ao_frame, text=('Contrast File ' + str(self.ao_rows) +  ':'), bg='#ECECEC')
        file_box = tk.Entry(self.ao_frame, textvariable=ao_file, width=25, highlightthickness=0, bd=1)
        filter_menu = tk.OptionMenu(self.ao_frame, self.__filter_str, 'Filter', 'J', 'H', 'K', 'G','Bp','Rp', 'R', 'I', 'L', 'LL', 'M')
        filter_str.set('Filter')
        # Place Widgets
        label.grid(column=0, row=self.ao_rows)
        file_box.grid(column=1, row=self.ao_rows)
        filter_menu.grid(column=2, row=self.ao_rows)
        # Add Widgets to list
        self.ao_file_boxes.append(file_box)
        self.ao_filter_menus.append(filter_str)

    def check_completeness(self):
        # When the user presses run need to check if all required fields are filled.
        # file names, coordinates, mass, prefix, num generated all need to be filled

        allow_run = True
        message = 'The following problems were detected:\n'

        if not self.__ao_check.get() and not self.__rv_check.get() and not self.__ruwe_check.get() and not self.__gaia_check.get():
            message = message + '- No analysis is selected. Please select at least one type of analysis.\n'
            allow_run = False
        # Analysis Options
        if self.__ao_check.get():
            # The AO test requires an input file and a filter selection
            if not  self.__ao_file.get():
                self.__ao_file_box.config(bg='lightcoral')
                allow_run = False
                message = message +'- No file containing AO data has been given.\n'
            if self.__filter_str.get() == 'Filter':
                allow_run = False
                message = message + '- No filter is selected. Please choose the filter that the AO data is in.\n'
        if self.__rv_check.get():
            # The RV test requires an input file and a resolution
            if not self.__rv_file.get():
                self.__rv_file_box.config(bg='lightcoral')
                allow_run = False
                message = message +'- No file containing RV data has been given.\n'
            if self.__resolution == 0.:
                self.__resolution_box.config(bg='lightcoral')
                allow_run = False
                message = message +'- No RV Resolution provided.\n'
        # Coordinates
        if not self.ra_str or self.ra_str == '00h00m00.00s' or not self.validate_ra(self.ra_str):
            self.ra_box.config(bg='lightcoral')
            allow_run = False
        if not self.dec_str or self.dec_str == '00d00m00.0s' or not self.validate_dec(self.dec_str):
            self.dec_box.config(bg='lightcoral')
            allow_run=False
        # Mass, Age, Number Generated
        if self.mass == 0 or not self.validate_mass(self.mass):
            self.mass_box.config(bg='lightcoral')
            allow_run=False
        if self.num_generated == 0 or not self.validate_num(self.num_generated):
            self.number_box.config(bg='lightcoral')
            allow_run = False
        # Prefix
        if self.__prefix.get() == '':
            self.prefix_box.config(bg='lightcoral')
            allow_run = False

        if not allow_run:
            message = message + '\nA required field is missing or invalid. Please check any fields marked in red and try again.\n'
            tk.messagebox.showinfo('Run Message', message)
            return -1
        else:
            return 0

    def gui_print(self, new_text):
        # need to change the text in label to include this message on the end
        text_width = self.display_box.cget('width')
        if new_text == 'clc':
            self.display_box.delete('1.0', 'end')
        else:
            wrap_new_text = wrap(new_text, text_width)
            new_message = ''
            for i in range(0, len(wrap_new_text)):
                new_message += wrap_new_text[i] + '\n'
            self.display_box.insert('end', new_message)
        self.root.update()
        return

    def update_status(self, new_status):
        self.status_box_text.set(('Status: ' + new_status))
        self.gui_print(('\nStatus: '+ new_status))
        print(('\nStatus: '+ new_status))
        self.root.update()

    # Analysis Options validation functions

    def validate_resolution(self, new_text):
        # Check that number generated is an integer greater than zero
        if not new_text:
            # box cleared
            self.__resolution = 0.
            self.__resolution_box.config(bg='white')
            return True
        try:
            n = float(new_text)
            if n > 0.:
                self.__resolution = n
                self.__resolution_box.config(bg='white')
                return True
            else:
                self.__resolution_box.config(bg='lightcoral')
                return False
        except:
            self.__resolution_box.config(bg='lightcoral')
            return False

    # Stellar Info validation functions

    def validate_ra(self, new_text):
        # Check that RA makes sense
        if not new_text:
            # box cleared
            self.ra_str= '00h00m00.00s'
            self.ra_box.config(bg='white')
            return True
        try:
            hours = new_text[0:2]
            minutes = new_text[3:5]
            seconds = new_text[6:-2]
            if new_text[2] == 'h' and new_text[5] == 'm' and new_text[-1] == 's':
                hours = int(hours)
                minutes = int(minutes)
                seconds = float(seconds)
                if hours < 24 and minutes < 60 and seconds < 60:
                    self.ra_str = new_text
                    self.ra_box.config(bg='white')
                    return True
                else:
                    self.ra_str = new_text
                    self.ra_box.config(bg='lightcoral')
                    return False
            else:
                self.ra_str = new_text
                self.ra_box.config(bg='lightcoral')
                return False
        except:
            self.ra_str = new_text
            self.ra_box.config(bg='lightcoral')
            return False

    def validate_dec(self, new_text):
        # Check that DEC makes sense
        if not new_text:
            # box cleared
            self.dec_str = '00d00m00.0s'
            self.dec_box.config(bg='white')
            return True
        try:
            sign = new_text[0]
            degrees = new_text[1:3]
            minutes = new_text[4:6]
            seconds = new_text[7:-2]
            if new_text[3] == 'd' and new_text[6] == 'm' and new_text[-1] == 's':
                degrees = int(degrees)
                minutes = int(minutes)
                seconds = float(seconds)
                if degrees <= 90 and minutes < 60 and seconds < 60:
                    if sign == '+' or sign == '-':
                        self.dec_str = new_text
                        self.dec_box.config(bg='white')
                        return True
                    else:
                        self.dec_box.config(bg='lightcoral')
                        return False
                else:
                    self.dec_str = new_text
                    self.dec_box.config(bg='lightcoral')
                    return False
            else:
                self.dec_str = new_text
                self.dec_box.config(bg='lightcoral')
                return False
        except:
            self.dec_str = new_text
            self.dec_box.config(bg='lightcoral')
            return False

    def validate_num(self, new_text):
        # Check that number generated is an integer greater than zero
        if not new_text:
            # box cleared
            self.num_generated = 0
            self.number_box.config(bg='white')
            return True
        try:
            n = int(new_text)
            if n > 0:
                self.num_generated = n
                self.number_box.config(bg='white')
                return True
            else:
                self.number_box.config(bg='lightcoral')
                return False
        except:
            self.number_box.config(bg='lightcoral')
            return False

    def validate_mass(self, new_text):
        # Check that mass is a number greater than zero
        if not new_text:
            # box cleared
            self.mass = 0
            self.mass_box.config(bg='white')
            return True
        try:
            m = float(new_text)
            if m > 0:
                self.mass = m
                self.mass_box.config(bg='white')
                return True
            else:
                self.mass_box.config(bg='lightcoral')
                return False
        except:
            self.mass_box.config(bg='lightcoral')
            return False

    def validate_age(self, new_text):
        # Check that age is a number greater than zero
        if not new_text:
            # box cleared
            self.age = 5
            self.age_box.insert(-1, '5')
            self.age_box.config(bg='white', fg='gray')
            return True
        try:
            a = float(new_text)
            if a > 0:
                self.age = a
                self.age_box.config(bg='white', fg='black')
                return True
            else:
                self.age_box.config(bg='lightcoral', fg='black')
                return False
        except:
            self.age_box.config(bg='lightcoral', fg='black')
            return False

    def validate_jitter(self, new_text):
        # Check that jitter is a number greater than zero
        if not new_text:
            # box cleared
            self.added_jitter = 20
            self.added_jitter_box.config(bg='white', fg='gray')
            self.added_jitter_box.insert(-1, '20')
            return True
        try:
            j = float(new_text)
            if j >= 0:
                self.added_jitter = j
                self.added_jitter_box.config(bg='white', fg='black')
                return True
            else:
                self.added_jitter_box.config(bg='lightcoral', fg='black')
                return False
        except:
            self.added_jitter_box.config(bg='lightcoral', fg='black')
            return False

    def validate_floor(self, new_text):
        # Check that jitter is a number greater than zero
        if not new_text:
            # box cleared
            self.rv_floor = 20
            self.rv_floor_box.config(bg='white', fg='gray')
            self.rv_floor_box.insert(-1, '20')
            return True
        try:
            j = float(new_text)
            if j >= 0:
                self.rv_floor = j
                self.rv_floor_box.config(bg='white', fg='black')
                return True
            else:
                self.rv_floor_box.config(bg='lightcoral', fg='black')
                return False
        except:
            self.rv_floor_box.config(bg='lightcoral', fg='black')
            return False

    # Advanced Options validation functions
    def validate_P(self, new_text, box_name):
        # Period must be greater than, or equal to 0.1, if accepted, the other period boxes should be grayed out
        if not new_text:
            # Text box cleared
            if box_name.endswith('entry'):
                self.__P_fixed = None
                self.P1_box.config(bg='white')
                self.P2_box.config(state='normal')
                self.P3_box.config(state='normal')
            elif box_name.endswith('entry2'):
                self.__P_min = None
                self.P2_box.config(bg='white')
                if self.__P_max is None:
                    self.P1_box.config(state='normal')
            elif box_name.endswith('entry3'):
                self.__P_max = None
                self.P3_box.config(bg='white')
                if self.__P_min is None:
                    self.P1_box.config(state='normal')
            return True
        try:
            p = float(new_text)
            if p >= 0.1:  # If valid set the desired limit to the given value, and disable other boxes as needed
                if box_name.endswith('entry'):
                    self.__P_fixed = p
                    self.P1_box.config(bg='white')
                    self.P2_box.config(state='disabled')
                    self.P3_box.config(state='disabled')
                elif box_name.endswith('entry2'):
                    self.__P_min = p
                    self.P2_box.config(bg='white')
                    self.P1_box.config(state='disabled')
                elif box_name.endswith('entry3'):
                    self.__P_max = p
                    self.P3_box.config(bg='white')
                    self.P1_box.config(state='disabled')
                return True
            else:  # If invalid set the desired limit to its default value, color and disable other boxes as needed
                if box_name.endswith('entry'):
                    self.P1_box.config(bg='lightcoral')
                    self.P2_box.config(state='disabled')
                    self.P3_box.config(state='disabled')
                elif box_name.endswith('entry2'):
                    self.__P_min = 0.1
                    self.P2_box.config(bg='lightcoral')
                    self.P1_box.config(state='disabled')
                elif box_name.endswith('entry3'):
                    self.__P_max = float('inf')
                    self.P3_box.config(bg='lightcoral')
                    self.P1_box.config(state='disabled')
                return False
        except:  # If very invalid set the desired limit to its default value, color and disable other boxes as needed
            if box_name.endswith('entry'):
                self.P1_box.config(bg='lightcoral')
                self.P2_box.config(state='disabled')
                self.P3_box.config(state='disabled')
            elif box_name.endswith('entry2'):
                self.__P_min = 0.1
                self.P2_box.config(bg='lightcoral')
                self.P1_box.config(state='disabled')
            elif box_name.endswith('entry3'):
                self.__P_max = float('inf')
                self.P3_box.config(bg='lightcoral')
                self.P1_box.config(state='disabled')
            return False

    def validate_i(self, new_text, box_name):
        # TODO: the GUI box says *cos*i, not sin(i), and the paper 
        # refers to cosi running from 0 to 1. This may need to be fixed.
        # For now I am assuming the parameter is cosi
        # sin(i) must be greater than -1 and less than 1, if accepted, the other boxes should be grayed out
        # a limit value of 'transit' can also be accepted
        if not new_text:
            # box cleared
            if box_name.endswith('entry4'):
                self.__i_fixed = None
                self.i1_box.config(bg='white')
                self.i2_box.config(state='normal')
                self.i3_box.config(state='normal')
            elif box_name.endswith('entry5'):
                self.__i_min = None
                self.i2_box.config(bg='white')
                if self.__i_max is None:
                    self.i1_box.config(state='normal')
            elif box_name.endswith('entry6'):
                self.__i_max = None
                self.i3_box.config(bg='white')
                if self.__i_min is None:
                    self.i1_box.config(state='normal')
            return True
        try:
            i = float(new_text)
            if -1 <= i <= 1:
                if box_name.endswith('entry4'):
                    self.__i_fixed = i
                    self.i1_box.config(bg='white')
                    self.i2_box.config(state='disabled')
                    self.i3_box.config(state='disabled')
                elif box_name.endswith('entry5'):
                    self.__i_min = i
                    self.i2_box.config(bg='white')
                    self.i1_box.config(state='disabled')
                elif box_name.endswith('entry6'):
                    self.__i_max = i
                    self.i3_box.config(bg='white')
                    self.i1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry4'):
                    self.i1_box.config(bg='lightcoral')
                    self.i2_box.config(state='disabled')
                    self.i3_box.config(state='disabled')
                elif box_name.endswith('entry5'):
                    self.__i_min = -1.
                    self.i2_box.config(bg='lightcoral')
                    self.i1_box.config(state='disabled')
                elif box_name.endswith('entry6'):
                    self.__i_max = 1.
                    self.i3_box.config(bg='lightcoral')
                    self.i1_box.config(state='disabled')
                return False
        except:
            if new_text == 'transit':
                if box_name.endswith('entry4'):
                    self.__i_fixed = new_text
                    self.i1_box.config(bg='white')
                    self.i2_box.config(state='disabled')
                    self.i3_box.config(state='disabled')
                elif box_name.endswith('entry5'):
                    self.__i_min = new_text
                    self.i2_box.config(bg='white')
                    self.i1_box.config(state='disabled')
                elif box_name.endswith('entry6'):
                    self.__i_max = new_text
                    self.i3_box.config(bg='white')
                    self.i1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry4'):
                    self.i1_box.config(bg='lightcoral')
                    self.i2_box.config(state='disabled')
                    self.i3_box.config(state='disabled')
                elif box_name.endswith('entry5'):
                    self.__i_min = -1
                    self.i2_box.config(bg='lightcoral')
                    self.i1_box.config(state='disabled')
                elif box_name.endswith('entry6'):
                    self.__i_max = 1
                    self.i3_box.config(bg='lightcoral')
                    self.i1_box.config(state='disabled')
                return False

    def validate_e(self, new_text, box_name):
        # e must be greater than 0 and less than 1, if accepted, the other boxes should be grayed out
        if not new_text:
            # Text box cleared
            if box_name.endswith('entry7'):
                self.__e_fixed = None
                self.e1_box.config(bg='white')
                self.e2_box.config(state='normal')
                self.e3_box.config(state='normal')
            elif box_name.endswith('entry8'):
                self.__e_min = None
                self.e2_box.config(bg='white')
                if self.__e_max is None:
                    self.e1_box.config(state='normal')
            elif box_name.endswith('entry9'):
                self.__e_max = None
                self.e3_box.config(bg='white')
                if self.__e_min is None:
                    self.e1_box.config(state='normal')
            return True
        try:
            e = float(new_text)
            if 0 <= e <= 1:  # If valid set the desired limit to the given value, and disable other boxes as needed
                if box_name.endswith('entry7'):
                    self.__e_fixed = e
                    self.e1_box.config(bg='white')
                    self.e2_box.config(state='disabled')
                    self.e3_box.config(state='disabled')
                elif box_name.endswith('entry8'):
                    self.__e_min = e
                    self.e2_box.config(bg='white')
                    self.e1_box.config(state='disabled')
                elif box_name.endswith('entry9'):
                    self.__e_max = e
                    self.e3_box.config(bg='white')
                    self.e1_box.config(state='disabled')
                return True
            else:  # If invalid set the desired limit to its default value, color and disable other boxes as needed
                if box_name.endswith('entry7'):
                    self.e1_box.config(bg='lightcoral')
                    self.e2_box.config(state='disabled')
                    self.e3_box.config(state='disabled')
                elif box_name.endswith('entry8'):
                    self.__e_min = 0.
                    self.e2_box.config(bg='lightcoral')
                    self.e1_box.config(state='disabled')
                elif box_name.endswith('entry9'):
                    self.__e_max = 1.
                    self.e3_box.config(bg='lightcoral')
                    self.e1_box.config(state='disabled')
                return False
        except: # If invalid set the desired limit to its default value, color and disable other boxes as needed
            if box_name.endswith('entry7'):
                self.e1_box.config(bg='lightcoral')
                self.e2_box.config(state='disabled')
                self.e3_box.config(state='disabled')
            elif box_name.endswith('entry8'):
                self.__e_min = 0.
                self.e2_box.config(bg='lightcoral')
                self.e1_box.config(state='disabled')
            elif box_name.endswith('entry9'):
                self.__e_max = 1.
                self.e3_box.config(bg='lightcoral')
                self.e1_box.config(state='disabled')
            return False

    def validate_arg_peri(self, new_text, box_name):
        # arg peri must be greater than 0 and less than pi, if accepted, the other boxes should be grayed out
        if not new_text:
            # box cleared
            if box_name.endswith('entry10'):
                self.__arg_peri_fixed = None
                self.w1_box.config(bg='white')
                self.w2_box.config(state='normal')
                self.w3_box.config(state='normal')
            elif box_name.endswith('entry11'):
                self.__arg_peri_min = None
                self.w2_box.config(bg='white')
                if self.__arg_peri_max is None:
                    self.w1_box.config(state='normal')
            elif box_name.endswith('entry12'):
                self.__arg_peri_max = None
                self.w3_box.config(bg='white')
                if self.__arg_peri_min is None:
                    self.w1_box.config(state='normal')
            return True
        try:
            w = float(new_text)
            if 0 <= w <= np.pi:
                if box_name.endswith('entry10'):
                    self.__arg_peri_fixed = w
                    self.w1_box.config(bg='white')
                    self.w2_box.config(state='disabled')
                    self.w3_box.config(state='disabled')
                elif box_name.endswith('entry11'):
                    self.__arg_peri_min = w
                    self.w2_box.config(bg='white')
                    self.w1_box.config(state='disabled')
                elif box_name.endswith('entry12'):
                    self.__arg_peri_max = w
                    self.w3_box.config(bg='white')
                    self.w1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry10'):
                    self.w1_box.config(bg='lightcoral')
                    self.w2_box.config(state='disabled')
                    self.w3_box.config(state='disabled')
                elif box_name.endswith('entry11'):
                    self.__arg_peri_min = 0.
                    self.w2_box.config(bg='lightcoral')
                    self.w1_box.config(state='disabled')
                elif box_name.endswith('entry12'):
                    self.__arg_peri_max = np.pi
                    self.w3_box.config(bg='lightcoral')
                    self.w1_box.config(state='disabled')
                return False
        except:
            if box_name.endswith('entry10'):
                self.w1_box.config(bg='lightcoral')
                self.w2_box.config(state='disabled')
                self.w3_box.config(state='disabled')
            elif box_name.endswith('entry11'):
                self.__arg_peri_min = 0.
                self.w2_box.config(bg='lightcoral')
                self.w1_box.config(state='disabled')
            elif box_name.endswith('entry12'):
                self.__arg_peri_max = np.pi
                self.w3_box.config(bg='lightcoral')
                self.w1_box.config(state='disabled')
            return False

    def validate_m(self, new_text, box_name):
        # *mass ratio q* must be greater than 0 and less than 1, 
        # if accepted, the other boxes should be grayed out
        if not new_text:
            # box cleared
            if box_name.endswith('entry13'):
                self.__m_fixed = None
                self.m1_box.config(bg='white')
                self.m2_box.config(state='normal')
                self.m3_box.config(state='normal')
            elif box_name.endswith('entry14'):
                self.__m_min = None
                self.m2_box.config(bg='white')
                if self.__m_max is None:
                    self.m1_box.config(state='normal')
            elif box_name.endswith('entry15'):
                self.__m_max = None
                self.m3_box.config(bg='white')
                if self.__m_min is None:
                    self.m1_box.config(state='normal')
            return True
        try:
            m = float(new_text)
            if 0 <= m <= 1:
                if box_name.endswith('entry13'):
                    self.__m_fixed = m
                    self.m1_box.config(bg='white')
                    self.m2_box.config(state='disabled')
                    self.m3_box.config(state='disabled')
                elif box_name.endswith('entry14'):
                    self.__m_min = m
                    self.m2_box.config(bg='white')
                    self.m1_box.config(state='disabled')
                elif box_name.endswith('entry15'):
                    self.__m_max = m
                    self.m3_box.config(bg='white')
                    self.m1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry13'):
                    self.m1_box.config(bg='lightcoral')
                    self.m2_box.config(state='disabled')
                    self.m3_box.config(state='disabled')
                elif box_name.endswith('entry14'):
                    self.__m_min = 0.
                    self.m2_box.config(bg='lightcoral')
                    self.m1_box.config(state='disabled')
                elif box_name.endswith('entry15'):
                    self.__m_max = 1.
                    self.m3_box.config(bg='lightcoral')
                    self.m1_box.config(state='disabled')
                return False
        except:
            if box_name.endswith('entry13'):
                self.m1_box.config(bg='lightcoral')
                self.m2_box.config(state='disabled')
                self.m3_box.config(state='disabled')
            elif box_name.endswith('entry14'):
                self.__m_min = 0.
                self.m2_box.config(bg='lightcoral')
                self.m1_box.config(state='disabled')
            elif box_name.endswith('entry15'):
                self.__m_max = 1.
                self.m3_box.config(bg='lightcoral')
                self.m1_box.config(state='disabled')
            return False

    def validate_a(self, new_text, box_name):
        # a must be greater than 0, if accepted, the other boxes should be grayed out
        if not new_text:
            # box cleared
            if box_name.endswith('entry16'):
                self.__a_fixed = None
                self.a1_box.config(bg='white')
                self.a2_box.config(state='normal')
                self.a3_box.config(state='normal')
            elif box_name.endswith('entry17'):
                self.__a_min = None
                self.a2_box.config(bg='white')
                if self.__a_max is None:
                    self.a1_box.config(state='normal')
            elif box_name.endswith('entry18'):
                self.__a_max = None
                self.a3_box.config(bg='white')
                if self.__a_min is None:
                    self.a1_box.config(state='normal')
            return True
        try:
            a = float(new_text)
            if a >= 0:
                if box_name.endswith('entry16'):
                    self.__a_fixed = a
                    self.a1_box.config(bg='white')
                    self.a2_box.config(state='disabled')
                    self.a3_box.config(state='disabled')
                elif box_name.endswith('entry17'):
                    self.__a_min = a
                    self.a2_box.config(bg='white')
                    self.a1_box.config(state='disabled')
                elif box_name.endswith('entry18'):
                    self.__a_max = a
                    self.a3_box.config(bg='white')
                    self.a1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry16'):
                    self.a1_box.config(bg='lightcoral')
                    self.a2_box.config(state='disabled')
                    self.a3_box.config(state='disabled')
                elif box_name.endswith('entry17'):
                    self.__a_min = 0.
                    self.a2_box.config(bg='lightcoral')
                    self.a1_box.config(state='disabled')
                elif box_name.endswith('entry18'):
                    self.__a_max = float('inf')
                    self.a3_box.config(bg='lightcoral')
                    self.a1_box.config(state='disabled')
                return False
        except:
            if box_name.endswith('entry16'):
                self.a1_box.config(bg='lightcoral')
                self.a2_box.config(state='disabled')
                self.a3_box.config(state='disabled')
            elif box_name.endswith('entry17'):
                self.__a_min = 0.
                self.a2_box.config(bg='lightcoral')
                self.a1_box.config(state='disabled')
            elif box_name.endswith('entry18'):
                self.__a_max = float('inf')
                self.a3_box.config(bg='lightcoral')
                self.a1_box.config(state='disabled')
            return False

    def validate_phi(self, new_text, box_name):
        # phi must be greater than 0 and less than 2pi, if accepted, the other boxes should be grayed out
        if not new_text:
            # box cleared
            if box_name.endswith('entry19'):
                self.__phi_fixed = None
                self.phi1_box.config(bg='white')
                self.phi2_box.config(state='normal')
                self.phi3_box.config(state='normal')
            elif box_name.endswith('entry20'):
                self.__phi_min = None
                self.phi2_box.config(bg='white')
                if self.__phi_max is None:
                    self.phi1_box.config(state='normal')
            elif box_name.endswith('entry21'):
                self.__phi_max = None
                self.phi3_box.config(bg='white')
                if self.__phi_min is None:
                    self.phi1_box.config(state='normal')
            return True
        try:
            phi = float(new_text)
            if 0 <= phi <= 2*np.pi:
                if box_name.endswith('entry19'):
                    self.__phi_fixed = phi
                    self.phi1_box.config(bg='white')
                    self.phi2_box.config(state='disabled')
                    self.phi3_box.config(state='disabled')
                elif box_name.endswith('entry20'):
                    self.__phi_min = phi
                    self.phi2_box.config(bg='white')
                    self.phi1_box.config(state='disabled')
                elif box_name.endswith('entry21'):
                    self.__phi_max = phi
                    self.phi3_box.config(bg='white')
                    self.phi1_box.config(state='disabled')
                return True
            else:
                if box_name.endswith('entry19'):
                    self.__phi_fixed = phi
                    self.phi1_box.config(bg='lightcoral')
                    self.phi2_box.config(state='disabled')
                    self.phi3_box.config(state='disabled')
                elif box_name.endswith('entry20'):
                    self.__phi_min = phi
                    self.phi2_box.config(bg='lightcoral')
                    self.phi1_box.config(state='disabled')
                elif box_name.endswith('entry21'):
                    self.__phi_max = phi
                    self.phi3_box.config(bg='lightcoral')
                    self.phi1_box.config(state='disabled')
                return False
        except:
            if box_name.endswith('entry19'):
                self.phi1_box.config(bg='lightcoral')
                self.phi2_box.config(state='disabled')
                self.phi3_box.config(state='disabled')
            elif box_name.endswith('entry20'):
                self.__phi_min = 0.
                self.phi2_box.config(bg='lightcoral')
                self.phi1_box.config(state='disabled')
            elif box_name.endswith('entry21'):
                self.__phi_max = 2*np.pi
                self.phi3_box.config(bg='lightcoral')
                self.phi1_box.config(state='disabled')
            return False

    def validate_pd_mu(self, new_text):
        # Check that mu is a number greater than zero
        if not new_text:
            # box cleared
            self.__pd_mu = 5.03
            self.PD_mu_box.insert(-1, '5.03')
            self.PD_mu_box.config(bg='white', fg='gray')
            return True
        try:
            j = float(new_text)
            if j >= 0:
                self.__pd_mu = j
                self.PD_mu_box.config(bg='white', fg='black')
                return True
            else:
                self.PD_mu_box.config(bg='lightcoral', fg='black')
                return False
        except:
            self.PD_mu_box.config(bg='lightcoral', fg='black')
            return False

    def validate_pd_sig(self, new_text):
        # Check that sigma is a number greater than zero
        if not new_text:
            # box cleared
            self.__pd_sig = 2.28
            self.PD_sig_box.config(bg='white', fg='gray')
            self.PD_sig_box.insert(-1, '2.28')
            return True
        try:
            j = float(new_text)
            if j >= 0:
                self.__pd_sig = j
                self.PD_sig_box.config(bg='white', fg='black')
                return True
            else:
                self.PD_sig_box.config(bg='lightcoral', fg='black')
                return False
        except:
            self.PD_sig_box.config(bg='lightcoral', fg='black')
            return False

    def validate_mass_exp(self, new_text):
        # Check that gamma is a number
        if not new_text:
            # box cleared
            self.__q_exp = 0.
            self.mass_exp_box.config(bg='white', fg='gray')
            self.mass_exp_box.insert(-1, '0.0')
            return True
        try:
            j = float(new_text)
            self.__q_exp = j
            self.mass_exp_box.config(bg='white', fg='black')
            return True
        except:
            self.mass_exp_box.config(bg='lightcoral', fg='black')
            return False

    # Help button functions
    @ staticmethod
    def help_age():
        message = """Optional. Estimate of the age of the star. This will effect the modeled magnitude of the stars. 
        If none is entered an age of 5 Gyr will be assumed."""
        tk.messagebox.showinfo('Star Age Help', message)

    @staticmethod
    def help_coord():
        message = 'The right ascension (RA) or declination (Dec) of the star. Use the format 00h00m00.00s for RA and +00d00m00.0s for Dec'
        tk.messagebox.showinfo('Coordinates Help', message)

    @staticmethod
    def help_added_jitter():
        message = 'Optional. A jitter term to be added in quadrature to the measurement error.'
        tk.messagebox.showinfo('Added Jitter Help', message)

    @staticmethod
    def help_floor():
        message = """Optional. A lower limit on the radial velocity semi-amplitude. All generated companions with semi-amplitudes lower than this floor will not be rejected by the RV test. See documentation for more detail"""
        tk.messagebox.showinfo('RV Floor Help', message)

    @staticmethod
    def help_advanced():
        message = 'Alter possible limits of orbit parameters for generated companions. All values must be numbers, ' \
                  'except for limits on cos(i), any one of which may be "transit", to limit to only transiting companions'
        tk.messagebox.showinfo('Advanced Help', message)

    @staticmethod
    def help_period():
        message = 'The default period distribution is a Log-Normal distribution with \u03BC=5.03 and \u03C3=2.28. '\
                  'To use a different log normal distribution enter the mean and standard deviation here.'
        tk.messagebox.showinfo('Period Distribution Help', message)

    @staticmethod
    def help_mass():
        message = 'The default mass ratio distribution is a uniform distribution, a.k.a., a power-law distribution'\
                'with \u03B3=0.0. To use a different power law distribution enter the desired exponent here.'
        tk.messagebox.showinfo('Mass Ratio Distribution Help', message)

    # Run and Close Functions
    def run_code(self):
        self.run_button.focus()
        go = self.check_completeness()
        if go == 0:
            self.parent.run()

    def close(self):
        self.root.quit()

    # Access Functions
    def get_ao_filename(self):
        if self.__ao_check.get():
            file_list = [x.get() for x in self.ao_file_boxes]
            return file_list
        else:
            return ['']

    def get_rv_filename(self):
        if self.__rv_check.get()==1:
            return self.__rv_file.get()
        else:
            return ''

    def get_ruwe(self):
        return self.__ruwe_check.get()

    def get_gaia(self):
        return self.__gaia_check.get()

    def get_filter(self):
        filter_list = [x.get() for x in self.ao_filter_menus]
        return filter_list

    def get_resolution(self):
        return self.__resolution

    def get_prefix(self):
        return self.__prefix.get()

    def get_extra(self):
        if self.__extra.get() == 1:
            return True
        else:
            return False

    def get_all_out_bool(self):
        if self.__all_out.get() == 1:
            return True
        else:
            return False

    def get_limits(self):
        # puts all the limit information together in a dictionary
        limits = {}
        limits["P"] = OrderedDict()
        limits["P"]["fixed"] = self.__P_fixed
        limits["P"]["min"] = self.__P_min
        limits["P"]["max"] = self.__P_max

        limits["cos_i"] = OrderedDict()
        limits["cos_i"]["fixed"] = self.__i_fixed
        limits["cos_i"]["min"] = self.__i_min
        limits["cos_i"]["max"] = self.__i_max
        
        limits["ecc"] = OrderedDict()
        limits["ecc"]["fixed"] = self.__e_fixed
        limits["ecc"]["min"] = self.__e_min
        limits["ecc"]["max"] = self.__e_max

        limits["arg_peri"] = OrderedDict()
        limits["arg_peri"]["fixed"] = self.__arg_peri_fixed
        limits["arg_peri"]["min"] = self.__arg_peri_min
        limits["arg_peri"]["max"] = self.__arg_peri_max

        limits["q"] = OrderedDict()
        limits["q"]["fixed"] = self.__m_fixed
        limits["q"]["min"] = self.__m_min
        limits["q"]["max"] = self.__m_max

        limits["a"] = OrderedDict()
        limits["a"]["fixed"] = self.__a_fixed
        limits["a"]["min"] = self.__a_min
        limits["a"]["max"] = self.__a_max

        limits["phase"] = OrderedDict()
        limits["phase"]["fixed"] = self.__phi_fixed
        limits["phase"]["min"] = self.__phi_min
        limits["phase"]["max"] = self.__phi_max

        return limits

    def get_p_dist(self):
        return self.__pd_mu, self.__pd_sig

    def get_q_dist(self):
        return self.__q_exp

    def get_gaia_limit(self):
        return self.__gaia_limit.get()