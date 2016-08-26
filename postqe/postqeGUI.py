#!/usr/bin/env python3
#encoding: UTF-8


import wx
from wxPlot1DDialog import Plot1DDialog, Plot2DDialog
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt

import time, sys, os
import numpy as np
from readutils import read_line, read_n_real_numbers,\
read_charge_file_iotk, read_charge_file_hdf5, read_wavefunction_file_iotk,\
read_wavefunction_file_hdf5, read_pp_out_file, write_charge, create_header
from compute_vs import compute_v_bare, compute_v_h, compute_v_xc, compute_G
from celldm import calcola_celldm
from plot import plot1D, plot2D
from postqe import get_from_xml
import settings


class RedirectText:
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl
 
    def write(self,string):
        self.out.WriteText(string)


class MyApp(wx.App):
    def OnInit(self):
        frame = MyFrame("postqe", (50, 60), (450, 340))
        frame.Show()
        self.SetTopWindow(frame)
        return True
    
class MyFrame(wx.Frame):
    def __init__(self, title, pos, size):
                
        # Initialize some flags
        self.IsXmlRead = False
        self.IsChargeRead = False
        self.IsWavefunctionRead = False
        
        wx.Frame.__init__(self, None, -1, title, pos, size)                
        panel = wx.Panel(self, -1)
        panel.SetBackgroundColour("White")
        
        text = wx.TextCtrl(self, wx.ID_ANY, size=size, style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(text, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.Fit()
        
        self.Bind(wx.EVT_CLOSE, self.OnQuit)                             
        self.createMenuBar()       
        self.CreateStatusBar()        
        self.SetStatusText("Welcome to QE postprocessing!")
        
        # Redirect stout to the TextCtrl in the main panel
        self.redir=RedirectText(text)
        sys.stdout=self.redir
        sys.stderr=self.redir

    ############################################################################
    #
    # Menu functions
    #
    
    def menuData(self):
        """
        This function returns the menu items, including shortkeys and help strings.
        To modify the menus of the application you need to change them here.
        """
        return (("&File",
            ("&Open xml data file...\tctrl+X", "Load charge density data", self.OnOpenxml),
            ("Open &binary charge density file...", "Load charge density data\
            in binary format (platform dependent)", self.OnOpenBinaryCharge),
            ("Open HDF5 &charge density file...\tctrl+C", "Load charge density data\
            in HDF5 format", self.OnOpenHDF5Charge),
            ("Open binary &wavefunctions file...", "Load wavefunction\
            data in binary format (platform dependent)", self.OnOpenBinaryWavefunction),
            ("Open HDF5 &wavefunctions file...\tctrl+W", "Load wavefunction\
            data in HDF5 format", self.OnOpenHDF5Wavefunction),
            ("Save charge as text file\tctrl+S", "Save charge\
            data in text format", self.OnSaveCharge),
            ("Save last calculated quantity as text file\tctrl+R", "Save last\
            calculated quantity in text format", self.OnSaveLast),
            ("", "", ""),
            ("Quit...\tctrl+Q", "Exit the program", self.OnQuit)),
                ("&Compute",
            ("Compute &V_bare\tctrl+V", "Compute the bare potential from the\
            electronic charge", self.OnComputeVbare),
            ("Compute V_bare+V_&H\tctrl+X", "Compute the bare plus Hartree\
            potential from the electronic charge", self.OnComputeVbare_VH),
            ("Compute V_bare+V_&H+V_ex\tctrl+X", "Compute the bare plus Hartree\
            plus exchange-correlation potential from the electronic charge", self.OnComputeVbare_VH_Vxc)),
                ("&Plot",
            ("&1D Plot\tctrl+P", "Generate a 1D plot of the charge or other quantity", self.OnPlot1D),
            ("&2D Plot\tctrl+O", "Generate a 2D plot of the charge or other quantity", self.OnPlot2D)),
                ("&Help",
            ("&Help contents\tctrl+H", "Compute the bare potential from the electronic charge", self.OnHelp),
            ("&About\tctrl+A", "Compute the bare plus the Hartree potential from the electronic charge", self.OnAbout)))
        
        
    def createMenuBar(self):
        menuBar = wx.MenuBar()
        for eachMenuData in self.menuData():
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1:]
            menuBar.Append(self.createMenu(menuItems), menuLabel)
        self.SetMenuBar(menuBar)
        
        
    def createMenu(self, menuData):
        menu = wx.Menu()
        for eachLabel, eachStatus, eachHandler in menuData:
            if not eachLabel:
                menu.AppendSeparator()
                continue
            menuItem = menu.Append(-1, eachLabel, eachStatus)
            self.Bind(wx.EVT_MENU, eachHandler, menuItem)
        return menu


    ############################################################################
    #
    # Handling methods, rather self explaining
    #
    
    def OnOpenxml(self, event):
        wildcard =  "Xml file (*.xml)|*.xml|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Choose an xml QE file", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                self.ecutwfc, self.ecutrho, self.ibrav, self.alat, self.a, self.b,\
                self.functional, self.atomic_positions, self.atomic_species,\
                self.nat, self.ntyp = get_from_xml(fname,settings.schema)    
                self.celldms = calcola_celldm(self.alat,self.a[0],self.a[1],self.a[2],self.ibrav)
                self.IsXmlRead = True
                self.SetStatusText("Read xml file: "+fname)
            except:
                wx.MessageBox("Something went wrong while opening the xml file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
        
    def OnOpenBinaryCharge(self, event):
        wildcard =  "dat file (*.dat)|*.dat|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Choose a binary charge file (QE format)", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                self.charge = read_charge_file_iotk(fname)
                self.nr = self.charge.shape
                self.IsChargeRead = True
                self.SetStatusText("Read charge file: "+fname)
            except:
                wx.MessageBox("Something went wrong while opening the charge file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
        
    def OnOpenHDF5Charge(self, event):
        wildcard =  "HDF5 files (*.hdf5)|*.hdf5|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Choose an HDF5 charge file", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                print (fname)
                self.charge = read_charge_file_hdf5(fname)
                self.nr = self.charge.shape
                self.IsChargeRead = True
                self.SetStatusText("Read charge file: "+fname)
            except:
                wx.MessageBox("Something went wrong while opening the charge file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
        
    def OnOpenBinaryWavefunction(self, event):
        wx.MessageBox("This command is not implemented yet! Sorry.",
        "Not Implemented", wx.OK | wx.ICON_INFORMATION, self)
        
        
    def OnOpenHDF5Wavefunction(self, event):
        wildcard =  "HDF5 files (*.hdf5)|*.hdf5|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Choose a HDF5 wavefunctions file", os.getcwd(),"", wildcard)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetPath()
                self.wavefunctions = read_wavefunction_file_hdf5(fname)
                self.IsChargeRead = True
                self.SetStatusText("Read wavefunctions file: "+fname)
            except:
                wx.MessageBox("Something went wrong while opening the wavefunction file... not loaded.",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()

        
    def OnSaveCharge(self, event):
        if (not self.IsChargeRead):
            wx.MessageBox("No charge file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
            
        wildcard =  "txt file (*.txt)|*.txt|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Save charge file in text format...", os.getcwd(),"", wildcard, wx.FD_SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetFilename()
                dirname = dialog.GetDirectory()
                
                if (self.IsXmlRead):
                    header = create_header("",self.nr,self.ibrav,self.celldms,self.nat,\
                    self.ntyp,self.atomic_species,self.atomic_positions)
                else:
                    header = ""
                
                oldpath = os.getcwd()
                os.chdir(dirname)
                write_charge(fname,self.charge,header)
                os.chdir(oldpath)   
                
                self.SetStatusText("Charge written on text file: "+fname)
            except:
                wx.MessageBox("Something wrong while writing the charge file...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
    
    
    def OnSaveLast(self, event):
        try:
            self.vtot
        except:
            wx.MessageBox("No quantity calculated!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
            
        wildcard =  "txt file (*.txt)|*.txt|" \
                    "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Save last calculated quantity in text format...", os.getcwd(),"", wildcard, wx.FD_SAVE)
        if dialog.ShowModal() == wx.ID_OK:
            try:
                fname = dialog.GetFilename()
                dirname = dialog.GetDirectory()
                
                if (self.IsXmlRead):
                    header = create_header("",self.nr,self.ibrav,self.celldms,self.nat,\
                    self.ntyp,self.atomic_species,self.atomic_positions)
                else:
                    header = ""
                
                oldpath = os.getcwd()
                os.chdir(dirname)
                write_charge(fname,self.vtot,header)
                os.chdir(oldpath)   
                
                self.SetStatusText("Last calculated quantity written on text file: "+fname)
            except:
                wx.MessageBox("Something wrong while writing the file...",
                "", wx.OK | wx.ICON_EXCLAMATION, self)  

        dialog.Destroy()
        
        
    def OnComputeVbare(self, event):
        if (not self.IsXmlRead):
            wx.MessageBox("No xml file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        if (not self.IsChargeRead):
            wx.MessageBox("No charge file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
                
        try:
            self.v_bare = compute_v_bare(self.ecutrho, self.alat, self.a[0],\
            self.a[1], self.a[2], self.nr, self.atomic_positions, self.atomic_species,\
            settings.pseudodir)  
            self.vtot = self.v_bare
            self.SetStatusText("V bare calculated!")
        except:
            wx.MessageBox("Something wrong while computing the bare potential...",
            "", wx.OK | wx.ICON_EXCLAMATION, self) 
    
        
    def OnComputeVbare_VH(self, event):
        if (not self.IsXmlRead):
            wx.MessageBox("No xml file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        if (not self.IsChargeRead):
            wx.MessageBox("No charge file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
                
        try:
            self.v_bare = compute_v_bare(self.ecutrho, self.alat, self.a[0],\
            self.a[1], self.a[2], self.nr, self.atomic_positions, self.atomic_species,\
            settings.pseudodir)              
            self.v_h =  compute_v_h(self.charge,self.ecutrho,self.alat,self.b)
            self.vtot = self.v_bare + self.v_h
            self.SetStatusText("Bare and Hartree potentials calculated!")
        except:
            wx.MessageBox("Something wrong while computing the bare and Hartree potentials...",
            "", wx.OK | wx.ICON_EXCLAMATION, self) 
    
        
    def OnComputeVbare_VH_Vxc(self, event):
        if (not self.IsXmlRead):
            wx.MessageBox("No xml file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        if (not self.IsChargeRead):
            wx.MessageBox("No charge file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
                
        try:
            self.v_bare = compute_v_bare(self.ecutrho, self.alat, self.a[0],\
            self.a[1], self.a[2], self.nr, self.atomic_positions, self.atomic_species,\
            settings.pseudodir)              
            self.v_h =  compute_v_h(self.charge,self.ecutrho,self.alat,self.b)
            self.charge_core = np.zeros(self.charge.shape)    # only for now, later in input
            self.v_xc = compute_v_xc(self.charge,self.charge_core,str(self.functional))
            self.vtot = self.v_bare + self.v_h + self.v_xc
            self.SetStatusText("Bare, Hartree and xc potentials calculated!")
        except:
            wx.MessageBox("Something wrong while computing the bare, Hartree and xc potentials...",
            "", wx.OK | wx.ICON_EXCLAMATION, self) 
    
    
    def OnPlot1D(self, event):
               
        if (not self.IsXmlRead):
            wx.MessageBox("No xml file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        if (not self.IsChargeRead):
            wx.MessageBox("No charge file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        try:
            # Create the Plot1D Dialog and get the x0,e1,nx values
            dlg = Plot1DDialog()
            result = dlg.ShowModal()
            dlg.Destroy()
            
            if result == wx.ID_OK:
                selection = int(dlg.radiobox.GetSelection())
                x0 = np.array([float(dlg.x0_0.GetValue()), float(dlg.x0_1.GetValue()), float(dlg.x0_2.GetValue())])
                e1 = np.array([float(dlg.e1_0.GetValue()), float(dlg.e1_1.GetValue()),float(dlg.e1_2.GetValue())])
                nx = int(dlg.nx.GetValue())

                # Compute the G vectors if needed
                self.SetStatusText("Generating the g vectors...")
                try:                
                    self.G
                except:
                    self.G = compute_G(self.b,self.charge.shape,self.ecutrho,self.alat)
                
                
                if (selection==0):
                    try:
                        toplot = self.charge
                    except:
                        wx.MessageBox("No charge file loaded!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==1):
                    try:
                        toplot = self.v_bare
                    except:
                        wx.MessageBox("No Vbare calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==2):
                    try:
                        toplot = self.v_bare + self.v_h  
                    except:
                        wx.MessageBox("No Vbare+VH calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==3):
                    try:
                        toplot = self.v_bare + self.v_h + self.v_xc 
                    except:
                        wx.MessageBox("No Vbare+VH calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                
                self.SetStatusText("Plotting...")
                self.plot1D = plot1D(toplot,self.G,x0,e1,nx)
                self.plot1D.show()
                self.SetStatusText("Plotting... done!")
        except:
            wx.MessageBox("Something wrong while plotting the charge...",
            "", wx.OK | wx.ICON_EXCLAMATION, self) 
            self.SetStatusText("")
    
        
    def OnPlot2D(self, event):
        if (not self.IsXmlRead):
            wx.MessageBox("No xml file loaded!",
                "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
            return
        
        try:
            # Create the Plot1D Dialog and get the x0,e1,nx values
            dlg = Plot2DDialog()
            result = dlg.ShowModal()
            dlg.Destroy()
            
            if result == wx.ID_OK:
                selection = int(dlg.radiobox.GetSelection())
                x0 = np.array([float(dlg.x0_0.GetValue()), float(dlg.x0_1.GetValue()), float(dlg.x0_2.GetValue())])
                e1 = np.array([float(dlg.e1_0.GetValue()), float(dlg.e1_1.GetValue()),float(dlg.e1_2.GetValue())])
                nx = int(dlg.nx.GetValue())
                e2 = np.array([float(dlg.e2_0.GetValue()), float(dlg.e2_1.GetValue()),float(dlg.e2_2.GetValue())])
                ny = int(dlg.ny.GetValue())
                dlg.Destroy()
                
                # Compute the G vectors if needed
                self.SetStatusText("Generating the g vectors...")
                try:                
                    self.G
                except:
                    self.G = compute_G(self.b,self.charge.shape,self.ecutrho,self.alat)
                
                if (selection==0):
                    try:
                        toplot = self.charge
                    except:
                        wx.MessageBox("No charge file loaded!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==1):
                    try:
                        toplot = self.v_bare
                    except:
                        wx.MessageBox("No Vbare calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==2):
                    try:
                        toplot = self.v_bare + self.v_h  
                    except:
                        wx.MessageBox("No Vbare+VH calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                elif (selection==3):
                    try:
                        toplot = self.v_bare + self.v_h + self.v_xc 
                    except:
                        wx.MessageBox("No Vbare+VH calculated!",
                        "Warning!", wx.OK | wx.ICON_EXCLAMATION, self)  
                        return
                                       
                self.SetStatusText("Plotting...")
                self.plot2D = plot2D(toplot,self.G,x0,e1,e2,nx,ny,True)
                self.plot2D.show()
                self.SetStatusText("Plotting... done!")
        except:
            wx.MessageBox("Something wrong while plotting the charge...",
            "", wx.OK | wx.ICON_EXCLAMATION, self) 
            self.SetStatusText("")
        
    
    def OnQuit(self, event):
        plt.close("all")  # close all matplotlib figures
        self.Destroy()

                
    def OnHelp(self, event):
        wx.MessageBox("This is the post processing tool of Quantum Espresso",
        "Help Contents", wx.OK | wx.ICON_INFORMATION, self)


    def OnAbout(self, event):
        wx.MessageBox("This is the post processing tool of Quantum Espresso",
        "About postqe", wx.OK | wx.ICON_INFORMATION, self)      
        

if __name__ == "__main__":
    app = MyApp(False)
    app.MainLoop()
