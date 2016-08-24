#!/usr/bin/env python3
#encoding: UTF-8

import wx

class PosIntValidator(wx.PyValidator):
    def __init__(self):
        wx.PyValidator.__init__(self)

    def Clone(self):
        """
        Note that every validator must implement the Clone() method.
        """
        return PosIntValidator()

    def Validate(self, win):
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()
        try:
            int(text)
            if (int(text)<1):
                raise
            textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        except:
            wx.MessageBox("This field must contain a valid integer number!", "Error")
            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True


class FloatValidator(wx.PyValidator):
    def __init__(self):
        wx.PyValidator.__init__(self)

    def Clone(self):
        """
        Note that every validator must implement the Clone() method.
        """
        return FloatValidator()

    def Validate(self, win):
        textCtrl = self.GetWindow()
        text = textCtrl.GetValue()
        try:
            float(text)
            textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
        except:
            wx.MessageBox("This field must contain a valid real number!", "Error")
            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False

    def TransferToWindow(self):
        return True

    def TransferFromWindow(self):
        return True
    

class Plot1DDialog(wx.Dialog):
    def __init__(self,texttitle="Plot 1D Charge"):
        wx.Dialog.__init__(self, None, -1, texttitle)
        
        # Plot what?
        choicesList = ['charge','Vbare','Vbare+VH','Vtot']
        self.radiobox = wx.RadioBox(self, -1, "", (450, 10), wx.DefaultSize, choicesList, 4, wx.RA_SPECIFY_COLS)
        
        # Create the text controls
        self.notes = wx.StaticText(self, -1, "Please enter the input paramters below: ")
        self.x0_l = wx.StaticText(self, -1, "Starting point (x,y,z): ")
        self.x0_0 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.x0_1 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.x0_2 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.e1_l = wx.StaticText(self, -1, "Direction (x,y,z): ")
        self.e1_0 = wx.TextCtrl(self, -1, "1.0", validator = FloatValidator())
        self.e1_1 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.e1_2 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.nx_l = wx.StaticText(self, -1, "Number of points along above direction (nx): ")
        self.nx = wx.TextCtrl(self, -1, "20", validator = PosIntValidator())

        # Use standard button IDs
        self.okay = wx.Button(self, wx.ID_OK)
        self.okay.SetDefault()
        self.cancel = wx.Button(self, wx.ID_CANCEL)

        # Layout with sizers
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.radiobox, 0, wx.ALL, 5)
        self.sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.ALL, 5)
        self.sizer.Add(self.notes, 0, wx.ALL, 5)
        self.sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.ALL, 5)

        self.fgs = wx.FlexGridSizer(4, 4, 5, 5)
        self.fgs.Add(self.x0_l, 0, wx.ALIGN_RIGHT)
        self.fgs.Add(self.x0_0, 0, wx.EXPAND)
        self.fgs.Add(self.x0_1, 0, wx.EXPAND)
        self.fgs.Add(self.x0_2, 0, wx.EXPAND)
        self.fgs.Add(self.e1_l, 0, wx.ALIGN_RIGHT)
        self.fgs.Add(self.e1_0, 0, wx.EXPAND)
        self.fgs.Add(self.e1_1, 0, wx.EXPAND)
        self.fgs.Add(self.e1_2, 0, wx.EXPAND)
        self.fgs.Add(self.nx_l, 0, wx.ALIGN_RIGHT)
        self.fgs.Add(self.nx, 0, wx.EXPAND)
        self.fgs.AddGrowableCol(1)
        
        self.sizer.Add(self.fgs, 0, wx.EXPAND|wx.ALL, 5)

        btns = wx.StdDialogButtonSizer()
        btns.AddButton(self.okay)
        btns.AddButton(self.cancel)
        btns.Realize()
        
        self.sizer.Add(btns, 0, wx.EXPAND|wx.ALL, 5)
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)
    

class Plot2DDialog(Plot1DDialog):
    def __init__(self,texttitle="Plot 2D Charge"):
        Plot1DDialog.__init__(self, texttitle)
        
        self.e2_l = wx.StaticText(self, -1, "Direction 2 (x,y,z): ")
        self.e2_0 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.e2_1 = wx.TextCtrl(self, -1, "1.0", validator = FloatValidator())
        self.e2_2 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.ny_l = wx.StaticText(self, -1, "Number of points along above direction (nx): ")
        self.ny = wx.TextCtrl(self, -1, "20", validator = PosIntValidator())
        
        self.fgs2 = wx.FlexGridSizer(2, 4, 5, 5)
        self.fgs2.Add(self.e2_l, 0, wx.ALIGN_RIGHT)
        self.fgs2.Add(self.e2_0, 0, wx.EXPAND)
        self.fgs2.Add(self.e2_1, 0, wx.EXPAND)
        self.fgs2.Add(self.e2_2, 0, wx.EXPAND)
        self.fgs2.Add(self.ny_l, 0, wx.ALIGN_RIGHT)
        self.fgs2.Add(self.ny, 0, wx.EXPAND)
        self.fgs2.AddGrowableCol(1)    

        self.sizer.Add(self.fgs2, 0, wx.EXPAND|wx.ALL, 5)
        
        btns = wx.StdDialogButtonSizer()
        btns.AddButton(self.okay)
        btns.AddButton(self.cancel)
        btns.Realize()
        
        self.sizer.Add(btns, 0, wx.EXPAND|wx.ALL, 5)
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)