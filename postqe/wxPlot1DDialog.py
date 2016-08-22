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
    def __init__(self):
        wx.Dialog.__init__(self, None, -1, "Plot 1D Charge")
        
        # Create the text controls
        notes = wx.StaticText(self, -1, "Please enter the input paramters below: ")
        x0_l = wx.StaticText(self, -1, "Starting point (x,y,z): ")
        self.x0_0 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.x0_1 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.x0_2 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        e1_l = wx.StaticText(self, -1, "Direction (x,y,z): ")
        self.e1_0 = wx.TextCtrl(self, -1, "1.0", validator = FloatValidator())
        self.e1_1 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        self.e1_2 = wx.TextCtrl(self, -1, "0.0", validator = FloatValidator())
        nx_l = wx.StaticText(self, -1, "Number of points along above direction (nx): ")
        self.nx = wx.TextCtrl(self, -1, "20", validator = PosIntValidator())

        # Use standard button IDs
        okay = wx.Button(self, wx.ID_OK)
        okay.SetDefault()
        cancel = wx.Button(self, wx.ID_CANCEL)

        # Layout with sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(notes, 0, wx.ALL, 5)
        sizer.Add(wx.StaticLine(self), 0, wx.EXPAND|wx.ALL, 5)

        fgs = wx.FlexGridSizer(4, 4, 5, 5)
        fgs.Add(x0_l, 0, wx.ALIGN_RIGHT)
        fgs.Add(self.x0_0, 0, wx.EXPAND)
        fgs.Add(self.x0_1, 0, wx.EXPAND)
        fgs.Add(self.x0_2, 0, wx.EXPAND)
        fgs.Add(e1_l, 0, wx.ALIGN_RIGHT)
        fgs.Add(self.e1_0, 0, wx.EXPAND)
        fgs.Add(self.e1_1, 0, wx.EXPAND)
        fgs.Add(self.e1_2, 0, wx.EXPAND)
        fgs.Add(nx_l, 0, wx.ALIGN_RIGHT)
        fgs.Add(self.nx, 0, wx.EXPAND)
        fgs.AddGrowableCol(1)
        
        sizer.Add(fgs, 0, wx.EXPAND|wx.ALL, 5)

        btns = wx.StdDialogButtonSizer()
        btns.AddButton(okay)
        btns.AddButton(cancel)
        btns.Realize()
        
        sizer.Add(btns, 0, wx.EXPAND|wx.ALL, 5)
        self.SetSizer(sizer)
        sizer.Fit(self)
    

