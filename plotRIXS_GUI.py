"""
This demo demonstrates how to draw a dynamic mpl (matplotlib) 
plot in a wxPython application.
It allows "live" plotting as well as manual zooming to specific
regions.
Both X and Y axes allow "auto" or "manual" settings. For Y, auto
mode sets the scaling of the graph to see all the data points.
For X, auto mode makes the graph "follow" the data. Set it X min
to manual 0 to always see the whole data from the beginning.
Note: press Enter in the 'manual' text box to make a new value 
affect the plot.
Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 31.07.2008
"""
import os
import pprint
import random
import sys
import wx

# The recommended way to use wx with mpl is with the WXAgg
# backend. 
#
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import numpy as np
import pylab


class DataGen(object):
    """ A silly class that generates pseudo-random data for
        display in the plot.
    """
    def __init__(self, init=50):
        self.data = self.init = init
        
    def next(self):
        self._recalc_data()
        return self.data
    
    def _recalc_data(self):
        delta = random.uniform(-0.5, 0.5)
        r = random.random()

        if r > 0.9:
            self.data += delta * 15
        elif r > 0.8: 
            # attraction to the initial value
            delta += (0.5 if self.init > self.data else -0.5)
            self.data += delta
        else:
            self.data += delta


class BoundControlBox(wx.Panel):
    """ A static box with a couple of radio buttons and a text
        box. Allows to switch between an automatic mode and a 
        manual mode with an associated value.
    """
    def __init__(self, parent, ID, label, initval):
        wx.Panel.__init__(self, parent, ID)
        
        self.value = initval
        
        box = wx.StaticBox(self, -1, label)
        sizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        
        self.radio_auto = wx.RadioButton(self, -1, 
            label="Auto", style=wx.RB_GROUP)
        self.radio_manual = wx.RadioButton(self, -1,
            label="Manual")
        self.manual_text = wx.TextCtrl(self, -1, 
            size=(45,-1),
            value=str(initval),
            style=wx.TE_PROCESS_ENTER)
        
        self.Bind(wx.EVT_UPDATE_UI, self.on_update_manual_text, self.manual_text)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.manual_text)
        
        manual_box = wx.BoxSizer(wx.HORIZONTAL)
        manual_box.Add(self.radio_manual, flag=wx.ALIGN_CENTER_VERTICAL)
        manual_box.Add(self.manual_text, flag=wx.ALIGN_CENTER_VERTICAL)
        
        sizer.Add(self.radio_auto, 0, wx.ALL, 10)
        sizer.Add(manual_box, 0, wx.ALL, 10)
        
        self.SetSizer(sizer)
        sizer.Fit(self)
    
    def on_update_manual_text(self, event):
        self.manual_text.Enable(self.radio_manual.GetValue())
    
    def on_text_enter(self, event):
        self.value = self.manual_text.GetValue()
    
    def is_auto(self):
        return self.radio_auto.GetValue()
        
    def manual_value(self):
        return self.value


class GraphFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Demo: dynamic matplotlib graph'
    
    def __init__(self, seedname):
        wx.Frame.__init__(self, None, -1, self.title)
        
        #self.datagen = DataGen()
        #self.data = [self.datagen.next()]
        self.paused = True
        self.seedname = seedname
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)        
        self.redraw_timer.Start(100)

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
                
        self.menubar.Append(menu_file, "&File")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        self.panel = wx.Panel(self)

        self.init_plot()
        self.canvas = FigCanvas(self.panel, -1, self.fig)

        self.xmin_control = BoundControlBox(self.panel, -1, "X min", self.xmin)
        self.xmax_control = BoundControlBox(self.panel, -1, "X max", self.xmax)
        self.ymin_control = BoundControlBox(self.panel, -1, "Y min", self.ymin)
        self.ymax_control = BoundControlBox(self.panel, -1, "Y max", self.ymax)
        
        #self.pause_button = wx.Button(self.panel, -1, "Pause")
        #self.Bind(wx.EVT_BUTTON, self.on_pause_button, self.pause_button)
        #self.Bind(wx.EVT_UPDATE_UI, self.on_update_pause_button, self.pause_button)
        
        self.cb_grid = wx.CheckBox(self.panel, -1, 
            "Show Grid",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_grid, self.cb_grid)
        self.cb_grid.SetValue(True)
        
        self.cb_ylab = wx.CheckBox(self.panel, -1, 
            "Show Y labels",
            style=wx.ALIGN_RIGHT)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_ylab, self.cb_ylab)        
        self.cb_ylab.SetValue(True)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        #self.hbox1.Add(self.pause_button, border=5, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_grid, border=5, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.hbox1.AddSpacer(10)
        self.hbox1.Add(self.cb_ylab, border=5, flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.xmin_control, border=5, flag=wx.ALL)
        self.hbox2.Add(self.xmax_control, border=5, flag=wx.ALL)
        self.hbox2.AddSpacer(24)
        self.hbox2.Add(self.ymin_control, border=5, flag=wx.ALL)
        self.hbox2.Add(self.ymax_control, border=5, flag=wx.ALL)
        
        #self.pause_button = wx.Button(self.panel, -1, "Pause")
        #self.Bind(wx.EVT_BUTTON, self.on_pause_button, self.pause_button)
        #self.Bind(wx.EVT_UPDATE_UI, self.on_update_pause_button, self.pause_button)
        self.hbox3 = wx.BoxSizer(wx.HORIZONTAL)
        self.close_button = wx.Button(self.panel, 1, 'Close')
        self.Bind(wx.EVT_BUTTON, self.OnClose, id=1)
        self.hbox3.Add(self.close_button)

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, flag=wx.LEFT | wx.TOP | wx.GROW)        
        self.vbox.Add(self.hbox1, 0, flag=wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox2, 0, flag=wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox3, 0, flag=wx.ALIGN_CENTER | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def OnClose(self, event):
        self.Close()
    
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def init_plot(self):
        self.dpi = 100
        self.fig = Figure((3.0, 4.0), dpi=self.dpi)

        self.axes = self.fig.add_subplot(111)
        self.axes.set_axis_bgcolor('white')
        self.axes.set_title('L-edge RIXS Spectrum', size=12)
        self.axes.set_ylabel('Energy Loss (eV)', size=10)
        self.axes.set_xlabel('Incoming X-ray Energy (eV)', size=10)
        
        pylab.setp(self.axes.get_xticklabels(), fontsize=8)
        pylab.setp(self.axes.get_yticklabels(), fontsize=8)

        # plot the data as a line series, and save the reference 
        # to the plotted line series
        #

        self.data_nsp = np.loadtxt(self.seedname + '.RIXS-nsf.dat')
        self.data_sp = np.loadtxt(self.seedname + '.RIXS-sf.dat')

        self.xmin = min(self.data_nsp[:,0]) 
        self.xmax = max(self.data_nsp[:,0]) 
        self.ymin = min(self.data_nsp[:,1])
        self.ymax = max(self.data_nsp[:,1]) 
        print 'xmin and xmax ',self.xmin, self.xmax 
        print 'ymin and ymax ',self.ymin, self.ymax 

        divX = len(np.unique(self.data_nsp[:,0])) #incoming energy
        divY = len(np.unique(self.data_nsp[:,1])) #energy loss
        self.Z1 = self.data_nsp[:,2].reshape(divX,divY)
        self.Z2 = self.data_sp[:,2].reshape(divX,divY)
        self.Ztol = self.Z1 + self.Z2
        self.plot_data = self.axes.imshow(np.transpose(self.Ztol), interpolation='bilinear', \
                         extent=[self.xmin, self.xmax, self.ymin, self.ymax],
                         origin='lower',vmax=self.Ztol.max()/2.0, vmin=0, cmap = 'coolwarm')
        #plt.colorbar(orientation ='vertical')
        #plt.title('Total')

        #self.plot_data = self.axes.plot(
        #    self.dat[:,0], self.dat[:,1], 
        #    linewidth=1,
        #    color='blue',
        #    )[0]
        self.axes.set_xbound(lower=self.xmin, upper=self.xmax)
        self.axes.set_ybound(lower=self.ymin, upper=self.ymax)
        #self.plot_data = self.axes.plot(
        #    self.data, 
        #    linewidth=1,
        #    color=(1, 1, 0),
        #    )[0]

    def draw_plot(self):
        """ Redraws the plot
        """
        # when xmin is on auto, it "follows" xmax to produce a 
        # sliding window effect. therefore, xmin is assigned after
        # xmax.
        #
        if self.xmax_control.is_auto():
            xmax = self.xmax
        else:
            xmax = float(self.xmax_control.manual_value())
            
        if self.xmin_control.is_auto():            
            xmin = self.xmin
        else:
            xmin = float(self.xmin_control.manual_value())

        # for ymin and ymax, find the minimal and maximal values
        # in the data set and add a mininal margin.
        # 
        # note that it's easy to change this scheme to the 
        # minimal/maximal value in the current display, and not
        # the whole data set.
        # 
        if self.ymin_control.is_auto():
            ymin = self.ymin
        else:
            ymin = float(self.ymin_control.manual_value())
        
        if self.ymax_control.is_auto():
            ymax = self.ymax
        else:
            ymax = float(self.ymax_control.manual_value())

        self.axes.set_xbound(lower=xmin, upper=xmax)
        self.axes.set_ybound(lower=ymin, upper=ymax)
        
        # anecdote: axes.grid assumes b=True if any other flag is
        # given even if b is set to False.
        # so just passing the flag into the first statement won't
        # work.
        #
        if self.cb_grid.IsChecked():
            self.axes.grid(True, color='gray')
        else:
            self.axes.grid(False)

        # Using setp here is convenient, because get_xticklabels
        # returns a list over which one needs to explicitly 
        # iterate, and setp already handles this.
        #  
        pylab.setp(self.axes.get_yticklabels(), 
            visible=self.cb_ylab.IsChecked())
        
        self.plot_data.set_extent(extent=[self.xmin, self.xmax, self.ymin, self.ymax])
        
        self.canvas.draw()
    
    def on_pause_button(self, event):
        self.paused = not self.paused
    
    def on_update_pause_button(self, event):
        label = "Resume" if self.paused else "Pause"
        self.pause_button.SetLabel(label)
    
    def on_cb_grid(self, event):
        self.draw_plot()
    
    def on_cb_ylab(self, event):
        self.draw_plot()
    
    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)
    
    def on_redraw_timer(self, event):
        # if paused do not add data, but still redraw the plot
        # (to respond to scale modifications, grid change, etc.)
        #
        if not self.paused:
            self.data.append(self.datagen.next())
        
        self.draw_plot()
    
    def on_exit(self, event):
        self.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')


if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = GraphFrame()
    app.frame.Show()
    app.MainLoop()

