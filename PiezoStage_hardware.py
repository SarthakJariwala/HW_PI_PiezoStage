import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import pickle
import sys
import seabreeze.spectrometers as sb
from pipython import GCSDevice
import time
from ScopeFoundry import HardwareComponent

class PiezoStageHW(HardwareComponent):
    
    ## Define name of this hardware plug-in

    
    def setup(self):
        # Define your hardware settings here.
        # These settings will be displayed in the GUI and auto-saved with data files
        self.name = 'piezostage'
        self.settings.New('x_pos', dtype=float, unit='um')
        self.settings.New('y_pos', dtype=float, unit='um')

    def connect(self):
        # Open connection to the device:
        self.pi_device = GCSDevice("E-710")	# Creates a Controller instant
        self.pi_device.ConnectNIgpib(board=0,device=4) # Connect to GPIB board
        #self.ui.status_textBrowser.append('Connected: {}'.format(self.pi_device.qIDN().strip()))
#        print('connected: {}'.format(self.pi_device.qIDN().strip()))
        
        self.axes = self.pi_device.axes[0:2] # selecting x and y axis of the stage
        
        self.pi_device.INI()
        self.pi_device.REF(axes=self.axes)
        
        self.pi_device.SVO(axes=self.axes, values=[True,True])	# Turn on servo control for both axes
        #self.ui.status_textBrowser.append("Current Stage Position:\n{}".format(self.pi_device.qPOS(axes=self.axes)))
#        print(self.pi_device.qPOS(axes=self.axes))

        self.settings['x_pos'] = self.pi_device.qPOS(axes=self.axes)['1']
        self.settings['y_pos'] = self.pi_device.qPOS(axes=self.axes)['2']
        #Connect settings to hardware:
        #self.settings.x_pos.connect_to_hardware(self.x_axis)
        #self.settings.y_pos.connect_to_hardware(self.y_axis)

    
        #Take an initial sample of the data.
        self.read_from_hardware()
        
    def disconnect(self):
        #Disconnect the device and remove connections from settings
        self.settings.disconnect_all_from_hardware()
        if hasattr(self, 'pi_device'):
            self.pi_device.close()
            del self.pi_device
            self.pi_device = None