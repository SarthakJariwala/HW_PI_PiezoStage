from __future__ import print_function, absolute_import, division
from pi_stage_slowscan import PIStage2DSlowScan
#from ScopeFoundryHW.mcl_stage import MCLStage2DSlowScan

import numpy as np
import time
import pyqtgraph as pg
#from ScopeFoundryHW.mcl_stage.mcl_stage_slowscan import MCLStage3DStackSlowScan

class OceanOptics_PI_2DSlowScan(PIStage2DSlowScan):
    
    name = 'oceanoptics_pi_2DSlowScan'

    def pre_scan_setup(self):
        #hardware 
        self.spec_hw = self.app.hardware['oceanoptics']
        self.pi_device_hw = self.app.hardware['piezostage']
        #ph = self.picoharp_hw.picoharp # low level hardware
        
        #scan specific setup
        
        # create data arrays
        
        #self.num_hist_chans = self.picoharp_hw.calc_num_hist_chans()

        #time_trace_map_shape = self.scan_shape + (self.num_hist_chans,)
        # self.time_trace_map = np.zeros(self.scan_shape, dtype=float)
        
        
        # self.time_trace_map_h5 = self.h5_meas_group.create_dataset('time_trace_map', 
        #                                                            shape=time_trace_map_shape,
        #                                                            dtype=float, 
        #                                                            compression='gzip')
        
        #self.time_array = self.h5_meas_group['time_array'] = ph.time_array[0:self.num_hist_chans]*1e-3
        #self.elapsed_time = self.h5_meas_group['elapsed_time'] = np.zeros(self.scan_shape, dtype=float)
        
        #self.app.settings_auto_save()
        

        # pyqt graph
        self.initial_scan_setup_plotting = True


        # set up experiment
        # experimental parameters already connected via LoggedQuantities
        
        # open shutter 
        # self.gui.shutter_servo_hc.shutter_open.update_value(True)
        # time.sleep(0.5)
        

        
    def post_scan_cleanup(self):
        # close shutter 
        #self.gui.shutter_servo_hc.shutter_open.update_value(False)
        pass
    
    def collect_pixel(self, pixel_num, k, j, i):
        
        # collect data
        print('OceanOptics_PI_2DSlowScan', 'collect_pixel', pixel_num, k, j, i)
        t0 = time.time()

        #hist_data, elapsed_time = self.read_picoharp_histogram()
        
        #ph = self.picoharp_hw.picoharp
        spec = self.spec_hw.spec
        pi_device = self.pi_devce.spec

        # #ph.start_histogram()
        # while not ph.check_done_scanning():
        #     self.picoharp_hw.settings.count_rate0.read_from_hardware()
        #     self.picoharp_hw.settings.count_rate1.read_from_hardware()
        #     if self.picoharp_hw.settings['Tacq'] > 0.2:
        #         ph.read_histogram_data()
        #     time.sleep(0.005) #self.sleep_time)  
        # ph.stop_histogram()
        # #ta = time.time()
        # ph.read_histogram_data()

        # hist_data = ph.histogram_data
        self._read_spectrometer()
        self.save_array[:,1] = self.y
        self.save_array[:,0] = self.spec.wavelengths()

        t1 = time.time()
        elapsed_time = t1-t0
        
        # store in arrays
        #self.time_trace_map[k,j,i, :] = hist_data[0:self.num_hist_chans]
        #self.time_trace_map_h5[k,j,i, :] = hist_data[0:self.num_hist_chans]
        
        self.elapsed_time[k,j,i] = elapsed_time

        # display count-rate
        self.display_image_map[k,j,i] = self.save_array[:,1].sum() * 1.0/elapsed_time
        
        import datetime
        print('pixel',  datetime.timedelta(seconds=(self.Npixels - pixel_num)*elapsed_time*1e-3), 'left')
        
        print( 'pixel done' )
    
    # def read_picoharp_histogram(self):
    #     print("asdf")

    #     ph = self.picoharp_hw.picoharp

    #     ph.start_histogram()

    #     while not ph.check_done_scanning():
    #         if self.picoharp_hw.settings['Tacq'] > 200:
    #             ph.read_histogram_data()
    #         time.sleep(0.005) #self.sleep_time)  
    #     ph.stop_histogram()
    #     #ta = time.time()
    #     ph.read_histogram_data()

    #     return ph.histogram_data, ph.read_elapsed_meas_time()


    def _read_spectrometer(self):
        if hasattr(self, 'spec'):
            intg_time_ms = 1 #for now, hardcode for intg time of 1ms #self.spec_hw.settings['intg_time']
            self.spec.integration_time_micros(intg_time_ms*1e3)
            
            scans_to_avg = 1 #for now, hardcode as 1 scan to average #self.settings['scans_to_avg']
            Int_array = np.zeros(shape=(2048,scans_to_avg))
            
            for i in range(scans_to_avg): #software average
                data = self.spec.spectrum(correct_dark_counts=True)#self.spec_hw.settings['correct_dark_counts'])#ui.correct_dark_counts_checkBox.isChecked())
                Int_array[:,i] = data[1]
                self.y = np.mean(Int_array, axis=-1)

        
    def update_display(self):
        PIStage2DSlowScan.update_display(self)
        
        # setup lifetime window
        if not hasattr(self, 'lifetime_graph_layout'):
            self.lifetime_graph_layout = pg.GraphicsLayoutWidget()
            self.lifetime_plot = self.lifetime_graph_layout.addPlot()
            self.lifetime_plotdata = self.lifetime_plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
            self.lifetime_plot.setLogMode(False, True)
        self.lifetime_graph_layout.show()
        pg.QtGui.QApplication.processEvents()

        
        kk, jj, ii = self.current_scan_index
        #ph = self.picoharp_hw.picoharp
        #self.lifetime_plotdata.setData(self.time_array,  1+ph.histogram_data[0:self.num_hist_chans])