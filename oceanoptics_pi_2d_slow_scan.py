from __future__ import print_function, absolute_import, division
#from ScopeFoundryHW.mcl_stage import MCLStage2DSlowScan
from pi_stage_slowscan import PIStage2DSlowScan
import numpy as np
import time
import pyqtgraph as pg
#from ScopeFoundryHW.mcl_stage.mcl_stage_slowscan import MCLStage3DStackSlowScan

class OceanOptics_PI_2DSlowScan(PIStage2DSlowScan):
    
    name = 'OceanOptics_PI_2DSlowScan'  
    
    def pre_scan_setup(self):
        #hardware 
        #self.spec_hw = self.app.hardware['oceanoptics']
        #spec = self.spec_hw.spec # low level hardware
        
        #scan specific setup
        
        # create data arrays
        
        #self.num_hist_chans = self.picoharp_hw.calc_num_hist_chans()

        #time_trace_map_shape = self.scan_shape + (self.num_hist_chans,)
        #self.time_trace_map = np.zeros(time_trace_map_shape, dtype=float)
        
        
        #self.time_trace_map_h5 = self.h5_meas_group.create_dataset('time_trace_map', 
         #                                                          shape=time_trace_map_shape,
          #                                                         dtype=float, 
           #                                                        compression='gzip')
        
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
    
    # def collect_pixel(self, pixel_num, k, j, i):
        
    #     # collect data
    #     print('Picoharp_MCL_2DSlowScan', 'collect_pixel', pixel_num, k, j, i)
    #     t0 = time.time()

    #     #hist_data, elapsed_time = self.read_picoharp_histogram()
        
    #     ph = self.picoharp_hw.picoharp
    #     ph.start_histogram()
    #     while not ph.check_done_scanning():
    #         self.picoharp_hw.settings.count_rate0.read_from_hardware()
    #         self.picoharp_hw.settings.count_rate1.read_from_hardware()
    #         if self.picoharp_hw.settings['Tacq'] > 0.2:
    #             ph.read_histogram_data()
    #         time.sleep(0.005) #self.sleep_time)  
    #     ph.stop_histogram()
    #     #ta = time.time()
    #     ph.read_histogram_data()

    #     hist_data = ph.histogram_data
    #     elapsed_time = ph.read_elapsed_meas_time()
        
    #     # store in arrays
    #     self.time_trace_map[k,j,i, :] = hist_data[0:self.num_hist_chans]
    #     self.time_trace_map_h5[k,j,i, :] = hist_data[0:self.num_hist_chans]
        
    #     self.elapsed_time[k,j,i] = elapsed_time

    #     # display count-rate
    #     self.display_image_map[k,j,i] = hist_data[0:self.num_hist_chans].sum() * 1.0/elapsed_time
        
    #     import datetime
    #     print('pixel',  datetime.timedelta(seconds=(self.Npixels - pixel_num)*elapsed_time*1e-3), 'left')
        
    #     print( 'pixel done' )

    def _read_spectrometer(self):
        if hasattr(self, 'spec'):
            intg_time_ms = self.spec_hw.settings['intg_time']
            self.spec.integration_time_micros(intg_time_ms*1e3)
            
            scans_to_avg = self.spec_measure.settings['scans_to_avg']
            Int_array = np.zeros(shape=(2048,scans_to_avg))
            
            for i in range(scans_to_avg): #software average
                data = self.spec.spectrum(correct_dark_counts=self.spec_hw.settings['correct_dark_counts'])#ui.correct_dark_counts_checkBox.isChecked())
                Int_array[:,i] = data[1]
                self.y = np.mean(Int_array, axis=-1)

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

        
    def update_display(self):
        PIStage2DSlowScan.update_display(self)
        
        # setup lifetime window
        # if not hasattr(self, 'lifetime_graph_layout'):
        #     self.lifetime_graph_layout = pg.GraphicsLayoutWidget()
        #     self.lifetime_plot = self.lifetime_graph_layout.addPlot()
        #     self.lifetime_plotdata = self.lifetime_plot.plot()
        #     self.lifetime_plot.setLogMode(False, True)
        # self.lifetime_graph_layout.show()
        if hasattr(self, 'spec') and hasattr(self, 'pi_device') and hasattr(self, 'y'):
            self.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
            pg.QtGui.QApplication.processEvents()
        
        # kk, jj, ii = self.current_scan_index
        # spec = self.spec_hw.spec
        # self.lifetime_plotdata.setData(self.time_array,  1+ph.histogram_data[0:self.num_hist_chans])

    def run(self):
        """
        Runs when measurement is started. Runs in a separate thread from GUI.
        It should not update the graphical interface directly, and should only
        focus on data acquisition.
        """

        self.pi_device = self.pi_device_hw.pi_device
        self.spec = self.spec_hw.spec
        self.axes = self.pi_device_hw.axes

        x_start = self.settings['h0']
        y_start = self.settings['v0']
        
        x_scan_size = self.settings['h1'] - self.settings['h0']
        y_scan_size = self.settings['v1'] - self.settings['v0']
        
        x_step = self.settings['dh']
        y_step = self.settings['dv']
            
        y_range = self.settings['Nv'] #int(np.ceil(y_scan_size/y_step))
        x_range = self.settings['Nh'] #int(np.ceil(x_scan_size/x_step))
        
        # Define empty array for saving intensities
        data_array = np.zeros(shape=(x_range*y_range,2048))
        
        # Move to the starting position
        self.pi_device.MOV(axes=self.axes, values=[x_start,y_start])
        

        k = 0
        for i in range(y_range):
            for j in range(x_range):
                if self.interrupt_measurement_called:
                    break
                self._read_spectrometer()
                data_array[k,:] = self.y
                self.pi_device.MVR(axes=self.axes[0], values=[x_step])
                #self.ui.progressBar.setValue(np.floor(100*((k+1)/(x_range*y_range))))
                print(100*((k+1)/(x_range*y_range)))
                self.pi_device_hw.read_from_hardware()
                k+=1

            # TODO
            # if statement needs to be modified to keep the stage at the finish y-pos for line scans in x, and same for y
            if i == y_range-1: # this if statement is there to keep the stage at the finish position (in x) and not bring it back like we were doing during the scan 
                self.pi_device.MVR(axes=self.axes[1], values=[y_step])
                self.pi_device_hw.read_from_hardware()
            else:
                self.pi_device.MVR(axes=self.axes, values=[-x_scan_size, y_step])
                self.pi_device_hw.read_from_hardware()

            if self.interrupt_measurement_called:
                break

        save_dict = {"Wavelengths": self.spec.wavelengths(), "Intensities": data_array,
                 "Scan Parameters":{"X scan start (um)": x_start, "Y scan start (um)": y_start,
                                    "X scan size (um)": x_scan_size, "Y scan size (um)": y_scan_size,
                                    "X step size (um)": x_step, "Y step size (um)": y_step},
                                    "OceanOptics Parameters":{"Integration Time (ms)": self.spec_hw.settings['intg_time'],
                                                              "Scans Averages": self.spec_measure.settings['scans_to_avg'],
                                                              "Correct Dark Counts": self.spec_hw.settings['correct_dark_counts']}
                 }
    
        pickle.dump(save_dict, open(self.app.settings['save_dir']+"/"+self.app.settings['sample']+"_raw_PL_spectra_data.pkl", "wb"))


        
# class Picoharp_MCL_3DSlowScan(MCLStage3DStackSlowScan):
    
#     name = 'Picoharp_MCL_3DSlowScan'
    
#     def pre_scan_setup(self):
#         #hardware 
#         self.picoharp_hw = self.app.hardware['picoharp']
#         ph = self.picoharp_hw.picoharp # low level hardware
        
#         #scan specific setup
        
#         # create data arrays
        
#         cr0 = self.picoharp_hw.settings.count_rate0.read_from_hardware()
#         rep_period_s = 1.0/cr0
#         time_bin_resolution = self.picoharp_hw.settings['Resolution']*1e-12
#         self.num_hist_chans = int(np.ceil(rep_period_s/time_bin_resolution))

#         time_trace_map_shape = self.scan_shape + (self.num_hist_chans,)
#         self.time_trace_map = np.zeros(time_trace_map_shape, dtype=float)
        
        
#         self.time_trace_map_h5 = self.create_h5_framed_dataset(name='time_trace_map', single_frame_map=self.time_trace_map)
        
#         self.time_array = self.h5_meas_group['time_array'] = ph.time_array[0:self.num_hist_chans]*1e-3
        
#         self.elapsed_time = np.zeros(self.scan_shape, dtype=float)
#         self.elapsed_time_h5 = self.create_h5_framed_dataset('elasped_time', self.elapsed_time)
        
#         #self.app.settings_auto_save()
        
#         # pyqt graph
#         #self.initial_scan_setup_plotting = True

        
#     def post_scan_cleanup(self):
#         # close shutter 
#         #self.gui.shutter_servo_hc.shutter_open.update_value(False)
#         pass
    
#     def collect_pixel(self, pixel_num, frame_i, k, j, i):
        
#         # collect data
#         print(pixel_num, frame_i, k, j, i)
#         t0 = time.time()
        
#         hist_data, elapsed_time = Picoharp_MCL_2DSlowScan.read_picoharp_histogram(self)

#         # store in arrays
#         self.time_trace_map[k,j,i, :] = hist_data[0:self.num_hist_chans]
#         self.time_trace_map_h5[frame_i,k,j,i, :] = hist_data[0:self.num_hist_chans]
        
#         self.elapsed_time[k,j,i] = elapsed_time
#         self.elapsed_time_h5[frame_i, k,j,i] = elapsed_time

#         # display count-rate
#         self.display_image_map[k,j,i] = hist_data[0:self.num_hist_chans].sum() * 1.0/elapsed_time
        
#         print( 'pixel done' )
        
#     def on_new_frame(self, frame_i):
#         MCLStage3DStackSlowScan.on_new_frame(self, frame_i)
#         self.extend_h5_framed_dataset(self.time_trace_map_h5, frame_i)
#         self.extend_h5_framed_dataset(self.elapsed_time_h5, frame_i)
        
#     def on_end_frame(self, frame_i):
#         MCLStage3DStackSlowScan.on_end_frame(self, frame_i)
#         self.h5_file.flush()
        
#     def update_display(self):
#         PIStage2DSlowScan.update_display(self)
#     