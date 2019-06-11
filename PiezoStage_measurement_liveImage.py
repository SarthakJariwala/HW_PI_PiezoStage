from ScopeFoundry import Measurement
from ScopeFoundry.helper_funcs import sibling_path, load_qt_ui_file
from ScopeFoundry import h5_io
from OceanOptics_hardware import OceanOpticsHW
from OceanOptics_measurement import OceanOpticsMeasure
import pyqtgraph as pg
import numpy as np
import time
import pickle
fromt qtpy import QtCore

class PiezoStageMeasureLive(Measurement):
	
	# this is the name of the measurement that ScopeFoundry uses 
	# when displaying your measurement and saving data related to it    
	name = "oceanoptics_scan"
	
	def setup(self):
		"""
		Runs once during App initialization.
		This is the place to load a user interface file,
		define settings, and set up data structures. 
		"""
		
		# Define ui file to be used as a graphical interface
		# This file can be edited graphically with Qt Creator
		# sibling_path function allows python to find a file in the same folder
		# as this python module
		self.ui_filename = sibling_path(__file__, "spec_pz.ui")
		
		#Load ui file and convert it to a live QWidget of the user interface
		self.ui = load_qt_ui_file(self.ui_filename)

		# Measurement Specific Settings
		# This setting allows the option to save data to an h5 data file during a run
		# All settings are automatically added to the Microscope user interface

		self.settings.New('x_start', dtype=float, unit='um', vmin=0)
		self.settings.New('y_start', dtype=float, unit='um', vmin=0)
	
		self.settings.New('x_size', dtype=float, initial=1, unit='um', vmin=1)
		self.settings.New('y_size', dtype=float, initial=1, unit='um', vmin=1)

		self.settings.New('x_step', dtype=float, initial=1, unit='um', vmin=1)
		self.settings.New('y_step', dtype=float, initial=1, unit='um', vmin=1)


		
		# Create empty numpy array to serve as a buffer for the acquired data
		self.buffer = np.zeros(120, dtype=float)
		
		# Define how often to update display during a run
		self.display_update_period = 0.1 
		
		# Convenient reference to the hardware used in the measurement
		#self.func_gen = self.app.hardware['virtual_function_gen']
		self.spec_hw = self.app.hardware['oceanoptics']
		self.pi_device_hw = self.app.hardware['piezostage']
		
		self.spec_measure = self.app.measurements['oceanoptics_measure']

		self.initial_scan_setup_plotting = False

	def setup_figure(self):
		"""
		Runs once during App initialization, after setup()
		This is the place to make all graphical interface initializations,
		build plots, etc.
		"""
		
		# connect ui widgets to measurement/hardware settings or functions
		self.ui.start_scan_pushButton.clicked.connect(self.start)
		self.ui.interrupt_scan_pushButton.clicked.connect(self.interrupt)

		self.pi_device_hw.settings.x_position.connect_to_widget(self.ui.x_pos_doubleSpinBox)
		self.pi_device_hw.settings.y_position.connect_to_widget(self.ui.y_pos_doubleSpinBox)

		#TODO: connect to ROI widget so that these changes are reflected.
		#currently, there is only a one way connection w/ ROI > spinboxes
		self.settings.x_start.connect_to_widget(self.ui.x_start_doubleSpinBox)
		self.settings.y_start.connect_to_widget(self.ui.y_start_doubleSpinBox)
		self.settings.x_size.connect_to_widget(self.ui.x_size_doubleSpinBox)
		self.settings.y_size.connect_to_widget(self.ui.y_size_doubleSpinBox)
		self.settings.x_step.connect_to_widget(self.ui.x_step_doubleSpinBox)
		self.settings.y_step.connect_to_widget(self.ui.y_step_doubleSpinBox)
		
		self.spec_hw.settings.intg_time.connect_to_widget(self.ui.intg_time_doubleSpinBox)
		self.spec_hw.settings.correct_dark_counts.connect_to_widget(self.ui.correct_dark_counts_checkBox)
		self.spec_measure.settings.scans_to_avg.connect_to_widget(self.ui.scans_to_avg_spinBox)


		# Set up pyqtgraph graph_layout in the UI
		self.graph_layout=pg.GraphicsLayoutWidget()
		self.ui.plot_groupBox.layout().addWidget(self.graph_layout)

		# # Create PlotItem object (a set of axes)  
		self.plot = self.graph_layout.addPlot(title="Spectrometer Live Reading")
		self.plot.setLabel('left', 'Intensity', unit='a.u.')
		self.plot.setLabel('bottom', 'Wavelength', unit='nm')

		#stage ui base
		self.stage_layout=pg.GraphicsLayoutWidget()
		self.ui.stage_groupBox.layout().addWidget(self.stage_layout)
		self.stage_plot = self.stage_layout.addPlot(title="Stage view")
		self.stage_plot.setXRange(0, 100)
		self.stage_plot.setYRange(0, 100)
		self.stage_plot.setLimits(xMin=0, xMax=100, yMin=0, yMax=100) 

		#region of interest
		self.scan_roi = pg.ROI([0,0],[25, 25], movable=True)
		self.scan_roi.addScaleHandle([1, 1], [0, 0])
		self.scan_roi.addScaleHandle([0, 0], [1, 1])		
		self.scan_roi.sigRegionChangeFinished.connect(self.mouse_update_scan_roi)
		self.stage_plot.addItem(self.scan_roi)

		#image display container
		self.img_items = []				
		self.img_item = pg.ImageItem()
		self.img_items.append(self.img_item)
		self.stage_plot.addItem(self.img_item)
		self.stage_plot.showGrid(x=True, y=True)
		self.stage_plot.setAspectLocked(lock=True, ratio=1)

		#histogram
		self.hist_lut = pg.HistogramLUTItem()
		self.stage_layout.addItem(self.hist_lut)
		
		#arrow
		self.current_stage_pos_arrow = pg.ArrowItem()
		self.current_stage_pos_arrow.setZValue(100)
		self.stage_plot.addItem(self.current_stage_pos_arrow)
		if hasattr(self, 'stage'):
			self.pi_device_hw.settings.x_position.updated_value.connect(self.update_arrow_pos, QtCore.Qt.UniqueConnection)
			self.pi_device_hw.settings.y_position.updated_value.connect(self.update_arrow_pos, QtCore.Qt.UniqueConnection)

		# for lqname in 'x_start x_size y_start y_size'.split():
		# 	self.settings.as_dict()[lqname].updated_value.connect(self.update_scan_roi)
		# # Create PlotDataItem object ( a scatter plot on the axes )
		self.optimize_plot_line = self.plot.plot([0])

	def mouse_update_scan_roi(self):
		x0,y0 =  self.scan_roi.pos()
		w, h =  self.scan_roi.size()
		print(w)
		print(h)
		#self.h_center.update_value(x0 + w/2)
		#self.v_center.update_value(y0 + h/2)
		#self.h_span.update_value(w-self.dh.val)
		#self.v_span.update_value(h-self.dv.val)
		self.settings['x_start'] = x0
		self.settings['y_start'] = y0
		self.settings['x_size'] = w
		self.settings['y_size'] = h
		#self.compute_scan_params()
		#self.update_scan_roi()

	def update_arrow_pos(self):
		x = self.pi_device_hw.settings['x_position']
		y = self.pi_device_hw.settings['y_position']
		self.current_stage_pos_arrow.setPos(x,y)

	#DONE?
	def update_display(self):
		"""
		Displays (plots) the numpy array self.buffer. 
		This function runs repeatedly and automatically during the measurement run.
		its update frequency is defined by self.display_update_period
		"""
		if self.initial_scan_setup_plotting:
			self.img_item = pg.ImageItem()
			self.img_items.append(self.img_item)
			self.stage_plot.addItem(self.img_item)
			self.hist_lut.setImageItem(self.img_item)
	
			self.img_item.setImage(self.display_image_map[0,:])
			#self.imshow_extent
			#self.log.debug('update_display set bounds {} {} {} {}'.format(x0, x1, y0, y1))
			self.img_item_rect = QtCore.QRectF(self.settings.x_start, self.settings.y_start,
				self.settings.x_size, self.settings.y_size)
			self.img_item.setRect(self.img_item_rect)
			#self.log.debug('update_display set bounds {}'.format(self.img_item_rect))
			
			self.initial_scan_setup_plotting = False
		else:
			#if self.settings.scan_type.val in ['raster']
			kk, jj = self.scan_index
			self.disp_img = self.display_image_map[kk,:].T
			self.img_item.setImage(self.disp_img, autoRange=False, autoLevels=True)
			self.img_item.setRect(self.img_item_rect) # Important to set rectangle after setImage for non-square pixels
			#self.update_LUT()

		if hasattr(self, 'spec') and hasattr(self, 'pi_device') and hasattr(self, 'y'):
			self.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
			pg.QtGui.QApplication.processEvents()

	#refer to slow scan code
	def run(self):
		"""
		Runs when measurement is started. Runs in a separate thread from GUI.
		It should not update the graphical interface directly, and should only
		focus on data acquisition.
		"""
		
		
		self.pi_device = self.pi_device_hw.pi_device
		self.spec = self.spec_hw.spec
		self.axes = self.pi_device_hw.axes

		x_start = self.settings['x_start']
		y_start = self.settings['y_start']
		
		x_scan_size = self.settings['x_size']
		y_scan_size = self.settings['y_size']
		
		x_step = self.settings['x_step']
		y_step = self.settings['y_step']
			
		y_range = int(np.ceil(y_scan_size/y_step))
		x_range = int(np.ceil(x_scan_size/x_step))

		self.display_image_map = np.zeros((y_range, x_range), dtype=float)
		
		# Define empty array for saving intensities
		data_array = np.zeros(shape=(x_range*y_range,2048))

		self.initial_scan_setup_plotting = True
		
		# Move to the starting position
		self.pi_device.MOV(axes=self.axes, values=[x_start,y_start])
		

		self.scan_index = (0, 0)
		for i in range(y_range):
			for j in range(x_range):
				if self.interrupt_measurement_called:
					break
				self._read_spectrometer()
				data_array[self.scan_index,:] = self.y
				self.collect_pixel(i*j, j, i)
				self.pi_device.MVR(axes=self.axes[0], values=[x_step])
				#self.ui.progressBar.setValue(np.floor(100*((k+1)/(x_range*y_range))))
				#print(100*((self.scan_index+1)/(x_range*y_range)))
				self.pi_device_hw.read_from_hardware()
				self.scan_index = (j, i)

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

	def collect_pixel(self, k, j): #got rid of pixel_num
		
		# collect data
		#print('OceanOptics_PI_2DSlowScan', 'collect_pixel', pixel_num, k, j, i)
		t0 = time.time()

		#hist_data, elapsed_time = self.read_picoharp_histogram()
		

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
		current_spec_array = self.get_single_spec()
		current_intensities = current_spec_array[:,1] 

		t1 = time.time()
		#elapsed_time = t1-t0
		
		# store in arrays
		#self.time_trace_map[k,j,i, :] = hist_data[0:self.num_hist_chans]
		#self.time_trace_map_h5[k,j,i, :] = hist_data[0:self.num_hist_chans]
		
		#self.elapsed_time[k,j] = elapsed_time

		# display count-rate
		self.display_image_map[k,j] = current_intensities.sum() #* 1.0/elapsed_time
		
		#import datetime
		#print('pixel',  datetime.timedelta(seconds=(self.Npixels - pixel_num)*elapsed_time*1e-3), 'left')
		
		print( 'pixel done' )

	def get_single_spec(self):
		save_array = np.zeros(shape=(2048,2))
		save_array[:,1] = self.y
		save_array[:,0] = self.spec.wavelengths()

		return save_array
		#np.savetxt(self.app.settings['save_dir']+"/"+self.app.settings['sample']+".txt", save_array, fmt = '%.5f', 
		#		   header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')