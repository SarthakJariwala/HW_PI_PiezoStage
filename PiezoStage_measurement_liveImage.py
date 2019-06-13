from ScopeFoundry import Measurement
from ScopeFoundry.helper_funcs import sibling_path, load_qt_ui_file
import pyqtgraph as pg
import numpy as np
import time
import pickle
import os.path

class PiezoStageMeasureLive(Measurement):

	
	def setup(self):
		"""
		Runs once during App initialization.
		This is the place to load a user interface file,
		define settings, and set up data structures. 
		"""
		self.name = "oceanoptics_scan_liveupdate"
		
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
	
		self.settings.New('x_size', dtype=float, initial=1, unit='um', vmin=0)
		self.settings.New('y_size', dtype=float, initial=1, unit='um', vmin=0)

		self.settings.New('x_step', dtype=float, initial=1, unit='um', vmin=.001)
		self.settings.New('y_step', dtype=float, initial=1, unit='um', vmin=.001)



		self.settings.New('x_clicked', dtype=float, initial=0, unit='um', vmin=0, vmax=100, ro=True)
		self.settings.New('y_clicked', dtype=float, initial=0, unit='um', vmin=0, vmax=100, ro=True)
		
		
		# Define how often to update display during a run
		self.display_update_period = 0.1 
		
		# Convenient reference to the hardware used in the measurement
		self.spec_hw = self.app.hardware['oceanoptics']
		self.pi_device_hw = self.app.hardware['piezostage']
		
		self.spec_measure = self.app.measurements['oceanoptics_measure']

		self.scan_complete = False

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
		self.settings.x_start.connect_to_widget(self.ui.x_start_doubleSpinBox)
		self.settings.y_start.connect_to_widget(self.ui.y_start_doubleSpinBox)	
		self.settings.x_size.connect_to_widget(self.ui.x_size_doubleSpinBox)
		self.settings.y_size.connect_to_widget(self.ui.y_size_doubleSpinBox)
		self.settings.x_step.connect_to_widget(self.ui.x_step_doubleSpinBox)
		self.settings.y_step.connect_to_widget(self.ui.y_step_doubleSpinBox)
		self.settings.x_clicked.connect_to_widget(self.ui.x_clicked_doubleSpinBox)
		self.settings.y_clicked.connect_to_widget(self.ui.y_clicked_doubleSpinBox)
		self.ui.move_to_selected_pushButton.connect(self.move_to_selected)

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
		
		# # Create PlotDataItem object ( a scatter plot on the axes )
		self.optimize_plot_line = self.plot.plot([0])

		#stage ui base
		self.stage_layout=pg.GraphicsLayoutWidget()
		self.ui.stage_groupBox.layout().addWidget(self.stage_layout)
		self.stage_plot = self.stage_layout.addPlot(title="Stage view")
		self.stage_plot.setXRange(0, 100)
		self.stage_plot.setYRange(0, 100)
		self.stage_plot.setLimits(xMin=0, xMax=100, yMin=0, yMax=100) 

		#region of interest - allows user to select scan area
		self.scan_roi = pg.ROI([0,0],[25, 25], movable=True)
		self.scan_roi.addScaleHandle([1, 1], [0, 0])
		self.scan_roi.addScaleHandle([0, 0], [1, 1])		
		self.scan_roi.sigRegionChangeFinished.connect(self.mouse_update_scan_roi)
		self.stage_plot.addItem(self.scan_roi)

		self.ui.x_start_doubleSpinBox.valueChanged.connect(self.update_roi_start)
		self.ui.y_start_doubleSpinBox.valueChanged.connect(self.update_roi_start)
		self.ui.x_size_doubleSpinBox.valueChanged.connect(self.update_roi_size)
		self.ui.y_size_doubleSpinBox.valueChanged.connect(self.update_roi_size)

		#image display container
		# self.imv = pg.ImageView()
		# self.imv.getView().setAspectLocked(lock=False, ratio=1)
		# self.imv.getView().setMouseEnabled(x=True, y=True)

		#image on stage plot
		self.img_item = pg.ImageItem()
		self.stage_plot.addItem(img_item)

		#arrow showing stage location
		self.current_stage_pos_arrow = pg.ArrowItem()
		self.current_stage_pos_arrow.setZValue(100)
		self.stage_plot.addItem(self.current_stage_pos_arrow)
		self.pi_device_hw.settings.x_position.updated_value.connect(self.update_arrow_pos, QtCore.Qt.UniqueConnection)
		self.pi_device_hw.settings.y_position.updated_value.connect(self.update_arrow_pos, QtCore.Qt.UniqueConnection)

		#Define crosshairs that will show up after scan, event handling.
		self.vLine = pg.InfiniteLine(angle=90, movable=False, pen='w')
		self.hLine = pg.InfiniteLine(angle=0, movable=False, pen='w')
		proxy = pg.SignalProxy(self.stage_plot.scene().sigMouseMoved, rateLimit=60, slot=self.chMove) #connect plot item to mouse moved, which handles crosshair movement
		self.stage_plot.scene().sigMouseClicked.connect(self.chClick)

	def chClick(self, event):
		'''
		Handle crosshair clicking, which toggles movement on and off.
		'''
		items = self.stage_plot.scene().items(event.scenePos()) #get items at clicked position
		if (self.vLine in items or hLine in items): #if crosshair is clicked, toggle movement
			self.move_ch = not self.move_ch
		if not self.move_ch: #if crosshair has been dropped, update movement
			ch_pos = self.stage_plot.vb.mapSceneToView(pos)
			self.settings['x_clicked'] = ch_pos[0]
			self.settings['y_clicked'] = ch_pos[1]

	def chMove(self, event):
		'''
		Handle crosshair movement. Crosshair will only move with mouse if toggled on.
		'''
		pos = event[0]
		if self.stage_plot.sceneBoundingRect().contains(pos): #check if mouse within bounds
			mousePoint = self.stage_plot.vb.mapSceneToView(pos) #convert device coordinates to scene coordinates
			if self.move_ch: #move crosshair only if toggled on by clicking
				self.vLine.setPos(mousePoint.x())
				self.hLine.setPos(mousePoint.y())

	def move_to_selected(self):
		'''
		Move stage to position selected by crosshairs.
		'''
		if self.scan_complete:
			x = self.settings['x_clicked']
			y = self.settings['y_clicked']
			self.pi_device.MOV(axes=self.axes, values=[x, y])
			self.pi_device_hw.read_from_hardware()

	
	def mouse_update_scan_roi(self):
		'''
		Update settings to reflect region of interest.
		'''
		x0,y0 =  self.scan_roi.pos()
		w, h =  self.scan_roi.size()
		self.settings['x_start'] = x0
		self.settings['y_start'] = y0
		self.settings['x_size'] = w
		self.settings['y_size'] = h

	def update_roi_start(self):
		'''
		Update region of interest start position according to spinboxes
		'''
		self.scan_roi.setPos((self.settings['x_start'], self.settings['y_start']))

	def update_roi_size(self):
		'''
		Update region of interest size according to spinboxes
		'''
		self.scan_roi.setSize((self.settings['x_size'], self.settings['y_size']))

	def update_arrow_pos(self):
		'''
		Update arrow position on image to stage position
		'''
		x = self.pi_device_hw.settings['x_position']
		y = self.pi_device_hw.settings['y_position']
		self.current_stage_pos_arrow.setPos(x,y)


	def update_display(self):
		"""
		Displays (plots) the numpy array self.buffer. 
		This function runs repeatedly and automatically during the measurement run.
		its update frequency is defined by self.display_update_period
		"""
		self.img_item.setPos(self.settings['x_start'], self.settings['y_start'])
		if hasattr(self, 'spec') and hasattr(self, 'pi_device') and hasattr(self, 'y'): #first, check if setup has happened
			#plot wavelengths vs intensity
			self.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True) #plot wavelength vs intensity
			pg.QtGui.QApplication.processEvents()

			disp_img = self.display_image_map #transpose to use for setImage, which takes 3d array (x, y, intensity)
			self.img_item.setImage(image=disp_img, autoLevels=True, autoRange=False)
			#self.imv.setImage(img=disp_img, autoRange=False, autoLevels=True)
			#self.imv.show()

		if self.scan_complete:
			stage_plot.addItem(self.hLine)
			stage_plot.addItem(self.vLine)

			middle_x = self.settings['x_start'] + (self.settings['x_size']/2)
			middle_y = self.settings['y_start'] + (self.settings['y_size']/2)
			self.hLine.setPos(middle_y)
			self.vLine.setPos(middle_x)



	def run(self):
		"""
		Runs when measurement is started. Runs in a separate thread from GUI.
		It should not update the graphical interface directly, and should only
		focus on data acquisition.

		Runs until scan is completed or interrupted.
		"""
		self.check_filename()
		
		self.pi_device = self.pi_device_hw.pi_device
		self.spec = self.spec_hw.spec
		self.axes = self.pi_device_hw.axes

		x_start = self.settings['x_start']
		y_start = self.settings['y_start']
		
		x_scan_size = self.settings['x_size']
		y_scan_size = self.settings['y_size']
		
		x_step = self.settings['x_step']
		y_step = self.settings['y_step']
		
		if y_scan_size == 0:
			y_scan_size = 1#self.settings['y_size'] = 1
			y_step = 1#self.settings['y_step'] = 1
        
		if x_scan_size == 0:
			x_scan_size = 1#self.settings['x_size'] = 1
			x_step = 1#self.settings['x_step'] = 1
        
		if y_step == 0:
			y_step = 1#self.settings['y_step'] = 1
            
		if x_step == 0:
			x_step = 1#self.settings['x_step'] = 1
			
		#number of scans in x and y
		self.y_range = int(np.ceil(y_scan_size/y_step))
		self.x_range = int(np.ceil(x_scan_size/x_step))
		
		# Define empty array for saving intensities
		data_array = np.zeros(shape=(self.x_range*self.y_range,2048))

		# Define empty array for image map
		#self.display_image_map = np.zeros((y_range, x_range), dtype=float)

		#Store spectrum for each pixel
		self.display_image_map = np.zeros((2048, self.x_range, self.y_range), dtype=float)
		
		# Move to the starting position
		self.pi_device.MOV(axes=self.axes, values=[x_start,y_start])
		self.pi_device_hw.read_from_hardware()
		

		k = 0 #keep track of scan/'pixel' number
		for i in range(self.y_range):
			for j in range(self.x_range):
				if self.interrupt_measurement_called:
					break
				self._read_spectrometer()
				data_array[k,:] = self.y
				#intensities_sum = data_array[k,:].sum()
				self.display_image_map[:, j, i] = self.y#intensities_sum
				self.pi_device.MVR(axes=self.axes[0], values=[x_step])
				#self.ui.progressBar.setValue(np.floor(100*((k+1)/(x_range*y_range))))
				print(100*((k+1)/(self.x_range*self.y_range)))
				self.pi_device_hw.read_from_hardware()
				k+=1
			# TODO
			# if statement needs to be modified to keep the stage at the finish y-pos for line scans in x, and same for y
			if i == self.y_range-1: # this if statement is there to keep the stage at the finish position (in x) and not bring it back like we were doing during the scan 
				self.pi_device.MVR(axes=self.axes[1], values=[y_step])
				self.pi_device_hw.read_from_hardware()
			else:				
				self.pi_device.MVR(axes=self.axes[1], values=[y_step])
				self.pi_device.MOV(axes=self.axes[0], values=[x_start])
				self.pi_device_hw.read_from_hardware()
				#self.pi_device_hw.read_from_hardware()

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
		self.scan_complete = True;

	def _read_spectrometer(self):
		'''
		Read spectrometer according to settings and update self.y (intensities array)
		'''
		if hasattr(self, 'spec'):
			intg_time_ms = self.spec_hw.settings['intg_time']
			self.spec.integration_time_micros(intg_time_ms*1e3) #seabreeze error checking
			
			scans_to_avg = self.spec_measure.settings['scans_to_avg']
			Int_array = np.zeros(shape=(2048,scans_to_avg))
			
			for i in range(scans_to_avg): #software average
				data = self.spec.spectrum(correct_dark_counts=self.spec_hw.settings['correct_dark_counts'])#acquire wavelengths and intensities from spec
				Int_array[:,i] = data[1]
				self.y = np.mean(Int_array, axis=-1)

	def check_filename(self):
		'''
		If no sample name given or duplicate sample name given, fix the problem by appending a unique number.
		'''
		samplename = self.app.settings['sample']
		filename = samplename + "_raw_PL_spectra_data.pkl"
		directory = self.app.settings['save_dir']
		if samplename == "":
			self.app.settings['sample'] = int(time.time())
		if (os.path.exists(directory+"/"+filename)):
			self.app.settings['sample'] = samplename + str(int(time.time()))
