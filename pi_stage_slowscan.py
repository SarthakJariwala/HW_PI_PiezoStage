from __future__ import division, print_function
#import numpy as np
from ScopeFoundry.scanning import BaseRaster2DSlowScan, BaseRaster2DFrameSlowScan
#from ScopeFoundry import Measurement, LQRange
import time

class PIStage2DSlowScan(BaseRaster2DSlowScan):
    
    name = "PIStage2DSlowScan"
    def __init__(self, app):
        BaseRaster2DSlowScan.__init__(self, app, h_limits=(1,100), v_limits=(1,100),
                                      h_spinbox_step = 0.1, v_spinbox_step=0.1,
                                      h_unit="um", v_unit="um")        
    
    def setup(self):
        BaseRaster2DSlowScan.setup(self)
        
        self.settings.New("h_axis", initial="X", dtype=str, choices=("X", "Y"))
        self.settings.New("v_axis", initial="Y", dtype=str, choices=("X", "Y"))
        
        self.ax_map = dict(X=0, Y=1)

        #Hardware
        self.stage = self.app.hardware['piezostage']
        
        # self.settings.h_axis.add_listener(self.on_new_stage_limits)
        # self.settings.v_axis.add_listener(self.on_new_stage_limits)
        #self.stage.settings.x_max.add_listener(self.on_new_stage_limits)
        
    # def on_new_stage_limits(self):
    #     h_axis = self.settings['h_axis'].lower()
    #     v_axis = self.settings['v_axis'].lower()
    #     h_max = self.stage.settings[h_axis + '_max']
    #     v_max = self.stage.settings[v_axis + '_max']
        
    #     self.set_h_limits(0.1, h_max-0.1)
    #     self.set_v_limits(0.1, v_max-0.1)
        
    def setup_figure(self):
        BaseRaster2DSlowScan.setup_figure(self)
        self.set_details_widget(widget=self.settings.New_UI(include=['h_axis', 'v_axis']))
        
    def move_position_start(self, h,v):
        #self.stage.y_position.update_value(x)
        #self.stage.y_position.update_value(y)
        
        S = self.settings
        
        coords = [None, None, None]
        coords[self.ax_map[S['h_axis']]] = h
        coords[self.ax_map[S['v_axis']]] = v
        
        #self.stage.move_pos_slow(x,y,None)
        self.stage.abs_mov(h, v)
    
    # def move_position_slow(self, h,v, dh,dv):
    #     self.move_position_start(h, v)

    # def move_position_fast(self,  h,v, dh,dv):
    #     #self.stage.x_position.update_value(x)
    #     S = self.settings        
    #     coords = [None, None, None]
    #     coords[self.ax_map[S['h_axis']]] = h
    #     coords[self.ax_map[S['v_axis']]] = v
    #     self.stage.move_pos_fast(*coords)
    #     #self.stage.move_pos_fast(x, y, None)
    #     #self.current_stage_pos_arrow.setPos(x, y)
    #     self.stage.settings.x_pos.read_from_hardware()
    #     self.stage.settings.y_pos.read_from_hardware()

    def move_position(self, h, v, dh, dv):
        self.stage.abs_mov(h, v)
        
    
class PIStage2DFrameSlowScan(BaseRaster2DFrameSlowScan):
    
    name = "PIStage2DFrameSlowScan"
    
    def __init__(self, app):
        BaseRaster2DFrameSlowScan.__init__(self, app, h_limits=(0,75), v_limits=(0,75), h_unit="um", v_unit="um")        
    
    def setup(self):
        PIStage2DSlowScan.setup(self)

    def move_position_start(self, h,v):
        PIStage2DSlowScan.move_position_start(self, h, v)
    
    # def move_position_slow(self, h,v, dh,dv):
    #     PIStage2DSlowScan.move_position_slow(self, h,v, dh,dv)
        
    # def move_position_fast(self,  h,v, dh,dv):
    #     PIStage2DSlowScan.move_position_fast(self,  h,v, dh,dv)
    
    # def on_new_stage_limits(self):
    #     PIStage2DSlowScan.on_new_stage_limits(self)

class Delay_PI_2DSlowScan(PIStage2DSlowScan):
    
    name = 'Delay_PI_2DSlowScan'

    def setup_figure(self):
        PIStage2DSlowScan.setup_figure(self)
        self.set_details_widget(widget=self.settings.New_UI(include=['h_axis', 'v_axis', 'pixel_time', 'frame_time']))


    def scan_specific_setup(self):
        #self.settings['pixel_time'] = 0.01
        self.settings.pixel_time.change_readonly(False)
        
    def collect_pixel(self, pixel_num, k, j, i):
        time.sleep(self.settings['pixel_time'])
        
    def post_scan_cleanup(self):
        pass
        
    def update_display(self):
        #MCLStage2DSlowScan.update_display(self)
        self.stage.read_from_hardware()
        self.stage.read_from_hardware()