# =============================================================================
# Top-6-Class MetaStudio (Top6Meta)
# Version: 1.0
#
# A Python HPC framework for modeling architected materials and metastructures.
#
# Authors:
#   - Agyapal Singh [1]
#   - Georgios Mermigkis [2]
#   - Panagiotis Hadjidoukas [2]
#   - Nikolaos Karathanasopoulos [1]
#
# Affiliations:
#   [1] New York University, Department of Engineering,
#       Abu Dhabi, United Arab Emirates
#   [2] Computer Engineering and Informatics Department,
#       University of Patras, Greece
#
# © 2026 The Authors
#
# License: MIT License
# =============================================================================
import sys
import numpy as np
import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QLabel, QPushButton, QComboBox,
                             QButtonGroup, QStackedWidget, QFormLayout, 
                             QLineEdit, QGridLayout, QCheckBox, QSlider,
                             QMessageBox, QDialog, QGroupBox, QFileDialog,
                             QProgressDialog, QFrame, QColorDialog, QDialogButtonBox,
                             QTableWidget, QTableWidgetItem, QHeaderView, QAction, QSplitter)
from PyQt5.QtGui import QFont, QPalette, QColor, QRegExpValidator, QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt, QRegExp, QTimer, QThread, pyqtSignal, QObject

sys.path.append('../python')
import tpms_core            # import tpms_core.py
import tpms_moments         # import tpms_moments.py
np.bool = np.bool_          # fix the bool type error (conda env problems)

from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.StlAPI import StlAPI_Reader
from OCC.Core.TopoDS import TopoDS_Shape
import os
import csv

import subprocess
import pickle
import tempfile
import webbrowser

# global variables for the TPMS interface 
# All the variables get default vales for each type using the function "update_configuration" in the TPMSInterface class
type = "TPMS"
form = "Lattice"
form_shape = "Cubic"                                        # Cubic or Cylindrical for TPMS, Spinodal
general_topo = "Sheet"
specific_topo = "Gyroid (GY)"
general_topo_for_hybrid = []                                # Hybrid has 1 topological for each layer
specific_topo_for_hybrid = []                               # Hybrid has 1 topological for each layer
x_rep_for_hybrid = []                                       # Hybrid has 1 repetition for each layer
y_rep_for_hybrid = []                                       # Hybrid has 1 repetition for each layer
z_repo_for_hybrid = []                                      # Hybrid has 1 repetition for each layer
vol_fraction_for_hybrid = []                                # Hybrid has 1 volume fraction for each layer
waves_number = 1000                                         # for spinodal  (waves range: 100-10000) - (dont have decimal)
number_of_layers = 2                                        # for hybrid and layered (hybrid 2 or 3, layered 2..9)
transition_quality = 5                                     # for hybrid (transition quality range: 1-40 - translate to "Low", "Medium", "High")
transition_location = "0.5"                                 # for hybrid - is a value between 0.2-0.8 with default value 0.5
layer_densities = []                                        # for layered - each layer has a different density, they have to sum up to 100 
                                                            # layer_density for each layer = volume_fraction for each layer
volume_fraction = 30                                        # for TPMS, Spinodal types (the range is 10-60% for all the cases)
grading = False                                             # for TPMS, Spinodal types (Hybrid has rading but it is not a boolean)
grading_for_hybrid = "Linear"                               # for Hybrid
grading_direction = "X"                                     # this & the following will be in a popup if grading is selected
grading_type = "Linear"                                     # 1=Linear, 2=cosine, 3=Sinusoidal, 4=VType Min at Centre, 5=VType Max at Centre     
grading_starting_volume = 30
grading_ending_volume = 40
length = 10                                                 # dimensions range: 0.2-1000 (radius too)
height = 10
width = 10
radius = 5                                                  # if form is "Cylindrical"
x_repetitions = 2                                           # repetitions range: 1-1000 (dont have decimal)
y_repetitions = 2
z_repetitions = 2
x_stretching = 1                                            # stretching range: 0.5-2
y_stretching = 1
z_stretching = 1
x_rotation = 0                                              # rotation range: 0-180 deg
y_rotation = 0
z_rotation = 0
resolution_points = 50                                      # resolution range: 50-1000
bottom_face_thickness = 4                                   # if form is "Sandwich"                                                 (btf range: 0-height/2)
top_face_thickness = 3                                      # if form is "Sandwich"                                                 (tft range: 0-height/2)
left_support_thick = 10                                     # if structure is "Beam"
right_support_thick = 10
bending_radius = 5                                          # StrutD: Grading = 0
stretching_radius = 5                                       
vertical_radius = 5                                         
joint_radius = 5                                            
min_bending_radius = 5                                      # StrutD: Grading = [1, 2, 3, 4, 5]
min_stretching_radius = 5                                  
min_vertical_radius = 5         
min_joint_radius = 5                                        
max_bending_radius = 10
max_stretching_radius = 10
max_vertical_radius = 10
max_joint_radius = 10
volume_fraction_strut = 20                                  # VolFracBased: Grading = 0
min_volume_fraction_strut = 20                              # VolFracBased: Grading = [1, 2, 3, 4, 5]
max_volume_fraction_strut = 40                                 
cylindrical_hybrid_type = "Triple-axis cylindrical (y, x, z)"   # For only Cylindrical in Hybrid: Type 1~ Cyl hybrid along one axis %2~ Cyl hybrid along 2 axes %%3~ Cyl hybrid along all three axiss
advanced_options = False                                    # Used in the second page to show/hide the advanced options
                                                            # Need a global variable to check if the values are filled when needed to
design_type = "Volume Fraction"                             # for strut: "VolFracBased" or "StrutDiaSizeBased", default is first   
model_color = (0.5059, 0.9686, 0.5490)
run_parallel = True
model_ipc_structure = False                                 # global variable for the Model IPC Structure checkbox

class TPMSInterface(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Top-6-Class Modeling Studio")
        self.setMinimumSize(1000, 800)  # Set a minimum size for the window because we don't have it fixed now
        #self.setGeometry(100, 100, 1400, 1000)
        
        # main central widget
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
        
        # title + stay updated button
        header_layout = QGridLayout()
        title_label = QLabel("Architected Materials and Structures (Top6Class)")
        title_label.setFont(QFont('Arial', 16))
        title_label.setAlignment(Qt.AlignCenter)

        stay_updated_btn = QPushButton("Stay updated →")
        stay_updated_btn.setFlat(True)
        stay_updated_btn.setCursor(Qt.PointingHandCursor)
        stay_updated_btn.setStyleSheet("""
            QPushButton {
                color: #ff8c00;
                border: none;
                background: transparent;
                font-size: 12px;
            }
            QPushButton:hover {
                color: #1e90ff;
                text-decoration: underline;
            }
        """)
        stay_updated_btn.clicked.connect(
            lambda: webbrowser.open("https://forms.gle/cArVMLCCcdhhZWzL9")
        )

        # grid placement
        header_layout.addWidget(title_label, 0, 1, alignment=Qt.AlignCenter)
        header_layout.addWidget(stay_updated_btn, 0, 2, alignment=Qt.AlignRight)

        # column stretch to force symmetry
        header_layout.setColumnStretch(0, 1)
        header_layout.setColumnStretch(1, 0)
        header_layout.setColumnStretch(2, 1)

        main_layout.addLayout(header_layout)

        # category Selection
        category_layout = QHBoxLayout()
        self.category_buttons = QButtonGroup()
        categories = ['Strut', 'TPMS', 'Spinodal', 'Hybrid', 'IPCs']
        
        self.category_buttons_dict = {}
        for category in categories:
            button = QPushButton(category)
            button.setCheckable(True)
            button.setStyleSheet("""
                QPushButton { 
                    background-color: #f0f8ff; 
                    border: 2px solid #87cefa; 
                    color: #1e90ff; 
                    padding: 10px;
                }
                QPushButton:checked {
                    background-color: #1e90ff;
                    color: white;
                }
            """)
            self.category_buttons.addButton(button)
            category_layout.addWidget(button)
            self.category_buttons_dict[category] = button
            button.clicked.connect(lambda checked, cat=category: self.update_configuration(cat))
        
        main_layout.addLayout(category_layout)
        
        # Box for colour selector + parallel checkbox
        top_box = QHBoxLayout()
        
        # colour selector
        self.color_selector = ColorSelectorWidget(self)
        top_box.addWidget(self.color_selector)

        # parallel checkbox
        self.parallel_checkbox = QCheckBox("Run in Parallel")
        self.parallel_checkbox.setChecked(run_parallel)
        self.parallel_checkbox.toggled.connect(self.toggle_parallel)
        top_box.addWidget(self.parallel_checkbox)

        main_layout.addLayout(top_box)

        # 3D Visualization Area and load the initial mesh
        self.plotter = pvqt.QtInteractor()
        initial_mesh = pv.read('defaults/tpms.stl')
        self.current_mesh = initial_mesh
        main_layout.addWidget(self.plotter.interactor)
        self.plotter.add_mesh(self.current_mesh, color=model_color, opacity=1, show_edges=True)
        self.plotter.add_axes()
        self.plotter.reset_camera()
        
        # configuration Area based on the selected category
        self.config_widget = QWidget()
        self.config_layout = QVBoxLayout()
        self.config_widget.setLayout(self.config_layout)
        main_layout.addWidget(self.config_widget)
        
        # proceed Button
        self.proceed_button = QPushButton("Proceed to Modeling")
        self.proceed_button.setStyleSheet("""
            background-color: #4682b4; 
            color: white; 
            padding: 15px; 
            font-size: 16px;
        """)
        self.proceed_button.clicked.connect(self.second_page)
        self.proceed_button.setEnabled(False)
        main_layout.addWidget(self.proceed_button)
        
        self.category_buttons_dict['TPMS'].setChecked(True)
        self.update_configuration('TPMS')

        # helper menus
        self._create_menus()
    
    
    # this function is used to update the configuration based on the selected category
    def update_configuration(self, category):
        
        # update the global "type" variable
        global type
        
        # original name was "Layered" but needed to change to "IPCs" for the button
        if category == "IPCs":
            category = "Layered"
            type = "Layered"
        else:
            type = category
        
        # also update the global variables for the form, general topo, specific topo (they need to be reset)
        global form, general_topo, specific_topo, number_of_layers, volume_fraction, waves_number
        global transition_quality, grading_for_hybrid, length, height, width, x_repetitions, y_repetitions, z_repetitions
        global grading, grading_direction, grading_type, grading_starting_volume, grading_ending_volume
        global bottom_face_thickness, top_face_thickness, resolution_points, radius
        global x_stretching, y_stretching, z_stretching, x_rotation, y_rotation, z_rotation
        global bending_radius, stretching_radius, vertical_radius, joint_radius, number_of_layers
        global general_topo_for_hybrid, specific_topo_for_hybrid, layer_densities, transition_location
        global x_rep_for_hybrid, y_rep_for_hybrid, z_repo_for_hybrid, vol_fraction_for_hybrid
        global form_shape, design_type, volume_fraction_strut, min_volume_fraction_strut, max_volume_fraction_strut
        global cylindrical_hybrid_type, model_ipc_structure
        global min_bending_radius, min_stretching_radius, min_vertical_radius, min_joint_radius
        global max_bending_radius, max_stretching_radius, max_vertical_radius, max_joint_radius
        global left_support_thick, right_support_thick
        
        form = "Lattice"
        form_shape = "Cubic"
        model_ipc_structure = False
        
        if category == 'Strut':
            general_topo = "Skeletal"
            design_type = "Volume Fraction"
            specific_topo = "Standard Octet Lattice"
            volume_fraction_strut = 20
            min_volume_fraction_strut = 20
            max_volume_fraction_strut = 40
            min_bending_radius = 5
            min_stretching_radius = 5
            min_vertical_radius = 5
            min_joint_radius = 5
            max_bending_radius = 10
            max_stretching_radius = 10
            max_vertical_radius = 10
            max_joint_radius = 10
            length = 10
            height = 10
            width = 10
            grading = False
            grading_direction = "Constant"
            grading_type = "Constant"
            bottom_face_thickness = 2
            top_face_thickness = 1
            left_support_thick = 10
            right_support_thick = 10
            resolution_points = 50
            bending_radius = 5
            stretching_radius = 5
            vertical_radius = 5
            joint_radius = 5
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
        elif category == "TPMS":
            general_topo = "Sheet"
            form_shape = "Cubic"
            specific_topo = "Gyroid (GY)"
            volume_fraction = 30
            grading = False
            grading_direction = "X"
            grading_type = "Linear"
            grading_starting_volume = 30
            grading_ending_volume = 40
            bottom_face_thickness = 4
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            length = 10
            height = 10
            width = 10
            radius = 5
            resolution_points = 50
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            x_stretching = 1
            y_stretching = 1
            z_stretching = 1
            x_rotation = 0
            y_rotation = 0
            z_rotation = 0
        elif category == "Spinodal":
            general_topo = "Sheet"
            form_shape = "Cubic"
            volume_fraction = 30
            waves_number = 1000
            grading = False
            grading_direction = "X"
            grading_type = "Linear"
            grading_starting_volume = 30
            grading_ending_volume = 40
            bottom_face_thickness = 4
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            length = 10
            height = 10
            width = 10
            radius = 5
            resolution_points = 121
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            x_stretching = 1
            y_stretching = 1
            z_stretching = 1
            x_rotation = 0
            y_rotation = 0
            z_rotation = 0
        elif category == "Hybrid":
            number_of_layers = 2
            transition_quality = 5
            grading_for_hybrid = "Linear"
            resolution_points = 50
            cylindrical_hybrid_type = "Triple-axis cylindrical (y, x, z)"
            transition_location = "0.5"
            length = 10
            height = 20
            width = 10
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            number_of_layers = 2
            general_topo_for_hybrid = ["Sheet", "Sheet"]
            specific_topo_for_hybrid = ["Gyroid (GY)", "Gyroid (GY)"]
            x_rep_for_hybrid = [1, 1]
            y_rep_for_hybrid = [2, 2]
            z_repo_for_hybrid = [2, 2]
            vol_fraction_for_hybrid = [30, 30]
        elif category == "Layered":
            general_topo = "Sheet"
            specific_topo = "Gyroid (GY)"
            number_of_layers = 2
            length = 10
            height = 10
            width = 10
            radius = 5
            x_repetitions = 3
            y_repetitions = 3
            z_repetitions = 3
            bottom_face_thickness = 4
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            resolution_points = 61
            number_of_layers = 2
            layer_densities = [70, 30]
        
        # clear previous configuration completely
        while self.config_layout.count():
            item = self.config_layout.takeAt(0)
            if item.widget():
                widget = item.widget()
                self.config_layout.removeWidget(widget)
                widget.deleteLater()
            elif item.layout():
                layout = item.layout()
                while layout.count():
                    sub_item = layout.takeAt(0)
                    if sub_item.widget():
                        sub_widget = sub_item.widget()
                        layout.removeWidget(sub_widget)
                        sub_widget.deleteLater()
                    elif sub_item.layout():
                        nested_layout = sub_item.layout()
                        while nested_layout.count():
                            nested_sub_item = nested_layout.takeAt(0)
                            if nested_sub_item.widget():
                                nested_sub_widget = nested_sub_item.widget()
                                nested_layout.removeWidget(nested_sub_widget)
                                nested_sub_widget.deleteLater()
        
        # explicitly reset the sliding-bars because they wouldnt go away
        if hasattr(self, 'volume_fraction_slider'):
            del self.volume_fraction_slider
        if hasattr(self, 'volume_fraction_label'):
            del self.volume_fraction_label
        if hasattr(self, 'grading_check'):
            del self.grading_check
        
        # clear the plot so I can load the new model based on the category
        self.plotter.clear()
        
        # disable all category buttons and enable the selected one
        for cat, button in self.category_buttons_dict.items():
            button.setChecked(cat == category)
        
        # create configuration based on selected category
        if category == 'Strut':
            mesh = pv.read('defaults/strut.stl')
            self.create_strut_configuration()
        elif category == 'TPMS':
            mesh = pv.read('defaults/tpms.stl')
            self.create_tpms_configuration()
        elif category == 'Spinodal':
            mesh = pv.read('defaults/spinodal.stl')
            self.create_spinodal_configuration()
        elif category == 'Hybrid':
            mesh = pv.read('defaults/hybrid.stl')
            self.create_hybrid_configuration()
        elif category == 'Layered':
            mesh1 = pv.read('defaults/layered_0.stl')
            mesh2 = pv.read('defaults/layered_1.stl')
            self.current_mesh = [mesh1, mesh2]
            
            # Add each mesh separately with different colors
            self.plotter.add_mesh(mesh1, color=model_color, opacity=1, show_edges=True, name='layer_1')
            darker_model_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
            self.plotter.add_mesh(mesh2, color=darker_model_color, opacity=1, show_edges=True, name='layer_2')
            self.create_layered_configuration()
                    
        # add the mesh to the plotter - in layered we don't do that because they are already added
        if category != 'Layered':
            self.current_mesh = mesh
            self.plotter.add_mesh(mesh, color=model_color, opacity=1, show_edges=True)
        self.plotter.reset_camera()
        
        # enable proceed button
        self.proceed_button.setEnabled(True)
    
    # ------------------------ Configuration Functions (setters) --------------------------------
    def update_form(self, text):
        global form
        form = text
        
        # based on the form selection, update the dimensions (if we need radius or not)
        if type == "Strut" or type == "Layered":
            self.setup_dimensions(text)
    
    def update_form_shape(self, text):
        global form_shape
        form_shape = text
    
    def update_general_topo(self, text):
        global general_topo
        general_topo = text
    
    def update_design_type(self, text):
        global design_type
        design_type = text
    
    def update_specific_topo(self, text):
        global specific_topo
        specific_topo = text
    
    def update_layers(self, text):
        global number_of_layers
        number_of_layers = text
    
    def update_volume_fraction(self, value):
        global volume_fraction
        volume_fraction = value
    
    def update_waves_number(self, text):
        global waves_number
        waves_number = text
    
    def update_distribution_hybrid(self, text):
        global grading_for_hybrid
        grading_for_hybrid = text
    
    def update_length(self, text):
        global length
        length = text
    
    def update_height(self, text):
        global height
        height = text
    
    def update_width(self, text):
        global width
        width = text
    
    def update_radius(self, text):
        global radius
        try:
            radius = int(text) if text else 0
        except ValueError:
            radius = 0
            
    def update_x_repetitions(self, text):
        global x_repetitions
        x_repetitions = text
    
    def update_y_repetitions(self, text):
        global y_repetitions
        y_repetitions = text
        
    def update_z_repetitions(self, text):
        global z_repetitions
        z_repetitions = text
        
    def is_grading_checked(self):
        global grading
        grading = not grading
        
        if grading:
            self.show_grading_popup()
            
        # if grading is enabled, we need to hide the volume fraction slider and label
        if hasattr(self, 'volume_fraction_slider'):
            self.volume_fraction_slider.setVisible(not grading)
            self.volume_fraction_label.setVisible(not grading)
    
    def is_model_ipc_checked(self):
        global model_ipc_structure
        model_ipc_structure = not model_ipc_structure
        
    def toggle_parallel(self):
        global run_parallel
        run_parallel = not run_parallel
    
    
    # the form selection is common for all categories
    def create_common_form_selection(self):
        form_layout = QHBoxLayout()
        if type == "Layered":
            form_label = QLabel("Design Specifications:")

        elif type == "Spinodal" or type == "TPMS":
            form_label = QLabel("Specimen Form:")
        else:
            form_label = QLabel("Form:")
        self.form_combo = QComboBox()
        self.form_combo.addItems([
            "Lattice", 
            # "Cylindrical",            # TPMS, SPIN seems to be using it but needs a different variable=SHAPE
                                        # Strut, Layred are not using it till now
            "Sandwich",
            "Beam or Plate"
        ])
        self.form_combo.currentTextChanged.connect(lambda text: self.update_form(text))
        form_layout.addWidget(form_label)
        form_layout.addWidget(self.form_combo)
        return form_layout
    
    def create_strut_configuration(self):
        
        # label
        strut_label = QLabel("<b>Specifications</b>")
        strut_label.setFont(QFont('Arial', 15))
        self.config_layout.addWidget(strut_label)
        
        # form Selection
        #self.config_layout.addLayout(self.create_common_form_selection())
        
        # topological Selection
        topo_layout = QGridLayout()
        
        # # general Topological Selection
        # general_topo_label = QLabel("General Topological Selection:")
        # self.general_topo_combo = QComboBox()
        # self.general_topo_combo.addItems(["Skeletal"])
        # self.general_topo_combo.activated.connect(lambda: self.update_general_topo(self.general_topo_combo.currentText()))
        # topo_layout.addWidget(general_topo_label, 0, 0)
        # topo_layout.addWidget(self.general_topo_combo, 0, 1)
        
        # Design Type Selection
        design_type_label = QLabel("Design Type:")
        self.design_type_combo = QComboBox()
        self.design_type_combo.addItems([
            "Strut Diameter",
            "Volume Fraction"
        ])
        self.design_type_combo.activated.connect(lambda: self.update_design_type(self.design_type_combo.currentText()))
        self.design_type_combo.setCurrentText("Volume Fraction")
        topo_layout.addWidget(design_type_label, 0, 0)
        topo_layout.addWidget(self.design_type_combo, 0, 1)
        
        # specific Topological Selection
        specific_topo_label = QLabel("Architecture:")
        self.specific_topo_combo = QComboBox()
        self.specific_topo_combo.addItems([
            "Standard Octet Lattice", 
            "Reinforced Octet", 
            "Octahedral", 
            "Reinforced Octahedral",
            "Circular Octahedral",
            "BCC"
        ])
        self.specific_topo_combo.activated.connect(lambda: self.update_specific_topo(self.specific_topo_combo.currentText()))
        topo_layout.addWidget(specific_topo_label, 1, 0)
        topo_layout.addWidget(self.specific_topo_combo, 1, 1)
        
        self.config_layout.addLayout(topo_layout)
        
        # dimensions are given in the first page for the strut
        # we create the dimension container here to add the dimensions because if form == "Cylindrical" we need radius
        # so we need to update the container based on the form selection - on the fly
        self.dim_container = QWidget()
        self.dim_layout = QFormLayout(self.dim_container)
        self.config_layout.addWidget(self.dim_container)
        
        
        # Horizontal layout for Material Grading and Model IPC Structure checkboxes
        grading_ipc_layout = QHBoxLayout()
        grading_ipc_layout.setSpacing(10)
        
        # Material Grading Checkbox
        self.grading_check = QCheckBox("Material Grading")
        self.grading_check.stateChanged.connect(lambda: self.is_grading_checked())
        grading_ipc_layout.addWidget(self.grading_check)
        
        # Model IPC Structure Checkbox
        self.model_ipc_check = QCheckBox("Model IPC Structure")
        self.model_ipc_check.stateChanged.connect(self.is_model_ipc_checked) 
        grading_ipc_layout.addWidget(self.model_ipc_check)
        # This prevents the checkboxes from filling the entire width, keeping them close.
        grading_ipc_layout.addStretch(1) 
        self.config_layout.addLayout(grading_ipc_layout)
        
        # initial dimension setup
        self.setup_dimensions("Lattice")  
    
    def create_tpms_configuration(self):
        
        # label
        tpms_label = QLabel("<b>Specifications</b>")
        tpms_label.setFont(QFont('Arial', 15))
        self.config_layout.addWidget(tpms_label)
        
        # form Selection
        self.config_layout.addLayout(self.create_common_form_selection())
        
        # form Shape Selection
        form_shape_layout = QHBoxLayout()
        form_shape_label = QLabel("Specimen Shape:")
        self.form_shape_combo = QComboBox()
        self.form_shape_combo.addItems(["Cubic", "Cylindrical"])
        self.form_shape_combo.currentTextChanged.connect(lambda text: self.update_form_shape(text))
        form_shape_layout.addWidget(form_shape_label)
        form_shape_layout.addWidget(self.form_shape_combo)
        self.config_layout.addLayout(form_shape_layout)
        
        # topological Selection
        topo_layout = QGridLayout()
        
        # general Topological Selection
        general_topo_label = QLabel("Topology Type:")
        self.general_topo_combo = QComboBox()
        self.general_topo_combo.addItems(["Sheet", "Skeletal", "Exoskeletal"])
        self.general_topo_combo.activated.connect(lambda: self.update_general_topo(self.general_topo_combo.currentText()))
        topo_layout.addWidget(general_topo_label, 0, 0)
        topo_layout.addWidget(self.general_topo_combo, 0, 1)
        
        # Specific Topological Selection
        specific_topo_label = QLabel("Topology Name:")
        self.specific_topo_combo = QComboBox()
        self.specific_topo_combo.addItems(["Gyroid (GY)", "IWP", "Pcell (SPC)", "FischerkochS (FKS)", "Schwarz_Diamond (SCD)", "Neovius (NE)", "Lidinoid (LD)",
                                            "SchwarzHexagonal (SCH)", "SplitP (SLP)", "I2Y (I2Y)", "Fisher–KochC(S) (FKCS)", "F-RD (FRD)"])
        self.specific_topo_combo.activated.connect(lambda: self.update_specific_topo(self.specific_topo_combo.currentText()))
        topo_layout.addWidget(specific_topo_label, 1, 0)
        topo_layout.addWidget(self.specific_topo_combo, 1, 1)
        
        self.config_layout.addLayout(topo_layout)
        
        # material Distribution
        dist_layout = QVBoxLayout()
        
        # Volume Fraction Slider and Label (FormLayout)
        volume_layout = QFormLayout()
        self.volume_fraction_slider = QSlider(Qt.Horizontal)
        self.volume_fraction_slider.setMinimum(10)
        self.volume_fraction_slider.setMaximum(60)
        self.volume_fraction_label = QLabel("Volume Fraction: 30%")
        self.volume_fraction_slider.setValue(30)
        
        self.volume_fraction_slider.valueChanged.connect(
            lambda: [self.volume_fraction_label.setText(f"Volume Fraction: {self.volume_fraction_slider.value()}%"),
            self.update_volume_fraction(self.volume_fraction_slider.value())]
        )
        
        volume_layout.addRow(self.volume_fraction_label, self.volume_fraction_slider)

        # Material Grading and Model IPC Structure Checkboxes (QHBoxLayout)
        grading_ipc_layout = QHBoxLayout()
        grading_ipc_layout.setSpacing(10)
        
        self.grading_check = QCheckBox("Material Grading")
        self.grading_check.stateChanged.connect(lambda: self.is_grading_checked())
        grading_ipc_layout.addWidget(self.grading_check)
        
        self.model_ipc_check = QCheckBox("Model IPC Structure")
        self.model_ipc_check.stateChanged.connect(self.is_model_ipc_checked) 
        grading_ipc_layout.addWidget(self.model_ipc_check)
        
        grading_ipc_layout.addStretch(1)

        # Add layouts in the desired order:
        self.config_layout.addLayout(volume_layout)
        self.config_layout.addLayout(grading_ipc_layout)
        self.config_layout.addLayout(dist_layout)
    
    def create_spinodal_configuration(self):
        
        # add a label
        spinodal_label = QLabel("<b>Specifications</b>")
        spinodal_label.setFont(QFont('Arial', 15))
        self.config_layout.addWidget(spinodal_label)

        # form Selection
        self.config_layout.addLayout(self.create_common_form_selection())
        
        # form Shape Selection
        form_shape_layout = QHBoxLayout()
        form_shape_label = QLabel("Specimen Shape:")
        self.form_shape_combo = QComboBox()
        self.form_shape_combo.addItems(["Cubic", "Cylindrical"])
        self.form_shape_combo.currentTextChanged.connect(lambda text: self.update_form_shape(text))
        form_shape_layout.addWidget(form_shape_label)
        form_shape_layout.addWidget(self.form_shape_combo)
        self.config_layout.addLayout(form_shape_layout)
        
        # topological Selection
        topo_layout = QHBoxLayout()
        topo_label = QLabel("Topology Type:")
        self.topo_combo = QComboBox()
        self.topo_combo.addItems(["Sheet", "Skeletal", "Exoskeletal"])
        self.topo_combo.activated.connect(lambda: self.update_general_topo(self.topo_combo.currentText()))
        
        waves_layout = QFormLayout()
        self.waves_input = QLineEdit()
        self.waves_input.setText("1000")
        self.waves_input.setValidator(QRegExpValidator(QRegExp("[0-9]*")))
        self.waves_input.textChanged.connect(lambda text: self.update_waves_number(text))
        waves_layout.addRow("Waves Number:", self.waves_input)
        
        topo_layout.addWidget(topo_label)
        topo_layout.addWidget(self.topo_combo)
        
        self.config_layout.addLayout(topo_layout)
        self.config_layout.addLayout(waves_layout)
        
        # Material Distribution
        dist_layout = QVBoxLayout()
        
        # Volume Fraction Slider and Label (FormLayout)
        volume_layout = QFormLayout()
        self.volume_fraction_slider = QSlider(Qt.Horizontal)
        self.volume_fraction_slider.setMinimum(10)
        self.volume_fraction_slider.setMaximum(60)
        self.volume_fraction_label = QLabel("Volume Fraction: 30%")
        self.volume_fraction_slider.setValue(30)
        
        self.volume_fraction_slider.valueChanged.connect(
            lambda: [self.volume_fraction_label.setText(f"Volume Fraction: {self.volume_fraction_slider.value()}%"),
            self.update_volume_fraction(self.volume_fraction_slider.value())]
        )
        
        volume_layout.addRow(self.volume_fraction_label, self.volume_fraction_slider)
        
        # Horizontal layout for Material Grading and Model IPC Structure checkboxes
        grading_ipc_layout = QHBoxLayout()
        grading_ipc_layout.setSpacing(10)
        
        self.grading_check = QCheckBox("Material Grading")
        self.grading_check.stateChanged.connect(lambda: self.is_grading_checked())
        grading_ipc_layout.addWidget(self.grading_check)
        
        self.model_ipc_check = QCheckBox("Model IPC Structure")
        self.model_ipc_check.stateChanged.connect(self.is_model_ipc_checked) 
        grading_ipc_layout.addWidget(self.model_ipc_check)
        
        grading_ipc_layout.addStretch(1)
        
        # Add layouts in the desired order (volume first, then checkboxes)
        self.config_layout.addLayout(volume_layout)
        self.config_layout.addLayout(grading_ipc_layout)
        
        # The original dist_layout is added last.
        self.config_layout.addLayout(dist_layout)
    
    def create_hybrid_configuration(self):
        
        # topological Selection
        topo_layout = QFormLayout()
        
        # For now we only support 2 layers, that's why we have it commented out
        # self.layers_input = QLineEdit()
        # self.layers_input.setText("2")
        # self.layers_input.setValidator(QRegExpValidator(QRegExp("^[2-3]$")))
        # self.layers_input.textChanged.connect(lambda text: self.update_layers(text))
        
        # topo_layout.addRow("Number of Layers (2-3):", self.layers_input)
        # self.config_layout.addLayout(topo_layout)
        
        # dimensions are given in the first page for the hybrid (that was missing in the first page)
        # we create the dimension container here to add the dimensions because if form == "Cylindrical" we need radius
        # so we need to update the container based on the form selection - on the fly
        # we dont actually need the Cylindrical option here, but we already have it
        self.dim_container = QWidget()
        self.dim_layout = QFormLayout(self.dim_container)
        self.config_layout.addWidget(self.dim_container)
        
        # initial dimension setup
        self.setup_dimensions("Lattice")  
        
        # material Distribution
        dist_layout = QHBoxLayout()
        dist_layout.setSpacing(6)

        dist_label = QLabel("Material Distribution:")
        dist_label.setSizePolicy(QLabel().sizePolicy().horizontalPolicy(), dist_label.sizePolicy().verticalPolicy())
        dist_label.setMinimumWidth(130)

        self.dist_combo = QComboBox()
        self.dist_combo.addItems(["Linear", "Cylindrical", "Spherical"])
        self.dist_combo.setMinimumWidth(150)
        self.dist_combo.activated.connect(lambda: self.update_distribution_hybrid(self.dist_combo.currentText()))

        # Add left-side elements
        dist_layout.addWidget(dist_label)
        dist_layout.addWidget(self.dist_combo)
        
        # stretch to push the element to the right
        dist_layout.addStretch(1) 
        
        # "Model IPC Structure" checkbox to the right
        self.model_ipc_check = QCheckBox("Model IPC Structure")
        self.model_ipc_check.stateChanged.connect(self.is_model_ipc_checked) 
        dist_layout.addWidget(self.model_ipc_check)
        
        self.config_layout.addLayout(dist_layout)
    
    def create_layered_configuration(self):
        # form Selection
        #self.config_layout.addLayout(self.create_common_form_selection())
        
        # topological Selection
        topo_layout = QGridLayout()
        
        # general Topological Selection
        general_topo_label = QLabel("Topology Type:")
        self.general_topo_combo = QComboBox()
        self.general_topo_combo.addItems(["Sheet", "Skeletal", "Exoskeletal"])
        self.general_topo_combo.activated.connect(lambda: self.update_general_topo(self.general_topo_combo.currentText()))
        topo_layout.addWidget(general_topo_label, 0, 0)
        topo_layout.addWidget(self.general_topo_combo, 0, 1)
        
        # Specific Topological Selection
        specific_topo_label = QLabel("Topology Name:")
        self.specific_topo_combo = QComboBox()
        self.specific_topo_combo.addItems(["Gyroid (GY)", "IWP", "Pcell (SPC)", "FischerkochS (FKS)", "Schwarz_Diamond (SCD)", "Neovius (NE)",
                                            "Lidinoid (LD)", "SchwarzHexagonal (SCH)", "SplitP (SLP)", "I2Y (I2Y)", "Fisher–KochC(S) (FKCS)", "F-RD (FRD)"])
        self.specific_topo_combo.activated.connect(lambda: self.update_specific_topo(self.specific_topo_combo.currentText()))
        topo_layout.addWidget(specific_topo_label, 1, 0)
        topo_layout.addWidget(self.specific_topo_combo, 1, 1)
        
        # Number of Layers - We stick with 2 layers for now
        # layers_label = QLabel("Number of Layers (2-9):")
        # self.layers_input = QLineEdit()
        # self.layers_input.setText("2")
        # self.layers_input.setValidator(QRegExpValidator(QRegExp("^[2-9]$")))
        # self.layers_input.textChanged.connect(lambda text: self.update_layers(text))
        # topo_layout.addWidget(layers_label, 2, 0)
        # topo_layout.addWidget(self.layers_input, 2, 1)
        
        self.config_layout.addLayout(topo_layout)
        
        # in the layered configuration, dimensions are given in the first page        
        # we create the dimension container here to add the dimensions because if form == "Cylindrical" we need radius
        # so we need to update the container based on the form selection - on the fly
        self.dim_container = QWidget()
        self.dim_layout = QFormLayout(self.dim_container)
        self.config_layout.addWidget(self.dim_container)
        
        # initial dimension setup
        self.setup_dimensions("Lattice")  
    
    
    # ------------- Proceed Function, Popup for empty values, Popup for grading ---------------
    # this is used because if form == "Cylindrical" we need to add radius to the dimensions
    def setup_dimensions(self, form_type):
        # clear existing layout
        while self.dim_layout.count():
            child = self.dim_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
                
        # we need to default the values for the dimensions again
        # otherwise when going from one form to another, the values will be kept and some may be empty in the new form
        global length, height, width, radius

        if form_type == "Cylindrical":
            length = 10
            height = 10
            width = 10
            radius = 5
            
            self.radius_input = QLineEdit()
            self.radius_input.setText("5")
            self.radius_input.setValidator(QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.radius_input.textChanged.connect(lambda text: self.update_radius(text))
            
            self.height_input = QLineEdit()
            self.height_input.setText("10")
            self.height_input.setValidator(QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))
            
            self.dim_layout.addRow("Radius (mm):", self.radius_input)
            self.dim_layout.addRow("Height (mm):", self.height_input)
            
        else:

            if type == "Hybrid" or type == "Layered" or type == "Strut":
                # add label
                layers_label = QLabel("<b>Total Specimen Dimensions</b>")
                layers_label.setFont(QFont('Arial', 15))
                self.dim_layout.addRow(layers_label)
            
            length=10
            width = 10

            self.length_input = QLineEdit()
            self.length_input.setText("10")
            self.length_input.setValidator(QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.length_input.textChanged.connect(lambda text: self.update_length(text))
            
            self.height_input = QLineEdit()
            # Hybrid has 20 as default height
            if type == "Hybrid":
                height = 20
                self.height_input.setText("20")
            else:
                height = 10
                self.height_input.setText("10")

            self.height_input.setValidator(QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))
            
            self.width_input = QLineEdit()
            self.width_input.setText("10")
            self.width_input.setValidator(QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.width_input.textChanged.connect(lambda text: self.update_width(text))
            
            self.dim_layout.addRow("Length (mm):", self.length_input)
            self.dim_layout.addRow("Width (mm):", self.width_input)
            self.dim_layout.addRow("Height (mm):", self.height_input)
    
    def popup_empty_values(self):
        msg = QMessageBox()
        msg.setWindowTitle("Empty Values")
        msg.setText("Please fill in all the required values.")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
    
    def popup_wrong_range_values(self, variable):
        msg = QMessageBox()
        msg.setWindowTitle("Invalid Range")
        
        if variable == "WavesNumber":
            msg.setText(f"The value for {variable} is not within the valid range (100-10000).")
        elif variable == "Height":
            msg.setText(f"The value of Height, Width, Length, Radius should be between 0.2 and 1000.")
        
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()

    def popup_beam_cylindrical_not_allowed(self):
        msg = QMessageBox()
        msg.setWindowTitle("Invalid Configuration")
        msg.setText("The combination of 'Beam or Plate' form with 'Cylindrical' shape is not allowed for Spinodal structures.")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
    
    def show_grading_popup(self):
        
        # Grading is available in TPMS, SPIN and Strut.
        dialog = GradingDialog(parent=self, global_var=type)
        if dialog.exec_() == QDialog.Accepted:
            values = dialog.get_values()
            global grading_direction, grading_type, grading_starting_volume, grading_ending_volume
            
            if type == "TPMS":
                grading_direction = values['direction']
                grading_type = values['type']
                grading_starting_volume = values['starting_volume']
                grading_ending_volume = values['ending_volume']
            elif type == "Spinodal":
                grading_direction = values['direction']
                grading_type = values['type']
                grading_starting_volume = values['starting_volume']
                grading_ending_volume = values['ending_volume']
            elif type == "Strut":
                grading_direction = values['direction']
                grading_type = values['type']
    
    def second_page(self):
        
        global length, height, width, radius, bottom_face_thickness, top_face_thickness, left_support_thick, right_support_thick
        global resolution_points, bending_radius, stretching_radius, vertical_radius, joint_radius, x_repetitions, y_repetitions, z_repetitions
        global x_stretching, y_stretching, z_stretching, x_rotation, y_rotation, z_rotation, transition_location, transition_quality
        global waves_number, number_of_layers
        global volume_fraction_strut, min_volume_fraction_strut, max_volume_fraction_strut
        global min_bending_radius, min_stretching_radius, min_vertical_radius, min_joint_radius
        global max_bending_radius, max_stretching_radius, max_vertical_radius, max_joint_radius
        
        # check if all the required values are filled in and also check their range
        # set the default values of the second page variables, because if we leave 1 value empty and we go back to the first page
        # and then to the second page, the value will be kept and the user will not know that it is empty
        if type == "TPMS":
            length = 10
            height = 10
            width = 10
            radius = 5
            bottom_face_thickness = 4
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            resolution_points = 50
            x_stretching = 1
            y_stretching = 1
            z_stretching = 1
            x_rotation = 0
            y_rotation = 0
            z_rotation = 0

            if form == "Beam or Plate" and form_shape == "Cylindrical":
                self.popup_beam_cylindrical_not_allowed()
                return
            
        elif type == "Spinodal":
            
            length = 10
            height = 10
            width = 10
            radius = 5
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            resolution_points = 121
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            bottom_face_thickness = 4
            
            if not waves_number:
                self.popup_empty_values()
                return
            
            if int(waves_number) < 100 or int(waves_number) > 10000:
                self.popup_wrong_range_values("WavesNumber")
                return

            if form == "Beam or Plate" and form_shape == "Cylindrical":
                self.popup_beam_cylindrical_not_allowed()
                return
                    
        elif type == "Hybrid":
            
            resolution_points = 50
            transition_location = "0.5"
            transition_quality = 5
            
            if not number_of_layers or not length or not height or not width:
                self.popup_empty_values()
                return
            
            if float(height) < 0.2 or float(height) > 1000 or float(width) < 0.2 or float(width) > 1000 or float(length) < 0.2 or float(length) > 1000:
                self.popup_wrong_range_values("Height")
                return
            
        elif type == "Layered":
            
            x_repetitions = 3
            y_repetitions = 3
            z_repetitions = 3
            resolution_points = 61
            bottom_face_thickness = 4
            top_face_thickness = 3
            left_support_thick = 10
            right_support_thick = 10
            
            if form == "Cylindrical":
                if not number_of_layers or radius == '' or height == '': 
                    self.popup_empty_values()
                    return
                
                # use the inputs as checking otherwise we would get 1 strange error once in a while
                # height covers all height, width, length, radius - for the pop-up
                height = float(self.height_input.text().strip())
                radius = float(self.radius_input.text().strip())
                
                if height < 0.2 or height > 1000 or radius < 0.2 or radius > 1000:
                    self.popup_wrong_range_values("Height")
                    return
                
            else:
                if not number_of_layers or not length or not height or not width:
                    self.popup_empty_values()
                    return
                
                height = float(self.height_input.text().strip())
                width = float(self.width_input.text().strip())
                length = float(self.length_input.text().strip())
                
                if height < 0.2 or height > 1000 or width < 0.2 or width > 1000 or length < 0.2 or length > 1000:
                    self.popup_wrong_range_values("Height")
                    return
             
        elif type == "Strut":
            
            # set the default values again for the variables of the second page
            bottom_face_thickness = 2
            top_face_thickness = 1
            left_support_thick = 10
            right_support_thick = 10
            resolution_points = 50
            bending_radius = 5
            stretching_radius = 5
            vertical_radius = 5
            joint_radius = 5
            x_repetitions = 2
            y_repetitions = 2
            z_repetitions = 2
            volume_fraction_strut = 20
            min_volume_fraction_strut = 20
            max_volume_fraction_strut = 40
            min_bending_radius = 5
            min_stretching_radius = 5                                  
            min_vertical_radius = 5         
            min_joint_radius = 5                                        
            max_bending_radius = 10
            max_stretching_radius = 10
            max_vertical_radius = 10
            max_joint_radius = 10 
            
            
            if form == "Cylindrical":
                if radius == '' or height == '':
                    self.popup_empty_values()
                    return
                
                # use the inputs as checking otherwise we would get 1 strange error once in a while
                # height covers all height, width, length, radius - for the pop-up
                height = float(self.height_input.text().strip())
                radius = float(self.radius_input.text().strip())
                
                if height < 0.2 or height > 1000 or radius < 0.2 or radius > 1000:
                    self.popup_wrong_range_values("Height")
                    return

            else:
                if not length or not height or not width:
                    self.popup_empty_values()
                    return
                
                height = float(self.height_input.text().strip())
                width = float(self.width_input.text().strip())
                length = float(self.length_input.text().strip())
                
                if length < 0.2 or length > 1000 or height < 0.2 or height > 1000 or width < 0.2 or width > 1000:
                    self.popup_wrong_range_values("Height")
                    return
        
                        
        # open the second page
        # if a second window exists, close and delete it
        try:
            # create a fresh instance
            self._second_window = SecondPage()
            self._second_window.show()
            self.hide()
        except Exception as e:
            print(f"Error creating second window: {e}")
        #self.hide()

    # helper menus
    def _create_menus(self):
        menubar = self.menuBar()

        # Application menu
        app_menu = menubar.addMenu("Top-6-Class MetaStudio")

        about_action = QAction("About Top-6-Class MetaStudio", self)
        about_action.triggered.connect(self.show_about)

        app_menu.addAction(about_action)
    
    def show_about(self):
        dialog = AboutDialog(self)
        dialog.exec_()


    # are you sure you want to close the application?
    def closeEvent(self, event):
        # create custom message box
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle("Exit Application")
        msg_box.setText("Are you sure you want to close the application?")
        
        # set "warning" icon
        msg_box.setIcon(QMessageBox.Warning)
        
        # create custom buttons
        yes_button = msg_box.addButton("Yes, Exit", QMessageBox.YesRole)
        no_button = msg_box.addButton("Cancel", QMessageBox.NoRole)
        
        # apply custom styling
        msg_box.setStyleSheet("""
            QMessageBox {
                background-color: #f0f0f0;
            }
            QMessageBox QLabel {
                color: #2c3e50;
                font-size: 12pt;
                padding: 10px;
            }
            QPushButton {
                background-color: #3498db;
                color: white;
                padding: 8px 16px;
                border-radius: 4px;
                min-width: 100px;
                font-size: 10pt;
            }
            QPushButton:hover {
                background-color: #2980b9;
            }
            QPushButton[text="Cancel"] {
                background-color: #95a5a6;
            }
            QPushButton[text="Cancel"]:hover {
                background-color: #7f8c8d;
            }
        """)
        
        # show the message box and handle response
        msg_box.exec_()
        
        if msg_box.clickedButton() == yes_button:
            event.accept()
        else:
            event.ignore()

    # Called when the window is shown (including when returning from second page)
    def showEvent(self, event):
        
        super().showEvent(event)
        # Refresh the color picker display and plots when window is shown
        if hasattr(self, 'color_selector'):
            self.color_selector.refresh_display()



class KeybindsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("PyVista Keybinds")
        self.setMinimumSize(450, 550)
        self.setFixedSize(450, 550)
        layout = QVBoxLayout()
        
        # Create table for keybinds
        table = QTableWidget()
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Key", "Action"])
        
        # PyVista default keybinds
        keybinds = [
            ("q", "Close the rendering window"),
            ("f", "Focus on point"),
            ("v", "Isometric camera view"),
            ("w", "Switch to wireframe representation"),
            ("s", "Switch to surface representation"),
            ("r", "Reset camera to view all actors"),
            ("Shift+Click", "Rubber band selection"),
            ("p", "Pick a point/cell"),
            ("Left Click + Drag", "Rotate camera"),
            ("Shift + Left Click + Drag", "Pan camera"),
            ("Middle Click + Drag", "Pan camera"),
            ("Right Click + Drag", "Zoom camera"),
            ("Mouse Wheel", "Zoom camera"),
            ("Ctrl+Left Click", "Pick point"),
            ("Up/Down Arrow", "Zoom in/out"),
            ("Left/Right Arrow", "Rotate camera"),
        ]
        
        table.setRowCount(len(keybinds))
        
        for row, (key, action) in enumerate(keybinds):
            key_item = QTableWidgetItem(key)
            action_item = QTableWidgetItem(action)
            table.setItem(row, 0, key_item)
            table.setItem(row, 1, action_item)
        
        # Style the table
        table.setStyleSheet("""
            QTableWidget {
                background-color: white;
                border: 2px solid #87cefa;
                border-radius: 5px;
                gridline-color: #e0e0e0;
            }
            QTableWidget::item {
                padding: 8px;
                color: #333;
            }
            QTableWidget::item:selected {
                background-color: #1e90ff;
                color: white;
            }
            QHeaderView::section {
                background-color: #4682b4;
                color: white;
                padding: 8px;
                border: none;
                font-weight: bold;
                font-size: 14px;
            }
        """)
        
        # Adjust column widths
        header = table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        table.setSelectionBehavior(QTableWidget.SelectRows)
        
        # Disable scrolling
        table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        
        layout.addWidget(table)
        self.setLayout(layout)
        
        # Style the dialog
        self.setStyleSheet("""
            QDialog {
                background-color: #f0f8ff;
            }
        """)

class ColorSelectorWidget(QWidget):
    # Class variable to store the last selected color name across all instances
    _last_selected_color = 'Light Blue'
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.colors = {
            'Light Blue': (0.53, 0.81, 0.98),
            'Ocean Blue': (0.27, 0.51, 0.71),
            'Forest Green': (0.13, 0.55, 0.13),
            'Coral Red': (1.0, 0.5, 0.31),
            'Royal Purple': (0.5, 0.0, 0.5),
            'Golden Yellow': (1.0, 0.84, 0.0),
            'Silver Gray': (0.75, 0.75, 0.75),
            'Rose Pink': (1.0, 0.41, 0.71),
            'Custom Color': None                    # Special entry for color picker
        }
        
        # Track if we have a custom color active
        self.has_custom_color = False
        self.custom_color_value = None
        
        # Update the class variable with current color state
        ColorSelectorWidget._last_selected_color = self.get_last_selected_colour()
        self.init_ui()
    
    def get_last_selected_colour(self):
        """Find which color name corresponds to the current global model_color"""
        global model_color
        
        # Check if current model_color matches any predefined color
        for color_name, color_value in self.colors.items():
            if color_value is not None:  # Skip 'Custom Color'
                # Compare with small tolerance for floating point comparison
                if (abs(model_color[0] - color_value[0]) < 0.001 and 
                    abs(model_color[1] - color_value[1]) < 0.001 and 
                    abs(model_color[2] - color_value[2]) < 0.001):
                    return color_name
        
        # If no match found, it's a custom color
        self.has_custom_color = True
        self.custom_color_value = model_color
        return 'Custom Color'
    
    def init_ui(self):
        layout = QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)
        
        color_label = QLabel("Model Color:")
        color_label.setStyleSheet("color: #1e90ff;")
        
        self.color_combo = QComboBox()
        self.color_combo.addItems(self.colors.keys())
        
        # Set current selection
        current_color_name = self.get_last_selected_colour()
        self.color_combo.setCurrentText(current_color_name)
        
        self.color_combo.setStyleSheet("""
            QComboBox {
                background-color: #f0f8ff;
                border: 2px solid #87cefa;
                border-radius: 5px;
                padding: 5px;
                min-width: 150px;
                color: #1e90ff;
            }
            QComboBox:hover {
                border: 2px solid #1e90ff;
            }
            QComboBox::drop-down {
                border: none;
            }
            QComboBox::down-arrow {
                image: none;
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid #1e90ff;
                margin-right: 5px;
            }
            QComboBox QAbstractItemView {
                background-color: white;
                border: 2px solid #87cefa;
                selection-background-color: #1e90ff;
                selection-color: white;
            }
        """)
        
        self.color_combo.currentTextChanged.connect(self.on_color_selection_changed)
        # Also connect to activated signal to catch when the same item is selected again
        self.color_combo.activated.connect(self.on_color_activated)
        
        layout.addWidget(color_label)
        layout.addWidget(self.color_combo)
        layout.addStretch()
        self.setLayout(layout)
    
    def open_color_picker(self):
        """Open color picker dialog and return RGB tuple (0-1 range)"""
        global model_color
        color_dialog = QColorDialog(self)
        color_dialog.setOption(QColorDialog.ShowAlphaChannel, False)
        
        # Set initial color to current model_color
        initial_color = QColor(
            int(model_color[0] * 255),
            int(model_color[1] * 255),
            int(model_color[2] * 255)
        )
        color_dialog.setCurrentColor(initial_color)
        
        if color_dialog.exec_() == QColorDialog.Accepted:
            color = color_dialog.selectedColor()
            # Convert from 0-255 range to 0-1 range for PyVista
            rgb_tuple = (color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0)
            return rgb_tuple
        return None
    
    def on_color_activated(self, index):
        """Handle when user clicks on an item (even if it's the same as current)"""
        color_name = self.color_combo.itemText(index)
        if color_name == "Custom Color":
            self.handle_custom_color()
    
    def on_color_selection_changed(self, color_name):
        """Handle when the selection actually changes"""
        if color_name != "Custom Color":
            self.has_custom_color = False
            self.custom_color_value = None
            ColorSelectorWidget._last_selected_color = color_name
            self.update_color(color_name)
        else:
            # Don't call handle_custom_color() to avoid the double-dialog issue.
            # The activated signal handler already does this.
            if self.custom_color_value:
                self.update_mesh_color(self.custom_color_value)
    
    def handle_custom_color(self):
        """Handle custom color selection"""
        custom_color = self.open_color_picker()
        
        if custom_color:
            self.has_custom_color = True
            self.custom_color_value = custom_color
            ColorSelectorWidget._last_selected_color = 'Custom Color'
            self.update_mesh_color(custom_color)
            
            # Ensure combo shows "Custom Color"
            self.color_combo.blockSignals(True)
            self.color_combo.setCurrentText("Custom Color")
            self.color_combo.blockSignals(False)
        else:
            # User cancelled - revert combo box
            self.refresh_combo_display()
    
    def update_color(self, color_name):
        """Update the colour based on the selected option"""
        global model_color
        model_color = self.colors[color_name]
        self.update_mesh_color(model_color)
    
    def refresh_display(self):
        """Refresh the display to match current global model_color"""
        global model_color
        self.update_mesh_color(model_color)
        self.refresh_combo_display()
    
    def refresh_combo_display(self):
        """Refresh the combo box to show the correct current selection"""
        current_color_name = self.get_last_selected_colour()
        
        self.color_combo.blockSignals(True)
        self.color_combo.setCurrentText(current_color_name)
        self.color_combo.blockSignals(False)
        
        ColorSelectorWidget._last_selected_color = current_color_name
    
    def update_mesh_color(self, color):
        """Update the mesh color in the plotter"""
        global model_color, type, model_ipc_structure
        model_color = color
        
        # find the main window by traversing up the parent hierarchy
        widget = self
        while widget.parent() is not None:
            widget = widget.parent()
            if isinstance(widget, TPMSInterface):
                if hasattr(widget, 'plotter'):
                    if hasattr(widget, 'current_mesh'):
                        widget.plotter.clear()
                        
                        # Check if we're dealing with layered meshes (list of meshes)
                        # This now includes both Layered type AND IPC structure cases
                        if (type == "Layered" or model_ipc_structure == True) and isinstance(widget.current_mesh, list):
                            base_color = color
                            
                            # Create variations of the base color for each layer
                            colors = [
                                base_color,  
                                (base_color[0] * 0.7, base_color[1] * 0.7, base_color[2] * 0.7) # darker
                            ]
                            
                            for i, mesh in enumerate(widget.current_mesh):
                                layer_color = colors[i] if i < len(colors) else base_color
                                widget.plotter.add_mesh(
                                    mesh,
                                    color=layer_color,
                                    opacity=1,
                                    show_edges=True,
                                    name=f'layer_{i}'
                                )
                        else:
                            # Handle single mesh
                            widget.plotter.add_mesh(
                                widget.current_mesh,
                                color=color,
                                opacity=1,
                                show_edges=True
                            )
                        
                        widget.plotter.add_axes()
                        widget.plotter.reset_camera()
                        return
                    else:
                        print("No current_mesh found in TPMSInterface")
                else:
                    print("No plotter found in TPMSInterface")
                break
            elif isinstance(widget, SecondPage):
                if hasattr(widget, 'plotter'):
                    if hasattr(widget, 'current_mesh'):
                        widget.plotter.clear()
                        
                        # Check if we're dealing with layered meshes (list of meshes)
                        # This now includes both Layered type AND IPC structure cases
                        if (type == "Layered" or model_ipc_structure == True) and isinstance(widget.current_mesh, list):
                            base_color = color
                            
                            # Create variations of the base color for each layer
                            colors = [
                                base_color,
                                (base_color[0] * 0.7, base_color[1] * 0.7, base_color[2] * 0.7)  # darker
                            ]
                            
                            for i, mesh in enumerate(widget.current_mesh):
                                layer_color = colors[i] if i < len(colors) else base_color
                                widget.plotter.add_mesh(
                                    mesh,
                                    color=layer_color,
                                    opacity=1,
                                    show_edges=True,
                                    name=f'layer_{i}'
                                )
                            
                            # Update individual layer plotters if they exist
                            if hasattr(widget, 'plotter0') and hasattr(widget, 'plotter1'):
                                if hasattr(widget, 'current_mesh_reinf') and hasattr(widget, 'current_mesh_compl'):
                                    # Clear and update plotter0 (reinforcement layer)
                                    widget.plotter0.clear()
                                    widget.plotter0.add_mesh(
                                        widget.current_mesh_reinf,
                                        color=base_color,
                                        opacity=1,
                                        show_edges=True
                                    )
                                    widget.plotter0.add_axes()
                                    widget.plotter0.reset_camera()
                                    
                                    # Clear and update plotter1 (compliant layer)
                                    widget.plotter1.clear()
                                    darker_color = (base_color[0] * 0.7, base_color[1] * 0.7, base_color[2] * 0.7)
                                    widget.plotter1.add_mesh(
                                        widget.current_mesh_compl,
                                        color=darker_color,
                                        opacity=1,
                                        show_edges=True
                                    )
                                    widget.plotter1.add_axes()
                                    widget.plotter1.reset_camera()
                        else:
                            # Handle single mesh
                            widget.plotter.add_mesh(
                                widget.current_mesh,
                                color=color,
                                opacity=1,
                                show_edges=True
                            )
                        
                        widget.plotter.add_axes()
                        widget.plotter.reset_camera()
                        return
                    else:
                        print("No current_mesh found in SecondPage")
                else:
                    print("No plotter found in SecondPage")
                break

class GradingDialog(QDialog):
    
    # Grading is different of TPMS(==Spin) and Strut
    def __init__(self, parent=None, global_var=None):
        super().__init__(parent)
        self.setWindowTitle("Material Grading Parameters")
        self.setModal(True)
        
        # initialize default values
        self.grading_direction = "X"
        self.grading_type = "Linear"
        self.grading_starting_volume = "30"
        self.grading_ending_volume = "40"
        
        # create the form layout
        form_layout = QFormLayout()
        
        # different UI based on type
        if global_var == "TPMS" or global_var == "Spinodal":
            
            # Direction ComboBox
            self.direction_combo = QComboBox()
            self.direction_combo.addItems(["X", "Y", "Z"])
            self.direction_combo.setCurrentText(self.grading_direction)
            form_layout.addRow("Grading Direction:", self.direction_combo)
            
            # Type ComboBox
            self.type_combo = QComboBox()
            self.type_combo.addItems(["Linear", "Cosine", "Sinusoidal", "VType Min at Centre", "VType Max at Centre"])
            self.type_combo.setCurrentText(self.grading_type)
            form_layout.addRow("Grading Type:", self.type_combo)
            
            # Starting Volume Entry
            self.starting_volume = QLineEdit(self.grading_starting_volume)
            self.starting_volume.setValidator(QRegExpValidator(QRegExp("[0-9]*")))
            form_layout.addRow("Starting Volume (%):", self.starting_volume)
            
            # Ending Volume Entry
            self.ending_volume = QLineEdit(self.grading_ending_volume)
            self.ending_volume.setValidator(QRegExpValidator(QRegExp("[0-9]*")))
            form_layout.addRow("Ending Volume (%):", self.ending_volume)
            
        elif global_var == "Strut":
            # Only Direction and Type for Strut
            
            # Direction ComboBox with different options
            self.direction_combo = QComboBox()
            self.direction_combo.addItems(["X", "Y", "Z"])
            self.direction_combo.setCurrentText("X")  # Default for Strut
            form_layout.addRow("Direction:", self.direction_combo)
            
            # Type ComboBox with numerical descriptions
            self.type_combo = QComboBox()
            self.type_combo.addItems(["Linear", "Cosine", "Sinusoidal", "VType Min at Centre", "VType Max at Centre"])
            self.type_combo.setCurrentText("Linear")  # Default for Strut
            form_layout.addRow("Type:", self.type_combo)
            
        # buttons
        button_layout = QHBoxLayout()
        self.ok_button = QPushButton("OK")
        self.cancel_button = QPushButton("Cancel")
        button_layout.addWidget(self.ok_button)
        button_layout.addWidget(self.cancel_button)
        
        # main layout
        main_layout = QVBoxLayout()
        main_layout.addLayout(form_layout)
        main_layout.addLayout(button_layout)
        self.setLayout(main_layout)
        
        # connect buttons
        self.ok_button.clicked.connect(self.validate_and_accept)
        self.cancel_button.clicked.connect(self.reject)
        
        # set size - adjust size for Strut which has fewer fields
        if global_var == "Strut":
            self.setFixedSize(300, 150)
        else:
            self.setFixedSize(300, 200)
        
        # style
        self.setStyleSheet("""
            QDialog {
                background-color: rgb(240, 248, 255);
            }
            QLabel {
                color: rgb(70, 130, 180);
            }
            QLineEdit, QComboBox {
                background-color: white;
                border: 1px solid rgb(70, 130, 180);
                border-radius: 3px;
                padding: 2px;
            }
            QPushButton {
                background-color: rgb(70, 130, 180);
                color: white;
                border: none;
                border-radius: 3px;
                padding: 5px 15px;
            }
            QPushButton:hover {
                background-color: rgb(100, 149, 237);
            }
        """)
        
    def validate_and_accept(self):
        # For Strut, no validation needed as there are no text inputs
        # For TPMS/Spinodal, validate the volume entries
        if hasattr(self, 'starting_volume') and hasattr(self, 'ending_volume'):
            if not self.starting_volume.text() or not self.ending_volume.text():
                QMessageBox.warning(self, "Input Error", "Please fill in all the required values.")
                return
        
        self.accept()
        
    def get_values(self):
        # Return different values based on what fields we have
        if hasattr(self, 'starting_volume') and hasattr(self, 'ending_volume'):
            # TPMS or Spinodal case
            return {
                'direction': self.direction_combo.currentText(),
                'type': self.type_combo.currentText(),
                'starting_volume': self.starting_volume.text(),
                'ending_volume': self.ending_volume.text()
            }
        else:
            # Strut case
            return {
                'direction': self.direction_combo.currentText(),
                'type': self.type_combo.currentText()
            }

class MomentsDialog(QDialog):
    def __init__(self, moments=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Averaged Cross-Sectional Properties")
        self.setModal(True)
        self.setFixedSize(800, 500)
        
        # main layout with margins for better spacing
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(20, 10, 20, 20)
        main_layout.setSpacing(10)
        self.setLayout(main_layout)
        
        # add title
        title_label = QLabel("Effective Property Computations")
        title_label.setStyleSheet("""
            font-size: 18px;
            font-weight: bold;
            color: #2c3e50;
            padding: 10px;
        """)
        title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(title_label)
        
        # create a frame for the grid
        frame = QFrame()
        frame.setFrameStyle(QFrame.StyledPanel | QFrame.Raised)
        frame.setStyleSheet("""
            QFrame {
                background-color: white;
                border: 1px solid #dcdde1;
                border-radius: 5px;
            }
        """)
        frame_layout = QVBoxLayout()
        frame.setLayout(frame_layout)
        
        # create grid for displaying moment data
        grid_widget = QWidget()
        grid_layout = QGridLayout()
        grid_layout.setSpacing(1)
        grid_widget.setLayout(grid_layout)
        
        # styles
        header_style = """
            font-weight: bold;
            font-size: 14px;
            padding: 12px;
            background-color: #4682b4;
            color: white;
            border: 1px solid #3873a3;
            border-radius: 3px;
        """
        
        value_style_even = """
            padding: 10px;
            background-color: #f8f9fa;
            border: 1px solid #e9ecef;
            border-radius: 3px;
        """
        
        value_style_odd = """
            padding: 10px;
            background-color: #ffffff;
            border: 1px solid #e9ecef;
            border-radius: 3px;
        """
        
        # add headers
        headers = [
            'Plane',
            'Area (mm²)',
            'Ixx (mm⁴)',
            'Iyy (mm⁴)',
            'Izz (mm⁴)',
            'Ixy/xz/yz (mm⁴)'
        ]
        
        for col, header in enumerate(headers):
            label = QLabel(header)
            label.setStyleSheet(header_style)
            label.setAlignment(Qt.AlignCenter)
            grid_layout.addWidget(label, 0, col)
        
        # add data rows
        if moments:
            planes = {'xy': 1, 'xz': 2, 'yz': 3}
            for plane, row in planes.items():
                # choose style based on row number
                current_style = value_style_odd if row % 2 else value_style_even
                
                # plane name
                plane_label = QLabel(plane.upper())
                plane_label.setStyleSheet(current_style)
                plane_label.setAlignment(Qt.AlignCenter)
                grid_layout.addWidget(plane_label, row, 0)
                
                # area
                area = moments[plane]['Area']
                area_label = QLabel(f"{area:.3f}")
                area_label.setStyleSheet(current_style)
                area_label.setAlignment(Qt.AlignRight)
                grid_layout.addWidget(area_label, row, 1)
                
                # moments of inertia
                if plane == 'xy':
                    ixx, iyy = moments[plane]['Ixx'], moments[plane]['Iyy']
                    izz, ixy = "-", moments[plane]['Ixy']
                elif plane == 'xz':
                    ixx, iyy = moments[plane]['Ixx'], "-"
                    izz, ixy = moments[plane]['Izz'], moments[plane]['Ixz']
                else:  # yz plane
                    ixx, iyy = "-", moments[plane]['Iyy']
                    izz, ixy = moments[plane]['Izz'], moments[plane]['Iyz']
                
                # add moment values
                for col, value in enumerate([ixx, iyy, izz, ixy], start=2):
                    if isinstance(value, (float, np.float64)):
                        text = f"{value:.3f}"
                    else:
                        text = value
                    label = QLabel(text)
                    label.setStyleSheet(current_style)
                    label.setAlignment(Qt.AlignRight)
                    grid_layout.addWidget(label, row, col)
        
        # add grid to frame
        frame_layout.addWidget(grid_widget)
        main_layout.addWidget(frame)
        
        # add description
        desc_label = QLabel("Values shown in mm² for areas and mm⁴ for moments of inertia")
        desc_label.setStyleSheet("""
            color: #666;
            font-style: italic;
            padding: 5px;
        """)
        desc_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(desc_label)

        # add asterisk note regarding voxel element analysis
        asterisk_label = QLabel("<i>* Computations use voxel element analysis</i>")
        asterisk_label.setStyleSheet("""
            color: #666;
            font-style: italic;
            padding: 2px;
        """)
        asterisk_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(asterisk_label)
        
        # add close button
        close_button = QPushButton("Close")
        close_button.setFixedWidth(120)
        close_button.setStyleSheet("""
            QPushButton {
                background-color: #4682b4;
                color: white;
                padding: 12px;
                font-size: 14px;
                border-radius: 5px;
                border: none;
            }
            QPushButton:hover {
                background-color: #3873a3;
            }
            QPushButton:pressed {
                background-color: #2f6491;
            }
        """)
        close_button.clicked.connect(self.accept)
        
        # add save data button
        save_button = QPushButton("Save Data")
        save_button.setFixedWidth(120)
        save_button.setStyleSheet("""
            QPushButton {
                background-color: #4682b4;
                color: white;
                padding: 12px;
                font-size: 14px;
                border-radius: 5px;
                border: none;
            }
            QPushButton:hover {
                background-color: #3873a3;
            }
            QPushButton:pressed {
                background-color: #2f6491;
            }
        """)
        save_button.clicked.connect(lambda: self.save_data(moments))
        
        # create button container for center alignment
        button_container = QHBoxLayout()
        button_container.addStretch()
        button_container.addWidget(close_button)
        button_container.addWidget(save_button)
        button_container.addStretch()
        main_layout.addLayout(button_container)
    
    
    # save the moments data to a .txt file
    def save_data(self, moments):
        print("Saving data...")

        file_path, _ = QFileDialog.getSaveFileName(self, "Save Moments Data", "", "Text Files (*.txt)")
        if file_path:
            with open(file_path, 'w') as file:
                file.write("Plane, Area (mm²), Ixx (mm⁴), Iyy (mm⁴), Izz (mm⁴), Ixy/xz/yz (mm⁴)\n")
                for plane, data in moments.items():
                    print(plane, data)
                    area = f"{data.get('Area', 0):.6f}"
                    ixx = f"{data.get('Ixx', '-'):.6f}" if isinstance(data.get("Ixx"), (int, float)) else "-"
                    iyy = f"{data.get('Iyy', '-'):.6f}" if isinstance(data.get("Iyy"), (int, float)) else "-"
                    izz = f"{data.get('Izz', '-'):.6f}" if isinstance(data.get("Izz"), (int, float)) else "-"

                    if plane == "xy":
                        ixy = f"{data.get('Ixy', '-'):.6f}" if isinstance(data.get("Ixy"), (int, float)) else "-"
                    elif plane == "xz":
                        ixy = f"{data.get('Ixz', '-'):.6f}" if isinstance(data.get("Ixz"), (int, float)) else "-"
                    elif plane == "yz":
                        ixy = f"{data.get('Iyz', '-'):.6f}" if isinstance(data.get("Iyz"), (int, float)) else "-"
                    else:
                        ixy = "-"

                    print(f"{plane}, {area}, {ixx}, {iyy}, {izz}, {ixy}")
                    file.write(f"{plane}, {area}, {ixx}, {iyy}, {izz}, {ixy}\n")

            print(f"Data saved to {file_path}")
            QMessageBox.information(self, "Data Saved", f"Data saved to {file_path}")

class GenerationProcess(QObject):
    #Drop-in replacement for GenerationThread using subprocess instead of threading.
    
    finished = pyqtSignal(object)  # emit a tuple or dict with results (F,V ...)
    error = pyqtSignal(str)
    stopped = pyqtSignal()
    
    def __init__(self, option, params):
        super().__init__()
        self.option = option
        self.params = params
        self.process = None
        self.stop_requested = False
        self.temp_dir = None
        self.check_timer = QTimer()
        self.check_timer.timeout.connect(self._check_process)
        
    def start(self):
        print("Starting generation process...")
        
        # Create temporary directory for communication
        self.temp_dir = tempfile.mkdtemp(prefix="generation_")
        
        # Serialize parameters
        params_file = os.path.join(self.temp_dir, "params.pkl")
        with open(params_file, 'wb') as f:
            pickle.dump((self.option, self.params), f)
        
        # Start subprocess
        script_path = os.path.join(os.path.dirname(__file__), "generation_worker.py")
        self.process = subprocess.Popen([
            sys.executable, script_path, 
            params_file, 
            self.temp_dir
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Start checking process status
        self.check_timer.start(100)  # Check every 100ms
        
    def stop(self):
        print("Stopping generation process...")
        self.stop_requested = True
        if self.process and self.process.poll() is None:
            self.process.terminate()
            # Force kill if it doesn't terminate nicely
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
        self._cleanup()
        
    def _check_process(self):
        # Check if process has finished and handle results
        if self.process is None:
            return
            
        if self.process.poll() is not None:  # Process has finished
            self.check_timer.stop()
            
            if self.stop_requested:
                self.stopped.emit()
                self._cleanup()
                return
            
            # check for results
            result_file = os.path.join(self.temp_dir, "result.pkl")
            error_file = os.path.join(self.temp_dir, "error.txt")
            
            if os.path.exists(result_file):
                try:
                    with open(result_file, 'rb') as f:
                        result = pickle.load(f)
                    self.finished.emit(result)
                except Exception as e:
                    self.error.emit(f"Failed to load result: {str(e)}")
            elif os.path.exists(error_file):
                try:
                    with open(error_file, 'r') as f:
                        error_msg = f.read()
                    self.error.emit(error_msg)
                except Exception as e:
                    self.error.emit(f"Process failed: {str(e)}")
            else:
                # check subprocess stderr
                stderr = self.process.stderr.read().decode() if self.process.stderr else ""
                if stderr:
                    self.error.emit(f"Process error: {stderr}")
                else:
                    self.error.emit("Process completed but no result found")
            
            self._cleanup()
    
    def _cleanup(self):
        # Clean up temporary files and resources
        if self.temp_dir and os.path.exists(self.temp_dir):
            try:
                import shutil
                shutil.rmtree(self.temp_dir)
            except:
                pass
        self.temp_dir = None
        self.process = None
        
    def isRunning(self):
        # check if the process is still running
        return self.process is not None and self.process.poll() is None


class SecondPage(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Architected Material Modeling - Further Configuration")

        if type == "Hybrid" and model_ipc_structure == True:
            self.setMinimumSize(1400, 920)      # larger because we need the expansion for cyl option + IPC
        elif type == "Strut" and model_ipc_structure == True and design_type == "Strut Diameter" and grading == True:
            self.setMinimumSize(1400, 920)
        else:
            self.setMinimumSize(1400, 800)      # min size because we dont have it fixed and visual problems can occur
        #self.setFixedSize(1400, 820)
        
        # main central widget
        main_widget = QWidget()
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
        
        # status bar used for thread messages when running the generation
        self.status_bar = self.statusBar()
        self.status_bar.setStyleSheet("background-color: #f0f0f0; color: #333; font-size: 16px;")
        
        # split layout into left (configuration) and right (visualization)
        split_layout = QHBoxLayout()
        
        # left side - Configuration
        config_widget = QWidget()
        self.config_layout = QVBoxLayout()
        config_widget.setLayout(self.config_layout)
        config_widget.setMaximumWidth(600)
        
        # create a wrapper widget/layout for the config side to handle the stretch
        # this stretch only affects the left side
        config_wrapper = QVBoxLayout()
        config_wrapper.addWidget(config_widget)
        config_wrapper.addStretch()  
        
        # right side - 3D Visualization
        self.plotter = pvqt.QtInteractor()
        self.current_mesh = None  # will hold the current mesh for color updates
        self.current_mesh_reinf = None
        self.current_mesh_compl = None

        # right side - Colour picker + label for the main plotter
        top_right_layout = QHBoxLayout()
        self.main_plotter_label = QLabel("Main 3D Model Visualization")
        self.main_plotter_label.setStyleSheet("color: #1e90ff; font-size: 16px;")
        top_right_layout.addWidget(self.main_plotter_label)
        
        # Small help button right next to the label
        self.questions_icon = QPushButton("Help (?)")
        self.questions_icon.setStyleSheet("""
            QPushButton {
                color: #1e90ff;
                font-size: 11px;
                font-weight: bold;
                background-color: #f0f8ff;
                border: 2px solid #87cefa;
                border-radius: 4px;
                padding: 3px 8px;
            }
            QPushButton:hover {
                background-color: #e6f3ff;
                border: 2px solid #1e90ff;
            }
            QPushButton:pressed {
                background-color: #d0e8ff;
                border: 2px solid #4682b4;
            }
        """)
        self.questions_icon.setToolTip("View keyboard shortcuts")
        self.questions_icon.clicked.connect(self.show_keybinds_dialog)
        
        top_right_layout.addWidget(self.questions_icon)
        top_right_layout.addStretch()  # Push color selector to the right
        
        self.color_selector = ColorSelectorWidget(self)
        self.color_selector.setFixedWidth(self.color_selector.sizeHint().width())
        top_right_layout.addWidget(self.color_selector, alignment=Qt.AlignRight)
        
        # labels for total surface, total vol fraction
        right_side_layout = QVBoxLayout()
        right_side_layout.addLayout(top_right_layout)
        right_side_layout.addWidget(self.plotter.interactor, stretch=1)
        
        # tolerance badge at the top right corner
        tolerance_badge = QLabel("Volume Tolerance: 1e-2")
        tolerance_badge.setStyleSheet("""
            background-color: rgba(70, 130, 180, 15%);
            color: #4682b4;
            border: 1px solid #4682b4;
            border-radius: 8px;
            padding: 4px 8px;
            font-size: 12px;
            font-weight: bold;
        """)
        tolerance_badge.setAlignment(Qt.AlignCenter)
        right_side_layout.addWidget(tolerance_badge, alignment=Qt.AlignRight | Qt.AlignTop)

        # info layout which is part of the right side (surface, vol fraction, moments button)
        info_layout = QHBoxLayout()
        
        # info panel: left side (total surface, total vol fraction)
        labels_layout = QVBoxLayout()
        self.total_surface_label = QLabel("Total Surface Area: 0.0 mm²")
        self.total_vol_fraction_label = QLabel("Total Volume Fraction: 0.0%")
        labels_layout.addWidget(self.total_surface_label)
        labels_layout.addWidget(self.total_vol_fraction_label)
        
        # info panel: right side (moments button) - at first have it disabled because we dont have moments yet
        button_layout = QVBoxLayout()
        
        # check if we need to add the moments button (Layred or IPC_Y)
        is_standard_layered = (type == "Layered" and model_ipc_structure == False)
        is_excluded_ipc_type = (
            (type == "TPMS" or type == "Strut" or type == "Spinodal" or type == "Hybrid") and model_ipc_structure == True
        )
        if not (is_standard_layered or is_excluded_ipc_type):
            self.moments_button = QPushButton("Effective Cross-Sectional Properties")
            self.moments_button.setFixedWidth(260)  
            self.moments_button.setStyleSheet("""
                QPushButton:enabled {
                    background-color: #4682b4;
                    color: white;
                    padding: 10px;
                    font-size: 14px;
                }
                QPushButton:disabled {
                    background-color: rgba(70, 130, 180, 40%);
                    color: rgba(255, 255, 255, 70%);
                    padding: 10px;
                    font-size: 14px;
                }
            """)
            self.moments_button.setToolTip("Enable this by completing the generation process first")
            self.moments_button.clicked.connect(lambda: self.show_moments())
            self.moments_button.setEnabled(False)
            
            # variable that will store the moments and will be used to pass them to the moments dialog
            self.moments = None
        
            button_layout.addWidget(self.moments_button)
        button_layout.addStretch() 
        
        # add the labels and button to the info layout, and then add the info layout to the right side layout
        info_layout.addLayout(labels_layout)
        info_layout.addStretch()  
        info_layout.addLayout(button_layout)
        right_side_layout.addLayout(info_layout)
        
        # add the wrapped config layout and plotter to the split layout
        split_layout.addLayout(config_wrapper)
        split_layout.addLayout(right_side_layout)
        main_layout.addLayout(split_layout)
        
        # bottom buttons
        button_layout = QHBoxLayout()
        self.generate_button = QPushButton("Generate")
        self.generate_button.clicked.connect(self.generate)
        self.save_button = QPushButton("Save CAD")
        self.save_button.clicked.connect(self.save_cardfile)
        
        self.generate_button.setStyleSheet("""
            background-color: #4682b4; 
            color: white; 
            padding: 15px; 
            font-size: 16px;
        """)
        self.save_button.setStyleSheet("""
            background-color: #4682b4; 
            color: white; 
            padding: 15px; 
            font-size: 16px;
        """)
        
        button_layout.addWidget(self.generate_button)
        button_layout.addWidget(self.save_button)
        
        main_layout.addLayout(button_layout)
        
        # initialize configuration based on type
        self.setup_configuration()
        
        # load initial STL model
        self.load_initial_model()


                
    # this function is used to load the initial 3D model
    def load_initial_model(self):
        if type == "TPMS":
            model_path = "defaults/tpms.stl"
            grid = pv.read(model_path)
            self.current_mesh = grid  # store current mesh for color updates
            self.plotter.add_mesh(grid, color=model_color, opacity=1, show_edges=True)
        elif type == "Strut":
            model_path = "defaults/strut.stl"
            grid = pv.read(model_path)
            self.current_mesh = grid  # store current mesh for color updates
            self.plotter.add_mesh(grid, color=model_color, opacity=1, show_edges=True)
        elif type == "Spinodal":
            model_path = "defaults/spinodal.stl"
            grid = pv.read(model_path)
            self.current_mesh = grid  # store current mesh for color updates
            self.plotter.add_mesh(grid, color=model_color, opacity=1, show_edges=True)
        elif type == "Hybrid":
            model_path = "defaults/hybrid.stl"
            grid = pv.read(model_path)
            self.current_mesh = grid  # store current mesh for color updates
            self.plotter.add_mesh(grid, color=model_color, opacity=1, show_edges=True)
        elif type == "Layered":
            model_path_0 = "defaults/layered_0.stl"
            model_path_1 = "defaults/layered_1.stl"
            
            # Load both meshes
            grid_0 = pv.read(model_path_0)
            grid_1 = pv.read(model_path_1)
            
            # store current meshes as a list for color updates
            self.current_mesh = [grid_0, grid_1]
            self.current_mesh_reinf = grid_0
            self.current_mesh_compl = grid_1
            
            # Add each mesh with different colors to main plotter
            self.plotter.add_mesh(grid_0, color=model_color, opacity=1, show_edges=True, name='layer_0')
            darker_model_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
            self.plotter.add_mesh(grid_1, color=darker_model_color, opacity=1, show_edges=True, name='layer_1')

            # Also load individual layers if plotters exist
            if hasattr(self, 'plotter0') and hasattr(self, 'plotter1'):
                self.load_individual_layers()
                
        # Only add axes and reset camera for main plotter if not layered (since layered handles it above)
        self.plotter.add_axes()
        self.plotter.reset_camera()
    
    # this function is used to load the 2 layers in layered
    def load_individual_layers(self, model_path_0="defaults/layered_0.stl", model_path_1="defaults/layered_1.stl"):
        if hasattr(self, 'plotter0') and hasattr(self, 'plotter1'):
            # Load layer 0 into plotter0
            grid_0 = pv.read(model_path_0)
            self.plotter0.add_mesh(grid_0, color=model_color, opacity=1, show_edges=True)
            self.plotter0.reset_camera()
            
            # Load layer 1 into plotter1 with darker color
            grid_1 = pv.read(model_path_1)
            darker_model_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
            self.plotter1.add_mesh(grid_1, color=darker_model_color, opacity=1, show_edges=True)
            self.plotter1.reset_camera()
    
    # this function is used to setup the configuration based on the selected type
    def setup_configuration(self):
        
        while self.config_layout.count():
            item = self.config_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        
        # add title
        if type == "Layered":
            title_label = QLabel("IPCs Configuration")
        else:
            title_label = QLabel(f"{type} Configuration")
        title_label.setFont(QFont('Arial', 17))
        title_label.setAlignment(Qt.AlignCenter)
        self.config_layout.addWidget(title_label)
        
                    
        # if form == "Sandwich"
        if form == "Sandwich":
            self.sandwitch_config()
            self.config_layout.addSpacing(10)
            
        if form == "Beam or Plate":
            self.beam_plate_config()
            self.config_layout.addSpacing(10)
        
        # add appropriate configuration based on type
        if type == 'TPMS':
            self.setup_tpms_config()
        elif type == 'Strut':
            self.setup_strut_config()
        elif type == 'Spinodal':
            self.setup_spinodal_config()
        elif type == 'Hybrid':
            self.setup_hybrid_config()
        elif type == 'Layered':
            self.setup_layered_config()

    # this function is used to create an input field with a label
    def create_input_field(self, label, default_value="", validator=None):
        layout = QHBoxLayout()
        label_widget = QLabel(label)
        input_widget = QLineEdit()
        input_widget.setText(str(default_value))
        if validator:
            input_widget.setValidator(validator)
        layout.addWidget(label_widget)
        layout.addWidget(input_widget)
        return layout, input_widget
    
    # this function is used to create a section label
    def create_section_label(self, text):
        label = QLabel(text)
        label.setFont(QFont('Arial', 16))
        label.setStyleSheet("color: #4682b4;")
        return label

    # this function is used to show the keybinds dialog
    def show_keybinds_dialog(self):
        dialog = KeybindsDialog(self)
        dialog.exec_()
        self.questions_icon.clearFocus()

    # this function takes the V: vertices, F: faces, after the generation and creates a PyVista mesh
    # that is then used to display the mesh in the PyVista plotter of the second page
    def create_pyvista_mesh(self, F, V):
        
        vertices = np.array(V, dtype=np.float32)
        faces = np.array(F, dtype=np.int32)
        
        # PyVista expects faces in format [n, id1, id2, id3] where n is number of points (3 for triangles)
        # we need to insert the count (3) at the start of each face
        faces_with_count = np.column_stack((np.full(len(faces), 3), faces))
        faces_with_count = faces_with_count.flatten()
        
        # create the PyVista mesh
        mesh = pv.PolyData(vertices, faces_with_count)
        return mesh
    
    # this function is used to update the mesh in the GUI after generation
    def update_mesh(self, F0, V0, F1=None, V1=None):
        
        # clear existing mesh
        self.plotter.clear_actors()
            
        # IPCs, so F0,V0 and F1,V1
        if model_ipc_structure == True:
        
            # create two new meshes from F0,V0 and F1,V1
            mesh1 = self.create_pyvista_mesh(F0, V0)
            mesh2 = self.create_pyvista_mesh(F1, V1)
            self.current_mesh = [mesh1, mesh2]
            self.current_mesh_compl = mesh1
            self.current_mesh_reinf = mesh2
            
            # add both meshes to the plotter with different colors
            self.plotter.add_mesh(mesh1, color=model_color,  opacity=1, show_edges=True)
            darker_model_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
            self.plotter.add_mesh(mesh2, color=darker_model_color,  opacity=1, show_edges=True)
            self.plotter.reset_camera()
            self.plotter.render()
            
            if hasattr(self, "plotters_container"):
                self.plotters_container.setVisible(True)
        
            if hasattr(self, 'plotter0') and hasattr(self, 'plotter1'):
                self.plotter0.clear_actors()
                self.plotter1.clear_actors()
                
                self.plotter0.add_mesh(mesh1, color=model_color, opacity=1, show_edges=True)
                self.plotter0.reset_camera()
                self.plotter0.render()
                
                self.plotter1.add_mesh(mesh2, color=darker_model_color, opacity=1, show_edges=True)
                self.plotter1.reset_camera()
                self.plotter1.render()
        
        # no IPCs, so single F,V
        else:
            
            # create new mesh from F, V (F,V come from generation)
            new_mesh = self.create_pyvista_mesh(F0, V0)
            self.current_mesh = new_mesh
            
            # add the new mesh to the plotter
            self.plotter.add_mesh(new_mesh, color=model_color,  opacity=1, show_edges=True)
            self.plotter.reset_camera()
            self.plotter.render()
    
    # this function is used to update the mesh in the GUI after generation for layered
    def update_mesh_layered(self, F1, V1, F2, V2):
        
        #clear existing meshes
        self.plotter.clear_actors()
        
        # create two new meshes from F1,V1 and F2,V2
        mesh1 = self.create_pyvista_mesh(F1, V1)
        mesh2 = self.create_pyvista_mesh(F2, V2)
        self.current_mesh = [mesh1, mesh2]
        self.current_mesh_compl = mesh1
        self.current_mesh_reinf = mesh2
        
        # add both meshes to the plotter with different colors
        self.plotter.add_mesh(mesh1, color=model_color,  opacity=1, show_edges=True)
        darker_model_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
        self.plotter.add_mesh(mesh2, color=darker_model_color,  opacity=1, show_edges=True)
        self.plotter.reset_camera()
        self.plotter.render()
        
        # update the individual layer plotters if they exist
        if hasattr(self, 'plotter0') and hasattr(self, 'plotter1'):
            self.plotter0.clear_actors()
            self.plotter1.clear_actors()
            
            self.plotter0.add_mesh(mesh1, color=model_color, opacity=1, show_edges=True)
            self.plotter0.reset_camera()
            self.plotter0.render()
            
            self.plotter1.add_mesh(mesh2, color=darker_model_color, opacity=1, show_edges=True)
            self.plotter1.reset_camera()
            self.plotter1.render()

    
    # ------------------------ Update Functions for variables (setters) --------------------------------
    def update_bottom_face_thickness(self, text):
        global bottom_face_thickness
        bottom_face_thickness = text    

    def update_top_face_thickness(self, text):
        global top_face_thickness
        top_face_thickness = text
    
    def update_left_support_thickness(self, text):
        global left_support_thick
        left_support_thick = text
        
    def update_right_support_thickness(self, text):
        global right_support_thick
        right_support_thick = text
    
    def update_length(self, text):
        global length
        length = text
    
    def update_height(self, text):
        global height
        height = text
    
    def update_width(self, text):
        global width
        width = text
    
    def update_radius(self, text):
        global radius
        radius = text
    
    def update_resolution(self, text):
        global resolution_points
        resolution_points = text
    
    def update_cylindrical_hybrid_type(self, text):
        global cylindrical_hybrid_type
        cylindrical_hybrid_type = text
        
    def update_x_repetitions(self, text):
        global x_repetitions
        x_repetitions = text
        
    def update_y_repetitions(self, text):
        global y_repetitions
        y_repetitions = text
        
    def update_z_repetitions(self, text):
        global z_repetitions
        z_repetitions = text
    
    def update_x_stretch(self, text):
        global x_stretching
        x_stretching = text
    
    def update_y_stretch(self, text):
        global y_stretching
        y_stretching = text
        
    def update_z_stretch(self, text):
        global z_stretching
        z_stretching = text
        
    def update_x_rotation(self, text):
        global x_rotation
        x_rotation = text
    
    def update_y_rotation(self, text):
        global y_rotation
        y_rotation = text
    
    def update_z_rotation(self, text):
        global z_rotation
        z_rotation = text
    
    def update_volume_fraction_strut(self, text):
        global volume_fraction_strut
        volume_fraction_strut = text
    
    def update_min_volume_fraction_for_strut(self, text):
        global min_volume_fraction_strut
        min_volume_fraction_strut = text
    
    def update_max_volume_fraction_for_strut(self, text):
        global max_volume_fraction_strut
        max_volume_fraction_strut = text
    
    def update_bending_radius(self, text):
        global bending_radius
        bending_radius = text
    
    def update_bending_radius_min(self, text):
        global min_bending_radius
        min_bending_radius = text
    
    def update_bending_radius_max(self, text):
        global max_bending_radius
        max_bending_radius = text
        
    def update_stretching_radius(self, text):
        global stretching_radius
        stretching_radius = text
    
    def update_stretching_radius_min(self, text):
        global min_stretching_radius
        min_stretching_radius = text
    
    def update_stretching_radius_max(self, text):
        global max_stretching_radius
        max_stretching_radius = text
    
    def update_vertical_radius(self, text):
        global vertical_radius
        vertical_radius = text
    
    def update_vertical_radius_min(self, text):
        global min_vertical_radius
        min_vertical_radius = text
        
    def update_vertical_radius_max(self, text):
        global max_vertical_radius
        max_vertical_radius = text
    
    def update_joint_radius(self, text):
        global joint_radius
        joint_radius = text
    
    def update_joint_radius_min(self, text):
        global min_joint_radius
        min_joint_radius = text
    
    def update_joint_radius_max(self, text):
        global max_joint_radius
        max_joint_radius = text
        
    def update_general_topo_for_hybrid(self, text, index):
        global general_topo_for_hybrid
        general_topo_for_hybrid[index-1] = text
    
    def update_specific_topo_for_hybrid(self, text, index):
        global specific_topo_for_hybrid
        specific_topo_for_hybrid[index-1] = text
        
    def update_x_rep_for_hybrid(self, text, index):
        global x_rep_for_hybrid
        x_rep_for_hybrid[index-1] = text
        
    def update_y_rep_for_hybrid(self, text, index):
        global y_rep_for_hybrid
        y_rep_for_hybrid[index-1] = text
    
    def update_z_rep_for_hybrid(self, text, index):
        global z_rep_for_hybrid
        z_rep_for_hybrid[index-1] = text
    
    def update_vol_fraction_for_hybrid(self, text, index):
        global vol_fraction_for_hybrid
        vol_fraction_for_hybrid[index-1] = text
    
    def update_transition_location(self, value):
        global transition_location
        transition_location = value
    
    def update_layer_density(self, text, index):
        global layer_densities
        try:
            # convert the text to integer before storing
            value = int(text) if text else 0
            layer_densities[index - 1] = value
        except ValueError:
            layer_densities[index - 1] = 0
    
    def update_transition_quality(self, value):
        global transition_quality
        transition_quality = value
    
    def update_transition_quality_input_from_slider(self):
        """Update the input field when slider changes"""
        # Block signals to prevent infinite loop
        self.transition_quality_input.blockSignals(True)
        self.transition_quality_input.setText(str(self.transition_quality_slider.value()))
        self.transition_quality_input.blockSignals(False)

    def update_transition_quality_slider_from_input(self):
        """Update the slider when input field changes"""
        try:
            value = int(self.transition_quality_input.text())
            if 1 <= value <= 40:  # Validate bounds
                # Block signals to prevent infinite loop
                self.transition_quality_slider.blockSignals(True)
                self.transition_quality_slider.setValue(value)
                self.transition_quality_slider.blockSignals(False)
                
                # Update the label and call the main update function
                self.transition_quality_label.setText(f"Quality (k): {value}")
                self.update_transition_quality(value)
            else:
                # If value is out of bounds, revert to slider value
                self.transition_quality_input.blockSignals(True)
                self.transition_quality_input.setText(str(self.transition_quality_slider.value()))
                self.transition_quality_input.blockSignals(False)
        except ValueError:
            pass
        
    
    # ------------------------ Configuration Functions based on type --------------------------------
    def sandwitch_config(self):
        self.config_layout.addWidget(self.create_section_label("Faceplate Configuration"))
        
        # strut has different starting default values
        if type == "Strut":
            bottom_face_thickness, self.bottom_face_thickness_input = self.create_input_field("Bottom Face Thickness (mm):", "2", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            top_face_thickness, self.top_face_thickness_input = self.create_input_field("Top Face Thickness (mm):", "1", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
        else:
            bottom_face_thickness, self.bottom_face_thickness_input = self.create_input_field("Bottom Face Thickness (mm):", "4", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            top_face_thickness, self.top_face_thickness_input = self.create_input_field("Top Face Thickness (mm):", "3", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
        
        # connect we setters
        self.bottom_face_thickness_input.textChanged.connect(lambda text: self.update_bottom_face_thickness(text))
        self.top_face_thickness_input.textChanged.connect(lambda text: self.update_top_face_thickness(text))
        
        self.config_layout.addLayout(bottom_face_thickness)
        self.config_layout.addLayout(top_face_thickness)
    
    def beam_plate_config(self):
        self.config_layout.addWidget(self.create_section_label("Thickness Configuration"))
        
        left_support_thickness, self.left_support_thickness_input = self.create_input_field("Left Support Thickness (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
        right_support_thickness, self.right_support_thickness_input = self.create_input_field("Right Support Thickness (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
        
        # connect we setters
        self.left_support_thickness_input.textChanged.connect(lambda text: self.update_left_support_thickness(text))
        self.right_support_thickness_input.textChanged.connect(lambda text: self.update_right_support_thickness(text))
        
        self.config_layout.addLayout(left_support_thickness)
        self.config_layout.addLayout(right_support_thickness)
    
    def setup_tpms_config(self):
        # dimensions section
        self.config_layout.addWidget(self.create_section_label("Dimensions"))
        if form_shape == "Cylindrical":
            
            # radius & connect with update function
            radius_layout, self.radius_input = self.create_input_field("Radius (R, mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.radius_input.textChanged.connect(lambda text: self.update_radius(text))
            
            # height & connect with update function
            height_layout, self.height_input = self.create_input_field("Height (Y, mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))
            
            self.config_layout.addLayout(radius_layout)
            self.config_layout.addLayout(height_layout)
        else:
            
            # length, height, width & connect with update functions
            length_layout, self.length_input = self.create_input_field("Length (X, mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.length_input.textChanged.connect(lambda text: self.update_length(text))
            height_layout, self.height_input = self.create_input_field("Width (Y, mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))    # Note: don't get confused when seeing "width" -> update_height func call
            width_layout, self.width_input = self.create_input_field("Height (Z, mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.width_input.textChanged.connect(lambda text: self.update_width(text))
            
            self.config_layout.addLayout(length_layout)
            self.config_layout.addLayout(height_layout)
            self.config_layout.addLayout(width_layout)
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # repetitions & connect with update functions (not in advanced options)
        self.config_layout.addWidget(self.create_section_label("Repetitions"))
        x_rep_layout, self.x_rep_input = self.create_input_field("X-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.x_rep_input.textChanged.connect(lambda text: self.update_x_repetitions(text))
        y_rep_layout, self.y_rep_input = self.create_input_field("Y-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.y_rep_input.textChanged.connect(lambda text: self.update_y_repetitions(text))
        z_rep_layout, self.z_rep_input = self.create_input_field("Z-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.z_rep_input.textChanged.connect(lambda text: self.update_z_repetitions(text))
        
        self.config_layout.addLayout(x_rep_layout)
        self.config_layout.addLayout(y_rep_layout)
        self.config_layout.addLayout(z_rep_layout)
        self.config_layout.addSpacing(10)
        
        
        # resolution Section
        self.config_layout.addWidget(self.create_section_label("Resolution"))
        points_layout, self.points_input = self.create_input_field("Points:", "50", QRegExpValidator(QRegExp("[0-9]*")))
        self.points_input.textChanged.connect(lambda text: self.update_resolution(text))
        self.config_layout.addLayout(points_layout)
        
        # label for button
        self.config_layout.addSpacing(10)
        advanced_header_layout = QHBoxLayout()
        advanced_header_layout.addWidget(self.create_section_label("Advanced Options"))
        
        # config button
        self.advanced_button = QPushButton("Configure")
        self.advanced_button.setMaximumWidth(100)
        self.advanced_button.clicked.connect(self.popup_advanced_options_adv_tpms)
        
        advanced_header_layout.addWidget(self.advanced_button)
        advanced_header_layout.addStretch()
        self.config_layout.addLayout(advanced_header_layout)
        
        # Initialize advanced options as false
        global advanced_options
        advanced_options = False
        
        if model_ipc_structure == True:
            # 1 plot for reinforcement preview, 1 plot for complementary preview
            self.config_layout.addSpacing(10)
            
            # Track which view is currently shown in main plotter
            self.main_plotter_view = "main"  # Can be "main", "reinforcement", or "complementary"
            
            # Create container for the two additional plotters on the left side
            self.plotters_container = QWidget()
            plotters_layout = QHBoxLayout()
            self.plotters_container.setLayout(plotters_layout)
            
            # Create plotter0 for layer 0
            plot0 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer0_header = QHBoxLayout()
            layer0_label = self.create_section_label("Reinforcement Layer")
            layer0_header.addWidget(layer0_label)
            
            # Arrow button to preview in main plotter
            self.reinf_preview_btn = QPushButton("↗")
            self.reinf_preview_btn.setMaximumWidth(30)
            self.reinf_preview_btn.setMaximumHeight(25)
            self.reinf_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.reinf_preview_btn.setToolTip("Preview in main plotter")
            self.reinf_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("reinforcement"))
            layer0_header.addWidget(self.reinf_preview_btn)
            layer0_header.addStretch()
            
            plot0.addLayout(layer0_header)
            
            self.plotter0 = pvqt.QtInteractor()
            self.plotter0.setMaximumHeight(230)  
            self.plotter0.setMinimumHeight(220)
            plotter0_widget = QWidget()
            plotter0_layout = QVBoxLayout()
            plotter0_layout.addWidget(self.plotter0.interactor)
            plotter0_widget.setLayout(plotter0_layout)
            plot0.addWidget(plotter0_widget)
            
            # Create plotter1 for layer 1
            plot1 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer1_header = QHBoxLayout()
            layer1_label = self.create_section_label("Complementary Layer")
            layer1_header.addWidget(layer1_label)
            
            # Arrow button to preview in main plotter
            self.compl_preview_btn = QPushButton("↗")
            self.compl_preview_btn.setMaximumWidth(30)
            self.compl_preview_btn.setMaximumHeight(25)
            self.compl_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.compl_preview_btn.setToolTip("Preview in main plotter")
            self.compl_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("complementary"))
            layer1_header.addWidget(self.compl_preview_btn)
            layer1_header.addStretch()
            
            plot1.addLayout(layer1_header)
            
            self.plotter1 = pvqt.QtInteractor()
            self.plotter1.setMaximumHeight(230)
            self.plotter1.setMinimumHeight(220)
            plotter1_widget = QWidget()
            plotter1_layout = QVBoxLayout()
            plotter1_layout.addWidget(self.plotter1.interactor)
            plotter1_widget.setLayout(plotter1_layout)
            plot1.addWidget(plotter1_widget)
            
            # Add both plotters to the container horizontally
            plotters_layout.addLayout(plot0)
            plotters_layout.addLayout(plot1)
            
            # Add the plotters container to the left side configuration layout
            self.config_layout.addWidget(self.plotters_container)
            self.plotters_container.setVisible(False)

    def setup_strut_config(self):
        # resolution Section
        self.config_layout.addWidget(self.create_section_label("Resolution"))
        points_layout, self.points_input = self.create_input_field("Points:", "50", QRegExpValidator(QRegExp("[0-9]*")))
        self.points_input.textChanged.connect(lambda text: self.update_resolution(text))
        self.config_layout.addLayout(points_layout)
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # radius Section & connect with update functions
        global design_type
        if design_type == "Strut Diameter":
            self.config_layout.addWidget(self.create_section_label("Radius"))
        elif design_type == "Volume Fraction":
            self.config_layout.addWidget(self.create_section_label("Volume Fraction (%)"))
        
        # no grading, if StrutD we need 4 vals, if VolFrac we need 1 val
        global grading
        if grading == False or grading == "False":
            if design_type == "Strut Diameter":
                bending_layout, self.bending_input = self.create_input_field("Bending Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.bending_input.textChanged.connect(lambda text: self.update_bending_radius(text))
                
                stretching_layout, self.stretching_input = self.create_input_field("Stretching Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.stretching_input.textChanged.connect(lambda text: self.update_stretching_radius(text))
                
                vertical_layout, self.vertical_input = self.create_input_field("Vertical Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.vertical_input.textChanged.connect(lambda text: self.update_vertical_radius(text))
                
                joint_layout, self.joint_input = self.create_input_field("Joint Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.joint_input.textChanged.connect(lambda text: self.update_joint_radius(text))
                
                self.config_layout.addLayout(bending_layout)
                self.config_layout.addLayout(stretching_layout)
                self.config_layout.addLayout(vertical_layout)
                self.config_layout.addLayout(joint_layout)
            elif design_type == "Volume Fraction":
                volume_fraction_for_strut_layout, self.volume_fraction_for_strut_input = self.create_input_field("Volume Fraction (%):", "20", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.volume_fraction_for_strut_input.textChanged.connect(lambda text: self.update_volume_fraction_strut(text))
                self.config_layout.addLayout(volume_fraction_for_strut_layout)
            
        else:  # grading is True, for StrutD we need 8 vals, for VolFrac we need 2 vals
            if design_type == "Strut Diameter":
                bending_layout, self.bending_input = self.create_input_field("Min Bending Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.bending_input.textChanged.connect(lambda text: self.update_bending_radius_min(text))
                stretching_layout, self.stretching_input = self.create_input_field("Min Stretching Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.stretching_input.textChanged.connect(lambda text: self.update_stretching_radius_min(text))
                vertical_layout, self.vertical_input = self.create_input_field("Min Vertical Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.vertical_input.textChanged.connect(lambda text: self.update_vertical_radius_min(text))
                joint_layout, self.joint_input = self.create_input_field("Min Joint Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.joint_input.textChanged.connect(lambda text: self.update_joint_radius_min(text))
                
                bending_layout_max, self.bending_input_max = self.create_input_field("Max Bending Radius (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.bending_input_max.textChanged.connect(lambda text: self.update_bending_radius_max(text))
                stretching_layout_max, self.stretching_input_max = self.create_input_field("Max Stretching Radius (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.stretching_input_max.textChanged.connect(lambda text: self.update_stretching_radius_max(text))
                vertical_layout_max, self.vertical_input_max = self.create_input_field("Max Vertical Radius (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.vertical_input_max.textChanged.connect(lambda text: self.update_vertical_radius_max(text))
                joint_layout_max, self.joint_input_max = self.create_input_field("Max Joint Radius (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.joint_input_max.textChanged.connect(lambda text: self.update_joint_radius_max(text))
                
                self.config_layout.addLayout(bending_layout)
                self.config_layout.addLayout(stretching_layout)
                self.config_layout.addLayout(vertical_layout)
                self.config_layout.addLayout(joint_layout)
                self.config_layout.addLayout(bending_layout_max)
                self.config_layout.addLayout(stretching_layout_max)
                self.config_layout.addLayout(vertical_layout_max)
                self.config_layout.addLayout(joint_layout_max)
            elif design_type == "Volume Fraction":
                min_volume_fraction_for_strut_layout, self.min_volume_fraction_for_strut_input = self.create_input_field("Min Volume Fraction (%):", "20", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.min_volume_fraction_for_strut_input.textChanged.connect(lambda text: self.update_min_volume_fraction_for_strut(text))
                max_volume_fraction_for_strut_layout, self.max_volume_fraction_for_strut_input = self.create_input_field("Max Volume Fraction (%):", "40", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
                self.max_volume_fraction_for_strut_input.textChanged.connect(lambda text: self.update_max_volume_fraction_for_strut(text))
                
                self.config_layout.addLayout(min_volume_fraction_for_strut_layout)
                self.config_layout.addLayout(max_volume_fraction_for_strut_layout)
        
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # repetitions & connect with update functions (not in advanced options)
        self.config_layout.addWidget(self.create_section_label("Repetitions"))
        x_rep_layout, self.x_rep_input = self.create_input_field("X-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.x_rep_input.textChanged.connect(lambda text: self.update_x_repetitions(text))
        y_rep_layout, self.y_rep_input = self.create_input_field("Y-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.y_rep_input.textChanged.connect(lambda text: self.update_y_repetitions(text))
        z_rep_layout, self.z_rep_input = self.create_input_field("Z-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.z_rep_input.textChanged.connect(lambda text: self.update_z_repetitions(text))
        
        self.config_layout.addLayout(x_rep_layout)
        self.config_layout.addLayout(y_rep_layout)
        self.config_layout.addLayout(z_rep_layout)
        self.config_layout.addSpacing(10)
        
        # no advanced options for strut anymore
        global advanced_options
        advanced_options = False
        
        if model_ipc_structure == True:
            # 1 plot for reinforcement preview, 1 plot for complementary preview
            self.config_layout.addSpacing(10)
            
            # Track which view is currently shown in main plotter
            self.main_plotter_view = "main"  # Can be "main", "reinforcement", or "complementary"
            
            # Create container for the two additional plotters on the left side
            self.plotters_container = QWidget()
            plotters_layout = QHBoxLayout()
            self.plotters_container.setLayout(plotters_layout)
            
            # Create plotter0 for layer 0
            plot0 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer0_header = QHBoxLayout()
            layer0_label = self.create_section_label("Reinforcement Layer")
            layer0_header.addWidget(layer0_label)
            
            # Arrow button to preview in main plotter
            self.reinf_preview_btn = QPushButton("↗")
            self.reinf_preview_btn.setMaximumWidth(30)
            self.reinf_preview_btn.setMaximumHeight(25)
            self.reinf_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.reinf_preview_btn.setToolTip("Preview in main plotter")
            self.reinf_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("reinforcement"))
            layer0_header.addWidget(self.reinf_preview_btn)
            layer0_header.addStretch()
            
            plot0.addLayout(layer0_header)
            
            self.plotter0 = pvqt.QtInteractor()
            self.plotter0.setMaximumHeight(230)  
            self.plotter0.setMinimumHeight(220)
            plotter0_widget = QWidget()
            plotter0_layout = QVBoxLayout()
            plotter0_layout.addWidget(self.plotter0.interactor)
            plotter0_widget.setLayout(plotter0_layout)
            plot0.addWidget(plotter0_widget)
            
            # Create plotter1 for layer 1
            plot1 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer1_header = QHBoxLayout()
            layer1_label = self.create_section_label("Complementary Layer")
            layer1_header.addWidget(layer1_label)
            
            # Arrow button to preview in main plotter
            self.compl_preview_btn = QPushButton("↗")
            self.compl_preview_btn.setMaximumWidth(30)
            self.compl_preview_btn.setMaximumHeight(25)
            self.compl_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.compl_preview_btn.setToolTip("Preview in main plotter")
            self.compl_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("complementary"))
            layer1_header.addWidget(self.compl_preview_btn)
            layer1_header.addStretch()
            
            plot1.addLayout(layer1_header)
            
            self.plotter1 = pvqt.QtInteractor()
            self.plotter1.setMaximumHeight(230)
            self.plotter1.setMinimumHeight(220)
            plotter1_widget = QWidget()
            plotter1_layout = QVBoxLayout()
            plotter1_layout.addWidget(self.plotter1.interactor)
            plotter1_widget.setLayout(plotter1_layout)
            plot1.addWidget(plotter1_widget)
            
            # Add both plotters to the container horizontally
            plotters_layout.addLayout(plot0)
            plotters_layout.addLayout(plot1)
            
            # Add the plotters container to the left side configuration layout
            self.config_layout.addWidget(self.plotters_container)
            self.plotters_container.setVisible(False)
    
    def setup_spinodal_config(self):
        # dimensions section & connect with update functions
        self.config_layout.addWidget(self.create_section_label("Dimensions"))
        if form_shape == "Cylindrical":
            radius_layout, self.radius_input = self.create_input_field("Radius (mm):", "5", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.radius_input.textChanged.connect(lambda text: self.update_radius(text))
            height_layout, self.height_input = self.create_input_field("Height (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))
            
            self.config_layout.addLayout(radius_layout)
            self.config_layout.addLayout(height_layout)
        else:
            length_layout, self.length_input = self.create_input_field("Length (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.length_input.textChanged.connect(lambda text: self.update_length(text))
            height_layout, self.height_input = self.create_input_field("Height (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.height_input.textChanged.connect(lambda text: self.update_height(text))
            width_layout, self.width_input = self.create_input_field("Width (mm):", "10", QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*")))
            self.width_input.textChanged.connect(lambda text: self.update_width(text))
            
            self.config_layout.addLayout(length_layout)
            self.config_layout.addLayout(height_layout)
            self.config_layout.addLayout(width_layout)
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # repetitions & connect with update functions (not in advanced options)
        self.config_layout.addWidget(self.create_section_label("Repetitions"))
        x_rep_layout, self.x_rep_input = self.create_input_field("X-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.x_rep_input.textChanged.connect(lambda text: self.update_x_repetitions(text))
        y_rep_layout, self.y_rep_input = self.create_input_field("Y-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.y_rep_input.textChanged.connect(lambda text: self.update_y_repetitions(text))
        z_rep_layout, self.z_rep_input = self.create_input_field("Z-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
        self.z_rep_input.textChanged.connect(lambda text: self.update_z_repetitions(text))
        
        self.config_layout.addLayout(x_rep_layout)
        self.config_layout.addLayout(y_rep_layout)
        self.config_layout.addLayout(z_rep_layout)
        self.config_layout.addSpacing(10)
        
        # resolution Section & connect with update functions
        self.config_layout.addWidget(self.create_section_label("Resolution"))
        points_layout, self.points_input = self.create_input_field("Points:", "121", QRegExpValidator(QRegExp("[0-9]*")))
        self.points_input.textChanged.connect(lambda text: self.update_resolution(text))
        self.config_layout.addLayout(points_layout)
        
        # set initial state of advanced options (we don't have advanced options for spinodal)
        global advanced_options
        advanced_options = False
        
        if model_ipc_structure == True:
            # 1 plot for reinforcement preview, 1 plot for complementary preview
            self.config_layout.addSpacing(10)
            
            # Track which view is currently shown in main plotter
            self.main_plotter_view = "main"  # Can be "main", "reinforcement", or "complementary"
            
            # Create container for the two additional plotters on the left side
            self.plotters_container = QWidget()
            plotters_layout = QHBoxLayout()
            self.plotters_container.setLayout(plotters_layout)
            
            # Create plotter0 for layer 0
            plot0 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer0_header = QHBoxLayout()
            layer0_label = self.create_section_label("Reinforcement Layer")
            layer0_header.addWidget(layer0_label)
            
            # Arrow button to preview in main plotter
            self.reinf_preview_btn = QPushButton("↗")
            self.reinf_preview_btn.setMaximumWidth(30)
            self.reinf_preview_btn.setMaximumHeight(25)
            self.reinf_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.reinf_preview_btn.setToolTip("Preview in main plotter")
            self.reinf_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("reinforcement"))
            layer0_header.addWidget(self.reinf_preview_btn)
            layer0_header.addStretch()
            
            plot0.addLayout(layer0_header)
            
            self.plotter0 = pvqt.QtInteractor()
            self.plotter0.setMaximumHeight(230)  
            self.plotter0.setMinimumHeight(220)
            plotter0_widget = QWidget()
            plotter0_layout = QVBoxLayout()
            plotter0_layout.addWidget(self.plotter0.interactor)
            plotter0_widget.setLayout(plotter0_layout)
            plot0.addWidget(plotter0_widget)
            
            # Create plotter1 for layer 1
            plot1 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer1_header = QHBoxLayout()
            layer1_label = self.create_section_label("Complementary Layer")
            layer1_header.addWidget(layer1_label)
            
            # Arrow button to preview in main plotter
            self.compl_preview_btn = QPushButton("↗")
            self.compl_preview_btn.setMaximumWidth(30)
            self.compl_preview_btn.setMaximumHeight(25)
            self.compl_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.compl_preview_btn.setToolTip("Preview in main plotter")
            self.compl_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("complementary"))
            layer1_header.addWidget(self.compl_preview_btn)
            layer1_header.addStretch()
            
            plot1.addLayout(layer1_header)
            
            self.plotter1 = pvqt.QtInteractor()
            self.plotter1.setMaximumHeight(230)
            self.plotter1.setMinimumHeight(220)
            plotter1_widget = QWidget()
            plotter1_layout = QVBoxLayout()
            plotter1_layout.addWidget(self.plotter1.interactor)
            plotter1_widget.setLayout(plotter1_layout)
            plot1.addWidget(plotter1_widget)
            
            # Add both plotters to the container horizontally
            plotters_layout.addLayout(plot0)
            plotters_layout.addLayout(plot1)
            
            # Add the plotters container to the left side configuration layout
            self.config_layout.addWidget(self.plotters_container)
            self.plotters_container.setVisible(False)
    
    def setup_hybrid_config(self):
        # topological Selection Section
        #self.config_layout.addWidget(self.create_section_label("Topological Selection"))
        
        # order of specific topologies that are shown in the combo boxes (2 cases for now)
        specific_topo_defaults = ["Gyroid (GY)", "FischerkochS (FKS)", "Pcell (SPC)"]  
        
        # based on the number of layers, initialize the global variables with their default values
        global general_topo_for_hybrid, specific_topo_for_hybrid, x_rep_for_hybrid, y_rep_for_hybrid, z_rep_for_hybrid, vol_fraction_for_hybrid
        general_topo_for_hybrid = ["Sheet"] * int(number_of_layers)
        specific_topo_for_hybrid = [
            specific_topo_defaults[i % len(specific_topo_defaults)]
            for i in range(int(number_of_layers))
        ]
        x_rep_for_hybrid = ["2"] * int(number_of_layers)                        
        y_rep_for_hybrid = ["2"] * int(number_of_layers)
        z_rep_for_hybrid = ["2"] * int(number_of_layers)
        vol_fraction_for_hybrid = ["30"] * int(number_of_layers)
        
    
        # Layer Selection based on the number of layers (topological, reps, vol fraction)
        for i in range(1, int(number_of_layers) + 1):
            layer_group = QGroupBox(f"Layer {i}")
            layer_layout = QVBoxLayout()
            
            # general Topology for Layer
            general_topo_layout = QHBoxLayout()
            general_topo_label = QLabel("General Topology:")
            general_topo_combo = QComboBox()
            general_topo_combo.addItems(["Sheet", "Skeletal", "Exoskeletal"])
            
            # Store the combo box in an instance variable or list to maintain reference
            # Store the layer index for the update function because otherwise only the last layer index will be used
            general_topo_combo.setProperty('layer_index', i)  
            general_topo_combo.currentIndexChanged.connect(
                lambda _, combo=general_topo_combo: self.update_general_topo_for_hybrid(combo.currentText(), combo.property('layer_index'))
            )
            
            general_topo_layout.addWidget(general_topo_label)
            general_topo_layout.addWidget(general_topo_combo)
            
            # specific Topology for Layer
            specific_topo_layout = QHBoxLayout()
            specific_topo_label = QLabel("Specific Topology:")
            specific_topo_combo = QComboBox()
            items = [
                "Gyroid (GY)", "IWP", "Pcell (SPC)", "FischerkochS (FKS)", 
                "Schwarz_Diamond (SCD)", "Neovius (NE)", "Lidinoid (LD)",
                "SchwarzHexagonal (SCH)", "SplitP (SLP)", "I2Y (I2Y)", 
                "Fisher–KochC(S) (FKCS)", "F-RD (FRD)"
            ]
            specific_topo_combo.addItems(items)
            default_value = specific_topo_for_hybrid[i-1]
            specific_topo_combo.setCurrentIndex(items.index(default_value))

            # Store the combo box in an instance variable or list to maintain reference
            # Store the layer index for the update function because otherwise only the last layer index will be used
            specific_topo_combo.setProperty('layer_index', i)  
            specific_topo_combo.currentIndexChanged.connect(
                lambda _, combo=specific_topo_combo: self.update_specific_topo_for_hybrid(combo.currentText(), combo.property('layer_index'))
            )
            
            specific_topo_layout.addWidget(specific_topo_label)
            specific_topo_layout.addWidget(specific_topo_combo)
            
            
            # repetitions x,y,z for Layer
            rep_layout = QHBoxLayout()
            x_rep_layout, x_rep_input = self.create_input_field("X-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
            x_rep_input.setProperty('layer_index', i)
            x_rep_input.textChanged.connect(
                lambda text, input=x_rep_input: self.update_x_rep_for_hybrid(text, input.property('layer_index'))
            )
            
            y_rep_layout, y_rep_input = self.create_input_field("Y-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
            y_rep_input.setProperty('layer_index', i)
            y_rep_input.textChanged.connect(
                lambda text, input=y_rep_input: self.update_y_rep_for_hybrid(text, input.property('layer_index'))
            )
            
            z_rep_layout, z_rep_input = self.create_input_field("Z-axis:", "2", QRegExpValidator(QRegExp("[0-9]*")))
            z_rep_input.setProperty('layer_index', i)
            z_rep_input.textChanged.connect(
                lambda text, input=z_rep_input: self.update_z_rep_for_hybrid(text, input.property('layer_index'))
            )
            
            rep_layout.addLayout(x_rep_layout)
            rep_layout.addLayout(y_rep_layout)
            rep_layout.addLayout(z_rep_layout)
            
            # volume fraction for Layer
            volume_layout = QFormLayout()
            volume_fraction_slider = QSlider(Qt.Horizontal)  
            volume_fraction_slider.setMinimum(10)
            volume_fraction_slider.setMaximum(60)
            volume_fraction_label = QLabel(f"Volume Fraction of Layer {i}: 30%")  
            volume_fraction_slider.setValue(30)

            # store the layer index
            volume_fraction_slider.setProperty('layer_index', i)
            volume_fraction_label.setProperty('layer_index', i)

            # update the connection to use local variables
            volume_fraction_slider.valueChanged.connect(
                lambda value, slider=volume_fraction_slider, label=volume_fraction_label: [
                    label.setText(f"Volume Fraction of Layer {slider.property('layer_index')}: {value}%"),
                    self.update_vol_fraction_for_hybrid(value, slider.property('layer_index'))
                ]
            )
            volume_layout.addRow(volume_fraction_label, volume_fraction_slider)
            
            
            
            layer_layout.addLayout(general_topo_layout)
            layer_layout.addLayout(specific_topo_layout)
            layer_layout.addLayout(rep_layout)
            layer_layout.addLayout(volume_layout)
            layer_group.setLayout(layer_layout)
            self.config_layout.addWidget(layer_group)
        
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # resolution Section & connect with update functions
        self.config_layout.addWidget(self.create_section_label("Resolution"))
        points_layout, self.points_input = self.create_input_field("Points:", "50", QRegExpValidator(QRegExp("[0-9]*")))
        self.points_input.textChanged.connect(lambda text: self.update_resolution(text))
        self.config_layout.addLayout(points_layout)
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # Cylindrical Hybrid Type - CylHyType in backend (for Cylindrical only)
        global grading_for_hybrid
        if grading_for_hybrid == "Cylindrical":
            self.config_layout.addWidget(self.create_section_label("Cylindrical Hybrid Type"))
            cylhy_layout = QHBoxLayout()
            self.cylhy_combo = QComboBox()
            self.cylhy_combo.addItems(["Single-axis cylindrical (y)", "Double-axis cylindrical (y, x)", "Triple-axis cylindrical (y, x, z)"])
            self.cylhy_combo.currentIndexChanged.connect(lambda index: self.update_cylindrical_hybrid_type(index))
            self.cylhy_combo.setCurrentIndex(2)  # default to Triple-axis cylindrical
            cylhy_layout.addWidget(self.cylhy_combo)
            cylhy_layout.addStretch()
            self.config_layout.addLayout(cylhy_layout)
            
            # add some spacing between the sections
            self.config_layout.addSpacing(10)
            
        # Transition Section
        transition_header_layout = QHBoxLayout()
        transition_header_layout.addWidget(self.create_section_label("Transition Options"))
        
        # label for button
        self.transition_conf_button = QPushButton("Configure")
        self.transition_conf_button.setMaximumWidth(100)
        self.transition_conf_button.clicked.connect(self.popup_transition_options_hybrid)
        transition_header_layout.addWidget(self.transition_conf_button)
        transition_header_layout.addStretch()
        self.config_layout.addLayout(transition_header_layout)
        
        # Hybrid does not have advanced options anymore (reps are in each layer section)
        global advanced_options
        advanced_options = False
        
        if model_ipc_structure == True:
            # 1 plot for reinforcement preview, 1 plot for complementary preview
            self.config_layout.addSpacing(10)
            
            # Track which view is currently shown in main plotter
            self.main_plotter_view = "main"  # Can be "main", "reinforcement", or "complementary"
            
            # Create container for the two additional plotters on the left side
            self.plotters_container = QWidget()
            plotters_layout = QHBoxLayout()
            self.plotters_container.setLayout(plotters_layout)
            
            # Create plotter0 for layer 0
            plot0 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer0_header = QHBoxLayout()
            layer0_label = self.create_section_label("Reinforcement Layer")
            layer0_header.addWidget(layer0_label)
            
            # Arrow button to preview in main plotter
            self.reinf_preview_btn = QPushButton("↗")
            self.reinf_preview_btn.setMaximumWidth(30)
            self.reinf_preview_btn.setMaximumHeight(25)
            self.reinf_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.reinf_preview_btn.setToolTip("Preview in main plotter")
            self.reinf_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("reinforcement"))
            layer0_header.addWidget(self.reinf_preview_btn)
            layer0_header.addStretch()
            
            plot0.addLayout(layer0_header)
            
            self.plotter0 = pvqt.QtInteractor()
            self.plotter0.setMaximumHeight(230)  
            self.plotter0.setMinimumHeight(220)
            plotter0_widget = QWidget()
            plotter0_layout = QVBoxLayout()
            plotter0_layout.addWidget(self.plotter0.interactor)
            plotter0_widget.setLayout(plotter0_layout)
            plot0.addWidget(plotter0_widget)
            
            # Create plotter1 for layer 1
            plot1 = QVBoxLayout()
            
            # Header with label and preview arrow button
            layer1_header = QHBoxLayout()
            layer1_label = self.create_section_label("Complementary Layer")
            layer1_header.addWidget(layer1_label)
            
            # Arrow button to preview in main plotter
            self.compl_preview_btn = QPushButton("↗")
            self.compl_preview_btn.setMaximumWidth(30)
            self.compl_preview_btn.setMaximumHeight(25)
            self.compl_preview_btn.setStyleSheet("""
                QPushButton {
                    background-color: #4682b4;
                    color: white;
                    border-radius: 4px;
                    font-size: 16px;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #5a9bd4;
                }
                QPushButton:pressed {
                    background-color: #3a6a94;
                }
            """)
            self.compl_preview_btn.setToolTip("Preview in main plotter")
            self.compl_preview_btn.clicked.connect(lambda: self.preview_in_main_plotter("complementary"))
            layer1_header.addWidget(self.compl_preview_btn)
            layer1_header.addStretch()
            
            plot1.addLayout(layer1_header)
            
            self.plotter1 = pvqt.QtInteractor()
            self.plotter1.setMaximumHeight(230)
            self.plotter1.setMinimumHeight(220)
            plotter1_widget = QWidget()
            plotter1_layout = QVBoxLayout()
            plotter1_layout.addWidget(self.plotter1.interactor)
            plotter1_widget.setLayout(plotter1_layout)
            plot1.addWidget(plotter1_widget)
            
            # Add both plotters to the container horizontally
            plotters_layout.addLayout(plot0)
            plotters_layout.addLayout(plot1)
            
            # Add the plotters container to the left side configuration layout
            self.config_layout.addWidget(self.plotters_container)
            self.plotters_container.setVisible(False)
        
    def setup_layered_config(self):
        
        # repetitions & connect with update functions (not in advanced options)
        self.config_layout.addWidget(self.create_section_label("Repetitions"))
        x_rep_layout, self.x_rep_input = self.create_input_field("X-axis:", "3", QRegExpValidator(QRegExp("[0-9]*")))
        self.x_rep_input.textChanged.connect(lambda text: self.update_x_repetitions(text))
        y_rep_layout, self.y_rep_input = self.create_input_field("Y-axis:", "3", QRegExpValidator(QRegExp("[0-9]*")))
        self.y_rep_input.textChanged.connect(lambda text: self.update_y_repetitions(text))
        z_rep_layout, self.z_rep_input = self.create_input_field("Z-axis:", "3", QRegExpValidator(QRegExp("[0-9]*")))
        self.z_rep_input.textChanged.connect(lambda text: self.update_z_repetitions(text))
        
        self.config_layout.addLayout(x_rep_layout)
        self.config_layout.addLayout(y_rep_layout)
        self.config_layout.addLayout(z_rep_layout)
        self.config_layout.addSpacing(10)
        
        # resolution Section & connect with update functions
        self.config_layout.addWidget(self.create_section_label("Resolution"))
        points_layout, self.points_input = self.create_input_field("Points:", "61", QRegExpValidator(QRegExp("[0-9]*")))
        self.points_input.textChanged.connect(lambda text: self.update_resolution(text))
        self.config_layout.addLayout(points_layout)
        
        # add some spacing between the sections
        self.config_layout.addSpacing(10)
        
        # based on the number of layers, initialize the global variables
        # they need to sum up to 100
        # for naming: layer densities == layer volume fractions
        global layer_densities
        portion = 100 // int(number_of_layers)
        layer_densities = [portion] * int(number_of_layers)
        
        # density Section for each layer 
        # if the number of layers is 2, the default values are 70 and 30
        # otherwise, the default values are equal portions
        self.config_layout.addWidget(self.create_section_label("Layer Volume Fractions"))
        if number_of_layers == 2:
            layer_densities = [70, 30]
            
            layer_layout_1, layer_input_1 = self.create_input_field("Reinforcement Phase VF (%):", "70", QRegExpValidator(QRegExp("[0-9]*")))
            layer_input_1.setProperty('layer_index', 1)
            layer_input_1.textChanged.connect(
                lambda text, input=layer_input_1: self.update_layer_density(text, input.property('layer_index'))
            )
            self.config_layout.addLayout(layer_layout_1)
            
            layer_layout_2, layer_input_2 = self.create_input_field("Complementary Phase VF (%):", "30", QRegExpValidator(QRegExp("[0-9]*")))
            layer_input_2.setProperty('layer_index', 2)
            layer_input_2.textChanged.connect(
                lambda text, input=layer_input_2: self.update_layer_density(text, input.property('layer_index'))
            )
            self.config_layout.addLayout(layer_layout_2)
            
        else:
            portion = 100 // int(number_of_layers)
            layer_densities = [portion] * int(number_of_layers)
        
        
            for i in range(1, int(number_of_layers) + 1):
                layer_layout, layer_input = self.create_input_field(f"Layer {i} Volume Fraction (%):", str(layer_densities[i-1]), QRegExpValidator(QRegExp("[0-9]*")))
                
                # Store the layer index as a property of the input field
                layer_input.setProperty('layer_index', i)
                layer_input.textChanged.connect(
                lambda text, input=layer_input: self.update_layer_density(text, input.property('layer_index'))
                )
                self.config_layout.addLayout(layer_layout)
            
            
        # 1 plot for reinforcement preview, 1 plot for complementary preview
        self.config_layout.addSpacing(30)
        
        # Create container for the two additional plotters on the left side
        plotters_container = QWidget()
        plotters_layout = QHBoxLayout()
        plotters_container.setLayout(plotters_layout)
        
        # Create plotter0 for layer 0
        plot0 = QVBoxLayout()
        layer0_label = self.create_section_label("Reinforcement Layer")
        plot0.addWidget(layer0_label)
        
        self.plotter0 = pvqt.QtInteractor()
        self.plotter0.setMaximumHeight(230)  
        self.plotter0.setMinimumHeight(220)
        plotter0_widget = QWidget()
        plotter0_layout = QVBoxLayout()
        plotter0_layout.addWidget(self.plotter0.interactor)
        plotter0_widget.setLayout(plotter0_layout)
        plot0.addWidget(plotter0_widget)
        
        # Create plotter1 for layer 1
        plot1 = QVBoxLayout()
        layer1_label = self.create_section_label("Complementary Layer")
        plot1.addWidget(layer1_label)
        self.plotter1 = pvqt.QtInteractor()
        self.plotter1.setMaximumHeight(230)
        self.plotter1.setMinimumHeight(220)
        plotter1_widget = QWidget()
        plotter1_layout = QVBoxLayout()
        plotter1_layout.addWidget(self.plotter1.interactor)
        plotter1_widget.setLayout(plotter1_layout)
        plot1.addWidget(plotter1_widget)
        
        # Add both plotters to the container horizontally
        plotters_layout.addLayout(plot0)
        plotters_layout.addLayout(plot1)
        
        # Add the plotters container to the left side configuration layout
        self.config_layout.addWidget(plotters_container)
        self.load_individual_layers()


    # ------------- Helper Function - ADVANCED TPMS --------------------------------
    def popup_advanced_options_adv_tpms(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Advanced Options")
        dialog.setMinimumWidth(300)
        dialog_layout = QVBoxLayout(dialog)

        # Store input widgets for validation
        self.stretch_inputs = {}
        self.rotation_inputs = {}
        
        # Store temporary values (will only update globals on OK)
        self.temp_values = {}

        # Stretching Group
        stretch_group = QGroupBox("Stretching")
        stretch_layout = QVBoxLayout()
        
        global x_stretching, y_stretching, z_stretching
        
        # Initialize temporary values from current globals
        self.temp_values['x_stretching'] = x_stretching
        self.temp_values['y_stretching'] = y_stretching
        self.temp_values['z_stretching'] = z_stretching
        
        # X-Stretch
        x_stretch_layout, x_stretch_input = self.create_input_field(
            "X-axis:", str(x_stretching), QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*"))
        )
        self.stretch_inputs['x'] = x_stretch_input
        x_stretch_input.textChanged.connect(lambda text: self.update_temp_value('x_stretching', text))
        x_stretch_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        stretch_layout.addLayout(x_stretch_layout)

        # Y-Stretch
        y_stretch_layout, y_stretch_input = self.create_input_field(
            "Y-axis:", str(y_stretching), QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*"))
        )
        self.stretch_inputs['y'] = y_stretch_input
        y_stretch_input.textChanged.connect(lambda text: self.update_temp_value('y_stretching', text))
        y_stretch_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        stretch_layout.addLayout(y_stretch_layout)
        
        # Z-Stretch
        z_stretch_layout, z_stretch_input = self.create_input_field(
            "Z-axis:", str(z_stretching), QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*"))
        )
        self.stretch_inputs['z'] = z_stretch_input
        z_stretch_input.textChanged.connect(lambda text: self.update_temp_value('z_stretching', text))
        z_stretch_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        stretch_layout.addLayout(z_stretch_layout)
        
        stretch_group.setLayout(stretch_layout)
        dialog_layout.addWidget(stretch_group)

        # Rotation Group
        rot_group = QGroupBox("Rotation")
        rot_layout = QVBoxLayout()
        
        global x_rotation, y_rotation, z_rotation
        
        # Initialize temporary values from current globals
        self.temp_values['x_rotation'] = x_rotation
        self.temp_values['y_rotation'] = y_rotation
        self.temp_values['z_rotation'] = z_rotation
        
        # X-Rotation
        x_rot_layout, x_rot_input = self.create_input_field(
            "X-axis (°):", str(x_rotation), QDoubleValidator()
        )
        self.rotation_inputs['x'] = x_rot_input
        x_rot_input.textChanged.connect(lambda text: self.update_temp_value('x_rotation', text))
        x_rot_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        rot_layout.addLayout(x_rot_layout)
        
        # Y-Rotation
        y_rot_layout, y_rot_input = self.create_input_field(
            "Y-axis (°):", str(y_rotation), QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*"))
        )
        self.rotation_inputs['y'] = y_rot_input
        y_rot_input.textChanged.connect(lambda text: self.update_temp_value('y_rotation', text))
        y_rot_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        rot_layout.addLayout(y_rot_layout)
        
        # Z-Rotation
        z_rot_layout, z_rot_input = self.create_input_field(
            "Z-axis (°):", str(z_rotation), QRegExpValidator(QRegExp(r"[0-9]*\.?[0-9]*"))
        )
        self.rotation_inputs['z'] = z_rot_input
        z_rot_input.textChanged.connect(lambda text: self.update_temp_value('z_rotation', text))
        z_rot_input.textChanged.connect(self.validate_dialog_inputs_adv_tpms)
        rot_layout.addLayout(z_rot_layout)
        
        rot_group.setLayout(rot_layout)
        dialog_layout.addWidget(rot_group)
        
        # Button box with Reset option
        button_box = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel | QDialogButtonBox.Reset
        )
        
        # Store button reference for enabling/disabling
        self.dialog_ok_button = button_box.button(QDialogButtonBox.Ok)
        
        button_box.accepted.connect(lambda: self.accept_advanced_options_adv_tpms(dialog))
        button_box.rejected.connect(dialog.reject)
        button_box.button(QDialogButtonBox.Reset).clicked.connect(self.reset_advanced_options_adv_tpms)
        
        dialog_layout.addWidget(button_box)
        
        # Validate initial state
        self.validate_dialog_inputs_adv_tpms()
        
        # Show the Dialog
        dialog.exec()

    def update_temp_value(self, key, text):
        # Update temporary value without affecting globals
        if text.strip():
            try:
                self.temp_values[key] = float(text)
            except ValueError:
                pass

    def accept_advanced_options_adv_tpms(self, dialog):
        # Only commit temporary values to globals when OK is clicked
        global x_stretching, y_stretching, z_stretching
        global x_rotation, y_rotation, z_rotation
        global advanced_options
        
        # Commit temporary values to globals
        x_stretching = self.temp_values.get('x_stretching', 1)
        y_stretching = self.temp_values.get('y_stretching', 1)
        z_stretching = self.temp_values.get('z_stretching', 1)
        x_rotation = self.temp_values.get('x_rotation', 0)
        y_rotation = self.temp_values.get('y_rotation', 0)
        z_rotation = self.temp_values.get('z_rotation', 0)
        
        # Update your actual global update functions if needed
        self.update_x_stretch(str(x_stretching))
        self.update_y_stretch(str(y_stretching))
        self.update_z_stretch(str(z_stretching))
        self.update_x_rotation(str(x_rotation))
        self.update_y_rotation(str(y_rotation))
        self.update_z_rotation(str(z_rotation))
        
        # Check if any values differ from defaults
        has_advanced = (x_stretching != 1 or y_stretching != 1 or z_stretching != 1 or
                        x_rotation != 0 or y_rotation != 0 or z_rotation != 0)
        
        advanced_options = has_advanced
        
        dialog.accept()

    def reset_advanced_options_adv_tpms(self):
        # Reset all fields to defaults in the dialog
        # Reset temporary values
        self.temp_values['x_stretching'] = 1
        self.temp_values['y_stretching'] = 1
        self.temp_values['z_stretching'] = 1
        self.temp_values['x_rotation'] = 0
        self.temp_values['y_rotation'] = 0
        self.temp_values['z_rotation'] = 0
        
        # Update all input fields
        for key, input_widget in self.stretch_inputs.items():
            input_widget.setText("1")
        
        for key, input_widget in self.rotation_inputs.items():
            input_widget.setText("0")

    def validate_dialog_inputs_adv_tpms(self):
        # Check if all inputs are valid and enable/disable OK button
        all_valid = True
        
        # Check all stretch inputs
        for input_widget in self.stretch_inputs.values():
            if input_widget.text().strip() == "":
                all_valid = False
                break
        
        # Check all rotation inputs
        if all_valid:
            for input_widget in self.rotation_inputs.values():
                if input_widget.text().strip() == "":
                    all_valid = False
                    break
        
        # Enable/disable OK button
        if hasattr(self, 'dialog_ok_button'):
            self.dialog_ok_button.setEnabled(all_valid)
    
    # ------------- Helper Function - Transition HYBRID --------------------------------
    def popup_transition_options_hybrid(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Transition Options")
        dialog.setMinimumWidth(300)
        dialog_layout = QVBoxLayout(dialog)

        global transition_location, transition_quality
        transition_location = float(transition_location)
        transition_quality = int(transition_quality)
        
        # Transition Location
        transition_layout = QVBoxLayout()
        self.transition_slider = QSlider(Qt.Horizontal)
        self.transition_slider.setMinimum(20)
        self.transition_slider.setMaximum(80)
        self.transition_slider.setValue(int(transition_location * 100))
        self.transition_slider.setTickInterval(1)
        self.transition_slider.setTickPosition(QSlider.TicksBelow)
        self.transition_label = QLabel(f"Location: {transition_location:.2f}")
        
        # connect the slider with the update function & also update the labels
        self.transition_slider.valueChanged.connect(
            lambda: [self.transition_label.setText(f"Location: {self.transition_slider.value() / 100:.2f}"),
                    self.update_transition_location(self.transition_slider.value() / 100)]
        )
        
        transition_layout.addWidget(self.transition_label)
        transition_layout.addWidget(self.transition_slider)
        dialog_layout.addLayout(transition_layout)
        
        # Transition Quality
        dialog_layout.addWidget(self.create_section_label("Transition Quality"))
        transition_quality_layout = QVBoxLayout()
        self.transition_quality_slider = QSlider(Qt.Horizontal)
        self.transition_quality_slider.setMinimum(1)            # Actual value range
        self.transition_quality_slider.setMaximum(40)
        self.transition_quality_slider.setValue(int(transition_quality))
        self.transition_quality_slider.setTickInterval(5)
        self.transition_quality_slider.setTickPosition(QSlider.TicksBelow)
        self.transition_quality_label = QLabel(f"Quality (k): {int(transition_quality)}")
        label_container = QWidget()
        label_layout = QHBoxLayout()
        label_layout.addWidget(QLabel("High"))
        label_layout.addStretch()
        label_layout.addWidget(QLabel("Medium"))
        label_layout.addStretch()
        label_layout.addWidget(QLabel("Low"))
        label_container.setLayout(label_layout)
        transition_quality_layout.addWidget(self.transition_quality_label)
        transition_quality_layout.addWidget(self.transition_quality_slider)
        transition_quality_layout.addWidget(label_container)    # Add Low/Medium/High labels
        dialog_layout.addLayout(transition_quality_layout)
        
        # connect slider changes - updates both label and calls update function
        self.transition_quality_slider.valueChanged.connect(
            lambda: [
                self.transition_quality_label.setText(f"Quality (k): {self.transition_quality_slider.value()}"),
                self.update_transition_quality(self.transition_quality_slider.value())
            ]
        )
        
        

        # Button box
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(dialog.accept)
        button_box.rejected.connect(dialog.reject)
        dialog_layout.addWidget(button_box)

        # Show the Dialog
        dialog.exec()
    
    
    # -------- Helper Function - Empty Vals, Generate, ranges TPMS, Swap main plot with subplots  ----------------
    def preview_in_main_plotter(self, layer_type):

        # If clicking the same layer that's already showing, switch back to main
        if self.main_plotter_view == layer_type:
            layer_type = "main"
        
        # Clear the main plotter
        self.plotter.clear()
        
        # Get current color from global variable
        global model_color
        
        # Update the main plotter label to show current view
        if layer_type == "reinforcement":
            if self.current_mesh_reinf is not None:
                self.plotter.add_mesh(self.current_mesh_compl, color=model_color, opacity=1, show_edges=True)
                self.main_plotter_label.setText("Main 3D Model Visualization - Reinforcement Layer Preview")
                self.main_plotter_label.setStyleSheet("color: #ff6b6b; font-size: 16px; font-weight: bold;")
                
                # Update button appearance
                self.reinf_preview_btn.setText("↙")
                self.reinf_preview_btn.setToolTip("Return to main view")
                self.color_selector.setDisabled(True)  # Disable color selector in layer preview mode
                if hasattr(self, 'compl_preview_btn'):
                    self.compl_preview_btn.setText("↗")
                    self.compl_preview_btn.setToolTip("Preview in main plotter")
                    
            
        elif layer_type == "complementary":
            if self.current_mesh_compl is not None:
                # Use darker color for complementary layer (consistent with ColorSelectorWidget)
                darker_color = (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
                self.plotter.add_mesh(self.current_mesh_reinf, color=darker_color, opacity=1, show_edges=True)
                self.main_plotter_label.setText("Main 3D Model Visualization - Complementary Layer Preview")
                self.main_plotter_label.setStyleSheet("color: #4ecdc4; font-size: 16px; font-weight: bold;")
                
                # Update button appearance
                self.compl_preview_btn.setText("↙")
                self.compl_preview_btn.setToolTip("Return to main view")
                self.color_selector.setDisabled(True)  # Disable color selector in layer preview mode
                if hasattr(self, 'reinf_preview_btn'):
                    self.reinf_preview_btn.setText("↗")
                    self.reinf_preview_btn.setToolTip("Preview in main plotter")
            
        else:  # main view
            if self.current_mesh is not None:
                # Handle both single mesh and list of meshes
                if isinstance(self.current_mesh, list):
                    # Multiple layers - use color variations
                    colors = [
                        model_color,
                        (model_color[0] * 0.7, model_color[1] * 0.7, model_color[2] * 0.7)
                    ]
                    for i, mesh in enumerate(self.current_mesh):
                        layer_color = colors[i] if i < len(colors) else model_color
                        self.plotter.add_mesh(mesh, color=layer_color, opacity=1, show_edges=True, name=f'layer_{i}')
                else:
                    # Single mesh
                    self.plotter.add_mesh(self.current_mesh, color=model_color, opacity=1, show_edges=True)
                self.main_plotter_label.setText("Main 3D Model Visualization")
                self.main_plotter_label.setStyleSheet("color: #1e90ff; font-size: 16px;")
                
                # Reset all button appearances
                if hasattr(self, 'reinf_preview_btn'):
                    self.reinf_preview_btn.setText("↗")
                    self.reinf_preview_btn.setToolTip("Preview in main plotter")
                if hasattr(self, 'compl_preview_btn'):
                    self.compl_preview_btn.setText("↗")
                    self.compl_preview_btn.setToolTip("Preview in main plotter")
                
                # re-enable colour picker
                self.color_selector.setDisabled(False)
        
        # Update current view tracker
        self.main_plotter_view = layer_type
        
        # Reset camera and render
        self.plotter.reset_camera()
        self.plotter.render()
    
    def popup_empty_values(self):
        msg = QMessageBox()
        msg.setWindowTitle("Empty Values")
        msg.setText("Please fill in all the required values.")
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
    
    def popup_wrong_range_values(self, variable):
        msg = QMessageBox()
        msg.setWindowTitle("Invalid Range")
        
        if variable == "Sandwich":
            msg.setText(f"Face thickness should be between 0 and half of the height: {height}.")
        elif variable == "Resolution":
            msg.setText("Resolution points should be between 50 and 1000.")
        elif variable == "Repetitions":
            msg.setText("Repetitions (x, y, z) should be between 2 and 50.")
        elif variable == "Dimensions":
            msg.setText("Dimensions should be between 0.2 and 1000.")
        elif variable == "Stretching":
            msg.setText("Stretching should be between 0.5 and 2.")
        elif variable == "Rotation":
            msg.setText("Rotation should be between 0 and 180.")
        elif variable == "Grading_Strut":
            msg.setText("When using grading for Strut, repetitions (x, y, z) should be > 2 and < 6.")
        elif variable == "Volume_Fraction_Strut":
            msg.setText("Volume Fraction range for Strut is 1 - 40%.")
        elif variable == "Volume_Fraction_Strut_Gradient":
            msg.setText("Min Volume Fraction for Strut should be less than Max Volume Fraction.")
        elif variable == "Strut_Diameter_Strut":
            msg.setText("For Strut Diameter grading, Min values should be less than Max values.")
        elif variable == "Support_Thickness":
            msg.setText("Support Thickness can't be negative.")
        elif variable == "Repetitions_Hybrid_equal":
            msg.setText("Repetitions must be equal when using Cylindrical or Spherical grading for Hybrid structures.")
        
        msg.setIcon(QMessageBox.Warning)
        msg.setStandardButtons(QMessageBox.Ok)
        msg.exec_()
    
    
    # function to show the progress popup during the saving process (for .stp files)
    def show_step_saving_popup(self):
        self.progress_popup = QMessageBox()
        self.progress_popup.setWindowTitle("Saving...")
        self.progress_popup.setText("This might take a while. Please wait while the model is being saved.")
        self.progress_popup.setIcon(QMessageBox.Information)
        self.progress_popup.setStandardButtons(QMessageBox.NoButton)
        self.progress_popup.show()
    
    def check_range_resolution(self):
        # we return 1 to stop the generation process
        if int(resolution_points) < 50 or int(resolution_points) > 1000:
            self.popup_wrong_range_values("Resolution")
            return 1
    
    def check_range_reps(self):
        if int(x_repetitions) < 2 or int(x_repetitions) > 50 or int(y_repetitions) < 2 or int(y_repetitions) > 50 or int(z_repetitions) < 2 or int(z_repetitions) > 50:
            self.popup_wrong_range_values("Repetitions")
            return 1
        
    def check_volume_fraction_strut_ranges(self):
        if float(volume_fraction_strut) > 40:
            self.popup_wrong_range_values("Volume_Fraction_Strut")
            return 1

        if grading == True:
            if float(min_volume_fraction_strut) < 1 or float(min_volume_fraction_strut) > 40 or float(max_volume_fraction_strut) < 1 or float(max_volume_fraction_strut) > 40:
                self.popup_wrong_range_values("Volume_Fraction_Strut")
                return 1
            if float(min_volume_fraction_strut) > float(max_volume_fraction_strut):
                self.popup_wrong_range_values("Volume_Fraction_Strut_Gradient")
                return 1 
    
    def check_strut_diameter_ranges(self):
        global grading
        global min_bending_radius, max_bending_radius, min_stretching_radius, max_stretching_radius
        global min_vertical_radius, max_vertical_radius, min_joint_radius, max_joint_radius
        
        if grading == True:
            if not (float(min_bending_radius) < float(max_bending_radius)):
                self.popup_wrong_range_values("Strut_Diameter_Strut")
                return 1
            if not (float(min_stretching_radius) < float(max_stretching_radius)):
                self.popup_wrong_range_values("Strut_Diameter_Strut")
                return 1
            if not (float(min_vertical_radius) < float(max_vertical_radius)):
                self.popup_wrong_range_values("Strut_Diameter_Strut")
                return 1
            if not (float(min_joint_radius) < float(max_joint_radius)):
                self.popup_wrong_range_values("Strut_Diameter_Strut")
                return 1
    
    def check_range_dims(self):
        if form == "Cylindrical":
            if float(height) < 0.2 or float(height) > 1000 or float(radius) < 0.2 or float(radius) > 1000:
                self.popup_wrong_range_values("Dimensions")
                return 1
            
        else:
            if float(height) < 0.2 or float(height) > 1000 or float(length) < 0.2 or float(length) > 1000 or float(width) < 0.2 or float(width) > 1000:
                self.popup_wrong_range_values("Dimensions")
                return 1
        
    def check_range_strch(self):
        if float(x_stretching) < 0.5 or float(x_stretching) > 2 or float(y_stretching) < 0.5 or float(y_stretching) > 2 or float(z_stretching) < 0.5 or float(z_stretching) > 2:
            self.popup_wrong_range_values("Stretching")
            return 1
    
    def check_range_rot(self):
        if float(x_rotation) < 0 or float(x_rotation) > 180 or float(y_rotation) < 0 or float(y_rotation) > 180 or float(z_rotation) < 0 or float(z_rotation) > 180:
            self.popup_wrong_range_values("Rotation")
            return 1
    
    def check_range_reps_for_hybrid(self):
        for i in range(1, int(number_of_layers) + 1):
            if int(x_rep_for_hybrid[i-1]) < 2 or int(x_rep_for_hybrid[i-1]) > 50 or int(y_rep_for_hybrid[i-1]) < 2 or int(y_rep_for_hybrid[i-1]) > 50 or int(z_rep_for_hybrid[i-1]) < 2 or int(z_rep_for_hybrid[i-1]) > 50:
                self.popup_wrong_range_values("Repetitions")
                return 1
        
        if grading_for_hybrid == "Cylindrical" or grading_for_hybrid == "Spherical":
            for i in range(1, int(number_of_layers) + 1):
                x_rep = int(x_rep_for_hybrid[i-1])
                y_rep = int(y_rep_for_hybrid[i-1])
                z_rep = int(z_rep_for_hybrid[i-1])
                if x_rep != y_rep or y_rep != z_rep:
                    self.popup_wrong_range_values("Repetitions_Hybrid_equal")
                    return 1



    # this function is used to load the Moments of Inertia dialog
    def show_moments(self):
        print("Show Moments of Inertia")
        dialog = MomentsDialog(self.moments, parent=self)
        dialog.exec_()
            
    def generate(self):
        
        # check if all the required values are filled in (some inputs need to be checked only if the advanced options are enabled)
        if type == "TPMS" or type == "Spinodal":
            if not radius or not height or not length or not width or not resolution_points or x_repetitions is None or y_repetitions is None or z_repetitions is None or x_repetitions == '' or y_repetitions == '' or z_repetitions == '': 
                self.popup_empty_values()
                return
            
            if not left_support_thick or not right_support_thick:
                self.popup_empty_values()
                return

            if left_support_thick < 0 or right_support_thick < 0:
                self.popup_wrong_range_values("Support_Thickness")
                return
            
            if advanced_options == True:
                if (x_stretching is None or y_stretching is None or z_stretching is None or
                    x_rotation is None or y_rotation is None or z_rotation is None or
                    x_stretching == '' or y_stretching == '' or z_stretching == '' or
                    x_rotation == '' or y_rotation == '' or z_rotation == ''):
                    self.popup_empty_values()
                    return
            
            # check the range of resolution points
            resol_check = self.check_range_resolution()
            if resol_check == 1:
                return
            
            # check the range of repetitions
            reps_check = self.check_range_reps()
            if reps_check == 1:
                return
            
            # check the range of dimensions
            dims_check = self.check_range_dims()
            if dims_check == 1:
                return
            
            if type == "TPMS" and advanced_options == True:
                # check the range of stretching
                strch_check = self.check_range_strch()
                if strch_check == 1:
                    return
                
                # check the range of rotation
                rot_check = self.check_range_rot()
                if rot_check == 1:
                    return
                   
        elif type == "Strut":
            if not bending_radius or not stretching_radius or not vertical_radius or not joint_radius or not resolution_points or not x_repetitions or not y_repetitions or not z_repetitions or not volume_fraction_strut or not min_volume_fraction_strut or not max_volume_fraction_strut \
                or not min_bending_radius or not max_bending_radius or not min_stretching_radius or not max_stretching_radius or not min_vertical_radius or not max_vertical_radius or not min_joint_radius or not max_joint_radius:
                self.popup_empty_values()
                return

            # check the range of resolution points
            resol_check = self.check_range_resolution()
            if resol_check == 1:
                return
            
            # check the range of repetitions
            reps_check = self.check_range_reps()
            if reps_check == 1:
                return
            
            # if Volume Fraction design is applied, check the ranges
            if design_type == "Volume Fraction":
                vol_frac_check = self.check_volume_fraction_strut_ranges()
                if vol_frac_check == 1:
                    return
            
            elif design_type == "Strut Diameter":
                rad_check = self.check_strut_diameter_ranges()
                if rad_check == 1:
                    return
                
        elif type == "Hybrid":
            if not resolution_points:
                self.popup_empty_values()
                return
            
            # for hybrid we check values for each layer
            for i in range(1, int(number_of_layers) + 1):
                if not general_topo_for_hybrid[i-1] or not specific_topo_for_hybrid[i-1] or not x_rep_for_hybrid[i-1] or not y_rep_for_hybrid[i-1] or not z_rep_for_hybrid[i-1]:
                    self.popup_empty_values()
                    return
                
            # check the range of resolution points
            resol_check = self.check_range_resolution()
            if resol_check == 1:
                return
            
            # check the range of repetitions
            reps_check = self.check_range_reps_for_hybrid()
            if reps_check == 1:
                return
                
        elif type == "Layered":
            if not resolution_points or not all(layer_densities) or x_repetitions is None or y_repetitions is None or z_repetitions is None or x_repetitions == '' or y_repetitions == '' or z_repetitions == '':
                self.popup_empty_values()
                return
            
            # check the range of resolution points
            resol_check = self.check_range_resolution()
            if resol_check == 1:
                return
            
            # check the range of repetitions
            reps_check = self.check_range_reps()
            if reps_check == 1:
                return
        
        if form == "Sandwich":
            if not bottom_face_thickness or not top_face_thickness:
                self.popup_empty_values()
                return
            
            # range = 0...height/2
            if float(bottom_face_thickness) < 0 or float(bottom_face_thickness) > float(height)/2 or float(top_face_thickness) < 0 or float(top_face_thickness) > float(height)/2:
                self.popup_wrong_range_values("Sandwich")
                return
            
        # if type == "Layered" all the layer densities should sum up to 100
        if type == "Layered" and sum(layer_densities) != 100:
            msg = QMessageBox()
            msg.setWindowTitle("Invalid Layer Densities")
            msg.setText("Layer densities should sum up to 100.")
            msg.setIcon(QMessageBox.Warning)
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return
        
        # let the check happen and show in the status bar that the generation is running
        self.status_bar.showMessage("Running Generation...")
        self.status_bar.setStyleSheet("background-color: #cce5ff; color: #004085; font-size: 16px;")
        
        # print 
        if form == "Sandwich":
            print(f"Bottom Face Thickness: {bottom_face_thickness}")
            print(f"Top Face Thickness: {top_face_thickness}")
            
        if type == "TPMS" or type == "Spinodal":
            print(f"Radius: {radius}")
            print(f"Height: {height}")
            print(f"Length: {length}")
            print(f"Width: {width}")
            print(f"Resolution Points: {resolution_points}")
            print(f"X Repetitions: {x_repetitions}")
            print(f"Y Repetitions: {y_repetitions}")
            print(f"Z Repetitions: {z_repetitions}")
            
            if type == "TPMS":
                print(f"X Stretching: {x_stretching}")
                print(f"Y Stretching: {y_stretching}")
                print(f"Z Stretching: {z_stretching}")
                print(f"X Rotation: {x_rotation}")
                print(f"Y Rotation: {y_rotation}")
                print(f"Z Rotation: {z_rotation}")
            
        elif type == "Strut":
            print(f"Bending Radius: {bending_radius}")
            print(f"Stretching Radius: {stretching_radius}")
            print(f"Vertical Radius: {vertical_radius}")
            print(f"Joint Radius: {joint_radius}")
            print(f"Resolution Points: {resolution_points}")
            print(f"X Repetitions: {x_repetitions}")
            print(f"Y Repetitions: {y_repetitions}")
            print(f"Z Repetitions: {z_repetitions}")
        
        elif type == "Hybrid":
            print(f"Resolution Points: {resolution_points}")
            print(f"X Repetitions: {x_rep_for_hybrid}")
            print(f"Y Repetitions: {y_rep_for_hybrid}")
            print(f"Z Repetitions: {z_rep_for_hybrid}")
            print(f"Volume Fractions: {vol_fraction_for_hybrid}")
            print(f"Transition Location: {transition_location}")
            print(f"Number of Layers: {number_of_layers}")
            print(f"General Topology for Hybrid: {general_topo_for_hybrid}")
            print(f"Specific Topology for Hybrid: {specific_topo_for_hybrid}")
        
        elif type == "Layered":
            print(f"Resolution Points: {resolution_points}")
            print(f"Layer Densities: {layer_densities}")

        
        # BACKEND based on the type
        if type == "TPMS":
            
            params = dict()

            if specific_topo == "Gyroid (GY)":
                params['Archi'] = "GY"
            elif specific_topo == "IWP":
                params['Archi'] = "IWP"
            elif specific_topo == "Pcell (SPC)":
                params['Archi'] = "SPC"
            elif specific_topo == "FischerkochS (FKS)":
                params['Archi'] = "FKS"
            elif specific_topo == "Schwarz_Diamond (SCD)":
                params['Archi'] = "SCD"
            elif specific_topo == "Neovius (NE)":
                params['Archi'] = "NE"
            elif specific_topo == "Lidinoid (LD)":
                params['Archi'] = "LD"
            elif specific_topo == "SchwarzHexagonal (SCH)":
                params['Archi'] = "SCH"
            elif specific_topo == "SplitP (SLP)":
                params['Archi'] = "SLP"
            elif specific_topo == "I2Y (I2Y)":
                params['Archi'] = "I2Y"
            elif specific_topo == "Fisher–KochC(S) (FKCS)":
                params['Archi'] = "FKCS"
            elif specific_topo == "F-RD (FRD)":
                params['Archi'] = "FRD"
                
            params['Type'] = general_topo
            if general_topo == "Sheet":
                params['Type'] = "TPSF"
            elif general_topo == "Skeletal":
                params['Type'] = "TPSN"
            elif general_topo == "Exoskeletal":
                params['Type'] = "TPSX"
            
            # Lattice, Sandwich, Beam or Plate
            if form == "Beam or Plate":
                params['Structure'] = "Beam"
                
                params['LSuprt_tkns'] = float(left_support_thick)
                params['RSuprt_tkns'] = float(right_support_thick)
            else:
                params['Structure'] = form
            
            params['Shape'] = form_shape
                
            if form_shape=="Cubic":
                params['a'] = float(length)
                params['b'] = float(height)
                params['c'] = float(width)
            elif form_shape=="Cylindrical":
                params['rad'] = float(radius)
                params['heit'] = float(height)
            
            params['Bfcplt_tkns'] = float(bottom_face_thickness)
            params['Tfcplt_tkns'] = float(top_face_thickness)
            
            params['nx'] = int(x_repetitions)
            params['ny'] = int(y_repetitions)
            params['nz'] = int(z_repetitions)
            
            params['Volu_Tol'] = 0.01       # with 0.001 we get some strange infinite loops
            
            params['Volume_Fraction'] = float(volume_fraction)
            
            # because the range check applies only when advanced options are enabled, a user may insert a value outside the range
            # but then disable the advanced options and the value will be used - so we need to check the range here
            if advanced_options == True:
                params['XD_Roll'] = float(x_rotation)
                params['YD_pitch'] = float(y_rotation)
                params['ZD_yaw'] = float(z_rotation)
            else:
                params['XD_Roll'] = 0
                params['YD_pitch'] = 0
                params['ZD_yaw'] = 0
            
            # because the range check applies only when advanced options are enabled, a user may insert a value outside the range
            # but then disable the advanced options and the value will be used - so we need to check the range here
            if advanced_options == True:
                params['Alph_StrchX'] = float(x_stretching)
                params['Beta_StrchY'] = float(y_stretching)
                params['Gamma_StrchZ'] = float(z_stretching)
            else:
                params['Alph_StrchX'] = 1
                params['Beta_StrchY'] = 1
                params['Gamma_StrchZ'] = 1
            
            params['MDP'] = int(resolution_points)
            
            if model_ipc_structure == True:
                params['IPC'] = "IPC_Y"
            else:
                params['IPC'] = "IPC_N"
            
            # gradation parameters
            if grading == False:
                params['Gradation'] = 'Constant'
            
            else:
                params['Gradation'] = 'Graded'
                
                if grading_direction == "X":
                    params['GradationDirection'] = 'XGraded'
                elif grading_direction == "Y":
                    params['GradationDirection'] = 'YGraded'
                elif grading_direction == "Z":
                    params['GradationDirection'] = 'ZGraded'
                
                if grading_type == "Linear":
                    params['FGType'] = 1
                elif grading_type == "Cosine":
                    params['FGType'] = 2
                elif grading_type == "Sinusoidal":
                    params['FGType'] = 3
                elif grading_type == "VType Min at Centre":
                    params['FGType'] = 4
                elif grading_type == "VType Max at Centre":
                    params['FGType'] = 5
                    
                params['Volume_Fraction1'] = float(grading_starting_volume)
                params['Volume_Fraction2'] = float(grading_ending_volume)
            
            # disable generation button
            self.generate_button.setEnabled(False)
            
            # thread connection
            params['run_parallel'] = run_parallel
            
            # we did the checks in the frontend, don't do them again in the backend
            params['from_GUI'] = True

            print(params)
            self.exec_thread = GenerationProcess("TPMS", params)
            self.exec_thread.finished.connect(self.on_generation_finished)
            self.exec_thread.error.connect(self.on_generation_error)
            self.exec_thread.start()
            
            #F,V, Final_Vol_Frac, Final_Surface = tpms_core.generate_tpms(**params)
        
        elif type == "Spinodal":
            
            params = dict()
            
            # Spinodal has fixed specific topology = SPIN 
            # Calling generate_tpms with Archi=SPIN
            params['Archi'] = "SPIN"
            
            params['Type'] = general_topo
            if general_topo == "Sheet":
                params['Type'] = "TPSF"
            elif general_topo == "Skeletal":
                params['Type'] = "TPSN"
            elif general_topo == "Exoskeletal":
                params['Type'] = "TPSX"
            
            # Lattice, Sandwich, Beam or Plate
            if form == "Beam or Plate":
                params['Structure'] = "Beam"
            else:
                params['Structure'] = form
                
            params['Shape'] = form_shape
                
            if form_shape=="Cubic":
                params['a'] = float(length)
                params['b'] = float(width)
                params['c'] = float(height)
            elif form_shape=="Cylindrical":
                params['rad'] = float(radius)
                params['heit'] = float(height)
            
            params['Bfcplt_tkns'] = bottom_face_thickness
            params['Tfcplt_tkns'] = top_face_thickness
            
            params['nx'] = int(x_repetitions)
            params['ny'] = int(y_repetitions)
            params['nz'] = int(z_repetitions)
            
            params['Volu_Tol'] = 0.01
            
            params['Volume_Fraction'] = float(volume_fraction)
            
            # roation and stretching are not used in Spinodal (just sending default values)
            params['XD_Roll'] = x_rotation
            params['YD_pitch'] = y_rotation
            params['ZD_yaw'] = z_rotation
            
            params['Alph_StrchX'] = x_stretching
            params['Beta_StrchY'] = y_stretching
            params['Gamma_StrchZ'] = z_stretching
            
            params['MDP'] = int(resolution_points)
            params['W_Tnum'] = int(waves_number)
            
            if model_ipc_structure == True:
                params['IPC'] = "IPC_Y"
            else:
                params['IPC'] = "IPC_N"    
            
            # gradation parameters
            if grading == False:
                params['Gradation'] = 'Constant'
            
            else:
                params['Gradation'] = 'Graded'
                
                if grading_direction == "X":
                    params['GradationDirection'] = 'XGraded'
                elif grading_direction == "Y":
                    params['GradationDirection'] = 'YGraded'
                elif grading_direction == "Z":
                    params['GradationDirection'] = 'ZGraded'
                
                if grading_type == "Linear":
                    params['FGType'] = 1
                elif grading_type == "Cosine":
                    params['FGType'] = 2
                elif grading_type == "Sinusoidal":
                    params['FGType'] = 3
                elif grading_type == "VType Min at Centre":
                    params['FGType'] = 4
                elif grading_type == "VType Max at Centre":
                    params['FGType'] = 5
                    
                params['Volume_Fraction1'] = float(grading_starting_volume)
                params['Volume_Fraction2'] = float(grading_ending_volume)
            
            # disable generation button
            self.generate_button.setEnabled(False)
            
            params['run_parallel'] = run_parallel

            # we did the checks in the frontend, don't do them again in the backend
            params['from_GUI'] = True

            print(params)
            self.exec_thread = GenerationProcess("Spinodal", params)
            self.exec_thread.finished.connect(self.on_generation_finished)
            self.exec_thread.error.connect(self.on_generation_error)
            self.exec_thread.start()
            
            #F,V, Final_Vol_Frac, Final_Surface = tpms_core.generate_spin(**params)
             
        elif type == "Strut":
            
            params = dict()
            
            # the specific topo is not used in Strut
            params['Archi'] = ""
            
            # Design Type
            if design_type == "Strut Diameter":
                params['DesignType'] = "StrutD"
            elif design_type == "Volume Fraction":
                params['DesignType'] = "VolFracBased"
            
            params['Type'] = general_topo
            if general_topo == "Sheet":
                params['Type'] = "TPSF"
            elif general_topo == "Skeletal":
                params['Type'] = "TPSN"
            elif general_topo == "Exoskeletal":
                params['Type'] = "TPSX"
            
            params['Structure'] = form
            if form == "Lattice":
                params['Structure'] = "Cuboid"
            
            # some standard values for STRUT
            params['BOX'] = 0
            params['BOXWIDTH'] = 1
            params['SANDWICH'] = 0
            
            params['BLayerTckns'] = bottom_face_thickness
            params['TLayerTckns'] = top_face_thickness
            
            # lattice family based on the specific topo
            if specific_topo == "Standard Octer Lattice":
                params['latticeFamily'] = 5
            elif specific_topo == "Reinforced Octet":
                params['latticeFamily'] = 6
            elif specific_topo == "Octahedral":
                params['latticeFamily'] = 7
            elif specific_topo == "Reinforced Octahedral":
                params['latticeFamily'] = 8
            elif specific_topo == "Circular Octahedral":
                params['latticeFamily'] = 9
            elif specific_topo == "BCC":
                params['latticeFamily'] = 10
            elif specific_topo == "FCC":
                params['latticeFamily'] = 11
            elif specific_topo == "Dodecahedron":
                params['latticeFamily'] = 12
                
            
            # fixed for STRUT
            params['sampleName'] = "StrutSample"
            
            params['xRep'] = int(y_repetitions)
            params['yRep'] = int(x_repetitions)
            params['zRep'] = int(z_repetitions)
            
            params['sx'] = float(width)
            params['sy'] = float(length)
            params['sz'] = float(height)
            
            params['MDP'] = int(resolution_points)
            params['finalLatticeRes'] = int(resolution_points)
            
            params['strutStretch'] = 1
            
            # grading is different for STRUT
            if grading == False:
                params['GradationDirection'] = 'Constant'
                params['FGType'] = 0 # constant
                
                if design_type == "Strut Diameter":
                    # these 4 are needed if we don't have grading
                    params['bendingStrutRadius'] = float(bending_radius) / 100.0
                    params['stretchingStrutRadius'] = float(stretching_radius) / 100.0
                    params['verticalStrutRadius'] = float(vertical_radius) / 100.0
                    params['jointRadius'] = float(joint_radius) / 100.0
                
                elif design_type == "Volume Fraction":
                    params['volumeFraction'] = float(volume_fraction_strut)
            
            elif grading == True:
                if grading_direction == "X":
                    params['GradationDirection'] = 'XGraded'
                elif grading_direction == "Y":
                    params['GradationDirection'] = 'YGraded'
                elif grading_direction == "Z":
                    params['GradationDirection'] = 'ZGraded'
                
                if grading_type == "Linear":
                    params['FGType'] = 1
                elif grading_type == "Cosine":
                    params['FGType'] = 2
                elif grading_type == "Sinusoidal":
                    params['FGType'] = 3
                elif grading_type == "VType Min at Centre":
                    params['FGType'] = 4
                elif grading_type == "VType Max at Centre":
                    params['FGType'] = 5
                
                if design_type == "Strut Diameter":
                    # these 8 types are needed if we have grading
                    params['bendingStrutRadius_min'] = float(min_bending_radius) / 100.0
                    params['bendingStrutRadius_max'] = float(max_bending_radius) / 100.0
                    params['stretchingStrutRadius_min'] = float(min_stretching_radius) / 100.0
                    params['stretchingStrutRadius_max'] = float(max_stretching_radius) / 100.0
                    params['verticalStrutRadius_min'] = float(min_vertical_radius) / 100.0
                    params['verticalStrutRadius_max'] = float(max_vertical_radius) / 100.0
                    params['jointRadius_min'] = float(min_joint_radius) / 100.0
                    params['jointRadius_max'] = float(max_joint_radius) / 100.0
                elif design_type == "Volume Fraction":
                    params['volumeFraction_min'] = float(min_volume_fraction_strut)
                    params['volumeFraction_max'] = float(max_volume_fraction_strut)
            
            if model_ipc_structure == True:
                params['IPC'] = "IPC_Y"
            else:
                params['IPC'] = "IPC_N"
                
            params['run_parallel'] = run_parallel

            # we did the checks in the frontend, don't do them again in the backend
            params['from_GUI'] = True
            
            # disable generation button
            self.generate_button.setEnabled(False)
            
            print(params)
            self.exec_thread = GenerationProcess("Strut", params)
            self.exec_thread.finished.connect(self.on_generation_finished)
            self.exec_thread.error.connect(self.on_generation_error)
            self.exec_thread.start()
        
            #F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_strut(**params)
            
        elif type == "Hybrid":
            
            params = dict()
            
            params['Hybrid_layers'] = int(number_of_layers)
            
            # this is standard for HYBRID
            params['Base_class'] = ['TPMS'] * int(number_of_layers)

            params['Archi'] = specific_topo_for_hybrid
            for i, topo in enumerate(specific_topo_for_hybrid):
                if topo == "Gyroid (GY)":
                    params['Archi'][i] = "GY"
                elif topo == "IWP":
                    params['Archi'][i] = "IWP"
                elif topo == "Pcell (SPC)":
                    params['Archi'][i] = "SPC"
                elif topo == "FischerkochS (FKS)":
                    params['Archi'][i] = "FKS"
                elif topo == "Schwarz_Diamond (SCD)":
                    params['Archi'][i] = "SCD"
                elif topo == "Neovius (NE)":
                    params['Archi'][i] = "NE"
                elif topo == "Lidinoid (LD)":
                    params['Archi'][i] = "LD"
                elif topo == "SchwarzHexagonal (SCH)":
                    params['Archi'][i] = "SCH"
                elif topo == "SplitP (SLP)":
                    params['Archi'][i] = "SLP"
                elif topo == "I2Y (I2Y)":
                    params['Archi'][i] = "I2Y"
                elif topo == "Fisher–KochC(S) (FKCS)":
                    params['Archi'][i] = "FKCS"
                elif topo == "F-RD (FRD)":
                    params['Archi'][i] = "FRD"
            
            params['Type'] = general_topo_for_hybrid
            for i, topo in enumerate(general_topo_for_hybrid):
                if topo == "Sheet":
                    params['Type'][i] = "TPSF"
                elif topo == "Skeletal":
                    params['Type'][i] = "TPSN"
                elif topo == "Exoskeletal":
                    params['Type'][i] = "TPSX"
                    
            if grading_for_hybrid == "Linear":
                params['HybridType'] = 1
            elif grading_for_hybrid == "Cylindrical":
                params['HybridType'] = 2
                
                # if "Cylindrical" we also need to specify the type of hybridization
                if cylindrical_hybrid_type == 0:
                    params['CylHyType'] = 1
                elif cylindrical_hybrid_type == 1:
                    params['CylHyType'] = 2
                elif cylindrical_hybrid_type == 2:
                    params['CylHyType'] = 3
                else: 
                    print(f"Unknown cylindrical hybrid type: {cylindrical_hybrid_type}")
                
            elif grading_for_hybrid == "Spherical":
                params['HybridType'] = 3
            
            
            
            params['a'] = float(length)
            params['b'] = float(width)
            params['c'] = float(height)
            
            # 1 value of x,y,z for each layer (we need them as integers)      
            params['nx'] = [int(x) for x in x_rep_for_hybrid]
            params['ny'] = [int(y) for y in y_rep_for_hybrid]
            params['nz'] = [int(z) for z in z_rep_for_hybrid]

            params['Volu_Tol'] = 0.01
            
            # volume fraction for each layer
            vol_fraction_in = [float(i) for i in vol_fraction_for_hybrid]
            params['Volume_Fraction'] = vol_fraction_in
            
            params['MDP'] = int(resolution_points)
            params['W_Tnum'] = 1000
            
            params['trans_quality'] = transition_quality
            
            # for the transition location we need to have a value for each layer except the last one
            #params['trans_position'] = float(transition_location) * int(length)
            # params['trans'] = float(transition_location) * float(length)
            params['trans'] = float(transition_location) * float(height)
            
            # disable generation button
            self.generate_button.setEnabled(False)
            
            if model_ipc_structure == True:
                params['IPC'] = "IPC_Y"
            else:
                params['IPC'] = "IPC_N"
            params['run_parallel'] = run_parallel

            # we did the checks in the frontend, don't do them again in the backend
            params['from_GUI'] = True
    
            print(params)
            self.exec_thread = GenerationProcess("Hybrid", params)
            self.exec_thread.finished.connect(self.on_generation_finished)
            self.exec_thread.error.connect(self.on_generation_error)
            self.exec_thread.start()
            
            #F, V, Final_Vol_Frac, Final_Surface = tpms_core.generate_hybrid(**params)

        elif type == "Layered":
        
            params = dict()

            params['Num_layers'] = int(number_of_layers)
            
            # this is standard for LAYERED (we dont need it multiplied by the number of layers)
            params['Base_class'] = 'TPMS'

            # Layered has 1 specific topo for all layers
            if specific_topo == "Gyroid (GY)":
                params['Archi'] = "GY"
            elif specific_topo == "IWP":
                params['Archi'] = "IWP"
            elif specific_topo == "Pcell (SPC)":
                params['Archi'] = "SPC"
            elif specific_topo == "FischerkochS (FKS)":
                params['Archi'] = "FKS"
            elif specific_topo == "Schwarz_Diamond (SCD)":
                params['Archi'] = "SCD"
            elif specific_topo == "Neovius (NE)":
                params['Archi'] = "NE"
            elif specific_topo == "Lidinoid (LD)":
                params['Archi'] = "LD"
            elif specific_topo == "SchwarzHexagonal (SCH)":
                params['Archi'] = "SCH"
            elif specific_topo == "SplitP (SLP)":
                params['Archi'] = "SLP"
            elif specific_topo == "I2Y (I2Y)":
                params['Archi'] = "I2Y"
            elif specific_topo == "Fisher–KochC(S) (FKCS)":
                params['Archi'] = "FKCS"
            elif specific_topo == "F-RD (FRD)":
                params['Archi'] = "FRD"
            
            params['Type'] = general_topo
            if general_topo == "Sheet":
                params['Type'] = "TPSF"
            elif general_topo == "Skeletal":
                params['Type'] = "TPSN"
            elif general_topo == "Exoskeletal":
                params['Type'] = "TPSX"
            
            # params['Structure'] = form
            # if form == "Lattice":
            #     params['Structure'] = "Cuboid"
            
            
            params['a'] = float(length)
            params['b'] = float(width)
            params['c'] = float(height)
            
            # params['Bfcplt_tkns'] = bottom_face_thickness
            # params['Tfcplt_tkns'] = top_face_thickness
            
            params['nx'] = int(x_repetitions)
            params['ny'] = int(y_repetitions)
            params['nz'] = int(z_repetitions)
            
            params['Volu_Tol'] = 0.01

            # we actually don't need a separate value in the GUI for the volume fraction
            params['Volume_Fraction'] = float(volume_fraction)
            
            # for "Layered" use default values for rotation, stretching
            params['XD_Roll'] = x_rotation
            params['YD_pitch'] = y_rotation
            params['ZD_yaw'] = z_rotation
            
            params['Alph_StrchX'] = x_stretching
            params['Beta_StrchY'] = y_stretching
            params['Gamma_StrchZ'] = z_stretching
            
            params['MDP'] = int(resolution_points)
            params['Layer_density'] = layer_densities
            params['IPC'] = "IPC_N"
            
            # fixed for LAYERED
            params['W_Tnum'] = 1000
            
            # disable generation button
            self.generate_button.setEnabled(False)
            
            params['run_parallel'] = run_parallel

            # we did the checks in the frontend, don't do them again in the backend
            params['from_GUI'] = True

            print(params)
            self.exec_thread = GenerationProcess("Layered", params)
            self.exec_thread.finished.connect(self.on_generation_finished)
            self.exec_thread.error.connect(self.on_generation_error)
            self.exec_thread.start()
            
            #F0, V0, F1, V1, Final_Vol_Frac, Final_Surface = tpms_core.generate_layered(**params)
            
    def save_cardfile(self):
        # ask the user what type of file they want to save & the path for it
        file_path, _ = QFileDialog.getSaveFileName(
            self, 
            "Save File", 
            "", 
            "STL Files (*.stl);;STEP Files (*.stp *.step);;All Files (*)"
        )
        
        if not file_path:
            return
            
        print(file_path)
        
        # handle Layered case separately because it has two meshes
        try:
            if hasattr(self, 'F0') and hasattr(self, 'V0') and hasattr(self, 'F1') and hasattr(self, 'V1'):
                if file_path.endswith(".stl"):
                    base_path = file_path[:-4] if file_path.endswith('.stl') else file_path
                    file0 = f"{base_path}_0.stl"
                    file1 = f"{base_path}_1.stl"
                    
                    # use your existing create_pyvista_mesh function to ensure consistency
                    mesh0 = self.create_pyvista_mesh(self.F0, self.V0)
                    mesh1 = self.create_pyvista_mesh(self.F1, self.V1)
                    
                    mesh0.save(file0)
                    mesh1.save(file1)
                    
                    QMessageBox.information(
                        self, 
                        "Success", 
                        f"Layered model successfully saved to:\n{file0}\n{file1}"
                    )
                    return
                    
                elif file_path.endswith(".stp") or file_path.endswith(".step"):
                    # create temporary directory for intermediate files
                    import tempfile
                    temp_dir = tempfile.mkdtemp()
                    
                    try:
                        # create both meshes
                        mesh0 = self.create_pyvista_mesh(self.F0, self.V0)
                        mesh1 = self.create_pyvista_mesh(self.F1, self.V1)
                        
                        # save as temporary STL files
                        stl_file0 = os.path.join(temp_dir, "temp0.stl")
                        stl_file1 = os.path.join(temp_dir, "temp1.stl")
                        mesh0.save(stl_file0)
                        mesh1.save(stl_file1)
                        
                        # initialize STEP writer
                        step_writer = STEPControl_Writer()
                        shape0 = TopoDS_Shape()
                        shape1 = TopoDS_Shape()
                        stl_reader = StlAPI_Reader()
                        
                        # read both STL files
                        if not stl_reader.Read(shape0, stl_file0):
                            raise Exception("Failed to read first layer STL file")
                        if not stl_reader.Read(shape1, stl_file1):
                            raise Exception("Failed to read second layer STL file")
                        
                        # transfer both shapes to STEP
                        if not step_writer.Transfer(shape0, STEPControl_AsIs):
                            raise Exception("Failed to transfer first layer to STEP")
                        if not step_writer.Transfer(shape1, STEPControl_AsIs):
                            raise Exception("Failed to transfer second layer to STEP")
                        
                        # show progress popup
                        self.show_step_saving_popup()
                        QApplication.processEvents()
                        
                        # write combined STEP file
                        if not step_writer.Write(file_path):
                            raise Exception("Failed to write STEP file")
                        
                        QMessageBox.information(
                            self, 
                            "Success", 
                            f"Layered model successfully saved to {file_path}"
                        )
                        
                    except Exception as e:
                        QMessageBox.critical(
                            self, 
                            "Error", 
                            f"Failed to save STEP file: {str(e)}"
                        )
                        return
                        
                    finally:
                        # clean up temporary files
                        if os.path.exists(stl_file0):
                            os.remove(stl_file0)
                        if os.path.exists(stl_file1):
                            os.remove(stl_file1)
                        if os.path.exists(temp_dir):
                            os.rmdir(temp_dir)
                        if hasattr(self, 'progress_popup'):
                            self.progress_popup.done(0)
                    
                    return
                    
            # original handling for non-layered cases
            dataset = self.plotter.renderer.GetActors().GetLastActor().GetMapper().GetInput()
            mesh = pv.wrap(dataset)
                
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save file: {str(e)}")
            return
        
        # handle all other cases (TPMS, Spinodal, Strut, Hybrid), Layered don't reach here
        if file_path.endswith(".stl"):
            if mesh is not None:
                mesh.save(file_path)
                QMessageBox.information(self, "Success", f"Model successfully saved to {file_path}")
            else:
                QMessageBox.critical(self, "Error", "Failed to save file: No mesh data available.")
            return
        
        elif file_path.endswith(".stp") or file_path.endswith(".step"):
            if hasattr(self, 'type') and self.type == "Strut":
                QMessageBox.critical(self, "Error", "Failed to save file: STEP format is not supported for Strut type.")
                return
                
            self.show_step_saving_popup()
            QApplication.processEvents()
        
            stl_file = file_path.replace(".stp", ".stl").replace(".step", ".stl")
            mesh.save(stl_file)
            
            shape = TopoDS_Shape()
            stl_reader = StlAPI_Reader()
            
            if mesh is not None:
                status = stl_reader.Read(shape, stl_file)
            else:
                QMessageBox.critical(self, "Error", "Failed to save file: No mesh data available.")
                return
            
            step_writer = STEPControl_Writer()
            status = step_writer.Transfer(shape, STEPControl_AsIs)
            
            if status > 0:
                try: 
                    status = step_writer.Write(file_path)
                    QMessageBox.information(self, "Success", f"Model successfully saved to {file_path}")
                
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to save file: {str(e)}")
                    return
                
                finally:
                    os.remove(stl_file)
                    if hasattr(self, 'progress_popup'):
                        self.progress_popup.done(0)
            else:
                raise Exception("Failed to transfer shape to STEP format")
            
            return
    
    # function used by the gen thread to update the self F,V .. to be used in the GUI
    def on_generation_finished(self, result):
        
        # enable the generation button
        self.generate_button.setEnabled(True)
        
        if result is not None:
            try:
                if type == "TPMS":
                    
                    if model_ipc_structure == True:
                        self.F0 = result["F_reinf"]
                        self.V0 = result["V_reinf"]
                        # self.F = self.F0    # utilize them for the moments out the if/else
                        # self.V = self.V0
                        self.F1 = result["F_compl"]
                        self.V1 = result["V_compl"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F0, self.V0, self.F1, self.V1)
                        
                    else:
                        self.F  = result["F"]
                        self.V  = result["V"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                    
                        # calculate the moments & enable the button, update surface and volume labels
                        try:
                            self.moments = tpms_moments.get_moments(self.F, self.V)
                            print("moments=", self.moments)
                            self.moments_button.setEnabled(True)    # it changes stylesheet too, because on init we set 2
                                                                    # styles based on the stage (enabled/disabled)
                        except Exception as e:
                            print("Error in moments calculation: ", e)
                        
                        # if IPC_Y this is the main right window mesh update (if IPC_N thats the only one we need)
                        self.update_mesh(self.F, self.V)
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                    
                    self.status_bar.showMessage("Generation completed successfully.")
                    self.status_bar.setStyleSheet("background-color: #d4edda; color: #155724; font-size: 16px;")
                    
                elif type == "Spinodal":
                    
                    if model_ipc_structure == True:
                        self.F0 = result["F_reinf"]
                        self.V0 = result["V_reinf"]
                        # self.F = self.F0    # utilize them for the moments out the if/else
                        # self.V = self.V0
                        self.F1 = result["F_compl"]
                        self.V1 = result["V_compl"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        self.update_mesh(self.F0, self.V0, self.F1, self.V1)
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                    
                    else:
                        self.F  = result["F"]
                        self.V  = result["V"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        # calculate the moments & enable the button, update surface and volume labels, update the mesh
                        try:
                            self.moments = tpms_moments.get_moments(self.F, self.V)
                            print("moments=", self.moments)
                            self.moments_button.setEnabled(True)
                        except Exception as e:
                            print("Error in moments calculation: ", e)
                    
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F, self.V)
                    
                    self.status_bar.showMessage("Generation completed successfully.")
                    self.status_bar.setStyleSheet("background-color: #d4edda; color: #155724; font-size: 16px;")
                    
                elif type == "Strut":
                    
                    if model_ipc_structure == True:
                        self.F0 = result["F_reinf"]
                        self.V0 = result["V_reinf"]
                        # self.F = self.F0    # utilize them for the moments out the if/else
                        # self.V = self.V0
                        self.F1 = result["F_compl"]
                        self.V1 = result["V_compl"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F0, self.V0, self.F1, self.V1)
                    else:
                        self.F  = result["F"]
                        self.V  = result["V"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                    
                        # calculate the moments & enable the button, update surface and volume labels, update the mesh
                        try:
                            self.moments = tpms_moments.get_moments(self.F, self.V)
                            print("moments=", self.moments)
                            self.moments_button.setEnabled(True)
                        except Exception as e:
                            print("Error in moments calculation: ", e)
                    
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F, self.V)
                    
                    self.status_bar.showMessage("Generation completed successfully.")
                    self.status_bar.setStyleSheet("background-color: #d4edda; color: #155724; font-size: 16px;")
                    
                elif type == "Hybrid":
                    
                    if model_ipc_structure == True:
                        self.F0 = result["F_reinf"]
                        self.V0 = result["V_reinf"]
                        # self.F = self.F0    # utilize them for the moments out the if/else
                        # self.V = self.V0
                        self.F1 = result["F_compl"]
                        self.V1 = result["V_compl"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F0, self.V0, self.F1, self.V1)
                    else:
                        self.F  = result["F"]
                        self.V  = result["V"]
                        self.Final_Vol_Frac = result["Final_Vol_Frac"]
                        self.Final_Surface = result["Final_Surface"]
                        
                        # calculate the moments & enable the button, update surface and volume labels, update the mesh
                        try:
                            self.moments = tpms_moments.get_moments(self.F, self.V)
                            print("moments=", self.moments)
                            self.moments_button.setEnabled(True)
                        except Exception as e:
                            print("Error in moments calculation: ", e)
                        
                        self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                        self.update_mesh(self.F, self.V)
                    
                    self.status_bar.showMessage("Generation completed successfully.")
                    self.status_bar.setStyleSheet("background-color: #d4edda; color: #155724; font-size: 16px;")
                
                elif type == "Layered":
                    self.F0  = result["F0"]
                    self.V0  = result["V0"]
                    self.F1  = result["F1"]
                    self.V1  = result["V1"]
                    self.Final_Vol_Frac = result["Final_Vol_Frac"]
                    self.Final_Surface = result["Final_Surface"]
                    
                    # surface, volume and mesh update
                    self.update_surface_and_volume(self.Final_Vol_Frac, self.Final_Surface)
                    self.update_mesh_layered(self.F0, self.V0, self.F1, self.V1)
                    
                    self.status_bar.showMessage("Generation completed successfully.")
                    self.status_bar.setStyleSheet("background-color: #d4edda; color: #155724; font-size: 16px;")
                    
            except KeyError as e:
                self.status_bar.showMessage(f"Missing key in result: {str(e)}")
                self.status_bar.setStyleSheet("background-color: #f8d7da; color: #721c24; font-size: 16px;")
        
        else:
            self.status_bar.showMessage("Generation returned no results.")
            self.status_bar.setStyleSheet("background-color: #f8d7da; color: #721c24; font-size: 16px;")
                
    # function used by the gen thread to return an error if it occurs
    def on_generation_error(self, error):
        self.status_bar.showMessage(f"Generation failed: {error}")
        self.status_bar.setStyleSheet("background-color: #f8d7da; color: #721c24; font-size: 16px;")
        
        print("Generation error:", error)
        # enable the generation button
        self.generate_button.setEnabled(True)
    
    def update_surface_and_volume(self, Final_Vol_Frac, Final_Surface):
        self.total_surface_label.setText(f"Surface Area: {Final_Surface:.2f} mm²")
        self.total_vol_fraction_label.setText(f"Volume: {Final_Vol_Frac:.2f} %")
        
    def closeEvent(self, event):
        # check if process exists and is running
        if hasattr(self, 'exec_thread') and self.exec_thread.process is not None and self.exec_thread.process.poll() is None:
            reply = QMessageBox.question(
                self,
                "Confirm Exit",
                "Generation is still in progress. Do you really want to close the application?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                event.ignore()
                return
            else:
                if hasattr(self, 'exec_thread'):
                    # Stop the process
                    self.exec_thread.stop()
                    print("Generation process terminated.")
        
        # CRITICAL: Clean up VTK widgets before accepting close event (need that specifically for windows)
        self.cleanup_vtk_widgets()
        
        # because we hide the first window, we now find the first window to bring it at front
        app = QApplication.instance()
        if app:
            for widget in app.allWidgets():
                if isinstance(widget, TPMSInterface) and widget != self:
                    widget.show()
                    widget.raise_()             # Bring to front
                    widget.activateWindow()     # Make it active
                    break
        
        event.accept()

    def cleanup_vtk_widgets(self):
        """Clean up all VTK widgets to prevent OpenGL context errors on Windows"""
        try:
            # Clean up the PyVista plotter
            if hasattr(self, 'plotter') and self.plotter is not None:
                try:
                    # PyVista specific cleanup
                    self.plotter.close()  # This is critical for PyVista
                    self.plotter = None
                    #print("PyVista plotter cleaned up")
                except Exception as e:
                    print(f"Error cleaning up PyVista plotter: {e}")
            
            # Find all VTK widgets recursively
            vtk_widgets = []
            
            # Common VTK widget class names to look for
            vtk_widget_types = [
                'QVTKRenderWindowInteractor',
                'QVTKOpenGLNativeWidget',
                'QVTKOpenGLWidget',
                'QVTKWidget',
                'QtInteractor'
            ]
            
            def find_vtk_widgets(widget):
                for child in widget.findChildren(QWidget):
                    class_name = child.__class__.__name__
                    if any(vtk_type in class_name for vtk_type in vtk_widget_types):
                        vtk_widgets.append(child)
                    # Also check if widget has VTK methods
                    elif hasattr(child, 'GetRenderWindow'):
                        vtk_widgets.append(child)
            
            find_vtk_widgets(self)
            
            # Clean up each VTK widget
            for widget in vtk_widgets:
                try:
                    # For PyVista QtInteractor widgets
                    if hasattr(widget, 'close'):
                        widget.close()
                        
                    if hasattr(widget, 'GetRenderWindow'):
                        render_window = widget.GetRenderWindow()
                        if render_window:
                            # Finalize the render window to clean up OpenGL context
                            render_window.Finalize()
                            
                    if hasattr(widget, 'GetInteractor'):
                        interactor = widget.GetInteractor()
                        if interactor:
                            interactor.TerminateApp()
                            
                except Exception as e:
                    print(f"Error cleaning up VTK widget: {e}")
                    
        except Exception as e:
            print(f"Error during VTK cleanup: {e}")


class AboutDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("About Top-6-Class MetaStudio")
        self.setFixedSize(420, 260)  # non-resizable
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowContextHelpButtonHint)

        layout = QVBoxLayout(self)

        label = QLabel(
            "<b>Top-6-Class MetaStudio (Top6Meta)</b><br>"
            "Version 1.0<br><br>"
            "A Python HPC framework for modeling architected materials "
            "and metastructures.<br><br>"
            "<b>Authors:</b><br>"
            "Agyapal Singh<br>"
            "Georgios Mermigkis<br>"
            "Panagiotis Hadjidoukas<br>"
            "Nikolaos Karathanasopoulos<br><br>"
            "© 2026"
        )
        label.setTextFormat(Qt.RichText)
        label.setAlignment(Qt.AlignLeft)
        label.setWordWrap(True)

        layout.addWidget(label)

        btn = QPushButton("OK")
        btn.setFixedWidth(80)
        btn.clicked.connect(self.accept)
        btn.setDefault(True)

        layout.addWidget(btn, alignment=Qt.AlignCenter)

def main():
    app = QApplication(sys.argv)
    
    # set fusion style for a modern look
    app.setStyle('Fusion')
    
    # light mode palette
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(240, 248, 255))        # AliceBlue
    palette.setColor(QPalette.WindowText, QColor(70, 130, 180))     # SteelBlue
    app.setPalette(palette)
    
    window = TPMSInterface()
    
    # make the window non-resizable - DISABLED for now
    #window.setFixedSize(1000, 750)
    
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()