import math
import numpy as np
import os

from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QPushButton,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QMessageBox,
    QComboBox
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import QApplication

import vtk
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from ui.resource_locate import resource_path


# --- Helper functions for color interpolation ---
def lerp(a, b, t):
    return a + (b - a) * t

def color_lerp(c1, c2, t):
    return [lerp(c1[i], c2[i], t) for i in range(3)]

def get_coolwarm_color(value):
    """
    value in [0,1].
    Interpolates between:
        0.0: [0.2298, 0.2987, 0.7537]  (cool/blue)
        0.5: [0.8650, 0.8650, 0.8650]  (light gray)
        1.0: [0.7050, 0.0150, 0.1499]  (warm/red)
    """
    v = max(0.0, min(1.0, value))

    cool_color = [0.2298, 0.2987, 0.7537]
    mid_color  = [0.8650, 0.8650, 0.8650]
    warm_color = [0.7050, 0.0150, 0.1499]

    if v <= 0.5:
        scale = (v - 0.0) / 0.5
        return color_lerp(cool_color, mid_color, scale)
    else:
        scale = (v - 0.5) / 0.5
        return color_lerp(mid_color, warm_color, scale)


class ProteinViewerPage(QWidget):
    """
    ProteinViewerPage now supports:
      - Two display modes ("Global Distance" vs. "Local Edge-Sum").
      - Distance-based label hiding when camera zooms out.
      - Locating a particular node via a drop-down + "Locate Node" button.
      - Black edges and gray edge-value text.
    """
    def __init__(self, mer_name, pdb_content, interactions, total_weight_sums, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.interactions = interactions
        self.total_weight_sums = total_weight_sums  # Distances from the best source
        self.on_back = on_back

        self.mers = {}
        self.positions = None
        self.mer_names = []

        # Store a dictionary for "edge sum" values
        self.edge_sum_by_node = {}

        # Current mode can be 'distance' or 'edge_sum'
        self.display_mode = 'distance'

        # We'll keep references to label actors so we can hide/show them dynamically
        self.node_text_actors = []
        self.edge_text_actors = []

        # Thresholds for distance-based label hiding:
        self.node_label_hide_dist = 3.0
        self.edge_label_hide_dist = 3.0

        # --- Main layout ---
        self.main_layout = QVBoxLayout(self)

        # --- Top Bar ---
        top_bar = QHBoxLayout()

        # Back button
        back_btn = QPushButton("Back")
        back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(back_btn, alignment=Qt.AlignLeft)

        # Current mode label
        self.mode_label = QLabel("Current Mode: Global Distance")
        self.mode_label.setFont(QFont("Arial", 11, QFont.Bold))
        top_bar.addWidget(self.mode_label, alignment=Qt.AlignCenter)

        # Info button
        info_btn = QPushButton("Info")
        info_btn.setToolTip("Learn about the display modes and visualization details.")
        info_btn.clicked.connect(self.show_info)
        top_bar.addWidget(info_btn, alignment=Qt.AlignRight)

        # Global Distance Mode button
        self.distance_btn = QPushButton("Global Distance Mode")
        self.distance_btn.clicked.connect(self.set_distance_mode)
        self.distance_btn.setToolTip("Node values = total distance from the best source Mer.")
        top_bar.addWidget(self.distance_btn, alignment=Qt.AlignRight)

        # Local Edge-Sum Mode button
        self.edge_sum_btn = QPushButton("Local Edge-Sum Mode")
        self.edge_sum_btn.clicked.connect(self.set_edge_sum_mode)
        self.edge_sum_btn.setToolTip("Node values = sum of connected edges' weights.")
        top_bar.addWidget(self.edge_sum_btn, alignment=Qt.AlignRight)

        # --- Node Locator: drop-down + button ---
        self.node_select_combobox = QComboBox()
        self.node_select_combobox.setEnabled(False)
        top_bar.addWidget(self.node_select_combobox)

        self.locate_btn = QPushButton("Locate Node")
        self.locate_btn.setEnabled(False)
        self.locate_btn.clicked.connect(self.locate_selected_node)
        top_bar.addWidget(self.locate_btn, alignment=Qt.AlignRight)

        top_bar.addStretch()
        self.main_layout.addLayout(top_bar)

        # --- Progress Label + Bar ---
        self.progress_label = QLabel("Processing graph...")
        self.progress_label.setAlignment(Qt.AlignCenter)
        self.progress_label.setFont(QFont("Arial", 14, QFont.Bold))
        self.main_layout.addWidget(self.progress_label, alignment=Qt.AlignCenter)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(self.progress_bar, alignment=Qt.AlignCenter)

        # --- Loading Image ---
        loading_image_path = resource_path(os.path.join("assets", "images", "protein.png"))
        self.loading_label = QLabel()
        if os.path.exists(loading_image_path):
            pixmap = QPixmap(loading_image_path)
            scaled_pixmap = pixmap.scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self.loading_label.setPixmap(scaled_pixmap)
        else:
            self.loading_label.setText("Loading...")
        self.loading_label.setAlignment(Qt.AlignCenter)
        self.main_layout.addWidget(self.loading_label, alignment=Qt.AlignCenter)

        # --- VTK Setup ---
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)  # White background
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()

        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        self.vtk_widget.hide()  # hidden until data is drawn

        # Schedule parsing and drawing after a short delay
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    def show_info(self):
        """ Show a popup describing the two modes and label-hiding logic. """
        info_text = (
            "This 3D viewer supports two display modes:\n\n"
            "1) Global Distance Mode:\n"
            "   Each node's value is the total distance (or weight) from the best source Mer.\n\n"
            "2) Local Edge-Sum Mode:\n"
            "   Each node's value is the sum of the weights of edges connected to it.\n\n"
            "DISTANCE-BASED LABEL HIDING:\n"
            "   Labels automatically hide when you zoom out, and show when you zoom in.\n\n"
            "LOCATE NODE:\n"
            "   Choose a node from the drop-down and click 'Locate Node' to center on it."
        )
        QMessageBox.information(self, "Display Modes", info_text)

    def set_distance_mode(self):
        """Switch to 'distance' display mode."""
        self.display_mode = 'distance'
        self.mode_label.setText("Current Mode: Global Distance")
        self.redraw_scene()

    def set_edge_sum_mode(self):
        """Switch to 'edge_sum' display mode."""
        self.display_mode = 'edge_sum'
        self.mode_label.setText("Current Mode: Local Edge-Sum")
        self.redraw_scene()

    def parse_and_draw(self):
        """Parses PDB content, initializes node data, draws the scene."""
        self.progress_bar.setValue(10)
        QApplication.processEvents()

        self.mers, _ = self.parse_pdb_content(self.pdb_content)
        self.mer_names = list(self.mers.keys())

        self.progress_bar.setValue(30)
        QApplication.processEvents()

        if not self.mers:
            self.progress_label.setText("No Mers found.")
            self.progress_bar.setValue(100)
            return

        # Normalize positions
        positions = np.array([self.mers[name]['position'] for name in self.mer_names])
        positions -= positions.mean(axis=0)
        scale = np.max(np.linalg.norm(positions, axis=1))
        if scale > 0:
            positions /= scale
        self.positions = positions

        # Edge sums
        self.edge_sum_by_node = {}
        for name in self.mer_names:
            self.edge_sum_by_node[name] = 0.0

        for interaction in self.interactions:
            from_mer = interaction.from_mer
            to_mer = interaction.to_mer
            weight = interaction.weight
            if from_mer in self.edge_sum_by_node:
                self.edge_sum_by_node[from_mer] += weight
            if to_mer in self.edge_sum_by_node:
                self.edge_sum_by_node[to_mer] += weight

        # Populate the combo box for "Locate Node" now that we know mer_names
        self.node_select_combobox.clear()
        for n in self.mer_names:
            self.node_select_combobox.addItem(n)
        self.node_select_combobox.setEnabled(True)
        self.locate_btn.setEnabled(True)

        self.draw_graph()

    def redraw_scene(self):
        """Remove existing visuals, then re-draw."""
        self.renderer.RemoveAllViewProps()
        self.node_text_actors.clear()
        self.edge_text_actors.clear()
        self.draw_graph()

    def draw_graph(self):
        """Build and display the 3D scene based on current display mode."""
        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # Determine node values
        if self.display_mode == 'distance':
            node_values = [self.total_weight_sums.get(name, 0.0) for name in self.mer_names]
        else:  # 'edge_sum'
            node_values = [self.edge_sum_by_node.get(name, 0.0) for name in self.mer_names]

        # Normalize for color mapping
        min_val = min(node_values) if node_values else 0.0
        max_val = max(node_values) if node_values else 1.0
        if abs(max_val - min_val) < 1e-12:
            normalized_vals = [0.5]*len(node_values)
        else:
            normalized_vals = [(val - min_val)/(max_val - min_val) for val in node_values]

        # Convert normalized values to colors
        colors_rgba = []
        for val in normalized_vals:
            r, g, b = get_coolwarm_color(val)
            a = 1.0
            colors_rgba.append([r, g, b, a])

        # Highlight best-source node as yellow
        try:
            source_index = self.mer_names.index(self.mer_name)
            colors_rgba[source_index] = [1.0, 1.0, 0.0, 1.0]
        except ValueError:
            pass

        # Create vtkPolyData for node positions
        vtk_points = vtk.vtkPoints()
        vtk_colors = vtk.vtkUnsignedCharArray()
        vtk_colors.SetNumberOfComponents(3)

        for i, pos in enumerate(self.positions):
            vtk_points.InsertNextPoint(pos[0], pos[1], pos[2])
            rgba_255 = [int(c * 255) for c in colors_rgba[i][:3]]
            vtk_colors.InsertNextTypedTuple(rgba_255)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        polydata.GetPointData().SetScalars(vtk_colors)

        # Sphere glyphs
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(0.02)

        glyph_mapper = vtk.vtkGlyph3DMapper()
        glyph_mapper.SetInputData(polydata)
        glyph_mapper.SetSourceConnection(sphere_source.GetOutputPort())
        glyph_mapper.SetColorModeToDirectScalars()
        glyph_mapper.SetScalarModeToUsePointData()
        glyph_mapper.SetScaleFactor(1.0)
        glyph_mapper.ScalingOff()

        glyph_actor = vtk.vtkActor()
        glyph_actor.SetMapper(glyph_mapper)
        self.renderer.AddActor(glyph_actor)

        self.progress_bar.setValue(70)
        QApplication.processEvents()

        
        bond_color = [0.3, 0.3, 0.3]
        line_polydata = vtk.vtkPolyData()
        line_points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()

        edge_midpoints = []
        edge_weights = []

        for interaction in self.interactions:
            from_mer = interaction.from_mer
            to_mer = interaction.to_mer
            weight = interaction.weight

            if from_mer in self.mer_names and to_mer in self.mer_names:
                i = self.mer_names.index(from_mer)
                j = self.mer_names.index(to_mer)

                p1_id = line_points.InsertNextPoint(self.positions[i])
                p2_id = line_points.InsertNextPoint(self.positions[j])

                line_cell = vtk.vtkLine()
                line_cell.GetPointIds().SetId(0, p1_id)
                line_cell.GetPointIds().SetId(1, p2_id)
                lines.InsertNextCell(line_cell)

                midpoint = (self.positions[i] + self.positions[j]) / 2.0
                edge_midpoints.append(midpoint)
                edge_weights.append(weight)

        line_polydata.SetPoints(line_points)
        line_polydata.SetLines(lines)

        line_mapper = vtk.vtkPolyDataMapper()
        line_mapper.SetInputData(line_polydata)

        line_actor = vtk.vtkActor()
        line_actor.SetMapper(line_mapper)
        line_actor.GetProperty().SetColor(bond_color)
        line_actor.GetProperty().SetOpacity(0.5)
        line_actor.GetProperty().SetLineWidth(1.0)
        self.renderer.AddActor(line_actor)

        self.progress_bar.setValue(85)
        QApplication.processEvents()

        # Node labels
        for idx, name in enumerate(self.mer_names):
            if self.display_mode == 'distance':
                node_val = self.total_weight_sums.get(name, 0.0)
            else:
                node_val = self.edge_sum_by_node.get(name, 0.0)

            label_text = f"{name} ({node_val:.2f})"

            vector_text = vtk.vtkVectorText()
            vector_text.SetText(label_text)

            text_mapper = vtk.vtkPolyDataMapper()
            text_mapper.SetInputConnection(vector_text.GetOutputPort())

            text_actor = vtk.vtkFollower()
            text_actor.SetMapper(text_mapper)
            text_actor.SetScale(0.006, 0.006, 0.006)
            text_actor.SetPosition(
                self.positions[idx][0],
                self.positions[idx][1],
                self.positions[idx][2] + 0.04
            )
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)  # black text
            self.renderer.AddActor(text_actor)

            self.node_text_actors.append(text_actor)

        # Edge labels (gray)
        for midpoint, weight in zip(edge_midpoints, edge_weights):
            if math.isinf(weight):
                continue
            label_text = f"{weight:.2f}"

            vector_text = vtk.vtkVectorText()
            vector_text.SetText(label_text)

            text_mapper = vtk.vtkPolyDataMapper()
            text_mapper.SetInputConnection(vector_text.GetOutputPort())

            text_actor = vtk.vtkFollower()
            text_actor.SetMapper(text_mapper)
            text_actor.SetScale(0.005, 0.005, 0.005)
            text_actor.SetPosition(midpoint[0], midpoint[1], midpoint[2])
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.5, 0.5, 0.5)  # gray text
            self.renderer.AddActor(text_actor)

            self.edge_text_actors.append(text_actor)

        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.RemoveAllObservers()
        camera.AddObserver("ModifiedEvent", self.on_camera_modified)

        self.progress_bar.setValue(100)
        self.progress_label.setText("Processing complete!")
        QApplication.processEvents()

        # Hide loading image and progress bar after short delay, show the VTK widget
        QTimer.singleShot(500, self.show_graph)

    def show_graph(self):
        """Hide loading UI, show the 3D scene."""
        self.progress_label.hide()
        self.progress_bar.hide()
        self.loading_label.hide()

        self.main_layout.addWidget(self.vtk_widget)
        self.vtk_widget.show()

        self.interactor.Initialize()
        self.interactor.Start()

    def on_camera_modified(self, caller, event):
        """Hide or show labels based on how far the camera is from (0,0,0)."""
        camera = self.renderer.GetActiveCamera()
        cx, cy, cz = camera.GetPosition()
        dist = math.sqrt(cx**2 + cy**2 + cz**2)

        # Node labels
        if dist > self.node_label_hide_dist:
            for actor in self.node_text_actors:
                actor.VisibilityOff()
        else:
            for actor in self.node_text_actors:
                actor.VisibilityOn()

        # Edge labels
        if dist > self.edge_label_hide_dist:
            for actor in self.edge_text_actors:
                actor.VisibilityOff()
        else:
            for actor in self.edge_text_actors:
                actor.VisibilityOn()

        self.vtk_widget.GetRenderWindow().Render()

    def locate_selected_node(self):
        """
        Center and zoom on the chosen node from the combo box.
        """
        node_name = self.node_select_combobox.currentText()
        if node_name not in self.mer_names:
            return  # safety check

        idx = self.mer_names.index(node_name)
        pos = self.positions[idx]

        camera = self.renderer.GetActiveCamera()
        # Make 'pos' the new focal point
        camera.SetFocalPoint(pos[0], pos[1], pos[2])

        # Position the camera so that we're looking at 'pos' from some offset in Z
        offset_distance = 0.3  # Adjust to zoom in/out more
        camera.SetPosition(pos[0], pos[1] - 0.1, pos[2] + offset_distance)
        camera.SetViewUp(0, 1, 0)
        self.renderer.ResetCameraClippingRange()

        self.vtk_widget.GetRenderWindow().Render()

    def parse_pdb_content(self, pdb_content):
        """
        Simple parse of PDB lines to extract a single (x,y,z) position per Mer name,
        effectively using the first encountered atom of each residue as that Mer's position.
        """
        mers = {}
        lines = pdb_content.split('\n')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    mer_name = (
                        line[17:20].strip() + "-" +
                        line[22:26].strip() + "(" +
                        line[21].strip() + ")"
                    )
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if mer_name not in mers:
                        # Use the first encountered atom as the Mer center
                        mers[mer_name] = {
                            'position': (x, y, z),
                            'bond_count': 0
                        }
                except ValueError:
                    pass
        return mers, None
