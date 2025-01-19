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
    QMessageBox
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
        0.5: [0.8650, 0.8650, 0.8650]  (white/light gray)
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
    ProteinViewerPage now supports two display modes:
    1) Global Distance Mode: Each node's value = distance from 'best source' Mer.
    2) Local Edge-Sum Mode: Each node's value = sum of the weights of edges connected to it.
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

        # Main layout
        self.main_layout = QVBoxLayout(self)

        # --- Top Bar with Back Button, Mode Switching, Info Button ---
        top_bar = QHBoxLayout()

        back_btn = QPushButton("Back")
        back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(back_btn, alignment=Qt.AlignLeft)

        self.mode_label = QLabel("Current Mode: Global Distance")  # Default text
        self.mode_label.setFont(QFont("Arial", 11, QFont.Bold))
        top_bar.addWidget(self.mode_label, alignment=Qt.AlignCenter)

        # Info button to explain each mode
        info_btn = QPushButton("Info")
        info_btn.setToolTip("Learn about the two display modes and how this 3D visualization is computed.")
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

        top_bar.addStretch()
        self.main_layout.addLayout(top_bar)

        # --- Progress Label and Bar ---
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

        # --- VTK Widget (initially hidden) ---
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)  # White background
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()

        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        self.vtk_widget.hide()  # Hidden until data is drawn

        # Schedule parsing and drawing after a short delay
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    def show_info(self):
        """
        Show a popup describing the two modes and the overall 3D visualization approach.
        """
        info_text = (
            "This 3D viewer supports two different ways of labeling and coloring each Mer (node):\n\n"
            "1) Global Distance Mode:\n"
            "   Each node's value is the total distance (or total weight) from the 'best source' Mer.\n"
            "   This reveals how centrally located or 'close' a Mer is to the chosen best source.\n\n"
            "2) Local Edge-Sum Mode:\n"
            "   Each node's value is the sum of the weights of edges connected to it.\n"
            "   This highlights how many strong connections a Mer has in its immediate neighborhood.\n\n"
            "In both modes, edges are shown with their individual weights, and node colors are\n"
            "scaled using a cool-to-warm color gradient. The node of the best source Mer is displayed in yellow.\n\n"
            "Use the buttons on the top-right to switch between these two display modes."
        )
        QMessageBox.information(self, "Display Modes", info_text)

    def set_distance_mode(self):
        """Switch to the mode displaying each node's distance from the 'best source' Mer."""
        self.display_mode = 'distance'
        self.mode_label.setText("Current Mode: Global Distance")
        self.redraw_scene()

    def set_edge_sum_mode(self):
        """
        Switch to the mode displaying each node's sum of connected edge weights
        (a local measure of node connectivity).
        """
        self.display_mode = 'edge_sum'
        self.mode_label.setText("Current Mode: Local Edge-Sum")
        self.redraw_scene()

    def parse_and_draw(self):
        """Initial data parsing + scene drawing, called once after a short delay."""
        self.progress_bar.setValue(10)
        QApplication.processEvents()

        # 1. Parse PDB content to get Mers and positions
        self.mers, _ = self.parse_pdb_content(self.pdb_content)
        self.mer_names = list(self.mers.keys())

        self.progress_bar.setValue(30)
        QApplication.processEvents()

        if not self.mers:
            self.progress_label.setText("No Mers found.")
            self.progress_bar.setValue(100)
            return

        # 2. Assign normalized positions for visualization
        positions = np.array([self.mers[name]['position'] for name in self.mer_names])
        positions -= positions.mean(axis=0)
        scale = np.max(np.linalg.norm(positions, axis=1))
        if scale > 0:
            positions /= scale
        self.positions = positions

        # 3. Prepare "edge sum" data: sum of all edges for each node
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

        # 4. Draw the scene in the current mode (default: distance)
        self.draw_graph()

    def redraw_scene(self):
        """
        Clear the existing 3D scene and re-draw it using the current display mode.
        """
        self.renderer.RemoveAllViewProps()
        self.draw_graph()

    def draw_graph(self):
        """Build and display the 3D scene according to self.display_mode."""
        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # 1. Determine node values
        if self.display_mode == 'distance':
            # Use self.total_weight_sums
            node_values = [self.total_weight_sums.get(name, 0.0) for name in self.mer_names]
        else:  # 'edge_sum'
            node_values = [self.edge_sum_by_node.get(name, 0.0) for name in self.mer_names]

        # 2. Normalize these values for color mapping
        min_val = min(node_values) if node_values else 0.0
        max_val = max(node_values) if node_values else 1.0
        if abs(max_val - min_val) < 1e-12:
            normalized_vals = [0.5]*len(node_values)
        else:
            normalized_vals = [(val - min_val)/(max_val - min_val) for val in node_values]

        # 3. Convert normalized values to cool-to-warm colors
        colors_rgba = []
        for val in normalized_vals:
            r, g, b = get_coolwarm_color(val)
            a = 1.0
            colors_rgba.append([r, g, b, a])

        # If the source Mer is in this list, color it yellow to highlight
        try:
            source_index = self.mer_names.index(self.mer_name)
            colors_rgba[source_index] = [1.0, 1.0, 0.0, 1.0]  # RGBA for yellow
        except ValueError:
            pass

        # 4. Create a vtkPolyData for the Mer positions
        vtk_points = vtk.vtkPoints()
        vtk_colors = vtk.vtkUnsignedCharArray()
        vtk_colors.SetNumberOfComponents(3)  # We'll store RGB in [0..255]

        for i, pos in enumerate(self.positions):
            vtk_points.InsertNextPoint(pos[0], pos[1], pos[2])
            rgba_255 = [int(c * 255) for c in colors_rgba[i][:3]]
            vtk_colors.InsertNextTypedTuple(rgba_255)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        polydata.GetPointData().SetScalars(vtk_colors)

        # 5. Sphere glyphs for each node
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(0.02)  # Adjust as needed

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

        # 6. Draw lines for interactions
        self.progress_bar.setValue(70)
        QApplication.processEvents()

        bond_color = [0, 0, 0]  # black
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

        # 7. Add 3D labels for each node
        self.progress_bar.setValue(85)
        QApplication.processEvents()

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
            text_actor.SetScale(0.01, 0.01, 0.01)  # Adjust as needed
            text_actor.SetPosition(
                self.positions[idx][0],
                self.positions[idx][1],
                self.positions[idx][2] + 0.04
            )
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)  # Black text
            self.renderer.AddActor(text_actor)

        # 8. Add labels for each edge
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
            text_actor.SetScale(0.006, 0.006, 0.006)
            text_actor.SetPosition(midpoint[0], midpoint[1], midpoint[2])
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)
            self.renderer.AddActor(text_actor)

        # Finish up
        self.progress_bar.setValue(95)
        QApplication.processEvents()

        self.renderer.ResetCamera()

        self.progress_bar.setValue(100)
        self.progress_label.setText("Processing complete!")
        QApplication.processEvents()

        # Hide loading image and progress bar after short delay, show the VTK widget
        QTimer.singleShot(500, self.show_graph)

    def show_graph(self):
        """Hide loading UI, show the 3D graph."""
        self.progress_label.hide()
        self.progress_bar.hide()
        self.loading_label.hide()

        self.main_layout.addWidget(self.vtk_widget)
        self.vtk_widget.show()

        self.interactor.Initialize()
        self.interactor.Start()

    def parse_pdb_content(self, pdb_content):
        """
        Simple parse of PDB lines to extract a single (x,y,z) position
        per Mer name. This lumps all atoms for a Mer into one 'center',
        based on the first encountered coordinate lines in the PDB.
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
                    mers[mer_name] = {
                        'position': (x, y, z),
                        'bond_count': 0
                    }
                except ValueError:
                    pass
        return mers, None
