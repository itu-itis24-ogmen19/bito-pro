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
from PySide6.QtCore import Qt, QTimer, QSize
from PySide6.QtGui import QFont, QPixmap, QIcon
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
    ProteinViewerPage now shows in "Show Protein Info":
      - The chosen PDB file name
      - The center chosen Mer (self.mer_name)
      - Total Mers
      - Total Edges
    """
    def __init__(
        self,
        mer_name,
        pdb_content,
        interactions,
        total_weight_sums,
        on_back,
        pdb_file_name=None  # <--- NEW optional parameter
    ):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.interactions = interactions
        self.total_weight_sums = total_weight_sums
        self.on_back = on_back

        # Store the file name for "Show Protein Info"
        self.pdb_file_name = pdb_file_name or "Unknown File"

        # Data structures
        self.mers = {}
        self.positions = None
        self.mer_names = []
        self.edge_sum_by_node = {}

        # Which display mode? 'distance' or 'edge_sum'
        self.display_mode = 'distance'

        # Label actors for dynamic hiding
        self.node_text_actors = []
        self.edge_text_actors = []

        # Distance thresholds for hiding labels
        self.node_label_hide_dist = 3.0
        self.edge_label_hide_dist = 3.0

        # --- Main layout ---
        self.main_layout = QVBoxLayout(self)

        # =========================
        # Top Bar
        # =========================
        top_bar = QHBoxLayout()

        back_btn = QPushButton("Back")
        back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(back_btn, alignment=Qt.AlignLeft)

        # Make the mode label smaller (9pt instead of 11)
        self.mode_label = QLabel("Current Mode: Global Distance")
        self.mode_label.setFont(QFont("Arial", 9, QFont.Bold))
        top_bar.addWidget(self.mode_label, alignment=Qt.AlignCenter)

        # Info button with an icon instead of text:
        info_btn = QPushButton()
        info_btn.setToolTip("Information about modes and features.")
        info_icon_path = resource_path(os.path.join("assets", "images", "info_icon.png"))
        if os.path.exists(info_icon_path):
            info_btn.setIcon(QIcon(info_icon_path))
            info_btn.setIconSize(QSize(24, 24))
        else:
            info_btn.setText("Info")
        info_btn.clicked.connect(self.show_info)
        top_bar.addWidget(info_btn, alignment=Qt.AlignRight)

        self.distance_btn = QPushButton("Global Distance Mode")
        self.distance_btn.clicked.connect(self.set_distance_mode)
        self.distance_btn.setToolTip("Node values = total distance from the best source Mer.")
        top_bar.addWidget(self.distance_btn, alignment=Qt.AlignRight)

        self.edge_sum_btn = QPushButton("Local Edge-Sum Mode")
        self.edge_sum_btn.clicked.connect(self.set_edge_sum_mode)
        self.edge_sum_btn.setToolTip("Node values = sum of connected edges' weights.")
        top_bar.addWidget(self.edge_sum_btn, alignment=Qt.AlignRight)

        # **Node selection**: For full-length names, we set size policy
        self.node_select_combobox = QComboBox()
        self.node_select_combobox.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.node_select_combobox.setMinimumWidth(250)
        self.node_select_combobox.setEnabled(False)
        top_bar.addWidget(self.node_select_combobox)

        self.locate_btn = QPushButton("Locate Node")
        self.locate_btn.setEnabled(False)
        self.locate_btn.clicked.connect(self.locate_selected_node)
        top_bar.addWidget(self.locate_btn, alignment=Qt.AlignRight)

        # **Show Node Info** button
        self.info_node_btn = QPushButton("Show Node Info")
        self.info_node_btn.setEnabled(False)
        self.info_node_btn.clicked.connect(self.show_node_info)
        top_bar.addWidget(self.info_node_btn, alignment=Qt.AlignRight)

        # **Show Protein Info** button
        self.protein_info_btn = QPushButton("Show Protein Info")
        self.protein_info_btn.setEnabled(False)
        self.protein_info_btn.clicked.connect(self.show_protein_info)
        top_bar.addWidget(self.protein_info_btn, alignment=Qt.AlignRight)

        top_bar.addStretch()
        self.main_layout.addLayout(top_bar)

        # =========================
        # Progress Label + Bar
        # =========================
        self.progress_label = QLabel("Processing graph...")
        self.progress_label.setAlignment(Qt.AlignCenter)
        self.progress_label.setFont(QFont("Arial", 14, QFont.Bold))
        self.main_layout.addWidget(self.progress_label, alignment=Qt.AlignCenter)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(self.progress_bar, alignment=Qt.AlignCenter)

        # =========================
        # Loading Image
        # =========================
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

        # =========================
        # VTK Setup
        # =========================
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()

        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        self.vtk_widget.hide()  # hidden until data is drawn

        # Kick off the process
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    # ----------------------------------------------------------------
    # Info / Help
    # ----------------------------------------------------------------
    def show_info(self):
        """
        Show a popup describing the two modes, distance-based hiding, node-locating, etc.
        """
        info_text = (
            "This 3D viewer supports two different ways of labeling each Mer (node):\n\n"
            "1) Global Distance Mode:\n"
            "   Each node's value is the total distance (or weight) from the best source Mer.\n\n"
            "2) Local Edge-Sum Mode:\n"
            "   Each node's value is the sum of the weights of edges connected to it.\n\n"
            "DISTANCE-BASED LABEL HIDING:\n"
            "   Labels automatically hide when you zoom out, and show when you zoom in.\n\n"
            "LOCATE NODE:\n"
            "   Pick a node from the drop-down, then click 'Locate Node' to center the camera.\n\n"
            "SHOW NODE INFO:\n"
            "   Displays the selected node's current value and all connections (edges).\n\n"
            "SHOW PROTEIN INFO:\n"
            "   Shows how many nodes and edges exist in this loaded protein data."
        )
        QMessageBox.information(self, "Display Modes & Features", info_text)

    # ----------------------------------------------------------------
    # Set Display Modes
    # ----------------------------------------------------------------
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

    # ----------------------------------------------------------------
    # Parse PDB / Construct Graph
    # ----------------------------------------------------------------
    def parse_and_draw(self):
        """Load PDB, build node data structures, then draw the scene."""
        self.progress_bar.setValue(10)
        QApplication.processEvents()

        self.mers, _ = self.parse_pdb_content(self.pdb_content)
        self.mer_names = list(self.mers.keys())

        self.progress_bar.setValue(30)
        QApplication.processEvents()

        if not self.mer_names:
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

        # Build edge sums
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

        # Populate combo box
        self.node_select_combobox.clear()
        for n in self.mer_names:
            self.node_select_combobox.addItem(n)
        self.node_select_combobox.setEnabled(True)
        self.locate_btn.setEnabled(True)
        self.info_node_btn.setEnabled(True)
        self.protein_info_btn.setEnabled(True)

        self.draw_graph()

    def redraw_scene(self):
        """Clear and re-draw the 3D scene according to current mode."""
        self.renderer.RemoveAllViewProps()
        self.node_text_actors.clear()
        self.edge_text_actors.clear()
        self.draw_graph()

    # ----------------------------------------------------------------
    # Drawing the 3D Scene
    # ----------------------------------------------------------------
    def draw_graph(self):
        """Create vtk actors for nodes/edges/labels based on the selected mode."""
        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # Node values
        if self.display_mode == 'distance':
            node_values = [self.total_weight_sums.get(name, 0.0) for name in self.mer_names]
        else:
            node_values = [self.edge_sum_by_node.get(name, 0.0) for name in self.mer_names]

        min_val = min(node_values) if node_values else 0.0
        max_val = max(node_values) if node_values else 1.0
        if abs(max_val - min_val) < 1e-12:
            normalized_vals = [0.5]*len(node_values)
        else:
            normalized_vals = [(val - min_val)/(max_val - min_val) for val in node_values]

        # Convert to color
        colors_rgba = []
        for val in normalized_vals:
            r, g, b = get_coolwarm_color(val)
            a = 1.0
            colors_rgba.append([r, g, b, a])

        # Highlight best-source node
        try:
            source_index = self.mer_names.index(self.mer_name)
            colors_rgba[source_index] = [1.0, 1.0, 0.0, 1.0]
        except ValueError:
            pass

        # Build vtkPolyData for positions
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

        # Node glyphs
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

        # Edges
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
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)
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
            text_actor.GetProperty().SetColor(0.5, 0.5, 0.5)
            self.renderer.AddActor(text_actor)

            self.edge_text_actors.append(text_actor)

        # Reset camera + attach observer for label hiding
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.RemoveAllObservers()
        camera.AddObserver("ModifiedEvent", self.on_camera_modified)

        self.progress_bar.setValue(100)
        self.progress_label.setText("Processing complete!")
        QApplication.processEvents()

        QTimer.singleShot(500, self.show_graph)

    def show_graph(self):
        """Hide the loading UI and show the VTK widget."""
        self.progress_label.hide()
        self.progress_bar.hide()
        self.loading_label.hide()

        self.main_layout.addWidget(self.vtk_widget)
        self.vtk_widget.show()

        self.interactor.Initialize()
        self.interactor.Start()

    # ----------------------------------------------------------------
    # Label Hiding Callback
    # ----------------------------------------------------------------
    def on_camera_modified(self, caller, event):
        """Hide or show labels based on camera distance from (0,0,0)."""
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

    # ----------------------------------------------------------------
    # Locate Node
    # ----------------------------------------------------------------
    def locate_selected_node(self):
        """
        Center and zoom the camera on the selected node in the combo box.
        """
        node_name = self.node_select_combobox.currentText()
        if node_name not in self.mer_names:
            return  # safety check

        idx = self.mer_names.index(node_name)
        pos = self.positions[idx]

        camera = self.renderer.GetActiveCamera()
        # Make 'pos' the new focal point
        camera.SetFocalPoint(pos[0], pos[1], pos[2])

        # Some offset to zoom in from above
        offset_distance = 0.3
        camera.SetPosition(pos[0], pos[1] - 0.1, pos[2] + offset_distance)
        camera.SetViewUp(0, 1, 0)
        self.renderer.ResetCameraClippingRange()

        self.vtk_widget.GetRenderWindow().Render()

    # ----------------------------------------------------------------
    # Show Node Info
    # ----------------------------------------------------------------
    def show_node_info(self):
        """
        Displays a popup with the selected node's:
          - Current mode value
          - All connections (edges) with their weights
        """
        node_name = self.node_select_combobox.currentText()
        if node_name not in self.mer_names:
            return

        # Gather node's current mode value
        if self.display_mode == 'distance':
            node_value = self.total_weight_sums.get(node_name, 0.0)
            mode_str = "Global Distance"
        else:
            node_value = self.edge_sum_by_node.get(node_name, 0.0)
            mode_str = "Local Edge-Sum"

        # Gather all connections
        connections = []
        for inter in self.interactions:
            if inter.from_mer == node_name:
                connections.append((inter.to_mer, inter.weight))
            elif inter.to_mer == node_name:
                connections.append((inter.from_mer, inter.weight))

        # Build a string
        info_str = []
        info_str.append(f"Node: {node_name}")
        info_str.append(f"Current Mode: {mode_str}")
        info_str.append(f"Node Value: {node_value:.3f}\n")

        if connections:
            info_str.append("Connections:")
            for (other, w) in connections:
                info_str.append(f"  • {node_name} <--> {other}, weight={w:.3f}")
        else:
            info_str.append("(No connections found)")

        full_text = "\n".join(info_str)

        QMessageBox.information(self, "Node Info", full_text)

    # ----------------------------------------------------------------
    # Show Protein Info
    # ----------------------------------------------------------------
    def show_protein_info(self):
        """
        Displays:
          - The name of the PDB file
          - The best (center) chosen Mer
          - How many nodes (mers) exist
          - How many edges (connections)
        """
        node_count = len(self.mer_names)
        edge_count = len(self.interactions)

        info_text = (
            f"Protein Info:\n\n"
            f"Name of the chosen PDB file: {self.pdb_file_name}\n"
            f"Center chosen Mer: {self.mer_name}\n"
            f"Total Mers (Nodes): {node_count}\n"
            f"Total Edges (Connections): {edge_count}\n"
        )
        QMessageBox.information(self, "Protein Info", info_text)

    # ----------------------------------------------------------------
    # PDB Parsing
    # ----------------------------------------------------------------
    def parse_pdb_content(self, pdb_content):
        """
        Return a dict of MerName -> {position: (x,y,z), bond_count: 0}, etc.
        We use the first encountered line for each mer as the position.
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
                    # If we haven't yet assigned a position for this mer_name, do so
                    if mer_name not in mers:
                        mers[mer_name] = {
                            'position': (x, y, z),
                            'bond_count': 0
                        }
                except ValueError:
                    pass
        return mers, None
