import math
import numpy as np
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QProgressBar
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import QApplication
import os
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
        0.5: [0.8650, 0.8650, 0.8650]  (whiteish/light gray center)
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


class ProteinVisPyViewerPage(QWidget):
    def __init__(self, mer_name, pdb_content, interactions, total_weight_sums, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.interactions = interactions  # List of Interaction objects
        self.total_weight_sums = total_weight_sums  # Mapping of Mer names to total weights
        self.on_back = on_back

        self.mers = {}
        self.positions = None
        self.mer_names = []

        self.main_layout = QVBoxLayout(self)

        # Top bar with back button
        top_bar = QHBoxLayout()
        back_btn = QPushButton("Back")
        back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(back_btn, alignment=Qt.AlignLeft)
        top_bar.addStretch()
        self.main_layout.addLayout(top_bar)

        # Progress label and bar
        self.progress_label = QLabel("Processing graph...")
        self.progress_label.setAlignment(Qt.AlignCenter)
        self.progress_label.setFont(QFont("Arial", 14, QFont.Bold))
        self.main_layout.addWidget(self.progress_label, alignment=Qt.AlignCenter)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.main_layout.addWidget(self.progress_bar, alignment=Qt.AlignCenter)

        # Loading image instead of showing partial graph
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

        # Prepare the VTK widget (initially hidden)
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)  # White background
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()

        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        self.vtk_widget.hide()  # Initially hidden

        # Schedule parsing and drawing after a short delay
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    def parse_and_draw(self):
        # Update progress
        self.progress_bar.setValue(10)
        QApplication.processEvents()

        # Parse PDB content to get mers and positions
        self.mers, _ = self.parse_pdb_content(self.pdb_content)
        self.mer_names = list(self.mers.keys())

        self.progress_bar.setValue(30)
        QApplication.processEvents()

        if not self.mers:
            self.progress_label.setText("No Mers found.")
            self.progress_bar.setValue(100)
            return

        # Assign positions for visualization
        positions = np.array([self.mers[name]['position'] for name in self.mer_names])

        # Normalize positions for visualization
        positions -= positions.mean(axis=0)
        scale = np.max(np.linalg.norm(positions, axis=1))
        if scale > 0:
            positions /= scale

        self.positions = positions  # Store for later use

        # Assign total weight sums to Mers
        weight_sums = np.array([self.total_weight_sums.get(name, 0.0) for name in self.mer_names])

        # Normalize weight sums for coloring
        min_weight = weight_sums.min()
        max_weight = weight_sums.max()
        if max_weight != min_weight:
            normalized_weight = (weight_sums - min_weight) / (max_weight - min_weight)
        else:
            normalized_weight = np.full_like(weight_sums, 0.5)

        # Assign colors to Mers based on normalized weight
        colors_rgba = []
        for val in normalized_weight:
            r, g, b = get_coolwarm_color(val)
            a = 1.0
            colors_rgba.append([r, g, b, a])

        # If the source Mer is found, color it yellow
        try:
            source_index = self.mer_names.index(self.mer_name)
            colors_rgba[source_index] = [1.0, 1.0, 0.0, 1.0]  # RGBA for yellow
        except ValueError:
            pass

        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # Create a vtkPolyData to store points
        vtk_points = vtk.vtkPoints()
        vtk_colors = vtk.vtkUnsignedCharArray()
        vtk_colors.SetNumberOfComponents(3)  # We'll store RGB in [0..255]

        for i, pos in enumerate(self.positions):
            vtk_points.InsertNextPoint(pos[0], pos[1], pos[2])
            # Convert float colors to [0..255] range
            rgba_255 = [int(c * 255) for c in colors_rgba[i][:3]]
            vtk_colors.InsertNextTypedTuple(rgba_255)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)
        polydata.GetPointData().SetScalars(vtk_colors)

        # Use a sphere source for glyphs
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(0.02)   # Adjust as needed

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

        # Draw lines for interactions (in black)
        bond_color = [0, 0, 0]  # black
        line_polydata = vtk.vtkPolyData()
        line_points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()

        # Store midpoints + weights for labeling
        edge_midpoints = []
        edge_weights = []

        for interaction in self.interactions:
            from_mer = interaction.from_mer
            to_mer = interaction.to_mer
            weight = interaction.weight

            try:
                i = self.mer_names.index(from_mer)
                j = self.mer_names.index(to_mer)
            except ValueError:
                continue

            # Insert points for each line
            p1_id = line_points.InsertNextPoint(self.positions[i])
            p2_id = line_points.InsertNextPoint(self.positions[j])

            line_cell = vtk.vtkLine()
            line_cell.GetPointIds().SetId(0, p1_id)
            line_cell.GetPointIds().SetId(1, p2_id)
            lines.InsertNextCell(line_cell)

            # Compute midpoint for weight label
            pos1 = self.positions[i]
            pos2 = self.positions[j]
            midpoint = (pos1 + pos2) / 2.0
            edge_midpoints.append(midpoint)
            edge_weights.append(weight)

        line_polydata.SetPoints(line_points)
        line_polydata.SetLines(lines)

        # Mapper/actor for lines
        line_mapper = vtk.vtkPolyDataMapper()
        line_mapper.SetInputData(line_polydata)
        line_actor = vtk.vtkActor()
        line_actor.SetMapper(line_mapper)
        line_actor.GetProperty().SetColor(bond_color)
        line_actor.GetProperty().SetOpacity(0.5)
        line_actor.GetProperty().SetLineWidth(1.0)
        self.renderer.AddActor(line_actor)

        self.progress_bar.setValue(80)
        QApplication.processEvents()

        # Add 3D labels for each Mer (mer_name + total_weight)
        for idx, name in enumerate(self.mer_names):
            total_weight = self.total_weight_sums.get(name, 0.0)
            label_text = f"{name} ({total_weight:.2f})"

            # Use vtkVectorText + vtkFollower
            vector_text = vtk.vtkVectorText()
            vector_text.SetText(label_text)

            text_mapper = vtk.vtkPolyDataMapper()
            text_mapper.SetInputConnection(vector_text.GetOutputPort())

            text_actor = vtk.vtkFollower()
            text_actor.SetMapper(text_mapper)
            text_actor.SetScale(0.01, 0.01, 0.01)  # Adjust
            # Slightly offset above the sphere
            text_actor.SetPosition(
                positions[idx][0],
                positions[idx][1],
                positions[idx][2] + 0.04
            )
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)  # Black text
            self.renderer.AddActor(text_actor)

        # Add 3D labels for each edge's weight
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
            text_actor.SetScale(0.006, 0.006, 0.006)  # Slightly smaller
            text_actor.SetPosition(midpoint[0], midpoint[1], midpoint[2])
            text_actor.SetCamera(self.renderer.GetActiveCamera())
            text_actor.GetProperty().SetColor(0.0, 0.0, 0.0)
            self.renderer.AddActor(text_actor)

        self.progress_bar.setValue(90)
        QApplication.processEvents()

        # Adjust camera
        self.renderer.ResetCamera()

        # Finalize progress
        self.progress_bar.setValue(100)
        self.progress_label.setText("Processing graph...")
        QApplication.processEvents()

        # Hide loading image and progress bar after short delay and show the graph
        QTimer.singleShot(500, self.show_graph)

    def show_graph(self):
        # Remove loading label and progress elements
        self.progress_label.hide()
        self.progress_bar.hide()
        self.loading_label.hide()

        # Add the VTK widget to the layout and show it
        self.main_layout.addWidget(self.vtk_widget)
        self.vtk_widget.show()

        # Initialize/Start the interactor
        self.interactor.Initialize()
        self.interactor.Start()

    def parse_pdb_content(self, pdb_content):
        mers = {}
        lines = pdb_content.split('\n')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    # Extract Mer name and position
                    mer_name = line[17:20].strip() + "-" + line[22:26].strip() + "(" + line[21].strip() + ")"
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    mers[mer_name] = {
                        'position': (x, y, z),
                        'bond_count': 0  # Initialize bond count
                    }
                except ValueError:
                    pass
        return mers, None