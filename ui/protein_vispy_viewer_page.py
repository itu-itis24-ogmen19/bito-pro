# /Users/m3_air/24-25 fall/Bitirme/bito-pro-main/ui/protein_vispy_viewer_page.py

import math
import numpy as np
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QProgressBar
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QPixmap
from vispy import scene
from vispy.color import get_colormap
from vispy.scene import visuals
from PySide6.QtWidgets import QApplication
import os

from ui.resource_locate import resource_path


class ProteinVisPyViewerPage(QWidget):
    def __init__(self, mer_name, pdb_content, interactions, total_weight_sums, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.interactions = interactions  # List of Interaction objects
        self.total_weight_sums = total_weight_sums  # dict of mer_name: total_weight_sum
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
            # Scale pixmap to a reasonable size if it's large (e.g., 200x200)
            scaled_pixmap = pixmap.scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            self.loading_label.setPixmap(scaled_pixmap)
        else:
            self.loading_label.setText("Loading...")
        self.loading_label.setAlignment(Qt.AlignCenter)
        self.main_layout.addWidget(self.loading_label, alignment=Qt.AlignCenter)

        # VisPy Canvas (initially hidden)
        self.canvas = scene.SceneCanvas(keys='interactive', show=False)
        self.canvas.bgcolor = 'white'
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable'

        # Schedule parsing and drawing after a short delay
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    def parse_and_draw(self):
        # Update progress to show we are starting
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

        # Assign total weight sums to Mers
        weight_sums = np.array([self.total_weight_sums.get(name, 0.0) for name in self.mer_names])

        # Normalize weight sums for coloring
        min_weight = weight_sums.min()
        max_weight = weight_sums.max()
        if max_weight != min_weight:
            normalized_weight = (weight_sums - min_weight) / (max_weight - min_weight)
        else:
            normalized_weight = np.full_like(weight_sums, 0.5)

        # Map normalized weight sums to colors (low weight → red, high weight → blue)
        cmap = get_colormap('coolwarm')  # Reversed 'coolwarm' colormap
        colors = cmap.map(normalized_weight)

        # Set the color of the source mer to yellow
        try:
            source_index = self.mer_names.index(self.mer_name)
            colors[source_index] = [1.0, 1.0, 0.0, 1.0]  # RGBA for yellow
        except ValueError:
            pass  # If source mer not found, skip

        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # Create scatter plot for Mers with weight sum-based colors
        scatter = scene.visuals.Markers(parent=self.view.scene)
        scatter.set_data(positions, face_color=colors, size=10)

        # Add labels for Mers
        for idx, name in enumerate(self.mer_names):
            label_text = name
            text = scene.visuals.Text(
                text=label_text,
                pos=positions[idx],
                color='black',
                font_size=10,  # Adjusted for readability
                anchor_x='center',
                anchor_y='center',
                parent=self.view.scene
            )
            # Offset the label slightly above the Mer
            text.pos = (positions[idx][0], positions[idx][1], positions[idx][2] + 0.05)

        self.progress_bar.setValue(70)
        QApplication.processEvents()

        # Collect all weight values for bond coloring
        weights = [interaction.weight for interaction in self.interactions]
        weights = np.array(weights)

        if len(weights) > 0:
            # Normalize weights
            min_weight_val = weights.min()
            max_weight_val = weights.max()
            if max_weight_val != min_weight_val:
                normalized_weight_val = (weights - min_weight_val) / (max_weight_val - min_weight_val)
            else:
                normalized_weight_val = np.full_like(weights, 0.5)

            # Choose a colormap for bonds
            cmap_bonds = get_colormap('autumn')  # 'plasma' for high contrast

            # Map normalized weight to colors
            colors_bonds = ['black'] * len(normalized_weight_val)

            # Create lines for bonds and add weight labels
            for idx, interaction in enumerate(self.interactions):
                from_mer = interaction.from_mer  # string
                to_mer = interaction.to_mer      # string
                weight = interaction.weight  # float

                # Get indices of the Mers
                try:
                    i = self.mer_names.index(from_mer)
                    j = self.mer_names.index(to_mer)
                except ValueError:
                    continue  # Skip if Mer not found

                # Get positions
                pos1 = positions[i]
                pos2 = positions[j]

                # Get color for this bond
                color = colors_bonds[idx]

                # Create a line between the two Mers with weight-based color
                line = scene.visuals.Line(
                    pos=np.array([pos1, pos2]),
                    color=color,
                    width=0.5,
                    connect='segments',
                    parent=self.view.scene
                )

                # Compute midpoint for weight label
                midpoint = (pos1 + pos2) / 2.0

                # Create a text label for weight
                weight_text = f"{weight:.2f}"
                text = scene.visuals.Text(
                    text=weight_text,
                    pos=midpoint,
                    color='black',  # Contrast with bond color
                    font_size=4,    # Adjusted for readability
                    anchor_x='center',
                    anchor_y='center',
                    parent=self.view.scene
                )
                # Offset the label slightly above the bond
                text.pos = (midpoint[0], midpoint[1], midpoint[2] + 0.02)

        self.progress_bar.setValue(90)
        QApplication.processEvents()

        # Zoom to fit
        self.view.camera.set_range()

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

        # Now add the canvas to the layout and show it
        self.main_layout.addWidget(self.canvas.native)
        self.canvas.show()

    def parse_pdb_content(self, pdb_content):
        mers = {}
        lines = pdb_content.split('\n')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    # Assuming mer_name is formatted as "ResidueName-ResidueNumber(ChainID)"
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
