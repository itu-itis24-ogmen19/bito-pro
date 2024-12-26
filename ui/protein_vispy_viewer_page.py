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
    def __init__(self, mer_name, pdb_content, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.on_back = on_back

        self.atoms = []
        self.positions = None
        self.bonds = []

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

        # We will create the canvas but not add it yet.
        self.canvas = scene.SceneCanvas(keys='interactive', show=False)  
        self.canvas.bgcolor = 'white'
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable'

        # Parse and draw after a slight delay so UI can show initial progress
        QTimer.singleShot(100, self.parse_and_draw)

        self.setLayout(self.main_layout)

    def parse_and_draw(self):
        # Update progress to show we are starting
        self.progress_bar.setValue(10)
        QApplication.processEvents()

        # Parse PDB content to get atoms and bonds
        atoms, bonds = self.parse_pdb_atoms_and_bonds(self.pdb_content)
        self.atoms = atoms
        self.bonds = bonds

        self.progress_bar.setValue(30)
        QApplication.processEvents()

        if not atoms:
            self.progress_label.setText("No atoms found.")
            self.progress_bar.setValue(100)
            return

        # Extract positions and temp factors
        self.positions = np.array([(atom['x'], atom['y'], atom['z']) for atom in atoms])
        temp_factors = np.array([atom['temp_factor'] for atom in atoms])

        # Normalize temp factors
        min_tf = temp_factors.min()
        max_tf = temp_factors.max()
        if max_tf != min_tf:
            normalized_tf = (temp_factors - min_tf) / (max_tf - min_tf)
        else:
            normalized_tf = np.full_like(temp_factors, 0.5)

        self.progress_bar.setValue(50)
        QApplication.processEvents()

        # Map colors
        cmap = get_colormap('coolwarm')
        colors = cmap.map(normalized_tf)

        # Create scatter plot
        scatter = scene.visuals.Markers(parent=self.view.scene)
        scatter.set_data(self.positions, face_color=colors, size=10)

        self.progress_bar.setValue(70)
        QApplication.processEvents()

        # Add text labels
        for i, atom in enumerate(atoms):
            label_text = f"{atom['full_name']} (TF={atom['temp_factor']:.2f})"
            text = scene.visuals.Text(
                text=label_text,
                pos=self.positions[i],
                color='black',
                font_size=120,
                anchor_x='left',
                anchor_y='bottom',
                parent=self.view.scene
            )
            text.pos = (self.positions[i][0] + 0.3, self.positions[i][1] + 0.3, self.positions[i][2])

        # Create lines for bonds and add labels
        for bond in bonds:
            i, j = bond
            pos = np.array([self.positions[i], self.positions[j]])
            line = scene.visuals.Line(pos, color='black', width=1, parent=self.view.scene)

            # Compute midpoint and distance
            midpoint = (pos[0] + pos[1]) / 2.0
            dx = pos[1][0] - pos[0][0]
            dy = pos[1][1] - pos[0][1]
            dz = pos[1][2] - pos[0][2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            affinity = 1.0/dist if dist != 0 else 0.0

            edge_label = scene.visuals.Text(
                text=f"Affinity={affinity:.2f}",
                pos=midpoint,
                color='blue',
                font_size=100,
                anchor_x='center',
                anchor_y='center',
                parent=self.view.scene
            )
            edge_label.pos = (midpoint[0], midpoint[1] + 0.2, midpoint[2])

        self.progress_bar.setValue(90)
        QApplication.processEvents()

        # Zoom to fit
        self.view.camera.set_range()

        # Now UI is ready, set progress to 100%
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

    def parse_pdb_atoms_and_bonds(self, pdb_content, bond_threshold=4.5):
        atoms = []
        lines = pdb_content.split('\n')
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    temp_factor = float(line[60:66])
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    residue_number = line[22:26].strip()
                    atom_name = line[12:16].strip()
                    full_name = f"{atom_name}-{residue_name}({chain_id}){residue_number}"
                    atoms.append({
                        'x': x,
                        'y': y,
                        'z': z,
                        'temp_factor': temp_factor,
                        'full_name': full_name
                    })
                except ValueError:
                    pass

        # Calculate bonds based on distance threshold
        bonds = []
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                dx = atoms[i]['x'] - atoms[j]['x']
                dy = atoms[i]['y'] - atoms[j]['y']
                dz = atoms[i]['z'] - atoms[j]['z']
                dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                if dist < bond_threshold:
                    bonds.append((i, j))
        return atoms, bonds