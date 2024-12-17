# ui/protein_vispy_viewer_page.py
import math
import numpy as np
from PySide6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QHBoxLayout
from PySide6.QtCore import Qt
from vispy import scene
from vispy.color import get_colormap
from vispy.scene import visuals

class ProteinVisPyViewerPage(QWidget):
    def __init__(self, mer_name, pdb_content, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.on_back = on_back

        self.atoms = []
        self.positions = None
        self.bonds = []

        layout = QVBoxLayout(self)

        # Top bar with back button
        top_bar = QHBoxLayout()
        back_btn = QPushButton("Back")
        back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(back_btn, alignment=Qt.AlignLeft)
        top_bar.addStretch()
        layout.addLayout(top_bar)

        # VisPy Canvas
        self.canvas = scene.SceneCanvas(keys='interactive', show=True)
        self.canvas.bgcolor = 'white'  # White background for text visibility
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = 'turntable'

        # Embed the VisPy canvas into the PySide6 widget
        layout.addWidget(self.canvas.native)
        self.setLayout(layout)

        self.parse_and_draw()

    def parse_and_draw(self):
        # Parse PDB content to get atoms and bonds
        atoms, bonds = self.parse_pdb_atoms_and_bonds(self.pdb_content)
        self.atoms = atoms
        self.bonds = bonds

        if not atoms:
            return

        # Extract positions and temp factors as NumPy arrays
        self.positions = np.array([(atom['x'], atom['y'], atom['z']) for atom in atoms])
        temp_factors = np.array([atom['temp_factor'] for atom in atoms])

        # Normalize temp factors
        min_tf = temp_factors.min()
        max_tf = temp_factors.max()
        if max_tf != min_tf:
            normalized_tf = (temp_factors - min_tf) / (max_tf - min_tf)
        else:
            normalized_tf = np.full_like(temp_factors, 0.5)

        # Map colors from blue (low) to red (high)
        cmap = get_colormap('coolwarm')
        colors = cmap.map(normalized_tf)

        # Create scatter plot for atoms (larger nodes)
        scatter = scene.visuals.Markers(parent=self.view.scene)
        scatter.set_data(self.positions, face_color=colors, size=10)

        # Add text labels for nodes with much larger font size
        # Show full name and temp factor
        for i, atom in enumerate(atoms):
            label_text = f"{atom['full_name']} (TF={atom['temp_factor']:.2f})"
            text = scene.visuals.Text(
                text=label_text,
                pos=self.positions[i],
                color='black',
                font_size=120,   # Increased from 12 to 120
                anchor_x='left',
                anchor_y='bottom',
                parent=self.view.scene
            )
            # Offset so labels don't overlap nodes
            text.pos = (self.positions[i][0] + 0.3, self.positions[i][1] + 0.3, self.positions[i][2])

        # Create lines for bonds and show affinity as 1.0 / distance
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

            affinity = 1.0 / dist if dist != 0 else 0.0

            edge_label = scene.visuals.Text(
                text=f"Affinity={affinity:.2f}",
                pos=midpoint,
                color='blue',
                font_size=100,  # Increased from 10 to 100
                anchor_x='center',
                anchor_y='center',
                parent=self.view.scene
            )
            # Slight vertical offset for clarity
            edge_label.pos = (midpoint[0], midpoint[1] + 0.2, midpoint[2])

        # Zoom to fit
        self.view.camera.set_range()

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
