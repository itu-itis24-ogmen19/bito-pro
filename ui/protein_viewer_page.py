# ui/protein_viewer_page.py
import math
from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsTextItem, QPushButton
from PySide6.QtCore import Qt, QRectF
from PySide6.QtGui import QBrush, QColor, QPen

class ProteinViewerPage(QWidget):
    def __init__(self, mer_name, pdb_content, on_back):
        super().__init__()
        self.mer_name = mer_name
        self.pdb_content = pdb_content
        self.on_back = on_back

        layout = QVBoxLayout(self)
        
        title = QLabel(f"Mers Graph - {self.mer_name}")
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        self.view = QGraphicsView()
        self.scene = QGraphicsScene(self)
        self.view.setScene(self.scene)
        layout.addWidget(self.view)

        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        layout.addWidget(self.back_btn)

        self.parse_and_draw()

    def parse_and_draw(self):
        mer_t_factors = self.parse_pdb_content(self.pdb_content)
        # Draw nodes in a circular layout
        self.draw_nodes(mer_t_factors)

    def parse_pdb_content(self, pdb_content):
        lines = pdb_content.split('\n')
        mer_t_factors = {}
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resName = line[17:20].strip()
                chainID = line[21:22].strip()
                resSeq = line[22:26].strip()
                tFactorStr = line[60:66].strip()
                try:
                    tFactor = float(tFactorStr)
                except:
                    tFactor = 0.0
                merIdentifier = f"{resName}-{resSeq}({chainID})"
                mer_t_factors.setdefault(merIdentifier, []).append(tFactor)

        # Average
        avg_factors = {}
        for k,v in mer_t_factors.items():
            if v:
                avg_f = sum(v)/len(v)
            else:
                avg_f = 0.0
            avg_factors[k] = avg_f

        return avg_factors

    def draw_nodes(self, mer_t_factors):
        if not mer_t_factors:
            return
        max_val = max(mer_t_factors.values()) if mer_t_factors.values() else 1.0

        # center node is the mer_name
        # Other nodes arranged in a circle
        mer_names = list(mer_t_factors.keys())
        if self.mer_name not in mer_names:
            mer_names.insert(0, self.mer_name)  # Ensure source mer is present (if not found)
        # Move source mer to front
        if self.mer_name in mer_names:
            mer_names.remove(self.mer_name)
            mer_names.insert(0, self.mer_name)
        
        center_x, center_y = 0, 0
        radius = 200
        num = len(mer_names)

        for i, mer in enumerate(mer_names):
            val = mer_t_factors.get(mer,0.0)
            normalized = val/max_val if max_val != 0 else 0
            # center node at the center
            if i == 0:
                x = center_x
                y = center_y
            else:
                angle = (2*math.pi/(num-1))*(i-1) if num>1 else 0
                x = center_x + radius * math.cos(angle)
                y = center_y + radius * math.sin(angle)

            self.draw_node(mer, val, x, y, is_source=(i==0))

        self.scene.setSceneRect(-300, -300, 600, 600)

    def draw_node(self, mer_name, tFactor, x, y, is_source=False):
        size = 40 if not is_source else 50
        rect = QRectF(x - size/2, y - size/2, size, size)
        color = QColor(255, 100, 100) if is_source else QColor(255,0,0)
        ellipse = QGraphicsEllipseItem(rect)
        ellipse.setBrush(QBrush(color))
        ellipse.setPen(QPen(Qt.black))
        self.scene.addItem(ellipse)

        text_item = QGraphicsTextItem(mer_name)
        text_item.setPos(x - size/2, y - size/2)
        text_item.setDefaultTextColor(Qt.white)
        self.scene.addItem(text_item)
