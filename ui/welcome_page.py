from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QPushButton
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtCore import Qt

import os

from ui.resource_locate import resource_path



class WelcomePage(QWidget):
    def __init__(self, on_get_started):
        super().__init__()
        layout = QVBoxLayout(self)
        self.setLayout(layout)
        
        # Load an image
        image_path = resource_path(os.path.join("assets", "images", "protein_icon.png"))
        pixmap = QPixmap(image_path)
        image_label = QLabel()
        image_label.setPixmap(pixmap.scaled(100, 100, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        image_label.setAlignment(Qt.AlignCenter)
        
        layout.addWidget(image_label, alignment=Qt.AlignCenter)

        title = QLabel("Protein Folding Structures")
        title.setFont(QFont("Montserrat", 24, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)

        subtitle = QLabel("Exploring the Complexities of Protein Structures\nIstanbul Technical University")
        subtitle.setAlignment(Qt.AlignCenter)

        layout.addWidget(title)
        layout.addWidget(subtitle)

        get_started_btn = QPushButton("Get Started")
        get_started_btn.clicked.connect(on_get_started)
        layout.addWidget(get_started_btn, alignment=Qt.AlignCenter)
