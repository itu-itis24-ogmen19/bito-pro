import os
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QListWidget, 
    QHeaderView, QTableWidget, QTableWidgetItem, QPushButton, QFileDialog
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QFont

from ui.resource_locate import resource_path



class MerListPage(QWidget):
    def __init__(self, best_source_mer, download_data, processed_files, on_back, on_view_mer):
        super().__init__()
        self.best_source_mer = best_source_mer
        self.download_data = download_data
        self.processed_files = processed_files
        self.on_back = on_back
        self.on_view_mer = on_view_mer

        main_layout = QHBoxLayout(self)
        self.setLayout(main_layout)

        # Sidebar for Previous Processes
        sidebar = QVBoxLayout()
        sidebar_label = QLabel("Previous Processes")
        sidebar_label.setFont(QFont("Arial", 10, QFont.Bold))
        sidebar.addWidget(sidebar_label)

        self.file_list = QListWidget()
        for pf in self.processed_files:
            self.file_list.addItem(f"{pf.file_name} - {pf.timestamp}")
        sidebar.addWidget(self.file_list)

        self.file_list.itemDoubleClicked.connect(self.load_previous_result)

        # Main area
        main_area = QVBoxLayout()

        # Top bar for back button
        top_bar = QHBoxLayout()
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(self.back_btn, alignment=Qt.AlignLeft)
        top_bar.addStretch()
        main_area.addLayout(top_bar)

        self.title_label = QLabel(f"Mers - Best Source: {self.best_source_mer}")
        self.title_label.setAlignment(Qt.AlignCenter)
        main_area.addWidget(self.title_label)

        self.mer_table = QTableWidget()
        self.mer_table.setColumnCount(3)
        self.mer_table.setHorizontalHeaderLabels(["Mer", "", "Download"])
        self.mer_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.mer_table.verticalHeader().setVisible(False)
        main_area.addWidget(self.mer_table)

        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: #888;") 
        self.status_label.setAlignment(Qt.AlignCenter)
        main_area.addWidget(self.status_label)

        main_layout.addLayout(sidebar, 1)
        main_layout.addLayout(main_area, 3)

        self.populate_mer_table(self.best_source_mer, self.download_data)

    def populate_mer_table(self, best_source_mer, download_data):
        self.download_data = download_data
        self.best_source_mer = best_source_mer

        mer_names = list(self.download_data.keys())
        self.mer_table.setRowCount(len(mer_names))

        check_icon_path = resource_path(os.path.join("assets", "images", "yellow_check.png"))
        check_icon = QIcon(check_icon_path) if os.path.exists(check_icon_path) else QIcon()

        for row, mer_name in enumerate(mer_names):
            name_item = QTableWidgetItem(mer_name)
            name_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)

            if mer_name == self.best_source_mer and not check_icon.isNull():
                name_item.setIcon(check_icon)

            self.mer_table.setItem(row, 0, name_item)

            view_btn = QPushButton("View")
            view_btn.setToolTip("View the Dijkstra calculation visualization for this Mer.")
            view_btn.clicked.connect(lambda _, mn=mer_name: self.view_mer(mn))
            self.mer_table.setCellWidget(row, 1, view_btn)

            download_btn = QPushButton("Download")
            download_btn.setToolTip("Download the enhanced PDB file for this Mer.")
            download_btn.clicked.connect(lambda _, mn=mer_name: self.download_pdb(mn))
            self.mer_table.setCellWidget(row, 2, download_btn)

        self.title_label.setText(f"Mers - Best Source: {self.best_source_mer}")
        self.status_label.setText(f"Loaded {len(mer_names)} Mers.")

    def view_mer(self, mer_name):
        pdb_content = self.download_data.get(mer_name, "")
        if not pdb_content:
            self.status_label.setText(f"No PDB data found for {mer_name}.")
            return
        self.status_label.setText(f"Viewing Mer '{mer_name}'...")
        self.on_view_mer(mer_name, pdb_content)

    def download_pdb(self, mer_name):
        pdb_content = self.download_data.get(mer_name, "")
        if not pdb_content:
            self.status_label.setText(f"No PDB data found for {mer_name}.")
            return

        save_path, _ = QFileDialog.getSaveFileName(self, "Save PDB File", f"{mer_name}.pdb", "PDB Files (*.pdb)")
        if save_path:
            with open(save_path, 'w') as f:
                f.write(pdb_content)
            self.status_label.setText(f"{mer_name}.pdb saved successfully.")

    def load_previous_result(self, item):
        selected_text = item.text()
        for pf in self.processed_files:
            entry_str = f"{pf.file_name} - {pf.timestamp}"
            if entry_str == selected_text:
                self.status_label.setText(f"Loading previous result for '{pf.file_name}'...")
                self.populate_mer_table(pf.best_source_mer, pf.download_data)
                return
        self.status_label.setText("Could not find the selected previous result.")
