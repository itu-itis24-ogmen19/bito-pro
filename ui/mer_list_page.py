import os
from PySide6.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QLabel, QListWidget,
    QHeaderView, QTableWidget, QTableWidgetItem, QPushButton, QFileDialog
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QFont

from ui.resource_locate import resource_path
from polmst import generate_enhanced_pdb  # Import the function for generating enhanced PDB content


class MerListPage(QWidget):
    def __init__(self, best_source_mer, download_data, processed_files, interactions,
                 total_weight_sums, on_back, on_view_mer):
        super().__init__()
        self.best_source_mer = best_source_mer
        self.download_data = download_data
        self.processed_files = processed_files
        self.interactions = interactions
        self.total_weight_sums = total_weight_sums
        self.on_back = on_back
        self.on_view_mer = on_view_mer

        self.best_source_row_index = None

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

        # Top bar for Back and global "View Original 3D" button
        top_bar = QHBoxLayout()
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(self.back_btn, alignment=Qt.AlignLeft)

        # Global button to view the original (unprocessed) 3D structure.
        self.view_original_btn = QPushButton("View Original 3D")
        self.view_original_btn.setToolTip("View the original 3D structure as uploaded.")
        self.view_original_btn.clicked.connect(self.view_original_global)
        top_bar.addWidget(self.view_original_btn, alignment=Qt.AlignRight)

        main_area.addLayout(top_bar)

        # A clickable label for best source
        self.title_label = QLabel()
        self.title_label.setAlignment(Qt.AlignCenter)
        self.title_label.setTextFormat(Qt.RichText)
        self.title_label.setOpenExternalLinks(False)
        self.title_label.linkActivated.connect(self.locate_best_source_in_table)
        main_area.addWidget(self.title_label)

        # Table with 3 columns: Center Chosen Mer, View 3D Graph, Download PDB
        self.mer_table = QTableWidget()
        self.mer_table.setColumnCount(3)
        self.mer_table.setHorizontalHeaderLabels([
            "Center Chosen Mer",   # was "Mer"
            "View 3D Graph",       # enhanced view
            "Download PDB"         # download enhanced PDB file
        ])
        self.mer_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.mer_table.verticalHeader().setVisible(False)
        main_area.addWidget(self.mer_table)

        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: #888;")
        self.status_label.setAlignment(Qt.AlignCenter)
        main_area.addWidget(self.status_label)

        main_layout.addLayout(sidebar, 1)
        main_layout.addLayout(main_area, 3)

        # Populate the table
        self.populate_mer_table(self.best_source_mer, self.download_data)

    def populate_mer_table(self, best_source_mer, download_data):
        self.download_data = download_data
        self.best_source_mer = best_source_mer

        mer_names = list(self.download_data.keys())
        self.mer_table.setRowCount(len(mer_names))

        arrow_icon_path = resource_path(os.path.join("assets", "images", "yellow_arrow.png"))
        arrow_icon = QIcon(arrow_icon_path) if os.path.exists(arrow_icon_path) else QIcon()

        magnify_icon_path = resource_path(os.path.join("assets", "images", "magnify_icon.png"))
        magnify_icon = QIcon(magnify_icon_path) if os.path.exists(magnify_icon_path) else QIcon()

        download_icon_path = resource_path(os.path.join("assets", "images", "download_icon.png"))
        download_icon = QIcon(download_icon_path) if os.path.exists(download_icon_path) else QIcon()

        self.best_source_row_index = None

        for row, mer_name in enumerate(mer_names):
            name_item = QTableWidgetItem(mer_name)
            name_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)

            # CASE-INSENSITIVE match, strip extra spaces
            if mer_name.strip().lower() == self.best_source_mer.strip().lower():
                print(f"[DEBUG] Found best source match: '{mer_name}' == '{self.best_source_mer}'")
                name_item.setIcon(arrow_icon)
                self.best_source_row_index = row

            self.mer_table.setItem(row, 0, name_item)

            # Create "View" button with magnifying glass icon for the enhanced view
            view_btn = QPushButton()
            if not magnify_icon.isNull():
                view_btn.setIcon(magnify_icon)
            else:
                view_btn.setText("View")
            view_btn.setToolTip("View the enhanced 3D Graph for this Mer.")
            view_btn.clicked.connect(lambda _, mn=mer_name: self.view_mer(mn))
            self.mer_table.setCellWidget(row, 1, view_btn)

            # Create "Download" button with download icon
            download_btn = QPushButton()
            if not download_icon.isNull():
                download_btn.setIcon(download_icon)
            else:
                download_btn.setText("Download")
            download_btn.setToolTip("Download the enhanced PDB file for this Mer.")
            download_btn.clicked.connect(lambda _, mn=mer_name: self.download_pdb(mn))
            self.mer_table.setCellWidget(row, 2, download_btn)

        # Make the best source clickable in the title
        self.title_label.setText(
            f"Mers - Best Source: <a href='#'>{self.best_source_mer}</a>"
        )

        self.status_label.setText(f"Loaded {len(mer_names)} Mers.")

    def view_mer(self, mer_name):
        """On-demand generation of enhanced PDB content for 'mer_name'."""
        pf = self.find_processed_file()
        if not pf:
            self.status_label.setText("No PDB data found.")
            return

        if mer_name in pf.all_distance_sums:
            dist_map = pf.all_distance_sums[mer_name]
            pdb_content = generate_enhanced_pdb(dist_map, mer_name, pf.file_path, pf.mers)
        else:
            pdb_content = ""

        if not pdb_content:
            self.status_label.setText(f"No PDB data found for {mer_name}.")
            return

        self.status_label.setText(f"Viewing enhanced 3D for '{mer_name}'...")
        self.on_view_mer(mer_name, pdb_content, pf.interactions, pf.total_weight_sums)

    def view_original_global(self):
        """
        Opens a 3D view using the original (unprocessed) PDB file content.
        Since the original structure is identical regardless of the chosen Mer,
        this global view is shared.
        """
        pf = self.find_processed_file()
        if not pf:
            self.status_label.setText("No PDB data found.")
            return

        if not hasattr(pf, "original_content") or not pf.original_content:
            self.status_label.setText("No original PDB data available.")
            return

        self.status_label.setText("Viewing original 3D structure...")
        # Pass a dummy name since the structure is global.
        self.on_view_mer("Original Structure", pf.original_content, pf.interactions, pf.total_weight_sums)

    def download_pdb(self, mer_name):
        """Generate and save a PDB file for the chosen Mer."""
        pf = self.find_processed_file()
        if not pf:
            self.status_label.setText("No PDB data found.")
            return

        if mer_name not in pf.all_distance_sums:
            self.status_label.setText(f"No PDB data found for {mer_name}.")
            return

        dist_map = pf.all_distance_sums[mer_name]
        pdb_content = generate_enhanced_pdb(dist_map, mer_name, pf.file_path, pf.mers)
        if not pdb_content:
            self.status_label.setText(f"No PDB data generated for {mer_name}.")
            return

        save_path, _ = QFileDialog.getSaveFileName(
            self, "Save PDB File", f"{mer_name}.pdb", "PDB Files (*.pdb)"
        )
        if save_path:
            with open(save_path, 'w') as f:
                f.write(pdb_content)
            self.status_label.setText(f"{mer_name}.pdb saved successfully.")

    def locate_best_source_in_table(self, _link=None):
        """
        Highlight/scroll to the best source Mer row if found,
        otherwise a status message says it wasn't located.
        """
        if self.best_source_row_index is not None:
            self.mer_table.selectRow(self.best_source_row_index)
            self.mer_table.scrollToItem(
                self.mer_table.item(self.best_source_row_index, 0)
            )
            self.status_label.setText("Best source Mer located in the list.")
        else:
            self.status_label.setText("Could not locate best source Mer in the table.")

    def load_previous_result(self, item):
        """Double-click on a previous result in the sidebar to reload its data."""
        selected_text = item.text()
        for pf in self.processed_files:
            entry_str = f"{pf.file_name} - {pf.timestamp}"
            if entry_str == selected_text:
                self.status_label.setText(f"Loading previous result for '{pf.file_name}'...")
                self.populate_mer_table(pf.best_source_mer, pf.download_data)
                self.on_view_mer(
                    pf.best_source_mer,
                    pf.download_data,
                    self.processed_files,
                    pf.interactions,
                    pf.total_weight_sums
                )
                return
        self.status_label.setText("Could not find the selected previous result.")

    def find_processed_file(self):
        """Find the ProcessedFile whose best_source_mer matches the current best_source_mer (case-insensitive)."""
        for pf in self.processed_files:
            if pf.best_source_mer.strip().lower() == self.best_source_mer.strip().lower():
                return pf
        return None
