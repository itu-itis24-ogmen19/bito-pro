# ui/file_upload_page.py
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog, QHBoxLayout, QListWidget
)
from PySide6.QtCore import Qt
import datetime
import os

from polmst import process_pdb_file

class ProcessedFile:
    def __init__(self, file_name, best_source_mer, download_data, timestamp):
        self.file_name = file_name
        self.best_source_mer = best_source_mer
        self.download_data = download_data
        self.timestamp = timestamp

class FileUploadPage(QWidget):
    def __init__(self, on_back, on_mer_list, processed_files):
        super().__init__()
        self.on_back = on_back
        self.on_mer_list = on_mer_list
        self.processed_files = processed_files  # Store the processed_files list here

        self.selected_file = None
        
        main_layout = QHBoxLayout(self)
        
        # Sidebar for previous processes
        sidebar = QVBoxLayout()
        sidebar_label = QLabel("Previous Processes")
        sidebar.addWidget(sidebar_label)
        self.file_list = QListWidget()
        sidebar.addWidget(self.file_list)
        main_layout.addLayout(sidebar, 1)
        
        # Connect double-click signal on previous processes list
        self.file_list.itemDoubleClicked.connect(self.open_previous_result)
        
        # Main area
        main_area = QVBoxLayout()

        # Top bar for back button
        top_bar = QHBoxLayout()
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(self.back_btn, alignment=Qt.AlignLeft)
        top_bar.addStretch()
        main_area.addLayout(top_bar)

        self.status_label = QLabel("No file selected")
        self.status_label.setStyleSheet("color: #888;")
        self.status_label.setAlignment(Qt.AlignCenter)
        main_area.addWidget(self.status_label)

        self.select_btn = QPushButton("Select File")
        self.select_btn.clicked.connect(self.pick_file)
        self.select_btn.setToolTip("Select a PDB file for processing.")
        main_area.addWidget(self.select_btn, alignment=Qt.AlignCenter)

        self.process_btn = QPushButton("Process File")
        self.process_btn.clicked.connect(self.process_file)
        self.process_btn.setToolTip("Run the analysis on the selected PDB file.")
        self.process_btn.setEnabled(False)  # Disabled until a file is selected
        main_area.addWidget(self.process_btn, alignment=Qt.AlignCenter)
        
        main_layout.addLayout(main_area, 3)

        # Populate the sidebar with any previously processed files
        for pf in self.processed_files:
            self.file_list.addItem(f"{pf.file_name} - {pf.timestamp}")

    def pick_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PDB File", "", "PDB Files (*.pdb)")
        if file_path:
            self.selected_file = file_path
            self.status_label.setText(f"Selected File: {os.path.basename(file_path)}")
            self.process_btn.setEnabled(True)  # Now enable process button since file is chosen
        else:
            self.status_label.setText("No file selected")
            self.process_btn.setEnabled(False)

    def process_file(self):
        if not self.selected_file:
            self.status_label.setText("No file selected. Please choose a PDB file first.")
            return

        # Show processing status
        self.status_label.setText("Processing file, please wait...")
        self.select_btn.setEnabled(False)
        self.process_btn.setEnabled(False)

        best_source_mer, mers, shortest_paths, all_enhanced_pdbs = process_pdb_file(self.selected_file)
        download_data = all_enhanced_pdbs

        processed = ProcessedFile(
            file_name=os.path.basename(self.selected_file),
            best_source_mer=best_source_mer,
            download_data=download_data,
            timestamp=datetime.datetime.now()
        )
        self.processed_files.append(processed)

        # Update sidebar
        self.file_list.addItem(f"{processed.file_name} - {processed.timestamp}")

        self.status_label.setText(f"Processing completed for '{processed.file_name}'!")
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)

        # Pass processed_files to on_mer_list
        self.on_mer_list(best_source_mer, download_data, self.processed_files)

    def open_previous_result(self, item):
        # Find the corresponding processed file based on the text
        selected_text = item.text()
        for pf in self.processed_files:
            entry_str = f"{pf.file_name} - {pf.timestamp}"
            if entry_str == selected_text:
                self.status_label.setText(f"Loading results for '{pf.file_name}'...")
                self.on_mer_list(pf.best_source_mer, pf.download_data, self.processed_files)
                return
        self.status_label.setText("Could not find the selected previous result.")
