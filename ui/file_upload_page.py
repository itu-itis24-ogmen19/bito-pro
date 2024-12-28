# /Users/m3_air/24-25 fall/Bitirme/bito-pro-main/ui/file_upload_page.py

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog, QHBoxLayout, QListWidget
)
from PySide6.QtCore import Qt
import datetime
import os

from polmst import process_pdb_file

class ProcessedFile:
    def __init__(self, file_name, best_source_mer, download_data, interactions, total_weight_sums, timestamp):
        """
        Represents a processed PDB file with all relevant data.

        Parameters:
            file_name (str): The name of the PDB file.
            best_source_mer (str): The name of the best source Mer.
            download_data (dict): Mapping of Mer names to their enhanced PDB content.
            interactions (list): List of Interaction objects.
            total_weight_sums (dict): Mapping of Mer names to their total weight sums.
            timestamp (datetime): The timestamp when the file was processed.
        """
        self.file_name = file_name
        self.best_source_mer = best_source_mer
        self.download_data = download_data
        self.interactions = interactions
        self.total_weight_sums = total_weight_sums
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
        """
        Opens a file dialog for the user to select a PDB file.
        """
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PDB File", "", "PDB Files (*.pdb)")
        if file_path:
            self.selected_file = file_path
            self.status_label.setText(f"Selected File: {os.path.basename(file_path)}")
            self.process_btn.setEnabled(True)  # Enable process button since file is chosen
        else:
            self.status_label.setText("No file selected")
            self.process_btn.setEnabled(False)

    def process_file(self):
        """
        Processes the selected PDB file and updates the UI accordingly.
        """
        if not self.selected_file:
            self.status_label.setText("No file selected. Please choose a PDB file first.")
            return

        # Show processing status
        self.status_label.setText("Processing file, please wait...")
        self.select_btn.setEnabled(False)
        self.process_btn.setEnabled(False)

        #try:
        best_source_mer, mers, total_weight_sums, all_enhanced_pdbs, interactions = process_pdb_file(self.selected_file)
        #except Exception as e:
        #    self.status_label.setText(f"Error processing file: {str(e)}")
        #    self.select_btn.setEnabled(True)
        #    self.process_btn.setEnabled(True)
        #    return

        download_data = all_enhanced_pdbs

        processed = ProcessedFile(
            file_name=os.path.basename(self.selected_file),
            best_source_mer=best_source_mer,
            download_data=download_data,
            interactions=interactions,
            total_weight_sums=total_weight_sums,
            timestamp=datetime.datetime.now()
        )
        self.processed_files.append(processed)

        # Update sidebar
        self.file_list.addItem(f"{processed.file_name} - {processed.timestamp}")

        self.status_label.setText(f"Processing completed for '{processed.file_name}'!")
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)

        # Pass processed_files to on_mer_list
        self.on_mer_list(best_source_mer, download_data, self.processed_files, interactions, total_weight_sums)

    def open_previous_result(self, item):
        """
        Loads a previously processed PDB file when selected from the sidebar.
        """
        # Find the corresponding processed file based on the text
        selected_text = item.text()
        for pf in self.processed_files:
            entry_str = f"{pf.file_name} - {pf.timestamp}"
            if entry_str == selected_text:
                self.status_label.setText(f"Loading previous result for '{pf.file_name}'...")
                self.populate_mer_table(pf.best_source_mer, pf.download_data)
                self.on_mer_list(pf.best_source_mer, pf.download_data, self.processed_files, pf.interactions, pf.total_weight_sums)
                return
        self.status_label.setText("Could not find the selected previous result.")
