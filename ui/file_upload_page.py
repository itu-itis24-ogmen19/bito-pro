from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog,
    QHBoxLayout, QListWidget, QSpacerItem, QSizePolicy, QProgressBar
)
from PySide6.QtCore import Qt, QThread, Signal
import datetime
import os

from polmst import process_pdb_file  # Assuming this is a blocking call

class ProcessedFile:
    def __init__(
        self,
        file_name,
        file_path,
        best_source_mer,
        mers,
        download_data,
        interactions,
        total_weight_sums,
        timestamp,
        all_distance_sums,
        original_content  # NEW: store the original file content
    ):
        self.file_name = file_name
        self.file_path = file_path
        self.best_source_mer = best_source_mer
        self.mers = mers
        self.download_data = download_data
        self.interactions = interactions
        self.total_weight_sums = total_weight_sums
        self.timestamp = timestamp
        self.all_distance_sums = all_distance_sums
        self.original_content = original_content  # NEW

class WorkerThread(QThread):
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

    def run(self):
        try:
            best_source_mer, mers, total_weight_sums, interactions, all_distance_sums = process_pdb_file(self.file_path)
            self.finished.emit({
                'best_source_mer': best_source_mer,
                'mers': mers,
                'total_weight_sums': total_weight_sums,
                'interactions': interactions,
                'all_distance_sums': all_distance_sums
            })
        except Exception as e:
            self.error.emit(str(e))

class FileUploadPage(QWidget):
    def __init__(self, on_back, on_mer_list, on_view_raw, processed_files):
        super().__init__()
        self.on_back = on_back
        self.on_mer_list = on_mer_list
        self.on_view_raw = on_view_raw   # New callback for raw view
        self.processed_files = processed_files

        self.selected_file = None
        self.worker = None

        self.init_ui()

    def init_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(20, 20, 20, 20)
        main_layout.setSpacing(20)

        # Sidebar for previous processes
        sidebar = QVBoxLayout()
        sidebar_label = QLabel("Previous Processes")
        sidebar_label.setStyleSheet("font-weight: bold; font-size: 16px;")
        sidebar.addWidget(sidebar_label)

        self.file_list = QListWidget()
        self.file_list.setStyleSheet("""
            QListWidget {
                border: 1px solid #ccc;
                border-radius: 5px;
                background-color: #2C2F33;
                color: white;
            }
            QListWidget::item {
                padding: 10px;
            }
            QListWidget::item:selected {
                background-color: #7289DA;
                color: white;
            }
        """)
        sidebar.addWidget(self.file_list)

        main_layout.addLayout(sidebar, 1)

        # Main area
        main_area = QVBoxLayout()

        # Top bar with a Back button
        top_bar = QHBoxLayout()
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        top_bar.addWidget(self.back_btn, alignment=Qt.AlignLeft)
        top_bar.addStretch()
        main_area.addLayout(top_bar)

        # Spacer
        main_area.addSpacerItem(QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Status Label
        self.status_label = QLabel("No file selected")
        self.status_label.setStyleSheet("color: #CCCCCC; font-size: 14px;")
        self.status_label.setAlignment(Qt.AlignCenter)
        main_area.addWidget(self.status_label)

        # Spacer
        main_area.addSpacerItem(QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Select File Button
        self.select_btn = QPushButton("Select File")
        self.select_btn.setStyleSheet("""
            QPushButton {
                background-color: #23272A;
                color: white;
                padding: 10px 30px;
                border: 2px solid #7289DA;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #2C2F33;
            }
            QPushButton:disabled {
                background-color: #555555;
                border: 2px solid #888888;
            }
        """)
        self.select_btn.clicked.connect(self.pick_file)
        self.select_btn.setToolTip("Select a PDB file for processing.")
        main_area.addWidget(self.select_btn, alignment=Qt.AlignCenter)

        # Process File Button
        self.process_btn = QPushButton("Process File")
        self.process_btn.setStyleSheet("""
            QPushButton {
                background-color: #23272A;
                color: white;
                padding: 10px 30px;
                border: 2px solid #7289DA;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #2C2F33;
            }
            QPushButton:disabled {
                background-color: #555555;
                border: 2px solid #888888;
            }
        """)
        self.process_btn.clicked.connect(self.process_file)
        self.process_btn.setToolTip("Run the analysis on the selected PDB file.")
        self.process_btn.setEnabled(False)
        main_area.addWidget(self.process_btn, alignment=Qt.AlignCenter)

        # View Raw PDB Button (New)
        self.view_raw_btn = QPushButton("View Raw PDB")
        self.view_raw_btn.setStyleSheet("""
            QPushButton {
                background-color: #23272A;
                color: white;
                padding: 10px 30px;
                border: 2px solid #7289DA;
                border-radius: 5px;
                font-size: 14px;
            }
            QPushButton:hover {
                background-color: #2C2F33;
            }
            QPushButton:disabled {
                background-color: #555555;
                border: 2px solid #888888;
            }
        """)
        self.view_raw_btn.setToolTip("Display the selected PDB file without processing it.")
        self.view_raw_btn.clicked.connect(self.view_raw_pdb)
        self.view_raw_btn.setEnabled(False)
        main_area.addWidget(self.view_raw_btn, alignment=Qt.AlignCenter)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setFixedHeight(20)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                border: 1px solid #ccc;
                border-radius: 5px;
                text-align: center;
                background-color: #2C2F33;
                color: white;
            }
            QProgressBar::chunk {
                background-color: #7289DA;
                width: 20px;
            }
        """)
        main_area.addWidget(self.progress_bar)

        # Spacer
        main_area.addSpacerItem(QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding))

        main_layout.addLayout(main_area, 3)

        # Populate sidebar with previously processed files
        for pf in self.processed_files:
            self.file_list.addItem(f"{pf.file_name} - {pf.timestamp.strftime('%Y-%m-%d %H:%M:%S')}")

        # Connect list double-click
        self.file_list.itemDoubleClicked.connect(self.open_previous_result)

    def pick_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Select PDB File", "", "PDB Files (*.pdb)"
        )
        if file_path:
            self.selected_file = file_path
            self.status_label.setText(f"Selected File: {os.path.basename(file_path)}")
            self.process_btn.setEnabled(True)
            self.view_raw_btn.setEnabled(True)  # Enable raw view button
        else:
            self.status_label.setText("No file selected")
            self.process_btn.setEnabled(False)
            self.view_raw_btn.setEnabled(False)

    def process_file(self):
        if not self.selected_file:
            self.status_label.setText("No file selected. Please choose a PDB file first.")
            return

        # Update UI to indicate processing
        self.status_label.setText(
            "PDB file is being processed, it may take a couple of minutes depending on the size of the file."
        )
        self.select_btn.setEnabled(False)
        self.process_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate progress

        # Start processing in a separate thread
        self.worker = WorkerThread(self.selected_file)
        self.worker.finished.connect(self.on_processing_finished)
        self.worker.error.connect(self.on_processing_error)
        self.worker.start()

    def view_raw_pdb(self):
        if not self.selected_file:
            self.status_label.setText("No file selected.")
            return
        with open(self.selected_file, 'r') as f:
            pdb_content = f.read()
        self.on_view_raw(pdb_content, self.selected_file)

    def on_processing_finished(self, result):
        # Stop the progress bar
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)

        # Re-enable buttons
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)

        # Read the original file content (unprocessed)
        with open(self.selected_file, 'r') as f:
            original_content = f.read()

        # Create ProcessedFile instance (note the new original_content parameter)
        processed = ProcessedFile(
            file_name=os.path.basename(self.selected_file),
            file_path=self.selected_file,
            best_source_mer=result['best_source_mer'],
            mers=result['mers'],
            download_data={mer_name: None for mer_name in result['mers'].keys()},
            interactions=result['interactions'],
            total_weight_sums=result['total_weight_sums'],
            timestamp=datetime.datetime.now(),
            all_distance_sums=result['all_distance_sums'],
            original_content=original_content  # NEW: store the original PDB content
        )
        self.processed_files.append(processed)
        self.file_list.addItem(f"{processed.file_name} - {processed.timestamp.strftime('%Y-%m-%d %H:%M:%S')}")

        self.status_label.setText(f"Processing completed for '{processed.file_name}'!")

        # Proceed to the next step
        self.on_mer_list(
            processed.best_source_mer,
            processed.download_data,  # now has keys == all mer names
            self.processed_files,
            processed.interactions,
            processed.total_weight_sums,
            processed.all_distance_sums
        )

    def on_processing_error(self, error_message):
        # Stop the progress bar
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)

        # Re-enable buttons
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)

        self.status_label.setText(f"Error processing file: {error_message}")

    def open_previous_result(self, item):
        selected_text = item.text()
        for pf in self.processed_files:
            entry_str = f"{pf.file_name} - {pf.timestamp.strftime('%Y-%m-%d %H:%M:%S')}"
            if entry_str == selected_text:
                self.status_label.setText(f"Loading previous result for '{pf.file_name}'...")
                self.on_mer_list(
                    pf.best_source_mer,
                    pf.download_data,  # may contain None for each mer_name
                    self.processed_files,
                    pf.interactions,
                    pf.total_weight_sums,
                    pf.all_distance_sums
                )
                return
        self.status_label.setText("Could not find the selected previous result.")
