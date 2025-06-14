from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPushButton, QFileDialog,
    QHBoxLayout, QListWidget, QSpacerItem, QSizePolicy,
    QProgressBar, QMessageBox       # ← added QMessageBox
)
from PySide6.QtCore import Qt, QThread, Signal
import datetime
import os

from polmst import process_pdb_file


# ────────────────────────────────────────────────────────────────
#  Data container for sidebar history
# ────────────────────────────────────────────────────────────────
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
        original_content,
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
        self.original_content = original_content


# ────────────────────────────────────────────────────────────────
#  Worker thread (calls process_pdb_file)
# ────────────────────────────────────────────────────────────────
class WorkerThread(QThread):
    finished = Signal(object)
    error = Signal(str)

    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

    def run(self):
        try:
            (
                best_source_mer,
                mers,
                total_weight_sums,
                interactions,
                all_distance_sums,
                islands_removed,          # ← NEW info
            ) = process_pdb_file(self.file_path)

            self.finished.emit(
                {
                    "best_source_mer": best_source_mer,
                    "mers": mers,
                    "total_weight_sums": total_weight_sums,
                    "interactions": interactions,
                    "all_distance_sums": all_distance_sums,
                    "islands_removed": islands_removed,  # ← NEW key
                }
            )
        except Exception as e:
            self.error.emit(str(e))


# ────────────────────────────────────────────────────────────────
#  Main upload page
# ────────────────────────────────────────────────────────────────
class FileUploadPage(QWidget):
    def __init__(self, on_back, on_mer_list, on_view_raw, processed_files):
        super().__init__()
        self.on_back = on_back
        self.on_mer_list = on_mer_list
        self.on_view_raw = on_view_raw
        self.processed_files = processed_files
        self.selected_file = None
        self.worker = None
        self.init_ui()

    # -----------------------------------------------------------------
    #  UI builder (unchanged except for extra imports)
    # -----------------------------------------------------------------
    def init_ui(self):
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(20, 20, 20, 20)

        # Sidebar (history)
        sidebar = QVBoxLayout()
        sidebar_label = QLabel("Previous Processes")
        sidebar_label.setStyleSheet("font-weight: bold; font-size: 16px;")
        sidebar.addWidget(sidebar_label)

        self.file_list = QListWidget()
        sidebar.addWidget(self.file_list)
        main_layout.addLayout(sidebar, 1)

        # Main area
        main = QVBoxLayout()
        top = QHBoxLayout()
        self.back_btn = QPushButton("Back")
        self.back_btn.clicked.connect(self.on_back)
        top.addWidget(self.back_btn, alignment=Qt.AlignLeft)
        top.addStretch()
        main.addLayout(top)

        main.addSpacerItem(QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.status_label = QLabel("No file selected")
        self.status_label.setStyleSheet("color: #CCCCCC; font-size: 14px;")
        self.status_label.setAlignment(Qt.AlignCenter)
        main.addWidget(self.status_label)

        main.addSpacerItem(QSpacerItem(20, 20, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.select_btn = QPushButton("Select File")
        self.select_btn.clicked.connect(self.pick_file)
        main.addWidget(self.select_btn, alignment=Qt.AlignCenter)

        self.process_btn = QPushButton("Process File")
        self.process_btn.setEnabled(False)
        self.process_btn.clicked.connect(self.process_file)
        main.addWidget(self.process_btn, alignment=Qt.AlignCenter)

        self.view_raw_btn = QPushButton("View Raw PDB")
        self.view_raw_btn.setEnabled(False)
        self.view_raw_btn.clicked.connect(self.view_raw_pdb)
        main.addWidget(self.view_raw_btn, alignment=Qt.AlignCenter)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main.addWidget(self.progress_bar)

        main_layout.addLayout(main, 3)

        # populate sidebar with history
        for pf in self.processed_files:
            self.file_list.addItem(
                f"{pf.file_name} - {pf.timestamp.strftime('%Y-%m-%d %H:%M:%S')}"
            )
        self.file_list.itemDoubleClicked.connect(self.open_previous_result)

    # -----------------------------------------------------------------
    #  File-choosing helpers
    # -----------------------------------------------------------------
    def pick_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PDB File", "", "PDB Files (*.pdb)")
        if file_path:
            self.selected_file = file_path
            self.status_label.setText(f"Selected File: {os.path.basename(file_path)}")
            self.process_btn.setEnabled(True)
            self.view_raw_btn.setEnabled(True)
        else:
            self.status_label.setText("No file selected")
            self.process_btn.setEnabled(False)
            self.view_raw_btn.setEnabled(False)

    def process_file(self):
        if not self.selected_file:
            return
        self.status_label.setText("Processing… please wait.")
        self.select_btn.setEnabled(False)
        self.process_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)

        self.worker = WorkerThread(self.selected_file)
        self.worker.finished.connect(self.on_processing_finished)
        self.worker.error.connect(self.on_processing_error)
        self.worker.start()

    def view_raw_pdb(self):
        if not self.selected_file:
            return
        with open(self.selected_file, "r") as f:
            self.on_view_raw(f.read(), self.selected_file)

    # -----------------------------------------------------------------
    #  Callbacks
    # -----------------------------------------------------------------
    def on_processing_finished(self, result):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)

        with open(self.selected_file, "r") as f:
            original_content = f.read()

        processed = ProcessedFile(
            file_name=os.path.basename(self.selected_file),
            file_path=self.selected_file,
            best_source_mer=result["best_source_mer"],
            mers=result["mers"],
            download_data={m: None for m in result["mers"]},
            interactions=result["interactions"],
            total_weight_sums=result["total_weight_sums"],
            timestamp=datetime.datetime.now(),
            all_distance_sums=result["all_distance_sums"],
            original_content=original_content,
        )
        self.processed_files.append(processed)
        self.file_list.addItem(
            f"{processed.file_name} - {processed.timestamp.strftime('%Y-%m-%d %H:%M:%S')}"
        )

        # ── NEW POP-UP ──────────────────────────────────────────────
        if result.get("islands_removed"):
            QMessageBox.warning(
                self,
                "Small Island Removed",
                "One or more small, disconnected components of the protein "
                "were removed. Only the largest connected component is kept "
                "so that distance calculations remain meaningful.",
            )
        # ────────────────────────────────────────────────────────────

        self.status_label.setText(f"Processing completed for '{processed.file_name}'")

        self.on_mer_list(
            processed.best_source_mer,
            processed.download_data,
            self.processed_files,
            processed.interactions,
            processed.total_weight_sums,
            processed.all_distance_sums,
        )

    def on_processing_error(self, msg):
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 100)
        self.select_btn.setEnabled(True)
        self.process_btn.setEnabled(True)
        self.status_label.setText(f"Error: {msg}")

    def open_previous_result(self, item):
        text = item.text()
        for pf in self.processed_files:
            if text.startswith(pf.file_name):
                self.on_mer_list(
                    pf.best_source_mer,
                    pf.download_data,
                    self.processed_files,
                    pf.interactions,
                    pf.total_weight_sums,
                    pf.all_distance_sums,
                )
                break
