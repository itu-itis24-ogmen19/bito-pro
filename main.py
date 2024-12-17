# main.py
import sys
import numpy as np

from PySide6.QtWidgets import QApplication, QMainWindow

from ui.welcome_page import WelcomePage
from ui.file_upload_page import FileUploadPage
from ui.mer_list_page import MerListPage
from ui.protein_vispy_viewer_page import ProteinVisPyViewerPage  # New import

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein Folding Structures")
        self.resize(1000, 600)

        self.processed_files = []

        def go_to_welcome():
            self.setCentralWidget(WelcomePage(on_get_started=go_to_file_upload))

        def go_to_file_upload():
            self.setCentralWidget(FileUploadPage(
                on_back=go_to_welcome,
                on_mer_list=go_to_mer_list,
                processed_files=self.processed_files
            )) 

        def go_to_mer_list(best_source_mer, download_data, processed_files):
            self.processed_files = processed_files
            self.current_best_mer = best_source_mer
            self.current_download_data = download_data
            self.setCentralWidget(MerListPage(
                best_source_mer=best_source_mer,
                download_data=download_data,
                processed_files=self.processed_files,
                on_back=go_to_file_upload,
                on_view_mer=go_to_protein_viewer
            ))

        def go_to_protein_viewer(mer_name, pdb_content):
            def on_back():
                go_to_mer_list(self.current_best_mer, self.current_download_data, self.processed_files)
            self.setCentralWidget(ProteinVisPyViewerPage(
                mer_name=mer_name,
                pdb_content=pdb_content,
                on_back=on_back
            ))

        go_to_welcome()

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
