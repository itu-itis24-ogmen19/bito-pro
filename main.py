import sys
import numpy as np
from PySide6.QtWidgets import QApplication, QMainWindow

from ui.welcome_page import WelcomePage
from ui.file_upload_page import FileUploadPage
from ui.mer_list_page import MerListPage
from ui.protein_vispy_viewer_page import ProteinVisPyViewerPage

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein Folding Structures")
        self.resize(1000, 600)

        self.processed_files = []
        self.current_best_mer = None
        self.current_download_data = None
        self.current_interactions = None
        self.current_total_weight_sums = None
        self.current_all_distance_sums = None

        def go_to_welcome():
            self.setCentralWidget(WelcomePage(on_get_started=go_to_file_upload))

        def go_to_file_upload():
            self.setCentralWidget(FileUploadPage(
                on_back=go_to_welcome,
                on_mer_list=go_to_mer_list,
                processed_files=self.processed_files
            ))

        def go_to_mer_list(best_source_mer, download_data, processed_files, interactions, total_weight_sums, all_distance_sums):
            self.processed_files = processed_files
            self.current_best_mer = best_source_mer
            self.current_download_data = download_data
            self.current_interactions = interactions
            self.current_total_weight_sums = total_weight_sums
            self.current_all_distance_sums = all_distance_sums
            self.setCentralWidget(MerListPage(
                best_source_mer=best_source_mer,
                download_data=download_data,
                processed_files=self.processed_files,
                interactions=interactions,
                total_weight_sums=total_weight_sums,
                on_back=go_to_file_upload,
                on_view_mer=go_to_protein_viewer
            ))

        def go_to_protein_viewer(mer_name, pdb_content, interactions, _unused):
            # We ignore the old total_weight_sums passed here, because
            # we want distances specifically from the selected mer_name
            def on_back():
                go_to_mer_list(
                    self.current_best_mer,
                    self.current_download_data,
                    self.processed_files,
                    self.current_interactions,
                    self.current_total_weight_sums,
                    self.current_all_distance_sums
                )

            # Look up the specific distance map for the chosen Mer
            if self.current_all_distance_sums and mer_name in self.current_all_distance_sums:
                distance_map_for_selected_mer = self.current_all_distance_sums[mer_name]
            else:
                distance_map_for_selected_mer = {}

            self.setCentralWidget(
                ProteinVisPyViewerPage(
                    mer_name=mer_name,
                    pdb_content=pdb_content,
                    interactions=interactions,
                    total_weight_sums=distance_map_for_selected_mer,
                    on_back=on_back
                )
            )

        go_to_welcome()

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
