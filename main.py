import sys
from PySide6.QtWidgets import QApplication, QMainWindow

from ui.welcome_page import WelcomePage
from ui.file_upload_page import FileUploadPage
from ui.mer_list_page import MerListPage
from ui.protein_viewer_page import ProteinViewerPage

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein Folding Structures")
        self.resize(1000, 600)

        # Maintain a global list of processed_files to keep track across navigation if desired
        # If you prefer, you can let FileUploadPage manage it, but this is more centralized.
        self.processed_files = []

        def go_to_welcome():
            # Welcome page leads to file upload when "Get Started" is clicked
            self.setCentralWidget(WelcomePage(on_get_started=go_to_file_upload))

        def go_to_file_upload():
            # When returning to file upload, pass the processed_files to maintain history
            self.setCentralWidget(FileUploadPage(
                on_back=go_to_welcome,
                on_mer_list=go_to_mer_list,
                processed_files=self.processed_files
            ))

        def go_to_mer_list(best_source_mer, download_data, processed_files):
            # Update the global processed_files to keep everything in sync
            self.processed_files = processed_files
            self.setCentralWidget(MerListPage(
                best_source_mer=best_source_mer,
                download_data=download_data,
                processed_files=self.processed_files,
                on_back=go_to_file_upload,
                on_view_mer=go_to_protein_viewer
            ))

        def go_to_protein_viewer(mer_name, pdb_content):
            # If ProteinViewerPage provides a back function returning to mer_list,
            # we must pass processed_files here as well. 
            # We'll assume we return to the MerListPage showing that same mer.
            
            # We need the best_source_mer and full download_data for that mer again.
            # If you need the full context, store it from MerListPage or pass it along.
            # For simplicity, let's assume we have minimal info here.
            # Typically, you'd have a reference to the current best_source_mer, download_data,
            # and processed_files from the MerListPage in closures or through instance attributes.
            
            # Let's store them on the window as attributes when go_to_mer_list is called:
            # Assume we store them as self.current_best_mer, self.current_download_data
            # whenever go_to_mer_list is called:
            # This requires a small addition to go_to_mer_list:
            #   self.current_best_mer = best_source_mer
            #   self.current_download_data = download_data

            # We'll add that now:
            # Modify go_to_mer_list:
            # def go_to_mer_list(best_source_mer, download_data, processed_files):
            #     self.processed_files = processed_files
            #     self.current_best_mer = best_source_mer
            #     self.current_download_data = download_data
            #     self.setCentralWidget(MerListPage(...))

            # Similarly here, on back from ProteinViewerPage:
            def on_back():
                # Return to MerListPage with the stored current context
                go_to_mer_list(self.current_best_mer, self.current_download_data, self.processed_files)

            self.setCentralWidget(ProteinViewerPage(
                mer_name=mer_name,
                pdb_content=pdb_content,
                on_back=on_back
            ))

        # Attach current best mer and download_data attributes to store context when mer_list is shown
        def _go_to_mer_list_wrapper(best_source_mer, download_data, processed_files):
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

        # Replace the original go_to_mer_list with this wrapper
        go_to_mer_list = _go_to_mer_list_wrapper

        # Start with welcome page
        go_to_welcome()

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
