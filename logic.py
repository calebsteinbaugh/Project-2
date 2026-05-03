from PyQt6.QtWidgets import QMainWindow
from sequence_analysis_gui import Ui_Sequence_analysis_window
from RNA_protein_Class import RNASequence


class LogicController:
    """
    Controller class responsible for handling user input, launching the
    analysis window, and coordinating RNA sequence analysis operations.
    """

    def __init__(self, ui, window: QMainWindow) -> None:
        """
        Initialize the LogicController.

        Args:
            ui: The main input window UI.
            window: The main application window.
        """
        self.ui = ui
        self.window: QMainWindow = window
        self.rna_sequence: RNASequence | None = None
        
        self.analysis_window: QMainWindow | None = None
        self.analysis_ui: Ui_Sequence_analysis_window | None = None
        self.current_results: str = ""
        
        self.connect_signals()
        
    def connect_signals(self) -> None:
        """
        Connect UI buttons to their corresponding controller methods.
        """
        self.ui.analyze_button.clicked.connect(self.open_analysis_window)
        self.ui.pushButton.clicked.connect(self.clear_input)
        
    def clear_input(self) -> None:
        """
        Clear the RNA input text field.
        """
        self.ui.RNA_input.clear()
        
    def open_analysis_window(self) -> None:
        """
        Validate RNA input and open the analysis window if valid.

        If the sequence is invalid, an error message is displayed
        in the input field.
        """
        sequence: str = self.ui.RNA_input.toPlainText().strip()
        
        try:
            self.rna_sequence = RNASequence(sequence)
        except ValueError as error:
            self.ui.RNA_input.setText(str(error))
            return 
        
        self.analysis_window = QMainWindow()
        self.analysis_ui = Ui_Sequence_analysis_window()
        self.analysis_ui.setupUi(self.analysis_window)

        self.analysis_ui.METHOD_RESULTS_txt.setText("")
        self.analysis_ui.feedback_txt.setText("") if hasattr(self.analysis_ui, "feedback_txt") else None

        self.analysis_ui.run_methods_button.clicked.connect(self.run_selected_methods)
        self.analysis_ui.return_button.clicked.connect(self.return_to_input)
        self.analysis_ui.save_results_button.clicked.connect(self.handle_save_results)

        self.analysis_window.show()
        self.window.hide()
