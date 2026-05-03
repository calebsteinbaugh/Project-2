from PyQt6.QtWidgets import QMainWindow
from sequence_analysis_gui import Ui_Sequence_analysis_window
from RNA_protein_Class import RNASequence


class LogicController:
    
    def __init__(self, ui, window):
        self.ui = ui
        self.window = window
        self.rna_sequence = None
        
        self.analysis_window = None
        self.analysis_ui = None
        self.current_results = ""
        
        self.connect_signals()
        
    def connect_signals(self):
        self.ui.analyze_button.clicked.connect(self.open_analysis_window)
        self.ui.pushButton.clicked.connect(self.clear_input)
        
    def clear_input(self):
        self.ui.RNA_input.clear()
        
    def open_analysis_window(self):
        sequence = self.ui.RNA_input.toPlainText().strip()
        
        try:
            self.rna_sequence = RNASequence(sequence)
        except ValueError as error:
            self.ui.RNA_input.setText(str(error))
            return 
        
        self.analysis_window = QMainWindow()
        self.analysis_ui = Ui_Sequence_analysis_window()
        self.analysis_ui.setupUi(self.analysis_window)
        
        self.analysis_ui.METHOD_RESULTS_txt.setText("")
        self.analysis_ui.feedback_txt.setText("")

        self.analysis_ui.run_methods_button.clicked.connect(self.run_selected_methods)
        self.analysis_ui.return_button.clicked.connect(self.return_to_input)
        self.analysis_ui.save_results.clicked.connect(self.handle_save_results)
        
        self.analysis_window.show()
        self.window.hide()
        
    def return_to_input(self):
        if self.analysis_window:
            self.analysis_window.close()
        self.window.show()
        
    def run_selected_methods(self):
        if not self.rna_sequence or not self.analysis_ui:
            return
        
        results = ""
        
        if self.analysis_ui.show_codon_checkbox.isChecked():
            codons = self.rna_sequence.get_codons()
            results += "Codons:\n"
            results += " | ".join(codons) + "\n\n"
            
        if self.analysis_ui.translate_codon_checkbox.isChecked():
            protein = self.rna_sequence.translate()
            results += "Translation:\n"
            results += "-".join(protein) + "\n\n"
            
        if self.analysis_ui.find_start_codon_checkbox.isChecked():
            start_codon = self.rna_sequence.get_start_codon()
            results += "Start Codon:\n"
            
            if start_codon == -1:
                results += "No Start Codon found.\n\n"
            else:
                results += f"AUG found at index {start_codon}.\n\n"

        if self.analysis_ui.find_stop_codon_checkbox.isChecked():
            stop_codon = self.rna_sequence.get_stop_codon()
            results += "Stop Codon:\n"

            if stop_codon == -1:
                results += "No Stop Codon found.\n\n"
            else:
                results += f"Stop codon found at index {stop_codon}.\n\n"
                
        if self.analysis_ui.amino_acid_count_checkbox.isChecked():
            amino_count = self.rna_sequence.amino_acids_count()
            results += "Amino Acid Counts:\n"
            
            if not amino_count:
                results += "No amino acids detected.\n\n"
            else:
                for amino, count in amino_count.items():
                    results += f"{amino}: {count}\n"
                results += "\n"
                
        if self.analysis_ui.product_prediction_checkbox.isChecked():
            prediction = self.rna_sequence.predict_protein_functionality()
            results += "Product Prediction:\n"
            results += f"{prediction}\n\n"
            
        if results == "":
            results = "Please select at least one analysis method."
            
        self.current_results = results
        self.analysis_ui.METHOD_RESULTS_txt.setText(results)
        self.analysis_ui.feedback_txt.setText("")
        
    def handle_save_results(self):
        if self.current_results == "":
            self.analysis_ui.feedback_txt.setText("No results to save.")
            return

        try:
            with open("analysis_results.txt", "a") as results_file:
                results_file.write("=== RNA Sequence Analysis ===\n")
                results_file.write(self.current_results)
                results_file.write("\n\n")

            self.analysis_ui.feedback_txt.setText("Results saved.")

        except OSError:
            self.analysis_ui.feedback_txt.setText("Could not save results.")