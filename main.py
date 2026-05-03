import sys
from PyQt6.QtWidgets import QApplication, QMainWindow
from sequence_input_gui import Ui_MainWindow
from logic import LogicController


def main():
    app = QApplication(sys.argv)

    window = QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(window)

    controller = LogicController(ui, window)

    window.show()
    app.exec()


if __name__ == "__main__":
    main()