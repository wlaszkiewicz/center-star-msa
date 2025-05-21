import sys
from PyQt5.QtWidgets import QApplication
from gui.main_window import CenterStarApp

if __name__ == '__main__':
    """
    Initializes and runs the PyQt5 application.
    """
    app = QApplication(sys.argv)
    main_window = CenterStarApp()
    main_window.show()
    sys.exit(app.exec_())
