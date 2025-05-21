import sys
from PyQt5.QtCore import Qt, QPoint, QRect
from PyQt5.QtGui import QPixmap, QPainter, QRegion
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTextEdit, QPushButton, QFileDialog, QMessageBox,
    QSpinBox, QGroupBox, QFormLayout, QTableWidget, QTableWidgetItem,
    QSizePolicy, QTabWidget, QListWidget, QListWidgetItem, QApplication, QDialog
)

from gui.dialogs import AddSequenceDialog
from core.alignment import get_center_sequence, center_star
from core.fasta_parser import parse_fasta
from core.consensus import generate_consensus
from core.statistics import calculate_msa_statistics


class CenterStarApp(QMainWindow):
    """
    Main application window for the Center Star MSA Tool.
    Handles user interactions, sequence input, running alignments,
    and displaying results.
    """

    def __init__(self):
        """
        Initializes the main application window and its UI components.
        """
        super().__init__()
        self.setWindowTitle("Center Star Algorithm MSA Tool")
        screen = QApplication.primaryScreen()
        screen_geometry = screen.availableGeometry()
        window_width, window_height = 1800, 1300
        x = (screen_geometry.width() - window_width) // 2
        y = (screen_geometry.height() - window_height) // 2
        self.setGeometry(x, y, window_width, window_height)

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)

        self.results_available = False

        # --- Input Group ---
        input_group = QGroupBox("Input Sequences")
        input_main_layout = QVBoxLayout()
        self.sequence_list_widget = QListWidget()
        self.sequence_list_widget.setSelectionMode(QListWidget.ExtendedSelection)
        self.sequence_list_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.sequence_list_widget.setMinimumHeight(200)
        input_main_layout.addWidget(self.sequence_list_widget)

        input_buttons_layout = QHBoxLayout()
        self.load_fasta_button = QPushButton("Load FASTA")
        self.load_fasta_button.clicked.connect(self.load_fasta_file)
        self.add_manual_button = QPushButton("Add Manually")
        self.add_manual_button.clicked.connect(self.add_manual_sequence)
        self.load_test_data_button = QPushButton("Load Test Data")
        self.load_test_data_button.clicked.connect(self.load_test_sequences)
        self.remove_sequence_button = QPushButton("Remove Selected")
        self.remove_sequence_button.clicked.connect(self.remove_selected_sequences)
        self.clear_sequences_button = QPushButton("Clear All")
        self.clear_sequences_button.clicked.connect(self.clear_all_sequences)

        input_buttons_layout.addWidget(self.load_fasta_button)
        input_buttons_layout.addWidget(self.add_manual_button)
        input_buttons_layout.addWidget(self.load_test_data_button)
        input_buttons_layout.addWidget(self.remove_sequence_button)
        input_buttons_layout.addWidget(self.clear_sequences_button)
        input_main_layout.addLayout(input_buttons_layout)
        input_group.setLayout(input_main_layout)
        self.main_layout.addWidget(input_group)

        # --- Scoring Group ---
        scoring_group = QGroupBox("Scoring Scheme")
        scoring_layout = QFormLayout()
        self.match_score_input = QSpinBox()
        self.match_score_input.setRange(0, 100)
        self.match_score_input.setValue(1)
        self.mismatch_penalty_input = QSpinBox()
        self.mismatch_penalty_input.setRange(-100, 0)
        self.mismatch_penalty_input.setValue(-1)
        self.gap_penalty_input = QSpinBox()
        self.gap_penalty_input.setRange(-100, 0)
        self.gap_penalty_input.setValue(-2)
        scoring_layout.addRow("Match Score:", self.match_score_input)
        scoring_layout.addRow("Mismatch Penalty:", self.mismatch_penalty_input)
        scoring_layout.addRow("Gap Penalty:", self.gap_penalty_input)
        scoring_group.setLayout(scoring_layout)
        self.main_layout.addWidget(scoring_group)

        self.run_button = QPushButton("Generate Alignment and Statistics")
        self.run_button.setStyleSheet("QPushButton {font-size: 30px; padding: 10px;}")
        self.run_button.clicked.connect(self.run_alignment)
        self.main_layout.addWidget(self.run_button)

        # --- Output Tabs ---
        self.output_tabs = QTabWidget()
        self.output_tabs.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Tab 1: Alignment & Summary
        self.tab1_widget = QWidget()
        self.tab1_layout = QVBoxLayout(self.tab1_widget)
        self.results_display = QTextEdit()
        self.results_display.setReadOnly(True)
        self.results_display.setFontFamily("Courier New")
        self.results_display.setLineWrapMode(QTextEdit.NoWrap)
        self.results_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.tab1_layout.addWidget(self.results_display)
        self.output_tabs.addTab(self.tab1_widget, "📊 Alignment & Summary")

        # Tab: Center Election Matrix
        self.tab_center_election_widget = QWidget()
        self.tab_center_election_layout = QVBoxLayout(self.tab_center_election_widget)
        self.center_election_matrix_table_widget = QTableWidget()
        self.center_election_matrix_table_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.center_election_matrix_table_widget.setAlternatingRowColors(True)
        self.tab_center_election_layout.addWidget(self.center_election_matrix_table_widget)
        self.output_tabs.addTab(self.tab_center_election_widget, "⭐ Center Election Matrix")

        # Tab 2: Pairwise Identity Matrix
        self.tab2_widget = QWidget()
        self.tab2_layout = QVBoxLayout(self.tab2_widget)
        self.identity_matrix_table_widget = QTableWidget()
        self.identity_matrix_table_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.identity_matrix_table_widget.setAlternatingRowColors(True)
        self.tab2_layout.addWidget(self.identity_matrix_table_widget)
        self.output_tabs.addTab(self.tab2_widget, "🆔 Pairwise Identity Matrix (%)")

        # Tab 3: Pairwise Distance Matrix
        self.tab3_widget = QWidget()
        self.tab3_layout = QVBoxLayout(self.tab3_widget)
        self.distance_matrix_table_widget = QTableWidget()
        self.distance_matrix_table_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.distance_matrix_table_widget.setAlternatingRowColors(True)
        self.tab3_layout.addWidget(self.distance_matrix_table_widget)
        self.output_tabs.addTab(self.tab3_widget, "📐 Pairwise Distance Matrix")
        self.main_layout.addWidget(self.output_tabs)

        # --- Save Buttons ---
        save_buttons_layout = QHBoxLayout()
        self.save_text_file_button = QPushButton("💾 Save All Output as TXT")
        self.save_text_file_button.clicked.connect(self.save_output_text_file)

        self.save_text_summary_image_button = QPushButton("🖼️ Save Summary (Tab 1) as Image")
        self.save_text_summary_image_button.clicked.connect(self.save_text_summary_as_image)
        self.save_center_election_matrix_image_button = QPushButton("🖼️ Save Center Election Matrix as Image")
        self.save_center_election_matrix_image_button.clicked.connect(self.save_center_election_matrix_as_image)
        self.save_identity_matrix_image_button = QPushButton("🖼️ Save Identity Matrix as Image")
        self.save_identity_matrix_image_button.clicked.connect(self.save_identity_matrix_as_image)
        self.save_distance_matrix_image_button = QPushButton("🖼️ Save Distance Matrix as Image")
        self.save_distance_matrix_image_button.clicked.connect(self.save_distance_matrix_as_image)

        save_buttons_layout.addWidget(self.save_text_file_button)
        save_buttons_layout.addWidget(self.save_text_summary_image_button)
        save_buttons_layout.addWidget(self.save_center_election_matrix_image_button)
        save_buttons_layout.addWidget(self.save_identity_matrix_image_button)
        save_buttons_layout.addWidget(self.save_distance_matrix_image_button)
        self.main_layout.addLayout(save_buttons_layout)


        self.output_tabs.currentChanged.connect(self._on_tab_changed)

        self._set_results_availability(False)


        self.main_layout.setStretchFactor(input_group, 2)
        self.main_layout.setStretchFactor(scoring_group, 1)
        self.main_layout.setStretchFactor(self.output_tabs, 7)

        self.current_sequences_with_ids = []
        self.full_text_output_for_saving = ""

    def _set_results_availability(self, available: bool):
        """
        Sets the availability of results and updates UI elements accordingly.

        Args:
            available (bool): True if results are available, False otherwise.
        """
        self.results_available = available
        if not available:
            self.results_display.clear()
            self.identity_matrix_table_widget.clearContents()
            self.identity_matrix_table_widget.setRowCount(0)
            self.identity_matrix_table_widget.setColumnCount(0)
            self.distance_matrix_table_widget.clearContents()
            self.distance_matrix_table_widget.setRowCount(0)
            self.distance_matrix_table_widget.setColumnCount(0)
            self.center_election_matrix_table_widget.clearContents()
            self.center_election_matrix_table_widget.setRowCount(0)
            self.center_election_matrix_table_widget.setColumnCount(0)
            self.full_text_output_for_saving = ""


        self.save_text_file_button.setEnabled(self.results_available)

        self._update_image_save_button_visibility()

    def _update_image_save_button_visibility(self):
        """
        Updates the visibility of specific "Save as Image" buttons
        based on the currently active tab and whether results are available.
        """
        show_summary_btn = False
        show_center_election_btn = False
        show_identity_btn = False
        show_distance_btn = False

        if self.results_available:
            current_index = self.output_tabs.currentIndex()
            if current_index == 0:
                show_summary_btn = True
            elif current_index == 1:
                show_center_election_btn = True
            elif current_index == 2:
                show_identity_btn = True
            elif current_index == 3:
                show_distance_btn = True

        self.save_text_summary_image_button.setVisible(show_summary_btn)
        self.save_center_election_matrix_image_button.setVisible(show_center_election_btn)
        self.save_identity_matrix_image_button.setVisible(show_identity_btn)
        self.save_distance_matrix_image_button.setVisible(show_distance_btn)

    def _on_tab_changed(self, index: int):
        """
        Slot for the QTabWidget's currentChanged signal.
        Updates save button visibility when the tab changes.

        Args:
            index (int): The index of the newly selected tab.
        """
        self._update_image_save_button_visibility()

    def _update_sequence_list_widget(self):
        """
        Refreshes the QListWidget displaying the current input sequences.
        """
        self.sequence_list_widget.clear()
        for seq_id, seq_data in self.current_sequences_with_ids:
            preview = seq_data[:30] + "..." if len(seq_data) > 30 else seq_data
            item_text = f"{seq_id} ({len(seq_data)}bp): {preview}"
            item = QListWidgetItem(item_text)
            item.setData(Qt.UserRole, (seq_id, seq_data))

            self.sequence_list_widget.addItem(item)

    def _get_existing_ids(self) -> list[str]:
        """
        Returns a list of IDs of currently loaded sequences.

        Returns:
            list[str]: A list of sequence IDs.
        """
        return [item[0] for item in self.current_sequences_with_ids]

    def load_fasta_file(self):
        """
        Opens a file dialog to load sequences from a FASTA file.
        Appends new sequences to the existing list and handles potential ID conflicts.
        """
        filepaths, _ = QFileDialog.getOpenFileNames(self, "Open FASTA File(s)", "",
                                                    "FASTA files (*.fasta *.fa *.fna *.faa);;All files (*)")
        if not filepaths:
            return

        sequences_loaded_count = 0
        sequences_appended_count = 0
        newly_parsed_sequences = []

        for filepath in filepaths:
            try:
                with open(filepath, 'r', encoding='utf-8') as f:
                    file_content = f.read()

                parsed_sequences_from_file = parse_fasta(file_content)

                if not parsed_sequences_from_file:
                    QMessageBox.warning(self, "Warning", f"No sequences found or parsed from {filepath}.")
                    continue

                sequences_loaded_count += len(parsed_sequences_from_file)

                existing_ids = self._get_existing_ids()
                temp_sequences_to_add = []
                ids_in_current_file_load = set()

                for seq_id, seq_data in parsed_sequences_from_file:
                    original_seq_id = seq_id
                    counter = 1
                    while seq_id in existing_ids or seq_id in ids_in_current_file_load:
                        seq_id = f"{original_seq_id}_{counter}"
                        counter += 1

                    if seq_id != original_seq_id:
                        QMessageBox.information(self, "ID Conflict Resolved",
                                                f"Sequence ID '{original_seq_id}' from {filepath} already exists or was duplicated in this load.\n"
                                                f"It has been renamed to '{seq_id}'.")

                    temp_sequences_to_add.append((seq_id, seq_data))
                    ids_in_current_file_load.add(seq_id)

                newly_parsed_sequences.extend(temp_sequences_to_add)

            except Exception as e:
                QMessageBox.critical(self, "Error", f"An error occurred while parsing {filepath}:\n{e}")


        if newly_parsed_sequences:
            self._set_results_availability(False)
            self.current_sequences_with_ids.extend(newly_parsed_sequences)
            self._update_sequence_list_widget()
            sequences_appended_count = len(newly_parsed_sequences)
            QMessageBox.information(self, "Success",
                                    f"Successfully appended {sequences_appended_count} new sequences from the selected file(s).\n"
                                    f"Total sequences processed from files: {sequences_loaded_count}.")
        elif sequences_loaded_count == 0 and filepaths:
            QMessageBox.information(self, "Info", "No new sequences were added from the selected file(s).")

    def load_test_sequences(self):
        """
        Loads a predefined set of test sequences.
        Asks for confirmation if sequences are already loaded.
        """
        if self.current_sequences_with_ids:
            reply = QMessageBox.question(self, "Confirm Load Test Data",
                                         "This will clear current sequences and load predefined test data. Continue?",
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return

        self._set_results_availability(False)
        test_sequences = [
            ("s1_test", "ATTGCCATT"),
            ("s2_test", "ATGGCCATT"),
            ("s3_test", "ATCCATTTTT"),
            ("s4_test", "ATCTTCTT"),
            ("s5_test", "ACTGACC")
        ]
        self.current_sequences_with_ids = list(test_sequences)
        self._update_sequence_list_widget()
        QMessageBox.information(self, "Test Data Loaded",
                                f"{len(self.current_sequences_with_ids)} test sequences have been loaded.")

    def add_manual_sequence(self):
        """
        Opens the AddSequenceDialog to manually add a new sequence.
        Invalidates current results if a sequence is added.
        """
        dialog = AddSequenceDialog(existing_ids=self._get_existing_ids(), parent=self)
        if dialog.exec_() == QDialog.Accepted:
            seq_id, seq_data = dialog.get_data()
            if seq_id and seq_data:
                self._set_results_availability(False)
                self.current_sequences_with_ids.append((seq_id, seq_data))
                self._update_sequence_list_widget()

                self.sequence_list_widget.setCurrentRow(self.sequence_list_widget.count() - 1)

    def remove_selected_sequences(self):
        """
        Removes selected sequences from the input list.
        Invalidates current results.
        """
        selected_items = self.sequence_list_widget.selectedItems()
        if not selected_items:
            QMessageBox.information(self, "Info", "No sequences selected to remove.")
            return

        ids_to_remove = set()
        for item in selected_items:
            data = item.data(Qt.UserRole)
            if data:
                ids_to_remove.add(data[0])


        self.current_sequences_with_ids = [
            item_tuple for item_tuple in self.current_sequences_with_ids
            if item_tuple[0] not in ids_to_remove
        ]
        self._update_sequence_list_widget()
        self._set_results_availability(False)

    def clear_all_sequences(self):
        """
        Clears all sequences from the input list after confirmation.
        Invalidates current results.
        """
        if not self.current_sequences_with_ids:
            QMessageBox.information(self, "Info", "Sequence list is already empty.")
            return
        reply = QMessageBox.question(self, "Confirm Clear",
                                     "Are you sure you want to clear all loaded sequences?",
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.current_sequences_with_ids = []
            self._update_sequence_list_widget()
            self._set_results_availability(False)

    def _format_center_election_matrix_for_text_save(self, seq_ids_ordered: list[str],
                                                     score_matrix: list[list[float]],
                                                     total_scores: list[float],
                                                     center_idx_val: int,
                                                     max_id_len: int) -> str:
        """Formats the center election matrix for inclusion in the text save file."""
        output = "\n\n--- Center Election Score Matrix ---\n"
        output += "This matrix shows the alignment scores between all pairs of your input sequences. "
        output += "The Center Star algorithm picks the sequence with the highest 'Sum' score as the 'center' to guide the multiple alignment.\n\n"

        num_seqs = len(seq_ids_ordered)
        if num_seqs == 0:
            return output + "  (No sequences for matrix)\n"

        score_example_len = 7
        max_name_len = max(len(sid) for sid in seq_ids_ordered) if seq_ids_ordered else 0
        col_width = max(score_example_len, max_name_len, len("Sum")) + 2

        header_str = f"{'':<{max_id_len + 2}}"
        for seq_id in seq_ids_ordered:
            header_str += f"{seq_id:<{col_width}}"
        header_str += f"| {'Sum':<{col_width}}\n"
        output += header_str
        output += "-" * len(header_str.strip()) + "\n"


        for i in range(num_seqs):
            row_str = f"{seq_ids_ordered[i]:<{max_id_len + 2}}"
            for j in range(num_seqs):
                if i == j:
                    row_str += f"{'-':<{col_width}}"
                else:
                    row_str += f"{score_matrix[i][j]:<{col_width}.1f}"
            row_str += f"| {total_scores[i]:<{col_width}.1f}\n"
            output += row_str

        output += f"\nChosen Center Sequence: {seq_ids_ordered[center_idx_val]} (Total Sum Score: {total_scores[center_idx_val]:.1f})\n"
        return output

    def populate_center_election_table(self, table_widget: QTableWidget,
                                       seq_ids_ordered: list[str],
                                       score_matrix: list[list[float]],
                                       total_scores: list[float],
                                       center_idx_val: int):
        """Populates the QTableWidget for the center election matrix."""
        num_seqs = len(seq_ids_ordered)
        if num_seqs == 0:
            table_widget.clearContents()
            table_widget.setRowCount(0)
            table_widget.setColumnCount(0)
            return

        table_widget.setRowCount(num_seqs)
        table_widget.setColumnCount(num_seqs + 1)

        header_labels = list(seq_ids_ordered) + ["Sum"]
        table_widget.setHorizontalHeaderLabels(header_labels)
        table_widget.setVerticalHeaderLabels(seq_ids_ordered)

        for i in range(num_seqs):
            for j in range(num_seqs):  # Pairwise scores
                item_text = ""
                if i == j:
                    item_text = "-"  # Diagonal
                else:
                    item_text = f"{score_matrix[i][j]:.1f}"
                item = QTableWidgetItem(item_text)
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                table_widget.setItem(i, j, item)


            sum_item = QTableWidgetItem(f"{total_scores[i]:.1f}")
            sum_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            sum_item.setFlags(sum_item.flags() & ~Qt.ItemIsEditable)
            if i == center_idx_val:  #  the chosen center's sum
                font = sum_item.font()
                font.setBold(True)
                sum_item.setFont(font)
            table_widget.setItem(i, num_seqs, sum_item)

        table_widget.resizeColumnsToContents()
        table_widget.resizeRowsToContents()

    def run_alignment(self):
        """
        Runs the Center Star alignment and statistics calculation process.
        Updates the GUI with the results, using plain text for the summary tab.
        """
        if not self.current_sequences_with_ids:
            QMessageBox.warning(self, "Input Error", "Please add at least one sequence to align.")
            return

        match_score = self.match_score_input.value()
        mismatch_penalty = self.mismatch_penalty_input.value()
        gap_penalty = self.gap_penalty_input.value()

        self.results_display.setText("⏳ Processing alignment and statistics, please wait...")
        self._set_results_availability(False)
        QApplication.processEvents()

        try:
            center_idx_val, center_item_tuple, _, all_seq_items_for_election, \
                pairwise_score_matrix_for_election, total_scores_for_election = \
                get_center_sequence(self.current_sequences_with_ids, match_score, mismatch_penalty, gap_penalty)

            center_seq_id_for_display = center_item_tuple[0] if center_item_tuple else "N/A"
            election_seq_ids = [item[0] for item in all_seq_items_for_election]
            max_id_len_election = max(len(sid) for sid in election_seq_ids) if election_seq_ids else 0

            if len(self.current_sequences_with_ids) > 1:
                self.populate_center_election_table(self.center_election_matrix_table_widget,
                                                    election_seq_ids, pairwise_score_matrix_for_election,
                                                    total_scores_for_election, center_idx_val)

            aligned_sequences_str = []
            if len(self.current_sequences_with_ids) == 1:
                aligned_sequences_str = [self.current_sequences_with_ids[0][1]]
            elif len(self.current_sequences_with_ids) > 1:
                aligned_sequences_str = center_star(self.current_sequences_with_ids, match_score, mismatch_penalty,
                                                    gap_penalty)
            else:
                self._set_results_availability(False);
                return

            if not aligned_sequences_str or \
                    not all(s is not None and isinstance(s, str) for s in aligned_sequences_str) or \
                    len(aligned_sequences_str) != len(self.current_sequences_with_ids):
                msg = "Alignment process failed or returned an unexpected result structure."
                self.results_display.setText(msg)
                QMessageBox.critical(self, "Alignment Error", msg)
                self._set_results_availability(False);
                return

            avg_pid, total_gaps_msa, cons_cols, sp_score, id_mat, dist_mat, num_aligned_seqs, aln_len = \
                calculate_msa_statistics(aligned_sequences_str, match_score, mismatch_penalty, gap_penalty)
            consensus_seq = generate_consensus(aligned_sequences_str) if aln_len > 0 else ""

            # --- Output for Tab 1 ---
            text_output_tab1 = "--- Alignment Parameters Used ---\n"
            text_output_tab1 += f"  - Match Score: {match_score}\n"
            text_output_tab1 += f"  - Mismatch Penalty: {mismatch_penalty}\n"
            text_output_tab1 += f"  - Gap Penalty: {gap_penalty}\n"
            if len(self.current_sequences_with_ids) > 1 and center_idx_val != -1:
                text_output_tab1 += f"\n-> Chosen Center Sequence (for MSA): {election_seq_ids[center_idx_val]} "
                text_output_tab1 += f"(Sum of pairwise scores: {total_scores_for_election[center_idx_val]:.1f})\n"
                text_output_tab1 += "   (See 'Center Election Matrix' tab for details)\n"

            text_output_tab1 += "\n\n--- Aligned Sequences & Per-Sequence Stats ---\n"
            if len(self.current_sequences_with_ids) > 1 and center_seq_id_for_display != "N/A":
                text_output_tab1 += f"The sequence '{center_seq_id_for_display}' (marked with '* ' below) was used as the initial center.\n\n"
            else:
                text_output_tab1 += "\n"

            current_sequence_ids_ordered = [item[0] for item in self.current_sequences_with_ids]
            max_id_len_msa = max(
                len(sid) for sid in current_sequence_ids_ordered) if current_sequence_ids_ordered else 10
            max_id_len_msa = max(max_id_len_msa, len("Consensus"))

            aln_display_len = aln_len if aln_len > 0 else 18
            header_parts = [
                f"{'ID':<{max_id_len_msa + 2}}", f"{'Aligned Sequence':<{aln_display_len + 2}}",
                f"{'Gaps':<6}", f"{'MatchCons':<10}", f"{'MismatchCons':<13}", f"{'ID%Cons':<8}"
            ]
            header_line = "  " + "  ".join(header_parts)
            separator_line = "  " + "-" * (len(header_line) - 2)

            text_output_tab1 += header_line + "\n"
            text_output_tab1 += separator_line + "\n"

            if consensus_seq:
                consensus_display_id = f"{'Consensus':<{max_id_len_msa + 2}}"
                consensus_gaps = consensus_seq.count('-')
                text_output_tab1 += f"{consensus_display_id}  {consensus_seq:<{aln_display_len + 2}} {consensus_gaps:<6} {'N/A':<10} {'N/A':<13} {'100.0%':<8}\n"
                text_output_tab1 += separator_line + "\n"

            for i, seq_str_aligned in enumerate(aligned_sequences_str):
                seq_id_original = current_sequence_ids_ordered[i]
                prefix = "* " if seq_id_original == center_seq_id_for_display and len(
                    self.current_sequences_with_ids) > 1 else "  "
                display_id = f"{prefix}{seq_id_original:<{max_id_len_msa - len(prefix) + len('  ')}}"
                gaps_in_seq = seq_str_aligned.count('-')
                matches_vs_cons, mismatches_vs_cons, id_vs_cons_val = 0, 0, 0.0
                if consensus_seq and aln_len > 0:
                    valid_pos_for_id_vs_cons = 0
                    for k in range(aln_len):
                        s_char, c_char = seq_str_aligned[k], consensus_seq[k]
                        if c_char != '-':
                            valid_pos_for_id_vs_cons += 1
                            if s_char == c_char:
                                matches_vs_cons += 1
                            elif s_char != '-':
                                mismatches_vs_cons += 1
                    if valid_pos_for_id_vs_cons > 0:
                        id_vs_cons_val = (matches_vs_cons / valid_pos_for_id_vs_cons) * 100
                line_parts_data = [
                    f"{display_id}", f"{seq_str_aligned:<{aln_display_len + 2}}",
                    f"{gaps_in_seq:<6}", f"{matches_vs_cons:<10}", f"{mismatches_vs_cons:<13}",
                    f"{id_vs_cons_val:<8.2f}%"
                ]
                text_output_tab1 += "  ".join(line_parts_data) + "\n"

            text_output_tab1 += "\nPer-Sequence Stats Explanation (vs Consensus):\n"
            text_output_tab1 += "  - Gaps: Number of '-' characters in this specific aligned sequence.\n"
            text_output_tab1 += "  - MatchCons: Number of characters matching 'Consensus' (Consensus not gap).\n"
            text_output_tab1 += "  - MismatchCons: Number of characters differing from 'Consensus' (neither are gaps, Consensus not gap).\n"
            text_output_tab1 += "  - ID%Cons: Percent identity to 'Consensus' (Consensus not gap positions).\n"

            text_output_tab1 += "\n\n--- Overall MSA Summary Statistics ---\n"
            text_output_tab1 += f"  - Number of Sequences in MSA: {num_aligned_seqs}\n"
            text_output_tab1 += f"  - Alignment Length (columns): {aln_len}\n"
            text_output_tab1 += f"  - Average Pairwise Identity (PID): {avg_pid:.2f}%\n"
            text_output_tab1 += f"  - Total Gaps in MSA: {total_gaps_msa}\n"
            text_output_tab1 += f"  - Fully Conserved Columns (non-gap): {cons_cols}\n"
            text_output_tab1 += f"  - Sum-of-Pairs Score (SP Score): {sp_score}\n"
            if consensus_seq:
                text_output_tab1 += f"  - Generated Consensus Sequence Length: {len(consensus_seq)}\n"

            self.results_display.setText(text_output_tab1)


            self.populate_matrix_table(self.identity_matrix_table_widget, id_mat, current_sequence_ids_ordered,
                                       "Identity")
            self.populate_matrix_table(self.distance_matrix_table_widget, dist_mat, current_sequence_ids_ordered,
                                       "Distance")

            self.full_text_output_for_saving = "--- Alignment Parameters Used ---\n"
            self.full_text_output_for_saving += f"Match Score: {match_score}, Mismatch Penalty: {mismatch_penalty}, Gap Penalty: {gap_penalty}\n"
            if len(self.current_sequences_with_ids) > 1:
                self.full_text_output_for_saving += self._format_center_election_matrix_for_text_save(
                    election_seq_ids, pairwise_score_matrix_for_election, total_scores_for_election, center_idx_val,
                    max_id_len_election)
            self.full_text_output_for_saving += "\n\n--- Aligned Sequences & Per-Sequence Stats ---\n"
            self.full_text_output_for_saving += header_line + "\n" + separator_line + "\n"
            if consensus_seq:
                consensus_display_id_plain = f"{'Consensus':<{max_id_len_msa + 2}}"
                consensus_gaps_plain = consensus_seq.count('-')
                self.full_text_output_for_saving += f"{consensus_display_id_plain}  {consensus_seq:<{aln_display_len + 2}} {consensus_gaps_plain:<6} {'N/A':<10} {'N/A':<13} {'100.0%':<8}\n"
                self.full_text_output_for_saving += separator_line + "\n"
            for i, seq_str_aligned in enumerate(aligned_sequences_str):
                seq_id_original = current_sequence_ids_ordered[i]
                prefix = "* " if seq_id_original == center_seq_id_for_display and len(
                    self.current_sequences_with_ids) > 1 else "  "
                display_id_plain = f"{prefix}{seq_id_original:<{max_id_len_msa - len(prefix) + len('  ')}}"
                gaps_in_seq_plain = seq_str_aligned.count('-')
                matches_vs_cons_plain, mismatches_vs_cons_plain, id_vs_cons_val_plain = 0, 0, 0.0
                if consensus_seq and aln_len > 0:
                    valid_pos_plain = 0
                    for k in range(aln_len):
                        s, c = seq_str_aligned[k], consensus_seq[k]
                        if c != '-': valid_pos_plain += 1
                        if s == c and c != '-':
                            matches_vs_cons_plain += 1
                        elif s != '-' and c != '-' and s != c:
                            mismatches_vs_cons_plain += 1
                    if valid_pos_plain > 0: id_vs_cons_val_plain = (matches_vs_cons_plain / valid_pos_plain) * 100
                line_parts_data_plain = [
                    f"{display_id_plain}", f"{seq_str_aligned:<{aln_display_len + 2}}",
                    f"{gaps_in_seq_plain:<6}", f"{matches_vs_cons_plain:<10}", f"{mismatches_vs_cons_plain:<13}",
                    f"{id_vs_cons_val_plain:<8.2f}%"
                ]
                self.full_text_output_for_saving += "  ".join(line_parts_data_plain) + "\n"

            self.full_text_output_for_saving += "\n\n--- Overall MSA Summary Statistics ---\n"
            self.full_text_output_for_saving += f"Number of Sequences in MSA: {num_aligned_seqs}\nAlignment Length: {aln_len}\n"
            self.full_text_output_for_saving += f"Average Pairwise Identity: {avg_pid:.2f}%\nTotal Gaps in MSA: {total_gaps_msa}\n"
            self.full_text_output_for_saving += f"Fully Conserved Columns: {cons_cols}\nSum-of-Pairs Score: {sp_score}\n"
            if consensus_seq: self.full_text_output_for_saving += f"Generated Consensus Length: {len(consensus_seq)}\n"
            self.full_text_output_for_saving += self._format_matrix_for_text_save(id_mat, current_sequence_ids_ordered,
                                                                                  "Pairwise Identity Matrix (%)",
                                                                                  max_id_len_msa)
            self.full_text_output_for_saving += self._format_matrix_for_text_save(dist_mat,
                                                                                  current_sequence_ids_ordered,
                                                                                  "Pairwise Distance Matrix",
                                                                                  max_id_len_msa)

            self.output_tabs.setCurrentIndex(0)
            self._set_results_availability(True)

        except Exception as e:
            self.results_display.setText(f"An error occurred during processing:\n{e}\n\n"
                                         "Please check your input sequences and scoring parameters.")
            import traceback
            traceback.print_exc()
            QMessageBox.critical(self, "Processing Error", f"An unexpected error occurred: {e}")
            self._set_results_availability(False)

    def _format_matrix_for_text_save(self, matrix_data: list[list[float]],
                                     header_labels: list[str], title: str,
                                     max_label_len: int) -> str:
        """
        Formats a given matrix (e.g., identity, distance) for text file saving.

        Args:
            matrix_data: The 2D list of floats representing the matrix.
            header_labels: List of sequence IDs for rows/columns.
            title: Title for this matrix section.
            max_label_len: Maximum length of a sequence ID for formatting.

        Returns:
            A string representation of the formatted matrix.
        """
        if not matrix_data or not header_labels or (matrix_data and not matrix_data[0]):
            return f"\n\n--- {title} ---\n  (Not applicable or empty)\n"

        num_seqs = len(header_labels)

        if num_seqs > 0 and (len(matrix_data) != num_seqs or len(matrix_data[0]) != num_seqs):
            return f"\n\n--- {title} ---\n  (Matrix data inconsistent with labels for {num_seqs} sequences)\n"
        elif num_seqs == 0 and matrix_data:
            return f"\n\n--- {title} ---\n  (No labels provided for matrix data)\n"

        output = f"\n\n--- {title} ---\n"

        val_example_len = 7
        col_width = max(val_example_len,
                        max(len(lbl) for lbl in header_labels) if header_labels else 0) + 2


        output += f"{'':<{max_label_len + 2}}"
        for lbl in header_labels:
            output += f"{lbl:<{col_width}}"
        output += "\n"
        output += "-" * ((max_label_len + 2) + num_seqs * col_width) + "\n"


        for i, row_label in enumerate(header_labels):
            row_str = f"{row_label:<{max_label_len + 2}}"
            for val_idx, val in enumerate(matrix_data[i]):
                row_str += f"{val:<{col_width - 1}.2f} "
            output += row_str.strip() + "\n"
        return output

    def populate_matrix_table(self, table_widget: QTableWidget,
                              matrix_data: list[list[float]],
                              header_labels: list[str],
                              matrix_type_str: str):
        """
        Populates a QTableWidget with matrix data (identity or distance).

        Args:
            table_widget: The QTableWidget to populate.
            matrix_data: The 2D list of floats (matrix values).
            header_labels: List of sequence IDs for headers.
            matrix_type_str: String like "Identity" or "Distance" for tooltips.
        """
        num_seqs = len(header_labels)


        if num_seqs == 0 or not matrix_data or \
                (matrix_data and not matrix_data[0] and num_seqs > 0) or \
                (matrix_data and matrix_data[0] and (len(matrix_data) != num_seqs or len(matrix_data[0]) != num_seqs)):
            table_widget.clearContents()
            table_widget.setRowCount(0)
            table_widget.setColumnCount(0)
            if num_seqs > 0:
                QMessageBox.warning(self, "Matrix Display Error",
                                    f"The {matrix_type_str} matrix data is inconsistent or empty, "
                                    f"cannot display for {num_seqs} sequences.")
            return

        table_widget.setRowCount(num_seqs)
        table_widget.setColumnCount(num_seqs)
        table_widget.setHorizontalHeaderLabels(header_labels)
        table_widget.setVerticalHeaderLabels(header_labels)

        for i in range(num_seqs):
            for j in range(num_seqs):
                value = matrix_data[i][j]
                item = QTableWidgetItem(f"{value:.2f}")
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                tooltip_suffix = "%" if matrix_type_str == "Identity" else ""
                item.setToolTip(
                    f"{matrix_type_str} between '{header_labels[i]}' and '{header_labels[j]}': {value:.2f}{tooltip_suffix}")
                table_widget.setItem(i, j, item)

        table_widget.resizeColumnsToContents()
        table_widget.resizeRowsToContents()

    def save_output_text_file(self):
        """
        Saves the consolidated text output (parameters, alignment, stats, matrices)
        to a user-specified TXT file.
        """
        if not self.full_text_output_for_saving:
            QMessageBox.warning(self, "No Output to Save", "Please run an alignment first to generate output.")
            return

        filepath, _ = QFileDialog.getSaveFileName(self, "Save All Output as TXT", "output/center_star_msa_output.txt",
                                                  "Text Files (*.txt);;All Files (*)")
        if filepath:
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write(self.full_text_output_for_saving)
                QMessageBox.information(self, "Success", f"Output successfully saved to {filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Error Saving TXT File", f"Could not save the file:\n{e}")

    def save_widget_as_image(self, widget: QWidget, dialog_title: str = "Save as Image",
                             default_filename: str = "output/output.png"):
        """
        Saves the visual content of a given QWidget (e.g., QTableWidget, QTextEdit) as an image file.

        Args:
            widget: The QWidget to capture.
            dialog_title: Title for the save file dialog.
            default_filename: Default filename for the saved image.
        """

        if isinstance(widget, QTableWidget) and (widget.rowCount() == 0 or widget.columnCount() == 0):
            QMessageBox.warning(self, "Cannot Save Image",
                                f"The table for '{dialog_title}' is empty or has no data to save.")
            return
        if isinstance(widget, QTextEdit) and not widget.toPlainText().strip():
            QMessageBox.warning(self, "Cannot Save Image",
                                f"The content for '{dialog_title}' is empty.")
            return

        filepath, _ = QFileDialog.getSaveFileName(self, dialog_title, default_filename,
                                                  "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg);;All Files (*)")
        if not filepath:
            return

        pixmap_to_save = None
        original_parent = None
        original_pos = None
        original_flags = None
        layout_management_disabled = False

        try:
            QApplication.processEvents()

            if isinstance(widget, QTableWidget):
                h_header = widget.horizontalHeader()
                v_header = widget.verticalHeader()

                render_width = sum(
                    widget.columnWidth(i) for i in range(widget.columnCount()) if not widget.isColumnHidden(i))
                if v_header.isVisible() and not v_header.isHidden():
                    render_width += v_header.width()

                render_height = sum(widget.rowHeight(i) for i in range(widget.rowCount()) if not widget.isRowHidden(i))
                if h_header.isVisible() and not h_header.isHidden():
                    render_height += h_header.height()

                render_width += widget.frameWidth() * 2
                render_height += widget.frameWidth() * 2

                if render_width <= 0 or render_height <= 0:
                    QMessageBox.warning(self, "Image Creation Error",
                                        f"Table '{dialog_title}' has zero or negative calculated dimensions for rendering.")
                    return

                original_geometry = widget.geometry()
                original_min_size = widget.minimumSize()
                original_max_size = widget.maximumSize()
                original_h_policy = widget.horizontalScrollBarPolicy()
                original_v_policy = widget.verticalScrollBarPolicy()
                original_size_policy = widget.sizePolicy()

                parent_layout_obj = widget.parentWidget().layout() if widget.parentWidget() else None
                if parent_layout_obj:

                    if hasattr(parent_layout_obj, 'setEnabled'):
                        parent_layout_obj.setEnabled(False)
                        layout_management_disabled = True

                widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
                widget.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

                widget.setMinimumSize(render_width, render_height)
                widget.setMaximumSize(render_width, render_height)
                widget.resize(render_width, render_height)
                widget.updateGeometry()

                QApplication.processEvents()
                widget.update()
                QApplication.processEvents()


                pixmap_to_save = QPixmap(render_width, render_height)
                pixmap_to_save.fill(widget.palette().window().color())


                painter = QPainter(pixmap_to_save)

                widget.render(painter, QPoint(), QRegion(widget.rect()),
                              QWidget.RenderFlags(QWidget.DrawWindowBackground | QWidget.DrawChildren))
                painter.end()


                widget.setMinimumSize(original_min_size)
                widget.setMaximumSize(original_max_size)
                widget.setSizePolicy(original_size_policy)
                widget.setHorizontalScrollBarPolicy(original_h_policy)
                widget.setVerticalScrollBarPolicy(original_v_policy)

                if layout_management_disabled and parent_layout_obj:
                    parent_layout_obj.setEnabled(True)

                else:
                    widget.setGeometry(original_geometry)

                widget.updateGeometry()
                QApplication.processEvents()
                widget.update()
                QApplication.processEvents()


            elif isinstance(widget, QTextEdit):
                doc = widget.document()
                QApplication.processEvents()

                doc_size = doc.size()
                image_width = int(doc_size.width())
                image_height = int(doc_size.height())

                if image_width <= 0 or image_height <= 0:
                    QMessageBox.warning(self, "Image Creation Error",
                                        f"Document content for '{dialog_title}' has no size. Ensure there is text.")
                    return

                pixmap_to_save = QPixmap(image_width, image_height)
                pixmap_to_save.fill(widget.palette().base().color())

                painter = QPainter(pixmap_to_save)
                doc.drawContents(painter)
                painter.end()

            else:
                if not widget.isVisible() or widget.width() <= 1 or widget.height() <= 1:
                    QMessageBox.warning(self, "Cannot Save Image",
                                        f"The content for '{dialog_title}' is not visible or has no size to grab.")
                    return
                pixmap_to_save = widget.grab()

            if pixmap_to_save is None or pixmap_to_save.isNull():
                QMessageBox.critical(self, "Image Creation Error",
                                     "Failed to create image from widget content. The pixmap is null.")
                return


            file_format_str = "PNG"
            if filepath.lower().endswith(".jpg") or filepath.lower().endswith(".jpeg"):
                file_format_str = "JPG"

            quality = -1 if file_format_str == "PNG" else 90

            if not pixmap_to_save.save(filepath, file_format_str, quality):
                QMessageBox.critical(self, "Image Save Error",
                                     f"Failed to save image as {file_format_str} to {filepath}.")
            else:
                QMessageBox.information(self, "Image Saved", f"Image successfully saved to {filepath}")

        except Exception as e:
            QMessageBox.critical(self, "Image Save Exception",
                                 f"An error occurred while preparing or saving the image:\n{e}")
            import traceback
            traceback.print_exc()
        finally:
            if layout_management_disabled and parent_layout_obj and hasattr(parent_layout_obj, 'setEnabled'):
                if not parent_layout_obj.isEnabled():
                    parent_layout_obj.setEnabled(True)
                    if widget.parentWidget(): widget.parentWidget().updateGeometry()
                    QApplication.processEvents()

    def save_text_summary_as_image(self):
        """Saves the content of the 'Alignment & Summary' QTextEdit (Tab 1) as an image."""
        if not self.results_display.toPlainText().strip():
            QMessageBox.warning(self, "No Text Output",
                                "The summary text output (Tab 1) is empty. Cannot save as image.")
            return
        self.save_widget_as_image(self.results_display, "Save Summary (Tab 1) as Image", "output/msa_summary_alignment.png")

    def save_center_election_matrix_as_image(self):
        """Saves the 'Center Election Matrix' QTableWidget as an image."""
        self.save_widget_as_image(self.center_election_matrix_table_widget, "Save Center Election Matrix as Image",
                                  "output/center_election_matrix.png")

    def save_identity_matrix_as_image(self):
        """Saves the 'Pairwise Identity Matrix' QTableWidget as an image."""
        self.save_widget_as_image(self.identity_matrix_table_widget, "Save Pairwise Identity Matrix as Image",
                                  "output/identity_matrix.png")

    def save_distance_matrix_as_image(self):
        """Saves the 'Pairwise Distance Matrix' QTableWidget as an image."""
        self.save_widget_as_image(self.distance_matrix_table_widget, "Save Pairwise Distance Matrix as Image",
                                  "output/distance_matrix.png")