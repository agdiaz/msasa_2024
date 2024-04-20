import numpy as np

class SequenceHandler:
    def __init__(self, gap_char=45):
        self.gap_char = gap_char
        self.padding = {0: self.pad_left, 1: self.pad_right}

    def pad_left(self, sequence):
        return np.insert(sequence, 0, self.gap_char)

    def pad_right(self, sequence):
        return np.append(sequence, self.gap_char)

    def generate_new_sequences(self, sequences, num_sequences, seq_length, addition_deletion_prob=0.5):
        direction = np.random.randint(0, 2)
        padding_fx = self.padding[direction]

        gap_positions = np.argwhere(sequences == self.gap_char)

        if np.random.rand() < addition_deletion_prob or len(gap_positions) == 0:
            new_sequences = np.zeros([num_sequences, seq_length + 1], dtype=np.uint8)

            seq_index_to_edit = np.random.randint(num_sequences)
            col_index_to_edit = np.random.randint(seq_length)

            for current_sequence_index in range(num_sequences):
                if current_sequence_index == seq_index_to_edit:
                    temp_seq = np.insert(sequences[current_sequence_index], col_index_to_edit, self.gap_char)
                else:
                    temp_seq = padding_fx(sequences[current_sequence_index])

                new_sequences[current_sequence_index] = temp_seq
        else:
            new_sequences = np.zeros([num_sequences, seq_length], dtype=np.uint8)
            seq_index_to_edit, col_index_to_edit = gap_positions[np.random.choice(len(gap_positions))]
            for current_sequence_index in range(num_sequences):
                if current_sequence_index == seq_index_to_edit:
                    temp_seq = np.delete(sequences[current_sequence_index], col_index_to_edit)
                    new_sequences[current_sequence_index] = padding_fx(temp_seq)
                else:
                    new_sequences[current_sequence_index] = sequences[current_sequence_index]

        while np.all(new_sequences[:, 0] == self.gap_char) or np.all(new_sequences[:, -1] == self.gap_char):
            if np.all(new_sequences[:, 0] == self.gap_char):
                new_sequences = new_sequences[:, 1:]
            if np.all(new_sequences[:, -1] == self.gap_char):
                new_sequences = new_sequences[:, :-1]

        return new_sequences
