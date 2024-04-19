from concurrent.futures import ThreadPoolExecutor
import numpy as np
from methodtools import lru_cache
from random import choice
from .objective_functions import ObjectiveFunction
from .sequences_generator import SequenceHandler

class SimulatedAnnealing:
    def __init__(
        self,
        sequences,
        logger,
        quality_function,
        extend=False,
        temp=100.0,
        cooling_rate=0.95,
        min_temp=0.001,
        no_changes_limit = 100,
        changes = 1,
        iteration_neighbors = 1,
    ):
        self.sequences = sequences
        self.extend_sequences = extend
        self.temp = temp
        self.cooling_rate = cooling_rate
        self.min_temp = min_temp
        self.no_changes_limit = no_changes_limit
        self.logger = logger
        self.changes = changes
        self.iteration_neighbors = iteration_neighbors

        self.quality_function: ObjectiveFunction = quality_function
        self.sequence_generator = SequenceHandler(ObjectiveFunction.GAP_SYMBOL)

    def log_async(
        self,
        iteration: int,
        temp: float,
        sequence_length: int,
        current_score: float,
        iteration_score: float,
        new_score: float,
        score_change: float,
        initial_score: float,
        best_score: float,
        total_accepted: int,
        total_rejected: int,
        accepted: bool,
        no_changes: int,
    ):
        try:
            if score_change > 0:
                acceptance = 1
            else:
                acceptance = np.exp(score_change / temp)

            historical_score_improvement: float = best_score - initial_score
            if initial_score != 0.0:
                historical_improvement_percentage: float = (historical_score_improvement / abs(initial_score)) * 100
                current_improvement_percentage: float = (score_change / abs(initial_score)) * 100
            else:
                historical_improvement_percentage: float = 0
                current_improvement_percentage: float = 0

            row_values = f"{iteration}\t{sequence_length}\t{temp:.10f}\t{current_score:+.10f}\t{iteration_score:+.10f}\t{new_score:+.10f}\t{score_change:+.10f}\t{current_improvement_percentage:+.3f}\t{best_score:+.10f}\t{historical_score_improvement:+.10f}\t{historical_improvement_percentage:+.3f}\t{total_accepted}\t{total_rejected}\t{accepted}\t{no_changes}\t{acceptance:.10f}"
            self.logger.info(row_values)
            print(row_values)
        except:
            pass

    @lru_cache(maxsize=204800)
    def cached_energy(self, current_sequences_bytes, num_sequences, seq_length) -> float:
        sequences = np.frombuffer(current_sequences_bytes, dtype=np.uint8).reshape((num_sequences, seq_length))
        current_score_columns = np.sort(sequences.transpose(), axis=1)

        scored_by_energy: float = self.quality_function.energy(
            current_score_columns.tobytes(),
            num_sequences,
            seq_length
        )

        return scored_by_energy

    def anneal(self):
        current_sequences, current_temp, iteration, no_changes, total_accepted, total_rejected = self.sequences.copy(), self.temp, 0, 0, 0, 0

        current_score = self.cached_energy(current_sequences.tobytes(), current_sequences.shape[0], current_sequences.shape[1])
        initial_score, best_score = current_score, current_score

        with ThreadPoolExecutor(max_workers=1) as executor:  # Using one worker for logging
            while (current_temp := current_temp * self.cooling_rate) > self.min_temp and no_changes < self.no_changes_limit:
                previous_score = current_score

                neighbors_best_score = None
                for _neighbor_index in np.arange(self.iteration_neighbors):
                    current_neighbor_sequences = current_sequences.copy()

                    for _change_index in np.arange(np.random.randint(1, self.changes)):
                        num_sequences, seq_length = current_neighbor_sequences.shape

                        current_neighbor_sequences = self.sequence_generator.generate_new_sequences(
                            current_neighbor_sequences.tobytes(),
                            num_sequences=num_sequences,
                            seq_length=seq_length,
                            addition_deletion_prob=0.5
                        )

                    current_neighbor_score = self.cached_energy(
                        current_neighbor_sequences.tobytes(),
                        current_neighbor_sequences.shape[0],
                        current_neighbor_sequences.shape[1]
                    )

                    if (neighbors_best_score is None) or (current_neighbor_score  >= neighbors_best_score):
                        neighbors_best_score = current_neighbor_score
                        new_score = current_neighbor_score
                        new_sequences = current_neighbor_sequences

                # Once the best new state is found, we proceed with the regular SA flow
                score_change = new_score - current_score

                if (accepted := self.should_accept(score_change, current_temp)):
                    current_sequences = new_sequences
                    current_score = new_score
                    total_accepted += 1

                    if new_score > best_score:
                        best_score = new_score
                else:
                    total_rejected += 1

                if previous_score == current_score:
                    no_changes += 1
                else:
                    no_changes = 0

                executor.submit(
                    self.log_async,
                    iteration,
                    current_temp,
                    current_sequences.shape[1],
                    current_score,
                    previous_score,
                    new_score,
                    score_change,
                    initial_score,
                    best_score,
                    total_accepted,
                    total_rejected,
                    accepted,
                    no_changes
                )

                iteration += 1

        print(f"Annealing finished using heuristic {self.quality_function}. Initial score={initial_score}; Final score={current_score}; Final temp={current_temp}; Final iteration={iteration-1}. Consecutive no-changes events: {no_changes}")
        print("SequencesEnergyCache:", self.cached_energy.cache_info())
        print("EnergyForSequencesCache:", self.quality_function.energy.cache_info())
        print("EnergyForColumnsCache:", self.quality_function.get_column_energy.cache_info())

        return (current_sequences, current_temp, initial_score, current_score)

    def should_accept(self, delta, current_temp) -> bool:
        return (delta > 0) or (np.random.rand() < np.exp(delta / current_temp))
