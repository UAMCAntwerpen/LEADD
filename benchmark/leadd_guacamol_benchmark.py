import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2" # Silence TensorFlow. We don't use it anyways.
import time
import json
import random
import pickle
import joblib
import argparse
from typing import List, Dict, Optional, Any
from collections import OrderedDict
import guacamol
from guacamol.benchmark_suites import goal_directed_benchmark_suite
from guacamol.goal_directed_generator import GoalDirectedGenerator
from guacamol.scoring_function import ScoringFunction
from guacamol.utils.data import get_time_string
import pyLEADD

class LEADD(GoalDirectedGenerator):
    def __init__(self, settings_file, output_directory, starting_population=None, benchmark_id=None, n_threads=1):
        self.pool = joblib.Parallel(n_jobs=n_threads)
        self.starting_population = starting_population
        assert(os.path.isfile(settings_file))
        assert(os.path.isdir(output_directory))
        self.settings = pyLEADD.LEADDSettings(settings_file)
        if benchmark_id is not None:
            self.configure_goal_directed_suite_v2_benchmark(benchmark_id)
        self.leadd = pyLEADD.LEADD(settings=self.settings, output_directory_path=output_directory)
        self.n_generations = 0
        self.n_scored_molecules = 0
        self.time_spend_designing = 0
        self.time_spend_scoring = 0

    def configure_goal_directed_suite_v2_benchmark(self, benchmark_id):
        """
        This function relieves the user from setting reasonable settings for each
        benchmark separately. Only the goal-directed suite v2 is supported.
        """
        assert(benchmark_id >= 1 and benchmark_id <= 20)
        benchmark = goal_directed_benchmark_suite(version_name="v2")[benchmark_id - 1]

        # Configure the population settings. Some benchmarks request very few
        # solutions, but LEADD needs a sizeable population of molecules to work.
        # We set the minimum size of the population to 100, and generate/keep the
        # the same number of molecules each generation.
        requested_n_molecules = max(benchmark.contribution_specification.top_counts)
        n_molecules = max([requested_n_molecules, 100])
        self.settings.SetNSeeds(n_molecules)
        self.settings.SetNChildrenPerGeneration(n_molecules)
        self.settings.SetNSurvivorsPerGeneration(n_molecules)

        # Configure the cyclicity settings. The mean number of cyclic atoms in
        # the molecules for each benchmark is hard-coded and, where possible,
        # derived from the structures towards which the scoring functions point.
        # The STDEV and Max of the distribution are defined relative to the mean.
        benchmark_mean_n_ring_atoms = { 1: 17,  2: 21,  3: 20,  4: 22,  5:  6,
                                        6: 17,  7:  6,  8:  6,  9:  6, 10: 23,
                                       11: 21, 12: 24, 13: 18, 14:  9, 15: 12,
                                       16: 15, 17: 15, 18: 17, 19: 20, 20: 20}
        mean_n_ring_atoms = benchmark_mean_n_ring_atoms[benchmark_id]
        self.settings.SetNRingAtomsMean(mean_n_ring_atoms)
        self.settings.SetNRingAtomsSTDEV(mean_n_ring_atoms/2)
        self.settings.SetMaxNRingAtoms(mean_n_ring_atoms*3)

        # If the present benchmark is a rediscovery benchmark only the best molecule
        # is considered in the score calculation. Hence, we use this maximum score
        # as termination criterion. Otherwise the average score is used.
        if benchmark_id <= 3:
            self.settings.SetAverageScoreAsTerminationCriterion(False)
        else:
            self.settings.SetAverageScoreAsTerminationCriterion(True)

        # If a starting population was specified, make sure that guided evolution
        # is disabled since both aren't compatible.
        if self.starting_population or benchmark.starting_population:
            if self.settings.UsingGuidedEvolution():
                self.settings.SetUsingGuidedEvolution(False)
                print("WARNING: Guided evolution isn't compatible with starting populations and was disabled.")

        # For reproducibility purposes set a random (but defined) PRNG seed.
        seed = random.randint(0, 2147483647)
        self.settings.SetPRNGSeed(seed)

    def score_smiles(self, smiles, scoring_function):
        scores = self.pool(joblib.delayed(scoring_function.score)(s) for s in smiles)
        return scores

    def score_population(self, scoring_function):
        smiles = [individual.GetSanitizedSMILES() for individual in self.leadd.GetPopulation() if individual.IsChild()]
        scores = self.pool(joblib.delayed(scoring_function.score)(s) for s in smiles)
        individual_idx = 0
        for individual in self.leadd.GetPopulation():
            if individual.IsChild():
                individual.SetScore(scores[individual_idx])
                individual_idx += 1
        assert(len(scores) == individual_idx)
        return len(scores)

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int, starting_population: Optional[List[str]] = None) -> List[str]:
        # Print the settings that will be used for the run. The population settings
        # were set before so we don't need to worry about them here.
        self.settings.Print()

        # If the benchmark specifies a starting population it has precedence over
        # the one specified by the user.
        if starting_population:
            self.starting_population = starting_population

        # If a starting population was specified, create and score it.
        if self.starting_population:
            population = pyLEADD.MakePopulationFromSMILES(self.starting_population)
            self.leadd.SetPopulation(population)
            self.score_population(scoring_function)

        # Design loop. Each iteration is a generation.
        while not self.leadd.TerminationCriteriaMet():
            # Expand the population with child molecules.
            start = time.time()
            self.leadd.GenerateChildren()
            self.time_spend_designing += time.time() - start
            # Score the children.
            start = time.time()
            n = self.score_population(scoring_function)
            leadd.IncreaseNScoringCalls(n)
            self.n_scored_molecules += n
            self.time_spend_scoring += time.time() - start
            # Keep the best individuals and wrap up the generation.
            start = time.time()
            self.leadd.SelectivePressure()
            self.time_spend_designing += time.time() - start
        self.n_generations += self.leadd.GetGenerationNumber()

        # Pass on the final designed molecules to the benchmark suite for analysis.
        final_population_smiles = [individual.GetSanitizedSMILES() for individual in self.leadd.GetPopulation()]
        return final_population_smiles[:number_molecules]

    def add_performance_data_to_guacamol_benchmark_result(self, benchmark_result):
        benchmark_result.metadata["n_generations"] = self.n_generations
        benchmark_result.metadata["n_scored_molecules"] = self.n_scored_molecules
        benchmark_result.metadata["time_spend_designing"] = self.time_spend_designing
        benchmark_result.metadata["time_spend_scoring"] = self.time_spend_scoring

    def clean_up(self):
        self.leadd.Cleanup()

def WriteGuacaMolBenchmarkResult(benchmark_result, output_json_path):
    results: Dict[str, Any] = OrderedDict()
    results["guacamol_version"] = guacamol.__version__
    results["benchmark_suite_version"] = "v2"
    results["timestamp"] = get_time_string()
    results["result"] = vars(benchmark_result)
    with open(output_json_path, "w") as file:
        file.write(json.dumps(results, indent=4))

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Evaluates LEADD's performance using one or more benchmarks from the GuacaMol goal-directed benchmark suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "starting_populations",
        help="Path to a pickled dictionary containing the starting populations of each benchmark.")
    parser.add_argument(
        "settings_file",
        help="Path to the the text file containing the base LEADD settings. Note that this script may opt to change some of them for certain benchmarks.")
    parser.add_argument(
        "output_directory",
        help="Path to the output directory.")
    parser.add_argument(
        "-r", "--random_start", action="store_true", default=False,
        help="Flag to use randomly generated molecules as a starting population instead of those in the pickled starting_populations.")
    parser.add_argument(
        "-n", "--n_threads", type=int, default=1,
        help="Number of threads to use for molecule scoring.")
    parser.add_argument(
        "-b", "--benchmark_ids", type=int, nargs="+", default=[id for id in range(1, 21)],
        help="IDs of the GuacaMol benchmarks to run.")
    args = parser.parse_args()
    assert(os.path.isfile(args.starting_populations))
    assert(os.path.isfile(args.settings_file))
    assert(os.path.isdir(args.output_directory))
    args.benchmark_ids = sorted(list(set(args.benchmark_ids)))
    return args

def Main():
    # Parse the command line arguments.
    args = ParseArgs()

    # Read the pickled starting populations.
    # These are supposed to be the top N molecules found by screening the GuacaMol
    # "all" dataset in SMILES format (using the same scoring functions as the
    # GuacaMol benchmark suite) stored in a dictionary, with the benchmark ID as
    # the key. In the original GuacaMol examples this step is done as part of the
    # benchmark, but since that implies repeating the virtual screen for every
    # replica we have separated it.
    with open(args.starting_populations, "rb") as file:
        starting_populations = pickle.load(file)

    # Iterate over the specified benchmarks.
    benchmarks = goal_directed_benchmark_suite(version_name="v2")
    for benchmark_id in args.benchmark_ids:
        benchmark = benchmarks[benchmark_id - 1]
        print(f"Current benchmark: {benchmark.name}")

        # Fetch the appropiate starting population for the benchmark.
        if args.random_start:
            starting_population = None
        else:
            starting_population = starting_populations[benchmark_id]

        # Create an output directory for the benchmark.
        benchmark_output_directory = os.path.join(args.output_directory, str(benchmark_id))
        os.mkdir(benchmark_output_directory)

        # Initialize a LEADD instance.
        leadd = LEADD(settings_file=args.settings_file,
                      output_directory=benchmark_output_directory,
                      starting_population=starting_population,
                      benchmark_id=benchmark_id,
                      n_threads=args.n_threads)

        # Use the benchmark to assess LEADD's performance.
        benchmark_result = benchmark.assess_model(leadd)
        # Add finer-grained computational performance data to the result's metadata.
        leadd.add_performance_data_to_guacamol_benchmark_result(benchmark_result)

        # Write the genetic operation frequencies to the report.
        leadd.leadd.WriteOperationFrequenciesToReport()

        # Release all LEADD resources.
        leadd.clean_up()

        # Write the results to a JSON file.
        output_json_path = os.path.join(benchmark_output_directory, "goal_directed_results.json")
        WriteGuacaMolBenchmarkResult(benchmark_result, output_json_path)

if __name__ == "__main__":
    Main()
