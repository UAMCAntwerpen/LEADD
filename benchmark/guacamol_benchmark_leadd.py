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
    def __init__(self, fragmentation_settings_file, reconstruction_settings_file, output_directory, starting_population=None, benchmark=None, n_threads=1):
        self.output_directory = output_directory
        self.pool = joblib.Parallel(n_jobs=n_threads)
        self.starting_population = starting_population
        self.fragmentation_settings = pyLEADD.FragmentationSettings(fragmentation_settings_file)
        self.reconstruction_settings = pyLEADD.ReconstructionSettings(reconstruction_settings_file)
        if benchmark is not None:
            self.configure_goal_directed_benchmark(benchmark)
        self.leadd = pyLEADD.LEADD(reconstruction_settings=self.reconstruction_settings, output_directory_path=self.output_directory)
        self.n_generations = 0
        self.n_scored_molecules = 0
        self.time_spend_designing = 0
        self.time_spend_scoring = 0

    def configure_goal_directed_benchmark(self, benchmark):
        """
        This function relieves the user from setting reasonable settings for each
        benchmark separately.
        """
        # Configure the population settings. Some benchmarks request very few
        # solutions, but LEADD needs a sizeable population of molecules to work.
        # We set the minimum size of the population to 100, and generate/keep the
        # the same number of molecules each generation.
        requested_n_molecules = max(benchmark.contribution_specification.top_counts)
        n_molecules = max([requested_n_molecules, 100])
        self.reconstruction_settings.SetNSeeds(n_molecules)
        self.reconstruction_settings.SetNChildrenPerGeneration(n_molecules)
        self.reconstruction_settings.SetNSurvivorsPerGeneration(n_molecules)

        # If the present benchmark is a rediscovery benchmark only the best molecule
        # is considered in the score calculation. Hence, we use this maximum score
        # as termination criterion. Otherwise the average score is used.
        if requested_n_molecules == 1:
            self.reconstruction_settings.SetAverageScoreAsTerminationCriterion(False)
        else:
            self.reconstruction_settings.SetAverageScoreAsTerminationCriterion(True)

        # Configure the cyclicity settings. The mean number of cyclic atoms in
        # the molecules for each benchmark is hard-coded and, where possible,
        # derived from the structures towards which the scoring functions point.
        # The STDEV and Max of the distribution are defined relative to the mean.
        benchmark_mean_n_ring_atoms = {
            "logP (target: -1.0)": 6,
            "logP (target: 8.0)": 12,
            "TPSA (target: 150.0)": 18,
            "CNS MPO": 6,
            "QED": 18,
            "C7H8N2O2":6,
            "Pioglitazone MPO": 17,
            "Celecoxib rediscovery": 17,
            "Troglitazone rediscovery": 21,
            "Thiothixene rediscovery": 20,
            "Aripiprazole similarity": 22,
            "Albuterol similarity": 6,
            "Mestranol similarity": 17,
            "C11H24": 6,
            "C9H10N2O2PF2Cl": 6,
            "Median molecules 1": 6,
            "Median molecules 2": 23,
            "Osimertinib MPO": 21,
            "Fexofenadine MPO": 24,
            "Ranolazine MPO": 18,
            "Perindopril MPO": 9,
            "Amlodipine MPO": 12,
            "Sitagliptin MPO": 15,
            "Zaleplon MPO": 15,
            "Valsartan SMARTS": 17,
            "Deco Hop": 20,
            "Scaffold Hop": 20}
        mean_n_ring_atoms = benchmark_mean_n_ring_atoms[benchmark.name]
        self.reconstruction_settings.SetNRingAtomsMean(mean_n_ring_atoms)
        self.reconstruction_settings.SetNRingAtomsSTDEV(mean_n_ring_atoms / 2)
        self.reconstruction_settings.SetMaxNRingAtoms(mean_n_ring_atoms * 3)

        # If a starting population was specified, make sure that guided evolution
        # is disabled since both aren't compatible.
        if self.starting_population or benchmark.starting_population:
            if self.reconstruction_settings.UsingGuidedEvolution():
                self.reconstruction_settings.SetUsingGuidedEvolution(False)
                print("WARNING: Guided evolution isn't compatible with starting populations and was disabled.")

        # For reproducibility purposes set a random (but defined) PRNG seed.
        if self.reconstruction_settings.GetPRNGSeed() == 0:
            seed = random.randint(0, 2147483647)
            self.reconstruction_settings.SetPRNGSeed(seed)

    def score_children(self, scoring_function):
        smiles = [individual.GetSanitizedSMILES() for individual in self.leadd.GetPopulation() if individual.IsChild()]
        scores = self.pool(joblib.delayed(scoring_function.score)(s) for s in smiles)
        individual_idx = 0
        for individual in self.leadd.GetPopulation():
            if individual.IsChild():
                individual.SetScore(scores[individual_idx])
                individual_idx += 1
        assert(len(scores) == individual_idx)
        return len(scores)

    def score_population(self, scoring_function):
        smiles = [individual.GetSanitizedSMILES() for individual in self.leadd.GetPopulation()]
        scores = self.pool(joblib.delayed(scoring_function.score)(s) for s in smiles)
        individual_idx = 0
        for individual in self.leadd.GetPopulation():
            individual.SetScore(scores[individual_idx])
            individual_idx += 1
        assert(len(scores) == individual_idx)
        return len(scores)

    def generate_optimized_molecules(self, scoring_function: ScoringFunction, number_molecules: int, starting_population: Optional[List[str]] = None) -> List[str]:
        # Print the settings that will be used for the run. The population settings
        # were set before so we don't need to worry about them here.
        self.reconstruction_settings.Print()

        # If the benchmark specifies a starting population it has precedence over
        # the one specified by the user.
        if starting_population:
            self.starting_population = starting_population

        # If a starting population was specified, create it.
        if self.starting_population:
            population = pyLEADD.MakePopulationFromSMILES(population_smiles=self.starting_population, fragmentation_settings=self.fragmentation_settings)
            self.leadd.SetPopulation(population)

        # Score the starting population.
        if self.reconstruction_settings.ScoreFirstPopulation():
            self.score_population(scoring_function)

        # Design loop. Each iteration is a generation.
        while not self.leadd.TerminationCriteriaMet():
            # Expand the population with child molecules.
            start = time.time()
            self.leadd.GenerateChildren()
            self.time_spend_designing += time.time() - start
            # Score the children.
            start = time.time()
            n = self.score_children(scoring_function)
            self.leadd.IncreaseNScoringCalls(n)
            self.n_scored_molecules += n
            self.time_spend_scoring += time.time() - start
            # Keep the best individuals and wrap up the generation.
            start = time.time()
            self.leadd.SelectivePressure()
            self.time_spend_designing += time.time() - start
        self.n_generations += self.leadd.GetGenerationNumber()

        # Store the serialized population.
        self.leadd.SavePopulation(os.path.join(self.output_directory, "population.rst"))

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

def LoadStartingPopulations(starting_populations_json_path):
    with open(starting_populations_json_path, "r") as file:
        benchmark_smiles_scores = json.load(file)
        starting_populations = {}
        for benchmark_name, smiles_scores_pairs in benchmark_smiles_scores.items():
            starting_populations[benchmark_name] = [p[0] for p in smiles_scores_pairs]
        return starting_populations

def WriteGuacaMolBenchmarkResult(benchmark_result, output_json_path):
    results: Dict[str, Any] = OrderedDict()
    results["guacamol_version"] = guacamol.__version__
    results["timestamp"] = get_time_string()
    results["result"] = vars(benchmark_result)
    with open(output_json_path, "w") as file:
        file.write(json.dumps(results, indent=4))

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Evaluates LEADD's performance using one or more benchmarks from the GuacaMol goal-directed benchmark suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "fragmentation_settings_file", type=str,
        help="Path to the the text file containing the base molecule fragmentation settings. These are used to convert benchmark- and user-specified starting_populations into LEADD compatible molecular representations.")
    parser.add_argument(
        "reconstruction_settings_file", type=str,
        help="Path to the the text file containing the base molecule design settings. Note that this script may opt to change some of them for certain benchmarks.")
    parser.add_argument(
        "output_directory", type=str,
        help="Path to the output directory.")
    parser.add_argument(
        "-s", "--starting_populations", type=str,
        help="Path to a JSON file containing the starting populations of each benchmark. If not provided randomly generated molecules will be used as a starting population instead.")
    parser.add_argument(
        "-n", "--n_threads", type=int, default=1,
        help="Number of threads to use for molecule scoring.")
    parser.add_argument(
        "-v", "--benchmark_suite_version", type=str, default="v2",
        help="GuacaMol goal-directed benchmark suite to run.")
    parser.add_argument(
        "-b", "--benchmark_ids", type=int, nargs="+",
        help="IDs of the GuacaMol suite benchmarks to run (1-based indexing). Defaults to the entire benchmark suite.")
    args = parser.parse_args()
    assert(os.path.isfile(args.fragmentation_settings_file))
    assert(os.path.isfile(args.reconstruction_settings_file))
    assert(os.path.isdir(args.output_directory))
    if args.starting_populations:
        assert(os.path.isfile(args.starting_populations))
    if args.benchmark_ids:
        args.benchmark_ids = sorted(list(set(args.benchmark_ids)))
    else:
        if args.benchmark_suite_version == "v2":
            args.benchmark_ids = list(range(1,21))
        elif args.benchmark_suite_version == "trivial":
            args.benchmark_ids = list(range(1,8))
    return args

def Main():
    # Parse the command line arguments.
    args = ParseArgs()

    # Read the starting populations JSON. These are supposed to be the top N
    # molecules found by screening the GuacaMol "all" dataset using the same scoring
    # functions as the GuacaMol benchmark suite. In the original GuacaMol examples
    # this step is done as part of the benchmark, but since that implies repeating
    # the virtual screen for every replica I've separated it.
    benchmark_starting_populations = None
    if args.starting_populations:
        benchmark_starting_populations = LoadStartingPopulations(args.starting_populations)

    # Iterate over the specified benchmarks.
    benchmarks = goal_directed_benchmark_suite(version_name=args.benchmark_suite_version)
    for benchmark_id in args.benchmark_ids:
        benchmark = benchmarks[benchmark_id - 1]
        print(f"Current benchmark: {benchmark.name}")

        # Fetch the appropiate starting population for the benchmark.
        starting_population = None
        if args.starting_populations:
            starting_population = benchmark_starting_populations[benchmark.name]

        # Create an output directory for the benchmark.
        benchmark_output_directory = os.path.join(args.output_directory, str(benchmark_id))
        os.mkdir(benchmark_output_directory)

        # Initialize a LEADD instance.
        leadd = LEADD(fragmentation_settings_file=args.fragmentation_settings_file,
                      reconstruction_settings_file=args.reconstruction_settings_file,
                      output_directory=benchmark_output_directory,
                      starting_population=starting_population,
                      benchmark=benchmark,
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
