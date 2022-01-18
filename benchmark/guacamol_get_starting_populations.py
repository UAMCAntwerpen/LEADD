import argparse
import joblib
import time
import json
from guacamol.utils.chemistry import canonicalize
from guacamol.benchmark_suites import goal_directed_benchmark_suite

def ParseArgs():
    parser = argparse.ArgumentParser(
        description="Retrieves the top N molecules, according to a GuacaMol goal-directed benchmark suite, from a SMILES file and stores them in a JSON file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "input", type=str,
        help="Path to the input SMILES file.")
    parser.add_argument(
        "output", type=str,
        help="Path to the output JSON file.")
    parser.add_argument(
        "-n", "--n_molecules", type=int, default=100,
        help="Minimum number of top N molecules to write to the JSON file. If this number is lower than the number of molecules requested by the benchmark it's upgraded to the latter.")
    parser.add_argument(
        "-v", "--benchmark_suite_version", type=str, default="v2",
        help="GuacaMol goal-directed benchmark suite to use.")
    parser.add_argument(
        "-t", "--n_threads", type=int, default=1,
        help="Number of threads to use for molecule scoring.")
    args = parser.parse_args()
    return args

def LoadSMILES(thread_pool, smiles_file):
    with open(smiles_file, "r") as file:
        return thread_pool(joblib.delayed(canonicalize)(smiles.strip()) for smiles in file)

def TopN(thread_pool, smiles, scoring_function, n):
    joblist = (joblib.delayed(scoring_function.score)(s) for s in smiles)
    scores = thread_pool(joblist)
    scored_smiles = list(zip(smiles, scores))
    scored_smiles = sorted(scored_smiles, key=lambda x: x[1], reverse=True)
    return scored_smiles[:n]

def Main():
    # Parse the command line arguments.
    args = ParseArgs()

    # Create a pool of worker threads.
    thread_pool = joblib.Parallel(n_jobs=args.n_threads)

    # Load the SMILES of the input file into memory.
    smiles = LoadSMILES(thread_pool, args.input)
    print(f"{len(smiles)} SMILES loaded.")

    # Loop over the individual benchmarks in the benchmark suite.
    starting_populations = {}
    benchmarks = goal_directed_benchmark_suite(version_name=args.benchmark_suite_version)
    for benchmark in benchmarks:
        # Determine the number of molecules to write out.
        n = args.n_molecules
        n_requested = max(benchmark.contribution_specification.top_counts)
        if n_requested > n:
            n = n_requested
        # Use the benchmark's scoring function to score the loaded SMILES.
        start_time = time.time()
        print(f"Scoring molecules with '{benchmark.name}' scoring function.")
        top_n_smiles = TopN(thread_pool, smiles, benchmark.objective, n)
        time_spend = time.time() - start_time
        print(f"{time_spend}s spent retrieving {len(top_n_smiles)} molecules.")
        # Store the best molecules.
        starting_populations[benchmark.name] = top_n_smiles

    # Write the molecules to the output JSON file.
    with open(args.output, "w") as file:
        json.dump(starting_populations, file, indent=4)

if __name__ == "__main__":
    Main()
