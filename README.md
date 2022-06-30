# LEADD
## About
LEADD (Lamarckian Evolutionary Algorithm for *de novo* Drug Design) is a tool for molecular design and optimization.

Molecules are represented as meta-graphs of molecular fragments. Fragments are extracted by fragmenting an input virtual library, and broken bonds are converted to labelled attachment points (a.k.a. connectors). Molecules are reassembled by genetic operators that combine fragments through said connectors. Knowledge-based connectivity rules, extracted from the same input library, are enforced throughout the process. A population of molecules is evolved stochastically through use of genetic operators, with the goal of optimizing a user-provided scoring function. A Lamarckian evolutionary mechanism adjusts the reproductive behaviour of the molecules based on the outcome of the previous generation.

A more detailed description of the algorithm can be found in the corresponding paper:
[Kerstjens, A., De Winter, H. LEADD: Lamarckian evolutionary algorithm for de novo drug design. J Cheminform 14, 3 (2022). https://doi.org/10.1186/s13321-022-00582-y](https://doi.org/10.1186/s13321-022-00582-y).

## Authors
* Alan Kerstjens (https://github.com/alankerstjens)

## Dependencies
* [CMake](https://cmake.org/) (>= 3.15)
* [C++17](https://isocpp.org/) and [OpenMP 4.0](https://www.openmp.org/) compliant compiler
* [Boost](https://www.boost.org/) (>= 1.73). Make sure [Boost.Program_options](https://www.boost.org/doc/libs/1_73_0/doc/html/program_options.html) and [Boost.Serialization](https://www.boost.org/doc/libs/1_73_0/libs/serialization/doc/index.html) are part of your Boost distribution!
* [RDKit](https://rdkit.org/) (>= 2021.03.2)
* [SQLite](https://sqlite.org/index.html) (>= 3.31.1)
* [HDF5](https://www.hdfgroup.org/) (>= 1.10.4)

It's recommended to also build the Python bindings. On top of the above, you will need the following:
* [Python](https://www.python.org/) (>= 3.8.8)
* [pybind11](https://github.com/pybind/pybind11) (>= 2.6.2). If you clone this repository it should be included in the distribution.

If you are interested in re-running the benchmarks you will also need:
* [GuacaMol](https://github.com/BenevolentAI/guacamol)

The versions listed above have been tested to work on both Ubuntu 20.04 and Windows 10. Older versions may also work, with the exception of the RDKit due to [a bug in older versions](https://github.com/rdkit/rdkit/issues/3994).

## Installation
The following sections assume you are installing LEADD on a GNU/Linux machine. If this isn't the case you may have to adapt these instructions slightly.

Clone this repository (along with the pybind11 submodule). The resulting directory will be referred to as `${LEADD}`.

```bash
git clone --recurse-submodules https://github.com/UAMCAntwerpen/LEADD.git
```

Install the rest of the dependencies. If you don't mind a bloated installation you can install them into an Anaconda environment using the provided `leadd_conda_env.yml` file.

```bash
cd ${LEADD}
conda env create -f leadd_conda_env.yml
conda activate LEADD
```

If you installed the RDKit through Anaconda you can build LEADD as follows:

```bash
mkdir build && cd build
cmake -DRDKit_INCLUDE_DIRS=${CONDA_PREFIX}/include/rdkit -DRDKit_LIBRARY_DIRS=${CONDA_PREFIX}/lib ..
make install clean
```

Alternatively, if you installed the RDKit from source you should have an environment variable `${RDBASE}` pointing to the RDKit's root directory:

```bash
mkdir build && cd build
cmake -DRDKit_INCLUDE_DIRS=${RDBASE}/Code -DRDKit_LIBRARY_DIRS=${RDBASE}/lib ..
make install clean
```

If you don't want to build the Python bindings add the `-DBUILD_PYTHON_BINDINGS=OFF` flag to the CMake commands. Otherwise you probably want to add `${LEADD}/lib` to your `${PYTHONPATH}` environment variable, for instance in your `.bashrc` file.

```bash
export PYTHONPATH="${PYTHONPATH}:${LEADD}/lib"
```

## Example workflow
If you'd like to follow along you can find most referenced files in the [example directory](example). Most executables are parallelized with OpenMP and will by default use all cores of your CPU. You can change this behaviour by setting the `OMP_NUM_THREADS` environment variable before running the command.

The first step is to create the fragments' SQLite3 database. In this example we will fragment [`PGK1_ligands.smi`](example/PGK1_ligands.smi), splitting acyclic regions into atomic fragments and defining connectors with MMFF94 atom types.

```bash
cd ${LEADD}/example
${LEADD}/bin/Fragment -i PGK1_ligands.smi -o fragments.db -s fragmentation_settings.txt
```

Thereafter, we precompute which fragments are compatible with each connector and store the output in a file.

```bash
${LEADD}/bin/PrecalculateConnectionQueryResults -i fragments.db -o fragments.cqr -t reconstruction_settings.txt
```

If you want to enable guided evolution you will need to generate a fragments similarity matrix. If your fragments database is very large, think twice before creating this matrix as it [can be very expensive!](#What-should-I-consider-when-selecting-input-molecules-to-create-the-fragments-database)

```bash
${LEADD}/bin/MakeFragmentSimilarityMatrix -i fragments.db -o fragments.h5
```

If you would like to use the [SAScore](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3225829/) filter or score heuristic you will also need a feature count library. For this you can use any collection of drug-like molecules, for instance [CHEMBL](https://www.ebi.ac.uk/chembl/).

```bash
${LEADD}/bin/MakeFeatureLibrary -i CHEMBL.smi -o CHEMBL.fl
```

Now you are ready to run LEADD itself. Make sure the absolute paths to the files we just created are within your `reconstruction_settings.txt` file. You'll have to [adapt](#Settings) the [provided settings file](example/reconstruction_settings.txt) a bit.

If you are using LEADD programatically with the Python bindings make sure your Python interpreter can find the LEADD shared library, for instance by adding `${LEADD}/lib` to your `${PYTHONPATH}`. The [provided example script](example/leadd_example.py) will attempt to rediscover [one of the ligands we just fragmented](example/CHEMBL3732863.smi):

```bash
python ${LEADD}/example/leadd_example.py reconstruction_settings.txt . -v
```

Alternatively, you can run LEADD as a standalone executable. Normally this requires [coupling a scoring function](#Scoring-settings). For this example we can use the included topological similarity scoring function to rediscover [the same ligand](example/CHEMBL3732863.smi):

```bash
${LEADD}/bin/StandaloneLEADD -s reconstruction_settings.txt -o . -v
```

Remember that LEADD is a stochastic algorithm. Try running it a couple of times to see what happens.

## Settings
This section covers how you can configure LEADD through its settings files. For an executable's options run it from the command line with the `--help` flag.

The settings files are plain text files where each line consists of a `KEY VALUE` pair separated by a single space. Most `KEY`s are optional. If an optional `KEY` is omitted or if the `VALUE` is blank LEADD will use a default `VALUE` instead (see source code). If the `KEY` is required by the executable LEADD will throw an error. Unless otherwise specified, a `KEY` should be assumed to be optional. It's recommended to specify all `KEY`s, including optional ones, to avoid unexpected behaviour. Empty lines and lines starting with `#` aren't parsed.

### Fragmentation settings
A template [fragmentation settings file](example/fragmentation_settings.txt) is provided with this distribution. For default setting values users should reference the [source code](source/FragmentationSettings.hpp).

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `ATOM_TYPING_SCHEME` | Enum | `DUMMY`, `ATOMIC_NUMBER`, `MMFF`, `MORGAN`, `HASHED_MORGAN` | Atom typing scheme used to define connections |
| `SYSTEMATIC_FRAGMENTATION_SCHEME` | Enum | `NONE`, `LINEAR`, `SUBGRAPH`, `CIRCULAR` | Fragmentation scheme for molecular regions treated systematically (normally acyclic regions) |
| `MIN_FRAGMENT_SIZE` | Integer | >= 0 | Minimum systematic fragment size (in number of heavy atoms). Only used when `SYSTEMATIC_FRAGMENTATION_SCHEME != NONE` |
| `MAX_FRAGMENT_SIZE` | Integer | >= 0 && >= `MIN_FRAGMENT_SIZE` | Maximum systematic fragment size (in number of heavy atoms). Only used when `SYSTEMATIC_FRAGMENTATION_SCHEME != NONE` |
| `FRAGMENT_RINGS` | Boolean | 0 or 1 | Flag to fragment cyclic regions systematically. *For testing purposes only* |

#### Morgan atom typing settings
| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `MORGAN_RADIUS` | Integer | >= 0 | Circular atomic environment radius used to define Morgan atom types. Only used when `ATOM_TYPING_SCHEME == MORGAN` or `ATOM_TYPING_SCHEME == HASHED_MORGAN` |
| `MORGAN_CONSIDER_CHIRALITY` | Boolean | 0 or 1 | Flag to use atom chirality as atomic invariant during Morgan atom type calculation. Only used when `ATOM_TYPING_SCHEME == MORGAN` or `ATOM_TYPING_SCHEME == HASHED_MORGAN` |
| `HASHED_MORGAN_N_BITS` | Integer | > 0 | Modulo used to collapse raw Morgan feature IDs into Morgan atom types. Only used when `ATOM_TYPING_SCHEME == HASHED_MORGAN` |


### Reconstruction settings
A template [reconstruction settings file](example/reconstruction_settings.txt) is provided with this distribution. For default setting values users should reference the [source code](source/ReconstructionSettings.hpp).

#### Input file settings
These settings are **required** by LEADD.

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `FRAGMENT_DATABASE_FILE` | String | Absolute paths | Path to the SQLite3 database generated by `Fragment`|
| `CONNECTION_QUERY_RESULTS_FILE` | String | Absolute paths | Path to the `.cqr` generated by `PrecalculateConnectionQueryResults`|

#### Fragment sampling settings

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `ACYCLIC_FREQUENCY_GAMMA` | Float | Any | Exponent applied to acyclic fragment frequencies when calculating their sampling weight |
| `ACYCLIC_LEVEL_GAMMA` | Float | Any | Exponent applied to acyclic fragment sizes when calculating their sampling weight |
| `RING_FREQUENCY_GAMMA` | Float | Any | Exponent applied to cyclic fragment frequencies when calculating their sampling weight |
| `RING_LEVEL_GAMMA` | Float | Any | Exponent applied to cyclic fragment sizes when calculating their sampling weight |
| `SCORE_GAMMA` | Float | Any | Exponent applied to molecules' scores when calculating fragment sampling weights in transfections and parent sampling weights |

#### Genetic operation probabilities
| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `PERIPHERAL_EXPANSION_PROBABILITY` | Integer | >= 0 | Sampling weight of the peripheral expansion operator |
| `PERIPHERAL_DELETION_PROBABILITY` | Integer | >= 0 | Sampling weight of the peripheral deletion operator |
| `PERIPHERAL_SUBSTITUTION_PROBABILITY` | Integer | >= 0 | Sampling weight of the peripheral substitution operator |
| `PERIPHERAL_TRANSFECTION_PROBABILITY` | Integer | >= 0 | Sampling weight of the peripheral transfection operator |
| `INTERNAL_EXPANSION_PROBABILITY` | Integer | >= 0 | Sampling weight of the internal expansion operator |
| `INTERNAL_DELETION_PROBABILITY` | Integer | >= 0 | Sampling weight of the internal deletion operator |
| `INTERNAL_SUBSTITUTION_PROBABILITY` | Integer | >= 0 | Sampling weight of the internal substitution operator |
| `INTERNAL_TRANSFECTION_PROBABILITY` | Integer | >= 0 | Sampling weight of the internal transfection operator |
| `TRANSLATION_PROBABILITY` | Integer | >= 0 | Sampling weight of the translation operator |
| `STEREO_FLIP_PROBABILITY` | Integer | >= 0 | Sampling weight of the stereo flip operator |

#### Designed molecule settings
| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `N_RING_ATOMS_MEAN` | Integer | > 0 | Mean number of ring atoms in designed molecules |
| `N_RING_ATOMS_STDEV` | Float | > 0 | Standard deviation in the number of ring atoms in designed molecules |
| `MAX_N_RING_ATOMS` | Integer | > 0 | Maximum allowed number of ring atoms |
| `PEAK_N_RING_ATOMS_PROBABILITY` | Float | 0 < v < 1 | Probability of keeping the number of ring atoms constant at `N_RING_ATOMS_MEAN` |
| `MIN_SEED_SIZE` | Integer | > 0 | Minimum number of heavy atoms in random starting molecules |
| `ASSIGN_UNSPECIFIED_STEREO` | Boolean | 0 or 1 | Flag to assign random stereochemistry to unspecified chiral centers and stereochemical double bonds |

#### Evolution settings
| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `PRNG_SEED` | Integer | >= 0 | Pseudo-random number generator seed. `0` is reserved for random seeds |
| `N_SEEDS` | Integer | > 0 | Number of random starting molecules |
| `N_GENERATIONS` | Integer | > 0 | Maximum number of generations |
| `N_CHILDREN_PER_GENERATION` | Integer | > 0 | Number of children bred each generation |
| `N_SURVIVORS_PER_GENERATION` | Integer | > 0 | Number of molecules advancing to the next generation |
| `MAX_CHILD_SIMILARITY` | Float | 0 < v < 1 | Maximum topological similarity between any two molecules of the population |
| `TERMINATION_SCORE` | Float | > 0 | Population score at which the process terminates |
| `AVERAGE_SCORE_AS_TERMINATION_CRITERION` | Boolean | 0 or 1 | Flag to use the average (instead of best) molecule score as `TERMINATION_SCORE` |
| `MAX_GENERATIONS_STUCK` | Integer | >= 0 | Maximum number of generations without improvement in the `TERMINATION_SCORE` |
| `MAX_ATTEMPTS_PER_GENERATION` | Integer | >= 0 | Maximum number of attempts at evolving a molecule before it's skipped |

#### Synthetic accessibility settings
The SAScore filter will be employed if the `FEATURE_LIBRARY_FILE` is specified and `MAX_SASCORE < 10`. The SAScore heuristic will be employed if `FEATURE_LIBRARY_FILE` is specified and `USE_SASCORE_HEURISTIC == 1`. Both can be used simultaneously.

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `FEATURE_LIBRARY_FILE` | String | Absolute paths | Path to the `.fl` generated by `MakeFeatureLibrary` |
| `MAX_SASCORE` | Float | 0 < v <= 10 | Maximum SAScore enforced by the SAScore filter |
| `SASCORE_HEURISTIC_MU` | Float | > 0 | Parameter of the SAScore heuristic Gaussian function |
| `SASCORE_HEURISTIC_SIGMA` | Float | > 0 | Parameter of the SAScore heuristic Gaussian function |
| `USE_SASCORE_HEURISTIC` | Boolean | 0 or 1 | Flag to enable the SAScore heuristic |

#### Guided evolution settings
Guided evolution will be enabled when `SIMILARITY_MATRIX_FILE` is specified and at least one of `ACYCLIC_POSITIVE_REINFORCEMENT`, `ACYCLIC_NEGATIVE_REINFORCEMENT`, `RING_POSITIVE_REINFORCEMENT` or `RING_NEGATIVE_REINFORCEMENT` are > 0.

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `SIMILARITY_MATRIX_FILE` | String | Absolute paths | Path to the HDF5 file generated by `MakeFragmentSimilarityMatrix` |
| `ACYCLIC_LEARNING_SIMILARITY_THRESHOLD` | Float | >= 0 | Minimum topological similarity between two acyclic fragments for one to trigger weight adjustments of the other |
| `RING_LEARNING_SIMILARITY_THRESHOLD` | Float | >= 0 | Minimum topological similarity between two cyclic fragments for one to trigger weight adjustments of the other |
| `ACYCLIC_POSITIVE_REINFORCEMENT` | Float | >= 0 | Reinforcement factor when increasing acyclic fragments sampling weights |
| `ACYCLIC_NEGATIVE_REINFORCEMENT` | Float | >= 0 | Reinforcement factor when decreasing acyclic fragments sampling weights |
| `RING_POSITIVE_REINFORCEMENT` | Float | >= 0 | Reinforcement factor when increasing cyclic fragments sampling weights |
| `RING_NEGATIVE_REINFORCEMENT` | Float | >= 0 | Reinforcement factor when decreasing cyclic fragments sampling weights |

#### Restart settings
By default LEADD generates a random starting population. Alternatively, users may specify their own starting populations. This can be done in two ways:
* Programatically with the `LEADD::SetPopulation` function. This can be done with both the [C++](source/LEADD.hpp) and [Python](source/pyLEADD.cpp) APIs.
* By specifying a `.rst` `RESTART_INPUT_FILE`.

Restart `.rst` files can be generated either by LEADD itself (if `RESTART_OUTPUT_FILE` is specified), by the `ConvertToReconstruction` executable or programatically with the `MakePopulationFromSMILES` function.

However, the user should be aware of the limitations of the latter two approaches. Since LEADD operates with meta-graphs, which are more information rich than regular molecular graphs, some compromises had to be made:
* Acylic regions will always be transformed into atomic fragments
* The generated connectors must also be present in the fragments database. Otherwise the molecule will be skipped.
* You won't be able to use guided evolution. This is because we can't guarantee that the generated fragments are present in the database nor calculate fragment similarities on the go in an efficient way.

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `RESTART_INPUT_FILE` | String | Absolute paths | Path to an input `.rst` file generated either by `ConvertToReconstruction` or LEADD |
| `RESET_WEIGHTS_ON_RESTART` | Boolean | 0 or 1 | Flag to reset connector weight arrays when restarting LEADD |
| `RESTART_OUTPUT_FILE` | String | Absolute paths | Path to an output `.rst` file written by LEADD |
| `N_GENERATIONS_PER_SAVE` | Integer | > 0 | Number of generations between writing a `.rst` file |

#### Scoring settings
LEADD comes out of the box with a topological similarity scoring function which will be employed if only `TEMPLATE_MOL_SMILES` is specified. This scoring function should only be used for testing purposes. For real drug discovery projects you should couple your own scoring function. The recommended way of doing this is programatically. A Python [example](example/leadd_example.py) of this can be found in this repository.

As a legacy way, you can also use the `StandaloneLEADD` executable. In this case, you must specify in the settings file `SCORING_FUNCTION_INPUT_FILE`, `SCORING_FUNCTION_CALL` and `SCORES_OUTPUT_FILE`, and make sure that your scoring function reads and writes files with the correct format. I strongly discourage this use because it can put some strain on the file system. More importantly, **the scoring function is evaluated as a shell command as is, posing a significant security risk! The user is responsible for making sure this isn't abused.**

| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `SCORING_FUNCTION_INPUT_FILE` | String | Absolute paths | Path to the space-tabulated SMILES file containing SMILES-ID pairs written by `StandaloneLEADD` and read by `SCORING_FUNCTION_CALL` |
| `SCORING_FUNCTION_CALL` | String | Any | Command ran by `StandaloneLEADD` to convert `SCORING_FUNCTION_INPUT_FILE` into `SCORES_OUTPUT_FILE` |
| `SCORES_OUTPUT_FILE` | String | Absolute paths | Path to the space-tabulated file containing ID-score pairs generated by `SCORING_FUNCTION_CALL` and read by `StandaloneLEADD` |
| `TEMPLATE_MOL_SMILES` | String | SMILES | SMILES of the reference molecule used by LEADD's topological similarity scoring function |
| `SCORE_FIRST_POPULATION` | Boolean | 0 or 1 | Flag to score the randomly generated/read population |

#### Reporting settings
LEADD can write out some statistics to CSVs and generate SVGs of the molecules over time. These settings are intended mostly for debugging purposes.
| Key | Type | Valid values | Description |
|:---:|:----:|:------------:|:-----------:|
| `WRITE_REPORT` | Boolean | 0 or 1 | Flag to write out evolution statistics to CSVs |
| `MONITOR_BEST_MOLECULE` | Boolean | 0 or 1 | Flag to write out the top molecule's connector weight arrays to CSV files and generate SVGs of its molecular graph |

## FAQ
### Does LEADD support all types of molecules?
LEADD was designed first and foremost with drug-like molecules in mind, and may or may not work well with:
* Very small or very large molecules
* Non-organic elements. I've only tested molecules with H, C, N, O, S, P and halogens.

### Can LEADD design ring systems?
LEADD's genetic operators work correctly with existing ring systems but won't create new ones. This is an arbitrary limitation due to the complexity of designing reasonable ring systems. If you are interested in ring design you could incorporate a suitable cyclization operator.

### Why do I need to specify the number of ring atoms in my designed molecules?
If you don't fragment ring systems (recommended) and you aren't using atomic fragments, LEADD's fragmentation scheme causes an overrepresentation of acyclic fragments. This could cause the design of very flexible non-drug-like molecules. As a solution LEADD samples acyclic and cyclic fragments separately. The [user settings](#Designed-molecule-settings) determine the likelihood of sampling each type of fragment given the current number of ring atoms. However, this is just a probabilistic model, and the number of ring atoms in your designed molecules will deviate from it should that lead to better fitness values. If you don't have a rough idea of how cyclic your designed molecules should be you can specify a very broad range. You will probably have to compensate for this by increasing the number of generations/population size.

### What should I consider when selecting input molecules to create the fragments database?
* Maximize the ring diversity. Since LEADD doesn't create new ring systems this is crucial.
* Avoid entirely cyclical molecules. By default these won't be fragmented and will be skipped. Note that if you decide to split them LEADD won't reassemble ring systems (see above).
* Be careful when using combinatorial chemistry libraries. You could face fragment frequency imbalances.
* If you intend on using guided evolution you should limit the number (*N*) of fragments. The time to calculate and space to store the fragments similarity matrix scale with *O(N^2)* complexity. LEADD's execution time and memory use will scale with *O(N)* complexity.

### Can I use a 3D scoring function?
The molecule objects used by LEADD are 2D. Naturally, you are free to generate 3D coordinates for them and score those instead. If you are going to use 3D scoring functions it's recommended to [enable the stereoflip operator](#Genetic-operation-probabilities). Please be aware that the conformations/orientations of a molecule can change significantly between genetic operations, which also may have implications during guided evolution.


## Troubleshooting
If you have any problems, questions or suggestions please open a GitHub Issue/Discussion.
