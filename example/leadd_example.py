import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import pyLEADD # LEADD must be imported after the RDKit!

class ECFPSimilarityScorer:
    def __init__(self, reference_molecule, radius):
        reference_molecule_noh = Chem.RemoveHs(reference_molecule)
        self.reference = AllChem.GetMorganFingerprint(reference_molecule_noh, radius)
        self.radius = radius

    def __call__(self, molecule):
        molecule_noh = Chem.RemoveHs(molecule)
        fingerprint = AllChem.GetMorganFingerprint(molecule_noh, self.radius)
        similarity = DataStructs.TanimotoSimilarity(self.reference, fingerprint)
        return similarity

def ParseArgs():
    parser = argparse.ArgumentParser(description="Sample program using LEADD's Python bindings and a ECFP similarity scoring function.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("settings", type=str,
        help="Path to a LEADD settings file.")
    parser.add_argument("output", type=str,
        help="Path to a directory where LEADD's output is stored.")
    parser.add_argument("-r", "--reference", type=str, default="CN1CCN(CCNc2nc3cc(Cl)c(cc3nc2NCCN2CCN(C)CC2)C#CCCO)CC1",
        help="SMILES of the reference molecule used for scoring.")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="Flag to print the progress of the evolution.")
    return parser.parse_args()

def Main():
    # Parse the command line arguments.
    args = ParseArgs()

    # Create a ECFP4 similarity scorer object using the user-provided reference molecule.
    reference = Chem.MolFromSmiles(args.reference)
    scorer = ECFPSimilarityScorer(reference, 2)

    # Read the LEADD settings from the user-specified file.
    settings = pyLEADD.LEADDSettings(settings_file_path=args.settings)
    settings.Print()

    # Initialize a LEADD designer object. This also initializes the starting population.
    leadd = pyLEADD.LEADD(settings=settings, output_directory_path=args.output)

    # Assign preliminary scores to the population members.
    # NOTE: to be as generic as possible, LEADD uses SMILES to communicate with
    # the scoring function. Since RDKit molecules are the base for LEADD molecules
    # it should be possible to use them directly. However, the RDKit uses
    # Boost.Python for its Python bindings whereas LEADD uses pybind11, making
    # direct compatibility non-trivial.
    if settings.ScoreFirstPopulation():
        for molecule in leadd.GetPopulation():
            smiles = molecule.GetSanitizedSMILES()
            score = scorer(Chem.MolFromSmiles(smiles))
            molecule.SetScore(score)

    # Design loop. Each iteration is a generation.
    while not leadd.TerminationCriteriaMet():
        # Expand the population with child molecules.
        leadd.GenerateChildren()
        # Score the children.
        for molecule in leadd.GetPopulation():
            if molecule.IsChild():
                smiles = molecule.GetSanitizedSMILES()
                score = scorer(Chem.MolFromSmiles(smiles))
                molecule.SetScore(score)
                leadd.IncreaseNScoringCalls(1)
        # Keep the best individuals and wrap up the generation.
        leadd.SelectivePressure()
        # Print some statistics.
        if args.verbose:
            print(f"Generation: {leadd.GetGenerationNumber()}, Score: {leadd.GetBestScore():.4f}")

    best_individual = leadd.GetBestIndividual()
    print(f"Best individual: {best_individual.GetSanitizedSMILES()}, Score: {best_individual.GetScore():.4f}")

    # Write out the designed molecules.
    output_molecules_path = os.path.join(args.output, "designed_molecules.smi")
    with open(output_molecules_path, "w") as file:
        for molecule in leadd.GetPopulation():
            file.write(molecule.GetSanitizedSMILES(), molecule.GetScore())

    # Release LEADD's resources.
    leadd.Cleanup()

if __name__ == "__main__":
    Main()
