#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "LEADD.cpp"

PYBIND11_MAKE_OPAQUE(std::list<ReconstructedMol>)

PYBIND11_MODULE(pyLEADD, module) {
  module.doc() = "Lamarckian Evolutionary Algorithm for de novo Drug Design";

  pybind11::class_<FragmentationSettings>(module, "FragmentationSettings")
    .def(pybind11::init<const std::string&>(),
         pybind11::arg("settings_file_path"))
    .def("SetAtomTyping", &FragmentationSettings::SetAtomTyping)
    .def("SetSystematicFragmentation", &FragmentationSettings::SetSystematicFragmentation)
    .def("SetMinFragmentSize", &FragmentationSettings::SetMinFragmentSize)
    .def("SetMaxFragmentSize", &FragmentationSettings::SetMaxFragmentSize)
    .def("SetFragmentRings", &FragmentationSettings::SetFragmentRings)
    .def("SetMorganRadius", &FragmentationSettings::SetMorganRadius)
    .def("SetMorganConsiderChirality", &FragmentationSettings::SetMorganConsiderChirality)
    .def("SetHashedMorganNBits", &FragmentationSettings::SetHashedMorganNBits)
    .def("Print", &FragmentationSettings::Print);

  pybind11::class_<ReconstructionSettings>(module, "ReconstructionSettings")
    .def(pybind11::init<const std::string&, bool>(),
         pybind11::arg("settings_file_path"), pybind11::arg("check_paths") = true)
    .def("SetPRNGSeed", &ReconstructionSettings::SetPRNGSeed)
    .def("SetNRingAtomsMean", &ReconstructionSettings::SetNRingAtomsMean)
    .def("SetNRingAtomsSTDEV", &ReconstructionSettings::SetNRingAtomsSTDEV)
    .def("SetMaxNRingAtoms", &ReconstructionSettings::SetMaxNRingAtoms)
    .def("SetPeakNRingAtomsProbability", &ReconstructionSettings::SetPeakNRingAtomsProbability)
    .def("SetNSeeds", &ReconstructionSettings::SetNSeeds)
    .def("SetNGenerations", &ReconstructionSettings::SetNGenerations)
    .def("SetNChildrenPerGeneration", &ReconstructionSettings::SetNChildrenPerGeneration)
    .def("SetNSurvivorsPerGeneration", &ReconstructionSettings::SetNSurvivorsPerGeneration)
    .def("SetAverageScoreAsTerminationCriterion", &ReconstructionSettings::SetAverageScoreAsTerminationCriterion)
    .def("SetUsingGuidedEvolution", &ReconstructionSettings::SetUsingGuidedEvolution)
    .def("UsingGuidedEvolution", &ReconstructionSettings::UsingGuidedEvolution)
    .def("ScoreFirstPopulation", &ReconstructionSettings::ScoreFirstPopulation)
    .def("WriteReport", &ReconstructionSettings::WriteReport)
    .def("Print", &ReconstructionSettings::Print);

  pybind11::class_<ReconstructedMol>(module, "ReconstructedMol")
    .def("SetScore", &ReconstructedMol::SetScore)
    .def("GetScore", &ReconstructedMol::GetScore)
    .def("IsChild", &ReconstructedMol::IsChild)
    .def("GetSanitizedSMILES", &ReconstructedMol::GetSanitizedSMILES,
        pybind11::return_value_policy::reference_internal);

  pybind11::class_<std::list<ReconstructedMol>>(module, "LEADDPopulation")
    .def(pybind11::init<>())
    .def("__iter__", [](std::list<ReconstructedMol>& p) {
      return pybind11::make_iterator(p.begin(), p.end());
    }, pybind11::keep_alive<0, 1>());

  module.def("MakePopulationFromSMILES", &MakePopulationFromSMILES,
             pybind11::arg("population_smiles"), pybind11::arg("fragmentation_settings"),
             pybind11::return_value_policy::move);

  pybind11::class_<LEADD>(module, "LEADD")
    .def(pybind11::init<ReconstructionSettings, const std::string&>(),
         pybind11::arg("reconstruction_settings"), pybind11::arg("output_directory_path"))
    .def("TerminationCriteriaMet", &LEADD::TerminationCriteriaMet)
    .def("GenerateChildren", &LEADD::GenerateChildren)
    .def("SelectivePressure", &LEADD::SelectivePressure)
    .def("IncreaseNScoringCalls", &LEADD::IncreaseNScoringCalls)
    .def("UpdateReport", &LEADD::UpdateReport)
    .def("WriteOperationFrequenciesToReport", &LEADD::WriteOperationFrequenciesToReport)
    .def("SetPopulation", &LEADD::SetPopulation,
         pybind11::arg("population"), pybind11::arg("reset_weights") = true)
    .def("SavePopulation",
         pybind11::overload_cast<const std::string&>(&LEADD::SavePopulation, pybind11::const_),
         pybind11::arg("output_file_path"))
    .def("GetPopulation", &LEADD::GetPopulation,
         pybind11::return_value_policy::reference_internal)
    .def("GetBestIndividual", &LEADD::GetBestIndividual,
         pybind11::return_value_policy::reference_internal)
    .def("GetGenerationNumber", &LEADD::GetGenerationNumber)
    .def("GetBestScore", &LEADD::GetBestScore)
    .def("GetSettings", &LEADD::GetSettings)
    .def("Cleanup", &LEADD::Cleanup);
};
