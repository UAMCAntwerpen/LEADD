#include "SAScore.hpp"

// Definition of the version integer of the FeatureLibrary class.
extern const unsigned feature_library_version = 20200604;


SAScoreHeuristic::SAScoreHeuristic() = default;
SAScoreHeuristic::SAScoreHeuristic(double mu, double sigma) :
	mu(mu), sigma(sigma) {};

double SAScoreHeuristic::operator()(double score, double sascore) const {
	double modifier = 0.0;
	if (sascore < mu) {
		modifier = 1.0;
	} else {
		double fraction = (sascore - mu) / sigma;
		modifier = exp(-0.5 * (fraction * fraction));
	};
	return score * modifier;
};


FeatureLibrary::FeatureLibrary() = default;

void FeatureLibrary::AddMolecule(const RDKit::ROMol& mol) {
	// Create a circular fingerprint (analogous to an ECFP_4 fingerprint) of the molecule.
	const RDKit::SparseIntVect<std::uint32_t>* fingerprint = RDKit::MorganFingerprints::getFingerprint(mol, 2);
	// Loop over the features in the fingerprint.
	FEATURE_COUNTS_MAP::iterator it, end_it = feature_counts.end();
	const std::map<std::uint32_t, int>& fingerprint_features = fingerprint->getNonzeroElements();
	for (const auto& fingerprint_feature : fingerprint_features) {
		std::uint32_t feature = fingerprint_feature.first;
		std::uint64_t count = fingerprint_feature.second;
		// Check if the feature has been previously recorded.
		it = feature_counts.find(feature);
		// If it has been recorded, increase its count.
		if (it != end_it) {
			it->second += count;
		// If not, insert it as a new feature in the counts table.
		} else {
			feature_counts.insert({ feature, count });
			end_it = feature_counts.end();
		};
	};
	// Delete the fingerprint to avoid memory leaks.
	delete fingerprint;
};

std::uint64_t FeatureLibrary::CalcPercentile() const {
	// Calculate the total count of all features.
	std::uint64_t total_count = 0;
	std::vector<std::uint32_t> features;
	features.reserve(feature_counts.size());
	for (const auto& feature_count : feature_counts) {
		features.push_back(feature_count.first);
		total_count += feature_count.second;
	};
	// Sort the features according to their counts in descending order.
	std::sort(features.begin(), features.end(),
		[&](std::uint32_t feature1, std::uint32_t feature2) {
			return feature_counts.at(feature1) > feature_counts.at(feature2);
		});
	// Find the 80th percentile count.
	double target_count = total_count * 0.8;
	std::uint64_t cumulative_count = 0, percentile = 0;
	for (std::uint32_t feature : features) {
		cumulative_count += feature_counts.at(feature);
		if (cumulative_count >= target_count) {
			percentile = feature_counts.at(feature);
			break;
		};
	};
	return percentile;
};

void FeatureLibrary::CalcFeatureScores() {
	// Calculates the scores associated with each feature within the library.
	feature_scores.clear();
	double percentile = static_cast<double>(CalcPercentile());
	for (const auto& feature_count : feature_counts) {
		std::uint32_t feature = feature_count.first;
		double count = static_cast<double>(feature_count.second);
		double score = log10(count / percentile);
		feature_scores.insert({ feature, score });
	};
};

double FeatureLibrary::GetScore(std::uint32_t feature) const {
	// Retrieves the score of a specific feature within the library.
	FEATURE_SCORES_MAP::const_iterator it = feature_scores.find(feature);
	if (it != feature_scores.end()) {
		return it->second;
	} else {
		return -4.0;
	};
};

double FeatureLibrary::CalcFeatureScore(const RDKit::ROMol& mol) const {
	// Calculates the feature score of a molecule.
	FEATURE_SCORES_MAP::const_iterator it, end_it = feature_scores.end();
	const RDKit::SparseIntVect<std::uint32_t>* fingerprint = RDKit::MorganFingerprints::getFingerprint(mol, 2);
	const std::map<std::uint32_t, int>& fingerprint_features = fingerprint->getNonzeroElements();
	double score = 0.0;
	std::uint64_t n_features = fingerprint_features.size();
	std::uint64_t feature_count = 0;
	for (const auto& fingerprint_feature : fingerprint_features) {
		std::uint32_t feature = fingerprint_feature.first;
		std::uint64_t count = fingerprint_feature.second;
		feature_count += count;
		it = feature_scores.find(feature);
		if (it != end_it) {
			score += it->second * count;
		} else {
			score += -4.0 * count;
		};
	};
	score /= static_cast<double>(feature_count);
	// Calculate a correction term for the fingerprint density that deems symmetrical
	// molecules easier to synthesize and viceversa.
	// NOTE: This term isn't part of the original publication, but was added by the
	// paper's author to the Python RDKit implementation.
	double symmetry_correction = 0.0;
	unsigned n_atoms = mol.getNumAtoms();
	if (n_atoms > n_features) {
		symmetry_correction = log(static_cast<double>(n_atoms) / static_cast<double>(n_features)) * 0.5;
	};
	delete fingerprint;
	return score + symmetry_correction;
};

double FeatureLibrary::CalcFeatureScore(const RDKit::SparseIntVect<std::uint32_t>& fingerprint, unsigned n_atoms) const {
	// Calculates the feature score of a molecule's fingerprint.
	FEATURE_SCORES_MAP::const_iterator it, end_it = feature_scores.end();
	const std::map<std::uint32_t, int>& fingerprint_features = fingerprint.getNonzeroElements();
	double score = 0.0;
	std::uint64_t n_unique_features = fingerprint_features.size();
	std::uint64_t feature_count = 0;
	for (const auto& fingerprint_feature : fingerprint_features) {
		std::uint32_t feature = fingerprint_feature.first;
		std::uint64_t count = fingerprint_feature.second;
		feature_count += count;
		it = feature_scores.find(feature);
		if (it != end_it) {
			score += it->second * count;
		} else {
			score += -4.0 * count;
		};
	};
	score /= static_cast<double>(feature_count);
	// Calculate a correction term for the fingerprint density that deems symmetrical
	// molecules easier to synthesize and viceversa.
	// NOTE: This term isn't part of the original publication, but was added by the
	// paper's author to the Python RDKit implementation.
	double symmetry_correction = 0.0;
	if (n_atoms > n_unique_features) {
		symmetry_correction = log(static_cast<double>(n_atoms) / static_cast<double>(n_unique_features)) * 0.5;
	};
	return score + symmetry_correction;
};

void FeatureLibrary::PrintCounts() const {
	for (const auto& feature_count : feature_counts) {
		std::cout << feature_count.first << ": " << feature_count.second << std::endl;
	};
};

void FeatureLibrary::PrintScores() const {
	for (const auto& feature_score : feature_scores) {
		std::cout << feature_score.first << ": " << feature_score.second << std::endl;
	};
};

unsigned CalcNumChiralCenters(const RDKit::ROMol& mol, bool detect_chirality) {
	unsigned n_chiral_centers = 0u;
	// If chirality wasn't specified previously:
	if (detect_chirality) {
		// Assign CIP ranks to the molcule's atoms.
		std::vector<unsigned> cip_ranks;
		std::vector<std::pair<int, int>> neighbors;
		RDKit::Chirality::assignAtomCIPRanks(mol, cip_ranks);
		// Loop over the atoms in the molecule.
		for (RDKit::ROMol::ConstAtomIterator ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
			// Use the CIP ranks to check if the atom is a chiral center.
			bool legal_center = false, has_dupes = true;
			std::tie(legal_center, has_dupes) = RDKit::Chirality::isAtomPotentialChiralCenter(*ai, mol, cip_ranks, neighbors);
			neighbors.clear();
			if (legal_center && !has_dupes) {
				++n_chiral_centers;
			};
		};
	// If chirality was already specified:
	} else {
		for (RDKit::ROMol::ConstAtomIterator ai = mol.beginAtoms(); ai != mol.endAtoms(); ++ai) {
			if ((*ai)->getChiralTag() != RDKit::Atom::ChiralType::CHI_UNSPECIFIED) {
				++n_chiral_centers;
			};
		};
	};
	return n_chiral_centers;
};

bool HasMacrocycle(const RDKit::ROMol& mol) {
	const RDKit::RingInfo* ring_info = mol.getRingInfo();
	const std::vector<std::vector<int>>& rings = ring_info->atomRings();
	for (const auto& ring : rings) {
		if (ring.size() > 8) {
			return true;
		};
	};
	return false;
};

double CalcComplexityScore(const RDKit::ROMol& mol, bool detect_chirality) {
	// Calculate the rest of the molecular complexity descriptors.
	double n_atoms = mol.getNumAtoms();
	double n_chiral_centers = CalcNumChiralCenters(mol, detect_chirality);
	double n_spiro_atoms = RDKit::Descriptors::calcNumSpiroAtoms(mol);
	double n_bridgehead_atoms = RDKit::Descriptors::calcNumBridgeheadAtoms(mol);
	// Use the molecular complexity descriptors to calculate the complexity penalties.
	double size_penalty = std::pow(n_atoms, 1.005) - n_atoms;
	double stereo_penalty = log10(n_chiral_centers + 1);
	double spiro_penalty = log10(n_spiro_atoms + 1);
	double bridgehead_penalty = log10(n_bridgehead_atoms + 1);
	double macrocycle_penalty = 0.0;
	if (HasMacrocycle(mol)) {
		macrocycle_penalty = log10(2);
	};
	return 0.0 - size_penalty - stereo_penalty - spiro_penalty - bridgehead_penalty - macrocycle_penalty;
};

double SAScore(const RDKit::ROMol& mol, const FeatureLibrary& feature_library, bool detect_chirality) {
	double feature_score = feature_library.CalcFeatureScore(mol);
	double complexity_score = CalcComplexityScore(mol, detect_chirality);
	double raw_sa_score = feature_score + complexity_score;
	// Scale the SAScore to a value between 1 and 10.
	double scaled_sa_score = 11.0 - (raw_sa_score + 5.0) * 9 / 6.5;
	if (scaled_sa_score > 8.0) {
		scaled_sa_score = 8.0 + log(scaled_sa_score - 8.0);
	};
	if (scaled_sa_score > 10.0) {
		scaled_sa_score = 10.0;
	} else if (scaled_sa_score < 1.0) {
		scaled_sa_score = 1.0;
	};
	return scaled_sa_score;
};

double SAScore(const RDKit::ROMol& mol, const RDKit::SparseIntVect<std::uint32_t>& fingerprint, const FeatureLibrary& feature_library, bool detect_chirality) {
	double feature_score = feature_library.CalcFeatureScore(fingerprint, mol.getNumAtoms());
	double complexity_score = CalcComplexityScore(mol, detect_chirality);
	double raw_sa_score = feature_score + complexity_score;
	// Scale the SAScore to a value between 1 and 10.
	double scaled_sa_score = 11.0 - (raw_sa_score + 5.0) * 9 / 6.5;
	if (scaled_sa_score > 8.0) {
		scaled_sa_score = 8.0 + log(scaled_sa_score - 8.0);
	};
	if (scaled_sa_score > 10.0) {
		scaled_sa_score = 10.0;
	} else if (scaled_sa_score < 1.0) {
		scaled_sa_score = 1.0;
	};
	return scaled_sa_score;
};
