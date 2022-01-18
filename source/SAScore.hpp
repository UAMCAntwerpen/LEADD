#pragma once
#ifndef _SASCORE_HPP_
#define _SASCORE_HPP_

#include <iostream>
#include <unordered_map>
#include <math.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

// Declaration of the version integer of the FeatureLibrary class.
extern const unsigned feature_library_version;

// Definition of a RDKit::Chirality function. This function is implemented in RDKit's Chirality.cpp,
// but for some reason isn't made available to the user in Chirality.h.
namespace RDKit {
	namespace Chirality {
		typedef std::pair<int, int> INT_PAIR;
		typedef std::vector<INT_PAIR> INT_PAIR_VECT;
		std::pair<bool, bool> isAtomPotentialChiralCenter(const Atom* atom, const ROMol& mol, const UINT_VECT& ranks, INT_PAIR_VECT& nbrs);
	};
};

typedef std::unordered_map<std::uint32_t, std::uint64_t> FEATURE_COUNTS_MAP;
typedef std::unordered_map<std::uint32_t, double> FEATURE_SCORES_MAP;

struct SAScoreComponents {
	double sascore = 0.0;
	double feature_score = 0.0;
	double complexity_score = 0.0;
	double size_penalty = 0.0;
	double stereo_penalty = 0.0;
	double spiro_penalty = 0.0;
	double bridgehead_penalty = 0.0;
	double macrocycle_penalty = 0.0;
};

class SAScoreHeuristic {
	double mu = 0.0;
	double sigma = 0.0;

	public:
		SAScoreHeuristic();
		SAScoreHeuristic(double mu, double sigma);

		double operator()(double score, double sascore) const;
};

class FeatureLibrary {
	FEATURE_COUNTS_MAP feature_counts;
	FEATURE_SCORES_MAP feature_scores;

	public:
		FeatureLibrary();
		void AddMolecule(const RDKit::ROMol& mol);

		void CalcFeatureScores();

		double GetScore(std::uint32_t feature) const;
		double CalcFeatureScore(const RDKit::ROMol& mol) const;
		double CalcFeatureScore(const RDKit::SparseIntVect<std::uint32_t>& fingerprint, unsigned n_atoms) const;

		void PrintCounts() const;
		void PrintScores() const;

	private:
		std::uint64_t CalcPercentile() const;

		template <class Archive>
		void serialize(Archive& ar, const unsigned version = feature_library_version) {
			ar & feature_counts & feature_scores;
		};

	friend class boost::serialization::access;
};

unsigned CalcNumChiralCenters(const RDKit::ROMol& mol, bool detect_chirality = true);
bool HasMacrocycle(const RDKit::ROMol& mol);
double CalcComplexityScore(const RDKit::ROMol& mol, bool detect_chirality = true);
double CalcComplexityScore(const RDKit::ROMol& mol, SAScoreComponents& sascore_components, bool detect_chirality = true);
double ScaleSAScore(double raw_sascore);
double SAScore(const RDKit::ROMol& mol, const FeatureLibrary& feature_library, bool detect_chirality = true);
double SAScore(const RDKit::ROMol& mol, const FeatureLibrary& feature_library, SAScoreComponents& sascore_components, bool detect_chirality = true);
double SAScore(const RDKit::ROMol& mol, const RDKit::SparseIntVect<std::uint32_t>& fingerprint, const FeatureLibrary& feature_library, bool detect_chirality = true);


#endif
