#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include "sample.hpp"

using namespace std;

ALIGNMENT_SCORE_TYPE SCORING_FUNCTION(FeatureAlignment &a, FeatureAlignment &b){
	if (a.charge != b.charge) return -1e9;
	MZ_TYPE mzDifference = (a.mz > b.mz) ? a.mz - b.mz : b.mz - a.mz;
	return pow(a.intensity * b.intensity, 0.5 * a.quality * b.quality) * (0.051 / (0.001 + mzDifference) - 1);
}

FEATURE_SCORE_TYPE FEATURE_SCORE_FUNCTION(Feature &a){
	return pow(log(a.intensity * a.quality), 3);
}

int main(int argc, char** argv) {
	string path = argv[1], option;
	ifstream fList("list.txt");
	ifstream fMap("map.txt");
	vector<string> paths;
	map<int, set<int> > dict;
	map<double, int> rank;
	set<int> matched;
	Sample sample(path);
	Sample sampleSummary(sample, 30);
	for(string s; getline(fList, s);) paths.push_back(s);
	for(int mz, num; fMap >> mz >> num;) dict[mz].insert(num);
	
	for (auto &f: sampleSummary.features){
		for (int i: dict[(int)(f.mz * 100)]){
			matched.insert(i);
		}
	}
	for (int i: matched){
		Sample other(paths[i] + ".summary");
		double score = sampleSummary.alignmentScore(other);
		rank[-score] = i;
	}
	int cnt = 0;
	cerr << "Top Matched samples in the database:\n";
	for (auto e: rank){
		if (cnt == 3) break;
		cnt++;
		Sample other(paths[e.second]);
		double score = sampleSummary.alignmentScore(other);
		cerr << paths[e.second] << "\t" << score << endl;
	}
	cerr << "\nDo you want to add current file to the database? (y/N) ";
	getline(cin, option);
	if (option.size() > 0 && (option[0] == 'Y' || option[0] == 'y')){
		fList.close();
		fMap.close();
		ofstream ofList("list.txt");
		ofstream ofMap("map.txt");
		for (string &s: paths) ofList << s << endl;
		ofList << path << endl;
		sampleSummary.write(path + ".summary");
		for (auto &mz: dict)
			for (int num: mz.second)
				ofMap << mz.first << "\t" << num << endl;
		for (auto &f: sampleSummary.features){
			ofMap << (int)(f.mz * 100) << "\t" << paths.size() << endl; 
		}
	}
	return 0;
}
