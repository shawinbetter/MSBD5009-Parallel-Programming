// ==============================================================================
// 
//  Write your code in "gensupermer_pthread.cpp"
//  Keep this file unchanged
// 
// ==============================================================================

#pragma once

#include <chrono>
#include <string>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
int gensupermers(char *reads, int *reads_offs, int K, int P, int num_of_reads, std::vector<std::string> &all_supermers, int num_threads);