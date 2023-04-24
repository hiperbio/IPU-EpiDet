#ifndef COMBINATION_H_   
#define COMBINATION_H_

unsigned long long binomial(unsigned int n, unsigned int r);

int next_pairing(int *pairing_conf, int numChunks);

int next_permutation(int *combination_table, int numSNPs);

int next_combination(int *combination_table, int numSNPs);

#endif
