
unsigned long long binomial(unsigned int n, unsigned int r) {
        unsigned long long nCk = 1;
        for(int i = 0; i < r; i++) {
                nCk *= (n - i);
                nCk /= (i + 1);
        }
        return nCk;
}

int next_pairing(int *pairing_conf, int numChunks) {

        if((pairing_conf[2] + 1) < numChunks) {
                pairing_conf[2]++;
                return 0;
        }
        if((pairing_conf[1] + 1) < numChunks) {
                pairing_conf[1]++;
                pairing_conf[2] = pairing_conf[1];
                return 0;
        }
        if((pairing_conf[0] + 1) < numChunks) {
                pairing_conf[0]++;
                pairing_conf[1] = pairing_conf[0];
                pairing_conf[2] = pairing_conf[0];
                return 0;
        }
	return 1;       // Problem
}

int next_permutation(int *combination_table, int numSNPs) {

                if((combination_table[2] + 1) < numSNPs) {
                        combination_table[2]++;
                        return 0;
                }
                if((combination_table[1] + 1) < numSNPs) {
                        combination_table[1]++;
                        combination_table[2] = 0;
                        return 0;
                }
                if((combination_table[0] + 1) < numSNPs) {
                        combination_table[0]++;
                        combination_table[1] = 0;
                        combination_table[2] = 0;
                        return 0;
                }
        return 1;       // Problem

}

int next_combination(int *combination_table, int numSNPs) {

                if((combination_table[2] + 1) < numSNPs) {
                        combination_table[2]++;
                        return 0;
                }
                if((combination_table[1] + 2) < numSNPs) {
                        combination_table[1]++;
                        combination_table[2] = combination_table[1] + 1;
                        return 0;
                }
                if((combination_table[0] + 3) < numSNPs) {
                        combination_table[0]++;
                        combination_table[1] = combination_table[0] + 1;
                        combination_table[2] = combination_table[0] + 2;
                        return 0;
                }
        return 1;       // Problem

}
