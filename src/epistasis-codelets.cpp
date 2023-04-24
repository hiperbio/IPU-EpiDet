#include <poplar/Vertex.hpp>

#include <float.h>
#include <math.h>

#define COMB_SIZE    	27

using namespace poplar;


/* MACROS used in the context of 3-way searches for deriving ...
 * ... genotype counts for 19 out of 27 possible genotypes. */

#define CALC_MACRO_X_Y(x, y, z, SNP_A_index, SNP_B_index, a, b, x1,y1,z1, x2,y2,z2); {\
	observedValues[(COMB_SIZE * 1) + x*9+y*3+z] = contTab2wayCases_1_2[(a*3+b)] - (observedValues[(COMB_SIZE * 1) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 1) + x2*9+y2*3+z2]);\
	observedValues[(COMB_SIZE * 0) + x*9+y*3+z] = contTab2wayControls_1_2[(a*3+b)] - (observedValues[(COMB_SIZE * 0) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 0) + x2*9+y2*3+z2]);\
}

#define CALC_MACRO_X_Z(x, y, z, SNP_A_index, SNP_B_index, a, b, x1,y1,z1, x2,y2,z2) {\
	observedValues[(COMB_SIZE * 1) + x*9+y*3+z] = contTab2wayCases_1_3[SNP_B_index * 9 + (a*3+b)] - (observedValues[(COMB_SIZE * 1) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 1) + x2*9+y2*3+z2]);\
	observedValues[(COMB_SIZE * 0) + x*9+y*3+z] = contTab2wayControls_1_3[SNP_B_index * 9 + (a*3+b)] - (observedValues[(COMB_SIZE * 0) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 0) + x2*9+y2*3+z2]);\
}

#define CALC_MACRO_Y_Z(x, y, z, SNP_A_index, SNP_B_index, a, b, x1,y1,z1, x2,y2,z2) {\
	observedValues[(COMB_SIZE * 1) + x*9+y*3+z] = contTab2wayCases_2_3[(SNP_A_index * numSNPs + SNP_B_index) * 4 + (a*2+b)] - (observedValues[(COMB_SIZE * 1) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 1) + x2*9+y2*3+z2]);\
	observedValues[(COMB_SIZE * 0) + x*9+y*3+z] = contTab2wayControls_2_3[(SNP_A_index * numSNPs + SNP_B_index) * 4 + (a*2+b)] - (observedValues[(COMB_SIZE * 0) + x1*9+y1*3+z1] + observedValues[(COMB_SIZE * 0) + x2*9+y2*3+z2]);\
}



class EpistasisVertex1 : public MultiVertex {
	public:

		Input<Vector<unsigned long long>> datasetCases_1;
		Input<Vector<unsigned long long>> datasetControls_1;
		Input<Vector<unsigned long long>> datasetCases_2;
		Input<Vector<unsigned long long>> datasetControls_2;
		Input<Vector<unsigned long long>> datasetCases_3;
		Input<Vector<unsigned long long>> datasetControls_3;
		Input<int> numSNPs;			
		Input<int> casesSize;
		Input<int> controlsSize;

		Output<Vector<int>> contTab1wayCases_1;
		Output<Vector<int>> contTab1wayControls_1;
		Output<Vector<int>> contTab1wayCases_2;
		Output<Vector<int>> contTab1wayControls_2;
		Output<Vector<int>> contTab1wayCases_3;
		Output<Vector<int>> contTab1wayControls_3;


		bool compute(unsigned workerId) {


			int casesSizePack64 = (int) ceil(((float) casesSize) / 64.0f);
			int controlsSizePack64 = (int) ceil(((float) controlsSize) / 64.0f);


			// for constructing mask for last iteration of loops iterating cases / controls
			int numBitsControls = controlsSize % 64;
			int numBitsCases = casesSize % 64;

			unsigned long long maskCases =  0 ^ ((1ull << numBitsCases) - 1);      
			unsigned long long maskControls = 0 ^ ((1ull << numBitsControls) - 1);

			int numPacks64 = (int) ceil(((float) SAMPLES_MAX) / 64.0f);

			for(int i=workerId; i < numSNPs; i+=6) {

				contTab1wayCases_1[i * 3 + 0] = 0;
				contTab1wayCases_1[i * 3 + 1] = 0;
				contTab1wayCases_1[i * 3 + 2] = 0;

				contTab1wayCases_2[i * 3 + 0] = 0;
				contTab1wayCases_2[i * 3 + 1] = 0;
				contTab1wayCases_2[i * 3 + 2] = 0;

				contTab1wayCases_3[i * 3 + 0] = 0;
				contTab1wayCases_3[i * 3 + 1] = 0;
				contTab1wayCases_3[i * 3 + 2] = 0;

				contTab1wayControls_1[i * 3 + 0] = 0;
				contTab1wayControls_1[i * 3 + 1] = 0;
				contTab1wayControls_1[i * 3 + 2] = 0;

				contTab1wayControls_2[i * 3 + 0] = 0;
				contTab1wayControls_2[i * 3 + 1] = 0;
				contTab1wayControls_2[i * 3 + 2] = 0;

				contTab1wayControls_3[i * 3 + 0] = 0;
				contTab1wayControls_3[i * 3 + 1] = 0;
				contTab1wayControls_3[i * 3 + 2] = 0;

				int cases_i, controls_i;

				// CASES

				for(cases_i = 0; cases_i < (casesSizePack64 - 1); cases_i++) {

					contTab1wayCases_1[i * 3 + 0] += __builtin_popcountll(datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_1[i * 3 + 1] += __builtin_popcountll(datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_1[i * 3 + 2] += __builtin_popcountll(~(datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));

					contTab1wayCases_2[i * 3 + 0] += __builtin_popcountll(datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_2[i * 3 + 1] += __builtin_popcountll(datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_2[i * 3 + 2] += __builtin_popcountll(~(datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));

					contTab1wayCases_3[i * 3 + 0] += __builtin_popcountll(datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_3[i * 3 + 1] += __builtin_popcountll(datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
					contTab1wayCases_3[i * 3 + 2] += __builtin_popcountll(~(datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));

				}

				contTab1wayCases_1[i * 3 + 0] += __builtin_popcountll(datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_1[i * 3 + 1] += __builtin_popcountll(datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_1[i * 3 + 2] += __builtin_popcountll(maskCases & ~(datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));

				contTab1wayCases_2[i * 3 + 0] += __builtin_popcountll(datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_2[i * 3 + 1] += __builtin_popcountll(datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_2[i * 3 + 2] += __builtin_popcountll(maskCases & ~(datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));

				contTab1wayCases_3[i * 3 + 0] += __builtin_popcountll(datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_3[i * 3 + 1] += __builtin_popcountll(datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]);
				contTab1wayCases_3[i * 3 + 2] += __builtin_popcountll(maskCases & ~(datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + i] | datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + i]));



				// CONTROLS

				for(controls_i = 0; controls_i < (controlsSizePack64 - 1); controls_i++) {

					contTab1wayControls_1[i * 3 + 0] += __builtin_popcountll(datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_1[i * 3 + 1] += __builtin_popcountll(datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_1[i * 3 + 2] += __builtin_popcountll(~(datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));

					contTab1wayControls_2[i * 3 + 0] += __builtin_popcountll(datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_2[i * 3 + 1] += __builtin_popcountll(datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_2[i * 3 + 2] += __builtin_popcountll(~(datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));

					contTab1wayControls_3[i * 3 + 0] += __builtin_popcountll(datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_3[i * 3 + 1] += __builtin_popcountll(datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
					contTab1wayControls_3[i * 3 + 2] += __builtin_popcountll(~(datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));

				}

				contTab1wayControls_1[i * 3 + 0] += __builtin_popcountll(datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_1[i * 3 + 1] += __builtin_popcountll(datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_1[i * 3 + 2] += __builtin_popcountll(maskControls & ~(datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));

				contTab1wayControls_2[i * 3 + 0] += __builtin_popcountll(datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_2[i * 3 + 1] += __builtin_popcountll(datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_2[i * 3 + 2] += __builtin_popcountll(maskControls & ~(datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));

				contTab1wayControls_3[i * 3 + 0] += __builtin_popcountll(datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_3[i * 3 + 1] += __builtin_popcountll(datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]);
				contTab1wayControls_3[i * 3 + 2] += __builtin_popcountll(maskControls & ~(datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + i] | datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + i]));


			}

			return true;
		}
};



class EpistasisVertex2 : public MultiVertex {
	public:

		Input<Vector<unsigned long long>> datasetCases_1;
		Input<Vector<unsigned long long>> datasetControls_1;
		Input<Vector<unsigned long long>> datasetCases_2;
		Input<Vector<unsigned long long>> datasetControls_2;
		Input<Vector<unsigned long long>> datasetCases_3;
		Input<Vector<unsigned long long>> datasetControls_3;
		Input<int> numSNPs;			
		Input<int> casesSize;
		Input<int> controlsSize;

		Output<Vector<int>> contTab2wayCases_2_3;
		Output<Vector<int>> contTab2wayControls_2_3;

		bool compute(unsigned workerId) {

			int casesSizePack64 = (int) ceil(((float) casesSize) / 64.0f);
			int controlsSizePack64 = (int) ceil(((float) controlsSize) / 64.0f);


			// for constructing mask for last iteration of loops iterating cases / controls
			int numBitsControls = controlsSize % 64;
			int numBitsCases = casesSize % 64;
			unsigned long long maskCases = (1ull << numBitsCases) - 1;      
			unsigned long long maskControls = (1ull << numBitsControls) - 1;

                        int numPacks64 = (int) ceil(((float) SAMPLES_MAX) / 64.0f);


			for(int i_j=workerId; i_j < numSNPs * numSNPs; i_j+=6) { 

				for(int idx=0;idx<4;idx++) {
					contTab2wayCases_2_3[i_j * 4 + idx] = 0;
					contTab2wayControls_2_3[i_j * 4 + idx] = 0;

				}

				int i = (int) (((float) i_j) / ((float) numSNPs));	
				int j = i_j % numSNPs;

				int cases_i, controls_i;

				unsigned long long casesArr[2 * 2];
				unsigned long long controlsArr[2 * 2];

				for(cases_i = 0; cases_i < casesSizePack64; cases_i++) {

					casesArr[0] = datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + i];
					casesArr[1] = datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + i];

					casesArr[2] = datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + j];
					casesArr[3] = datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + j];

					for(int a_i = 0; a_i < 2; a_i++) {
						for(int b_i = 0; b_i < 2; b_i++) {

							unsigned long long acc_2_3 = 0xFFFFFFFFFFFFFFFF;

							unsigned int comb_i = a_i * 2 + b_i;

							acc_2_3 = acc_2_3 & casesArr[0 * 2 + a_i];
							acc_2_3 = acc_2_3 & casesArr[1 * 2 + b_i];

							contTab2wayCases_2_3[(i * numSNPs + j) * 4 + comb_i] += __builtin_popcountll(acc_2_3);
						}
					}
				}


				for(controls_i = 0; controls_i < controlsSizePack64; controls_i++) {

					controlsArr[0] = datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + i];
					controlsArr[1] = datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + i];

					controlsArr[2] = datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + j];
					controlsArr[3] = datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + j];

					for(int a_i = 0; a_i < 2; a_i++) {
						for(int b_i = 0; b_i < 2; b_i++) {

							unsigned long long acc_2_3 = 0xFFFFFFFFFFFFFFFF;

							unsigned int comb_i = a_i * 2 + b_i;

							acc_2_3 = acc_2_3 & controlsArr[0 * 2 + a_i];
							acc_2_3 = acc_2_3 & controlsArr[1 * 2 + b_i];

							contTab2wayControls_2_3[(i * numSNPs + j) * 4 + comb_i] += __builtin_popcountll(acc_2_3);
						}
					}
				}

			}

			return true;
		}
};

class EpistasisVertex3 : public MultiVertex {
	public:

		Input<Vector<unsigned long long>> datasetCases_1;
		Input<Vector<unsigned long long>> datasetControls_1;
		Input<Vector<unsigned long long>> datasetCases_2;
		Input<Vector<unsigned long long>> datasetControls_2;
		Input<Vector<unsigned long long>> datasetCases_3;
		Input<Vector<unsigned long long>> datasetControls_3;
		Input<Vector<int>> contTab1wayCases_1;
		Input<Vector<int>> contTab1wayControls_1;
		Input<Vector<int>> contTab1wayCases_2;
		Input<Vector<int>> contTab1wayControls_2;
		Input<Vector<int>> contTab1wayCases_3;
		Input<Vector<int>> contTab1wayControls_3;
		Input<Vector<int>> contTab2wayCases_2_3;
		Input<Vector<int>> contTab2wayControls_2_3;
		Input<Vector<float>> lgammaPrecalc;	
		Input<int> numSNPs;			
		Input<int> casesSize;
		Input<int> controlsSize;

		InOut<int> numTilePhases;
		InOut<Vector<float>> output;				// score for best combination
		InOut<Vector<unsigned long long>> output_index;		// best combination SNP indexes	


		bool compute(unsigned workerId) {


			int casesSizePack64 = (int) ceil(((float) casesSize) / 64.0f);
			int controlsSizePack64 = (int) ceil(((float) controlsSize) / 64.0f);

			// for constructing mask for last iteration of loops iterating cases / controls
			int numBitsControls = controlsSize % 64;
			int numBitsCases = casesSize % 64;
			unsigned long long maskCases = (1ull << numBitsCases) - 1;      
			unsigned long long maskControls = (1ull << numBitsControls) - 1;

                        int numPacks64 = (int) ceil(((float) SAMPLES_MAX) / 64.0f);

			float bestScore = FLT_MAX;
			unsigned long long int bestCombination;		

			for(int snp_A=workerId; snp_A < numSNPs; snp_A+=6) {	// each worker processes a different SNP A

				// Calculate 2-way tables for 'numSNPs' SNP A x SNP C combinations

				int contTab2wayCases_1_3[BETA_PARAM * 9];
				int contTab2wayControls_1_3[BETA_PARAM * 9];

				for(int idx=0; idx < numSNPs; idx++) {

					contTab2wayCases_1_3[idx * 9 + 0 * 3 + 0] = 0;
					contTab2wayCases_1_3[idx * 9 + 0 * 3 + 1] = 0;
					contTab2wayCases_1_3[idx * 9 + 1 * 3 + 0] = 0;
					contTab2wayCases_1_3[idx * 9 + 1 * 3 + 1] = 0;

					contTab2wayControls_1_3[idx * 9 + 0 * 3 + 0] = 0;
					contTab2wayControls_1_3[idx * 9 + 0 * 3 + 1] = 0;
					contTab2wayControls_1_3[idx * 9 + 1 * 3 + 0] = 0;
					contTab2wayControls_1_3[idx * 9 + 1 * 3 + 1] = 0;
				}

				for(int snp_C=0; snp_C < numSNPs; snp_C++) {    

					int cases_i, controls_i;

					unsigned long long casesArr[2 * 2];
					unsigned long long controlsArr[2 * 2];

					for(cases_i = 0; cases_i < casesSizePack64; cases_i++) {

						casesArr[0] = datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + snp_A];
						casesArr[1] = datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + snp_A];

						casesArr[2] = datasetCases_3[0 * numSNPs * numPacks64 + cases_i * numSNPs + snp_C];
						casesArr[3] = datasetCases_3[1 * numSNPs * numPacks64 + cases_i * numSNPs + snp_C];

						for(int a_i = 0; a_i < 2; a_i++) {
							for(int b_i = 0; b_i < 2; b_i++) {

								unsigned long long acc_1_3 = 0xFFFFFFFFFFFFFFFF;

								unsigned int comb_i = a_i * 3 + b_i;

								acc_1_3 = acc_1_3 & casesArr[0 * 2 + a_i];
								acc_1_3 = acc_1_3 & casesArr[1 * 2 + b_i];

								contTab2wayCases_1_3[snp_C * 9 + comb_i] += __builtin_popcountll(acc_1_3);
							}
						}
					}

					// {0, 2} = {0, :} − ({0, 0} + {0, 1})
					contTab2wayCases_1_3[snp_C * 9 + 0 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 0] - (contTab2wayCases_1_3[snp_C * 9 + 0 * 3 + 0] + contTab2wayCases_1_3[snp_C * 9 + 0 * 3 + 1]);
					// {1, 2} = {1, :} − ({1, 0} + {1, 1})
					contTab2wayCases_1_3[snp_C * 9 + 1 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 1] - (contTab2wayCases_1_3[snp_C * 9 + 1 * 3 + 0] + contTab2wayCases_1_3[snp_C * 9 + 1 * 3 + 1]);
					// {2, 0} = {:, 0} − ({0, 0} + {1, 0})
					contTab2wayCases_1_3[snp_C * 9 + 2 * 3 + 0] = contTab1wayCases_3[snp_C * 3 + 0] - (contTab2wayCases_1_3[snp_C * 9 + 0 * 3 + 0] + contTab2wayCases_1_3[snp_C * 9 + 1 * 3 + 0]);
					// {2, 1} = {:, 1} − ({0, 1} + {1, 1}
					contTab2wayCases_1_3[snp_C * 9 + 2 * 3 + 1] = contTab1wayCases_3[snp_C * 3 + 1] - (contTab2wayCases_1_3[snp_C * 9 + 0 * 3 + 1] + contTab2wayCases_1_3[snp_C * 9 + 1 * 3 + 1]);
					// {2, 2} = {2, :} − ({2, 0} + {2, 1})
					contTab2wayCases_1_3[snp_C * 9 + 2 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 2] - (contTab2wayCases_1_3[snp_C * 9 + 2 * 3 + 0] + contTab2wayCases_1_3[snp_C * 9 + 2 * 3 + 1]);

					for(controls_i = 0; controls_i < controlsSizePack64; controls_i++) {

						controlsArr[0] = datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + snp_A];
						controlsArr[1] = datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + snp_A];

						controlsArr[2] = datasetControls_3[0 * numSNPs * numPacks64 + controls_i * numSNPs + snp_C];
						controlsArr[3] = datasetControls_3[1 * numSNPs * numPacks64 + controls_i * numSNPs + snp_C];

						for(int a_i = 0; a_i < 2; a_i++) {
							for(int b_i = 0; b_i < 2; b_i++) {

								unsigned long long acc_1_3 = 0xFFFFFFFFFFFFFFFF;

								unsigned int comb_i = a_i * 3 + b_i;

								acc_1_3 = acc_1_3 & controlsArr[0 * 2 + a_i];
								acc_1_3 = acc_1_3 & controlsArr[1 * 2 + b_i];

								contTab2wayControls_1_3[snp_C * 9 + comb_i] += __builtin_popcountll(acc_1_3);
							}
						}
					}

					// {0, 2} = {0, :} − ({0, 0} + {0, 1})
					contTab2wayControls_1_3[snp_C * 9 + 0 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 0] - (contTab2wayControls_1_3[snp_C * 9 + 0 * 3 + 0] + contTab2wayControls_1_3[snp_C * 9 + 0 * 3 + 1]);
					// {1, 2} = {1, :} − ({1, 0} + {1, 1})
					contTab2wayControls_1_3[snp_C * 9 + 1 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 1] - (contTab2wayControls_1_3[snp_C * 9 + 1 * 3 + 0] + contTab2wayControls_1_3[snp_C * 9 + 1 * 3 + 1]);
					// {2, 0} = {:, 0} − ({0, 0} + {1, 0})
					contTab2wayControls_1_3[snp_C * 9 + 2 * 3 + 0] = contTab1wayControls_3[snp_C * 3 + 0] - (contTab2wayControls_1_3[snp_C * 9 + 0 * 3 + 0] + contTab2wayControls_1_3[snp_C * 9 + 1 * 3 + 0]);
					// {2, 1} = {:, 1} − ({0, 1} + {1, 1})
					contTab2wayControls_1_3[snp_C * 9 + 2 * 3 + 1] = contTab1wayControls_3[snp_C * 3 + 1] - (contTab2wayControls_1_3[snp_C * 9 + 0 * 3 + 1] + contTab2wayControls_1_3[snp_C * 9 + 1 * 3 + 1]);
					// {2, 2} = {2, :} − ({2, 0} + {2, 1})
					contTab2wayControls_1_3[snp_C * 9 + 2 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 2] - (contTab2wayControls_1_3[snp_C * 9 + 2 * 3 + 0] + contTab2wayControls_1_3[snp_C * 9 + 2 * 3 + 1]);

				}


				for(int snp_B=0; snp_B < numSNPs; snp_B++) {

					// Calculates 2-way tables for single SNP A x SNP B combination
					// Also, precomputes AND between SNP A x SNP B for all samples and 4 genotype combinations (those combining only alleles of type 0 or type 1)

					int contTab2wayCases_1_2[9];
					int contTab2wayControls_1_2[9];

					contTab2wayCases_1_2[0 * 3 + 0] = 0;
					contTab2wayCases_1_2[0 * 3 + 1] = 0;
					contTab2wayCases_1_2[1 * 3 + 0] = 0;
					contTab2wayCases_1_2[1 * 3 + 1] = 0;
					contTab2wayControls_1_2[0 * 3 + 0] = 0;
					contTab2wayControls_1_2[0 * 3 + 1] = 0;
					contTab2wayControls_1_2[1 * 3 + 0] = 0;
					contTab2wayControls_1_2[1 * 3 + 1] = 0;

					// Stores up to 125 bitchunks of 64 samples of each genotype kind 
					unsigned long long precombined_AxB_cases[125 * 4];	// only 4 2-way genotypes are stored (those without allele 2)
					unsigned long long precombined_AxB_controls[125 * 4];

					int cases_i, controls_i;

					unsigned long long casesArr[2 * 3];
					unsigned long long controlsArr[2 * 3];

					for(cases_i = 0; cases_i < casesSizePack64; cases_i++) {

						casesArr[0] = datasetCases_1[0 * numSNPs * numPacks64 + cases_i * numSNPs + snp_A];
						casesArr[1] = datasetCases_1[1 * numSNPs * numPacks64 + cases_i * numSNPs + snp_A];

						casesArr[2] = datasetCases_2[0 * numSNPs * numPacks64 + cases_i * numSNPs + snp_B];
						casesArr[3] = datasetCases_2[1 * numSNPs * numPacks64 + cases_i * numSNPs + snp_B];

						for(int a_i = 0; a_i < 2; a_i++) {
							for(int b_i = 0; b_i < 2; b_i++) {

								unsigned long long acc_1_2 = 0xFFFFFFFFFFFFFFFF;

								unsigned int comb_i = a_i * 3 + b_i;

								acc_1_2 = acc_1_2 & casesArr[0 * 2 + a_i];
								acc_1_2 = acc_1_2 & casesArr[1 * 2 + b_i];

								precombined_AxB_cases[cases_i * 4 + (a_i * 2 + b_i)] = acc_1_2;

								contTab2wayCases_1_2[comb_i] += __builtin_popcountll(acc_1_2);
							}
						}
					}

					// {0, 2} = {0, :} − ({0, 0} + {0, 1})
					contTab2wayCases_1_2[0 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 0] - (contTab2wayCases_1_2[0 * 3 + 0] + contTab2wayCases_1_2[0 * 3 + 1]);

					// {1, 2} = {1, :} − ({1, 0} + {1, 1})
					contTab2wayCases_1_2[1 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 1] - (contTab2wayCases_1_2[1 * 3 + 0] + contTab2wayCases_1_2[1 * 3 + 1]);

					// {2, 0} = {:, 0} − ({0, 0} + {1, 0})
					contTab2wayCases_1_2[2 * 3 + 0] = contTab1wayCases_2[snp_B * 3 + 0] - (contTab2wayCases_1_2[0 * 3 + 0] + contTab2wayCases_1_2[1 * 3 + 0]);

					// {2, 1} = {:, 1} − ({0, 1} + {1, 1})
					contTab2wayCases_1_2[2 * 3 + 1] = contTab1wayCases_2[snp_B * 3 + 1] - (contTab2wayCases_1_2[0 * 3 + 1] + contTab2wayCases_1_2[1 * 3 + 1]);

					// {2, 2} = {2, :} − ({2, 0} + {2, 1})
					contTab2wayCases_1_2[2 * 3 + 2] = contTab1wayCases_1[snp_A * 3 + 2] - (contTab2wayCases_1_2[2 * 3 + 0] + contTab2wayCases_1_2[2 * 3 + 1]);

					
					for(controls_i = 0; controls_i < controlsSizePack64; controls_i++) {

						controlsArr[0] = datasetControls_1[0 * numSNPs * numPacks64 + controls_i * numSNPs + snp_A];
						controlsArr[1] = datasetControls_1[1 * numSNPs * numPacks64 + controls_i * numSNPs + snp_A];

						controlsArr[2] = datasetControls_2[0 * numSNPs * numPacks64 + controls_i * numSNPs + snp_B];
						controlsArr[3] = datasetControls_2[1 * numSNPs * numPacks64 + controls_i * numSNPs + snp_B];

						for(int a_i = 0; a_i < 2; a_i++) {
							for(int b_i = 0; b_i < 2; b_i++) {

								unsigned long long acc_1_2 = 0xFFFFFFFFFFFFFFFF;

								unsigned int comb_i = a_i * 3 + b_i;

								acc_1_2 = acc_1_2 & controlsArr[0 * 2 + a_i];
								acc_1_2 = acc_1_2 & controlsArr[1 * 2 + b_i];

								precombined_AxB_controls[controls_i * 4 + (a_i * 2 + b_i)] = acc_1_2;


								contTab2wayControls_1_2[comb_i] += __builtin_popcountll(acc_1_2);
							}
						}
					}


					// {0, 2} = {0, :} − ({0, 0} + {0, 1})
					contTab2wayControls_1_2[0 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 0] - (contTab2wayControls_1_2[0 * 3 + 0] + contTab2wayControls_1_2[0 * 3 + 1]);

					// {1, 2} = {1, :} − ({1, 0} + {1, 1})
					contTab2wayControls_1_2[1 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 1] - (contTab2wayControls_1_2[1 * 3 + 0] + contTab2wayControls_1_2[1 * 3 + 1]);

					// {2, 0} = {:, 0} − ({0, 0} + {1, 0})
					contTab2wayControls_1_2[2 * 3 + 0] = contTab1wayControls_2[snp_B * 3 + 0] - (contTab2wayControls_1_2[0 * 3 + 0] + contTab2wayControls_1_2[1 * 3 + 0]);

					// {2, 1} = {:, 1} − ({0, 1} + {1, 1})
					contTab2wayControls_1_2[2 * 3 + 1] = contTab1wayControls_2[snp_B * 3 + 1] - (contTab2wayControls_1_2[0 * 3 + 1] + contTab2wayControls_1_2[1 * 3 + 1]);

					// {2, 2} = {2, :} − ({2, 0} + {2, 1})
					contTab2wayControls_1_2[2 * 3 + 2] = contTab1wayControls_1[snp_A * 3 + 2] - (contTab2wayControls_1_2[2 * 3 + 0] + contTab2wayControls_1_2[2 * 3 + 1]);


					for(int snp_C=0; snp_C < numSNPs; snp_C++) {	

						int i;
						float score = 0;

						unsigned int observedValues[2 * COMB_SIZE];	
						unsigned int observedValues_8gen[2 * 8];

						const int pow_table[10] = {1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683};

						int combination_table[3];
						combination_table[0] = snp_A;
						combination_table[1] = snp_B;
						combination_table[2] = snp_C;

						for(i=0;i<COMB_SIZE;i++) {
							observedValues[i] = 0;
							observedValues[COMB_SIZE + i] = 0;
						}

						for(i=0;i<8;i++) {
							observedValues_8gen[i] = 0;
							observedValues_8gen[8 + i] = 0;
						}

						unsigned int comb_i;

						unsigned long long casesArr_3[2];
						unsigned long long controlsArr_3[2];

						const unsigned long long * datasetCases_3_ptr0 = &datasetCases_3[0 * numSNPs * numPacks64];
						const unsigned long long * datasetCases_3_ptr1 = &datasetCases_3[1 * numSNPs * numPacks64];

						for(cases_i = 0; cases_i < casesSizePack64; cases_i++) {
							casesArr_3[0] = datasetCases_3_ptr0[cases_i * numSNPs + combination_table[2]];
							casesArr_3[1] = datasetCases_3_ptr1[cases_i * numSNPs + combination_table[2]];

							comb_i = 8 * 1;

							for(int a_i_b_i = 0; a_i_b_i < 4; a_i_b_i++) {
								for(int c_i = 0; c_i < 2; c_i++) {
									observedValues_8gen[comb_i] += __builtin_popcountll(precombined_AxB_cases[cases_i * 4 + a_i_b_i] & casesArr_3[c_i]);
									comb_i++;
								}
							}
						}

						const unsigned long long * datasetControls_3_ptr0 = &datasetControls_3[0 * numSNPs * numPacks64];
						const unsigned long long * datasetControls_3_ptr1 = &datasetControls_3[1 * numSNPs * numPacks64];

						for(controls_i = 0; controls_i < controlsSizePack64; controls_i++) {
							controlsArr_3[0] = datasetControls_3_ptr0[controls_i * numSNPs + combination_table[2]];
							controlsArr_3[1] = datasetControls_3_ptr1[controls_i * numSNPs + combination_table[2]];

							comb_i = 8 * 0;

							for(int a_i_b_i = 0; a_i_b_i < 4; a_i_b_i++) {
								for(int c_i = 0; c_i < 2; c_i++) {


									observedValues_8gen[comb_i] += __builtin_popcountll(precombined_AxB_controls[controls_i * 4 + a_i_b_i] & controlsArr_3[c_i]);

									comb_i++;
								}
							}
						}


						// Copies from contiguous 8 genotypes (for controls and cases) to the correct contingency table positions

						observedValues[0] = observedValues_8gen[0];
						observedValues[COMB_SIZE + 0] = observedValues_8gen[8 + 0];

						observedValues[1] = observedValues_8gen[1];
						observedValues[COMB_SIZE + 1] = observedValues_8gen[8 + 1];

						observedValues[3] = observedValues_8gen[2];
						observedValues[COMB_SIZE + 3] = observedValues_8gen[8 + 2];

						observedValues[4] = observedValues_8gen[3];
						observedValues[COMB_SIZE + 4] = observedValues_8gen[8 + 3];

						observedValues[9] = observedValues_8gen[4];
						observedValues[COMB_SIZE + 9] = observedValues_8gen[8 + 4];

						observedValues[10] = observedValues_8gen[5];
						observedValues[COMB_SIZE + 10] = observedValues_8gen[8 + 5];

						observedValues[12] = observedValues_8gen[6];
						observedValues[COMB_SIZE + 12] = observedValues_8gen[8 + 6];

						observedValues[13] = observedValues_8gen[7];
						observedValues[COMB_SIZE + 13] = observedValues_8gen[8 + 7];



						/* The frequency counts for the following 19 genotypes are analytically derived with simple arithmetic operations. */
						// {0,0,2}, {0,1,2}, {0,2,0}, {0,2,1}, {0,2,2}, {1,0,2}, {1,1,2}, {1,2,0}, {1,2,1}, {1,2,2}, {2,0,0}, {2,0,1}, {2,0,2}, {2,1,0}, {2,1,1}, {2,1,2}, {2,2,0}, {2,2,1}, {2,2,2}


						// $\{0,0,2\}$ & $\{0,0,:\} - (\{0,0,0\} + \{0,0,1\})$
						CALC_MACRO_X_Y(0,0,2, combination_table[0], combination_table[1], 0,0, 0,0,0, 0,0,1);

						// $\{0,1,2\}$ & $\{0,1,:\} - (\{0,1,0\} + \{0,1,1\})$
						CALC_MACRO_X_Y(0,1,2, combination_table[0], combination_table[1], 0,1, 0,1,0, 0,1,1);

						// $\{0,2,0\}$ & $\{0,:,0\} - (\{0,0,0\} + \{0,1,0\})$
						CALC_MACRO_X_Z(0,2,0, combination_table[0], combination_table[2], 0,0, 0,0,0, 0,1,0);

						// $\{0,2,1\}$ & $\{0,:,1\} - (\{0,0,1\} + \{0,1,1\})$
						CALC_MACRO_X_Z(0,2,1, combination_table[0], combination_table[2], 0,1, 0,0,1, 0,1,1);

						// $\{0,2,2\}$ & $\{0,:,2\} - (\{0,0,2\} + \{0,1,2\})$
						CALC_MACRO_X_Z(0,2,2, combination_table[0], combination_table[2], 0,2, 0,0,2, 0,1,2);

						// $\{1,0,2\}$ & $\{1,0,:\} - (\{1,0,0\} + \{1,0,1\})$
						CALC_MACRO_X_Y(1,0,2, combination_table[0], combination_table[1], 1,0, 1,0,0, 1,0,1);

						// $\{1,1,2\}$ & $\{1,1,:\} - (\{1,1,0\} + \{1,1,1\})$
						CALC_MACRO_X_Y(1,1,2, combination_table[0], combination_table[1], 1,1, 1,1,0, 1,1,1);

						// $\{1,2,0\}$ & $\{1,:,0\} - (\{1,0,0\} + \{1,1,0\})$
						CALC_MACRO_X_Z(1,2,0, combination_table[0], combination_table[2], 1,0, 1,0,0, 1,1,0);

						// $\{1,2,1\}$ & $\{1,:,1\} - (\{1,0,1\} + \{1,1,1\})$
						CALC_MACRO_X_Z(1,2,1, combination_table[0], combination_table[2], 1,1, 1,0,1, 1,1,1);

						// $\{1,2,2\}$ & $\{1,:,2\} - (\{1,0,2\} + \{1,1,2\})$
						CALC_MACRO_X_Z(1,2,2, combination_table[0], combination_table[2], 1,2, 1,0,2, 1,1,2);

						// $\{2,0,0\}$ & $\{:,0,0\} - (\{0,0,0\} + \{1,0,0\})$
						CALC_MACRO_Y_Z(2,0,0, combination_table[1], combination_table[2], 0,0, 0,0,0, 1,0,0);

						// $\{2,0,1\}$ & $\{:,0,1\} - (\{0,0,1\} + \{1,0,1\})$
						CALC_MACRO_Y_Z(2,0,1, combination_table[1], combination_table[2], 0,1, 0,0,1, 1,0,1);

						// $\{2,0,2\}$ & $\{2,0,:\} - (\{2,0,0\} + \{2,0,1\})$
						CALC_MACRO_X_Y(2,0,2, combination_table[0], combination_table[1], 2,0, 2,0,0, 2,0,1);

						// $\{2,1,0\}$ & $\{:,1,0\} - (\{0,1,0\} + \{1,1,0\})$
						CALC_MACRO_Y_Z(2,1,0, combination_table[1], combination_table[2], 1,0, 0,1,0, 1,1,0);

						// $\{2,1,1\}$ & $\{:,1,1\} - (\{0,1,1\} + \{1,1,1\})$
						CALC_MACRO_Y_Z(2,1,1, combination_table[1], combination_table[2], 1,1, 0,1,1, 1,1,1);

						// $\{2,1,2\}$ & $\{2,1,:\} - (\{2,1,0\} + \{2,1,1\})$
						CALC_MACRO_X_Y(2,1,2, combination_table[0], combination_table[1], 2,1, 2,1,0, 2,1,1);

						// $\{2,2,0\}$ & $\{2,:,0\} - (\{2,0,0\} + \{2,1,0\})$
						CALC_MACRO_X_Z(2,2,0, combination_table[0], combination_table[2], 2,0, 2,0,0, 2,1,0);

						// $\{2,2,1\}$ & $\{2,:,1\} - (\{2,0,1\} + \{2,1,1\})$
						CALC_MACRO_X_Z(2,2,1, combination_table[0], combination_table[2], 2,1, 2,0,1, 2,1,1);

						// $\{2,2,2\}$ & $\{2,:,2\} - (\{2,0,2\} + \{2,1,2\})$
						CALC_MACRO_X_Z(2,2,2, combination_table[0], combination_table[2], 2,2, 2,0,2, 2,1,2);


						for (i=0; i< COMB_SIZE; i++) {
							unsigned int zerosCount = observedValues[COMB_SIZE * 0 + i];	
							unsigned int onesCount = observedValues[COMB_SIZE * 1 + i];

							score = score + lgammaPrecalc[zerosCount] + lgammaPrecalc[onesCount] - lgammaPrecalc[zerosCount + onesCount + 1];
						}
						score = fabs(score);

						if( (score <= bestScore) && (combination_table[0] != combination_table[1]) && (combination_table[0] != combination_table[2]) && (combination_table[1] != combination_table[2]) ) {
							bestScore = score;
							bestCombination =  ( ((unsigned long long int)combination_table[0]) << 0) | (((unsigned long long int)combination_table[1]) << 16) | (((unsigned long long int)combination_table[2] ) << 32);
							bestCombination = bestCombination | (((unsigned long long int)numTilePhases ) << 48);
						}

					}
				}
			}


			if((bestScore <= output[workerId])) {	
				output[workerId] = bestScore;
				output_index[workerId] = bestCombination;
			}


			if(workerId == 0) {
				numTilePhases++;
			}

			return true;
		}
};


