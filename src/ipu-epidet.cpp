/**
 *
 * ipu-epidet.cpp: Third-order exhaustive epistasis detection on Graphcore IPU-POD systems 
 *
 */

#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poputil/TileMapping.hpp>
#include <poplar/IPUModel.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <experimental/filesystem>
#include "combination.hpp"

#define MAX_LINE_SIZE 	134217728
#define MAX_NUM_LINES 	100000
#define NUM_TILES 	1472   	// number of tiles of each IPU 
#define SAMPLES_MAX 	4096	// maximum samples of each type 

using namespace poplar;
using namespace poplar::program;

// Sends data from the host to the IPU
struct HostToDeviceCallback : public LegacyStreamCallback {

	poplar::ArrayRef<unsigned long long> in;

	HostToDeviceCallback(poplar::ArrayRef<unsigned long long> in) : in{in} {}

	void fetch(void *p) override {
		std::copy(in.begin(), in.end(), static_cast<unsigned long long *>(p));
	}

};

// Sends data from the IPU to the host
struct DeviceToHostCallback : public LegacyStreamCallback {

	poplar::ArrayRef<float> out;

	DeviceToHostCallback(poplar::ArrayRef<float> out) : out{out} {}

	void fetch(void *p) override {
		std::memcpy((void*) out.begin(), p, (NUM_TILES * 6) * sizeof(float));
	}

};

// Sends data from the IPU to the host
struct DeviceToHostCallback_ull : public LegacyStreamCallback {

	poplar::ArrayRef<unsigned long long> out;

	DeviceToHostCallback_ull(poplar::ArrayRef<unsigned long long> out) : out{out} {}

	void fetch(void *p) override {
		std::memcpy((void*) out.begin(), p, (NUM_TILES * 6) * sizeof(unsigned long long));
	}

};

/* Loads a sample from dataset in 0's, 1's and 2's format */
int getValues(char* line, u_char* data, u_char* data_target)
{
	int num = 0;
	const char* tok;
	for (tok = strtok(line, ",");
			tok && *tok;
			tok = strtok(NULL, ",\n"))
	{
		if(data != NULL) {
			data[num] = atoi(tok);
		}
		num++;
	}
	if(data_target != NULL) {
		data_target[0] = data[num - 1];
	}
	return num;
}

/* Prepares compute set adding vertices and performing tile mapping */
ComputeSet buildEpistasisComputeSet1(Graph &graph, Tensor numTilePhases, Tensor datasetCases_1, Tensor datasetControls_1, Tensor datasetCases_2, Tensor datasetControls_2, Tensor datasetCases_3, Tensor datasetControls_3, Tensor contTab1wayCases_1, Tensor contTab1wayControls_1, Tensor contTab1wayCases_2, Tensor contTab1wayControls_2, Tensor contTab1wayCases_3, Tensor contTab1wayControls_3, Tensor contTab2wayCases_2_3, Tensor contTab2wayControls_2_3, Tensor lgammaPrecalc, Tensor numSNPs, Tensor casesSize, Tensor controlsSize, Tensor output, Tensor output_index) {
	ComputeSet epistasisCS1 = graph.addComputeSet("epistasisCS1");

	for (unsigned tile_i = 0; tile_i < NUM_TILES; tile_i++) {     

		auto v = graph.addVertex(epistasisCS1,              
				"EpistasisVertex1",
				{

				// Inputs
				{"numSNPs", numSNPs},   
				{"casesSize", casesSize},
				{"controlsSize", controlsSize},
				{"datasetCases_1", datasetCases_1[tile_i] },
				{"datasetControls_1", datasetControls_1[tile_i]},
				{"datasetCases_2", datasetCases_2[tile_i]},
				{"datasetControls_2", datasetControls_2[tile_i]},
				{"datasetCases_3", datasetCases_3[tile_i]},
				{"datasetControls_3", datasetControls_3[tile_i]},

				// Outputs
				{"contTab1wayCases_1", contTab1wayCases_1[tile_i]},
				{"contTab1wayControls_1", contTab1wayControls_1[tile_i]},
				{"contTab1wayCases_2", contTab1wayCases_2[tile_i]},
				{"contTab1wayControls_2", contTab1wayControls_2[tile_i]},
				{"contTab1wayCases_3", contTab1wayCases_3[tile_i]},
				{"contTab1wayControls_3", contTab1wayControls_3[tile_i]},

				});

		graph.setTileMapping(v, tile_i);   
		graph.setPerfEstimate(v, 0);
	}

	return epistasisCS1;
}

/* Prepares compute set adding vertices and performing tile mapping */
ComputeSet buildEpistasisComputeSet2(Graph &graph, Tensor numTilePhases, Tensor datasetCases_1, Tensor datasetControls_1, Tensor datasetCases_2, Tensor datasetControls_2, Tensor datasetCases_3, Tensor datasetControls_3, Tensor contTab1wayCases_1, Tensor contTab1wayControls_1, Tensor contTab1wayCases_2, Tensor contTab1wayControls_2, Tensor contTab1wayCases_3, Tensor contTab1wayControls_3, Tensor contTab2wayCases_2_3, Tensor contTab2wayControls_2_3, Tensor lgammaPrecalc, Tensor numSNPs, Tensor casesSize, Tensor controlsSize, Tensor output, Tensor output_index) {
	ComputeSet epistasisCS2 = graph.addComputeSet("epistasisCS2");

	for (unsigned tile_i = 0; tile_i < NUM_TILES; tile_i++) {     

		auto v = graph.addVertex(epistasisCS2,             
				"EpistasisVertex2",
				{

				// Inputs
				{"numSNPs", numSNPs}, 
				{"casesSize", casesSize},
				{"controlsSize", controlsSize},
				{"datasetCases_1", datasetCases_1[tile_i] },
				{"datasetControls_1", datasetControls_1[tile_i]},
				{"datasetCases_2", datasetCases_2[tile_i]},
				{"datasetControls_2", datasetControls_2[tile_i]},
				{"datasetCases_3", datasetCases_3[tile_i]},
				{"datasetControls_3", datasetControls_3[tile_i]},

				// Outputs
				{"contTab2wayCases_2_3", contTab2wayCases_2_3[tile_i]},         
				{"contTab2wayControls_2_3", contTab2wayControls_2_3[tile_i]},      

				});

		graph.setTileMapping(v, tile_i); 
                graph.setPerfEstimate(v, 0);
	}

	return epistasisCS2;
}

/* Prepares compute set adding vertices and performing tile mapping */
ComputeSet buildEpistasisComputeSet3(Graph &graph, Tensor numTilePhases, Tensor datasetCases_1, Tensor datasetControls_1, Tensor datasetCases_2, Tensor datasetControls_2, Tensor datasetCases_3, Tensor datasetControls_3, Tensor contTab1wayCases_1, Tensor contTab1wayControls_1, Tensor contTab1wayCases_2, Tensor contTab1wayControls_2, Tensor contTab1wayCases_3, Tensor contTab1wayControls_3, Tensor contTab2wayCases_2_3, Tensor contTab2wayControls_2_3, Tensor lgammaPrecalc, Tensor numSNPs, Tensor casesSize, Tensor controlsSize, Tensor output, Tensor output_index) {
	ComputeSet epistasisCS3 = graph.addComputeSet("epistasisCS3");

	for (unsigned tile_i = 0; tile_i < NUM_TILES; tile_i++) {     

		auto v = graph.addVertex(epistasisCS3,         
				"EpistasisVertex3",
				{

                                // Inputs
				{"numTilePhases", numTilePhases[tile_i]},
				{"numSNPs", numSNPs},   
				{"casesSize", casesSize},
				{"controlsSize", controlsSize},
				{"datasetCases_1", datasetCases_1[tile_i] },
				{"datasetControls_1", datasetControls_1[tile_i]},
				{"datasetCases_2", datasetCases_2[tile_i]},
				{"datasetControls_2", datasetControls_2[tile_i]},
				{"datasetCases_3", datasetCases_3[tile_i]},
				{"datasetControls_3", datasetControls_3[tile_i]},
				{"contTab1wayCases_1", contTab1wayCases_1[tile_i]},
				{"contTab1wayControls_1", contTab1wayControls_1[tile_i]},
				{"contTab1wayCases_2", contTab1wayCases_2[tile_i]},
				{"contTab1wayControls_2", contTab1wayControls_2[tile_i]},
				{"contTab1wayCases_3", contTab1wayCases_3[tile_i]},
				{"contTab1wayControls_3", contTab1wayControls_3[tile_i]},
				{"contTab2wayCases_2_3", contTab2wayCases_2_3[tile_i]},         
				{"contTab2wayControls_2_3", contTab2wayControls_2_3[tile_i]},      
				{"lgammaPrecalc", lgammaPrecalc[tile_i]},

				// Outputs
				{"output", output[tile_i]},                     
				{"output_index", output_index[tile_i]},         

				});

		graph.setTileMapping(v, tile_i);   
                graph.setPerfEstimate(v, 0);
	}

	return epistasisCS3;
}

int main(int argc, char **argv) {

	/* The only input to the compiled host binary is the dataset to process */
	if(argc < 2) {
		printf("USE: infile\n");
		return 1;
	}

	unsigned int numIPUs = NUM_IPUS;	 	
	unsigned int superchunkSize = ALPHA_PARAM;	
        int h_numSNPs_per_chunk = BETA_PARAM;   
        int numMiniChunks = superchunkSize / h_numSNPs_per_chunk;	

        int h_casesSize = (int) ceil(((float) SAMPLES_MAX) / 64.0f);
        int h_controlsSize = (int) ceil(((float) SAMPLES_MAX) / 64.0f);

        int lgammaPrecalcSize = SAMPLES_MAX * 2;	// times 2 to account for cases and controls

        poplar::Executable exe;


        struct timespec t_start, t_end;
	
	clock_gettime(CLOCK_MONOTONIC, &t_start);	// start time measurement

#if defined (USE_IPU_HW)

        /* Attach to IPU device if using real hardware */

	auto manager = DeviceManager::createDeviceManager();
	auto devices = manager.getDevices(poplar::TargetType::IPU, numIPUs);

	printf("Number of available IPU devices: %lu\n", devices.size());

	std::cout << "Trying to attach to IPU device ...\n";
	auto it = std::find_if(devices.begin(), devices.end(), [](Device &device) {
			return device.attach();
			});

	if (it == devices.end()) {
		std::cerr << "Error attaching to device\n";
		return -1;
	}

	auto device = std::move(*it);
	std::cout << "Attached to IPU device " << device.getId() << std::endl;

	std::experimental::filesystem::path exeName_cwd = std::experimental::filesystem::current_path() / ("ipu_program_pod" + std::to_string(numIPUs) + ".bin");
	const auto exeName = exeName_cwd.string();

	try {	// if IPU graph/program binary exists, load it from storage

		auto inf = std::ifstream(exeName);
		exe = poplar::Executable::deserialize(inf);

	} catch (const poplar::poplar_error& e) {
		std::cout << "IPU graph/program binary ('" << exeName << "') targeting IPU-POD" << std::to_string(numIPUs) << " has not been compiled yet. Will do it now...\n";

#else

	// Use the IPUModel emulator to run IPU code on CPU
	poplar::IPUModel ipuModel;
	ipuModel.numIPUs = numIPUs;     	// number of virtual IPUs to emulate on CPU
	ipuModel.tilesPerIPU = NUM_TILES;       // number of tiles per virtual IPU (e.g. 1472 for MK2/GC200 microarchitecture)
	auto device = ipuModel.createDevice();

#endif

		/* Construction of the IPU compute graph */
		
		auto target = device.getTarget();
		
		Graph graph(target, replication_factor(numIPUs));
		graph.addCodelets("src/epistasis-codelets.cpp", CodeletFileType :: Auto, "-DBETA_PARAM=36 -DSAMPLES_MAX=4096");        // Adds codelets to the compute graph  

		// Used for reseting the counters for the algorithm phases, between consequent executions of the program
		Tensor zeroConstant = graph.addConstant(INT, {}, 0);
		Tensor maxfloatConstant = graph.addConstant(FLOAT, {}, FLT_MAX);

		Tensor numTilePhases = graph.addVariable(INT, {NUM_TILES}, "numTilePhases");
		for(int i=0; i < NUM_TILES; i++) {
			graph.setInitialValue(numTilePhases[i], 0);
		}

		Tensor datasetCases_1 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_1");
		Tensor datasetControls_1 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_1");
		Tensor datasetCases_2 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_2");
		Tensor datasetControls_2 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_2");
		Tensor datasetCases_3 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_3");
		Tensor datasetControls_3 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_3");

		Tensor datasetCases_superChunk1 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_superChunk1");
		Tensor datasetControls_superChunk1 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_superChunk1");
		Tensor datasetCases_superChunk2 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_superChunk2");
		Tensor datasetControls_superChunk2 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_superChunk2");
		Tensor datasetCases_superChunk3 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_casesSize)}, "datasetCases_superChunk3");
		Tensor datasetControls_superChunk3 = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, (unsigned int) (2 * h_numSNPs_per_chunk * h_controlsSize)}, "datasetControls_superChunk3");

		Tensor contTab1wayCases_1 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayCases_1");
		Tensor contTab1wayControls_1 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayControls_1");

		Tensor contTab1wayCases_2 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayCases_2");
		Tensor contTab1wayControls_2 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayControls_2");

		Tensor contTab1wayCases_3 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayCases_3");
		Tensor contTab1wayControls_3 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (3 * h_numSNPs_per_chunk)}, "contTab2wayControls_3");

		// Only 4 genotypes are saved / used for SNPs B x SNPs C
		Tensor contTab2wayCases_2_3 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (4 * h_numSNPs_per_chunk * h_numSNPs_per_chunk)}, "contTab2wayCases_2_3");
		Tensor contTab2wayControls_2_3 = graph.addVariable(INT, {NUM_TILES, (unsigned int) (4 * h_numSNPs_per_chunk * h_numSNPs_per_chunk)}, "contTab2wayControls_2_3");

		Tensor lgammaPrecalc = graph.addVariable(FLOAT, {NUM_TILES, (unsigned int) (lgammaPrecalcSize + 1)}, "lgammaPrecalc");

		Tensor numSNPs = graph.addVariable(INT, {}, "numSNPs");
		Tensor casesSize = graph.addVariable(INT, {}, "casesSize");
		Tensor controlsSize = graph.addVariable(INT, {}, "controlsSize");
		Tensor output = graph.addVariable(FLOAT, {NUM_TILES, 6}, "output");
		Tensor output_index = graph.addVariable(UNSIGNED_LONGLONG, {NUM_TILES, 6}, "output_index");

		poputil::mapTensorLinearly(graph,numSNPs);
		poputil::mapTensorLinearly(graph,casesSize);
		poputil::mapTensorLinearly(graph,controlsSize);

		poputil::mapTensorLinearly(graph,zeroConstant);
		poputil::mapTensorLinearly(graph,maxfloatConstant);


		for(int i=0; i < NUM_TILES; i++) {

			graph.setTileMapping(numTilePhases[i], i);

			graph.setTileMapping(datasetCases_1[i], i);
			graph.setTileMapping(datasetControls_1[i], i);
			graph.setTileMapping(datasetCases_2[i], i);
			graph.setTileMapping(datasetControls_2[i], i);
			graph.setTileMapping(datasetCases_3[i], i);
			graph.setTileMapping(datasetControls_3[i], i);

			graph.setTileMapping(datasetCases_superChunk1[i], i);
			graph.setTileMapping(datasetControls_superChunk1[i], i);
			graph.setTileMapping(datasetCases_superChunk2[i], i);
			graph.setTileMapping(datasetControls_superChunk2[i], i);
			graph.setTileMapping(datasetCases_superChunk3[i], i);
			graph.setTileMapping(datasetControls_superChunk3[i], i);

			graph.setTileMapping(contTab1wayCases_1[i], i);
			graph.setTileMapping(contTab1wayControls_1[i], i);
			graph.setTileMapping(contTab1wayCases_2[i], i);
			graph.setTileMapping(contTab1wayControls_2[i], i);
			graph.setTileMapping(contTab1wayCases_3[i], i);
			graph.setTileMapping(contTab1wayControls_3[i], i);

			graph.setTileMapping(contTab2wayCases_2_3[i], i);
			graph.setTileMapping(contTab2wayControls_2_3[i], i);

			graph.setTileMapping(output[i], i);
			graph.setTileMapping(output_index[i], i);

			graph.setTileMapping(lgammaPrecalc[i], i);

		}


		int NumElemsToTransfer = 9 * h_numSNPs_per_chunk * h_numSNPs_per_chunk;
		int number2wayPermutationsOfChunks = numMiniChunks * numMiniChunks;


		DataStream inStream_datasetCases_superChunk1[numMiniChunks];
		DataStream inStream_datasetControls_superChunk1[numMiniChunks];
		DataStream inStream_datasetCases_superChunk2[numMiniChunks];
		DataStream inStream_datasetControls_superChunk2[numMiniChunks];
		DataStream inStream_datasetCases_superChunk3[numMiniChunks];
		DataStream inStream_datasetControls_superChunk3[numMiniChunks];

		for(int i=0; i < numMiniChunks; i++) {

			inStream_datasetCases_superChunk1[i] = graph.addHostToDeviceFIFO("datasetCases_superChunk1_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_casesSize, ReplicatedStreamMode::REPLICATE);
			inStream_datasetControls_superChunk1[i] = graph.addHostToDeviceFIFO("datasetControls_superChunk1_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_controlsSize, ReplicatedStreamMode::REPLICATE);

			inStream_datasetCases_superChunk2[i] = graph.addHostToDeviceFIFO("datasetCases_superChunk2_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_casesSize, ReplicatedStreamMode::REPLICATE);
			inStream_datasetControls_superChunk2[i] = graph.addHostToDeviceFIFO("datasetControls_superChunk2_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_controlsSize, ReplicatedStreamMode::REPLICATE);

			inStream_datasetCases_superChunk3[i] = graph.addHostToDeviceFIFO("datasetCases_superChunk3_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_casesSize, ReplicatedStreamMode::REPLICATE);
			inStream_datasetControls_superChunk3[i] = graph.addHostToDeviceFIFO("datasetControls_superChunk3_" + std::to_string(i), UNSIGNED_LONGLONG, 2 * h_numSNPs_per_chunk * h_controlsSize, ReplicatedStreamMode::REPLICATE);
		}



		auto inStream_lgammaPrecalc = graph.addHostToDeviceFIFO("lgammaPrecalc", FLOAT, lgammaPrecalcSize + 1, ReplicatedStreamMode::BROADCAST);
		auto inStream_numSNPs = graph.addHostToDeviceFIFO("numSNPs", INT, 1, ReplicatedStreamMode::BROADCAST);
		auto inStream_casesSize = graph.addHostToDeviceFIFO("casesSize", INT, 1, ReplicatedStreamMode::BROADCAST);
		auto inStream_controlsSize = graph.addHostToDeviceFIFO("controlsSize", INT, 1, ReplicatedStreamMode::BROADCAST);

		// Device to Host
		auto outStream_output = graph.addDeviceToHostFIFO("output", FLOAT, NUM_TILES * 6);
		auto outStream_output_index = graph.addDeviceToHostFIFO("output_index", UNSIGNED_LONGLONG, NUM_TILES * 6);





		/* Construction of IPU control program */

		ComputeSet cs1 = buildEpistasisComputeSet1(graph, numTilePhases, datasetCases_1, datasetControls_1, datasetCases_2, datasetControls_2, datasetCases_3, datasetControls_3, contTab1wayCases_1, contTab1wayControls_1, contTab1wayCases_2, contTab1wayControls_2, contTab1wayCases_3, contTab1wayControls_3, contTab2wayCases_2_3, contTab2wayControls_2_3, lgammaPrecalc, numSNPs, casesSize, controlsSize, output, output_index);

		ComputeSet cs2 = buildEpistasisComputeSet2(graph, numTilePhases, datasetCases_1, datasetControls_1, datasetCases_2, datasetControls_2, datasetCases_3, datasetControls_3, contTab1wayCases_1, contTab1wayControls_1, contTab1wayCases_2, contTab1wayControls_2, contTab1wayCases_3, contTab1wayControls_3, contTab2wayCases_2_3, contTab2wayControls_2_3, lgammaPrecalc, numSNPs, casesSize, controlsSize, output, output_index);

		ComputeSet cs3 = buildEpistasisComputeSet3(graph, numTilePhases, datasetCases_1, datasetControls_1, datasetCases_2, datasetControls_2, datasetCases_3, datasetControls_3, contTab1wayCases_1, contTab1wayControls_1, contTab1wayCases_2, contTab1wayControls_2, contTab1wayCases_3, contTab1wayControls_3, contTab2wayCases_2_3, contTab2wayControls_2_3, lgammaPrecalc, numSNPs, casesSize, controlsSize, output, output_index);


		auto progInstance = Sequence{};
		progInstance.add(Execute(cs1)); // 1-way counts
		progInstance.add(Execute(cs2)); // 2-way BxC counts
		progInstance.add(Execute(cs3)); // 2-way AxB and AxC counts + 3-way counts


		auto prog = Sequence({
				Copy(inStream_numSNPs, numSNPs),
				Copy(inStream_casesSize, casesSize),
				Copy(inStream_controlsSize, controlsSize),
				});

		prog.add(Copy(inStream_lgammaPrecalc, lgammaPrecalc[0]));
		for(int idx=1; idx < NUM_TILES; idx++) {
			prog.add(Copy(lgammaPrecalc[0], lgammaPrecalc[idx]));
		}

		for(int idx=0; idx < NUM_TILES; idx++) {
			prog.add(Copy(zeroConstant, numTilePhases[idx]));
			prog.add(Copy(maxfloatConstant, output[idx][0]));
			prog.add(Copy(maxfloatConstant, output[idx][1]));
			prog.add(Copy(maxfloatConstant, output[idx][2]));
			prog.add(Copy(maxfloatConstant, output[idx][3]));
			prog.add(Copy(maxfloatConstant, output[idx][4]));
			prog.add(Copy(maxfloatConstant, output[idx][5]));
		}

		for(int chunk_i=0; chunk_i < numMiniChunks; chunk_i++) {

			prog.add(Copy(inStream_datasetCases_superChunk1[chunk_i], datasetCases_superChunk1[chunk_i]));
			prog.add(Copy(inStream_datasetControls_superChunk1[chunk_i], datasetControls_superChunk1[chunk_i]));

			prog.add(Copy(inStream_datasetCases_superChunk2[chunk_i], datasetCases_superChunk2[chunk_i]));
			prog.add(Copy(inStream_datasetControls_superChunk2[chunk_i], datasetControls_superChunk2[chunk_i]));

			prog.add(Copy(inStream_datasetCases_superChunk3[chunk_i], datasetCases_superChunk3[chunk_i]));
			prog.add(Copy(inStream_datasetControls_superChunk3[chunk_i], datasetControls_superChunk3[chunk_i]));
		}

		unsigned int numPairingsOfChunks = binomial(numMiniChunks, 3) * (((((double)numMiniChunks)+1) * (((double)numMiniChunks)+2))/ ((((double)numMiniChunks) - 2) * (((double)numMiniChunks) - 1)));

		unsigned int numComputePhases = ceil(((double) numPairingsOfChunks) / ((double) NUM_TILES));    // depending on alpha and beta parameters, last compute phase might not make use of all tiles

		Sequence multiplePhasesSeq = Sequence({});


		int miniChunk_pairing_conf[3] = {0,0,0};

		for(unsigned computePhase_i = 0; computePhase_i < numComputePhases; computePhase_i++) {

			for(int tile_i=0; tile_i < NUM_TILES; tile_i++) {

				// Populate 3 dataset chunks inside each tile with data that is on the first alpha/beta tiles of each IPU
				multiplePhasesSeq.add(Copy(datasetCases_superChunk1[ miniChunk_pairing_conf[0] ], datasetCases_1[tile_i]));
				multiplePhasesSeq.add(Copy(datasetControls_superChunk1[ miniChunk_pairing_conf[0] ], datasetControls_1[tile_i]));
				multiplePhasesSeq.add(Copy(datasetCases_superChunk2[ miniChunk_pairing_conf[1] ], datasetCases_2[tile_i]));
				multiplePhasesSeq.add(Copy(datasetControls_superChunk2[ miniChunk_pairing_conf[1] ], datasetControls_2[tile_i]));
				multiplePhasesSeq.add(Copy(datasetCases_superChunk3[  miniChunk_pairing_conf[2] ], datasetCases_3[tile_i]));
				multiplePhasesSeq.add(Copy(datasetControls_superChunk3[ miniChunk_pairing_conf[2] ], datasetControls_3[tile_i]));

				int ret = next_pairing(miniChunk_pairing_conf, numMiniChunks);
			}

			multiplePhasesSeq.add(progInstance);
		}

		prog.add(multiplePhasesSeq);

		// Copies outputs to host
		prog.add(Copy(output, outStream_output));
		prog.add(Copy(output_index, outStream_output_index));


		/* Compilation of the IPU graph and programs */

		std::vector<poplar::program::Program> progs;
		progs.push_back(prog);
		exe = poplar::compileGraph(graph, progs);


#if defined (USE_IPU_HW)
		// Saving the IPU executable for later use, if the constructed graph is targeting real IPU devices
		auto outf = std::ofstream(exeName);
		exe.serialize(outf);
	}
#endif

	/* Load executable in Poplar engine */

	poplar::Engine engine(std::move(exe));
	engine.load(device);


	clock_gettime(CLOCK_MONOTONIC, &t_end);		// stop time measurement
	double timing_duration_compilingOrLoadingExe  = ((t_end.tv_sec + ((double) t_end.tv_nsec / 1000000000)) - (t_start.tv_sec + ((double) t_start.tv_nsec / 1000000000)));


	printf("Number of tiles to use per IPU processor: %d\n", NUM_TILES);

	int sampleSize;


	clock_gettime(CLOCK_MONOTONIC, &t_start);	// start time measurement

	/* Load dataset from storage */

	FILE* stream = fopen(argv[1], "r");	// file with dataset

	char * line = (char *) malloc(MAX_LINE_SIZE * sizeof(char));
	char * ret_fgets = fgets(line, MAX_LINE_SIZE, stream);
	int numCols = getValues(line, NULL, NULL);

	u_char * dataset = (u_char*) malloc(sizeof(u_char) * MAX_NUM_LINES * numCols);
	u_char * dataset_target = (u_char*) malloc(sizeof(u_char) * MAX_NUM_LINES);

	if(dataset == NULL) {
		printf("\nMemory allocation for dataset (genotype) failed.\n");
	}
	if(dataset_target == NULL) {
		printf("\nMemory allocation for dataset (phenotype) failed.\n");
	}

	sampleSize = 0;
	while (fgets(line, MAX_LINE_SIZE, stream))
	{
		getValues(line, dataset + (numCols * sampleSize), dataset_target + sampleSize);
		sampleSize++;
	}


	/* Counts the number of controls (0s) and cases (1s) */
	int numCases = 0;
	int numControls = 0;
	for(int i=0; i < sampleSize; i++) {
		if(dataset_target[i] == 1) {

			numCases++;
		}
		else {
			numControls++;
		}
	}

	printf("sample size: %d\n#cases: %d, #controls:%d\n", sampleSize, numCases, numControls);


	int h_numSNPs = numCols - 1;

	unsigned int numSuperChunks = ceil(((float) h_numSNPs) / superchunkSize);	


	auto hDatasetCases = std::vector<unsigned long long>(2 * h_numSNPs_per_chunk * numSuperChunks * numMiniChunks * h_casesSize);		// only 2 genotypes are represented in dataset internal representation
	auto hDatasetControls = std::vector<unsigned long long>(2 * h_numSNPs_per_chunk * numSuperChunks * numMiniChunks * h_controlsSize);


	int numSamplesOnes_64packed = h_casesSize;
	int numSamplesZeros_64packed = h_controlsSize;


	/* Binarizes dataset */
	for(int column_j=0; column_j < (numCols - 1); column_j++) {     
		int numSamples0Found = 0;
		int numSamples1Found = 0;
		for(int line_i=0; line_i < sampleSize; line_i++) {

			int datasetElement = dataset[line_i * numCols + column_j];


			if(dataset_target[line_i] == 1) {

				if(datasetElement == 2) {
					numSamples1Found++;
					continue;
				}

				int Ones_index = ((int) (column_j / h_numSNPs_per_chunk)) * h_numSNPs_per_chunk * 2 * numSamplesOnes_64packed + datasetElement * numSamplesOnes_64packed * h_numSNPs_per_chunk + ((int)(numSamples1Found / 64.0f)) * h_numSNPs_per_chunk + (column_j % h_numSNPs_per_chunk);
				hDatasetCases[Ones_index] = hDatasetCases[Ones_index] | (((unsigned long long) 1) << (numSamples1Found % 64));
				numSamples1Found++;
			}

			else {

				if(datasetElement == 2) {
					numSamples0Found++;
					continue;
				}

				int Zeros_index = ((int) (column_j / h_numSNPs_per_chunk)) * h_numSNPs_per_chunk * 2 * numSamplesZeros_64packed + datasetElement * numSamplesZeros_64packed * h_numSNPs_per_chunk + ((int)(numSamples0Found / 64.0f)) * h_numSNPs_per_chunk + (column_j % h_numSNPs_per_chunk);
				hDatasetControls[Zeros_index] = hDatasetControls[Zeros_index] | (((unsigned long long) 1) << (numSamples0Found % 64));
				numSamples0Found++;
			}

		}
	}


	clock_gettime(CLOCK_MONOTONIC, &t_end);		// stop time measurement
	double timing_duration_loaddata = ((t_end.tv_sec + ((double) t_end.tv_nsec / 1000000000)) - (t_start.tv_sec + ((double) t_start.tv_nsec / 1000000000)));


	/* Perform epistasis detection search */

	clock_gettime(CLOCK_MONOTONIC, &t_start);	// start time measurement

	printf("#SNPs: %d\n", h_numSNPs);
	unsigned long long h_numCombinations = binomial(h_numSNPs, 3); 
	printf("#combinations: %llu\n", h_numCombinations);

	auto h_lgammaPrecalc = std::vector<float>(lgammaPrecalcSize + 1);

	for(int i=1; i < (lgammaPrecalcSize + 2); i++) {	// precalculates lgamma() values
		h_lgammaPrecalc[i - 1] = lgamma((double)i);
	}



        float* h_output = (float*) malloc(sizeof(float) * NUM_TILES * 6 * numIPUs);
        unsigned long long *h_output_index = (unsigned long long *) malloc(sizeof(unsigned long long) * NUM_TILES * 6 * numIPUs);


	engine.connectStream("lgammaPrecalc", h_lgammaPrecalc.data());
	engine.connectStream("numSNPs", &h_numSNPs_per_chunk);  
	engine.connectStream("casesSize", &numCases);		
	engine.connectStream("controlsSize", &numControls);	


	for(int i=0; i < numIPUs; i++) {
		auto callback = std::make_unique<DeviceToHostCallback_ull>(poplar::ArrayRef<unsigned long long>(h_output_index  + (i * (NUM_TILES * 6)), NUM_TILES * 6));
		engine.connectStreamToCallback("output_index", i, std::move(callback));        	// middle parameter is the replicated_graph_id
	}


	for(int i=0; i < numIPUs; i++) {
		auto callback = std::make_unique<DeviceToHostCallback>(poplar::ArrayRef<float>(h_output  + (i * (NUM_TILES * 6)), NUM_TILES * 6));
		engine.connectStreamToCallback("output", i, std::move(callback));        	// middle parameter is the replicated_graph_id
	}


	float minScore = FLT_MAX;
	unsigned long long indexOfMinScore;



	/* Execute the IPU program */
	// std::cout << "Running graph program to perform epistasis detection run\n";

	int** superChunk_pairing_conf = (int**) malloc(sizeof(int*) * numIPUs);

	int initialPairingConf[3] = {0,0,0};

	for(int i=0; i<numIPUs; i++) {
		superChunk_pairing_conf[i] = (int*) malloc(sizeof(int) * 3);

		superChunk_pairing_conf[i][0] = initialPairingConf[0];
		superChunk_pairing_conf[i][1] = initialPairingConf[1];
		superChunk_pairing_conf[i][2] = initialPairingConf[2];

		next_permutation(initialPairingConf, numSuperChunks);
	}
	
	int best_superChunk_pairing_conf[3] = {-1, -1, -1};
	int best_tile = -1;


	int numSuperChunkCombsProcessed = 0;

	while(1) {

		/*
		for(int i=0; i < numIPUs; i++) {
			printf("%d: superchunks triplet: %d,%d,%d\n", i, superChunk_pairing_conf[i][0], superChunk_pairing_conf[i][1], superChunk_pairing_conf[i][2]);
		}
		*/

		for(int i=0; i < numMiniChunks; i++) {


			// number of elements to copy
			unsigned int numElements = 2 * h_numSNPs_per_chunk * h_casesSize;

			for(int ipu_i = 0; ipu_i < numIPUs; ipu_i++) {	// for each replica 'ipu_i'

				auto callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetCases.data() + (((superChunk_pairing_conf[ipu_i][0] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_casesSize)), numElements));
				engine.connectStreamToCallback("datasetCases_superChunk1_" + std::to_string(i), ipu_i, std::move(callback));	

				callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetControls.data() + (((superChunk_pairing_conf[ipu_i][0] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_controlsSize)), numElements));
				engine.connectStreamToCallback("datasetControls_superChunk1_" + std::to_string(i), ipu_i, std::move(callback));   


				callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetCases.data() + (((superChunk_pairing_conf[ipu_i][1] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_casesSize)), numElements));
				engine.connectStreamToCallback("datasetCases_superChunk2_" + std::to_string(i), ipu_i, std::move(callback));   

				callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetControls.data() + (((superChunk_pairing_conf[ipu_i][1] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_controlsSize)), numElements));
				engine.connectStreamToCallback("datasetControls_superChunk2_" + std::to_string(i), ipu_i, std::move(callback));   


				callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetCases.data() + (((superChunk_pairing_conf[ipu_i][2] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_casesSize)), numElements));
				engine.connectStreamToCallback("datasetCases_superChunk3_" + std::to_string(i), ipu_i, std::move(callback));   

				callback = std::make_unique<HostToDeviceCallback>(poplar::ArrayRef<unsigned long long>(hDatasetControls.data() + (((superChunk_pairing_conf[ipu_i][2] * numMiniChunks) + i) * (2 * h_numSNPs_per_chunk * h_controlsSize)), numElements));
				engine.connectStreamToCallback("datasetControls_superChunk3_" + std::to_string(i), ipu_i, std::move(callback));   

			}

		}

		engine.run();	// Executes the IPU program with inputs selected
		numSuperChunkCombsProcessed += numIPUs;	

		// Finds minimum score and its index
		for(int i=0; i<numIPUs; i++) {
			for(int j=0; j < (NUM_TILES * 6); j++) {

				// printf("score: %f\n", h_output[i * (NUM_TILES * 6) + j]);

				if(h_output[i * (NUM_TILES * 6) + j] <= minScore) {
					best_superChunk_pairing_conf[0] = superChunk_pairing_conf[i][0];
					best_superChunk_pairing_conf[1] = superChunk_pairing_conf[i][1];
					best_superChunk_pairing_conf[2] = superChunk_pairing_conf[i][2];

					best_tile = floor(((double) j) / 6);	// divided by 6 because each tile has 6 workers

					minScore = h_output[i * (NUM_TILES * 6) + j];
					indexOfMinScore = h_output_index[i * (NUM_TILES * 6) + j];  

				}

			}
		}

		// Determines next triplet of superchunks to be processed
		for(int i=0; i<numIPUs; i++) {

			int ret;
			for(int j=0; j < numIPUs; j++) {

				ret = next_permutation(superChunk_pairing_conf[i], numSuperChunks);  

			}
			if((i==0) && (ret != 0)) {
				goto end_of_loop;       // if Replica 0 has no next 3-way pairing of super-chunks to process, then no Replica has anything more to do

			}
		}
	}

	end_of_loop:

	int phaseNumber = indexOfMinScore >> 48 & 0xFFFF;

        int pairing_conf[3] = {0,0,0};

	// for retrieving the mini-chunk pairing
	// this executes very fast on the host, so no overhead from this
	for(unsigned computePhase_i = 0; computePhase_i <= phaseNumber; computePhase_i++) {

		for(int tile_i=0; tile_i < NUM_TILES; tile_i++) {

			if((computePhase_i == phaseNumber) && (tile_i == best_tile))
			{
				break;	// this break statement will also get out of the outter-loop, since we are at the last iteration of the outer-loop
			}
			else {
				int ret = next_pairing(pairing_conf, numMiniChunks);  
			}

		}
	}


	printf("Solution with best score\n");

	int snpIdxInsideMini_X = (indexOfMinScore >> 0 & 0xFFFF);	
	int snpIdxInsideMini_Y = (indexOfMinScore >> 16 & 0xFFFF);
	int snpIdxInsideMini_Z = (indexOfMinScore >> 32 & 0xFFFF);

	printf("%u ", (best_superChunk_pairing_conf[0] * superchunkSize) + (pairing_conf[0] * h_numSNPs_per_chunk) + snpIdxInsideMini_X);
	printf("%u ", (best_superChunk_pairing_conf[1] * superchunkSize) + (pairing_conf[1] * h_numSNPs_per_chunk) + snpIdxInsideMini_Y);
	printf("%u ", (best_superChunk_pairing_conf[2] * superchunkSize) + (pairing_conf[2] * h_numSNPs_per_chunk) + snpIdxInsideMini_Z);

	printf("%.2lf \n", minScore);


	clock_gettime(CLOCK_MONOTONIC, &t_end);
	double timing_duration_search = ((t_end.tv_sec + ((double) t_end.tv_nsec / 1000000000)) - (t_start.tv_sec + ((double) t_start.tv_nsec / 1000000000)));


        printf("Time to setup IPU graph/program:\t%0.3lf seconds\n", timing_duration_compilingOrLoadingExe);	
        printf("Time to load+preprocess dataset:\t%0.3lf seconds\n", timing_duration_loaddata);
	printf("Time to perform epistasis search:\t%0.3lf seconds\n", timing_duration_search); 

	std::cout << "Tera unique sets per sec. scaled to sample size: " << std::fixed << std::setprecision(6) << (((double) h_numCombinations * (double) (sampleSize) / (double)(timing_duration_search)) / 1e12) << std::endl;


	return 0;
}
