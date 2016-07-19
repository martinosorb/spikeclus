#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// Workaround to ease the translation from C# to C++
#define allocate2D(V, t, x, y) \
    V = new t*[x]; for (int i = 0; i < x; i++) { V[i] = new t[y]; }

namespace SpkDslowFilter {
class InterpDetection {
	int NChannels;// number of channels; is set when reading the data
	int** ChInd4a;
	int** ChInd4b;
	int** ChInd4c;
	//int** ChInd4d;
	int* ChInd4List;
	int* ChInd5;
	int** ChInd5a;
	int** ChInd5b;
	int** ChInd5c;
	//int** ChInd5d;
	int* ChInd5List;
	//Sampling rate and subsampling
	int Sampling;//=7000; set at runtime
	float FrameDt;//interval between frames. for gui.
	int sfi;// integer sampling rate
	bool HF;//frequency above 12kHz? Used for slow filter.
	bool Subsample;//whether only every second frame shall be used (used for slow filter).
	int* dfTI;//to pass parameters (increment and offsets) from the gui to the main function.
	int df;//increment for detection (used for backward detection subsampling), defined in initial estimation
	int dfAbs;//absolute value (of df)
	int dfSign;//sign (of df)
	long t0x;//where the first cutout starts
	int tx;//where to start the iteration initially (relative to cutout)
	int ti;//initial point to estimate (relative to cutout)
	int tf;//final point to estimate (relative to cutout)
	int tm;//where to start Aglobal estimation
	int tms;
	int tfiA;
	int ty;//where to end Aglobal estimation
	//slow filter
	int** Avgs1;
	int** Avgs2;
	int* Avgs3;
	int* vSmth;
	int tSmth;
	int tSmth1;
	int tdSmth;
	int tdSmth1;
	static const int TauFilt1 = 4;//i.e. 4*2=8 frames
	int TauFiltOffset;// = 48;//16 for smooth update (Taufilt=32) and 32 for overlapping, 64 for nonoverlapping. set at runtime.
	int NFblocks; //average over two, weighted; want ~5 ms take Sampling/5000; set at runtime if >20kHz
	static const int TauFilt=6*TauFilt1;//update timescale of slow offset voltage QmPreD in frames (here using 2*TauFilt1*TauFilt2)
	int FiltNorm;// = 18; set at runtime
	static const int dtTMx = 3;
	//global fluctuations
	int* Aglobal;
	int* Aglobaldiff;
	int* AglobalSmth;
	int* AglobalSdiff;
	//int Aglobaldiffold=0;
	int* Slice;
	//correlation with global fluctuations
	long SqIglobal;//sum of squared global increments
	long* SqIv;//sum of squared channel increments
	long** SIprod;//sum of product of global and channel voltage increments
	int* Vbias;
	int* FVbias;
	long* FVsbias;
	long* Vsqbias;
	bool ACF;//measure autocorrelation? (otherwise faster)
	//Variables for variance and mean
	int* Qd;//noise amplitude
	int* Qm;//median
	int* QmPre;//filtered raw
	int* QmPreD;//offset after slow timescale correction
	int** Qdiff;//signal amplitude
	int** Qmax;
	int* QmaxE;
	//Variables for the spike detection
	int* Sl4;//counter for spike length
	int* Z4next;//counter, skip frames when far from detection threshold
	int* Z5next;
	bool* Sl4x;//tag for removed spikes
	bool* AHP4;//tag for repolarizing current
	int* Amp4;//buffers spike amplitude
	int NSpikes4;//to count spikes
	//Variables for the spike detection
	int* Sl5;//counter for spike length
	bool* Sl5x;//tag for removed spikes
	bool* AHP5;//tag for repolarizing current
	int* Amp5;//buffers spike amplitude
	int NSpikes5;//to count spikes
	//cutouts
	int CutPre;
	int CutPost;
	int CutOffset;// = 13; set at runtime
	int tCut;//cutout interval=20; set at runtime
	int tCutLong;//=27; set at runtime
	int* tQm;
	int* tQmLong;
	int* tQmA;
	int* tQmX;//cutout indices
	int* tQmXLong;//cutout indices
	int* tShape;
	int tQmi;//initial cutout
	int tQmf;//final cutout
	int tQm0;//
	int tQmm;//median
	int Lspike;
	int Lsw;
	int Li;
	int Lf;
	int Lm;
	int* Lw;
	int Lmx;
	//Parameters for variance and mean updates
	static const int Tau_m0 = 4;//timescale for updating Qm (increment is Qd/Tau_m)
	static const int Qdmin=300;//set minimum value of Qd
	//Parameters for spike detection
	int threshold;// = 6*AmpScale;//threshold to detect spikes >6 is likely to be real spikes; set at runtime
	int AHPthr;// = 0;//signal should go below that threshold within Slmax-Slmin frames; set at runtime
	int Slmax;// = 8;//dead time in frames after peak, used for further testing; set at runtime
	int Sln0;//=2; set at runtime
	//Parameters for reading data
	static const int tInc = 192;//increment for reading data
	static const int Ascale0 = -192;//factor to multiply to raw traces to increase resolution; definition of ADC counts had been changed!
	int Ascale;
	int AscaleV;
	int AscaleG;
	static const int AmpScale=100;//want to use integer amplitudes
	//Parameters for recalibration events and artefact handling
	static const int artT =10;//to use after artefacts; to update Qm for 10 frames
	int* A;//control parameter for amplifier effects
	//instead of copying data around, want lookup variables that loop over one dimension of a matrix
	//2-loop
	int dt;
	static const int dtMx = 2;
	int dtPre;
	//(dtEMx-1)-loop --> still exists/necessary?
	int dtE;
	static 	const int dtEMx = 2;
	int dtEx;

	int ChInd4ListSize;
	int ChInd5ListSize;
	//Files to save the spikes etc.
	std::ofstream w; //for spikes
	std::ofstream wShapes; //for raw data
	std::ofstream wX; //for spikes
	std::ofstream wShapesX; //for raw data
	std::ofstream wInfo; //for other stuff
	std::ofstream wMean; //for avg. Voltage

	static const int recalibTrigger=1;
	int Acal;//for recalibration events
		
public:
	InterpDetection();
	~InterpDetection();
	int* SetInitialParams (long nFrames, double nSec, int sf, double sfd, int NCh, int* Indices);
	void openFiles(const std::string& name);
	void AvgVoltageDefault(unsigned short* vm, long t0, int t, int tInc); //want to compute an approximate 33 percentile
	void InitialEstimation(unsigned short* vm, long t0); //use this to get a better initial estimate of Qd. only fast transients.
	void StartDetection(unsigned short* vm, long t0, long nFrames, double nSec, double sfd, int* Indices);
	void skipLastReverse(int skipLast);
	void Iterate(unsigned short* vm, long t0, int tInc); // tInc shadowed to keep consistency with the vm array (stride size)
	void FinishDetection(unsigned short* vm, int skipLast, int tInc); // tInc is also shadowed here
};
};