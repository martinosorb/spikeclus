#include "SpkDonline.h"

namespace SpkDonline {
Detection::Detection() :
	// Set default parameters
	threshold(9), // threshold to detect spikes >11 is likely to be real spikes, but can and should be sorted afterwards
	AHPthr(0),    // signal should go below that threshold within MaxSl-Slmin frames
	MaxSl(8),     // dead time in frames after peak, used for further testing
	MinAvgAmp(5), // minimal avg. amplitude of peak (in units of Qd)
	MinSl(3)     // length considered for determining avg. spike amplitude
{}

void Detection::InitDetection(long nFrames, double nSec, int sf, int NCh,
                              long ti, long *Indices, unsigned int nCPU) {
  NChannels = NCh;
  tInc = ti;
  Qd = new int[NChannels];      // noise amplitude
  Qm = new int[NChannels];      // median
  Sl = new int[NChannels];      // counter for spike length
  AHP = new bool[NChannels];    // counter for repolarizing current
  Amp = new int[NChannels];     // buffers spike amplitude
  SpkArea = new int[NChannels]; // integrates over spike
  A = new int[NChannels];       // control parameter for amplifier effects
  ChInd = new int[NChannels];
  Slice = new int[NChannels];
  MaxSl = sf / 1000 + 1;
  MinSl = sf / 3000 + 2;
  Sampling = sf;
  Aglobal = new int[tInc];
  for (int i = 0; i < tInc; i++)
    Aglobal[i] = 0;
  for (int i = 0; i < NChannels; i++) {
    Qd[i] = 400;
    Qm[i] = Voffset * Ascale;
    Sl[i] = 0;
    AHP[i] = false;
    Amp[i] = 0;
    A[i] = artT; // start like after an out-of-linear-regime event
    SpkArea[i] = 0;
    ChInd[i] = Indices[i];
  }
  
  nthreads = nCPU;
  threads = new std::thread[nthreads];
}

void Detection::SetInitialParams(int thres, int maa, int ahpthr, int maxsl,
                                 int minsl) {
  // set the detection parameters
  threshold = thres;
  MinAvgAmp = maa;
  AHPthr = ahpthr;
  MaxSl = maxsl;
  MinSl = minsl;
}

// don't know how to make multiple threads write to the same file,
// maybe you'll have to buffer values during the detection (make Iterate a list
// instead of a void)
// and then write spikes for each block
void Detection::openSpikeFile(const char *name) {
  std::cout << "# Writing to: " << name << "\n";
  w.open(name);
  // fs = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
  // w = new StreamWriter(fs);
}

void Detection::MedianVoltage(unsigned short *vm) // easier to interpret, though
                                                  // it takes a little longer to
                                                  // run, but I'm sure there is
                                                  // a faster method in C++ for
                                                  // computing the median
{ // otherwise could try to take the mean also (have to ignore channels out of
  // the linear regime then) as signals are about 15% correlated
  for (int t = 0; t < tInc; t++) { // this function wastes most of the time
    for (int i = 0; i < NChannels; i++) { // loop across channels
      Slice[i] = vm[i + t*NChannels];        // vm [i] [t];
    }
    std::sort(Slice, Slice + NChannels);
    Aglobal[t] = Slice[NChannels / 2];
  }
}

void Detection::MeanVoltageThread(int threadID, unsigned short *vm) {

  int chunkSize = std::ceil( (float) tInc/ (float) nthreads);

  int n;
  int Vsum;
  for (int t = threadID*chunkSize; t < tInc and t < (threadID+1)*chunkSize; t++) { 
    n = 1; // constant offset doesn't matter, avoid zero division
    Vsum = 0;
    for (int i = 0; i < NChannels; i++) { // loop across channels
      if (((vm[i + t*NChannels] + 4) % NChannels) > 10) {
        Vsum += (vm[i + t* NChannels ]); // !!! Indexing changed
        n++;
      }
    }
    Aglobal[t] = Vsum / n;
  }

}

void Detection::MeanVoltage(unsigned short *vm) {
  for (int threadID = 0; threadID < nthreads; threadID++) {
    threads[threadID] = std::thread( [=] { MeanVoltageThread(threadID, vm); });
  }
  for (int threadID = 0; threadID < nthreads; threadID++) { 
    threads[threadID].join();
  }
}

void Detection::IterateThread(int threadID, unsigned short *vm, long t0) {
  
  int a; // buffer for Iterate() now is thread dependant
  
  // Loop over data (offline algorithm)
  for (int t = 0; t < tInc; t++) { 

  // Number of channels associated to a thread
  int chunkSize = std::ceil( (float) NChannels/ (float) nthreads);

  // Loop accross all channels associated to this thread
  for (int i = threadID*chunkSize; i < NChannels and i < (threadID+1)*chunkSize; i++) { 

    // CHANNEL OUT OF LINEAR REGIME
    if (((vm[i + t*NChannels] + 4) % NChannels) < 10) {
      if (A[i] < artT) { // reset only when it starts leaving the linear regime
        Sl[i] = 0;
        A[i] = artT;
      }
    }
    // DEFAULT OPERATIONS
    else if (A[i] == 0) {
      // Difference between ADC counts and Qm
      a = (vm[i + t*NChannels] - Aglobal[t]) * Ascale - Qm[i]; 
      // UPDATE Qm and Qd
      if (a > 0) {
        if (a > Qd[i]) {
          Qm[i] += Qd[i] / Tau_m0;
          if (a < (5 * Qd[i])) {
            Qd[i]++;
          } else if ((Qd[i] > Qdmin) & (a > (6 * Qd[i]))) {
            Qd[i]--;
          }
        } else if (Qd[i] > Qdmin) { // set a minimum level for Qd
          Qd[i]--;
        }
      } else if (a < -Qd[i]) {
        Qm[i] -= Qd[i] / Tau_m0 / 2;
      }
      // TREATMENT OF THRESHOLD CROSSINGS
      if (Sl[i] > 0) { // Sl frames after peak value
        // default
        Sl[i] = (Sl[i] + 1) % (MaxSl + 1); // increment Sl[i]
        if (Sl[i] < MinSl) { // calculate area under first and second frame
                              // after spike
          SpkArea[i] += a;
        }
        // check whether it does repolarize
        else if (a < (AHPthr * Qd[i])) {
          AHP[i] = true;
        }
        // accept spikes after MaxSl frames if...
        if ((Sl[i] == MaxSl) & (AHP[i])) {
          if ((2 * SpkArea[i]) > (MinSl * MinAvgAmp * Qd[i])) {

            output_mtx.lock();

            w << ChInd[i] << " " << t0 + t - MaxSl + 1 << " "
              << -Amp[i] * Ascale / Qd[i] << "\n";

            output_mtx.unlock();
            
          }
          Sl[i] = 0;
        }
        // check whether current ADC count is higher
        else if (Amp[i] < a) {
          Sl[i] = 1; // reset peak value
          Amp[i] = a;
          AHP[i] = false;  // reset AHP
          SpkArea[i] += a; // not resetting this one (anyway don't need to
                            // care if the spike is wide)
        }
      }
      // check for threshold crossings
      else if (a > ((threshold * Qd[i]) / 2)) {
        Sl[i] = 1;
        Amp[i] = a;
        AHP[i] = false;
        SpkArea[i] = a;
      }
    }
    // AFTER CHANNEL WAS OUT OF LINEAR REGIME
    else {
      Qm[i] = (2 * Qm[i] + (vm[i + t*NChannels] - Aglobal[t]) * Ascale +
                2 * Qd[i]) /
              3; // update Qm
      A[i]--;
    }
  }
  
  }
}

void Detection::Iterate(unsigned short *vm, long t0) {
  // SPIKE DETECTION  
  for (int threadID = 0; threadID < nthreads; threadID++) {
    threads[threadID] = std::thread( [=] { IterateThread(threadID, vm, t0); });
  }
  for (int threadID = 0; threadID < nthreads; threadID++) { 
    threads[threadID].join();
  }
}

void Detection::FinishDetection() // write spikes in interval after last
                                  // recalibration; close file
{
  w.close();
}
}
