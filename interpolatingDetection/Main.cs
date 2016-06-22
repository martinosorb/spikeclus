using System;
using System.Runtime.InteropServices;
using BW;
using System.IO;
using System.Collections;
using System.Windows.Forms;
using System.Diagnostics;
using System.Threading;
using System.Drawing;

namespace SpkDonline {
	class Detection {
		int NChannels;// number of channels; is set when reading the data
		int[] ChInd;//indices for parallelization
		//Variables for variance and mean
		int[] Qd;//noise amplitude
		int[] Qm;//median
		//Variables for the spike detection
		int[] Sl;//counter for spike length
		bool[] AHP;//counter for repolarizing current
		int[] Amp;//buffers spike amplitude
		int[] SpkArea;//integrates over spike
		//Parameters for variance and mean updates
		const int Tau_m0 = 4;//timescale for updating Qm (increment is Qd/Tau_m)
		const int Qdmin=200;//set minimum value of Qd
		//Parameters for spike detection
		int threshold = 9;//threshold to detect spikes >11 is likely to be real spikes, but can and should be sorted afterwards
		int AHPthr = 0;//signal should go below that threshold within MaxSl-Slmin frames
		int MaxSl = 8;//dead time in frames after peak, used for further testing
		int MinAvgAmp = 5;//minimal avg. amplitude of peak (in units of Qd)
		int MinSl =3;//length considered for determining avg. spike amplitude
		//Parameters for reading data
		const int tInc = 100;//increment for reading data, has to be changed in main program as well
		const int Ascale = -64;//factor to multiply to raw traces to increase resolution; definition of ADC counts had been changed!
		const int Voffset =0;//mean ADC counts, as initial value for Qm
		//Parameters for recalibration events and artefact handling
		const int artT =10;//to use after artefacts; to update Qm for 10 frames
		int[] A;//control parameter for amplifier effects
		//Files to save the spikes etc.
		int Sampling=7000;
		static FileStream fs;//for spikes
		StreamWriter w;
		int[] Aglobal= new int[tInc];
		int[] Slice;


		public void InitDetection (long nFrames, double nSec, int sf, int NCh, int[] Indices)
		{
			NChannels = NCh;
			Qd = new int[NChannels];//noise amplitude
			Qm = new int[NChannels];//median
			Sl = new int[NChannels];//counter for spike length
			AHP = new bool[NChannels];//counter for repolarizing current
			Amp = new int[NChannels];//buffers spike amplitude
			SpkArea = new int[NChannels];//integrates over spike
			A = new int[NChannels];//control parameter for amplifier effects
			ChInd = new int[NChannels];
			Slice = new int[NChannels];
			MaxSl = sf / 1000 + 1;
			MinSl = sf / 3000 + 2;
			Sampling = sf;
			for (int i=0;i<NChannels;i++){
				Qd[i]=400;
				Qm[i]=Voffset*Ascale;
				Sl[i]=0;
				AHP[i]=false;
				Amp[i]=0;
				A[i]=artT;//start like after an out-of-linear-regime event
				SpkArea[i]=0;
				ChInd[i]=Indices[i];
			}
			w.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
		}

		public void SetInitialParams ()
		{
			// ask for some paramters
			Form form1 = new Form ();
			Button button1 = new Button ();
			button1.Text = "Go!";
			button1.Location = new Point (10, 10);
			button1.DialogResult = DialogResult.OK;
			form1.Text = "Detection Parameters";
			form1.FormBorderStyle = FormBorderStyle.FixedDialog;
			form1.AcceptButton = button1;
			form1.StartPosition = FormStartPosition.CenterScreen;
			form1.Controls.Add (button1);
			form1.SetBounds(0,0,300,400);
			NumericUpDown numericUpDown1 = new NumericUpDown ();
			numericUpDown1.Parent = form1;
			numericUpDown1.Maximum = 20;
			numericUpDown1.Minimum = 9;
			numericUpDown1.Value = 12;
			numericUpDown1.Location = new Point (10, 80);
			Label l1 = new Label ();
			l1.Text = "Detection Threshold";
			l1.AutoSize = true;
			l1.Location = new Point (numericUpDown1.Left, numericUpDown1.Bottom);
			l1.Size = numericUpDown1.Size;
			l1.Parent = form1;
			NumericUpDown numericUpDown2 = new NumericUpDown ();
			numericUpDown2.Parent = form1;
			numericUpDown2.Maximum = 15;//should be lower than the threshold
			numericUpDown2.Minimum = 0;
			numericUpDown2.Increment = 1;
			numericUpDown2.Value = 7;
			numericUpDown2.Location = new Point (10, 140);
			Label l2 = new Label ();
			l2.Text = "minimum avg. Depolarization Amplitude (current)";//same units as threshold, but integrated
			l2.AutoSize = true;
			l2.Location = new Point (numericUpDown2.Left, numericUpDown2.Bottom);
			l2.Size = numericUpDown2.Size;
			l2.Parent = form1;
			NumericUpDown numericUpDown3 = new NumericUpDown ();//ideally want to have multiples of frames here (converted to ms)
			numericUpDown3.Parent = form1;
			numericUpDown3.Maximum = 1.0m;
			numericUpDown3.Minimum = 0.0m;
			numericUpDown3.Increment = 0.1m;
			numericUpDown3.Value = 0.4m;
			numericUpDown3.DecimalPlaces=1;
			numericUpDown3.Location = new Point (10, 200);
			Label l3 = new Label ();
			l3.Text = "minimum Depolarization Width (ms)";
			l3.AutoSize = true;
			l3.Location = new Point (numericUpDown3.Left, numericUpDown3.Bottom);
			l3.Size = numericUpDown3.Size;
			l3.Parent = form1;
			NumericUpDown numericUpDown4 = new NumericUpDown ();//ideally want to have multiples of frames here (converted to ms)
			numericUpDown4.Parent = form1;
			numericUpDown4.Maximum = 2.0m;
			numericUpDown4.Minimum = 0.5m;
			numericUpDown4.Increment = 0.1m;
			numericUpDown4.Value = 1.0m;
			numericUpDown4.DecimalPlaces=1;
			numericUpDown4.Location = new Point (10, 260);
			Label l4 = new Label ();
			l4.Text = "maximum Depolarization Width (ms)";//actually only measure from peak here
			l4.AutoSize = true;
			l4.Location = new Point (numericUpDown4.Left, numericUpDown4.Bottom);
			l4.Size = numericUpDown4.Size;
			l4.Parent = form1;
			NumericUpDown numericUpDown5 = new NumericUpDown ();
			numericUpDown5.Parent = form1;
			numericUpDown5.Maximum = 10;//will include voltage steps if I set it that high...
			numericUpDown5.Minimum = -5;
			numericUpDown5.Value = 0;
			numericUpDown5.Location = new Point (10, 320);
			Label l5 = new Label ();
			l5.Text = "Repolarization threshold";//again same units as threshold
			l5.AutoSize = true;
			l5.Location = new Point (numericUpDown5.Left, numericUpDown5.Bottom);
			l5.Size = numericUpDown5.Size;
			l5.Parent = form1;
			form1.ShowDialog ();
			// set the parameters
			threshold = (int)numericUpDown1.Value;
			MinAvgAmp = (int)numericUpDown2.Value;
			AHPthr = (int)numericUpDown5.Value;
			MaxSl = (int)(Sampling*numericUpDown4.Value/ 1000 +0.5m) + 1;
			MinSl = (int)(Sampling*numericUpDown3.Value/ 1000 +0.5m);
			form1.Dispose();
		}

		//don't know how to make multiple threads write to the same file, 
		//maybe you'll have to buffer values during the detection (make Iterate a list instead of a void)
		//and then write spikes for each block
		public void openSpikeFile (string name) {
			fs = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			w = new StreamWriter(fs);
		}

		public void MedianVoltage (short[][] vm)//easier to interpret, though it takes a little longer to run, but I'm sure there is a faster method in C++ for computing the median
		{//otherwise could try to take the mean also (have to ignore channels out of the linear regime then) as signals are about 15% correlated
			for (int t=0; t<tInc; t++) {//this function wastes most of the time
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice [i] = (vm [i] [t]);
				}
				Array.Sort (Slice);
				Aglobal[t] = Slice [NChannels / 2];
			}
		}

		public void MeanVoltage (short[][] vm)// if median takes too long...
		{
			int n;
			int Vsum;
			for (int t=0; t<tInc; t++) {
				n=1;//constant offset doesn't matter, avoid zero division
				Vsum=0;
				for (int i=1; i<NChannels; i++) {//loop across channels
					if (((vm[i][t]+4)%4096)>10) {
						Vsum+= (vm [i] [t]);
						n++;
					}
				}
				Aglobal[t] = Vsum/n;
			}
		}

		public void Iterate (short[][] vm, long t0){
			int a;//to buffer the difference between ADC counts and Qm
			for (int t=0; t<tInc; t++) {//loop over data, will be removed for an online algorithm
				// SPIKE DETECTION
				for (int i=0; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm[i][t]+4)%4096)<10) {
						if (A[i]<artT) {//reset only when it starts leaving the linear regime
							Sl[i]=0;
							A[i]=artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A[i]==0) {
						a = (vm[i][t]-Aglobal[t])*Ascale-Qm[i];//difference between ADC counts and Qm
						//UPDATE Qm and Qd
						if (a>0) {
							if (a>Qd[i]) {
								Qm[i]+=Qd[i]/Tau_m0;
								if  (a<(5*Qd[i])) {
									Qd[i]++;
								}
								else if ((Qd[i]>Qdmin) & (a>(6*Qd[i]))) {
									Qd[i]--;
								}
							}
							else if (Qd[i]>Qdmin){//set a minimum level for Qd
								Qd[i]--;
							}
						}
						else if (a<-Qd[i]){
							Qm[i]-=Qd[i]/Tau_m0/2;

						}
						//TREATMENT OF THRESHOLD CROSSINGS
						if (Sl[i]>0) {//Sl frames after peak value
							//default
							Sl[i]=(Sl[i]+1)%(MaxSl+1);// increment Sl[i]
							if (Sl[i]<MinSl){//calculate area under first and second frame after spike
								SpkArea[i]+=a;
							}
							//check whether it does repolarize
							else if (a<(AHPthr*Qd[i])) {
								AHP[i]=true;
							}
							//accept spikes after MaxSl frames if...
							if ((Sl[i]==MaxSl) & (AHP[i])) {
								if ((2*SpkArea[i])>(MinSl*MinAvgAmp*Qd[i])){
									w.WriteLine("{0} {1} {2}", ChInd[i], t0+t-MaxSl+1, -Amp[i]*Ascale/Qd[i]);
								}
								Sl[i]=0;
							}
							//check whether current ADC count is higher
							else if (Amp[i]<a) {
								Sl[i]=1;//reset peak value
								Amp[i]=a;
								AHP[i]=false;//reset AHP
								SpkArea[i]+=a;//not resetting this one (anyway don't need to care if the spike is wide)
							}
						}
						//check for threshold crossings
						else if (a>((threshold*Qd[i])/2)) {
							Sl[i]=1;
							Amp[i]=a;
							AHP[i]=false;
							SpkArea[i]=a;
						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm[i]=(2*Qm[i]+(vm[i][t]-Aglobal[t])*Ascale+2*Qd[i])/3;//update Qm
						A[i]--;
					}
				}
			}
		}

		public void FinishDetection (){//write spikes in interval after last recalibration; close file
			w.Flush();
			w.Close();
			fs.Close();
		}
	}
	public class MainClass{
		[STAThread]
		static public void Main ()
		{
			// Initialize an instance of the BrwRdr class
			BrwRdr brwRdr = new BrwRdr ();
			// you must use a valid full file path

			// Create and display a fileChooserDialog
			OpenFileDialog openFileDialog1 = new OpenFileDialog ();
			openFileDialog1.Filter = "brw files (*.brw)|*.brw|All files (*.*)|*.*";
			openFileDialog1.FilterIndex = 1;
			openFileDialog1.RestoreDirectory = true;
			openFileDialog1.InitialDirectory = Directory.GetCurrentDirectory ();

			if (openFileDialog1.ShowDialog () == DialogResult.OK) {

				// Open the file for reading.
				brwRdr.Open (openFileDialog1.FileName);
				// the number of frames of the recording
				long nFrames = brwRdr.RecNFrames;
				// the duration in seconds of the recording
				double nSec = brwRdr.RecDuration / 1000;
				// the sampling frequency
				double sfd = (double) brwRdr.SamplingRate;
				int sf = (int) sfd;
				Console.WriteLine ("# Number of frames: {0}", nFrames);
				Console.WriteLine ("# Duration (s): {0}", nSec);
				Console.WriteLine ("# Sampling rate: {0}", sfd);

				SaveFileDialog saveFileDialog1 = new SaveFileDialog ();
				saveFileDialog1.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
				saveFileDialog1.FilterIndex = 1;
				saveFileDialog1.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_online_Spikes";
				saveFileDialog1.DefaultExt = "txt";
				saveFileDialog1.InitialDirectory = Directory.GetCurrentDirectory ();
				saveFileDialog1.Title = "Save Spikes As";
				if (saveFileDialog1.ShowDialog () == DialogResult.OK) {
					//list of channels
					ChCoord[] Channels = brwRdr.GetRecChs (StreamType.Raw);
					int[] Indices = new int[Channels.Length];
					for (int i=0; i<Channels.Length; i++) {
						Indices [i] = (Channels [i].Col - 1) + 64 * (Channels [i].Row - 1);
					}
					// how many frames to analyze
					long nDumpFrames = nFrames;
					Detection SpkD;
					SpkD = new Detection ();
					// open the output file
					SpkD.openSpikeFile (saveFileDialog1.FileName);
					saveFileDialog1.Dispose ();
					SpkD.InitDetection (nFrames, nSec, sf, Channels.Length, Indices);
					SpkD.SetInitialParams ();

					// measure execution time
					var sw = new Stopwatch ();
					sw.Start ();

					const long tInc = 100;//has to be changed in Detection class as well
					int tCut = sf/1000 + sf/1000 + 6;
					const int NChannels = 4096;
					for (long t0=0; t0<nDumpFrames-tInc; t0+=tInc-tCut) {
						if ((t0/tInc)%100==0){
							Console.WriteLine ("{0} sec", t0/sf);
						}
						short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
						SpkD.MedianVoltage(vm);
						SpkD.Iterate (vm, t0);
					}
					SpkD.FinishDetection ();//should be a different object when parallelizing, just closes files.
					sw.Stop ();
					Console.WriteLine ("Elapsed time: {0}", sw.Elapsed); // TimeSpan
					Console.WriteLine ("Milliseconds/frame: {0}", sw.Elapsed.TotalMilliseconds / nFrames); // TimeSpan
				}
			}
		}
	}
}



