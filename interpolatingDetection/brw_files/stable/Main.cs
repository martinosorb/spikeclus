using System;
using System.Runtime.InteropServices;
using BW;
using System.IO;
using System.Collections;
using System.Windows.Forms;
using System.Diagnostics;
using System.Threading;
using System.Drawing;

namespace SpkD45 {
	class Detection {
		int NChannels;// number of channels; is set when reading the data
		int[,] ChInd4;
		int[,] ChInd5;
		int[] W4 = {0,2,3,1};
		int[] W5 = {0,1,1,1};
		//Variables for variance and mean
		int[] Qd;//noise amplitude
		int[] Qm;//median
		int[,] Qdiff;//signal amplitude
		int[] Qmax;
		//correlation index with global fluctuations
		long SqIglobal;//sum of squared global increments
		long[] SqIv;//sum of squared channel increments
		long[,] SIprod;//sum of product of global and channel voltage increments
		int SqIg=0;
		//int SqIgroot=0;
		int[] SIp;
		int[] Vbias;
		int[] FVbias;
		long[] FVsbias;
		long[] Vsqbias;
		//Variables for the spike detection
		int[] Sl4;//counter for spike length
		bool[] Sl4x;//tag for removed spikes
		bool[] AHP4;//tag for repolarizing current
		int[] Amp4;//buffers spike amplitude
		int[] Qd4;//noise amplitude
		int[,] Sx4;
		int[,] Sx5;
		//Variables for the spike detection
		int[] Sl5;//counter for spike length
		bool[] Sl5x;//tag for removed spikes
		bool[] AHP5;//tag for repolarizing current
		int[] Amp5;//buffers spike amplitude
		int[] Qd5;//noise amplitude
		//Parameters for variance and mean updates
		const int Tau_m0 = 4;//timescale for updating Qm (increment is Qd/Tau_m)
		const int Qdmin=450;//set minimum value of Qd
		//Parameters for spike detection
		int threshold5 = 10;//threshold to detect spikes >11 is likely to be real spikes, but can and should be sorted afterwards
		int threshold4 = 10;//threshold to detect spikes >11 is likely to be real spikes, but can and should be sorted afterwards
		int AHPthr = 0;//signal should go below that threshold within Slmax-Slmin frames
		int Slmax = 8;//dead time in frames after peak, used for further testing
		//int Slmin =3;// number of frames that should be below Amp-AHPthr; need this now for excluding the detection of multiple spikes
		int Sln0=2;
		//Parameters for reading data
		//const int tInc = 128;//increment for reading data, has to be changed in main program as well
		const int Ascale = -64;//factor to multiply to raw traces to increase resolution; definition of ADC counts had been changed!
		//Parameters for recalibration events and artefact handling
		const int artT =10;//to use after artefacts; to update Qm for 10 frames
		int[] A;//control parameter for amplifier effects
		//Files to save the spikes etc.
		int CutOffset = 13;
		int CutAfter = 7;
		int CutAfterLong =14;//use for spikes that do not repolarize properly (compound events?)
		int Sampling=7000;
		static FileStream fs;//for spikes
		StreamWriter w;
		static FileStream fsShapes;//for raw data
		StreamWriter wShapes;
		static FileStream fsX;//for spikes
		StreamWriter wX;
		static FileStream fsShapesX;//for raw data
		StreamWriter wShapesX;
		static FileStream fsInfo;//for other stuff
		StreamWriter wInfo;
		static FileStream fsMean;//for avg. Voltage
		StreamWriter wMean;
		int Aglobal=0;
		int Aglobaldiff=0;
		int Aglobaldiffold=0;
		int recalibTrigger=1;
		int Acal=3000;//for recalibration events
		// gui elements
		Form logwindow;
		TextBox logbox;
		int[] Slice;
		int dt=0;
		
		public void SetInitialParams (long nFrames, double nSec, int sf, int NCh, int[] Indices)
		{
			NChannels = NCh;
			int[] SInd = new int[4096];
			int[] SInd4 = new int[4096];
			int[] SInd5 = new int[4096];
			ChInd4 = new int[NCh, 16];
			ChInd5 = new int[NCh, 13];
			Sx4 = new int[NCh,4];
			Sx5 = new int[NCh,4];
			for (int i=0; i<NCh; i++) {//fillvalues
				SInd[i]=-1;
				SInd4[i]=0;
				SInd5[i]=0;
				for (int j=0; j<4; j++) {
					Sx4[i,j]=0;
					Sx5[i,j]=0;
				}
				for (int j=0; j<13; j++) {
					ChInd4 [i, j] = -1;
					ChInd5 [i, j] = -1;
				}
				ChInd4 [i, 13] = -1;
				ChInd4 [i, 14] = -1;
				ChInd4 [i, 15] = -1;
			}
			for (int i=0; i<NCh; i++) {//find active channels and number of neighbors
				SInd [Indices [i]] = i;
				SInd4 [Indices [i]] += 1;
				SInd5 [Indices [i]] += 1;
				if ((Indices [i] % 64) >0) {//channel to the left
					SInd4 [Indices [i] - 1] += 1;
					SInd5 [Indices [i] - 1] += 1;
					if ((Indices [i] / 64) >0) {//up and left
						SInd4 [Indices [i] - 65] += 1;
					}
				}
				if ((Indices [i] % 64) <63) {//right
					SInd5 [Indices [i] + 1] += 1;
				}
				if ((Indices [i] / 64) >0) {//up
					SInd4 [Indices [i] - 64] += 1;
					SInd5 [Indices [i] - 64] += 1;
				}
				if ((Indices [i] / 64) <63) {//down
					SInd5 [Indices [i] + 64] += 1;
				}
			}//channels with (all) neighbors have SInd4=4 or SInd5=5
			for (int i=0; i<4096; i++) {
				if (SInd4 [i] == 4) {
					ChInd4 [SInd [i], 0] = SInd [i];
					ChInd4 [SInd [i], 1] = SInd [i + 1];
					ChInd4 [SInd [i], 2] = SInd [i + 65];
					ChInd4 [SInd [i], 3] = SInd [i + 64];
					if (SInd5 [i] == 5) {
						ChInd4 [SInd [i], 12] = SInd [i];
					}
					if (SInd5 [(i + 1)%4096] == 5) {
						ChInd4 [SInd [i], 13] = SInd [(i + 1) % 4096];
					}
					if (SInd5 [(i + 64) % 4096] == 5) {
						ChInd4 [SInd [i], 15] = SInd [(i +64) % 4096];
					}
					if (SInd5 [(i + 65) % 4096] == 5) {
						ChInd4 [SInd [i], 14] = SInd [(i + 65) % 4096];
					}
					if (SInd4 [(i + 4032) % 4096] == 4) {
						ChInd4 [SInd [i], 4] = SInd [(i + 4032) % 4096];
						ChInd4 [SInd [i], 5] = SInd [(i + 4033) % 4096];
					}
					if (SInd4 [(i + 1) % 4096] == 4) {
						ChInd4 [SInd [i], 6] = SInd [(i + 2) % 4096];
						ChInd4 [SInd [i], 7] = SInd [(i + 66) % 4096];
					}
					if (SInd4 [(i + 64) % 4096] == 4) {
						ChInd4 [SInd [i], 8] = SInd [(i + 129) % 4096];
						ChInd4 [SInd [i], 9] = SInd [(i + 128) % 4096];
					}
					if (SInd4 [(i + 4095) % 4096] == 4) {
						ChInd4 [SInd [i], 10] = SInd [(i + 63) % 4096];
						ChInd4 [SInd [i], 11] = SInd [(i + 4095) % 4096];
					}
				}
				if (SInd5 [i] == 5) {
					ChInd5 [SInd [i], 0] = SInd [i];
					ChInd5 [SInd [i], 1] = SInd [i - 64];
					ChInd5 [SInd [i], 2] = SInd [i + 1];
					ChInd5 [SInd [i], 3] = SInd [i + 64];
					ChInd5 [SInd [i], 4] = SInd [i - 1];
					if (SInd4 [i] == 4) {
						ChInd5 [SInd [i], 11] = SInd [i];
						ChInd5 [SInd [i], 7] = SInd [(i + 65) % 4096];
					}
					if (SInd4 [(i + 4095)%4096] == 4) {
						ChInd5 [SInd [i], 12] = SInd [(i + 4095) % 4096];
						ChInd5 [SInd [i], 8] = SInd [(i + 63) % 4096];
					}
					if (SInd4 [(i + 4032) % 4096] == 4) {
						ChInd5 [SInd [i], 10] = SInd [(i + 4032) % 4096];
						ChInd5 [SInd [i], 6] = SInd [(i + 4033) % 4096];
					}
					if (SInd4 [(i + 4031)%4096] == 4) {
						ChInd5 [SInd [i], 9] = SInd [(i + 4031)%4096];
						ChInd5 [SInd [i], 5] = SInd [(i + 4031)%4096];
					}
				}
			}
			Qd = new int[NChannels];//noise amplitude
			Qm = new int[NChannels];//median
			Qdiff = new int[NChannels,2];
			Qmax = new int[NChannels];
			SqIv = new long[NChannels];//sum of squared channel increments
			SIprod = new long[NChannels,13];//sum of product of global and channel voltage increments
			SIp = new int[NChannels];
			Vbias = new int[NChannels];
			FVbias = new int[NChannels];
			FVsbias = new long[NChannels];
			Vsqbias = new long[NChannels];
			A = new int[NChannels];//control parameter for amplifier effects
			Qd4 = new int[NChannels];//noise amplitude
			Sl4 = new int[NChannels];//counter for spike length
			Sl4x = new bool[NChannels];
			AHP4 = new bool[NChannels];//counter for repolarizing current
			Amp4 = new int[NChannels];//buffers spike amplitude
			Qd5 = new int[NChannels];//noise amplitude
			Sl5 = new int[NChannels];//counter for spike length
			Sl5x = new bool[NChannels];
			AHP5 = new bool[NChannels];//counter for repolarizing current
			Amp5 = new int[NChannels];//buffers spike amplitude
			Slice = new int[NChannels];
			Slmax = sf / 835;
			CutOffset = sf / 835 + 6;
			CutAfter = sf / 1002;
			CutAfterLong = sf / 501;
			Sln0= sf / 3340;
			Sampling = sf;
			SqIglobal = 0;
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
			NumericUpDown numericUpDown1 = new NumericUpDown ();
			numericUpDown1.Parent = form1;
			numericUpDown1.Maximum = 20;
			numericUpDown1.Minimum = 6;
			numericUpDown1.Value = 10;
			numericUpDown1.Location = new Point (10, 70);
			Label l1 = new Label ();
			l1.Text = "Detection threshold (5 channel interpolation)";
			l1.AutoSize = true;
			l1.Location = new Point (numericUpDown1.Left, numericUpDown1.Bottom);
			l1.Size = numericUpDown1.Size;
			l1.Parent = form1;
			NumericUpDown numericUpDown3 = new NumericUpDown ();
			numericUpDown3.Parent = form1;
			numericUpDown3.Maximum = 20;
			numericUpDown3.Minimum = 6;
			numericUpDown3.Value = 9;
			numericUpDown3.Location = new Point (10, 120);
			Label l3 = new Label ();
			l3.Text = "Detection threshold (4 channel interpolation)";
			l3.AutoSize = true;
			l3.Location = new Point (numericUpDown3.Left, numericUpDown3.Bottom);
			l3.Size = numericUpDown3.Size;
			l3.Parent = form1;
			NumericUpDown numericUpDown2 = new NumericUpDown ();
			numericUpDown2.Parent = form1;
			numericUpDown2.Maximum = 5;
			numericUpDown2.Minimum = -5;
			numericUpDown2.Value = 0;
			numericUpDown2.Location = new Point (10, 170);
			Label l2 = new Label ();
			l2.Text = "Repolarization threshold";
			l2.AutoSize = true;
			l2.Location = new Point (numericUpDown2.Left, numericUpDown2.Bottom);
			l2.Size = numericUpDown2.Size;
			l2.Parent = form1;
			if (Indices [0] == 0) {
				DomainUpDown upDown1 = new DomainUpDown ();
				upDown1.Parent = form1;
				upDown1.Items.Add ("yes");
				upDown1.Items.Add ("no");
				upDown1.SelectedIndex = 0;
				upDown1.Location = new Point (10, 220);
				Label l0 = new Label ();
				l0.Text = "Recalibration";
				l0.AutoSize = true;
				l0.Location = new Point (upDown1.Left, upDown1.Bottom);
				l0.Size = upDown1.Size;
				l0.Parent = form1;
				form1.ShowDialog ();
				recalibTrigger = upDown1.SelectedIndex;
			}
			else {
				form1.ShowDialog ();
			}
			// set the parameters
			threshold5 = (int)numericUpDown1.Value;
			threshold4 = (int)numericUpDown3.Value;
			AHPthr = (int)numericUpDown2.Value;
			form1.Dispose();

						
			// open a box that logs acticity during detection
			logwindow = new Form();
			logwindow.Size = new Size(210, 180);
			logbox = new TextBox();
        	logbox.Parent = logwindow;
        	logbox.Dock = DockStyle.Fill;
        	logbox.Multiline = true;
			logwindow.Visible = true;
			logbox.Visible = true;
			logbox.AppendText("Detecting spikes...\n");
			logbox.Update();
			logwindow.Update();


			for (int i=0;i<NChannels;i++){
		   		Qd[i]=800;
				Qm[i]=0;
				Qdiff[i,0]=0;
				Qdiff[i,1]=0;
				Qmax[i]=0;
				SqIv[i]=0;//sum of squared channel increments
				SIp[i]=0;
				for (int ii=0; ii<13; ii++) {
					SIprod[i,ii]=0;//sum of product of global and channel voltage increments
				}
				FVbias[i]=0;
				Vbias[i]=0;
				FVsbias[i]=0;
				Vsqbias[i]=0;
				A[i]=0;
				Sl4[i]=0;
				Sl4x[i]=false;
				AHP4[i]=false;
				Amp4[i]=0;
				Sl5[i]=0;
				Sl5x[i]=false;
				AHP5[i]=false;
				Amp5[i]=0;
			}
			w.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wShapes.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wX.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wShapesX.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wInfo.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			wMean.BaseStream.Seek(0, SeekOrigin.Begin);   // Set the file pointer to the start.
			//write some info
			wInfo.WriteLine("# Number of frames:\n{0}", nFrames);
			wInfo.WriteLine("# Duration (s):\n{0}", nSec);
			wInfo.WriteLine("# Sampling rate:\n{0}", sf);
			wInfo.WriteLine("# Detection threshold (5,4):\n{0} {1}", threshold5, threshold4);
			wInfo.WriteLine("# Repolarization threshold:\n{0}", AHPthr);
			wInfo.WriteLine("# Recalibration trigger:\n{0}", recalibTrigger);
			wInfo.WriteLine ("# Recording channels:");
			for (int i=0; i<Indices.Length; i++) {
				wInfo.WriteLine ("{0}", Indices [i]);
			}
			wInfo.WriteLine ("# Recording channels4:");
			for (int i=0; i<NChannels; i++) {
				for (int j=0; j<12;j++){
					wInfo.Write("{0} ", ChInd4[i,j]);
				}
				wInfo.WriteLine ();
			}
			wInfo.WriteLine ("# Recording channels5:");
			for (int i=0; i<NChannels; i++) {
				for (int j=0; j<9;j++){
					wInfo.Write("{0} ", ChInd5[i,j]);
				}
				wInfo.WriteLine ();
			}
			wInfo.WriteLine("# Recalibration events:");
		}
		
		public void openSpikeFile (string name) {
			fs = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			w = new StreamWriter(fs);
		}
		
		public void openShapeFile (string name) {
			fsShapes = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wShapes = new StreamWriter(fsShapes);
		}
		
		public void openSpikeXFile (string name) {
			fsX = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wX = new StreamWriter(fsX);
		}
		
		public void openShapeXFile (string name) {
			fsShapesX = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wShapesX = new StreamWriter(fsShapesX);
		}

		public void openInfoFile (string name) {
			fsInfo = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wInfo = new StreamWriter(fsInfo);
		}

		public void openMeanFile (string name) {
			fsMean = new FileStream(name, FileMode.OpenOrCreate, FileAccess.Write);
			wMean = new StreamWriter(fsMean);
		}

		public void InitialEstimation (short[][] vm, long t0, long tInc) {//use this to get a better initial estimate of Qd
			if (t0 == 0) {
				for (int t=0; t<CutOffset; t++) {
					dt = (dt + 1) % 2;
					for (int i=1; i<NChannels; i++) {//loop across channels
						Slice [i] = ((vm [i] [t]) % 4095 + (vm [i] [t + 1]) % 4095 + (vm [i] [t + 2]) % 4095);
					}
					Array.Sort (Slice);
					Aglobal = Slice [NChannels / 2];
					for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
						//CHANNEL OUT OF LINEAR REGIME
						if (((vm [i] [t + 2] + 4) % 4096) < 10) {
							if (A [i] < artT) {//reset only when it starts leaving the linear regime
								A [i] = artT;
							}
						} else {
							Qm [i] = (2 * Qm [i] + (vm [i] [t] + vm [i] [t + 1] + vm [i] [t + 2] - Aglobal) * Ascale + 2 * Qd [i]) / 3;//update Qm
						}
					}
				}
			}
			for (int t=CutOffset; t<tInc-1; t++) {//loop over data, will be removed for an online algorithm
				dt = (dt + 1) % 2;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice [i] = ((vm [i] [t - 1]) % 4095 + (vm [i] [t]) % 4095 + (vm [i] [t + 1]) % 4095);
				}
				Array.Sort (Slice);
				Aglobal = Slice [NChannels / 2];
				// SPIKE DETECTION
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm [i] [t + 1] + 4) % 4096) < 10) {
						if (A [i] < artT) {//reset only when it starts leaving the linear regime
							A [i] = artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A [i] == 0) {
						Qdiff [i, dt] = (vm [i] [t - 1] + vm [i] [t] + vm [i] [t + 1] - Aglobal) * Ascale - Qm [i];//difference between ADC counts and Qm
						//UPDATE Qm and Qd
						if (Qdiff [i, dt] > 0) {
							if (Qdiff [i, dt] > Qd [i]) {
								Qm [i] += Qd [i] / Tau_m0;
								if (Qdiff [i, dt] < (5 * Qd [i])) {
									Qd [i]++;
								} else if ((Qd [i] > Qdmin) & (Qdiff [i, dt] > (6 * Qd [i]))) {
									Qd [i]--;
								}
							} else if (Qd [i] > Qdmin) {//set a minimum level for Qd
								Qd [i]--;
							}
						} else if (Qdiff [i, dt] < -Qd [i]) {
							Qm [i] -= Qd [i] / Tau_m0 / 2;

						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm [i] = (2 * Qm [i] + (vm [i] [t - 1] + vm [i] [t] + vm [i] [t + 1] - Aglobal) * Ascale + 2 * Qd [i]) / 3;//update Qm
						A [i]--;
					}
				}
			}
			if (t0 > 99 * tInc - CutOffset-1) {
				for (int i=0;i<NChannels;i++){
					Qm[i]=0;
					Qdiff[i,0]=0;
					Qdiff[i,1]=0;
					A[i]=0;
				}
			}
		}

		public void Iterate (short[][] vm, long t0, long tInc){
			int a4;//to buffer the difference between ADC counts and Qm
			int a5;//to buffer the difference between ADC counts and Qm
			int[] ChS4 = new int[4];
			int[] ChS5 = new int[4];
			int[] ChS4Arg= new int[4];
			int[] ChS5Arg= new int[4];
			if (t0==0){
				for (int t=0; t<CutOffset; t++) {
					dt=(dt+1)%2;
					for (int i=1; i<NChannels; i++) {//loop across channels
						Slice[i]=((vm[i][t])%4095+(vm[i][t+1])%4095+(vm[i][t+2])%4095);
					}
					Array.Sort(Slice);
					Aglobaldiff=Slice[NChannels / 2]-Aglobal;
					Aglobal=Slice[NChannels / 2];
					wMean.WriteLine("{0}", Aglobal);
					for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
						//CHANNEL OUT OF LINEAR REGIME
						if (((vm[i][t+2]+4)%4096)<10) {
							if (A[i]<artT) {//reset only when it starts leaving the linear regime
								A[i]=artT;
							}
						}
						else {
							Qm[i]=(2*Qm[i]+(vm[i][t]+vm[i][t+1]+vm[i][t+2]-Aglobal)*Ascale+2*Qd[i])/3;//update Qm
						}
					}
				}
			}
			for (int t=CutOffset; t<tInc-CutAfterLong; t++) {//loop over data, will be removed for an online algorithm
				dt=(dt+1)%2;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice[i]=((vm[i][t-1])%4095+(vm[i][t])%4095+(vm[i][t+1])%4095);
				}
				Array.Sort(Slice);
				Aglobaldiffold=Aglobaldiff+0;
				Aglobaldiff=Slice[NChannels / 2]-Aglobal;
				Aglobal=Slice[NChannels / 2];
				SqIglobal+=Aglobaldiff*Aglobaldiff;
				SqIg*=(Sampling-1);
				SqIg/=Sampling;
				SqIg+=Aglobaldiff;
				wMean.WriteLine("{0}", Aglobal);
				// RECALIBRATION EVENTS
				if (recalibTrigger==0){
					if (vm[0][t]<2500) {//newfiles:<2500_oldfiles:<1500
						if (Acal > 2000) {
							wInfo.Write("{0}", t+t0);//write time of recalibration event
							for (int i=0; i<NChannels; i++) {//loop across channels
								if (A[i]==0) {
									wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
								}
								else {
									wInfo.Write (" 0");
								}
							}
							wInfo.WriteLine();
							Console.WriteLine ("{0} sec", (t+t0)/Sampling);// to monitor progress of spike detection
							logbox.AppendText(String.Format("{0} sec\n", (t+t0)/Sampling));
							logwindow.Update();
							logbox.Update();
							Acal=0;//to remember last recalibration event
                    	}
					}
					Acal++;
				}
				else if ((t0+t)%Sampling==0) {//write spikes after every second
					wInfo.Write("{0}", t+t0);//write time of recalibration event
					for (int i=0; i<NChannels; i++) {//loop across channels
						if (A[i]==0) {
							wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
						}
						else {
							wInfo.Write (" 0");
						}
					}
					wInfo.WriteLine();
					Console.WriteLine ("{0} sec", (t+t0)/Sampling);// to monitor progress of spike detection
					logbox.AppendText(String.Format("{0} sec\n", (t+t0)/Sampling));
					logwindow.Update();
					logbox.Update();
				}
				// SPIKE DETECTION
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm[i][t+1]+4)%4096)<10) {
						if (A[i]<artT) {//reset only when it starts leaving the linear regime
							for (int ii=0; ii<Math.Min(artT-A[i],6); ii++) {//is only for one step in the past...ignoring others
								SIprod[i,12-ii]-=Aglobaldiffold*(vm[i][t+6-ii]-vm[i][t+3-ii]);
							}
							Qdiff[i,0]=0;
							Qdiff[i,1]=0;
							A[i]=artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A[i]==0) {
						SqIv[i]+=(vm[i][t+1]-vm[i][t-2])*(vm[i][t+1]-vm[i][t-2]);
						for (int iii=-6; iii<7; iii++){
							SIprod[i,iii+6]+=Aglobaldiff*(vm[i][t+1+iii]-vm[i][t-2+iii]);
						}
						//need a correlation variable which I can increment/decrement depending on whether 
						//the signal is correlated with the global signal.Significance? only update when high deviation?
						//Vbias=cumsum(fcorr*Aglobaldiff)/Sampling
						//fcorr+=1 if Qdiffdiff*Aglobaldiff>const
						//fcorr-=1 if Qdiffdiff*Aglobaldiff<const
						//might want to use 0 instead of const., as I want to converge...
						//unfortunately, update speed depends on Aglobaldiff, further, signal might be asymmetric
						//(update median increments only, i.e. if it correlates most of the time, increase corrcoeff)
						//--> need a threshold to get more reliable (e.g. Qd*std(Aglobaldiff) or Qd)
						Vbias[i]=FVbias[i]*SqIg/200;//local deviation of global signal;have to divide this by some time constant (e.g. 10)
						Qdiff[i,dt] = (vm[i][t-1]+vm[i][t]+vm[i][t+1]-Aglobal)*Ascale-Qm[i]-Vbias[i];//difference between ADC counts and Qm
						if (SIp[i]>0) {
							FVbias[i]++;
						}
						else {
							FVbias[i]--;
						}
						SIp[i]*=99;
						SIp[i]/=100;
						SIp[i]+=Aglobaldiff*(Qdiff[i,dt]-Qdiff[i,(dt+1)%2]);
						FVsbias[i]+=FVbias[i]/200;//want to see whether this matches the correlation structure
						Vsqbias[i]+=FVbias[i]/200*FVbias[i]/200;
						//UPDATE Qm and Qd
						if (Qdiff[i,dt]>0) {
							if (Qdiff[i,dt]>Qd[i]) {
								Qm[i]+=Qd[i]/Tau_m0;
								if  (Qdiff[i,dt]<(5*Qd[i])) {
									Qd[i]++;
								}
								else if ((Qd[i]>Qdmin) & (Qdiff[i,dt]>(6*Qd[i]))) {
									Qd[i]--;
								}
							}
							else if (Qd[i]>Qdmin){//set a minimum level for Qd
								Qd[i]--;
							}
						}
						else if (Qdiff[i,dt]<-Qd[i]){
							Qm[i]-=Qd[i]/Tau_m0/2;

						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm[i]=(2*Qm[i]+(vm[i][t-1]+vm[i][t]+vm[i][t+1]-Aglobal)*Ascale+2*Qd[i])/3;//update Qm
						A[i]--;
					}
					Qmax[i]=Math.Max(Qdiff[i,(dt+1)%2],Qdiff[i,dt]);
				}
				//4-channel interpolation
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					if (Sl4[i]>0) {//Sl frames after peak value
						for (int ii=0; ii<4; ii++) {
							if (A[ChInd4[i,ii]]==0){
								ChS4[ii]=Qmax[ChInd4[i,ii]];
								ChS4Arg[ii]=ii;
							}
							else {
								ChS4[ii]=0;
								ChS5Arg[ii]=ii;
							}
						}
						Array.Sort(ChS4,ChS4Arg);
						a4=0;
						Qd4[i]=0;
						for (int ii=0; ii<4; ii++) {
							a4+=ChS4[ii]*W4[ii];
							Qd4[i]+=Qd[ChInd4[i,ChS4Arg[ii]]]*W4[ii];
						}
						//TREATMENT OF THRESHOLD CROSSINGS
						if ((Sl4[i]==1) & ((Amp4[i])<((threshold4*Qd4[i])/2))){
							Sl4[i]=0;
							if (a4>((threshold4*Qd4[i])/2)) {
								Sl4[i]=1;
								Sl4x[i]=false;
								Amp4[i]=a4;
								AHP4[i]=false;
								for (int ii=0; ii<4; ii++) {
									Sx4[i,ii]=0;
								}
							}
						}
						else {
							//default
							Sl4[i]=(Sl4[i]+1)%(Slmax+1);// increment Sl[i]
							//check whether it repolarizes
							if (a4>(AHPthr*Qd4[i])) {
								AHP4[i]=true;
							}
							if (Sl4[i]==(Slmax-Sln0)){
								for (int ii=1; ii<5;ii++) {
									if (ChInd5[i,ii]>-1){
										if (ChInd5[i,ii]<ChInd4[i,0]){//have updated Sl4 already
											if (Sl4[ChInd5[i,ii]]>(Slmax-2*Sln0)){
												if ((Amp4[ChInd5[i,ii]]*Qd4[i])<(Amp4[i]*Qd4[ChInd5[i,ii]])){
													Sl4x[ChInd5[i,ii]]=true;
													Sx4[i,ii-1]=1;
													Sx4[i,ii%4]=1;
												}
											}
										}
										else {
											if ((Sl4[ChInd5[i,ii]]>(Slmax-2*Sln0-1)) & (Sl4[ChInd5[i,ii]]<(Slmax-1))){
												if ((Amp4[ChInd5[i,ii]]*Qd4[i])<(Amp4[i]*Qd4[ChInd5[i,ii]])){
													Sl4x[ChInd5[i,ii]]=true;
													Sx4[i,ii-1]=1;
													Sx4[i,ii%4]=1;
												}
											}
										}
									}
								}
								for (int ii=12; ii<16;ii++) {
									if (ChInd4[i,ii]>-1){
										if ((Sl5[ChInd4[i,ii]]>(Slmax-2*Sln0-1)) & (Sl5[ChInd4[i,ii]]<(Slmax-1))){
											if ((Amp5[ChInd4[i,ii]]*Qd4[i])<(Amp4[i]*Qd5[ChInd4[i,ii]])){
												Sl5x[ChInd4[i,ii]]=true;
												Sx4[i,ii-12]=1;
												//Sx4[i,(ii-11)%4]=1;
											}
										}
									}
								}
							}
							//accept spikes after Slmax frames if...
							if ((Sl4[i]==Slmax) & (!Sl4x[i])) {// & (AHP4[i]<(Slmax-Slmin))
								if (AHP4[i]) {
									wX.Write("{0} {1} {2} {3} {4} {5}", i, t0+t-Slmax+1, Amp4[i], Qd4[i], 1, Acal-Slmax+1);
									for (int ii=0; ii<4; ii++) {
										wX.Write(" {0}", Sx4[i,ii]);
									}
									wX.WriteLine ();
									for (int ii=0; ii<4;ii++){
										wShapesX.Write("{0} {1} ", ChInd4[i,ii], (Qm[ChInd4[i,ii]]+Vbias[ChInd4[i,ii]])/3);
										for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
											wShapesX.Write("{0} ", vm[ChInd4[i,ii]][t-CutOffset+jj]);
										}
									}
									for (int ii=4; ii<12;ii++){
										wShapesX.Write("{0} ", ChInd4[i,ii]);
										if (ChInd4[i,ii]>-1){
											wShapesX.Write("{0} ", (Qm[ChInd4[i,ii]]+Vbias[ChInd4[i,ii]])/3);
											for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
												wShapesX.Write("{0} ", vm[ChInd4[i,ii]][t-CutOffset+jj]);
											}
										}
										else {
											for (int jj=0; jj<(CutOffset+CutAfter+2); jj++){
												wShapesX.Write("0 ");
											}
										}
									}
									for (int ii=0; ii<4;ii++){
										if (Sx4[i,ii]==1) {
											if (ChInd5[ChInd4[i,ii],ii+5]>-1){
												wShapesX.Write("{0} {1} ", ChInd5[ChInd4[i,ii],ii+5], (Qm[ChInd5[ChInd4[i,ii],ii+5]]+Vbias[ChInd5[ChInd4[i,ii],ii+5]])/3);
												for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
													wShapesX.Write("{0} ", vm[ChInd5[ChInd4[i,ii],ii+5]][t-CutOffset+jj]);
												}
											}
										}
									}
								}
								else {
									wX.Write("{0} {1} {2} {3} {4} {5}", i, t0+t-Slmax+1, Amp4[i], Qd4[i], 0, Acal-Slmax+1);
									for (int ii=0; ii<4; ii++) {
										wX.Write(" {0}", Sx4[i,ii]);
									}
									wX.WriteLine ();
									for (int ii=0; ii<4;ii++){
										wShapesX.Write("{0} {1} ", ChInd4[i,ii], (Qm[ChInd4[i,ii]]+Vbias[ChInd4[i,ii]])/3);
										for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
										wShapesX.Write("{0} ", vm[ChInd4[i,ii]][t-CutOffset+jj]);
										}
									}
									for (int ii=4; ii<12;ii++){
										wShapesX.Write("{0} ", ChInd4[i,ii]);
										if (ChInd4[i,ii]>-1){
											wShapesX.Write("{0} ", (Qm[ChInd4[i,ii]]+Vbias[ChInd4[i,ii]])/3);
											for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
												wShapesX.Write("{0} ", vm[ChInd4[i,ii]][t-CutOffset+jj]);
											}
										}
										else {
											for (int jj=0; jj<(CutOffset+CutAfterLong+2); jj++){
												wShapesX.Write("0 ");
											}
										}
									}
									for (int ii=0; ii<4;ii++){
										if (Sx4[i,ii]==1) {
											if (ChInd5[ChInd4[i,ii],ii+5]>-1){
												wShapesX.Write("{0} {1} ", ChInd5[ChInd4[i,ii],ii+5], (Qm[ChInd5[ChInd4[i,ii],ii+5]]+Vbias[ChInd5[ChInd4[i,ii],ii+5]])/3);
												for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
													wShapesX.Write("{0} ", vm[ChInd5[ChInd4[i,ii],ii+5]][t-CutOffset+jj]);
												}
											}
										}
									}
								}
								wShapesX.WriteLine();
								Sl4[i]=0;
							}
							//check whether current ADC count is higher
							else if (Amp4[i]<a4) {
								Sl4[i]=1;//reset peak value
								Sl4x[i]=false;
								Amp4[i]=a4;
								AHP4[i]=false;//reset AHP
							}
						}
					}
					else if (ChInd4[i,0]>-1){
						if (Math.Max(Math.Min(Qmax[ChInd4[i,0]],Qmax[ChInd4[i,3]]),Math.Min(Qmax[ChInd4[i,1]],Qmax[ChInd4[i,2]]))>(2000)){
							for (int ii=0; ii<4; ii++) {
								if (A[ChInd4[i,ii]]==0){
									ChS4[ii]=Qmax[ChInd4[i,ii]];
									ChS4Arg[ii]=ii;
								}
								else {
									ChS4[ii]=0;
									ChS5Arg[ii]=ii;
								}
							}
							Array.Sort(ChS4,ChS4Arg);
							a4=0;
							Qd4[i]=0;
							for (int ii=0; ii<4; ii++) {
								a4+=ChS4[ii]*W4[ii];
								Qd4[i]+=Qd[ChInd4[i,ChS4Arg[ii]]]*W4[ii];
							}
							//check for threshold crossings
							if (a4>((threshold4*Qd4[i])/2)) {
								Sl4[i]=1;
								Sl4x[i]=false;
								Amp4[i]=a4;
								AHP4[i]=false;
								for (int ii=0; ii<4; ii++) {
									Sx4[i,ii]=0;
								}
							}
						}
					}
				}
				//5 channel interpolation
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					if (Sl5[i]>0) {//Sl frames after peak value
						for (int ii=1; ii<5; ii++) {
							if (A[ChInd5[i,ii]]==0){
								ChS5[ii-1]=Qmax[ChInd5[i,ii]];
								ChS5Arg[ii-1]=ii;
							}
							else {
								ChS5[ii-1]=0;
								ChS5Arg[ii-1]=ii;
							}
						}
						Array.Sort(ChS5,ChS5Arg);
						a5=3*Qmax[ChInd5[i,0]];
						Qd5[i]=3*Qd[ChInd5[i,0]];
						for (int ii=0; ii<4; ii++) {
							a5+=ChS5[ii]*W5[ii];
							Qd5[i]+=Qd[ChInd5[i,ChS5Arg[ii]]]*W5[ii];
						}
						//TREATMENT OF THRESHOLD CROSSINGS
						if (Sl5[i]>0) {//Sl frames after peak value
							if ((Sl5[i]==1) & ((Amp5[i])<((threshold5*Qd5[i])/2))){
								Sl5[i]=0;
								if (a5>((threshold5*Qd5[i])/2)) {
									Sl5[i]=1;
									Sl5x[i]=false;
									Amp5[i]=a5;
									AHP5[i]=false;
									for (int ii=0; ii<4; ii++) {
										Sx5[i,ii]=0;
									}
								}
							}
							else {
								//default
								Sl5[i]=(Sl5[i]+1)%(Slmax+1);// increment Sl[i]
								//check whether it doesn't repolarize
								if (a5<(AHPthr*Qd5[i])) {
									AHP5[i]=true;
								}
								if ((Sl5[i]==(Slmax-Sln0))){
									for (int ii=1; ii<5;ii++) {
										if (ChInd5[i,ii]<ChInd5[i,0]){//have updated Sl5 already
											if (Sl5[ChInd5[i,ii]]>(Slmax-2*Sln0)){
												if ((Amp5[ChInd5[i,ii]]*Qd5[i])<(Amp5[i]*Qd5[ChInd5[i,ii]])){
													Sl5x[ChInd5[i,ii]]=true;
													Sx5[i,ii-1]=1;
													Sx5[i,(ii)%4]=1;
												}
											}
										}
										else {
											if ((Sl5[ChInd5[i,ii]]>(Slmax-2*Sln0-1)) & (Sl5[ChInd5[i,ii]]<(Slmax-1))){
												if ((Amp5[ChInd5[i,ii]]*Qd5[i])<(Amp5[i]*Qd5[ChInd5[i,ii]])){
													Sl5x[ChInd5[i,ii]]=true;
													Sx5[i,ii-1]=1;
													Sx5[i,(ii)%4]=1;
												}
											}
										}
									}
									for (int ii=9; ii<13;ii++) {
										if (ChInd5[i,ii]>-1) {
											if (Sl4[ChInd5[i,ii]]>(Slmax-2*Sln0)){//have updated Sl4 already
												if ((Amp4[ChInd5[i,ii]]*Qd5[i])<(Amp5[i]*Qd4[ChInd5[i,ii]])){
													Sl4x[ChInd5[i,ii]]=true;
													Sx5[i,ii-9]=1;
												}
											}
										}
									}
								}
								//accept spikes after Slmax frames if...
								if ((Sl5[i]==Slmax) & (!Sl5x[i])) {
									if (AHP5[i]) {
										w.Write("{0} {1} {2} {3} {4} {5}", i, t0+t-Slmax+1, Amp5[i], Qd5[i], 1, Acal-Slmax+1);
										for (int ii=0; ii<4; ii++) {//to be removed...?
											w.Write(" {0}", (Sx5[i,ii]+Sx5[i,(ii+1)%4])/2);
										}
										w.WriteLine ();
										for (int ii=0; ii<5;ii++){
											wShapes.Write("{0} {1} ", ChInd5[i,ii], (Qm[ChInd5[i,ii]]+Vbias[ChInd5[i,ii]])/3);
											for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
											wShapes.Write("{0} ", vm[ChInd5[i,ii]][t-CutOffset+jj]);
											}
										}
										for (int ii=0; ii<4;ii++){
											wShapes.Write("{0} ", ChInd5[i,ii+5]);
											if (ChInd5[i,ii+5]>-1){
												wShapes.Write("{0} ", (Qm[ChInd5[i,ii+5]]+Vbias[ChInd5[i,ii+5]])/3);
												for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
													wShapes.Write("{0} ", vm[ChInd5[i,ii+5]][t-CutOffset+jj]);
												}
											}
											else {
												for (int jj=0; jj<(CutOffset+CutAfter+2); jj++){
													wShapes.Write("0 ");
												}
											}
										}
										for (int ii=0; ii<4;ii++){
											if ((Sx5[i,ii]+Sx5[i,(ii+1)%4])/2==1) {
												if (ChInd5[ChInd5[i,ii+1],ii+5]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],ii+5]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],ii+5]]+Vbias[ChInd5[ChInd5[i,ii+1],ii+5]])/3);
													for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],ii+5]][t-CutOffset+jj]);
													}
												}
												if (ChInd5[ChInd5[i,ii+1],ii+1]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],ii+1]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],ii+1]]+Vbias[ChInd5[ChInd5[i,ii+1],ii+1]])/3);
													for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],ii+1]][t-CutOffset+jj]);
													}
												}
												if (ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]]+Vbias[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]])/3);
													for (int jj=0; jj<(CutOffset+CutAfter+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]][t-CutOffset+jj]);
													}
												}
											}
										}
									}
									else {
										w.Write("{0} {1} {2} {3} {4} {5}", i, t0+t-Slmax+1, Amp5[i], Qd5[i], 0, Acal-Slmax+1);
										for (int ii=0; ii<4; ii++) {
											w.Write(" {0}", (Sx5[i,ii]+Sx5[i,(ii+1)%4])/2);
										}
										w.WriteLine ();
										for (int ii=0; ii<5;ii++){
											wShapes.Write("{0} {1} ", ChInd5[i,ii], (Qm[ChInd5[i,ii]]+Vbias[ChInd5[i,ii]])/3);
											for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
											wShapes.Write("{0} ", vm[ChInd5[i,ii]][t-CutOffset+jj]);
											}
										}
										for (int ii=0; ii<4;ii++){
											wShapes.Write("{0} ", ChInd5[i,ii+5]);
											if (ChInd5[i,ii+5]>-1){
												wShapes.Write("{0} ", (Qm[ChInd5[i,ii+5]]+Vbias[ChInd5[i,ii+5]])/3);
												for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
													wShapes.Write("{0} ", vm[ChInd5[i,ii+5]][t-CutOffset+jj]);
												}
											}
											else {
												for (int jj=0; jj<(CutOffset+CutAfterLong+2); jj++){
													wShapes.Write("0 ");
												}
											}
										}
										for (int ii=0; ii<4;ii++){
											if ((Sx5[i,ii]+Sx5[i,(ii+1)%4])/2==1) {
												if (ChInd5[ChInd5[i,ii+1],ii+5]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],ii+5]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],ii+5]]+Vbias[ChInd5[ChInd5[i,ii+1],ii+5]])/3);
													for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],ii+5]][t-CutOffset+jj]);
													}
												}
												if (ChInd5[ChInd5[i,ii+1],ii+1]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],ii+1]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],ii+1]]+Vbias[ChInd5[ChInd5[i,ii+1],ii+1]])/3);
													for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],ii+1]][t-CutOffset+jj]);
													}
												}
												if (ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]>-1){
													wShapes.Write("{0} ", ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]);
													wShapes.Write("{0} ", (Qm[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]]+Vbias[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]])/3);
													for (int jj=0; jj<(CutOffset+CutAfterLong+1); jj++){
														wShapes.Write("{0} ", vm[ChInd5[ChInd5[i,ii+1],(ii+1)%4+5]][t-CutOffset+jj]);
													}
												}
											}
										}
									}
									wShapes.WriteLine();
									Sl5[i]=0;
								}
								//check whether current ADC count is higher
								else if (Amp5[i]<a5) {
									Sl5[i]=1;//reset peak value
									Sl5x[i]=false;
									Amp5[i]=a5;
									AHP5[i]=false;//reset AHP
								}
							}
						}
					}
					else if (ChInd5[i,0]>-1){
						if (Qmax[ChInd5[i,0]]>(4*Qd[i])){
							for (int ii=1; ii<5; ii++) {
								if (A[ChInd5[i,ii]]==0){
									ChS5[ii-1]=Qmax[ChInd5[i,ii]];
									ChS5Arg[ii-1]=ii;
								}
								else {
									ChS5[ii-1]=0;
									ChS5Arg[ii-1]=ii;
								}
							}
							Array.Sort(ChS5,ChS5Arg);
							a5=3*Qmax[ChInd5[i,0]];
							Qd5[i]=3*Qd[ChInd5[i,0]];
							for (int ii=0; ii<4; ii++) {
								a5+=ChS5[ii]*W5[ii];
								Qd5[i]+=Qd[ChInd5[i,ChS5Arg[ii]]]*W5[ii];
							}
							//check for threshold crossings
							if (a5>((threshold5*Qd5[i])/2)) {
								Sl5[i]=1;
								Sl5x[i]=false;
								Amp5[i]=a5;
								AHP5[i]=false;
								for (int ii=0; ii<4; ii++) {
									Sx5[i,ii]=0;
								}
							}
						}
					}
				}
			}
		}

		public void IterateLast (short[][] vm, long t0){
			//continue baseline estimation
			for (int t=CutOffset; t<CutOffset+CutAfterLong-1; t++) {//loop over data, will be removed for an online algorithm
				dt=(dt+1)%2;
				for (int i=1; i<NChannels; i++) {//loop across channels
					Slice[i]=((vm[i][t-1])%4095+(vm[i][t])%4095+(vm[i][t+1])%4095);
				}
				Array.Sort(Slice);
				Aglobaldiffold=Aglobaldiff+0;
				Aglobaldiff=Slice[NChannels / 2]-Aglobal;
				Aglobal=Slice[NChannels / 2];
				SqIglobal+=Aglobaldiff*Aglobaldiff;
				SqIg*=(Sampling-1);
				SqIg/=Sampling;
				SqIg+=Aglobaldiff;
				wMean.WriteLine("{0}", Aglobal);
				// RECALIBRATION EVENTS
				if (recalibTrigger==0){
					if (vm[0][t]<2500) {//newfiles:<2500_oldfiles:<1500
						if (Acal > 2000) {
							wInfo.Write("{0}", t+t0);//write time of recalibration event
							for (int i=0; i<NChannels; i++) {//loop across channels
								if (A[i]==0) {
									wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
								}
								else {
									wInfo.Write (" 0");
								}
							}
							wInfo.WriteLine();
							Console.WriteLine ("{0} sec", (t+t0)/Sampling);// to monitor progress of spike detection
							logbox.AppendText(String.Format("{0} sec\n", (t+t0)/Sampling));
							logwindow.Update();
							logbox.Update();
							Acal=0;//to remember last recalibration event
						}
					}
					Acal++;
				}
				else if ((t0+t)%Sampling==0) {
					wInfo.Write("{0}", t+t0);//write time of recalibration event
					for (int i=0; i<NChannels; i++) {//loop across channels
						if (A[i]==0) {
							wInfo.Write(" {0}", Qd[i]);//write variance of recalibration event
						}
						else {
							wInfo.Write (" 0");
						}
					}
					wInfo.WriteLine();
					Console.WriteLine ("{0} sec", (t+t0)/Sampling);// to monitor progress of spike detection
					logbox.AppendText(String.Format("{0} sec\n", (t+t0)/Sampling));
					logwindow.Update();
					logbox.Update();
				}
				// SPIKE DETECTION
				for (int i=1-recalibTrigger; i<NChannels; i++) {//loop across channels
					//CHANNEL OUT OF LINEAR REGIME
					if (((vm[i][t+1]+4)%4096)<10) {
						if (A[i]<artT) {//reset only when it starts leaving the linear regime
							//for (int ii=0; ii<Math.Min(artT-A[i],6); ii++) {//is only for one step in the past...ignoring others
							//	SIprod[i,12-ii]-=Aglobaldiffold*(vm[i][t+6-ii]-vm[i][t+3-ii]);
							//}
							Qdiff[i,0]=0;
							Qdiff[i,1]=0;
							A[i]=artT;
						}
					}
					//DEFAULT OPERATIONS
					else if (A[i]==0) {
						SqIv[i]+=(vm[i][t+1]-vm[i][t-2])*(vm[i][t+1]-vm[i][t-2]);
						//for (int iii=-6; iii<7; iii++){
						//	SIprod[i,iii+6]+=Aglobaldiff*(vm[i][t+1+iii]-vm[i][t-2+iii]);
						//}
						//need a correlation variable which I can increment/decrement depending on whether 
						//the signal is correlated with the global signal.Significance? only update when high deviation?
						//Vbias=cumsum(fcorr*Aglobaldiff)/Sampling
						//fcorr+=1 if Qdiffdiff*Aglobaldiff>const
						//fcorr-=1 if Qdiffdiff*Aglobaldiff<const
						//might want to use 0 instead of const., as I want to converge...
						//unfortunately, update speed depends on Aglobaldiff, further, signal might be asymmetric
						//(update median increments only, i.e. if it correlates most of the time, increase corrcoeff)
						//--> need a threshold to get more reliable (e.g. Qd*std(Aglobaldiff) or Qd)
						Vbias[i]=FVbias[i]*SqIg/200;//local deviation of global signal;have to divide this by some time constant (e.g. 10)
						Qdiff[i,dt] = (vm[i][t-1]+vm[i][t]+vm[i][t+1]-Aglobal)*Ascale-Qm[i]-Vbias[i];//difference between ADC counts and Qm
						if (SIp[i]>0) {
							FVbias[i]++;
						}
						else {
							FVbias[i]--;
						}
						SIp[i]*=99;
						SIp[i]/=100;
						SIp[i]+=Aglobaldiff*(Qdiff[i,dt]-Qdiff[i,(dt+1)%2]);
						FVsbias[i]+=FVbias[i]/200;//want to see whether this matches the correlation structure
						Vsqbias[i]+=FVbias[i]/200*FVbias[i]/200;
						//UPDATE Qm and Qd
						if (Qdiff[i,dt]>0) {
							if (Qdiff[i,dt]>Qd[i]) {
								Qm[i]+=Qd[i]/Tau_m0;
								if  (Qdiff[i,dt]<(5*Qd[i])) {
									Qd[i]++;
								}
								else if ((Qd[i]>Qdmin) & (Qdiff[i,dt]>(6*Qd[i]))) {
									Qd[i]--;
								}
							}
							else if (Qd[i]>Qdmin){//set a minimum level for Qd
								Qd[i]--;
							}
						}
						else if (Qdiff[i,dt]<-Qd[i]){
							Qm[i]-=Qd[i]/Tau_m0/2;

						}
					}
					//AFTER CHANNEL WAS OUT OF LINEAR REGIME
					else {
						Qm[i]=(2*Qm[i]+(vm[i][t-1]+vm[i][t]+vm[i][t+1]-Aglobal)*Ascale+2*Qd[i])/3;//update Qm
						A[i]--;
					}
					Qmax[i]=Math.Max(Qdiff[i,(dt+1)%2],Qdiff[i,dt]);
				}
			}
			wMean.WriteLine("{0}", Aglobal);
		}


		public void FinishDetection ()
		{//write spikes in interval after last recalibration; close file
			wInfo.WriteLine ("#Sum(squared global fluctuations):");
			wInfo.WriteLine ("{0}", SqIglobal);
			wInfo.WriteLine ("#Sum(squared channel fluctuations):");
			for (int i=0; i<NChannels; i++) {//loop across channels
				wInfo.Write ("{0} ", SqIv[i]);
			}
			wInfo.WriteLine();
			wInfo.WriteLine ("#Sum(product of channel and global fluctuations):");
			for (int ii=0; ii<13; ii++) {//loop across timelags
				for (int i=0; i<NChannels; i++) {//loop across channels
					wInfo.Write ("{0} ", SIprod[i,ii]);
				}
				wInfo.WriteLine();
			}
			wInfo.WriteLine ("#Sum(avg. deviations from global fluctuations):");
			for (int i=0; i<NChannels; i++) {//loop across channels
				wInfo.Write ("{0} ", FVsbias[i]);
			}
			wInfo.WriteLine();
			wInfo.WriteLine ("#Sum(avg. squared deviations from global fluctuations):");
			for (int i=0; i<NChannels; i++) {//loop across channels
				wInfo.Write ("{0} ", Vsqbias[i]);
			}
			wInfo.WriteLine();
			w.Flush();
			w.Close();
			fs.Close();
			wShapes.Flush();
			wShapes.Close();
			fsShapes.Close();
			wX.Flush();
			wX.Close();
			fsX.Close();
			wShapesX.Flush();
			wShapesX.Close();
			fsShapesX.Close();
			wInfo.Flush();
			wInfo.Close();
			fsInfo.Close();
			wMean.Flush();
			wMean.Close();
			fsMean.Close();
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
			//openFileDialog1.InitialDirectory = "c:\\";
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
				int sf = brwRdr.SamplingRate;
				Console.WriteLine ("# Number of frames: {0}", nFrames);
				Console.WriteLine ("# Duration (s): {0}", nSec);
				Console.WriteLine ("# Sampling rate: {0}", sf);

				SaveFileDialog saveFileDialog1 = new SaveFileDialog ();
				saveFileDialog1.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
				saveFileDialog1.FilterIndex = 1;
				saveFileDialog1.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Spikes";
				saveFileDialog1.DefaultExt = "txt";
				saveFileDialog1.InitialDirectory = Directory.GetCurrentDirectory ();
				saveFileDialog1.Title = "Save Spikes As";
				
				if (saveFileDialog1.ShowDialog () == DialogResult.OK) {
					SaveFileDialog saveFileDialog2 = new SaveFileDialog ();
					saveFileDialog2.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
					saveFileDialog2.FilterIndex = 1;
					saveFileDialog2.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Shapes";
					saveFileDialog2.DefaultExt = "txt";
					saveFileDialog2.Title = "Save Shapes As";

					if (saveFileDialog2.ShowDialog () == DialogResult.OK) {
						SaveFileDialog saveFileDialog5 = new SaveFileDialog ();
						saveFileDialog5.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
						saveFileDialog5.FilterIndex = 1;
						saveFileDialog5.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_SpikesX";
						saveFileDialog5.DefaultExt = "txt";
						saveFileDialog5.Title = "Save Interpolated Shapes As";

						if (saveFileDialog5.ShowDialog () == DialogResult.OK) {
							SaveFileDialog saveFileDialog6 = new SaveFileDialog ();
							saveFileDialog6.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
							saveFileDialog6.FilterIndex = 1;
							saveFileDialog6.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_ShapesX";
							saveFileDialog6.DefaultExt = "txt";
							saveFileDialog6.Title = "Save Interpolated Shapes As";
					
							if (saveFileDialog6.ShowDialog () == DialogResult.OK) {
								SaveFileDialog saveFileDialog3 = new SaveFileDialog ();
								saveFileDialog3.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
								saveFileDialog3.FilterIndex = 1;
								saveFileDialog3.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Info";
								saveFileDialog3.DefaultExt = "txt";
								saveFileDialog3.Title = "Save Further Information As";

								if (saveFileDialog3.ShowDialog () == DialogResult.OK) {
									SaveFileDialog saveFileDialog4 = new SaveFileDialog ();
									saveFileDialog4.Filter = "txt files (*.txt)|*.txt|All files (*.*)|*.*";
									saveFileDialog4.FilterIndex = 1;
									saveFileDialog4.FileName = Path.GetFileNameWithoutExtension (openFileDialog1.FileName) + "_Avg";
									saveFileDialog4.DefaultExt = "txt";
									saveFileDialog4.Title = "Save Average Voltage Increments As";
								
									if (saveFileDialog4.ShowDialog () == DialogResult.OK) {
										//list of channels
										ChCoord[] Channels = brwRdr.GetRecChs (StreamType.Raw);
										Console.WriteLine ("# Number of recorded channels: {0}", Channels.Length);
										int[] Indices = new int[Channels.Length];
										for (int i=0; i<Channels.Length; i++) {
											Indices [i] = (Channels [i].Col - 1) + 64 * (Channels [i].Row - 1);
											//Console.WriteLine ("{0}", Indices [i]);
										}
										// how many frames to analyze
										long nDumpFrames = nFrames;
										Detection SpkD;
										SpkD = new Detection ();
			
										// open the output file
										SpkD.openSpikeFile (saveFileDialog1.FileName);
										saveFileDialog1.Dispose ();
										SpkD.openShapeFile (saveFileDialog2.FileName);
										saveFileDialog2.Dispose ();
										SpkD.openSpikeXFile (saveFileDialog5.FileName);
										saveFileDialog5.Dispose ();
										SpkD.openShapeXFile (saveFileDialog6.FileName);
										saveFileDialog6.Dispose ();
										SpkD.openInfoFile (saveFileDialog3.FileName);
										saveFileDialog3.Dispose ();
										SpkD.openMeanFile (saveFileDialog4.FileName);
										saveFileDialog4.Dispose ();
										SpkD.SetInitialParams (nFrames, nSec, sf, Channels.Length, Indices);
								
										// measure execution time
										var sw = new Stopwatch ();
										sw.Start ();
									
										const long tInc = 128;
										int tCut = sf/501 + sf/835 + 6;
										int CutOffset = sf / 835 + 6;
										const int NChannels = 4096;
										for (long t0=0; t0<99*tInc; t0+=tInc - CutOffset-1) {
											short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
											SpkD.InitialEstimation (vm, t0, tInc);
										}
										for (long t0=0; t0<nDumpFrames-tInc; t0+=tInc-tCut) {
											short[] [] vm = brwRdr.GetRawDataADCCounts (Channels, t0, tInc);
											SpkD.Iterate (vm, t0, tInc);
										}
										long t0a=nDumpFrames-tCut-((nDumpFrames-tInc)%(tInc-tCut));
										if (t0a < (nDumpFrames - tCut)) {
											short[] [] vma = brwRdr.GetRawDataADCCounts (Channels, t0a, nDumpFrames - t0a);
											SpkD.Iterate (vma, t0a, nDumpFrames - t0a);
										}
										short[] [] vmf = brwRdr.GetRawDataADCCounts (Channels, nDumpFrames-tCut, tCut);
										SpkD.IterateLast (vmf, nDumpFrames-tCut);
										SpkD.FinishDetection ();
										sw.Stop ();
										Console.WriteLine ("Elapsed time: {0}", sw.Elapsed); // TimeSpan
										Console.WriteLine ("Milliseconds/frame: {0}", sw.Elapsed.TotalMilliseconds / nFrames); // TimeSpan
									}
								}
							}
						}
					}
				}
			}
		}
	}
}




