#include "stdlib.h"
#include <iostream>
#include "string.h"
#include <vector>
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"

double waveform_baseline(const std::vector<double>* waveform, const int N) {
    /* Find the baseline of waveform from the first N samples; 
     * @parameters: waveform as a vector and N as an integer.
     * @return:     baseline as a double.
     */
    double sample_sum = 0.;
    int count = 0;
    for (int k=1; k <20; k++) {
        sample_sum += waveform->at(k);
        count++;
    }
    double baseline = sample_sum/float(count);
    return baseline;
}

double pulse_charge(const std::vector<double>* waveform, const double baseline, const int center, const int width) {
    /* Find charge of pulse in "width" before and after "center"
     * @parameters: waveform, baseline, center of integration region, width of region.
     * @return:     value of integram in the region  
     */
    double charge = 0.;
    int kmin = center - width;
    int kmax = center + width; 
    if ((center - width) < 0) {kmin = 0;}
    if ((center + width)>waveform->size()) {kmax = waveform->size();}
    for (int k = kmin; k < kmax; k++)  charge += waveform->at(k);
    return charge;
}

double CFD_timing(const std::vector<double>* waveform, const double baseline, const int global_imin, const float startp, const float endp, const float percentage) {
    /* Find timing of pulse with CFD technique.
     * @parameters: waveform, baseline, global minimum, range start point, range end point, percentage of rise time.
     * @return:     value of CFD time
     */
    double y_min = waveform->at(global_imin)-baseline;
    double y_end = endp*y_min;
    double y_start = startp*y_min;
    double rise_amplitude = (endp-startp)*y_min;
    int j = global_imin;
    int j_end = 0;
    while (((waveform->at(j)-baseline) < y_end) && (j > 1)) {
        j--;
        j_end = j-1;
    }
    j = global_imin;
    int j_start = 0;
    while (((waveform->at(j)-baseline) < y_start) && (j > 1)) {
        j--;
        j_start = j+1;
    }
    double b = (y_end-y_start)/(double(j_end)-double(j_start));
    double a = y_start - j_start*b;
    double iCFD = (percentage*rise_amplitude-a)/b;
    return iCFD;
}

int global_timing(const std::vector<double>* waveform, const double baseline, const int start, const int end) {
    /* Find the timing of the waveform from STT technique.
     * @parameters: waveform vector, baseline value.
     * @return:     integer with sample of global minimum.
     */
    double ymin = 99999999.;
    int imin = 0;
    int nsamples = waveform->size();
    for (int k = start; k < end; k++){
        if (waveform->at(k)-baseline < ymin) {
            ymin = waveform->at(k)-baseline;
            imin = k;
        }
    }
    return imin;
}

std::vector<int> count_pulses(const std::vector<double>* waveform, const double baseline, const double threshold_ADC) {
    /* Count number of times waveform goes below a threshold.
     * @parameters: waveform samples, baseline threshold level.
     * @return: by reference, number and position of pulses.
    */
    double V_per_unit = 1./16384.*2.*1.e3;
    double threshold = threshold_ADC/V_per_unit;
    std::vector<int> pulses;
    
    double ymin = 9999999.;
    int imin = 0;
    bool inside_pulse = false;
    for (int k=0; k<waveform->size(); k++) {

        if (waveform->at(k)-baseline < ymin) {
            ymin = waveform->at(k)-baseline;
            imin = k;
        }

        if (waveform->at(k) < threshold && !inside_pulse) {
            inside_pulse = true;
        }

        if (waveform->at(k) >= threshold && inside_pulse) {
            inside_pulse = false;
            pulses.push_back(imin);
            imin=0; ymin=9999999;
        }
    }

    return pulses;
}

void display_waveform(const std::vector<double>* waveform, const double baseline, const int STT_time, const double CFD_time, std::vector<int> pulses, const int event, const int channel) {
    
    /* Plot waveform into a png file in 'waveforms' folder. 
     * @parameters: waveform, baseline, times and channel/event numbers.
     * @return:     void (save image only)
    */

    // 14-bit @ 500 MS/s
    double ns_per_bin = 1./500.e6;
    double V_per_unit = 1./16384.*2.*1.e3;

    // Plot base waveform
    TCanvas* waveform_canvas = new TCanvas(Form("waveformEv%dCh%d"), Form("Ev%dCh%d", event, channel), 1200, 600);
    waveform_canvas->SetTitle(Form("Event %d Channel %d", event, channel));
    TH1D* waveform_histogram = new TH1D(
                             Form("waveformEv%dCh%d"), 
                             Form("Waveform Ev%dCh%d"), 
                             waveform->size(), 
                             0, 
                             waveform->size()-1
                            );
    for (int i=0; i < waveform->size()-1; i++) {waveform_histogram->SetBinContent(i,waveform->at(i)*V_per_unit);}
    waveform_canvas->cd(0);
    waveform_histogram->Draw("HIST");

    // Plot all pulses
    TGraph* time_graph = new TGraph();
    if (pulses.size() > 0) {
        time_graph->Set(pulses.size());
        time_graph->SetMarkerStyle(8);
        time_graph->SetMarkerColor(2);
        for (int i=0; i < pulses.size(); i++) {
            time_graph->SetPoint(i+1, pulses.at(i), waveform->at(pulses.at(i))*V_per_unit);
        }
        time_graph->Draw("SAME P");
    } 

    // Save file
    if (pulses.size() > 0) {
        waveform_canvas->SaveAs(Form("../waveforms/waveform_%d_%d.png", event, channel));    
    }
    return;
}

void analyze_data_230719(){

    // 14-bit @ 500 MS/s
    double ns_per_bin = 1./500.e6; // ns per bin
    double V_per_unit = 1./16384.*2.*1.e3; // mV / ADC unit

    // Pulse-finding threshold
    double threshold = 200; // mV
    threshold = threshold/V_per_unit;

    // Start and end of pulse ranges
    int start[32];
    int end[32];
    memset(start, 0, 32*sizeof(int));
    memset(end, 210, 32*sizeof(int));

    // Add file names to analyze here
    std::vector<std::string> filenames;
    filenames.push_back("~/Dropbox/Work/CERN_2023/root_files/root_run_000292.root");
    double energy = 0.24;
    std::string polarity = "plus";

    // Open files
    TChain *digitizer0 = new TChain("midas_data_D300");
    TChain *digitizer1 = new TChain("midas_data_D301");
    TChain *digitizer2 = new TChain("midas_data_D302");
    TChain *digitizer3 = new TChain("midas_data_D303"); 
    for(int i=0; i<(int)filenames.size(); i++){
        digitizer0->AddFile(filenames[i].c_str());
        digitizer1->AddFile(filenames[i].c_str());
        digitizer2->AddFile(filenames[i].c_str());
        digitizer3->AddFile(filenames[i].c_str());
    }  

    // Prepare output file
    char output_name[1024];
    sprintf(output_name, "~/Dropbox/Work/CERN_2023/analysis_files/pulse_information-%f-%s.root",energy,polarity.c_str());
    TFile* output_file = new TFile(output_name,"RECREATE");
    TTree* output_tree = new TTree("pulse_information", "Pulse information");
    double output_STT_times[32];
    double output_CFD_times[32];
    double output_baselines[32];
    double output_amplitudes[32];
    double output_charges[32];
    std::vector<int> output_STT_pulses[32];
    for (int j=0; j<32; j++) {
        output_tree->Branch(Form("STT_times%d",j), &output_STT_times[j], Form("STT_times%d/D",j));
        output_tree->Branch(Form("CFD_times%d",j), &output_CFD_times[j], Form("CFD_times%d/D",j));
        output_tree->Branch(Form("amplitudes%d",j), &output_amplitudes[j], Form("amplitudes%d/D",j));
        output_tree->Branch(Form("charges%d",j), &output_charges[j], Form("charges%d/D",j));
        output_tree->Branch(Form("baselines%d",j), &output_baselines[j], Form("baselines%d/D",j));
        output_tree->Branch(Form("STT_multipeaks%d",j), &output_STT_pulses[j], Form("STT_multipeaks%d",j));
    }

    // Initialize vectors for waveforms
    std::vector<double>* waveforms[32];
    for (int i=0; i<32; i++) {waveforms[i]=NULL;}
    for(int i=0; i<8; i++){
        digitizer0->SetBranchAddress(Form("Channel%d",i), &(waveforms[i+0*8]));
        digitizer1->SetBranchAddress(Form("Channel%d",i), &(waveforms[i+1*8]));
        digitizer2->SetBranchAddress(Form("Channel%d",i), &(waveforms[i+2*8]));
        digitizer3->SetBranchAddress(Form("Channel%d",i), &(waveforms[i+3*8]));
    }
    
    std::cout << "**Pre-analysis summary** " << std::endl;
    std::cout << "Digitizer 0: " << digitizer0->GetEntries() << std::endl;
    std::cout << "Digitizer 1: " << digitizer1->GetEntries() << std::endl;
    std::cout << "Digitizer 2: " << digitizer2->GetEntries() << std::endl;
    std::cout << "Digitizer 3: " << digitizer3->GetEntries() << std::endl;

    // Start analysis
    int baseline_samples = 20;

    for (int i=0; i<digitizer0->GetEntries(); i++) {

        printf("Progress: %4.2f%%.\r", (double(i)+1)/double(digitizer0->GetEntries())*100.);

        // Load data
        digitizer0->GetEntry(i);
        digitizer1->GetEntry(i);
        digitizer2->GetEntry(i);
        digitizer3->GetEntry(i);

        // Initialize waveform information
        double baselines[32]; memset(baselines, 0., sizeof(double));
        int global_times[32]; memset(global_times, 0, sizeof(int));
        double STT_times[32]; memset(STT_times, 0., sizeof(double));
        double CFD_times[32]; memset(CFD_times, 0., sizeof(double));
        double charges[32]; memset(charges, 0., sizeof(double));
        double amplitudes[32]; memset(amplitudes, 0., sizeof(double));
        std::vector<int> pulses[32];
        for (int j=0; j<32; j++) {
            output_STT_pulses[j].clear();
            pulses[j].clear();
        }

        // Make analysis
        for (int j=0; j<32; j++) {
            if (waveforms[j]->size()>0) {
                baselines[j]  = waveform_baseline(waveforms[j], baseline_samples);
                STT_times[j]  = global_timing(waveforms[j], baselines[j], 0, 200);
                CFD_times[j]  = CFD_timing(waveforms[j], baselines[j], STT_times[j], 0, 200, 0.1);
                charges[j]    = pulse_charge(waveforms[j], baselines[j], STT_times[j], 20);
                amplitudes[j] = waveforms[j]->at(STT_times[j]);
                std::vector<int> temporary_pulses = count_pulses(waveforms[j], baselines[j], threshold);
                for (int k=0; k < temporary_pulses.size(); k++) {
                    pulses[j].push_back(temporary_pulses.at(k)); 
                }
            }
        }

    // Save pulse information to root file
        for (int j=0; j<32; j++) {
            output_baselines[j]  = baselines[j];
            output_STT_times[j]  = STT_times[j];
            output_CFD_times[j]  = CFD_times[j];
            output_charges[j]    = charges[j];
            output_amplitudes[j] = amplitudes[j];
            for (int k=0; k < pulses[j].size(); k++) {
                output_STT_pulses[j].push_back(pulses[j].at(k));
            }
        }
        output_tree->Fill();
    }
    output_file->Write();
    printf("Progress: 100.00%%\n");
}
