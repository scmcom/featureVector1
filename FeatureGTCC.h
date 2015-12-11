//FeatureGTCC.h

#pragma once

#ifndef __FeatureGTCC__
#define __FeatureGTCC__

#include "FeatureExtraction.h"

#define BW_CORRECTION       1.0190
#define VERY_SMALL_NUMBER   1e-200
#ifndef M_PI
#define M_PI                3.14159265358979323846
#endif

#define myMax(x,y)  ((x) > (y) ? (x) : (y))
#define myMod(x,y)  ((x) - (y) * floor((x) / (y) ))
#define erb(x)      (24.7 * (4.37e-3 * (x) + 1.0 ))


class FeatureGTCC : public FeatureExtraction
{
public:
    FeatureGTCC(std::vector<SAMPLE>* signal, std::vector<std::vector<SAMPLE>> magSpec, int fs, int winSize, int hopSize) : FeatureExtraction(signal, magSpec, fs, winSize, hopSize){
        set_fs(fs);
        set_winSize(winSize);
        set_hopSize(hopSize);
        calculate_feature(signal, magSpec);

        id = FD_GTCC;
    }

    virtual void calculate_feature(std::vector<SAMPLE>* signal, std::vector<std::vector<SAMPLE>> magSpec){

        m_feature.clear();

        int numCoeffs = 26; // 26 by default
        set_nCols(numCoeffs);
        set_nRows(magSpec.size());

        int numfBins = magSpec[0].size();
        int nCols = get_nCols();
        int nRows = get_nRows();
        int nFFT = get_winSize();
        int fs = get_fs();
        m_feature.resize(get_nRows(), std::vector<SAMPLE>(get_nCols(), 0));
        SAMPLE fBinSize = fs / (SAMPLE)nFFT;

        SAMPLE lowFreqLimit = 0;
        SAMPLE upperFreqLimit = fs / 2;



        // get the filter points in fBin index
        std::vector<int> filterPoints;
        filterPoints = createFilterPoints(lowFreqLimit, upperFreqLimit, numCoeffs, nFFT, fs);

        std::vector<std::vector<SAMPLE>> filterBanks;
        filterBanks = createFilterBanks(filterPoints, numCoeffs, numfBins);

        std::vector<std::vector<SAMPLE>> filtSpecSum;
        filtSpecSum.resize(nRows, std::vector<SAMPLE>(nCols, 0));

        SAMPLE fSum = 0.0;
        for (int rowIdx = 0; rowIdx < nRows; ++rowIdx){
            for (int colIdx = 0; colIdx < numCoeffs; ++colIdx){
                for (int k = 0; k < numfBins; ++k){
                    fSum += magSpec[rowIdx][k] * filterBanks[colIdx][k];
                }
                filtSpecSum[rowIdx][colIdx] = fSum;
                fSum = 0.0;
            }
        }

        std::vector<std::vector<SAMPLE>> logFiltSpecSum;
        logFiltSpecSum.resize(nRows, std::vector<SAMPLE>(nCols, 0));

        for (int rowIdx = 0; rowIdx < nRows; rowIdx++){
            for (int colIdx = 0; colIdx < numCoeffs; colIdx++){
                logFiltSpecSum[rowIdx][colIdx] = log(filtSpecSum[rowIdx][colIdx]);
            }
        }

        for (int frameIdx = 0; frameIdx < nRows; frameIdx++){
            m_feature[frameIdx] = dct(logFiltSpecSum[frameIdx]);
        }

        return;
    }

    SAMPLE freq_to_erb(SAMPLE freqVal){
        SAMPLE erbVal;
        erbVal = 21.4 * log(0.00437*(freqVal)+1);
        //erbVal = erb(freqVal);
        return erbVal;
    }

    SAMPLE erb_to_freq(SAMPLE erbVal){
        SAMPLE freqVal;
        freqVal = ((pow(10, erbVal/21.4) - 1) / 0.00437);
        return freqVal;
    }

    std::vector<int> createFilterPoints(SAMPLE lowFreqLimit, SAMPLE upperFreqLimit, int numCoeffs, int nFFT, int fs){ //range: 0 - fs/2
        int numPoints = numCoeffs + 2;
        std::vector<int> filterPoints(numPoints);
        SAMPLE lowerbLimit = freq_to_erb(lowFreqLimit);
        SAMPLE uppererbLimit = freq_to_erb(upperFreqLimit);
        SAMPLE erbWinSize = (uppererbLimit - lowerbLimit) / (SAMPLE)(numPoints - 1);
        SAMPLE erbVal;
        SAMPLE freqVal;

        for (int i = 0; i < numPoints; i++){
            erbVal = lowerbLimit + i*erbWinSize;
            freqVal = erb_to_freq(erbVal);
            filterPoints[i] = floor((nFFT + 1)*freqVal / (SAMPLE)fs);
        }
        return filterPoints;
    }

    std::vector<std::vector<SAMPLE>> createFilterBanks(std::vector<int> filterPoints, int numCoeffs, int numfBins){
        std::vector<std::vector<SAMPLE>> filterBanks;
        filterBanks.resize(numCoeffs, std::vector<SAMPLE>(numfBins, 0));

        for (int m = 1; m < numCoeffs + 1; m++){
            for (int k = 0; k < numfBins; k++){
                if (k < filterPoints[m - 1]){
                    filterBanks[m - 1][k] = 0;
                }
                else if (k >= filterPoints[m - 1] && k <= filterPoints[m]){
                    filterBanks[m - 1][k] = (k - filterPoints[m - 1]) / (filterPoints[m] - filterPoints[m - 1]);
                }
                else if (k >= filterPoints[m] && k <= filterPoints[m + 1]){
                    filterBanks[m - 1][k] = (filterPoints[m + 1] - k) / (filterPoints[m + 1] - filterPoints[m]);
                }
                else if (k > filterPoints[m + 1]){
                    filterBanks[m - 1][k] = 0;
                }
            }
        }

        return filterBanks; //a m x k matrix of filter bank values
    }

    std::vector<SAMPLE> dct(std::vector<SAMPLE> logFiltSpecSum){
        int numCoeffs = logFiltSpecSum.size();
        std::vector<SAMPLE> retVal(numCoeffs);

        SAMPLE sum = 0.0;
        for (int k = 0; k < numCoeffs; k++){
            for (int n = 0; n < numCoeffs; n++){
                sum += logFiltSpecSum[n] * cos((M_PI*k / (SAMPLE)numCoeffs)*(n - .5));
            }
            retVal[k] = sum;
            sum = 0.0;
        }

        return retVal;
    }

};

#endif
