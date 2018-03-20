//
//  otsu_multi.h
//  DSC_segment
//
//  Created by Tuan Nguyen Trung on 6/14/17.
//  Copyright © 2017 Asger Nyman Christiansen. All rights reserved.
//

#ifndef otsu_multi_h
#define otsu_multi_h

#include <stdio.h>
#include <iostream>
#include <vector>

// Base on https://github.com/hipersayanX/MultiOtsuThreshold
// A TlresholdSelectionMethod fromGray-LevelHistograms
//      NOBUYUKI OTSU
//A Fast Algorithm for Multilevel Thresholding
//      PING-SUNG LIAO, TSE-SHENG CHEN* AND PAU-CHOO CHUNG+

inline std::vector<double> buildTables(std::vector<int> histogram)
{
    std::vector<long> P(histogram.size()+1);
    std::vector<long> S(histogram.size()+1);
    
    P[0] = 0;
    S[0] = 0;
    
    long sumP = 0;
    long sumS = 0;
    
    for (int i = 0; i < histogram.size(); i++) {
        sumP += long(histogram[i]);
        sumS += long(i * histogram[i]);
        P[i + 1] = sumP;
        S[i + 1] = sumS;
    }
    
    // Calculate the between-class variance for the interval u-v
    std::vector<double> H(histogram.size() * histogram.size(), 0.);
    
    for (int u = 0; u < histogram.size(); u++) {
        double * hLine = H.data() + u * histogram.size();
        
        for (int v = u + 1; v < histogram.size(); v++)
        {
            if(std::abs(P[v] - P[u]) > 0.000001)
                hLine[v] = (S[v] - S[u])*(S[v] - S[u]) / (P[v] - P[u]);
            else
            {
                
            }
        }
    }
    
    return H;
}

void for_loop(double *maxSum,
              std::vector<int> *thresholds,
              const std::vector<double> &H,
              int u,
              int vmax,
              int level,
              int levels,
              std::vector<int> *index)
{
    size_t classes = index->size() - 1;
    
    for (int i = u; i < vmax; i++) {
        (*index)[level] = i;
        
        if (level + 1 >= classes) {
            // Reached the end of the for loop.
            
            // Calculate the quadratic sum of al intervals.
            double sum = 0.;
            
            for (int c = 0; c < classes; c++) {
                int u = index->at(c);
                int v = index->at(c + 1);
                sum += H[v + u * levels];
            }
            
            if (*maxSum < sum) {
                // Return calculated threshold.
//                *thresholds = index->mid(1, thresholds->size());
                *thresholds = std::vector<int>(index->begin() + 1, index->begin()+1+thresholds->size());
                *maxSum = sum;
            }
        } else
            // Start a new for loop level, one position after current one.
            for_loop(maxSum,
                     thresholds,
                     H,
                     i + 1,
                     vmax + 1,
                     level + 1,
                     levels,
                     index);
    }
}

std::vector<int> otsu_muti(std::vector<int> histogram_in, int classes)
{
    std::vector<int> histogram = histogram_in;
    
    for (int i = 0; i < histogram.size(); i++) // Avoid exception
        histogram[i]++;
    
    double maxSum = 0;
    std::vector<int> thresholds(classes - 1, 0);
    std::vector<double> H = buildTables(histogram);
    
    std::vector<int> index(classes + 1);
    index[0] = 0;
    index[index.size() - 1] = histogram.size() - 1;
    
    for_loop(&maxSum,
             &thresholds,
             H,
             1,
             histogram.size() - classes + 1,
             1,
             histogram.size(), &index);
    
    return thresholds;
}

#endif /* otsu_multi_h */
