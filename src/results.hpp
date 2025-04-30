#pragma once

class Results{
public:
    double ratio;
    double totalTime;
    double baseTime;
    double errorBoundTime;
    double growthTime;
    double averageTightenTime;
    int numFP = -1;
    int numFN = -1;
    int numExpanded = -1;
    double writeToFileTime;
    double decompressionTime;

    double gtCtTime;
    double decompressedCtTime;

    int unsimplifiedGtCtSize = -1;
    int unsimplifiedDecompressedCtSize = -1;
    int simplifiedCTSize = -1; // should be the same for both :)
};