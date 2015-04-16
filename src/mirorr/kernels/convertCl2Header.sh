#!/bin/bash

echo "const std::string kMeanVarianceString =" > kMeanVariance.h
#sed -e 's/\\/\\\\/g;s/"/\\"/g;s/  /\\t/g;s/^/"/;s/$/\\n"/' kMeanVariance.cl >> kMeanVariance.h
sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/' kMeanVariance.cl >> kMeanVariance.h
echo ";" >> kMeanVariance.h

echo "const std::string kMatchString =" > kMatch.h 
#sed -e 's/\\/\\\\/g;s/"/\\"/g;s/  /\\t/g;s/^/"/;s/$/\\n"/' kMatch.cl >> kMatch.h
sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/' kMatch.cl >> kMatch.h  
echo ";" >> kMatch.h

# new kernel with moving block mean compute optimization
echo "const std::string newKmatchString =" > optKmatch.h 
sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/' optKmatch.cl >> optKmatch.h  
echo ";" >> optKmatch.h

# new kernel with moving block mean compute optimization - maks version
echo "const std::string newKmatchString =" > optKmatchMask.h 
sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/' optKmatchMask.cl >> optKmatchMask.h  
echo ";" >> optKmatchMask.h
