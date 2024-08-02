# Hourglass
Hourglass is an adaptive range filter that can adapt to recurring false positives to achieve optimal performance on skewed and adversarial range queries.

# Settings
- cmake version 3.28.3
- gcc version 13.2.0
- OpenSSL 3.0.13


# Build
``` 
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j16
```

# Run the example
We give a simple query example in example.cpp 
```
cd build
./hourglass
```
