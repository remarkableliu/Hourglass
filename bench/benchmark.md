# Build
``` 
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j16
```

# Generate workloads
```
bash ./bench/scripts/download_datasets.sh
bash ./bench/scripts/generate_datasets.sh ./build ./real_datasets
```

# Run benchmark
```
bash ./bench/scripts/execute_tests.sh ./build {YOUR_WORKLOAD_PATH}
```