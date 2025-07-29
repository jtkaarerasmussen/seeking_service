Code associated with the manuscript "Seeking Service: How Client Behavior Determines Cleaning Station Clustering" 

All simulations were done with python 3.12 and cargo 1.87

To run the example, first clone and navagate to the repository. Numpy and Matplotlib are required and can be installed with
```
pip3 install matplotlib numpy
```

Next, you have to the compile the rust simulation, which can be done with cargo as follows:
```
(cd sim && cargo build --release)
```
At this point, you can run the example with
```
python3 example_usage.py
```
