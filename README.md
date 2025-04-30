# TopoA: Topology Preserving Augmentation of Lossy Compressors

Related Publication: 
- Nathaniel Gorski, Xin Liang, Hanqi Guo, Lin Yan, and Bei Wang. A general framework for augmenting lossy compressors with topological guarantees. IEEE TVCG journal track at IEEE PacificVis 2025. 

[arxiv](https://arxiv.org/abs/2502.14022)

## Installation
TopoA uses VTK as a dependency. We suspect that the latest version will suffice, but we have not tested builds later than VTK 9.3.0.

```
cd TopoA
mkdir build && cd build
cmake -DVTK_DIR=[path_to_vtk] ..
make
```

If ```path_to_vtk``` is a relative path, then it should be relative to ```CMakeLists.txt``` rather, than the ```build``` folder.

Currently, we do not support Windows. Windows users should be able to use TopoA using WSL.

## Usage
Run ```TopoA -h``` for syntax options. In this section, we provide additional clarifications.

### Basic Usage

- TopoA currently can only compress 3D scalar fields in VTK image format (.vti).

- When compressing, you must specify a persistence threshold $\varepsilon$ and error bound $\xi$. Both values should be expressed as a percentage of the range of the scalar field, and should be within $(0,1]$. For example, a value of $\varepsilon = 0.04$ means that the persistence threshold should be 4% of the range. These variables are specified using ```-eps``` and ```-xi``` respectively.

- TopoA creates and deletes many temporary files. By default, this occurs in the CWD. You can specify the folder where these files are created using ```-of```.


### Specifying a Base Compressor

- We define a "base compressor" as the compressor that is augmented by TopoA. Currently, TopoA only supports five base comprssors: [ZFP](https://github.com/LLNL/zfp), [SZ3](https://github.com/szcompressor/SZ3), [TTHRESH](https://github.com/rballester/tthresh), cubic spline interpolation (CSI), and [Neurcomp](https://github.com/matthewberger/neurcomp?tab=readme-ov-file).
- The base compressor is specified by using the ```-bc``` flag followed by one of: ```ZFP```, ```SZ3```, ```TTHRESH```, ```CSI```, or ```Neurcomp```.
- The cubic spline interpolation model is built into TopoA and requires no additional setup.
- To augment one of ZFP, SZ3, or TTHRESH, you must place the binary executable for that compressor into the CWD.
	- You can alternatively use ```-bcFolder``` to specify the folder that contains the binary.
	- The binaries must have their original names i.e. ```sz3```, ```zfp```, and ```tthresh```.
- Neurcomp compresses data by training a neural network. As such, it may be desirable to run Neurcomp on a GPU. To enable the use of GPUs, Neurcomp must be run separately from TopoA. We provide steps below:
	1. Clone the Neurcomp repo into the same folder as TopoA, or clone it into the folder pointed to by ```-bcFolder```.
	2. Convert your volume into Numpy format (.npy) so that Neurcomp can compress it.
	3. Run Neurcomp to compress your scalar field (before running TopoA). Neurcomp will output two files: ```thenet.pth``` and ```thenet.json```.
	4. Place ```thenet.pth``` and ```thenet.json``` into the folder where TopoA will save its temporary files (see the previous section).
	5. Copy and paste the .npy file created in step ii into the same folder as ```thenet.pth``` and ```thenet.json```. Rename it to ```dummy.npy```
		- When decompressing, Neurcomp requires you to specify the original file in order to obtain the size of the original file and measure the reconstruction quality. Because we do not use the reconstruction quality outputted by Neurcomp, you could also make ```dummy.npy``` some other numpy array with the same dimensions as the input file, but with different values.
	6. Run TopoA as normal using ```-bc Neurcomp```.

### Experiment Flag

If you plan to benchmark TopoA, it may be convenient to use the ```-experiment``` flag. When used, this flag will cause TopoA to run the compressor and decompressor on a given input data file. It will collect the evaluation metrics from our publication and print them to the console. You can choose to write these files to a csv file using the ```-csv``` flag. If you specify a CSV file that already exists, the results will be appended onto the end of that file.

You do not need to specify filenames for the compressed or decompressed outputs. If you decide not to, then TopoA will use temporary files, and delete them when it is finished.

