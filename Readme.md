
#### Build instructions

###### Windows

Create a build directory, e.g. `C:\dev\git\Paclets\GUROBILink\build`
```
cd C:\dev\git\Paclets\GUROBILink
mkdir build
cd build
```
Configure the build
```
cmake -G "Visual Studio 15 2017 Win64" ^
      -DMATHEMATICA_LAYOUT="C:\Program Files\Wolfram Research\Mathematica\12.2"^
      -DINSTALL_PDB=ON ^
      ..\CSource
```
Open the solution in Visual Studio
```
cmake --open .
````
Build from the command line
```
cmake --build . --config Debug --target INSTALL
```
The built paclet will be assembled in `build\install` and can be enabled by evaluating

```
PacletDirectoryAdd["C:\\dev\\git\\Paclets\\GUROBILink\\build\\install\\GUROBILink"];
```

###### Mac

Create a build directory, e.g. `~/git/Paclets/GUROBILink/build`
```
cd ~/git/Paclets/GUROBILink
mkdir build
cd build
```
Configure the build
```
cmake -DMATHEMATICA_LAYOUT="/Applications/Mathematica 12.2.app/Contents" ../CSource
```
Build the debug configuration from the command line
```
cmake --build . --config Debug --target install
```
The built paclet will be assembled in `build/install` and can be enabled by evaluating

```
PacletDirectoryAdd["~/git/Paclets/GUROBILink/build/install/GUROBILink"];
```

