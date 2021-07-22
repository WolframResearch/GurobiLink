
# GurobiLink for Wolfram Language

GurobiLink is a package implemented in [Wolfram Language](https://www.wolfram.com/language/) and C++
using [LibraryLink](https://reference.wolfram.com/language/guide/LibraryLink.html) that provides an interface
to the [Gurobi](https://www.gurobi.com/) numerical optimization solver. Gurobi can be used
with the same ease as the built-in solvers and has excellent performance for particularly large or complicated
(e.g. mixed-integer) models. It supports continuous and mixed-integer linear programming (LP),
quadratic programming (QP), quadratically-constrained programming (QCP) and other classes of problems.

The package requires a Gurobi [license](https://reference.wolfram.com/language/workflow/GetALicenseForGUROBI.html)
(free academic licenses and evaluation licenses for commercial users are available, please see the
[Gurobi documentation](https://www.gurobi.com/documentation/quickstart.html) for license installation instructions).
It works with Gurobi version 9.0 and later and is bundled with Wolfram Language products (such as [Mathematica](https://www.wolfram.com/mathematica/)
or [Wolfram Desktop](https://www.wolfram.com/desktop/)) version [12.2](https://writings.stephenwolfram.com/2020/12/launching-version-12-2-of-wolfram-language-mathematica-228-new-functions-and-much-more/) and later.

GurobiLink makes the solver accessible as a plug-in through the Wolfram Language
[optimization method framework](https://reference.wolfram.com/language/OptimizationMethodFramework/tutorial/OptimizationMethodFramework.html)
and is used to implement [`Method -> "Gurobi"`](https://reference.wolfram.com/language/ref/method/GUROBI.html) in functons like [`NMinimize`](https://reference.wolfram.com/language/ref/NMinimize.html) or
[`ConvexOptimization`](https://reference.wolfram.com/language/ref/ConvexOptimization.html). This allows seamless integration of the Wolfram Language modeling capabilities
with the high performance of Gurobi, for example

```
In[1]:= ConvexOptimization[Total[x], Norm[x] <= 1, Element[x, Vectors[10, Integers]],
            Method -> "Gurobi"] // AbsoluteTiming

Out[1]= {0.004534, {x -> {0, 0, -1, 0, 0, 0, 0, 0, 0, 0}}}
```

## How to build

The build requires [CMake](https://cmake.org/) 3.15 or later, a C++ compiler with support for C++11, as well
as installations of Wolfram Language and the [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-optimizer-eula/).

The general steps for building the complete paclet are

```
git clone https://github.com/WolframResearch/GurobiLink.git GurobiLink
cd GurobiLink
mkdir build
cd build
cmake <options> ..
cmake --build . --target install
```

which will place the result by default in the `GurobiLink/build/install` directory.

The typical CMake options required for building are

* `WOLFRAM_INSTALLATION` -- Wolfram layout location, for example, `/usr/local/Wolfram/Mathematica/12.3`
* `GUROBI_VERSION` -- version number, for example, `9.1.2`
* `GUROBI_HOME` -- where Gurobi is installed, for example, `/opt/gurobi912/linux64`

It may not always be necessary to specify all of these, since the build system provides reasonable defaults.

The built paclet can then be enabled in a particular Wolfram Language session using
[`PacletDirectoryLoad`](https://reference.wolfram.com/language/ref/PacletDirectoryLoad.html)
or packed into a `.paclet` archive using
[`CreatePacletArchive`](https://reference.wolfram.com/language/ref/CreatePacletArchive.html).

The `.paclet` format is suitable for distribution or permanent installation using
[`PacletInstall`](https://reference.wolfram.com/language/ref/PacletInstall.html).

Some platform-specific examples follow:

#### Windows

Create a build directory, e.g. `C:\dev\git\Paclets\GurobiLink\build`
```
cd C:\dev\git\Paclets\GurobiLink
mkdir build
cd build
```
Configure the build using the 64-bit Visual Studio 2017 generator (see `cmake --help` for
other possible CMake generators)
```
cmake -G "Visual Studio 15 2017 Win64" ^
      -DWOLFRAM_INSTALLATION="C:\Program Files\Wolfram Research\Mathematica\12.3"^
      -DGUROBI_VERSION=9.1.2^
      -DGUROBI_HOME="C:\gurobi912\win64"^
      -DINSTALL_PDB=ON ^
      ..
```
The generated solution can be opened in the Visual Studio IDE using
```
cmake --open .
```
To build, select the appropriate configuration (for example, Debug), right-click
the desired target (for example, INSTALL) and choose 'Build' from the context menu.

Alternatively, to build the debug configuration from the command line use
```
cmake --build . --config Debug --target INSTALL
```
The development version of the paclet will be assembled in `build\install` and can be enabled in a
Wolfram Language session by

```
PacletDirectoryLoad["C:\\dev\\git\\Paclets\\GurobiLink\\build\\install\\GurobiLink"];
```
Sanity-check the build by running the included basic test file:

```
TestReport["C:\\dev\\git\\Paclets\\GurobiLink\\Tests\\GurobiLink_Basic.wlt"]
```

#### macOS

Create a build directory, e.g. `~/git/Paclets/GurobiLink/build`
```
cd ~/git/Paclets/GurobiLink
mkdir build
cd build
```
Configure the build using the default Unix Makefiles generator (see `cmake --help` for
other possible CMake generators)
```
cmake -DWOLFRAM_INSTALLATION="/Applications/Mathematica12.3.app/Contents" \
      -DGUROBI_VERSION=9.1.2 \
      -DGUROBI_HOME=/Library/gurobi912/mac64 \
      ..
```
Build the debug configuration from the command line
```
cmake --build . --config Debug --target install
```
The development version of the paclet will be assembled in `build/install` and can be enabled in a Wolfram Language
session by

```
PacletDirectoryLoad["~/git/Paclets/GurobiLink/build/install/GurobiLink"];
```
Sanity-check the build by running the included basic test file:

```
TestReport["~/git/Paclets/GurobiLink/Tests/GurobiLink_Basic.wlt"]
```

#### Linux

See the macOS instructions, except the build configuration step would be

```
cmake -DWOLFRAM_INSTALLATION="/usr/local/Wolfram/Mathematica/12.3" \
      -DGUROBI_VERSION=9.1.2 \
      -DGUROBI_HOME=/opt/gurobi912/linux64 \
      ..
```

### See also 
 * [LICENSE](LICENSE) - GurobiLink license
 * [CONTRIBUTING.md](CONTRIBUTING.md) - Guidelines for contributing to GurobiLink
