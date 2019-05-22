# tth-analyze-mem

Matrix element method for ttH, H->tautau channel

## Building instructions

Clone the library into `$CMSSW_BASE/src/tthAnalysis/tthMEM`:

```bash
git clone https://github.com/HEP-KBFI/tth-analyze-mem.git $CMSSW_BASE/src/tthAnalysis/tthMEM
```

In order to build the whole project, you also need to install numerical integration library VAMP by following these steps:

```bash
mkdir $CMSSW_BASE/VAMP
wget http://whizard.hepforge.org/oldsrc/vamp-2.3.0.tar.gz -P $CMSSW_BASE/VAMP
tar zxvf $CMSSW_BASE/VAMP/vamp-2.3.0.tar.gz -C $CMSSW_BASE/VAMP
mkdir -p $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/share/doc/vamp
cd $CMSSW_BASE/VAMP/vamp-2.3.0
./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.3.0/prefix
make -j4
make install
cp $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH
cp $CMSSW_BASE/src/tthAnalysis/tthMEM/test/vamp.xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/vamp.xml
```
Next up try to set up with `scram setup vamp` and verify the installation with `scram tool info vamp`. You should see output along the lines of

```bash
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Name : vamp
Version : 2.3.0
++++++++++++++++++++
SCRAM_PROJECT=no
LIB=vamp
LIBDIR=.../VAMP/vamp-2.3.0/prefix/lib
INCLUDE=.../VAMP/vamp-2.3.0/prefix/include/vamp
USE=f77compiler
```
You can check other versions of VAMP [here](http://whizard.hepforge.org/vamp.html).

Finally, build it with `scram b -j8` as usual.

## Running an example

At the time of writing, the project includes one executable, `runMEM`, which takes a Python configuration file as its argument. The program needs `.root` file(s) containing selected events; an example configuration file is provided in `data/runMEM_3l1tau_cfg.py` and example data file in `data` subdirectory. The output file is the original file + signal and background probabilities found with the MEM.

In order to run the MEM via SLURM, one has to call `python createJobs_3l1tau.py` in any directory and follow instructions on screen. Note that this assumes you have run the tth analysis: https://github.com/HEP-KBFI/tth-htt/ (with `select_root_output` set to `True`).

The project also includes some unit tests, which can be run with `scram b -j8 runtests` after building the project.

## Disabling logging

If you want to gain in speed, you are advised to build the project w/o no logging whatsoever. This vcan be achieved w/ the following command
```bash
USER_CPPFLAGS="-DDISABLE_LOGGING" scram b -j8
```

## Valgrind commands

Here we provide some example valgrind commands which help to debug the application. Of course, in order to benefit from the Valgrind bundle, one should definitely build the project w/ debug symbols enabled:
```bash
scram b -j8 USER_CXXFLAGS="-g -Wuninitialized"
```

<details>
<summary>Commands</summary>

Memory leak detection:
```bash
valgrind --tool=memcheck `cmsvgsupp`                               \
--suppressions=$CMSSW_BASE/src/tthAnalysis/tthMEM/data/tthMEM.supp \
--leak-check=yes                                                   \
--show-reachable=yes                                               \
--num-callers=20                                                   \
--track-fds=yes                                                    \
--track-origins=yes                                                \
--log-file="valgrind.log"                                          \
runMEM_3l1tau python/runMEM_3l1tau_2016_cfg.py                     \
&> out.log
```
The second line suppresses leaks from 3rd party libraries linked to `runMEM_3l1tau`. If you only see Python-related leaks, then you're good.

Memory consumption:
```bash
valgrind --tool=massif                         \
--depth=40                                     \
--time-stamp=yes                               \
--time-unit=ms                                 \
--threshold=0.1                                \
runMEM_3l1tau python/runMEM_3l1tau_2016_cfg.py \
&> out.log
```

Callgraph:
```bash
valgrind --tool=callgrind                      \
runMEM_3l1tau python/runMEM_3l1tau_2016_cfg.py \
&> out.log
```
</details>
