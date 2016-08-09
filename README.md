# tth-analyze-mem

Matrix element method for ttH, H->tautau channel

## Building instructions

### Dependencies

In order to build the whole project, you first need to install numerical integration library VAMP by following these steps:

```bash
mkdir $CMSSW_BASE/VAMP
wget http://www.hepforge.org/archive/whizard/vamp-2.3.0.tar.gz -P $CMSSW_BASE/VAMP
tar zxvf $CMSSW_BASE/VAMP/vamp-2.3.0.tar.gz -C $CMSSW_BASE/VAMP
mkdir -p $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/share/doc/vamp
cd $CMSSW_BASE/VAMP/vamp-2.3.0
./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.3.0/prefix
make -j4
make install
cp $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH
```

Create a new file with the following contents

```xml
<tool name="vamp" version="2.3.0">
  <lib name="vamp"/>
  <client>
    <environment name="LIBDIR" default="<your lib dir>"/>
    <environment name="INCLUDE" default="<your include dir>"/>
  </client>
  <use name="f77compiler"/>
</tool>
```

where the full paths are

```bash
<your lib dir>     = $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/lib
<your include dir> = $CMSSW_BASE/VAMP/vamp-2.3.0/prefix/include/vamp
```

Save it to `$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/vamp.xml`. Next up try to set up with `scram setup vamp` and verify the installation with `scram tool info vamp`. You should see output along the lines of

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
You can check other versions of VAMP [here](https://www.hepforge.org/downloads/whizard) (scroll down).

### The library itself

Clone the library into `$CMSSW_BASE/src/tthAnalysis/tthMEM`:

```bash
git clone git@github.com:HEP-KBFI/tth-analyze-mem.git $CMSSW_BASE/src/tthAnalysis/tthMEM
```

and build it with `scram b -j8` as usual.

## Running an example

At the time of writing, the project includes one executable, `runMEM`, which takes a Python configuration file as its argument. The program needs `.root` file(s) containing selected events; an example configuration file is provided in `data/runMEM_3l1tau_cfg.py` and example data file in `data` subdirectory. The output file is the original file + signal and background probabilities found with the MEM.

In order to run the MEM via SLURM, one has to call `python createJobs_3l1tau.py` in any directory and follow instructions on screen. Note that this assumes you have run the tth analysis: https://github.com/HEP-KBFI/tth-htt/ (with `select_root_output` set to `True`).

The project also includes some unit tests, which can be run with `scram b -j8 runtests` after building the project.
