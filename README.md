# tth-analyze-mem

Matrix element method for ttH, H->tautau channel

## Building instructions

In order to build the whole project, you first need to install numerical integration library VAMP by following these steps:

```bash
mkdir $CMSSW_BASE/VAMP
wget http://www.hepforge.org/archive/whizard/vamp-2.2.8.tar.gz -P $CMSSW_BASE/VAMP
tar zxvf $CMSSW_BASE/VAMP/vamp-2.2.8.tar.gz -C $CMSSW_BASE/VAMP
mkdir -p $CMSSW_BASE/VAMP/vamp-2.2.8/prefix/share/doc/vamp
./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.2.8/prefix
make -j4
make install
cp $CMSSW_BASE/VAMP/vamp-2.2.8/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH
```

Create a new file with the following contents

```xml
<tool name="vamp" version="2.2.8">
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
<your lib dir>     = $CMSSW_BASE/VAMP/vamp-2.2.8/prefix/lib
<your include dir> = $CMSSW_BASE/VAMP/vamp-2.2.8/prefix/include/vamp
```

Save it to `$CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/vamp.xml`.
Next up try to set up with `scram setup vamp` and verify the installation with `scram tool info vamp`.
You should see output along the lines of
```bash
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Name : vamp
Version : 2.2.8
++++++++++++++++++++
SCRAM_PROJECT=no
LIB=vamp
LIBDIR=.../VAMP/vamp-2.2.8/prefix/lib
INCLUDE=.../VAMP/vamp-2.2.8/prefix/include/vamp
USE=f77compiler
```

