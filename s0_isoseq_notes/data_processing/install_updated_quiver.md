# Installing an Updated version of Quiver (Genomic Consensus) in a virtual env
Brian Uapinyoying
Updated: 09/26/2016

There is a bug that can sometimes cause a `RuntimeError: Invalid input!` error to occur in the IceQuiver step that results in a corrupt cluster file. This stops the quiver step from continuing.  The version of GenomicConsensus/Quiver that ships with smrtAnalysis 2.3.0 bundle is outdated and the bug has been addressed in the new GenomicConsensus version, but there hasn't been a new bundle released by PacBio that updates the code.  Here are instructions on how to work around the issue.

### Example of the error output

**IceQuiver.log**
```
...
[DEBUG] [IceQuiver] 2016-09-09 12:05:56: [IceQuiver] Concatenating sam files between 25914 and 26033.
[DEBUG] [IceQuiver] 2016-09-09 12:05:59: [IceQuiver] Concatenation done
[DEBUG] [IceQuiver] 2016-09-09 12:05:59: [IceQuiver] Creating quiver cmds for c25914 to c26033
[DEBUG] [IceQuiver] 2016-09-09 12:05:59: [IceQuiver] Creating quiver bash script /home/brian/3-6kb_grp1_edl_8cell/cluster_out/quivered/c25914to26033.sh for c25914 to c26033.
[DEBUG] [IceQuiver] 2016-09-09 12:06:00: [IceQuiver] Submitting todo quiver jobs.
[DEBUG] [IceQuiver] 2016-09-09 12:06:00: [IceQuiver] Running CMD: bash /home/brian/3-6kb_grp1_edl_8cell/cluster_out/quivered/c25914to26033.sh 1>/home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.olog 2>/home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.elog
[ERROR] [IceQuiver] 2016-09-09 12:12:55: [IceQuiver] CMD exited with a non-zero code: bash /home/brian/cluster_out/quivered/c25914to26033.sh 1>/home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.olog 2>/home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.elog, 
Failed to run Quiver
Error log: /home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.elog
Out log: /home/brian/3-6kb_grp1_edl_8cell/cluster_out/log/quivered/c25914to26033.sh.olog
```

**c25914to26033.sh.elog**
```
...
[INFO] Quiver operating on c25937:2995-3505
[INFO] Quiver operating on c25937:3495-3716
[INFO] Quiver operating on c25942:61-505
Process QuiverWorkerProcess-14:
Traceback (most recent call last):
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/multiprocessing/process.py", line 258, in _bootstrap
    self.run()
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/Worker.py", line 91, in run
    self._run()
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/Worker.py", line 77, in _run
    result = self.onChunk(datum)
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/quiver/quiver.py", line 196, in onChunk
    refContig, options.coverage, self.quiverConfig)
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/quiver/quiver.py", line 116, in consensusAndVariantsForWindow
    quiverConfig)
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/quiver/utils.py", line 383, in consensusForAlignments
    refineDinucleotideRepeats(mms)
  File "/c1/apps/smrt_analysis/2.3.0/current/redist/python2.7/lib/python2.7/site-packages/GenomicConsensus/quiver/utils.py", line 164, in refineDinucleotideRepeats
    return cc.RefineDinucleotideRepeats(mms)
RuntimeError: Invalid input!
```

#### Setup a virtualenv

Copy the code below into a shell script and run at command line in the folder you would like to create it.

Warning:
* Make sure you unload smrtnalaysis 2.3.0 from your enviornment
* Do not load the SmrtShell or you will call the original quiver after installation.

```
!# /bin/bash
set -vex
export VENV_GC=${PWD}/VENV_GC
wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz
tar zxf virtualenv-1.11.6.tar.gz -C /tmp/
python /tmp/virtualenv-1.11.6/virtualenv.py --system-site-packages $VENV_GC
source ${VENV_GC}/bin/activate
```

#### Attempt to install GenomicConsensus
First try to install GenomicConsensus directly. Likely there will be an error asking for ConsensusCore as a dependency.

Download GenomicConsensus, go into folder and attempt to install:

* `$ git clone https://github.com/PacificBiosciences/GenomicConsensus.git`
* `$ cd GenomicConsensus`
* `$ python setup.py install`

#### Installing ConsensusCore
I had trouble with prereqs/dependencies for installing ConsensusCore. Specifically, boost version 1.61.0.  Best to just download and build the specified old boost 1.47.0 version. The latest verison 1.61.0 didn't seem to have the same directory tree after the build process (`./include` dir no longer there). Swig version 3.0.10 works fine, so I assume new versions are no problem, but I compiled it too just in case. You can test the one on your system first 

1. Download ConsensusCore: `$ git clone https://github.com/PacificBiosciences/ConsensusCore.git`

2. Manually download [boost version 1.47.0](http://sourceforge.net/projects/boost/files/boost/1.47.0/) & extract `$ tar -zxvf boost_1_47_0.tar.gz`

3. Go into downloaded boost folder and run bootstrap, then install script with custom prefix. Plan a custom boost install directory (that you have install permissions). Directory names with underscores won't work so name your new boost directory with dashes and dots only. I decided to place into the `./ConsensusCore/tools/boost-1.47.0`
    * `$ cd boost_1_47_0`
    * `$ ./bootstrap.sh`       
    * `$ ./b2 install --prefix=/home/brian/vGC/ConsensusCore/tools/boost-1.47.0`

4. Manually [download Swig](http://www.swig.org/download.html). Latest version works just fine. I used version 3.0.10. Make sure to add the `--prefix` option with your intended install directory.
    * `$ tar -zxvf swig-3.0.10.tar.gz`
    * `$ cd swig-3.0.10`
    * `$ ./configure --prefix=/home/brian/vGC/ConsensusCore/tools/swig-3.0.10`
    * `$ make`
    * `$ make install`

5. Installing ConsensusCore itself. Once finished building swig and boost, go back into ConsensusCore directory
    * `$ cd ../ConsensusCore`
    
6. For the install command you *must* specify the path to the boost parent folder `/<path>/<to>/boost-1.47.0/include`, this folder does not seem to exist anymore in the latest version of boost and will cause errors if you use say version 1.60.0. In the same command, you will also need to specify the actual swig binary path `/<path>/<to>/swig-3.0.10/bin/swig`. Run the setup.py with boost and swig paths provided in the same command e.g.
    *   `$ python setup.py install --boost=/home/brian/vGC/ConsensusCore/tools/boost-1.47.0/include --swig=/home/brianu/vGC/ConsensusCore/tools/swig-3.0.10/bin/swig`

7. Finally return to the GenomicConsensus folder and install
    * `$ cd ../GenomicConsensus`
    * `$ python setup.py install`
    * Test `$ quiver --version` you should get `2.1.0`
