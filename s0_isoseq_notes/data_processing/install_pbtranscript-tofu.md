##  Installing pbtranscript-ToFU
---

Brian Uapinyoying, updated 3/18/16

[Official cDNA_primer page](https://github.com/PacificBiosciences/cDNA_primer)

It's best to use the install script outlined in the [instructions](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu). If you try to use newer versions of the tools outlined in the instructions such as virtualenv, setuptools or bxpython it will likely not work and you won't know why.  I learned the hard way.

The provided script on the instruction page works best. However I had to tweak the script a bit to get it to work on my Ubuntu x64 14.04 LTS setup. 

##### 1.  **Create new folder in home directory**. Place the shell script below in the directory. 

**Note**: the header should be `#!</path/to/smrtanalysis>/smrtcmds/bin/smrtshell`
```bash
!#/opt/smrtanalysis/smrtcmds/bin/smrtshell
set -vex
#module load ccache
#module load zlib
export VENV_TOFU=${PWD}/VENV_TOFU
wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz
tar zxf virtualenv-1.11.6.tar.gz -C /tmp/
python /tmp/virtualenv-1.11.6/virtualenv.py --system-site-packages $VENV_TOFU
source ${VENV_TOFU}/bin/activate

pip install urllib3 # added to script
pip install hmac # added to script

pip install --upgrade pip
pip install --upgrade setuptools
pip install Cython
pip install numpy
pip install bx-python
pip install pysam
pip install pbcore
```

##### 2.  **Run the script**
```bash
bash <script_name.sh>
```

##### 3.  **Exit smrtshell and then run git clone to download install files**
```bash
git clone https://github.com/PacificBiosciences/cDNA_primer.git
```

##### 4.  **Go back into the VENV_TOFU and smrtshell**
```bash
source ${VENV_TOFU}/bin/activate
/opt/smrtanalysis/smrtcmds/bin/smrtshell
```

##### 5. **Create a second shell script with below in it and run**

```bash
set -vex
export VENV_TOFU=${PWD}/VENV_TOFU
source ${VENV_TOFU}/bin/activate  

cd cDNA_primer/pbtranscript-tofu/pbtranscript/pbtools/pbtranscript/branch/C/
python setup2.py build_ext --inplace
cd ../../../../
make

cd ../external_daligner/
unzip  DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod.zip
cd  DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod/
make
cp HPCdaligner HPCmapper LA4Ice LAcat LAcheck LAmerge LAshow LAsort LAsplit daligner $VENV_TOFU/bin

cd ../
unzip DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa.zip
cd  DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa/
make
cp Catrack DAM2fasta DB2fasta DB2quiva DBdust DBrm DBshow DBsplit DBstats fasta2DAM fasta2DB quiva2DB simulator ${VENV_TOFU}/bin
```


##### 6.  Test install and troubleshooting. There is a [known issue](https://github.com/PacificBiosciences/cDNA_primer/issues/3) post-installation

After TOFU installation, load the virtual enviornment and smrtShell
```bash
source ${VENV_TOFU}/bin/activate
/opt/smrtanalysis/smrtcmds/bin/smrtshell
```

Test if it works by trying:
```bash
collapse_isoforms_by_sam.py --help
```

If you get this error: `ImportError: No module named modified_bx_intervals.intersection_unique`

Go into the modified_bx_intervals directory
```bash
cd ${VENV_ENV}/lib/python2.7/site-packages/pbtools.pbtranscript-0.4.0-py2.7-linux-x86_64.egg/pbtools/pbtranscript/modified_bx_intervals/
```

See if you have these two files:
```bash
intersection_unique.so
__init__.py
```

If you are only missing the `__init__.py`, just create it in the folder using the touch command
```bash
touch __init__.py
```

Rerun the command with help to see if the error is gone

Worst case senario, you are missing 'intersection_unique.so'. Youwill need to reinstall pbtranscript-TOFU

Other installation issues see [github issues for cDNA_primer](https://github.com/PacificBiosciences/cDNA_primer/issues/17)
