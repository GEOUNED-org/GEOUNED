# Installation

Recommended method of installing is to make use of Conda / Mamba.
This results in installation of GEOUNED with your system Python and the most straightforward usage.

There are other options which are also documented.
The main difference between the installation options is the manner in which the main dependency FreeCAD is installed and this has consequences for how the code is used.

With some options (Conda) offering integration of FreeCAD Python API into your system Python and other options requiring GEOUNED scripts to be run with freecadcmd or freecad.cmd and integration of GEOUNED and FreeCAD performed with system path appending or freecad.pip 

## Linux (Ubuntu)

### Conda


In principle, installing any Conda/Mamba distribution will work. A few Conda/Mamba options are:
- [Miniforge](https://github.com/conda-forge/miniforge)
- [Anaconda](https://www.anaconda.com/download)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

In this example we will install Miniforge (which includes Mamba)
```bash
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Create a new environment, I've chosen Python 3.10 here but newer versions are
also supported.
```bash
mamba create --name new_env python=3.10 -y
```

Activate the environment
```bash
mamba activate new_env
```

Then you can install the cad_to_dagmc package
```bash
mamba install -y -c conda-forge geouned
```

### AppImage

### Apt-get

### Snap



## Mac OS

### Conda

## Windows

### Conda

### Portable FreeCAD

