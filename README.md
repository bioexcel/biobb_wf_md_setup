# Protein MD Setup tutorial using BioExcel Building Blocks (biobb)

Based on the official [Gromacs tutorial](http://www.mdtutorials.com/gmx/lysozyme/index.html).

This tutorial aims to illustrate the process of **setting up a simulation** system containing a **protein**, step by step, using the **BioExcel Building Blocks library (biobb)**. The particular example used is the **Lysozyme** protein (PDB code 1AKI).

### Biobb modules used

* [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
* [biobb_model](https://github.com/bioexcel/biobb_model): Tools to model macromolecular structures.
* [biobb_md](https://github.com/bioexcel/biobb_md): Tools to setup and run Molecular Dynamics simulations.
* [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.

### Auxiliar libraries used

* [nglview](http://nglviewer.org/#nglview): Jupyter/IPython widget to interactively view molecular structures and trajectories in notebooks.
* [ipywidgets](https://github.com/jupyter-widgets/ipywidgets): Interactive HTML widgets for Jupyter notebooks and the IPython kernel.
* [plotly](https://plot.ly/python/offline/): Python interactive graphing library integrated in Jupyter notebooks.

### Conda Installation

```console
conda install -c bioconda biobb_MD_setup //// NOTE: this is not yet available ////
```

## Pipeline steps

1. Input Parameters
2. Fetching PDB Structure
3. Fix Protein Structure
4. Create Protein System Topology
5. Create Solvent Box
6. Fill the Box with Water Molecules
7. Adding Ions
8. Energetically Minimize the System
9. Equilibrate the System (NVT)
11. Equilibrate the System (NPT)
12. Free Molecular Dynamics Simulation
13. Post-processing and Visualizing Resulting 3D Trajectory
14. Output Files

### Version
June 2019 Release

### Copyright & Licensing
This software has been developed in the [MMB group](http://mmb.irbbarcelona.org) at the [BSC](http://www.bsc.es/) & [IRB](https://www.irbbarcelona.org/) for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU H2020 [823830](http://cordis.europa.eu/projects/823830), EU H2020 [675728](http://cordis.europa.eu/projects/675728)).

* (c) 2015-2019 [Barcelona Supercomputing Center](https://www.bsc.es/)
* (c) 2015-2019 [Institute for Research in Biomedicine](https://www.irbbarcelona.org/)

Licensed under the
[Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0), see the file LICENSE for details.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")
