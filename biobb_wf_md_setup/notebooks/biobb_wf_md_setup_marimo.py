# GENERATED FILE, do not edit by hand: edit biobb_wf_md_setup.ipynb (+ biobb_wf_md_setup_marimo_overlay.py) and regenerate
# this file with ipynb_to_marimo.py (Claude Code skill: ipynb-to-marimo).

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")

with app.setup:
    # Like Jupyter, save all output files in the notebook's folder
    import marimo as _mo
    import os as _os

    _os.chdir(_mo.notebook_dir())


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Protein MD Setup tutorial using BioExcel Building Blocks (biobb)
    **Based on the official GROMACS tutorial:** [http://www.mdtutorials.com/gmx/lysozyme/index.html](http://www.mdtutorials.com/gmx/lysozyme/index.html)
    ***
    This tutorial aims to illustrate the process of **setting up a simulation system** containing a **protein**, step by step, using the **BioExcel Building Blocks library (biobb)**. The particular example used is the **Lysozyme** protein (PDB code 1AKI, [https://doi.org/10.2210/pdb1AKI/pdb](https://doi.org/10.2210/pdb1AKI/pdb)).
    ***

    ## Settings

    ### Biobb modules used

    - [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
    - [biobb_model](https://github.com/bioexcel/biobb_model): Tools to model macromolecular structures.
    - [biobb_gromacs](https://github.com/bioexcel/biobb_gromacs): Tools to setup and run Molecular Dynamics simulations.
    - [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.

    ### Auxiliary libraries used

    * [marimo](https://marimo.io/): Reactive Python notebooks, stored as plain Python files.
    * [molviewspec](https://molstar.org/mol-view-spec/): Python toolkit to describe molecular scenes and display them with the [Mol*](https://molstar.org/) viewer in notebooks.
    * [plotly](https://plot.ly/python/offline/): Python interactive graphing library integrated in Jupyter notebooks.


    ### Conda Installation and Launch

    ```console
    git clone https://github.com/bioexcel/biobb_wf_md_setup.git
    cd biobb_wf_md_setup
    conda env create -f conda_env/environment.yml
    conda activate biobb_wf_md_setup
    marimo edit biobb_wf_md_setup/notebooks/biobb_wf_md_setup_marimo.py
    ```

    ***
    ## Pipeline steps
    1. [Input Parameters](#input)
    2. [Fetching PDB Structure](#fetch)
    3. [Fix Protein Structure](#fix)
    4. [Create Protein System Topology](#top)
    5. [Create Solvent Box](#box)
    6. [Fill the Box with Water Molecules](#water)
    7. [Adding Ions](#ions)
    8. [Energetically Minimize the System](#min)
    9. [Equilibrate the System (NVT)](#nvt)
    10. [Equilibrate the System (NPT)](#npt)
    11. [Free Molecular Dynamics Simulation](#free)
    12. [Post-processing and Visualizing Resulting 3D Trajectory](#post)
    13. [Output Files](#output)
    14. [Questions & Comments](#questions)

    ***
    <img src="https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png" alt="Bioexcel2 logo"
    	title="Bioexcel2 logo" width="400" />
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="input"></a>
    ## Input parameters
    **Input parameters** needed:
    - **pdbCode**: PDB code of the protein structure (e.g. 1AKI, [https://doi.org/10.2210/pdb1AKI/pdb](https://doi.org/10.2210/pdb1AKI/pdb))
    """)
    return


@app.cell
def _(mo):
    import re
    from pathlib import Path
    import molviewspec as mvs

    # Type a PDB code and click "Run setup": the whole pipeline runs, and it runs
    # again, step by step, every time a new code is submitted.
    pdb_form = mo.ui.text(value="1AKI", label="**PDB code**").form(
        submit_button_label="Run setup",
        validate=lambda code: None if re.fullmatch(r"\s*[0-9][A-Za-z0-9]{3}\s*", code or "")
        else "A PDB code has 4 characters, e.g. 1AKI",
    )
    pdb_form
    return Path, mvs, pdb_form


@app.cell
def _(mo, pdb_form):
    # Nothing below this cell runs until the form above is submitted
    mo.stop(pdb_form.value is None, mo.md("**Enter a PDB code and click _Run setup_ to start the pipeline.**"))
    pdbCode = pdb_form.value.strip().upper()
    mo.md(f"Setting up an MD simulation for **{pdbCode}**")
    return (pdbCode,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="fetch"></a>
    ***
    ## Fetching PDB structure
    Downloading **PDB structure** with the **protein molecule** from the RCSB PDB database.<br>
    Alternatively, a **PDB file** can be used as starting structure. <br>
    ***
    **Building Blocks** used:
    - [Pdb](https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb) from **biobb_io.api.pdb**
    ***
    """)
    return


@app.cell
def _(pdbCode):
    # Downloading desired PDB file 
    # Import module
    from biobb_io.api.pdb import pdb

    # Create properties dict and inputs/outputs
    downloaded_pdb = pdbCode+'.pdb'
    _prop = {
        'pdb_code': pdbCode
    }

    #Create and launch bb
    pdb(output_pdb_path=downloaded_pdb,
        properties=_prop)
    return (downloaded_pdb,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="vis3D"></a>
    ### Visualizing 3D structure
    Visualizing the downloaded/given **PDB structure** using **Mol\***:
    """)
    return


@app.cell
def _(Path, downloaded_pdb, mo, mvs):
    # Show protein
    _builder = mvs.create_builder()
    _structure = _builder.download(url=downloaded_pdb).parse(format='pdb').model_structure()
    _structure.component(selector='all').representation(type='ball_and_stick').color(custom={
        'molstar_color_theme_name': 'element-symbol',
        # Grey carbons (Mol* colours them by chain by default)
        'molstar_color_theme_params': {'carbonColor': {'name': 'element-symbol', 'params': {}}}
    })
    mo.iframe(_builder.molstar_html(data={downloaded_pdb: Path(downloaded_pdb).read_bytes()}), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="fix"></a>
    ***
    ## Fix protein structure
    **Checking** and **fixing** (if needed) the protein structure:<br>
    - **Modeling** **missing side-chain atoms**, modifying incorrect **amide assignments**, choosing **alternative locations**.<br>
    - **Checking** for missing **backbone atoms**, **heteroatoms**, **modified residues** and possible **atomic clashes**.

    ***
    **Building Blocks** used:
    - [FixSideChain](https://biobb-model.readthedocs.io/en/latest/model.html#module-model.fix_side_chain) from **biobb_model.model.fix_side_chain**
    ***
    """)
    return


@app.cell
def _(downloaded_pdb, pdbCode):
    # Check & Fix PDB
    # Import module
    from biobb_model.model.fix_side_chain import fix_side_chain

    # Create prop dict and inputs/outputs
    fixed_pdb = pdbCode + '_fixed.pdb'

    # Create and launch bb
    fix_side_chain(input_pdb_path=downloaded_pdb, 
                 output_pdb_path=fixed_pdb)
    return (fixed_pdb,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualizing 3D structure
    Visualizing the fixed **PDB structure** using **Mol\***. In this particular example, the checking step didn't find any issue to be solved, so there is no difference between the original structure and the fixed one.
    """)
    return


@app.cell
def _(Path, fixed_pdb, mo, mvs):
    # Show protein
    _builder = mvs.create_builder()
    _structure = _builder.download(url=fixed_pdb).parse(format='pdb').model_structure()
    _structure.component(selector='all').representation(type='ball_and_stick').color(custom={
        'molstar_color_theme_name': 'element-symbol',
        # Grey carbons (Mol* colours them by chain by default)
        'molstar_color_theme_params': {'carbonColor': {'name': 'element-symbol', 'params': {}}}
    })
    mo.iframe(_builder.molstar_html(data={fixed_pdb: Path(fixed_pdb).read_bytes()}), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="top"></a>
    ***
    ## Create protein system topology
    **Building GROMACS topology** corresponding to the protein structure.<br>
    Force field used in this tutorial is [**amber99sb-ildn**](https://dx.doi.org/10.1002%2Fprot.22711): AMBER **parm99** force field with **corrections on backbone** (sb) and **side-chain torsion potentials** (ildn). Water molecules type used in this tutorial is [**spc/e**](https://pubs.acs.org/doi/abs/10.1021/j100308a038).<br>
    Adding **hydrogen atoms** if missing. Automatically identifying **disulfide bridges**. <br>

    Generating two output files:
    - **GROMACS structure** (gro file)
    - **GROMACS topology** ZIP compressed file containing:
        - *GROMACS topology top file* (top file)
        - *GROMACS position restraint file/s* (itp file/s)
    ***
    **Building Blocks** used:
    - [Pdb2gmx](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.pdb2gmx) from **biobb_gromacs.gromacs.pdb2gmx**
    ***
    """)
    return


@app.cell
def _(fixed_pdb, pdbCode):
    # Create system topology
    # Import module
    from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx

    # Create inputs/outputs
    output_pdb2gmx_gro = pdbCode+'_pdb2gmx.gro'
    output_pdb2gmx_top_zip = pdbCode+'_pdb2gmx_top.zip'

    # Create and launch bb
    pdb2gmx(input_pdb_path=fixed_pdb, 
            output_gro_path=output_pdb2gmx_gro, 
            output_top_zip_path=output_pdb2gmx_top_zip)
    return output_pdb2gmx_gro, output_pdb2gmx_top_zip


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualizing 3D structure
    Visualizing the generated **GRO structure** using **Mol\***. Note that **hydrogen atoms** were added to the structure by the **pdb2gmx GROMACS tool** when generating the **topology**.
    """)
    return


@app.cell
def _(Path, mo, mvs, output_pdb2gmx_gro):
    # Show protein
    _builder = mvs.create_builder()
    _structure = _builder.download(url=output_pdb2gmx_gro).parse(format='gro').model_structure()
    _structure.component(selector='all').representation(type='ball_and_stick').color(custom={
        'molstar_color_theme_name': 'element-symbol',
        # Grey carbons (Mol* colours them by chain by default)
        'molstar_color_theme_params': {'carbonColor': {'name': 'element-symbol', 'params': {}}}
    })
    mo.iframe(_builder.molstar_html(data={output_pdb2gmx_gro: Path(output_pdb2gmx_gro).read_bytes()}), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="box"></a>
    ***
    ## Create solvent box
    Define the unit cell for the **protein structure MD system** to fill it with water molecules.<br>
    A **cubic box** is used to define the unit cell, with a **distance from the protein to the box edge of 1.0 nm**. The protein is **centered in the box**.

    ***
    **Building Blocks** used:
    - [Editconf](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.editconf) from **biobb_gromacs.gromacs.editconf**
    ***
    """)
    return


@app.cell
def _(output_pdb2gmx_gro, pdbCode):
    # Editconf: Create solvent box
    # Import module
    from biobb_gromacs.gromacs.editconf import editconf

    # Create prop dict and inputs/outputs
    output_editconf_gro = pdbCode+'_editconf.gro'

    _prop = {
        'box_type': 'cubic',
        'distance_to_molecule': 1.0
    }

    #Create and launch bb
    editconf(input_gro_path=output_pdb2gmx_gro, 
             output_gro_path=output_editconf_gro,
             properties=_prop)
    return (output_editconf_gro,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="water"></a>
    ***
    ## Fill the box with water molecules
    Fill the unit cell for the **protein structure system** with water molecules.<br>
    The solvent type used is the default **Simple Point Charge water (SPC)**, a generic equilibrated 3-point solvent model.

    ***
    **Building Blocks** used:
    - [Solvate](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.solvate) from **biobb_gromacs.gromacs.solvate**
    ***
    """)
    return


@app.cell
def _(output_editconf_gro, output_pdb2gmx_top_zip, pdbCode):
    # Solvate: Fill the box with water molecules
    from biobb_gromacs.gromacs.solvate import solvate

    # Create prop dict and inputs/outputs
    output_solvate_gro = pdbCode+'_solvate.gro'
    output_solvate_top_zip = pdbCode+'_solvate_top.zip'

    # Create and launch bb
    solvate(input_solute_gro_path=output_editconf_gro, 
            output_gro_path=output_solvate_gro, 
            input_top_zip_path=output_pdb2gmx_top_zip, 
            output_top_zip_path=output_solvate_top_zip)
    return output_solvate_gro, output_solvate_top_zip


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualizing 3D structure
    Visualizing the **protein system** with the newly added **solvent box** using **Mol\***.<br> Note the **cubic box** filled with **water molecules** surrounding the **protein structure**, which is **centered** right in the middle of the cube.
    """)
    return


@app.cell
def _(Path, mo, mvs, output_solvate_gro):
    # Show protein
    _builder = mvs.create_builder()
    _structure = _builder.download(url=output_solvate_gro).parse(format='gro').model_structure()
    _structure.component(selector='protein').representation(type='cartoon').color(color='green')
    _structure.component(selector=mvs.ComponentExpression(label_comp_id='SOL')).representation(type='ball_and_stick').color(custom={'molstar_color_theme_name': 'element-symbol'})
    mo.iframe(_builder.molstar_html(data={output_solvate_gro: Path(output_solvate_gro).read_bytes()}), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ions"></a>
    ***
    ## Adding ions
    Add ions to neutralize the **protein structure** charge
    - [Step 1](#ionsStep1): Creating portable binary run file for ion generation
    - [Step 2](#ionsStep2): Adding ions to **neutralize** the system
    ***
    **Building Blocks** used:
    - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp**
    - [Genion](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.genion) from **biobb_gromacs.gromacs.genion**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ionsStep1"></a>
    ### Step 1: Creating portable binary run file for ion generation
    A simple **energy minimization** molecular dynamics parameters (mdp) properties will be used to generate the portable binary run file for **ion generation**, although **any legitimate combination of parameters** could be used in this step.
    """)
    return


@app.cell
def _(output_solvate_gro, output_solvate_top_zip, pdbCode):
    # Grompp: Creating portable binary run file for ion generation
    from biobb_gromacs.gromacs.grompp import grompp

    # Create prop dict and inputs/outputs
    output_gppion_tpr = pdbCode+'_gppion.tpr'
    _prop = {
        'simulation_type': 'ions',
        'maxwarn': 1
    }

    # Create and launch bb
    grompp(input_gro_path=output_solvate_gro, 
           input_top_zip_path=output_solvate_top_zip, 
           output_tpr_path=output_gppion_tpr,  
           properties=_prop)
    return grompp, output_gppion_tpr


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ionsStep2"></a>
    ### Step 2: Adding ions to neutralize the system
    Replace **solvent molecules** with **ions** to **neutralize** the system.
    """)
    return


@app.cell
def _(output_gppion_tpr, output_solvate_top_zip, pdbCode):
    # Genion: Adding ions to neutralize the system
    from biobb_gromacs.gromacs.genion import genion

    # Create prop dict and inputs/outputs
    output_genion_gro = pdbCode+'_genion.gro'
    output_genion_top_zip = pdbCode+'_genion_top.zip'
    _prop={
        'neutral':True,
        'concentration':0
    }

    # Create and launch bb
    genion(input_tpr_path=output_gppion_tpr, 
           output_gro_path=output_genion_gro, 
           input_top_zip_path=output_solvate_top_zip, 
           output_top_zip_path=output_genion_top_zip, 
           properties=_prop)
    return output_genion_gro, output_genion_top_zip


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Visualizing 3D structure
    Visualizing the **neutralized protein system** with the newly added **ions** using **Mol\***
    """)
    return


@app.cell
def _(Path, mo, mvs, output_genion_gro):
    # Show protein
    _builder = mvs.create_builder()
    _structure = _builder.download(url=output_genion_gro).parse(format='gro').model_structure()
    _structure.component(selector='protein').representation(type='cartoon').color(custom={'molstar_color_theme_name': 'secondary-structure'})
    ions = [mvs.ComponentExpression(label_comp_id='NA'), mvs.ComponentExpression(label_comp_id='CL')]
    _structure.component(selector=ions).representation(type='spacefill').color(custom={'molstar_color_theme_name': 'element-symbol'})
    mo.iframe(_builder.molstar_html(data={output_genion_gro: Path(output_genion_gro).read_bytes()}), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="min"></a>
    ***
    ## Energetically minimize the system
    Energetically minimize the **protein system** till reaching a desired potential energy.
    - [Step 1](#emStep1): Creating portable binary run file for energy minimization
    - [Step 2](#emStep2): Energetically minimize the **system** till reaching a force of 500 kJ/mol*nm.
    - [Step 3](#emStep3): Checking **energy minimization** results. Plotting energy by time during the **minimization** process.
    ***
    **Building Blocks** used:
    - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp**
    - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun**
    - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="emStep1"></a>
    ### Step 1: Creating portable binary run file for energy minimization
    The **minimization** type of the **molecular dynamics parameters (mdp) property** contains the main default parameters to run an **energy minimization**:

    -  integrator  = steep ; Algorithm (steep = steepest descent minimization)
    -  emtol       = 1000.0 ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
    -  emstep      = 0.01 ; Minimization step size (nm)
    -  nsteps      = 50000 ; Maximum number of (minimization) steps to perform

    In this particular example, the method used to run the **energy minimization** is the default **steepest descent**, but the **maximum force** is placed at **500 kJ/mol\*nm^2**, and the **maximum number of steps** to perform (if the maximum force is not reached) to **5,000 steps**.
    """)
    return


@app.cell
def _(grompp, output_genion_gro, output_genion_top_zip, pdbCode):
    # Grompp: Creating portable binary run file for mdrun

    # Create prop dict and inputs/outputs
    output_gppmin_tpr = pdbCode+'_gppmin.tpr'
    _prop = {
        'mdp':{
            'emtol':'500',
            'nsteps':'5000'
        },
        'simulation_type': 'minimization'
    }

    # Create and launch bb
    grompp(input_gro_path=output_genion_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppmin_tpr,  
           properties=_prop)
    return (output_gppmin_tpr,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="emStep2"></a>
    ### Step 2: Running Energy Minimization
    Running **energy minimization** using the **tpr file** generated in the previous step.
    """)
    return


@app.cell
def _(output_gppmin_tpr, pdbCode):
    # Mdrun: Running minimization
    from biobb_gromacs.gromacs.mdrun import mdrun

    # Create prop dict and inputs/outputs
    output_min_trr = pdbCode+'_min.trr'
    output_min_gro = pdbCode+'_min.gro'
    output_min_edr = pdbCode+'_min.edr'
    output_min_log = pdbCode+'_min.log'

    # Create and launch bb
    mdrun(input_tpr_path=output_gppmin_tpr, 
          output_trr_path=output_min_trr, 
          output_gro_path=output_min_gro, 
          output_edr_path=output_min_edr, 
          output_log_path=output_min_log)
    return mdrun, output_min_edr, output_min_gro


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="emStep3"></a>
    ### Step 3: Checking Energy Minimization results
    Checking **energy minimization** results. Plotting **potential energy** by time during the minimization process.
    """)
    return


@app.cell
def _(output_min_edr, pdbCode):
    # GMXEnergy: Getting system energy by time  
    from biobb_analysis.gromacs.gmx_energy import gmx_energy

    # Create prop dict and inputs/outputs
    output_min_ene_xvg = pdbCode+'_min_ene.xvg'
    _prop = {
        'terms':  ["Potential"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_min_edr, 
              output_xvg_path=output_min_ene_xvg, 
              properties=_prop)
    return gmx_energy, output_min_ene_xvg


@app.cell
def _(output_min_ene_xvg):
    import plotly.graph_objs as go

    # Read data from file and filter energy values higher than 1000 kJ/mol
    with open(output_min_ene_xvg, 'r') as energy_file:
        _x, _y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in energy_file
            if not line.startswith(("#", "@"))
            if float(line.split()[1]) < 1000
        ])

    # Create a scatter plot
    _fig = go.Figure(data=go.Scatter(x=_x, y=_y, mode='lines'))

    # Update layout
    _fig.update_layout(title="Energy Minimization",
                      xaxis_title="Energy Minimization Step",
                      yaxis_title="Potential Energy kJ/mol",
                      height=600)

    # Show the figure
    _fig
    return (go,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="nvt"></a>
    ***
    ## Equilibrate the system (NVT)
    Equilibrate the **protein system** in **NVT ensemble** (constant Number of particles, Volume and Temperature). Protein **heavy atoms** will be restrained using position restraining forces: movement is permitted, but only after overcoming a substantial energy penalty. The utility of position restraints is that they allow us to equilibrate our solvent around our protein, without the added variable of structural changes in the protein.

    - [Step 1](#eqNVTStep1): Creating portable binary run file for system equilibration
    - [Step 2](#eqNVTStep2): Equilibrate the **protein system** with **NVT** ensemble.
    - [Step 3](#eqNVTStep3): Checking **NVT Equilibration** results. Plotting **system temperature** by time during the **NVT equilibration** process.
    ***
    **Building Blocks** used:
    - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp**
    - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun**
    - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNVTStep1"></a>
    ### Step 1: Creating portable binary run file for system equilibration (NVT)
    The **nvt** type of the **molecular dynamics parameters (mdp) property** contains the main default parameters to run an **NVT equilibration** with **protein restraints** (see [GROMACS mdp options](http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)):

    -  Define                   = -DPOSRES
    -  integrator               = md
    -  dt                       = 0.002
    -  nsteps                   = 5000
    -  pcoupl                   = no
    -  gen_vel                  = yes
    -  gen_temp                 = 300
    -  gen_seed                 = -1

    In this particular example, the default parameters will be used: **md** integrator algorithm, a **step size** of **2fs**, **5,000 equilibration steps** with the protein **heavy atoms restrained**, and a temperature of **300K**.

    *Please note that for the sake of time this tutorial is only running 10ps of NVT equilibration, whereas in the [original example](http://www.mdtutorials.com/gmx/lysozyme/06_equil.html) the simulated time was 100ps.*
    """)
    return


@app.cell
def _(grompp, output_genion_top_zip, output_min_gro, pdbCode):
    # Grompp: Creating portable binary run file for NVT Equilibration

    # Create prop dict and inputs/outputs
    output_gppnvt_tpr = pdbCode+'_gppnvt.tpr'
    _prop = {
        'mdp':{
            'nsteps': 5000,
            'dt': 0.002,
            'Define': '-DPOSRES',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        },
        'simulation_type': 'nvt'
    }

    # Create and launch bb
    grompp(input_gro_path=output_min_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppnvt_tpr,  
           properties=_prop)
    return (output_gppnvt_tpr,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNVTStep2"></a>
    ### Step 2: Running NVT equilibration
    """)
    return


@app.cell
def _(mdrun, output_gppnvt_tpr, pdbCode):
    # Mdrun: Running Equilibration NVT

    # Create prop dict and inputs/outputs
    output_nvt_trr = pdbCode+'_nvt.trr'
    output_nvt_gro = pdbCode+'_nvt.gro'
    output_nvt_edr = pdbCode+'_nvt.edr'
    output_nvt_log = pdbCode+'_nvt.log'
    output_nvt_cpt = pdbCode+'_nvt.cpt'

    # Create and launch bb
    mdrun(input_tpr_path=output_gppnvt_tpr, 
          output_trr_path=output_nvt_trr, 
          output_gro_path=output_nvt_gro, 
          output_edr_path=output_nvt_edr, 
          output_log_path=output_nvt_log, 
          output_cpt_path=output_nvt_cpt)
    return output_nvt_cpt, output_nvt_edr, output_nvt_gro


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNVTStep3"></a>
    ### Step 3: Checking NVT Equilibration results
    Checking **NVT Equilibration** results. Plotting **system temperature** by time during the NVT equilibration process.
    """)
    return


@app.cell
def _(gmx_energy, output_nvt_edr, pdbCode):
    # GMXEnergy: Getting system temperature by time during NVT Equilibration  

    # Create prop dict and inputs/outputs
    output_nvt_temp_xvg = pdbCode+'_nvt_temp.xvg'
    _prop = {
        'terms':  ["Temperature"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_nvt_edr, 
              output_xvg_path=output_nvt_temp_xvg, 
              properties=_prop)
    return (output_nvt_temp_xvg,)


@app.cell
def _(go, output_nvt_temp_xvg):
    # Read temperature data from file
    with open(output_nvt_temp_xvg, 'r') as temperature_file:
        _x, _y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in temperature_file
            if not line.startswith(("#", "@"))
        ])

    # Create a scatter plot
    _fig = go.Figure(data=go.Scatter(x=_x, y=_y, mode='lines+markers'))

    # Update layout
    _fig.update_layout(title="Temperature during NVT Equilibration",
                      xaxis_title="Time (ps)",
                      yaxis_title="Temperature (K)",
                      height=600)

    # Show the figure
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="npt"></a>
    ***
    ## Equilibrate the system (NPT)
    Equilibrate the **protein system** in **NPT** ensemble (constant Number of particles, Pressure and Temperature).
    - [Step 1](#eqNPTStep1): Creating portable binary run file for system equilibration
    - [Step 2](#eqNPTStep2): Equilibrate the **protein system** with **NPT** ensemble.
    - [Step 3](#eqNPTStep3): Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NPT equilibration** process.
    ***
    **Building Blocks** used:
    - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp**
    - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun**
    - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNPTStep1"></a>
    ### Step 1: Creating portable binary run file for system equilibration (NPT)

    The **npt** type of the **molecular dynamics parameters (mdp) property** contains the main default parameters to run an **NPT equilibration** with **protein restraints** (see [GROMACS mdp options](http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)):

    -  Define                   = -DPOSRES
    -  integrator               = md
    -  dt                       = 0.002
    -  nsteps                   = 5000
    -  pcoupl = Parrinello-Rahman
    -  pcoupltype = isotropic
    -  tau_p = 1.0
    -  ref_p = 1.0
    -  compressibility = 4.5e-5
    -  refcoord_scaling = com
    -  gen_vel = no

    In this particular example, the default parameters will be used: **md** integrator algorithm, a **time step** of **2fs**, **5,000 equilibration steps** with the protein **heavy atoms restrained**, and a Parrinello-Rahman **pressure coupling** algorithm.

    *Please note that for the sake of time this tutorial is only running 10ps of NPT equilibration, whereas in the [original example](http://www.mdtutorials.com/gmx/lysozyme/07_equil2.html) the simulated time was 100ps.*
    """)
    return


@app.cell
def _(grompp, output_genion_top_zip, output_nvt_cpt, output_nvt_gro, pdbCode):
    # Grompp: Creating portable binary run file for NPT System Equilibration

    # Create prop dict and inputs/outputs
    output_gppnpt_tpr = pdbCode+'_gppnpt.tpr'
    _prop = {
        'mdp':{
            'nsteps':'5000',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        },
        'simulation_type': 'npt'
    }

    # Create and launch bb
    grompp(input_gro_path=output_nvt_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppnpt_tpr, 
           input_cpt_path=output_nvt_cpt,  
           properties=_prop)
    return (output_gppnpt_tpr,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNPTStep2"></a>
    ### Step 2: Running NPT equilibration
    """)
    return


@app.cell
def _(mdrun, output_gppnpt_tpr, pdbCode):
    # Mdrun: Running NPT System Equilibration

    # Create prop dict and inputs/outputs
    output_npt_trr = pdbCode+'_npt.trr'
    output_npt_gro = pdbCode+'_npt.gro'
    output_npt_edr = pdbCode+'_npt.edr'
    output_npt_log = pdbCode+'_npt.log'
    output_npt_cpt = pdbCode+'_npt.cpt'

    # Create and launch bb
    mdrun(input_tpr_path=output_gppnpt_tpr, 
          output_trr_path=output_npt_trr, 
          output_gro_path=output_npt_gro, 
          output_edr_path=output_npt_edr, 
          output_log_path=output_npt_log, 
          output_cpt_path=output_npt_cpt)
    return output_npt_cpt, output_npt_edr, output_npt_gro


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="eqNPTStep3"></a>
    ### Step 3: Checking NPT Equilibration results
    Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NPT equilibration** process.
    """)
    return


@app.cell
def _(gmx_energy, output_npt_edr, pdbCode):
    # GMXEnergy: Getting system pressure and density by time during NPT Equilibration  

    # Create prop dict and inputs/outputs
    output_npt_pd_xvg = pdbCode+'_npt_PD.xvg'
    _prop = {
        'terms':  ["Pressure","Density"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_npt_edr, 
              output_xvg_path=output_npt_pd_xvg, 
              properties=_prop)
    return (output_npt_pd_xvg,)


@app.cell
def _(go, output_npt_pd_xvg):
    from plotly.subplots import make_subplots

    # Read pressure and density data from file 
    with open(output_npt_pd_xvg,'r') as pd_file:
        _x, _y, z = zip(*[
            (float(line.split()[0]), float(line.split()[1]), float(line.split()[2]))
            for line in pd_file
            if not line.startswith(("#", "@"))
        ])

    trace1 = go.Scatter(
        x=_x,y=_y
    )
    trace2 = go.Scatter(
        x=_x,y=z
    )

    _fig = make_subplots(rows=1, cols=2, print_grid=False)
    _fig.append_trace(trace1, 1, 1)
    _fig.append_trace(trace2, 1, 2)

    _fig.update_layout(
        height=500,
        title='Pressure and Density during NPT Equilibration',
        showlegend=False,
        xaxis1_title='Time (ps)',
        yaxis1_title='Pressure (bar)',
        xaxis2_title='Time (ps)',
        yaxis2_title='Density (Kg*m^-3)'
    )

    # Show the figure
    _fig
    return (make_subplots,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Explore the energy files
    Choose a **simulation step** and an **energy term**: the term is extracted from that step's
    energy file with **GMXEnergy** and plotted right away. Change either selector and the plot updates.
    """)
    return


@app.cell
def _(mo, output_min_edr, output_npt_edr, output_nvt_edr):
    energy_step = mo.ui.dropdown(
        options={
            "Energy minimization": output_min_edr,
            "NVT equilibration": output_nvt_edr,
            "NPT equilibration": output_npt_edr,
        },
        value="NPT equilibration",
        label="Simulation step",
    )
    energy_term = mo.ui.dropdown(
        options=["Potential", "Kinetic-En.", "Total-Energy", "Temperature", "Pressure", "Density", "Volume"],
        value="Temperature",
        label="Energy term",
    )
    mo.hstack([energy_step, energy_term], justify="start", gap=2)
    return energy_step, energy_term


@app.cell
def _(Path, energy_step, energy_term, gmx_energy, go, mo, pdbCode):
    # Extract the chosen term (about a second) and plot it
    _xvg = pdbCode + '_explore_energy.xvg'
    Path(_xvg).unlink(missing_ok=True)
    gmx_energy(input_energy_path=energy_step.value,
               output_xvg_path=_xvg,
               properties={'terms': [energy_term.value], 'disable_logs': True})
    mo.stop(not Path(_xvg).exists(),
            mo.callout(mo.md(f"**{energy_term.value}** is not in the {energy_step.selected_key} energy file."), kind="warn"))

    with open(_xvg, 'r') as _xvg_file:
        _x, _y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in _xvg_file
            if not line.startswith(("#", "@"))
        ])

    _fig = go.Figure(data=go.Scatter(x=_x, y=_y, mode='lines'))
    _fig.update_layout(title=f"{energy_term.value} during {energy_step.selected_key}",
                       xaxis_title="Step" if energy_step.selected_key == "Energy minimization" else "Time (ps)",
                       yaxis_title=energy_term.value,
                       height=500)
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="free"></a>
    ***
    ## Free Molecular Dynamics Simulation
    Upon completion of the **two equilibration phases (NVT and NPT)**, the system is now well-equilibrated at the desired temperature and pressure. The **position restraints** can now be released. The last step of the **protein** MD setup is a short, **free MD simulation**, to ensure the robustness of the system.
    - [Step 1](#mdStep1): Creating portable binary run file to run a **free MD simulation**.
    - [Step 2](#mdStep2): Run short MD simulation of the **protein system**.
    - [Step 3](#mdStep3): Checking results for the final step of the setup process, the **free MD run**. Plotting **Root Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time during the **free MD run** step.
    ***
    **Building Blocks** used:
    - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp**
    - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun**
    - [GMXRms](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rms) from **biobb_analysis.gromacs.gmx_rms**
    - [GMXRgyr](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rgyr) from **biobb_analysis.gromacs.gmx_rgyr**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="mdStep1"></a>
    ### Step 1: Creating portable binary run file to run a free MD simulation

    The **free** type of the **molecular dynamics parameters (mdp) property** contains the main default parameters to run an **free MD simulation** (see [GROMACS mdp options](http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html)):

    -  integrator               = md
    -  dt                       = 0.002 (ps)
    -  nsteps                   = 50000

    In this particular example, the default parameters will be used: **md** integrator algorithm, a **time step** of **2fs**, and a total of **50,000 md steps** (100ps).

    *Please note that for the sake of time this tutorial is only running 100ps of free MD, whereas in the [original example](http://www.mdtutorials.com/gmx/lysozyme/08_MD.html) the simulated time was 1ns (1000ps).*
    """)
    return


@app.cell
def _(grompp, output_genion_top_zip, output_npt_cpt, output_npt_gro, pdbCode):
    # Grompp: Creating portable binary run file for mdrun

    # Create prop dict and inputs/outputs
    output_gppmd_tpr = pdbCode+'_gppmd.tpr'
    _prop = {
        'mdp':{
            'nsteps':'50000',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        },
        'simulation_type': 'free'
    }

    # Create and launch bb
    grompp(input_gro_path=output_npt_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppmd_tpr, 
           input_cpt_path=output_npt_cpt, 
           properties=_prop)
    return (output_gppmd_tpr,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="mdStep2"></a>
    ### Step 2: Running short free MD simulation
    """)
    return


@app.cell
def _(mdrun, output_gppmd_tpr, pdbCode):
    # Mdrun: Running free dynamics

    # Create prop dict and inputs/outputs
    output_md_trr = pdbCode+'_md.trr'
    output_md_gro = pdbCode+'_md.gro'
    output_md_edr = pdbCode+'_md.edr'
    output_md_log = pdbCode+'_md.log'
    output_md_cpt = pdbCode+'_md.cpt'

    # Create and launch bb
    mdrun(input_tpr_path=output_gppmd_tpr, 
          output_trr_path=output_md_trr, 
          output_gro_path=output_md_gro, 
          output_edr_path=output_md_edr, 
          output_log_path=output_md_log, 
          output_cpt_path=output_md_cpt)
    return output_md_gro, output_md_trr


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="mdStep3"></a>
    ### Step 3: Checking free MD simulation results
    Checking results for the final step of the setup process, the **free MD run**. Plotting **Root Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time during the **free MD run** step. **RMSd** against the **experimental structure** (input structure of the pipeline) and against the **minimized and equilibrated structure** (output structure of the NPT equilibration step).
    """)
    return


@app.cell
def _(output_gppmd_tpr, output_md_trr, pdbCode):
    # GMXRms: Computing Root Mean Square deviation to analyse structural stability 
    #         RMSd against minimized and equilibrated snapshot (backbone atoms)   

    from biobb_analysis.gromacs.gmx_rms import gmx_rms

    # Create prop dict and inputs/outputs
    output_rms_first = pdbCode+'_rms_first.xvg'
    _prop = {
        'selection':  'Backbone',
        #'selection': 'non-Water'
    }

    # Create and launch bb
    gmx_rms(input_structure_path=output_gppmd_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rms_first, 
              properties=_prop)
    return gmx_rms, output_rms_first


@app.cell
def _(gmx_rms, output_gppmin_tpr, output_md_trr, pdbCode):
    # GMXRms: Computing Root Mean Square deviation to analyse structural stability 
    #         RMSd against experimental structure (backbone atoms)   

    # Create prop dict and inputs/outputs
    output_rms_exp = pdbCode+'_rms_exp.xvg'
    _prop = {
        'selection':  'Backbone',
        #'selection': 'non-Water'
    }

    # Create and launch bb
    gmx_rms(input_structure_path=output_gppmin_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rms_exp, 
              properties=_prop)
    return (output_rms_exp,)


@app.cell
def _(go, make_subplots, output_rms_exp, output_rms_first):
    # Read RMS vs first snapshot data from file 
    with open(output_rms_first,'r') as rms_first_file:
        _x, _y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in rms_first_file
            if not line.startswith(("#", "@"))
        ])

    # Read RMS vs experimental structure data from file 
    with open(output_rms_exp,'r') as rms_exp_file:
        x2, y2 = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in rms_exp_file
            if not line.startswith(("#", "@"))
        ])

    _fig = make_subplots()
    _fig.add_trace(go.Scatter(x=_x, y=_y, mode="lines+markers", name="RMSd vs first"))
    _fig.add_trace(go.Scatter(x=_x, y=y2, mode="lines+markers", name="RMSd vs exp"))

    # Set layout including height
    _fig.update_layout(
        title="RMSd during free MD Simulation",
        xaxis=dict(title="Time (ps)"),
        yaxis=dict(title="RMSd (nm)"),
        height=600
    )

    # Show the figure
    _fig
    return


@app.cell
def _(output_gppmin_tpr, output_md_trr, pdbCode):
    # GMXRgyr: Computing Radius of Gyration to measure the protein compactness during the free MD simulation 

    from biobb_analysis.gromacs.gmx_rgyr import gmx_rgyr

    # Create prop dict and inputs/outputs
    output_rgyr = pdbCode+'_rgyr.xvg'
    _prop = {
        'selection':  'Backbone'
    }

    # Create and launch bb
    gmx_rgyr(input_structure_path=output_gppmin_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rgyr, 
              properties=_prop)
    return (output_rgyr,)


@app.cell
def _(go, output_rgyr):
    # Read Rgyr data from file
    with open(output_rgyr, 'r') as rgyr_file:
        _x, _y = zip(*[
            (float(line.split()[0]), float(line.split()[1]))
            for line in rgyr_file
            if not line.startswith(("#", "@"))
        ])

    # Create a scatter plot
    _fig = go.Figure(data=go.Scatter(x=_x, y=_y, mode='lines+markers'))

    # Update layout
    _fig.update_layout(title="Radius of Gyration",
                      xaxis_title="Time (ps)",
                      yaxis_title="Rgyr (nm)",
                      height=600)

    # Show the figure
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="post"></a>
    ***
    ## Post-processing and Visualizing resulting 3D trajectory
    Post-processing and Visualizing the **protein system** MD setup **resulting trajectory** using **Mol\***
    - [Step 1](#ppStep1): *Imaging* the resulting trajectory, **stripping out water molecules and ions** and **correcting periodicity issues**.
    - [Step 2](#ppStep2): Generating a *dry* structure, **removing water molecules and ions** from the final snapshot of the MD setup pipeline.
    - [Step 3](#ppStep3): Visualizing the *imaged* trajectory using the *dry* structure as a **topology**.
    ***
    **Building Blocks** used:
    - [GMXImage](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_image) from **biobb_analysis.gromacs.gmx_image**
    - [GMXTrjConvStr](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str) from **biobb_analysis.gromacs.gmx_trjconv_str**
    ***
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ppStep1"></a>
    ### Step 1: *Imaging* the resulting trajectory.
    Stripping out **water molecules and ions** and **correcting periodicity issues**
    """)
    return


@app.cell
def _(output_gppmd_tpr, output_md_trr, pdbCode):
    # GMXImage: "Imaging" the resulting trajectory
    #           Removing water molecules and ions from the resulting structure
    from biobb_analysis.gromacs.gmx_image import gmx_image

    # Create prop dict and inputs/outputs
    output_imaged_traj = pdbCode+'_imaged_traj.trr'
    _prop = {
        'center_selection':  'Protein',
        'output_selection': 'Protein',
        'pbc' : 'mol',
        'center' : True
    }

    # Create and launch bb
    gmx_image(input_traj_path=output_md_trr,
             input_top_path=output_gppmd_tpr,
             output_traj_path=output_imaged_traj, 
              properties=_prop)
    return (output_imaged_traj,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ppStep2"></a>
    ### Step 2: Generating the output *dry* structure.
    **Removing water molecules and ions** from the resulting structure
    """)
    return


@app.cell
def _(output_gppmd_tpr, output_md_gro, pdbCode):
    # GMXTrjConvStr: Converting and/or manipulating a structure
    #                Removing water molecules and ions from the resulting structure
    #                The "dry" structure will be used as a topology to visualize 
    #                the "imaged dry" trajectory generated in the previous step.
    from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str

    # Create prop dict and inputs/outputs
    output_dry_gro = pdbCode+'_md_dry.gro'
    _prop = {
        'selection':  'Protein'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_md_gro,
             input_top_path=output_gppmd_tpr,
             output_str_path=output_dry_gro, 
              properties=_prop)
    return (output_dry_gro,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="ppStep3"></a>
    ### Step 3: Visualizing the generated dehydrated trajectory.
    Using the **imaged trajectory** (output of the [Post-processing step 1](#ppStep1)) with the **dry structure** (output of the [Post-processing step 2](#ppStep2)) as a topology.
    The trajectory **plays automatically in a loop**. Use the **stop** button at the top left of the viewer to pause it and step through individual frames.
    """)
    return


@app.cell
def _(mo):
    traj_style = mo.ui.dropdown(
        options={"Cartoon": "cartoon", "Ball and stick": "ball_and_stick", "Spacefill": "spacefill", "Surface": "surface"},
        value="Cartoon",
        label="Representation",
    )
    traj_color = mo.ui.dropdown(
        options={
            "Secondary structure": "secondary-structure",
            "Residue index (rainbow)": "sequence-id",
            "Residue name": "residue-name",
            "Hydrophobicity": "hydrophobicity",
            "Element": "element-symbol",
        },
        value="Secondary structure",
        label="Colour",
    )
    # debounce: the viewer is rebuilt only when the slider is released, not at every step of a drag
    traj_speed = mo.ui.slider(start=20, stop=500, step=20, value=100, label="ms per frame", debounce=True)
    mo.hstack([traj_style, traj_color, traj_speed], justify="start", gap=2)
    return traj_color, traj_speed, traj_style


@app.cell
def _(
    Path,
    mo,
    mvs,
    output_dry_gro,
    output_imaged_traj,
    traj_color,
    traj_speed,
    traj_style,
):
    # Show trajectory (redrawn whenever one of the controls above changes)
    import subprocess

    # Number of frames of the imaged trajectory (from GROMACS gmx check)
    check = subprocess.run(['gmx', 'check', '-f', output_imaged_traj], capture_output=True, text=True)
    n_frames = int(next(line.split()[1] for line in check.stderr.splitlines() if line.startswith('Coords')))

    _builder = mvs.create_builder()
    # Imaged trajectory coordinates
    _builder.download(url=output_imaged_traj).parse(format='trr').coordinates(ref='trajectory')
    # Dry structure used as topology
    _structure = _builder.download(url=output_dry_gro).parse(format='gro').model_structure(coordinates_ref='trajectory', ref='structure')
    color = {'molstar_color_theme_name': traj_color.value}
    if traj_color.value == 'element-symbol':
        # Grey carbons (Mol* colours them by chain by default)
        color['molstar_color_theme_params'] = {'carbonColor': {'name': 'element-symbol', 'params': {}}}
    _structure.component(selector='protein').representation(type=traj_style.value).color(custom=color)
    # Play the trajectory in a loop
    _builder.animation(autoplay=True, loop=True).interpolate(
        kind='scalar', target_ref='structure', property='model_index',
        start=0, end=n_frames - 1, discrete=True, duration_ms=n_frames * traj_speed.value
    )
    data = {
        output_dry_gro: Path(output_dry_gro).read_bytes(),
        output_imaged_traj: Path(output_imaged_traj).read_bytes()
    }
    mo.iframe(_builder.molstar_html(data=data), height="600px")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    <a id="output"></a>
    ## Output files

    Important **Output files** generated:
    - **output_md_gro** (1aki_md.gro): **Final structure** (snapshot) of the MD setup protocol.
    - **output_md_trr** (1aki_md.trr): **Final trajectory** of the MD setup protocol.
    - **output_md_cpt** (1aki_md.cpt): **Final checkpoint file**, with information about the state of the simulation. It can be used to **restart** or **continue** a MD simulation.
    - **output_gppmd_tpr** (1aki_gppmd.tpr): **Final tpr file**, GROMACS portable binary run input file. This file contains the starting structure of the **MD setup free MD simulation step**, together with the molecular topology and all the simulation parameters. It can be used to **extend** the simulation.
    - **output_genion_top_zip** (1aki_genion_top.zip): **Final topology** of the MD system. It is a compressed zip file including a **topology file** (.top) and a set of auxiliary **include topology** files (.itp).

    **Analysis** (MD setup check) output files generated:
    - **output_rms_first** (1aki_rms_first.xvg): **Root Mean Square deviation (RMSd)** against **minimized and equilibrated structure** of the final **free MD run step**.
    - **output_rms_exp** (1aki_rms_exp.xvg): **Root Mean Square deviation (RMSd)** against **experimental structure** of the final **free MD run step**.
    - **output_rgyr** (1aki_rgyr.xvg): **Radius of Gyration** of the final **free MD run step** of the **setup pipeline**.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ***
    <a id="questions"></a>

    ## Questions & Comments

    Questions, issues, suggestions and comments are really welcome!

    * GitHub issues:
        * [https://github.com/bioexcel/biobb](https://github.com/bioexcel/biobb)

    * BioExcel forum:
        * [https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library](https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library)
    """)
    return


if __name__ == "__main__":
    app.run()
