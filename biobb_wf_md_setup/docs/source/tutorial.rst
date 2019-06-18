
Protein MD Setup tutorial using BioExcel Building Blocks (biobb)
================================================================

Based on the official Gromacs tutorial: http://www.mdtutorials.com/gmx/lysozyme/index.html
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

--------------

This tutorial aims to illustrate the process of **setting up a
simulation system** containing a **protein**, step by step, using the
**BioExcel Building Blocks library (biobb)**. The particular example
used is the **Lysozyme** protein (PDB code 1AKI). \*Biobb modules*\*
used:

-  `biobb_io <https://github.com/bioexcel/biobb_io>`__: Tools to fetch
   biomolecular data from public databases.
-  `biobb_model <https://github.com/bioexcel/biobb_model>`__: Tools to
   model macromolecular structures.
-  `biobb_md <https://github.com/bioexcel/biobb_md>`__: Tools to setup
   and run Molecular Dynamics simulations.
-  `biobb_analysis <https://github.com/bioexcel/biobb_analysis>`__:
   Tools to analyse Molecular Dynamics trajectories.

**Auxiliar libraries** used:

-  `nglview <http://nglviewer.org/#nglview>`__: Jupyter/IPython widget
   to interactively view molecular structures and trajectories in
   notebooks.
-  `ipywidgets <https://github.com/jupyter-widgets/ipywidgets>`__:
   Interactive HTML widgets for Jupyter notebooks and the IPython
   kernel.
-  `plotly <https://plot.ly/python/offline/>`__: Python interactive
   graphing library integrated in Jupyter notebooks.

Conda **Installation**:

-  **conda install -c bioconda biobb_MD_setup** //// *NOTE: this is not
   yet available* ////

--------------

Pipeline steps:
~~~~~~~~~~~~~~~

1.  `Input Parameters <#input>`__
2.  `Fetching PDB Structure <#fetch>`__
3.  `Fix Protein Structure <#fix>`__
4.  `Create Protein System Topology <#top>`__
5.  `Create Solvent Box <#box>`__
6.  `Fill the Box with Water Molecules <#water>`__
7.  `Adding Ions <#ions>`__
8.  `Energetically Minimize the System <#min>`__
9.  `Equilibrate the System (NVT) <#nvt>`__
10. `Equilibrate the System (NPT) <#npt>`__
11. `Free Molecular Dynamics Simulation <#free>`__
12. `Post-processing and Visualizing Resulting 3D Trajectory <#post>`__
13. `Output Files <#output>`__

--------------

\**\*

\*\ **## Input parameters**\ Input parameters\ **needed: -**\ pdbCode**:
PDB code of the protein structure (e.g. 1AKI)

.. parsed-literal::

    import nglview
    import ipywidgets
    
    pdbCode = "1AKI"

**## Fetching PDB structure Downloading**\ *\ PDB structure\ *\ **with
the**\ *\ protein molecule\ *\ **from the RCSB PDB database.
Alternatively, a**\ *\ PDB file\ *\ **can be used as starting structure.
** **Building Blocks** used: -
`Pdb <https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb>`__
from **biobb_io.api.pdb** \**\*

.. parsed-literal::

    # Downloading desired PDB file 
    # Import module
    from biobb_io.api.pdb import Pdb
    
    # Create properties dict and inputs/outputs
    downloaded_pdb = pdbCode+'.pdb'
    prop = {
        'pdb_code': pdbCode
    }
    
    #Create and launch bb
    Pdb(output_pdb_path=downloaded_pdb,
        properties=prop).launch()

### Visualizing 3D structure Visualizing the downloaded/given **PDB
structure** using **NGL**:

.. parsed-literal::

    # Show protein
    view = nglview.show_file(downloaded_pdb)
    view.add_representation(repr_type='ball+stick', selection='all')
    view._remote_call('setSize', target='Widget', args=['','600px'])
    view

.. parsed-literal::

    view.render_image()
    view.download_image(filename='image1.png')

.. parsed-literal::

    from IPython import display

.. parsed-literal::

    display.HTML("<img src='_static/image1.png'></img>")




.. raw:: html

    <img src='_static/image1.png'></img>



\*\ **## Fix protein structure**\ Checking\ **and**\ fixing\ **(if
needed) the protein structure: -**\ Modelingmissing side-chain
atoms\ **, modifying incorrect**\ amide assignments\ **,
choosing**\ alternative locations\ **. -**\ Checking\ **for
missing**\ backbone atoms\ **,**\ heteroatoms\ **,**\ modified
residues\ **and possible**\ atomic clashes**.

--------------

**Building Blocks** used: -
`FixSideChain <https://biobb-model.readthedocs.io/en/latest/model.html#module-model.fix_side_chain>`__
from **biobb_model.model.fix_side_chain** \**\*

.. parsed-literal::

    # Check & Fix PDB
    # Import module
    from biobb_model.model.fix_side_chain import FixSideChain
    
    # Create prop dict and inputs/outputs
    fixed_pdb = pdbCode + '_fixed.pdb'
    
    # Create and launch bb
    FixSideChain(input_pdb_path=downloaded_pdb, 
                 output_pdb_path=fixed_pdb).launch()

Visualizing 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~

Visualizing the fixed **PDB structure** using **NGL**. In this
particular example, the checking step didn’t find any issue to be
solved, so there is no difference between the original structure and the
fixed one.

.. parsed-literal::

    # Show protein
    view = nglview.show_file(fixed_pdb)
    view.add_representation(repr_type='ball+stick', selection='all')
    view._remote_call('setSize', target='Widget', args=['','600px'])
    view.camera='orthographic'
    view

.. parsed-literal::

    view.render_image()
    view.download_image(filename='image2.png')

.. parsed-literal::

    display.HTML("<img src='_static/image2.png'></img>")




.. raw:: html

    <img src='_static/image2.png'></img>



\*\ **## Create protein system topology**\ Building GROMACS
topology\ **corresponding to the protein structure. Force field used in
this tutorial
is**\ `amber99sb-ildn <https://dx.doi.org/10.1002%2Fprot.22711>`__\ **:
AMBER**\ parm99\ **force field with**\ corrections on backbone\ **(sb)
and**\ side-chain torsion potentials\ **(ildn). Water molecules type
used in this tutorial
is**\ `spc/e <https://pubs.acs.org/doi/abs/10.1021/j100308a038>`__\ **.
Adding**\ hydrogen atoms\ **if missing. Automatically
identifying**\ disulfide bridges**.

Generating two output files: - **GROMACS structure** (gro file) -
**GROMACS topology** ZIP compressed file containing: - *GROMACS topology
top file* (top file) - *GROMACS position restraint file/s* (itp file/s)
*Building Blocks\ *\ **used:
-**\ *\ *\ `Pdb2gmx <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.pdb2gmx>`__\ *\ *\ **from**\ *\ biobb_md.gromacs.pdb2gmx*

.. parsed-literal::

    # Create system topology
    # Import module
    from biobb_md.gromacs.pdb2gmx import Pdb2gmx
    
    # Create inputs/outputs
    output_pdb2gmx_gro = pdbCode+'_pdb2gmx.gro'
    output_pdb2gmx_top_zip = pdbCode+'_pdb2gmx_top.zip'
    
    # Create and launch bb
    Pdb2gmx(input_pdb_path=fixed_pdb, 
            output_gro_path=output_pdb2gmx_gro, 
            output_top_zip_path=output_pdb2gmx_top_zip).launch()

Visualizing 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~

Visualizing the generated **GRO structure** using **NGL**. Note that
**hydrogen atoms** were added to the structure by the **pdb2gmx GROMACS
tool** when generating the **topology**.

.. parsed-literal::

    # Show protein
    view = nglview.show_file(output_pdb2gmx_gro)
    view.add_representation(repr_type='ball+stick', selection='all')
    view._remote_call('setSize', target='Widget', args=['','600px'])
    view.camera='orthographic'
    view

.. parsed-literal::

    view.render_image()
    view.download_image(filename='image3.png')

.. parsed-literal::

    display.HTML("<img src='_static/image3.png'></img>")




.. raw:: html

    <img src='_static/image3.png'></img>



\*\ **## Create solvent box Define the unit cell for the**\ protein
structure MD system\ **to fill it with water molecules. A**\ cubic
box\ **is used to define the unit cell, with a**\ distance from the
protein to the box edge of 1.0 nm\ **. The protein is**\ centered in the
box**.

--------------

**Building Blocks** used: -
`Editconf <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.editconf>`__
from **biobb_md.gromacs.editconf** \**\*

.. parsed-literal::

    # Editconf: Create solvent box
    # Import module
    from biobb_md.gromacs.editconf import Editconf
    
    # Create prop dict and inputs/outputs
    output_editconf_gro = pdbCode+'_editconf.gro'
    
    prop = {
        'box_type': 'cubic',
        'distance_to_molecule': 1.0
    }
    
    #Create and launch bb
    Editconf(input_gro_path=output_pdb2gmx_gro, 
             output_gro_path=output_editconf_gro,
             properties=prop).launch()

\*\ **## Fill the box with water molecules Fill the unit cell for
the**\ protein structure system\ **with water molecules. The solvent
type used is the default**\ Simple Point Charge water (SPC)**, a generic
equilibrated 3-point solvent model.

--------------

**Building Blocks** used: -
`Solvate <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.solvate>`__
from **biobb_md.gromacs.solvate** \**\*

.. parsed-literal::

    # Solvate: Fill the box with water molecules
    from biobb_md.gromacs.solvate import Solvate
    
    # Create prop dict and inputs/outputs
    output_solvate_gro = pdbCode+'_solvate.gro'
    output_solvate_top_zip = pdbCode+'_solvate_top.zip'
    
    # Create and launch bb
    Solvate(input_solute_gro_path=output_editconf_gro, 
            output_gro_path=output_solvate_gro, 
            input_top_zip_path=output_pdb2gmx_top_zip, 
            output_top_zip_path=output_solvate_top_zip).launch()

Visualizing 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~

Visualizing the **protein system** with the newly added **solvent box**
using **NGL**. Note the **cubic box** filled with **water molecules**
surrounding the **protein structure**, which is **centered** right in
the middle of the cube.

.. parsed-literal::

    # Show protein
    view = nglview.show_file(output_solvate_gro)
    view.clear_representations()
    view.add_representation(repr_type='cartoon', selection='solute', color='green')
    view.add_representation(repr_type='ball+stick', selection='SOL')
    view._remote_call('setSize', target='Widget', args=['','600px'])
    view.camera='orthographic'
    view

.. parsed-literal::

    view.render_image()
    view.download_image(filename='image4.png')

.. parsed-literal::

    display.HTML("<img src='_static/image4.png'></img>")




.. raw:: html

    <img src='_static/image4.png'></img>



**## Adding ions Add ions to neutralize the**\ *\ protein
structure\ *\ **charge -**\ *\ *\ `Step 1 <#ionsStep1>`__\ *\ *\ **:
Creating portable binary run file for ion generation -**\ *\ *\ `Step
2 <#ionsStep2>`__\ *\ *\ **: Adding ions to**\ *\ neutralize\ *\ **the
system** **Building Blocks** used: -
`Grompp <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp>`__
from **biobb_md.gromacs.grompp** -
`Genion <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.genion>`__
from **biobb_md.gromacs.genion** \**\*

### Step 1: Creating portable binary run file for ion generation A
simple **energy minimization** molecular dynamics parameters (mdp)
properties will be used to generate the portable binary run file for
**ion generation**, although **any legitimate combination of
parameters** could be used in this step.

.. parsed-literal::

    # Grompp: Creating portable binary run file for ion generation
    from biobb_md.gromacs.grompp import Grompp
    
    # Create prop dict and inputs/outputs
    output_gppion_tpr = pdbCode+'_gppion.tpr'
    prop = {
        'mdp':{
            'type': 'minimization'
        }
    }
    
    # Create and launch bb
    Grompp(input_gro_path=output_solvate_gro, 
           input_top_zip_path=output_solvate_top_zip, 
           output_tpr_path=output_gppion_tpr,  
           properties=prop).launch()

### Step 2: Adding ions to neutralize the system Replace **solvent
molecules** with **ions** to **neutralize** the system.

.. parsed-literal::

    # Genion: Adding ions to neutralize the system
    from biobb_md.gromacs.genion import Genion
    
    # Create prop dict and inputs/outputs
    output_genion_gro = pdbCode+'_genion.gro'
    output_genion_top_zip = pdbCode+'_genion_top.zip'
    prop={
        'neutral':True
    }
    
    # Create and launch bb
    Genion(input_tpr_path=output_gppion_tpr, 
           output_gro_path=output_genion_gro, 
           input_top_zip_path=output_solvate_top_zip, 
           output_top_zip_path=output_genion_top_zip, 
           properties=prop).launch()

Visualizing 3D structure
~~~~~~~~~~~~~~~~~~~~~~~~

Visualizing the **neutralized protein system** with the newly added
**ions** using **NGL**

.. parsed-literal::

    # Show protein
    view = nglview.show_file(output_genion_gro)
    view.clear_representations()
    view.add_representation(repr_type='cartoon', selection='solute', color='sstruc')
    view.add_representation(repr_type='ball+stick', selection='NA')
    view.add_representation(repr_type='ball+stick', selection='CL')
    view._remote_call('setSize', target='Widget', args=['','600px'])
    view.camera='orthographic'
    view

.. parsed-literal::

    view.render_image()
    view.download_image(filename='image5.png')

.. parsed-literal::

    display.HTML("<img src='_static/image5.png'></img>")




.. raw:: html

    <img src='_static/image5.png'></img>



**## Energetically minimize the system Energetically minimize
the**\ *\ protein system\ *\ **till reaching a desired potential energy.
-**\ *\ *\ `Step 1 <#emStep1>`__\ *\ *\ **: Creating portable binary run
file for energy minimization -**\ *\ *\ `Step 2 <#emStep2>`__\ *\ *\ **:
Energetically minimize the**\ *\ system\ *\ **till reaching a force of
500 kJ mol-1 nm-1. -**\ *\ *\ `Step 3 <#emStep3>`__\ *\ *\ **:
Checking**\ *\ energy minimization\ *\ **results. Plotting energy by
time during the**\ *\ minimization\ *\ **process.** **Building Blocks**
used: -
`Grompp <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp>`__
from **biobb_md.gromacs.grompp** -
`Mdrun <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun>`__
from **biobb_md.gromacs.mdrun** -
`GMXEnergy <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy>`__
from **biobb_analysis.gromacs.gmx_energy** \**\*

### Step 1: Creating portable binary run file for energy minimization
The **minimization** type of the **molecular dynamics parameters (mdp)
property** contains the main default parameters to run an **energy
minimization**:

-  integrator = steep ; Algorithm (steep = steepest descent
   minimization)
-  emtol = 1000.0 ; Stop minimization when the maximum force < 1000.0
   kJ/mol/nm
-  emstep = 0.01 ; Minimization step size (nm)
-  nsteps = 50000 ; Maximum number of (minimization) steps to perform

In this particular example, the method used to run the **energy
minimization** is the default **steepest descent**, but the **maximum
force** is placed at **500 KJ/mol*nm^2**, and the **maximum number of
steps** to perform (if the maximum force is not reached) to **5,000
steps**.

.. parsed-literal::

    # Grompp: Creating portable binary run file for mdrun
    from biobb_md.gromacs.grompp import Grompp
    
    # Create prop dict and inputs/outputs
    output_gppmin_tpr = pdbCode+'_gppmin.tpr'
    prop = {
        'mdp':{
            'type': 'minimization',
            'emtol':'500',
            'nsteps':'5000'
        }
    }
    
    # Create and launch bb
    Grompp(input_gro_path=output_genion_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppmin_tpr,  
           properties=prop).launch()

### Step 2: Running Energy Minimization Running **energy minimization**
using the **tpr file** generated in the previous step.

.. parsed-literal::

    # Mdrun: Running minimization
    from biobb_md.gromacs.mdrun import Mdrun
    
    # Create prop dict and inputs/outputs
    output_min_trr = pdbCode+'_min.trr'
    output_min_gro = pdbCode+'_min.gro'
    output_min_edr = pdbCode+'_min.edr'
    output_min_log = pdbCode+'_min.log'
    
    # Create and launch bb
    Mdrun(input_tpr_path=output_gppmin_tpr, 
          output_trr_path=output_min_trr, 
          output_gro_path=output_min_gro, 
          output_edr_path=output_min_edr, 
          output_log_path=output_min_log).launch()

### Step 3: Checking Energy Minimization results Checking **energy
minimization** results. Plotting **potential energy** by time during the
minimization process.

.. parsed-literal::

    # GMXEnergy: Getting system energy by time  
    from biobb_analysis.gromacs.gmx_energy import GMXEnergy
    
    # Create prop dict and inputs/outputs
    output_min_ene_xvg = pdbCode+'_min_ene.xvg'
    prop = {
        'terms':  ["Potential"]
    }
    
    # Create and launch bb
    GMXEnergy(input_energy_path=output_min_edr, 
              output_xvg_path=output_min_ene_xvg, 
              properties=prop).launch()

.. parsed-literal::

    import plotly
    import plotly.graph_objs as go
    
    #Read data from file and filter energy values higher than 1000 Kj/mol^-1
    with open(output_min_ene_xvg,'r') as energy_file:
        x,y = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in energy_file 
                if not line.startswith(("#","@")) 
                if float(line.split()[1]) < 1000 
            ])
        )
    
    plotly.offline.init_notebook_mode(connected=True)
    
    plotly.offline.iplot({
        "data": [go.Scatter(x=x, y=y)],
        "layout": go.Layout(title="Energy Minimization",
                            xaxis=dict(title = "Energy Minimization Step"),
                            yaxis=dict(title = "Potential Energy KJ/mol-1")
                           )
    })



.. raw:: html

            <script type="text/javascript">
            window.PlotlyConfig = {MathJaxConfig: 'local'};
            if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
            if (typeof require !== 'undefined') {
            require.undef("plotly");
            requirejs.config({
                paths: {
                    'plotly': ['https://cdn.plot.ly/plotly-latest.min']
                }
            });
            require(['plotly'], function(Plotly) {
                window._Plotly = Plotly;
            });
            }
            </script>
            



.. raw:: html

    <div>
            
            
                <div id="64a55da9-6e3f-43f2-b541-0c0c57afb87e" class="plotly-graph-div" style="height:525px; width:100%;"></div>
                <script type="text/javascript">
                    require(["plotly"], function(Plotly) {
                        window.PLOTLYENV=window.PLOTLYENV || {};
                        window.PLOTLYENV.BASE_URL='https://plot.ly';
                        
                    if (document.getElementById("64a55da9-6e3f-43f2-b541-0c0c57afb87e")) {
                        Plotly.newPlot(
                            '64a55da9-6e3f-43f2-b541-0c0c57afb87e',
                            [{"type": "scatter", "uid": "ebf29f53-7894-4099-a70a-1f2494c559dc", "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 15.0, 16.0, 17.0, 18.0, 19.0, 21.0, 22.0, 23.0, 25.0, 26.0, 27.0, 28.0, 29.0, 31.0, 33.0, 34.0, 35.0, 36.0, 38.0, 39.0, 40.0, 42.0, 43.0, 44.0, 45.0, 46.0, 48.0, 49.0, 50.0, 52.0, 53.0, 54.0, 55.0, 56.0, 58.0, 59.0, 60.0, 62.0, 63.0, 64.0, 65.0, 67.0, 68.0, 69.0, 71.0, 72.0, 73.0, 74.0, 75.0, 77.0, 78.0, 79.0, 81.0, 82.0, 83.0, 84.0, 86.0, 87.0, 88.0, 90.0, 91.0, 92.0, 93.0, 95.0, 96.0, 97.0, 99.0, 100.0, 101.0, 102.0, 104.0, 105.0, 106.0, 108.0, 109.0, 110.0, 111.0, 112.0, 114.0, 115.0, 116.0, 118.0, 119.0, 120.0, 121.0, 123.0, 124.0, 125.0, 127.0, 128.0, 129.0, 130.0, 132.0, 133.0, 134.0, 136.0, 137.0, 138.0, 139.0, 140.0, 141.0, 143.0, 144.0, 145.0, 147.0, 148.0, 149.0, 150.0, 152.0, 153.0, 154.0, 156.0, 157.0, 158.0, 159.0, 160.0, 162.0, 163.0, 164.0, 166.0, 167.0, 168.0, 169.0, 171.0, 172.0, 173.0, 175.0, 176.0, 177.0, 178.0, 180.0, 181.0, 182.0, 184.0, 185.0, 186.0, 187.0, 189.0, 190.0, 191.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 200.0, 202.0, 203.0, 204.0, 205.0, 207.0, 208.0, 209.0, 211.0, 212.0, 213.0, 215.0, 216.0, 217.0, 218.0, 219.0, 221.0, 222.0, 223.0, 225.0, 226.0, 227.0, 228.0, 229.0, 231.0, 232.0, 233.0, 235.0, 236.0, 237.0, 238.0, 239.0, 241.0, 242.0, 243.0, 245.0, 246.0, 247.0, 248.0, 249.0, 251.0, 252.0, 253.0, 255.0, 256.0, 257.0, 258.0, 259.0, 261.0, 262.0, 263.0, 265.0, 266.0, 267.0, 268.0, 269.0, 271.0, 272.0, 273.0, 275.0, 276.0, 277.0, 278.0, 280.0, 281.0, 282.0, 284.0, 285.0, 286.0, 287.0, 289.0, 290.0, 291.0, 293.0, 294.0, 295.0, 296.0, 297.0, 299.0, 300.0, 301.0, 303.0, 304.0, 305.0, 306.0, 307.0, 309.0, 310.0, 311.0, 313.0, 314.0, 315.0, 316.0, 317.0, 319.0, 320.0, 321.0, 323.0, 324.0, 325.0, 326.0, 328.0, 329.0, 330.0, 332.0, 333.0, 334.0, 335.0, 337.0, 338.0, 339.0, 341.0, 342.0, 343.0, 344.0, 345.0, 347.0, 348.0, 349.0, 351.0, 352.0, 353.0, 354.0, 355.0, 357.0, 358.0, 359.0, 361.0, 362.0, 363.0, 364.0, 366.0, 367.0, 368.0, 370.0, 371.0, 372.0, 373.0, 375.0, 376.0, 377.0, 379.0, 380.0, 381.0, 382.0, 383.0, 384.0, 386.0, 387.0, 388.0, 390.0, 391.0, 392.0, 393.0, 395.0, 396.0, 397.0, 399.0, 400.0, 401.0, 402.0, 403.0, 405.0, 406.0, 407.0, 409.0, 410.0, 411.0, 412.0, 414.0, 415.0, 416.0, 418.0, 419.0, 420.0, 421.0, 422.0, 423.0, 425.0, 427.0, 428.0, 429.0, 430.0, 432.0, 433.0, 434.0, 436.0, 437.0, 438.0, 440.0, 441.0, 442.0, 443.0, 444.0, 446.0, 447.0, 448.0, 450.0, 451.0, 452.0, 453.0, 454.0, 455.0, 456.0, 458.0, 459.0, 460.0, 462.0, 463.0, 464.0, 465.0, 467.0, 468.0, 469.0, 471.0, 472.0, 473.0, 474.0, 475.0, 477.0, 478.0, 479.0, 481.0, 482.0, 483.0, 484.0, 486.0, 487.0, 488.0, 490.0, 491.0, 492.0, 493.0, 495.0, 496.0, 497.0, 499.0, 500.0, 501.0, 502.0, 503.0, 505.0, 506.0, 507.0, 509.0, 510.0, 511.0, 512.0, 513.0, 515.0, 516.0, 517.0, 519.0, 520.0, 521.0, 522.0, 524.0, 525.0, 526.0, 528.0, 529.0, 530.0, 531.0, 533.0, 534.0, 535.0, 537.0, 538.0, 539.0, 540.0, 541.0, 543.0, 544.0, 545.0, 547.0, 548.0, 549.0, 550.0, 551.0, 553.0, 554.0, 555.0, 557.0, 558.0, 559.0, 560.0, 561.0, 563.0, 564.0, 565.0, 567.0, 568.0, 569.0, 570.0, 572.0, 573.0, 574.0, 576.0, 577.0, 578.0, 579.0, 581.0, 582.0, 583.0, 585.0, 586.0, 587.0, 588.0, 589.0, 591.0, 592.0, 593.0, 595.0, 596.0, 597.0, 598.0, 599.0, 601.0, 602.0, 603.0, 605.0, 606.0, 607.0, 608.0, 610.0, 611.0, 612.0, 614.0, 615.0, 616.0, 617.0, 619.0, 620.0, 621.0, 623.0, 624.0, 625.0, 626.0, 627.0, 629.0, 630.0, 631.0, 633.0, 634.0, 635.0, 636.0, 637.0, 639.0, 640.0, 641.0, 643.0, 644.0, 645.0, 646.0, 648.0, 649.0, 650.0, 652.0, 653.0, 654.0, 655.0, 657.0, 658.0, 659.0, 661.0, 662.0, 663.0, 664.0, 665.0, 666.0, 668.0, 669.0, 670.0, 672.0, 673.0, 675.0, 676.0, 677.0, 679.0, 680.0, 681.0, 682.0, 683.0, 685.0, 686.0, 687.0, 689.0, 690.0, 691.0, 692.0, 693.0, 694.0, 695.0, 697.0, 699.0, 700.0, 701.0, 702.0, 704.0, 705.0, 706.0, 708.0, 709.0, 710.0, 711.0, 712.0, 714.0, 715.0, 716.0, 718.0, 719.0, 720.0, 721.0, 722.0, 724.0, 725.0, 726.0, 728.0, 729.0, 730.0, 731.0, 732.0, 733.0, 735.0, 736.0, 737.0, 739.0, 740.0, 741.0, 742.0, 744.0, 745.0, 746.0, 748.0, 749.0, 750.0, 751.0, 753.0, 754.0, 755.0, 757.0, 758.0, 759.0, 760.0, 761.0, 763.0, 764.0, 765.0, 767.0, 768.0, 769.0, 770.0, 772.0, 773.0, 774.0, 776.0, 777.0, 778.0, 779.0, 781.0, 782.0, 783.0, 785.0, 786.0, 787.0, 788.0, 789.0, 790.0, 792.0, 793.0, 794.0, 796.0, 797.0, 799.0, 800.0, 801.0, 803.0, 804.0, 805.0, 806.0, 807.0, 808.0, 809.0, 811.0, 812.0, 813.0, 815.0, 816.0, 817.0, 818.0, 819.0, 821.0, 822.0, 823.0, 825.0, 826.0, 827.0, 828.0, 830.0, 831.0, 832.0, 834.0, 835.0, 836.0, 837.0, 839.0, 840.0, 841.0, 843.0, 844.0, 845.0, 846.0, 847.0, 849.0, 850.0, 851.0, 853.0, 854.0, 855.0, 856.0, 858.0, 859.0, 860.0, 862.0, 863.0, 864.0, 865.0, 867.0, 868.0, 869.0, 871.0, 872.0, 873.0, 874.0, 875.0, 877.0, 878.0, 879.0, 881.0, 882.0, 883.0, 884.0, 885.0, 886.0, 888.0, 889.0, 890.0, 892.0, 893.0, 895.0, 896.0, 897.0, 899.0, 900.0, 901.0, 902.0, 903.0, 905.0, 906.0, 907.0, 909.0, 910.0, 911.0, 912.0, 913.0, 914.0, 915.0, 917.0, 919.0, 920.0, 921.0, 922.0, 924.0, 925.0, 926.0, 928.0, 929.0, 930.0, 931.0, 932.0, 934.0, 935.0, 936.0, 938.0, 939.0, 940.0, 941.0, 942.0, 944.0, 945.0, 946.0, 948.0, 949.0, 950.0, 951.0, 952.0, 953.0, 955.0, 956.0, 957.0, 959.0, 960.0, 961.0, 962.0, 964.0, 965.0, 966.0, 968.0, 969.0, 970.0, 971.0, 973.0, 974.0, 975.0, 977.0, 978.0, 979.0, 980.0, 981.0, 983.0, 984.0, 985.0, 987.0, 988.0, 989.0, 990.0, 991.0, 993.0, 994.0, 995.0, 997.0, 998.0, 999.0, 1000.0, 1002.0, 1003.0, 1004.0, 1006.0, 1007.0, 1008.0, 1009.0, 1011.0, 1012.0, 1013.0, 1015.0, 1016.0, 1017.0, 1018.0, 1019.0, 1021.0, 1022.0, 1023.0, 1025.0, 1026.0, 1027.0, 1028.0, 1030.0, 1031.0, 1032.0, 1034.0, 1035.0, 1036.0, 1037.0, 1039.0, 1040.0, 1041.0, 1043.0, 1044.0, 1045.0, 1046.0, 1047.0, 1048.0, 1049.0, 1051.0, 1053.0, 1054.0, 1055.0, 1056.0, 1058.0, 1059.0, 1060.0, 1062.0, 1063.0, 1064.0, 1065.0, 1066.0, 1068.0, 1069.0, 1070.0, 1072.0, 1073.0, 1074.0, 1075.0, 1076.0, 1078.0, 1079.0, 1080.0, 1082.0, 1083.0, 1084.0, 1085.0, 1086.0, 1087.0, 1089.0, 1091.0, 1092.0, 1093.0, 1094.0, 1096.0, 1097.0, 1098.0, 1100.0, 1101.0, 1102.0, 1103.0, 1104.0, 1106.0, 1107.0, 1108.0, 1110.0, 1111.0, 1112.0, 1113.0, 1114.0, 1115.0, 1116.0, 1118.0, 1120.0, 1121.0, 1122.0, 1123.0, 1125.0, 1126.0, 1127.0, 1129.0, 1130.0, 1131.0, 1133.0, 1134.0, 1135.0, 1136.0, 1137.0, 1139.0, 1140.0, 1141.0, 1143.0, 1144.0, 1145.0, 1146.0, 1147.0, 1148.0, 1149.0, 1151.0, 1152.0, 1153.0, 1155.0, 1156.0, 1158.0, 1159.0, 1160.0, 1162.0, 1163.0, 1164.0, 1165.0, 1166.0, 1168.0, 1169.0, 1170.0, 1172.0, 1173.0, 1174.0, 1175.0, 1176.0, 1178.0, 1179.0, 1180.0, 1182.0, 1183.0, 1184.0, 1185.0, 1186.0, 1187.0, 1188.0, 1190.0, 1192.0, 1193.0, 1194.0, 1195.0, 1197.0, 1198.0, 1199.0, 1201.0, 1202.0, 1203.0, 1204.0, 1205.0, 1207.0, 1208.0, 1209.0, 1211.0, 1212.0, 1213.0, 1214.0, 1215.0, 1217.0, 1218.0, 1219.0, 1221.0, 1222.0, 1223.0, 1224.0, 1226.0, 1227.0, 1228.0, 1230.0, 1231.0, 1232.0, 1233.0, 1234.0, 1235.0, 1237.0, 1238.0, 1239.0, 1241.0, 1242.0, 1244.0, 1245.0, 1246.0, 1248.0, 1249.0, 1250.0, 1251.0, 1252.0, 1254.0, 1255.0, 1256.0, 1258.0, 1259.0, 1260.0, 1261.0, 1262.0, 1264.0, 1265.0, 1266.0, 1268.0, 1269.0, 1270.0, 1271.0, 1272.0, 1274.0, 1275.0, 1276.0, 1278.0, 1279.0, 1280.0, 1281.0, 1282.0, 1283.0, 1285.0, 1286.0, 1287.0, 1289.0, 1290.0, 1292.0, 1293.0, 1294.0, 1296.0, 1297.0, 1298.0, 1299.0, 1300.0, 1302.0, 1303.0, 1304.0, 1306.0, 1307.0, 1308.0, 1309.0, 1310.0, 1312.0, 1313.0, 1314.0, 1316.0, 1317.0, 1318.0, 1319.0, 1320.0, 1322.0, 1323.0, 1324.0, 1326.0, 1327.0, 1328.0, 1329.0, 1330.0, 1331.0, 1333.0, 1334.0, 1335.0, 1337.0, 1338.0, 1340.0, 1341.0, 1342.0, 1344.0, 1345.0, 1346.0, 1347.0, 1348.0, 1350.0, 1351.0, 1352.0, 1354.0, 1355.0, 1356.0, 1357.0, 1358.0, 1360.0, 1361.0, 1362.0, 1364.0, 1365.0, 1366.0, 1367.0, 1368.0, 1370.0, 1371.0, 1372.0, 1374.0, 1375.0, 1376.0, 1377.0, 1378.0, 1379.0, 1381.0, 1382.0, 1383.0, 1385.0, 1386.0, 1388.0, 1389.0, 1390.0, 1392.0, 1393.0, 1394.0, 1395.0, 1396.0, 1398.0, 1399.0, 1400.0, 1402.0, 1403.0, 1404.0, 1405.0, 1406.0, 1408.0, 1409.0, 1410.0, 1412.0, 1413.0, 1414.0, 1415.0, 1416.0, 1418.0, 1419.0, 1420.0, 1422.0, 1423.0, 1424.0, 1425.0, 1426.0, 1427.0, 1429.0, 1430.0, 1431.0, 1433.0, 1434.0, 1436.0, 1437.0, 1438.0, 1440.0, 1441.0, 1442.0, 1443.0, 1444.0, 1446.0, 1447.0, 1448.0, 1450.0, 1451.0, 1452.0, 1453.0, 1454.0, 1456.0, 1457.0, 1458.0, 1460.0, 1461.0, 1462.0, 1463.0, 1464.0, 1466.0, 1467.0, 1468.0, 1470.0, 1471.0, 1472.0, 1473.0, 1474.0, 1475.0, 1477.0, 1478.0, 1479.0, 1481.0, 1482.0, 1484.0, 1485.0, 1486.0, 1488.0, 1489.0, 1490.0, 1491.0, 1492.0, 1494.0, 1495.0, 1496.0, 1498.0, 1499.0, 1500.0, 1501.0, 1502.0, 1504.0, 1505.0, 1506.0, 1508.0, 1509.0, 1510.0, 1511.0, 1512.0, 1514.0, 1515.0, 1516.0, 1518.0, 1519.0, 1520.0, 1521.0, 1522.0, 1523.0, 1525.0, 1526.0, 1527.0, 1529.0, 1530.0, 1532.0, 1533.0, 1534.0, 1536.0, 1537.0, 1538.0, 1539.0, 1540.0, 1542.0, 1543.0, 1544.0, 1546.0, 1547.0, 1548.0, 1549.0, 1550.0, 1552.0, 1553.0, 1554.0, 1556.0, 1557.0, 1558.0, 1559.0, 1560.0, 1561.0, 1563.0, 1564.0, 1565.0, 1567.0, 1568.0, 1570.0, 1571.0, 1572.0, 1574.0, 1575.0, 1576.0, 1577.0, 1578.0, 1580.0, 1581.0, 1582.0, 1584.0, 1585.0, 1586.0, 1587.0, 1588.0, 1590.0, 1591.0, 1592.0, 1594.0, 1595.0, 1596.0, 1597.0, 1598.0, 1600.0, 1601.0, 1602.0, 1604.0, 1605.0, 1606.0, 1607.0, 1608.0, 1609.0, 1611.0, 1612.0, 1613.0, 1615.0, 1616.0, 1618.0, 1619.0, 1620.0, 1622.0, 1623.0, 1624.0, 1625.0, 1626.0, 1628.0, 1629.0, 1630.0, 1632.0, 1633.0, 1634.0, 1635.0, 1636.0, 1638.0, 1639.0, 1640.0, 1642.0, 1643.0, 1644.0, 1645.0, 1646.0, 1648.0, 1649.0, 1650.0, 1652.0, 1653.0, 1654.0, 1655.0, 1656.0, 1657.0, 1659.0, 1660.0, 1661.0, 1663.0, 1664.0, 1666.0, 1667.0, 1668.0, 1670.0, 1671.0, 1672.0, 1673.0, 1674.0, 1676.0, 1677.0, 1678.0, 1680.0, 1681.0, 1682.0, 1683.0, 1684.0, 1686.0, 1687.0, 1688.0, 1690.0, 1691.0, 1692.0, 1693.0, 1694.0, 1696.0, 1697.0, 1698.0, 1700.0, 1701.0, 1702.0, 1703.0, 1704.0, 1705.0, 1707.0, 1708.0, 1709.0, 1711.0, 1712.0, 1714.0, 1715.0, 1716.0, 1718.0, 1719.0, 1720.0, 1721.0, 1722.0, 1724.0, 1725.0, 1726.0, 1728.0, 1729.0, 1730.0, 1731.0, 1732.0, 1734.0, 1735.0, 1736.0, 1738.0, 1739.0, 1740.0, 1741.0, 1742.0, 1744.0, 1745.0, 1746.0, 1748.0, 1749.0, 1750.0, 1751.0, 1752.0, 1753.0, 1755.0, 1757.0, 1758.0, 1759.0, 1760.0, 1762.0, 1763.0, 1764.0, 1766.0, 1767.0, 1768.0, 1770.0, 1771.0, 1772.0, 1773.0, 1774.0, 1775.0, 1776.0, 1778.0, 1779.0, 1780.0, 1782.0, 1783.0, 1784.0, 1785.0, 1787.0, 1788.0, 1789.0, 1791.0, 1792.0, 1793.0, 1794.0, 1795.0, 1796.0, 1798.0, 1799.0, 1800.0, 1802.0, 1803.0, 1805.0, 1806.0, 1807.0, 1809.0, 1810.0, 1811.0, 1812.0, 1813.0, 1815.0, 1816.0, 1817.0, 1819.0, 1820.0, 1821.0, 1822.0, 1823.0, 1825.0, 1826.0, 1827.0, 1829.0, 1830.0, 1831.0, 1832.0, 1833.0, 1834.0, 1836.0, 1837.0, 1838.0, 1840.0, 1841.0, 1843.0, 1844.0, 1845.0, 1847.0, 1848.0, 1849.0, 1850.0, 1851.0, 1853.0, 1854.0, 1855.0, 1857.0, 1858.0, 1859.0, 1860.0, 1861.0, 1863.0, 1864.0, 1865.0, 1867.0, 1868.0, 1869.0, 1870.0, 1871.0, 1872.0, 1873.0, 1875.0, 1877.0, 1878.0, 1879.0, 1880.0, 1882.0, 1883.0, 1884.0, 1886.0, 1887.0, 1888.0, 1889.0, 1890.0, 1892.0, 1893.0, 1894.0, 1896.0, 1897.0, 1898.0, 1899.0, 1900.0, 1902.0, 1903.0, 1904.0, 1906.0, 1907.0, 1908.0, 1909.0, 1910.0, 1911.0, 1913.0, 1914.0, 1915.0, 1917.0, 1918.0, 1919.0, 1920.0, 1922.0, 1923.0, 1924.0, 1926.0, 1927.0, 1928.0, 1929.0, 1931.0, 1932.0, 1933.0, 1935.0, 1936.0, 1937.0, 1938.0, 1939.0, 1941.0, 1942.0, 1943.0, 1945.0, 1946.0, 1947.0, 1948.0, 1950.0, 1951.0, 1952.0, 1954.0, 1955.0, 1956.0, 1957.0, 1959.0, 1960.0, 1961.0, 1963.0, 1964.0, 1965.0, 1966.0, 1967.0, 1968.0, 1970.0, 1971.0, 1972.0, 1974.0, 1975.0, 1977.0, 1978.0, 1979.0, 1981.0, 1982.0, 1983.0, 1984.0, 1985.0, 1987.0, 1988.0, 1989.0, 1991.0, 1992.0, 1993.0, 1994.0, 1995.0, 1997.0, 1998.0, 1999.0, 2001.0, 2002.0, 2003.0, 2004.0, 2005.0, 2006.0, 2007.0, 2009.0, 2010.0, 2011.0, 2013.0, 2014.0, 2015.0, 2016.0, 2018.0, 2019.0, 2020.0, 2022.0, 2023.0, 2024.0, 2025.0, 2027.0, 2028.0, 2029.0, 2031.0, 2032.0, 2033.0, 2034.0, 2035.0, 2037.0, 2038.0, 2039.0, 2041.0, 2042.0, 2043.0, 2044.0, 2046.0, 2047.0, 2048.0, 2050.0, 2051.0, 2052.0, 2053.0, 2055.0, 2056.0, 2057.0, 2059.0, 2060.0, 2061.0, 2062.0, 2063.0, 2064.0, 2066.0, 2067.0, 2068.0, 2070.0, 2071.0, 2072.0, 2073.0, 2075.0, 2076.0, 2077.0, 2079.0, 2080.0, 2081.0, 2082.0, 2083.0, 2085.0, 2086.0, 2087.0, 2089.0, 2090.0, 2091.0, 2092.0, 2094.0, 2095.0, 2096.0, 2098.0, 2099.0, 2100.0, 2101.0, 2103.0, 2104.0, 2105.0, 2107.0, 2108.0, 2109.0, 2110.0, 2111.0, 2113.0, 2114.0, 2115.0, 2117.0, 2118.0, 2119.0, 2120.0, 2121.0, 2122.0, 2124.0, 2125.0, 2126.0, 2128.0, 2129.0, 2131.0, 2132.0, 2133.0, 2135.0, 2136.0, 2137.0, 2138.0, 2139.0, 2141.0, 2142.0, 2143.0, 2145.0, 2146.0, 2147.0, 2148.0, 2149.0, 2150.0, 2152.0, 2153.0, 2154.0, 2156.0, 2157.0, 2158.0, 2159.0, 2161.0, 2162.0, 2163.0, 2165.0, 2166.0, 2167.0, 2168.0, 2169.0, 2171.0, 2172.0, 2173.0, 2175.0, 2176.0, 2177.0, 2178.0, 2179.0, 2181.0, 2182.0, 2183.0, 2185.0, 2186.0, 2187.0, 2188.0, 2190.0, 2191.0, 2192.0, 2194.0, 2195.0, 2196.0, 2197.0, 2198.0, 2200.0], "y": [-327612.125, -347306.34375, -371487.59375, -400072.0625, -427773.65625, -452044.53125, -457291.28125, -461900.9375, -463964.375, -467420.09375, -472346.75, -473178.28125, -479740.0, -482194.4375, -483604.21875, -485770.25, -486202.75, -487950.96875, -491239.84375, -491569.4375, -496079.96875, -497212.90625, -497915.25, -498762.6875, -498798.59375, -499044.6875, -503261.53125, -504590.9375, -505913.9375, -506153.25, -507574.0625, -508783.8125, -509316.40625, -510723.5625, -511569.15625, -512162.0, -513047.5625, -513196.71875, -514057.28125, -515344.9375, -515768.96875, -517365.65625, -517970.40625, -518392.53125, -518939.75, -519113.28125, -519502.125, -520976.59375, -521471.09375, -523377.875, -523769.65625, -524161.8125, -524393.8125, -524643.8125, -525843.375, -526626.4375, -528134.5, -528466.375, -528804.3125, -529018.9375, -529233.4375, -529242.3125, -530580.5625, -531358.5, -533036.125, -533256.9375, -533527.875, -533613.8125, -533751.8125, -534757.6875, -534980.8125, -536369.625, -536554.125, -536812.1875, -536875.625, -537043.8125, -537864.6875, -538098.6875, -539174.5625, -539362.5, -539596.75, -539682.875, -539870.0, -540481.625, -540768.8125, -541524.25, -541727.9375, -541908.1875, -542069.375, -542176.0, -542255.625, -542922.125, -543344.3125, -544164.875, -544323.3125, -544502.0, -544595.875, -544720.9375, -545209.8125, -545589.125, -546230.625, -546383.4375, -546560.0, -546650.125, -546800.875, -547170.375, -547469.5, -547933.25, -548091.0625, -548234.6875, -548365.5, -548466.1875, -548545.6875, -548559.3125, -549138.75, -549491.9375, -550260.625, -550351.75, -550511.875, -550521.625, -550657.0625, -551055.75, -551097.4375, -551613.5, -551739.8125, -551838.0625, -551935.25, -551975.4375, -552010.5625, -552502.375, -552619.375, -553261.8125, -553348.6875, -553464.625, -553493.125, -553570.1875, -553954.1875, -554041.8125, -554554.875, -554648.25, -554751.875, -554801.125, -554869.3125, -555198.0, -555378.8125, -555811.0, -555902.9375, -556006.1875, -556059.5625, -556136.0625, -556407.8125, -556605.9375, -556964.9375, -557057.875, -557153.875, -557216.4375, -557288.9375, -557297.125, -557313.5625, -557747.5, -557914.8125, -558078.1875, -558093.5, -558289.8125, -558423.4375, -558476.625, -558631.375, -558754.875, -558839.0, -558985.6875, -559091.25, -559186.375, -559305.0625, -559339.125, -559465.4375, -559617.25, -559699.75, -559889.875, -559985.875, -560059.0625, -560159.5625, -560189.625, -560290.5, -560469.5625, -560556.5, -560773.75, -560864.4375, -560929.875, -561021.1875, -561045.5625, -561132.0, -561328.0, -561404.125, -561644.5625, -561728.4375, -561787.9375, -561869.125, -561889.4375, -561959.0, -562176.4375, -562243.3125, -562512.0625, -562588.3125, -562645.375, -562713.25, -562733.0, -562781.75, -563025.125, -563086.5625, -563388.625, -563456.0, -563511.9375, -563564.5, -563586.375, -563610.1875, -563883.4375, -563945.9375, -564285.5625, -564342.6875, -564399.375, -564435.1875, -564460.0, -564685.9375, -564877.8125, -565181.3125, -565232.0625, -565313.875, -565327.0, -565408.75, -565557.0, -565603.6875, -565787.0, -565857.3125, -565904.75, -565975.0625, -565989.125, -566052.0, -566221.1875, -566269.75, -566481.5, -566545.0, -566591.5625, -566648.4375, -566665.9375, -566706.875, -566901.375, -566959.9375, -567204.6875, -567259.625, -567308.0, -567349.9375, -567373.1875, -567388.4375, -567612.0, -567696.0, -567976.4375, -568021.625, -568073.9375, -568098.375, -568128.9375, -568305.75, -568417.375, -568654.125, -568699.375, -568760.375, -568780.1875, -568836.0, -568965.9375, -569029.625, -569192.875, -569247.9375, -569292.5625, -569342.25, -569365.875, -569402.1875, -569553.5, -569639.25, -569830.875, -569879.375, -569925.8125, -569962.125, -569991.8125, -570005.125, -570181.3125, -570317.1875, -570539.9375, -570580.3125, -570629.9375, -570651.3125, -570688.0, -570821.9375, -570907.375, -571084.8125, -571128.9375, -571177.1875, -571205.0625, -571245.4375, -571355.0625, -571456.1875, -571595.3125, -571642.375, -571686.4375, -571725.0, -571757.9375, -571779.6875, -571788.0625, -571961.8125, -572065.875, -572296.375, -572327.875, -572379.5, -572386.0, -572430.75, -572554.8125, -572581.3125, -572741.0625, -572784.375, -572819.0625, -572854.0625, -572871.1875, -572887.8125, -573041.125, -573108.875, -573307.75, -573342.0, -573384.4375, -573399.9375, -573431.8125, -573552.125, -573609.5625, -573767.3125, -573805.4375, -573844.25, -573870.1875, -573898.0625, -573901.8125, -573905.1875, -574096.3125, -574171.125, -574242.0625, -574248.9375, -574334.75, -574391.5, -574415.75, -574480.3125, -574534.0625, -574573.875, -574638.0, -574683.125, -574726.25, -574776.375, -574795.875, -574848.6875, -574913.625, -574960.125, -575041.875, -575083.875, -575119.8125, -575161.5, -575183.75, -575220.5, -575220.75, -575241.75, -575388.625, -575388.875, -575585.0, -575612.5625, -575647.9375, -575657.375, -575681.5625, -575802.8125, -575834.25, -575990.1875, -576022.875, -576053.4375, -576076.375, -576093.25, -576096.5625, -576242.25, -576305.4375, -576492.5, -576517.875, -576552.9375, -576562.9375, -576585.1875, -576697.875, -576739.4375, -576892.6875, -576920.875, -576957.5, -576969.25, -577001.0, -577090.3125, -577128.4375, -577241.25, -577275.375, -577303.0, -577333.25, -577346.0, -577366.8125, -577471.75, -577513.375, -577646.1875, -577677.1875, -577704.625, -577728.0625, -577741.5, -577751.125, -577873.375, -577919.0625, -578073.8125, -578099.8125, -578128.375, -578144.25, -578158.5, -578260.5, -578340.625, -578481.0, -578505.5, -578544.1875, -578550.1875, -578589.125, -578660.5625, -578681.8125, -578769.8125, -578804.0, -578826.3125, -578861.3125, -578865.6875, -578898.0, -578980.9375, -578996.625, -579100.6875, -579132.6875, -579153.25, -579183.375, -579187.4375, -579210.4375, -579307.9375, -579319.1875, -579442.375, -579470.625, -579491.625, -579514.75, -579519.8125, -579531.75, -579646.4375, -579654.75, -579800.875, -579824.375, -579846.875, -579861.5625, -579869.1875, -579969.6875, -580049.4375, -580188.5, -580208.1875, -580244.5, -580245.25, -580282.5, -580349.9375, -580360.0625, -580444.75, -580475.875, -580493.8125, -580525.0, -580526.5625, -580553.875, -580635.0625, -580641.8125, -580745.4375, -580773.1875, -580791.5625, -580815.1875, -580819.0625, -580833.9375, -580931.625, -580941.5, -581067.3125, -581089.75, -581111.0, -581125.0, -581134.3125, -581221.3125, -581295.375, -581414.125, -581433.25, -581464.25, -581469.3125, -581499.5625, -581560.9375, -581578.8125, -581656.75, -581683.375, -581702.6875, -581727.75, -581735.0, -581754.4375, -581828.8125, -581852.6875, -581948.625, -581971.8125, -581992.5, -582009.9375, -582021.125, -582027.75, -582118.0625, -582161.9375, -582279.375, -582297.8125, -582321.75, -582330.125, -582347.0625, -582417.8125, -582451.8125, -582546.375, -582567.0, -582590.5625, -582602.6875, -582621.875, -582680.3125, -582722.5, -582797.8125, -582820.6875, -582842.125, -582860.125, -582875.3125, -582884.3125, -582884.9375, -582981.5, -583023.8125, -583153.875, -583166.875, -583194.5, -583232.75, -583265.3125, -583312.625, -583336.375, -583357.375, -583381.75, -583394.5625, -583416.6875, -583461.625, -583502.75, -583559.75, -583582.5, -583602.8125, -583623.25, -583637.4375, -583652.6875, -583654.4375, -583658.5625, -583756.0, -583793.5, -583830.0625, -583835.0, -583877.5625, -583908.5625, -583922.125, -583958.8125, -583985.875, -584006.125, -584037.0625, -584040.1875, -584073.5625, -584112.9375, -584126.0625, -584176.25, -584200.1875, -584216.625, -584240.6875, -584247.0625, -584268.5, -584319.5625, -584341.75, -584408.4375, -584429.125, -584447.25, -584463.1875, -584475.6875, -584482.3125, -584484.4375, -584571.375, -584590.5, -584703.875, -584717.5, -584737.125, -584740.625, -584753.75, -584821.125, -584834.0625, -584924.5, -584940.75, -584959.625, -584967.1875, -584981.625, -585037.9375, -585064.9375, -585138.0, -585155.875, -585173.3125, -585186.8125, -585198.0625, -585202.3125, -585270.4375, -585320.8125, -585409.375, -585424.0, -585443.875, -585450.6875, -585465.6875, -585518.0, -585544.625, -585614.75, -585631.75, -585649.9375, -585661.125, -585676.25, -585720.4375, -585760.8125, -585818.0625, -585835.875, -585853.4375, -585868.3125, -585881.0625, -585889.375, -585891.5625, -585963.1875, -586002.625, -586099.5, -586111.3125, -586133.5, -586161.9375, -586189.5, -586224.9375, -586244.4375, -586261.625, -586281.125, -586291.8125, -586310.5625, -586310.6875, -586325.4375, -586384.875, -586386.8125, -586464.375, -586480.3125, -586492.8125, -586504.75, -586509.4375, -586512.5, -586586.0625, -586597.75, -586693.25, -586705.375, -586721.6875, -586725.75, -586735.3125, -586794.25, -586809.4375, -586888.375, -586901.6875, -586919.0625, -586924.4375, -586939.125, -586985.8125, -587004.6875, -587065.3125, -587081.3125, -587095.6875, -587108.4375, -587116.25, -587122.5, -587179.625, -587211.375, -587285.4375, -587298.875, -587314.875, -587322.1875, -587333.375, -587379.3125, -587410.75, -587472.1875, -587486.75, -587503.5, -587512.0625, -587527.0625, -587564.1875, -587591.75, -587638.3125, -587654.75, -587669.1875, -587683.0625, -587692.75, -587702.3125, -587746.375, -587785.75, -587842.3125, -587856.8125, -587872.0625, -587882.75, -587894.1875, -587897.4375, -587899.125, -587966.875, -587977.5625, -588069.25, -588077.75, -588096.25, -588123.8125, -588146.9375, -588181.5, -588197.875, -588211.625, -588228.0, -588236.875, -588251.5, -588283.6875, -588312.5, -588353.375, -588369.0625, -588382.875, -588396.6875, -588406.4375, -588416.5, -588416.8125, -588419.625, -588489.125, -588514.1875, -588539.6875, -588542.375, -588572.0625, -588594.375, -588603.375, -588629.5, -588648.4375, -588661.8125, -588683.375, -588685.375, -588707.875, -588736.5, -588745.625, -588782.75, -588799.1875, -588810.5625, -588826.625, -588831.75, -588845.0, -588882.625, -588901.9375, -588952.1875, -588965.5, -588979.5, -588988.9375, -589000.125, -589001.625, -589006.125, -589067.0, -589067.75, -589146.875, -589157.125, -589169.625, -589173.875, -589179.9375, -589231.0, -589249.625, -589318.375, -589329.0, -589344.625, -589347.125, -589360.875, -589398.9375, -589409.3125, -589458.5, -589472.5, -589483.5, -589495.6875, -589500.25, -589507.625, -589554.0625, -589570.1875, -589630.4375, -589642.3125, -589654.4375, -589662.4375, -589669.5625, -589669.6875, -589726.0, -589757.375, -589830.25, -589839.6875, -589853.8125, -589856.875, -589866.875, -589909.1875, -589920.75, -589977.5625, -589988.875, -590002.125, -590008.6875, -590018.75, -590054.4375, -590077.1875, -590123.0625, -590135.625, -590147.875, -590158.1875, -590165.875, -590171.125, -590213.6875, -590248.25, -590303.375, -590314.4375, -590327.6875, -590334.4375, -590343.8125, -590378.0625, -590407.25, -590453.0625, -590465.3125, -590479.0625, -590486.625, -590498.6875, -590526.6875, -590550.9375, -590585.6875, -590599.3125, -590611.3125, -590623.1875, -590631.0625, -590640.1875, -590640.5625, -590643.0625, -590702.375, -590723.1875, -590744.875, -590747.5, -590772.5, -590791.5625, -590799.3125, -590822.3125, -590838.4375, -590849.75, -590868.0625, -590870.0625, -590889.0625, -590913.875, -590922.25, -590954.25, -590968.375, -590978.625, -590991.9375, -590997.0625, -591007.6875, -591040.25, -591059.8125, -591103.0, -591114.875, -591127.3125, -591135.0, -591145.5, -591146.0, -591151.125, -591202.0, -591221.875, -591241.0625, -591246.1875, -591267.9375, -591285.875, -591295.25, -591317.3125, -591331.875, -591344.0, -591360.0, -591364.1875, -591380.6875, -591403.75, -591416.75, -591446.125, -591459.1875, -591469.9375, -591482.125, -591489.375, -591498.5, -591498.9375, -591501.1875, -591556.25, -591577.5625, -591599.125, -591601.5625, -591627.4375, -591643.9375, -591651.4375, -591669.4375, -591686.0625, -591699.8125, -591719.9375, -591733.375, -591745.5, -591760.0, -591767.0625, -591781.0, -591801.3125, -591822.0, -591849.0625, -591861.25, -591873.4375, -591884.0, -591893.125, -591901.25, -591906.5625, -591907.8125, -591954.8125, -591979.75, -592044.5625, -592051.5, -592066.75, -592085.1875, -592101.125, -592123.875, -592136.4375, -592147.75, -592160.8125, -592167.625, -592180.5625, -592201.875, -592220.25, -592247.0, -592259.375, -592269.625, -592281.75, -592288.125, -592298.8125, -592324.0, -592346.0, -592377.6875, -592389.625, -592399.875, -592410.375, -592417.25, -592424.4375, -592424.6875, -592426.0625, -592480.125, -592498.875, -592518.0625, -592519.75, -592542.375, -592559.375, -592565.75, -592586.0625, -592600.625, -592610.625, -592626.5625, -592628.0, -592644.625, -592667.3125, -592673.8125, -592702.75, -592715.0625, -592723.4375, -592735.3125, -592739.4375, -592748.6875, -592778.25, -592793.25, -592833.0, -592842.625, -592854.0, -592860.3125, -592869.5, -592894.625, -592919.625, -592951.75, -592962.3125, -592972.8125, -592982.25, -592989.375, -592995.1875, -592996.125, -593037.75, -593069.6875, -593126.5, -593133.25, -593147.5625, -593163.3125, -593177.9375, -593197.4375, -593209.3125, -593220.0, -593232.5, -593238.75, -593251.1875, -593269.4375, -593285.5625, -593308.625, -593320.0, -593329.6875, -593340.875, -593346.8125, -593357.4375, -593379.0625, -593397.5, -593424.5, -593435.3125, -593445.125, -593455.5, -593461.1875, -593469.5, -593494.5625, -593518.125, -593549.8125, -593560.125, -593569.75, -593578.0625, -593584.9375, -593590.0625, -593590.6875, -593631.4375, -593658.25, -593713.75, -593720.1875, -593733.8125, -593749.75, -593763.0, -593782.375, -593793.375, -593803.375, -593814.875, -593820.4375, -593832.3125, -593850.875, -593865.5625, -593888.1875, -593898.4375, -593908.25, -593918.8125, -593924.4375, -593934.0625, -593955.5625, -593972.75, -593999.5625, -594009.875, -594019.0, -594028.1875, -594034.0625, -594041.25, -594066.375, -594088.5625, -594119.75, -594129.875, -594138.8125, -594147.0, -594153.25, -594157.625, -594157.9375, -594198.0625, -594219.0625, -594273.9375, -594279.5625, -594292.25, -594307.3125, -594320.3125, -594338.9375, -594349.6875, -594359.25, -594370.125, -594375.4375, -594386.625, -594404.5625, -594418.6875, -594440.625, -594450.875, -594459.625, -594469.75, -594474.75, -594484.0625, -594504.875, -594521.875, -594547.9375, -594558.125, -594566.0625, -594574.9375, -594580.625, -594587.125, -594611.75, -594633.3125, -594664.75, -594674.0, -594682.8125, -594689.6875, -594696.0, -594699.5, -594699.875, -594739.0, -594754.6875, -594808.625, -594813.5625, -594825.625, -594841.125, -594852.9375, -594871.75, -594881.875, -594890.75, -594901.0, -594905.75, -594916.0625, -594933.8125, -594947.3125, -594969.125, -594979.1875, -594987.3125, -594997.3125, -595001.8125, -595009.875, -595030.6875, -595046.8125, -595073.25, -595082.75, -595090.75, -595098.75, -595104.1875, -595109.6875, -595134.375, -595156.375, -595187.5625, -595196.1875, -595204.9375, -595211.0625, -595217.1875, -595219.5, -595220.3125, -595258.5625, -595268.75, -595321.4375, -595326.125, -595337.375, -595352.6875, -595364.125, -595382.6875, -595392.4375, -595400.6875, -595410.1875, -595415.0, -595424.5, -595442.125, -595455.5625, -595477.5, -595486.8125, -595494.6875, -595503.25, -595508.125, -595515.25, -595536.3125, -595552.75, -595578.9375, -595588.0, -595595.625, -595602.875, -595608.1875, -595612.9375, -595637.1875, -595659.9375, -595691.375, -595699.0625, -595707.6875, -595713.0, -595719.125, -595720.5, -595721.75, -595759.5625, -595764.5625, -595816.125, -595820.625, -595830.9375, -595845.625, -595857.4375, -595876.1875, -595885.375, -595893.0, -595902.125, -595906.5, -595915.1875, -595933.0, -595945.9375, -595968.25, -595976.75, -595984.3125, -595992.6875, -595997.3125, -596004.0, -596024.3125, -596041.125, -596067.4375, -596075.6875, -596083.25, -596090.125, -596095.6875, -596099.4375, -596123.3125, -596147.125, -596178.4375, -596185.75, -596194.0625, -596198.8125, -596205.3125, -596205.8125, -596207.0625, -596244.0, -596244.8125, -596295.3125, -596299.375, -596309.1875, -596324.25, -596335.375, -596353.875, -596362.8125, -596369.875, -596378.25, -596382.6875, -596390.9375, -596407.875, -596421.25, -596443.5, -596451.625, -596458.8125, -596466.375, -596470.875, -596476.625, -596497.4375, -596514.25, -596540.6875, -596548.25, -596555.375, -596561.625, -596567.125, -596569.875, -596570.1875, -596603.4375, -596618.0625, -596664.0, -596668.125, -596678.625, -596691.25, -596701.25, -596716.5625, -596725.4375, -596732.6875, -596741.6875, -596745.5625, -596754.5, -596769.5, -596779.9375, -596798.75, -596806.4375, -596813.75, -596821.9375, -596825.625, -596833.0, -596850.5625, -596863.0, -596884.75, -596892.625, -596899.8125, -596906.6875, -596911.125, -596915.9375, -596936.1875, -596953.0625, -596978.8125, -596986.4375, -596993.5625, -596999.125, -597004.25, -597006.4375, -597006.875, -597039.375, -597050.4375, -597095.5625, -597099.25, -597109.1875, -597121.6875, -597131.0625, -597146.6875, -597154.5, -597161.4375, -597170.1875, -597173.6875, -597182.8125, -597197.25, -597207.3125, -597225.375, -597233.4375, -597240.0, -597247.375, -597251.25, -597258.3125, -597275.6875, -597287.5625, -597309.5625, -597316.625, -597323.3125, -597330.125, -597334.25, -597338.6875, -597358.9375, -597375.8125, -597401.9375, -597408.6875, -597415.75, -597420.6875, -597425.5625, -597427.375, -597427.5, -597459.4375, -597467.1875, -597511.375, -597515.3125, -597524.4375, -597536.9375, -597546.125, -597561.5, -597569.1875, -597575.8125, -597583.8125, -597587.4375, -597595.5625, -597609.9375, -597620.0625, -597638.125, -597646.125, -597652.375, -597659.75, -597662.9375, -597669.375, -597686.25, -597699.125, -597720.9375, -597728.0, -597734.6875, -597740.5, -597744.5625, -597748.6875, -597769.0, -597786.6875, -597812.8125, -597818.8125, -597825.5625, -597830.375, -597835.5625, -597836.5, -597837.4375, -597868.5625, -597872.1875, -597915.375, -597918.8125, -597927.6875, -597940.3125, -597949.5, -597964.75, -597972.1875, -597978.6875, -597986.4375, -597989.75, -597997.3125, -598011.75, -598022.375, -598040.9375, -598048.0, -598054.0, -598060.875, -598064.6875, -598069.875, -598087.0625, -598100.625, -598123.25, -598129.4375, -598136.0625, -598141.5625, -598146.1875, -598149.0, -598169.6875, -598189.375, -598215.9375, -598222.0, -598229.125, -598233.0625, -598238.1875, -598238.5625, -598239.625, -598270.1875, -598282.75, -598294.375, -598296.125, -598310.875, -598319.4375, -598324.3125, -598334.0625, -598343.4375, -598351.25, -598363.1875, -598370.125, -598377.5, -598385.1875, -598389.375, -598396.9375, -598397.75, -598403.5, -598424.5, -598427.8125, -598455.5625, -598461.3125, -598467.4375, -598471.25, -598475.6875, -598494.1875, -598507.4375, -598531.3125, -598537.1875, -598543.5, -598548.375, -598552.75, -598554.375, -598554.8125, -598584.25, -598591.0625, -598631.4375, -598634.8125, -598643.5, -598654.6875, -598663.1875, -598677.25, -598684.6875, -598690.625, -598697.9375, -598701.25, -598708.3125, -598721.6875, -598731.375, -598748.25, -598755.0625, -598761.1875, -598767.75, -598771.3125, -598776.375, -598792.125, -598805.0625, -598825.4375, -598832.0625, -598837.9375, -598843.25, -598847.25, -598850.5625, -598850.875, -598877.25, -598893.75, -598929.9375, -598934.125, -598943.0, -598952.4375, -598961.125, -598972.9375, -598980.0625, -598986.4375, -598994.125, -598997.3125, -599004.875, -599016.75, -599025.875, -599040.6875, -599047.4375, -599053.4375, -599060.25, -599063.75, -599070.0625, -599083.6875, -599095.5, -599113.125, -599119.5, -599125.5, -599131.3125, -599135.4375, -599139.75, -599139.875, -599140.625, -599170.375, -599181.5, -599192.25, -599193.1875, -599206.125, -599215.3125, -599219.0625, -599230.125, -599238.5, -599244.0625, -599253.25, -599254.0625, -599263.875, -599276.125, -599279.625, -599295.375, -599302.25, -599307.3125, -599314.5625, -599316.375, -599322.375, -599338.0625, -599345.9375, -599367.1875, -599373.0, -599378.75, -599383.125, -599387.5625, -599388.9375, -599390.125, -599416.375, -599418.625, -599453.0625, -599457.375, -599462.875, -599464.5, -599468.125, -599489.0, -599495.0625, -599524.125, -599528.6875, -599535.1875, -599536.9375, -599542.375, -599559.1875, -599565.5, -599586.9375, -599592.8125, -599598.25, -599603.25, -599605.6875, -599608.25, -599628.75, -599640.0625, -599666.4375, -599671.25, -599677.4375, -599680.4375, -599684.5, -599701.0625, -599712.875, -599735.1875, -599740.5625, -599746.75, -599749.9375, -599755.625, -599768.8125, -599779.9375, -599796.5625, -599802.6875, -599808.0625, -599813.3125, -599817.25, -599820.6875, -599820.875, -599843.6875, -599864.875, -599896.0625, -599899.875, -599907.75, -599916.125, -599924.3125, -599935.0625, -599941.1875, -599947.0625, -599953.875, -599957.1875, -599964.0, -599974.25, -599983.375, -599995.9375, -600001.9375, -600007.625, -600013.75, -600017.375, -600023.1875, -600034.75, -600046.3125, -600061.1875, -600067.0, -600072.6875, -600078.125, -600082.0, -600085.9375, -600086.8125, -600087.5625, -600113.0625, -600113.8125, -600147.75, -600151.0625, -600156.5, -600157.125, -600160.625, -600181.25, -600184.1875, -600211.6875, -600216.0625, -600221.8125, -600223.4375, -600228.0, -600244.5, -600250.0625, -600271.875, -600276.8125, -600281.625, -600285.75, -600288.375, -600290.25, -600310.0, -600322.125, -600349.0625, -600353.1875, -600358.9375, -600360.75, -600365.25, -600380.9375, -600388.375, -600409.3125, -600414.4375, -600419.8125, -600422.875, -600427.0, -600440.3125, -600452.5625, -600469.8125, -600475.0625, -600480.5, -600485.125, -600488.625, -600490.75, -600491.625, -600513.25, -600523.1875, -600552.6875, -600556.0625, -600562.625, -600562.8125, -600569.0625, -600584.0, -600585.5625, -600604.9375, -600610.4375, -600614.5625, -600619.25, -600620.3125, -600623.25, -600642.1875, -600646.0625, -600670.4375, -600675.0625, -600679.5, -600682.4375, -600685.125, -600701.3125, -600713.8125, -600736.5625, -600740.9375, -600746.875, -600748.875, -600754.25, -600766.9375, -600773.0, -600789.375, -600794.4375, -600799.1875, -600803.875, -600806.6875, -600810.0, -600825.0625, -600835.6875, -600855.1875, -600859.9375, -600865.125, -600868.4375, -600872.0, -600872.9375, -600873.5625, -600897.125, -600898.3125, -600930.375, -600933.25, -600939.3125, -600949.125, -600956.625, -600968.75, -600974.1875, -600978.8125, -600984.4375, -600987.375, -600992.1875, -601003.5625, -601014.3125, -601029.125, -601034.375, -601039.1875, -601043.4375, -601047.25, -601050.25, -601050.8125, -601069.8125, -601083.9375, -601109.9375, -601113.375, -601119.5625, -601120.1875, -601126.25, -601139.4375, -601141.6875, -601158.5, -601163.8125, -601167.375, -601172.5, -601173.4375, -601177.1875, -601193.1875, -601197.5, -601218.5625, -601223.0625, -601227.5, -601230.375, -601233.125, -601233.3125, -601253.3125, -601264.5625, -601291.375, -601294.4375, -601300.0625, -601301.0, -601305.1875, -601320.0, -601323.375, -601343.125, -601347.25, -601352.125, -601355.125, -601358.1875, -601358.3125, -601377.125]}],
                            {"title": {"text": "Energy Minimization"}, "xaxis": {"title": {"text": "Energy Minimization Step"}}, "yaxis": {"title": {"text": "Potential Energy KJ/mol-1"}}},
                            {"showLink": false, "linkText": "Export to plot.ly", "plotlyServerURL": "https://plot.ly", "responsive": true}
                        ).then(function(){
                                
    var gd = document.getElementById('64a55da9-6e3f-43f2-b541-0c0c57afb87e');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })
                    };
                    });
                </script>
            </div>


\*\ **## Equilibrate the system (NVT) Equilibrate the**\ protein
system\ **in**\ NVT ensemble\ **(constant Number of particles, Volume
and Temperature). Protein**\ heavy atoms*\* will be restrained using
position restraining forces: movement is permitted, but only after
overcoming a substantial energy penalty. The utility of position
restraints is that they allow us to equilibrate our solvent around our
protein, without the added variable of structural changes in the
protein.

-  `Step 1 <#eqNVTStep1>`__: Creating portable binary run file for
   system equilibration
-  `Step 2 <#eqNVTStep2>`__: Equilibrate the **protein system** with
   **NVT** ensemble.
-  `Step 3 <#eqNVTStep3>`__: Checking **NVT Equilibration** results.
   Plotting **system temperature** by time during the **NVT
   equilibration** process. \*Building Blocks*\* used:
-  `Grompp <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp>`__
   from **biobb_md.gromacs.grompp**
-  `Mdrun <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun>`__
   from **biobb_md.gromacs.mdrun**
-  `GMXEnergy <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy>`__
   from **biobb_analysis.gromacs.gmx_energy** \**\*

### Step 1: Creating portable binary run file for system equilibration
(NVT) The **nvt** type of the **molecular dynamics parameters (mdp)
property** contains the main default parameters to run an **NVT
equilibration** with **protein restraints** (see `GROMACS mdp
options <http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html>`__):

-  Define = -DPOSRES
-  integrator = md
-  dt = 0.002
-  nsteps = 5000
-  pcoupl = no
-  gen_vel = yes
-  gen_temp = 300
-  gen_seed = -1

In this particular example, the default parameters will be used: **md**
integrator algorithm, a **step size** of **2fs**, **5,000 equilibration
steps** with the protein **heavy atoms restrained**, and a temperature
of **300K**.

*Please note that for the sake of time this tutorial is only running
10ps of NVT equilibration, whereas in the*\ `original
example <http://www.mdtutorials.com/gmx/lysozyme/06_equil.html>`__\ *the
simulated time was 100ps.*

.. parsed-literal::

    # Grompp: Creating portable binary run file for NVT Equilibration
    from biobb_md.gromacs.grompp import Grompp
    
    # Create prop dict and inputs/outputs
    output_gppnvt_tpr = pdbCode+'_gppnvt.tpr'
    prop = {
        'mdp':{
            'type': 'nvt',
            'nsteps': 5000,
            'dt': 0.002,
            'define': '-DPOSRES',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        }
    }
    
    # Create and launch bb
    Grompp(input_gro_path=output_min_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppnvt_tpr,  
           properties=prop).launch()

### Step 2: Running NVT equilibration

.. parsed-literal::

    # Mdrun: Running Equilibration NVT
    from biobb_md.gromacs.mdrun import Mdrun
    
    # Create prop dict and inputs/outputs
    output_nvt_trr = pdbCode+'_nvt.trr'
    output_nvt_gro = pdbCode+'_nvt.gro'
    output_nvt_edr = pdbCode+'_nvt.edr'
    output_nvt_log = pdbCode+'_nvt.log'
    output_nvt_cpt = pdbCode+'_nvt.cpt'
    
    # Create and launch bb
    Mdrun(input_tpr_path=output_gppnvt_tpr, 
          output_trr_path=output_nvt_trr, 
          output_gro_path=output_nvt_gro, 
          output_edr_path=output_nvt_edr, 
          output_log_path=output_nvt_log, 
          output_cpt_path=output_nvt_cpt).launch()

### Step 3: Checking NVT Equilibration results Checking **NVT
Equilibration** results. Plotting **system temperature** by time during
the NVT equilibration process.

.. parsed-literal::

    # GMXEnergy: Getting system temperature by time during NVT Equilibration  
    from biobb_analysis.gromacs.gmx_energy import GMXEnergy
    
    # Create prop dict and inputs/outputs
    output_nvt_temp_xvg = pdbCode+'_nvt_temp.xvg'
    prop = {
        'terms':  ["Temperature"]
    }
    
    # Create and launch bb
    GMXEnergy(input_energy_path=output_nvt_edr, 
              output_xvg_path=output_nvt_temp_xvg, 
              properties=prop).launch()

.. parsed-literal::

    import plotly
    import plotly.graph_objs as go
    
    # Read temperature data from file 
    with open(output_nvt_temp_xvg,'r') as temperature_file:
        x,y = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in temperature_file 
                if not line.startswith(("#","@")) 
            ])
        )
    
    plotly.offline.init_notebook_mode(connected=True)
    
    plotly.offline.iplot({
        "data": [go.Scatter(x=x, y=y)],
        "layout": go.Layout(title="Temperature during NVT Equilibration",
                            xaxis=dict(title = "Time (ps)"),
                            yaxis=dict(title = "Temperature (K)")
                           )
    })



.. raw:: html

            <script type="text/javascript">
            window.PlotlyConfig = {MathJaxConfig: 'local'};
            if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
            if (typeof require !== 'undefined') {
            require.undef("plotly");
            requirejs.config({
                paths: {
                    'plotly': ['https://cdn.plot.ly/plotly-latest.min']
                }
            });
            require(['plotly'], function(Plotly) {
                window._Plotly = Plotly;
            });
            }
            </script>
            



.. raw:: html

    <div>
            
            
                <div id="8140cc92-9d6b-4ba0-ba0e-7472fdfa447b" class="plotly-graph-div" style="height:525px; width:100%;"></div>
                <script type="text/javascript">
                    require(["plotly"], function(Plotly) {
                        window.PLOTLYENV=window.PLOTLYENV || {};
                        window.PLOTLYENV.BASE_URL='https://plot.ly';
                        
                    if (document.getElementById("8140cc92-9d6b-4ba0-ba0e-7472fdfa447b")) {
                        Plotly.newPlot(
                            '8140cc92-9d6b-4ba0-ba0e-7472fdfa447b',
                            [{"type": "scatter", "uid": "5321b266-9b40-4bfa-ae67-e22a09305fce", "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], "y": [300.741913, 299.490234, 300.505554, 302.478424, 301.352844, 299.399078, 302.126221, 300.414246, 299.876892, 300.086853, 302.289581]}],
                            {"title": {"text": "Temperature during NVT Equilibration"}, "xaxis": {"title": {"text": "Time (ps)"}}, "yaxis": {"title": {"text": "Temperature (K)"}}},
                            {"showLink": false, "linkText": "Export to plot.ly", "plotlyServerURL": "https://plot.ly", "responsive": true}
                        ).then(function(){
                                
    var gd = document.getElementById('8140cc92-9d6b-4ba0-ba0e-7472fdfa447b');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })
                    };
                    });
                </script>
            </div>


**## Equilibrate the system (NPT) Equilibrate the**\ *\ protein
system\ *\ **in**\ *\ NPT\ *\ **ensemble (constant Number of particles,
Pressure and Temperature). -**\ *\ *\ `Step
1 <#eqNPTStep1>`__\ *\ *\ **: Creating portable binary run file for
system equilibration -**\ *\ *\ `Step 2 <#eqNPTStep2>`__\ *\ *\ **:
Equilibrate the**\ *\ protein
system\ *\ **with**\ *\ NPT\ *\ **ensemble. -**\ *\ *\ `Step
3 <#eqNPTStep3>`__\ *\ *\ **: Checking**\ *\ NPT
Equilibration\ *\ **results. Plotting**\ *\ system pressure and
density\ *\ **by time during the**\ *\ NPT
equilibration\ *\ **process.** **Building Blocks** used: -
`Grompp <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp>`__
from **biobb_md.gromacs.grompp** -
`Mdrun <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun>`__
from **biobb_md.gromacs.mdrun** -
`GMXEnergy <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy>`__
from **biobb_analysis.gromacs.gmx_energy** \**\*

### Step 1: Creating portable binary run file for system equilibration
(NPT)

The **npt** type of the **molecular dynamics parameters (mdp) property**
contains the main default parameters to run an **NPT equilibration**
with **protein restraints** (see `GROMACS mdp
options <http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html>`__):

-  Define = -DPOSRES
-  integrator = md
-  dt = 0.002
-  nsteps = 5000
-  pcoupl = Parrinello-Rahman
-  pcoupltype = isotropic
-  tau_p = 1.0
-  ref_p = 1.0
-  compressibility = 4.5e-5
-  refcoord_scaling = com
-  gen_vel = no

In this particular example, the default parameters will be used: **md**
integrator algorithm, a **time step** of **2fs**, **5,000 equilibration
steps** with the protein **heavy atoms restrained**, and a
Parrinello-Rahman **pressure coupling** algorithm.

*Please note that for the sake of time this tutorial is only running
10ps of NPT equilibration, whereas in the*\ `original
example <http://www.mdtutorials.com/gmx/lysozyme/07_equil2.html>`__\ *the
simulated time was 100ps.*

.. parsed-literal::

    # Grompp: Creating portable binary run file for NPT System Equilibration
    from biobb_md.gromacs.grompp import Grompp
    
    # Create prop dict and inputs/outputs
    output_gppnpt_tpr = pdbCode+'_gppnpt.tpr'
    prop = {
        'mdp':{
            'type': 'npt',
            'nsteps':'5000',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        }
    }
    
    # Create and launch bb
    Grompp(input_gro_path=output_nvt_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppnpt_tpr, 
           input_cpt_path=output_nvt_cpt,  
           properties=prop).launch()

### Step 2: Running NPT equilibration

.. parsed-literal::

    # Mdrun: Running NPT System Equilibration
    from biobb_md.gromacs.mdrun import Mdrun
    
    # Create prop dict and inputs/outputs
    output_npt_trr = pdbCode+'_npt.trr'
    output_npt_gro = pdbCode+'_npt.gro'
    output_npt_edr = pdbCode+'_npt.edr'
    output_npt_log = pdbCode+'_npt.log'
    output_npt_cpt = pdbCode+'_npt.cpt'
    
    # Create and launch bb
    Mdrun(input_tpr_path=output_gppnpt_tpr, 
          output_trr_path=output_npt_trr, 
          output_gro_path=output_npt_gro, 
          output_edr_path=output_npt_edr, 
          output_log_path=output_npt_log, 
          output_cpt_path=output_npt_cpt).launch()

### Step 3: Checking NPT Equilibration results Checking **NPT
Equilibration** results. Plotting **system pressure and density** by
time during the **NPT equilibration** process.

.. parsed-literal::

    # GMXEnergy: Getting system pressure and density by time during NPT Equilibration  
    from biobb_analysis.gromacs.gmx_energy import GMXEnergy
    
    # Create prop dict and inputs/outputs
    output_npt_pd_xvg = pdbCode+'_npt_PD.xvg'
    prop = {
        'terms':  ["Pressure","Density"]
    }
    
    # Create and launch bb
    GMXEnergy(input_energy_path=output_npt_edr, 
              output_xvg_path=output_npt_pd_xvg, 
              properties=prop).launch()

.. parsed-literal::

    import plotly
    from plotly import tools
    import plotly.graph_objs as go
    
    # Read pressure and density data from file 
    with open(output_npt_pd_xvg,'r') as pd_file:
        x,y,z = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]),float(line.split()[2]))
                for line in pd_file 
                if not line.startswith(("#","@")) 
            ])
        )
    
    plotly.offline.init_notebook_mode(connected=True)
    
    trace1 = go.Scatter(
        x=x,y=y
    )
    trace2 = go.Scatter(
        x=x,y=z
    )
    
    fig = tools.make_subplots(rows=1, cols=2, print_grid=False)
    
    fig.append_trace(trace1, 1, 1)
    fig.append_trace(trace2, 1, 2)
    
    fig['layout']['xaxis1'].update(title='Time (ps)')
    fig['layout']['xaxis2'].update(title='Time (ps)')
    fig['layout']['yaxis1'].update(title='Pressure (bar)')
    fig['layout']['yaxis2'].update(title='Density (Kg*m^-3)')
    
    fig['layout'].update(title='Pressure and Density during NPT Equilibration')
    fig['layout'].update(showlegend=False)
    
    plotly.offline.iplot(fig)



.. raw:: html

            <script type="text/javascript">
            window.PlotlyConfig = {MathJaxConfig: 'local'};
            if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
            if (typeof require !== 'undefined') {
            require.undef("plotly");
            requirejs.config({
                paths: {
                    'plotly': ['https://cdn.plot.ly/plotly-latest.min']
                }
            });
            require(['plotly'], function(Plotly) {
                window._Plotly = Plotly;
            });
            }
            </script>
            



.. raw:: html

    <div>
            
            
                <div id="66848755-026a-48d7-a469-dc719e37fd88" class="plotly-graph-div" style="height:525px; width:100%;"></div>
                <script type="text/javascript">
                    require(["plotly"], function(Plotly) {
                        window.PLOTLYENV=window.PLOTLYENV || {};
                        window.PLOTLYENV.BASE_URL='https://plot.ly';
                        
                    if (document.getElementById("66848755-026a-48d7-a469-dc719e37fd88")) {
                        Plotly.newPlot(
                            '66848755-026a-48d7-a469-dc719e37fd88',
                            [{"type": "scatter", "uid": "c51c41f3-fbe3-4be4-b133-bb6f7f8fc524", "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], "xaxis": "x", "y": [-490.519562, -87.049194, -221.164536, -43.451157, 36.488144, -48.837311, -90.49791, 224.62384, 246.190063, 313.460724, -90.206894], "yaxis": "y"}, {"type": "scatter", "uid": "718cd668-716b-4a87-8446-d8786db3d158", "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], "xaxis": "x2", "y": [995.082031, 1016.396729, 1018.278259, 1025.212891, 1018.200439, 1018.467224, 1016.799316, 1022.625549, 1021.392273, 1020.431763, 1018.817871], "yaxis": "y2"}],
                            {"showlegend": false, "title": {"text": "Pressure and Density during NPT Equilibration"}, "xaxis": {"anchor": "y", "domain": [0.0, 0.45], "title": {"text": "Time (ps)"}}, "xaxis2": {"anchor": "y2", "domain": [0.55, 1.0], "title": {"text": "Time (ps)"}}, "yaxis": {"anchor": "x", "domain": [0.0, 1.0], "title": {"text": "Pressure (bar)"}}, "yaxis2": {"anchor": "x2", "domain": [0.0, 1.0], "title": {"text": "Density (Kg*m^-3)"}}},
                            {"showLink": false, "linkText": "Export to plot.ly", "plotlyServerURL": "https://plot.ly", "responsive": true}
                        ).then(function(){
                                
    var gd = document.getElementById('66848755-026a-48d7-a469-dc719e37fd88');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })
                    };
                    });
                </script>
            </div>


**## Free Molecular Dynamics Simulation Upon completion of the**\ *\ two
equilibration phases (NVT and NPT)\ *\ **, the system is now
well-equilibrated at the desired temperature and pressure.
The**\ *\ position restraints\ *\ **can now be released. The last step
of the**\ *\ protein\ *\ **MD setup is a short,**\ *\ free MD
simulation\ *\ **, to ensure the robustness of the system.
-**\ *\ *\ `Step 1 <#mdStep1>`__\ *\ *\ **: Creating portable binary run
file to run a**\ *\ free MD simulation\ *\ **. -**\ *\ *\ `Step
2 <#mdStep2>`__\ *\ *\ **: Run short MD simulation of the**\ *\ protein
system\ *\ **. -**\ *\ *\ `Step 3 <#mdStep3>`__\ *\ *\ **: Checking
results for the final step of the setup process, the**\ *\ free MD
run\ *\ **. Plotting**\ *\ Root Mean Square deviation
(RMSd)*\ **and**\ *\ Radius of Gyration (Rgyr)*\ **by time during
the**\ *\ free MD run\ *\ **step.** **Building Blocks** used: -
`Grompp <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp>`__
from **biobb_md.gromacs.grompp** -
`Mdrun <https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun>`__
from **biobb_md.gromacs.mdrun** -
`GMXRms <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rms>`__
from **biobb_analysis.gromacs.gmx_rms** -
`GMXRgyr <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rgyr>`__
from **biobb_analysis.gromacs.gmx_rgyr** \**\*

### Step 1: Creating portable binary run file to run a free MD
simulation

The **free** type of the **molecular dynamics parameters (mdp)
property** contains the main default parameters to run an **free MD
simulation** (see `GROMACS mdp
options <http://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html>`__):

-  integrator = md
-  dt = 0.002 (ps)
-  nsteps = 50000

In this particular example, the default parameters will be used: **md**
integrator algorithm, a **time step** of **2fs**, and a total of
**50,000 md steps** (100ps).

*Please note that for the sake of time this tutorial is only running
100ps of free MD, whereas in the*\ `original
example <http://www.mdtutorials.com/gmx/lysozyme/08_MD.html>`__\ *the
simulated time was 1ns (1000ps).*

.. parsed-literal::

    # Grompp: Creating portable binary run file for mdrun
    from biobb_md.gromacs.grompp import Grompp
    
    # Create prop dict and inputs/outputs
    output_gppmd_tpr = pdbCode+'_gppmd.tpr'
    prop = {
        'mdp':{
            'type': 'free',
            'nsteps':'50000',
            #'tc_grps': "DNA Water_and_ions" # NOTE: uncomment this line if working with DNA
        }
    }
    
    # Create and launch bb
    Grompp(input_gro_path=output_npt_gro, 
           input_top_zip_path=output_genion_top_zip, 
           output_tpr_path=output_gppmd_tpr, 
           input_cpt_path=output_npt_cpt, 
           properties=prop).launch()

### Step 2: Running short free MD simulation

.. parsed-literal::

    # Mdrun: Running free dynamics
    from biobb_md.gromacs.mdrun import Mdrun
    
    # Create prop dict and inputs/outputs
    output_md_trr = pdbCode+'_md.trr'
    output_md_gro = pdbCode+'_md.gro'
    output_md_edr = pdbCode+'_md.edr'
    output_md_log = pdbCode+'_md.log'
    output_md_cpt = pdbCode+'_md.cpt'
    
    # Create and launch bb
    Mdrun(input_tpr_path=output_gppmd_tpr, 
          output_trr_path=output_md_trr, 
          output_gro_path=output_md_gro, 
          output_edr_path=output_md_edr, 
          output_log_path=output_md_log, 
          output_cpt_path=output_md_cpt).launch()

### Step 3: Checking free MD simulation results Checking results for the
final step of the setup process, the **free MD run**. Plotting **Root
Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time
during the **free MD run** step. **RMSd** against the **experimental
structure** (input structure of the pipeline) and against the
**minimized and equilibrated structure** (output structure of the NPT
equilibration step).

.. parsed-literal::

    # GMXRms: Computing Root Mean Square deviation to analyse structural stability 
    #         RMSd against minimized and equilibrated snapshot (backbone atoms)   
    
    from biobb_analysis.gromacs.gmx_rms import GMXRms
    
    # Create prop dict and inputs/outputs
    output_rms_first = pdbCode+'_rms_first.xvg'
    prop = {
        'selection':  'Backbone',
        #'selection': 'non-Water'
    }
    
    # Create and launch bb
    GMXRms(input_structure_path=output_gppmd_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rms_first, 
              properties=prop).launch()

.. parsed-literal::

    # GMXRms: Computing Root Mean Square deviation to analyse structural stability 
    #         RMSd against experimental structure (backbone atoms)   
    
    from biobb_analysis.gromacs.gmx_rms import GMXRms
    
    # Create prop dict and inputs/outputs
    output_rms_exp = pdbCode+'_rms_exp.xvg'
    prop = {
        'selection':  'Backbone',
        #'selection': 'non-Water'
    }
    
    # Create and launch bb
    GMXRms(input_structure_path=output_gppmin_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rms_exp, 
              properties=prop).launch()

.. parsed-literal::

    import plotly
    import plotly.graph_objs as go
    
    # Read RMS vs first snapshot data from file 
    with open(output_rms_first,'r') as rms_first_file:
        x,y = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in rms_first_file 
                if not line.startswith(("#","@")) 
            ])
        )
    
    # Read RMS vs experimental structure data from file 
    with open(output_rms_exp,'r') as rms_exp_file:
        x2,y2 = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in rms_exp_file
                if not line.startswith(("#","@")) 
            ])
        )
        
    trace1 = go.Scatter(
        x = x,
        y = y,
        name = 'RMSd vs first'
    )
    
    trace2 = go.Scatter(
        x = x,
        y = y2,
        name = 'RMSd vs exp'
    )
    
    data = [trace1, trace2]
    
    plotly.offline.init_notebook_mode(connected=True)
    
    plotly.offline.iplot({
        "data": data,
        "layout": go.Layout(title="RMSd during free MD Simulation",
                            xaxis=dict(title = "Time (ps)"),
                            yaxis=dict(title = "RMSd (nm)")
                           )
    })




.. raw:: html

            <script type="text/javascript">
            window.PlotlyConfig = {MathJaxConfig: 'local'};
            if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
            if (typeof require !== 'undefined') {
            require.undef("plotly");
            requirejs.config({
                paths: {
                    'plotly': ['https://cdn.plot.ly/plotly-latest.min']
                }
            });
            require(['plotly'], function(Plotly) {
                window._Plotly = Plotly;
            });
            }
            </script>
            



.. raw:: html

    <div>
            
            
                <div id="0897032b-f64d-4a5b-9451-2dcd5a7f64c4" class="plotly-graph-div" style="height:525px; width:100%;"></div>
                <script type="text/javascript">
                    require(["plotly"], function(Plotly) {
                        window.PLOTLYENV=window.PLOTLYENV || {};
                        window.PLOTLYENV.BASE_URL='https://plot.ly';
                        
                    if (document.getElementById("0897032b-f64d-4a5b-9451-2dcd5a7f64c4")) {
                        Plotly.newPlot(
                            '0897032b-f64d-4a5b-9451-2dcd5a7f64c4',
                            [{"name": "RMSd vs first", "type": "scatter", "uid": "48c2de19-522a-457b-ac9e-090d67b90473", "x": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0], "y": [0.0, 0.0629591, 0.0567507, 0.0584769, 0.079135, 0.0661647, 0.0790737, 0.085431, 0.086772, 0.0745758, 0.0601999]}, {"name": "RMSd vs exp", "type": "scatter", "uid": "d4e8219a-b4de-42b0-95df-36b85b01ebf9", "x": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0], "y": [0.03315, 0.0680659, 0.0582172, 0.0617637, 0.0824452, 0.0671645, 0.0780765, 0.0852354, 0.0857494, 0.0728402, 0.0613163]}],
                            {"title": {"text": "RMSd during free MD Simulation"}, "xaxis": {"title": {"text": "Time (ps)"}}, "yaxis": {"title": {"text": "RMSd (nm)"}}},
                            {"showLink": false, "linkText": "Export to plot.ly", "plotlyServerURL": "https://plot.ly", "responsive": true}
                        ).then(function(){
                                
    var gd = document.getElementById('0897032b-f64d-4a5b-9451-2dcd5a7f64c4');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })
                    };
                    });
                </script>
            </div>


.. parsed-literal::

    # GMXRgyr: Computing Radius of Gyration to measure the protein compactness during the free MD simulation 
    
    from biobb_analysis.gromacs.gmx_rgyr import GMXRgyr
    
    # Create prop dict and inputs/outputs
    output_rgyr = pdbCode+'_rgyr.xvg'
    prop = {
        'selection':  'Backbone'
    }
    
    # Create and launch bb
    GMXRms(input_structure_path=output_gppmin_tpr,
             input_traj_path=output_md_trr,
             output_xvg_path=output_rgyr, 
              properties=prop).launch()

.. parsed-literal::

    import plotly
    import plotly.graph_objs as go
    
    # Read Rgyr data from file 
    with open(output_rgyr,'r') as rgyr_file:
        x,y = map(
            list,
            zip(*[
                (float(line.split()[0]),float(line.split()[1]))
                for line in rgyr_file 
                if not line.startswith(("#","@")) 
            ])
        )
    
    plotly.offline.init_notebook_mode(connected=True)
    
    plotly.offline.iplot({
        "data": [go.Scatter(x=x, y=y)],
        "layout": go.Layout(title="Radius of Gyration",
                            xaxis=dict(title = "Time (ps)"),
                            yaxis=dict(title = "Rgyr (nm)")
                           )
    })



.. raw:: html

            <script type="text/javascript">
            window.PlotlyConfig = {MathJaxConfig: 'local'};
            if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
            if (typeof require !== 'undefined') {
            require.undef("plotly");
            requirejs.config({
                paths: {
                    'plotly': ['https://cdn.plot.ly/plotly-latest.min']
                }
            });
            require(['plotly'], function(Plotly) {
                window._Plotly = Plotly;
            });
            }
            </script>
            



.. raw:: html

    <div>
            
            
                <div id="7d523280-f237-4cd1-8a21-a2f0bd864d77" class="plotly-graph-div" style="height:525px; width:100%;"></div>
                <script type="text/javascript">
                    require(["plotly"], function(Plotly) {
                        window.PLOTLYENV=window.PLOTLYENV || {};
                        window.PLOTLYENV.BASE_URL='https://plot.ly';
                        
                    if (document.getElementById("7d523280-f237-4cd1-8a21-a2f0bd864d77")) {
                        Plotly.newPlot(
                            '7d523280-f237-4cd1-8a21-a2f0bd864d77',
                            [{"type": "scatter", "uid": "0563512e-34d9-442e-b614-dde67a9cd002", "x": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0], "y": [0.03315, 0.0680659, 0.0582172, 0.0617637, 0.0824452, 0.0671645, 0.0780765, 0.0852354, 0.0857494, 0.0728402, 0.0613163]}],
                            {"title": {"text": "Radius of Gyration"}, "xaxis": {"title": {"text": "Time (ps)"}}, "yaxis": {"title": {"text": "Rgyr (nm)"}}},
                            {"showLink": false, "linkText": "Export to plot.ly", "plotlyServerURL": "https://plot.ly", "responsive": true}
                        ).then(function(){
                                
    var gd = document.getElementById('7d523280-f237-4cd1-8a21-a2f0bd864d77');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })
                    };
                    });
                </script>
            </div>


**## Post-processing and Visualizing resulting 3D trajectory
Post-processing and Visualizing the**\ *\ protein system\ *\ **MD
setup**\ *\ resulting
trajectory\ *\ **using**\ *\ NGL\ *\ **-**\ *\ *\ `Step
1 <#ppStep1>`__\ *\ *\ **: Imaging the resulting
trajectory,**\ *\ stripping out water molecules and
ions\ *\ **and**\ *\ correcting periodicity issues\ *\ **.
-**\ *\ *\ `Step 2 <#ppStep2>`__\ *\ *\ **: Generating a dry
structure,**\ *\ removing water molecules and ions\ *\ **from the final
snapshot of the MD setup pipeline. -**\ *\ *\ `Step
3 <#ppStep3>`__\ *\ *\ **: Visualizing the imaged trajectory using the
dry structure as a**\ *\ topology\ *\ **.** **Building Blocks** used: -
`GMXImage <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_image>`__
from **biobb_analysis.gromacs.gmx_image** -
`GMXTrjConvStr <https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str>`__
from **biobb_analysis.gromacs.gmx_trjconv_str** \**\*

### Step 1: *Imaging* the resulting trajectory. Stripping out **water
molecules and ions** and **correcting periodicity issues**

.. parsed-literal::

    # GMXImage: "Imaging" the resulting trajectory
    #           Removing water molecules and ions from the resulting structure
    from biobb_analysis.gromacs.gmx_image import GMXImage
    
    # Create prop dict and inputs/outputs
    output_imaged_traj = pdbCode+'_imaged_traj.trr'
    prop = {
        'center_selection':  'Protein',
        'output_selection': 'Protein',
        'pbc' : 'mol',
        'center' : True
    }
    
    # Create and launch bb
    GMXImage(input_traj_path=output_md_trr,
             input_top_path=output_gppmd_tpr,
             output_traj_path=output_imaged_traj, 
              properties=prop).launch()

### Step 2: Generating the output *dry* structure. **Removing water
molecules and ions** from the resulting structure

.. parsed-literal::

    # GMXTrjConvStr: Converting and/or manipulating a structure
    #                Removing water molecules and ions from the resulting structure
    #                The "dry" structure will be used as a topology to visualize 
    #                the "imaged dry" trajectory generated in the previous step.
    from biobb_analysis.gromacs.gmx_trjconv_str import GMXTrjConvStr
    
    # Create prop dict and inputs/outputs
    output_dry_gro = pdbCode+'_md_dry.gro'
    prop = {
        'selection':  'Protein'
    }
    
    # Create and launch bb
    GMXTrjConvStr(input_structure_path=output_md_gro,
             input_top_path=output_gppmd_tpr,
             output_str_path=output_dry_gro, 
              properties=prop).launch()

### Step 3: Visualizing the generated dehydrated trajectory. Using the
**imaged trajectory** (output of the `Post-processing step
1 <#ppStep1>`__) with the **dry structure** (output of the
`Post-processing step 2 <#ppStep2>`__) as a topology.

.. parsed-literal::

    # Show trajectory
    view = nglview.show_simpletraj(nglview.SimpletrajTrajectory(output_imaged_traj, output_dry_gro), gui=True)
    view

.. parsed-literal::

    from time import sleep
    # range number of frames for the animated gif trajectory
    for frame in range(0, 11):
        # set frame to update coordinates
        view.frame = frame
        # make sure to let NGL spending enough time to update coordinates
        sleep(0.5)
        view.download_image(filename='trj_image{}.png'.format(frame))
        # make sure to let NGL spending enough time to render before going to next frame
        sleep(2.0)

.. parsed-literal::

    import moviepy.editor as mpy

.. parsed-literal::

    # go to folder where the images are stored
    template = './trj_image{}.png'
    # get all (sorted) image files
    imagefiles = [template.format(str(i)) for i in range(0, 10, 1)]

.. parsed-literal::

    # make a gif file
    frame_per_second = 8
    im = mpy.ImageSequenceClip(imagefiles, fps=frame_per_second)
    im.write_gif('trajectory.gif', fps=frame_per_second)

.. parsed-literal::

    display.HTML("<img src='_static/trajectory.gif'></img>")




.. raw:: html

    <img src='_static/trajectory.gif'></img>



## Output files

Important **Output files** generated: - {{output_md_gro}}: **Final
structure** (snapshot) of the MD setup protocol. - {{output_md_trr}}:
**Final trajectory** of the MD setup protocol. - {{output_md_cpt}}:
**Final checkpoint file**, with information about the state of the
simulation. It can be used to **restart** or **continue** a MD
simulation. - {{output_gppmd_tpr}}: **Final tpr file**, GROMACS portable
binary run input file. This file contains the starting structure of the
**MD setup free MD simulation step**, together with the molecular
topology and all the simulation parameters. It can be used to **extend**
the simulation. - {{output_genion_top_zip}}: **Final topology** of the
MD system. It is a compressed zip file including a **topology file**
(.top) and a set of auxiliar **include topology** files (.itp).

**Analysis** (MD setup check) output files generated: -
{{output_rms_first}}: **Root Mean Square deviation (RMSd)** against
**minimized and equilibrated structure** of the final **free MD run
step**. - {{output_rms_exp}}: **Root Mean Square deviation (RMSd)**
against **experimental structure** of the final **free MD run step**. -
{{output_rgyr}}: **Radius of Gyration** of the final **free MD run
step** of the **setup pipeline**.
