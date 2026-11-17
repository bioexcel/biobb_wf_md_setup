# flake8: noqa
# ruff: noqa
# marimo-only cells for biobb_wf_md_setup_marimo.py, merged into the Jupyter notebook's
# cells by ipynb_to_marimo.py (Claude Code skill: ipynb-to-marimo). Not a runnable file:
# names like mo, gmx_energy or output_min_edr come from other cells of the notebook.
#
# Each "# %%" line starts a cell. A directive says where a group of cells goes:
#   # %% replace <hash>: <snippet>   replaces the code cell containing <snippet>
#   # %% after: <snippet> / before: <snippet>
#   # %%                             next cell of the same group
# <hash> guards the replaced cell: if it changes in the .ipynb, regeneration stops.


# %% replace 5dc11a29: pdbCode = "1AKI"
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

# %%
# Nothing below this cell runs until the form above is submitted
mo.stop(pdb_form.value is None, mo.md("**Enter a PDB code and click _Run setup_ to start the pipeline.**"))
pdbCode = pdb_form.value.strip().upper()
mo.md(f"Setting up an MD simulation for **{pdbCode}**")


# %% after: Pressure and Density during NPT Equilibration
mo.md(r"""
### Explore the energy files
Choose a **simulation step** and an **energy term**: the term is extracted from that step's
energy file with **GMXEnergy** and plotted right away. Change either selector and the plot updates.
""")

# %%
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

# %%
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


# %% replace dc8276cb: # Show trajectory
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

# %%
# Show trajectory (redrawn whenever one of the controls above changes)
import subprocess

# Number of frames of the imaged trajectory (from GROMACS gmx check)
check = subprocess.run(['gmx', 'check', '-f', output_imaged_traj], capture_output=True, text=True)
n_frames = int(next(line.split()[1] for line in check.stderr.splitlines() if line.startswith('Coords')))

builder = mvs.create_builder()
# Imaged trajectory coordinates
builder.download(url=output_imaged_traj).parse(format='trr').coordinates(ref='trajectory')
# Dry structure used as topology
structure = builder.download(url=output_dry_gro).parse(format='gro').model_structure(coordinates_ref='trajectory', ref='structure')
color = {'molstar_color_theme_name': traj_color.value}
if traj_color.value == 'element-symbol':
    # Grey carbons (Mol* colours them by chain by default)
    color['molstar_color_theme_params'] = {'carbonColor': {'name': 'element-symbol', 'params': {}}}
structure.component(selector='protein').representation(type=traj_style.value).color(custom=color)
# Play the trajectory in a loop
builder.animation(autoplay=True, loop=True).interpolate(
    kind='scalar', target_ref='structure', property='model_index',
    start=0, end=n_frames - 1, discrete=True, duration_ms=n_frames * traj_speed.value
)
data = {
    output_dry_gro: Path(output_dry_gro).read_bytes(),
    output_imaged_traj: Path(output_imaged_traj).read_bytes()
}
mo.iframe(builder.molstar_html(data=data), height="600px")
