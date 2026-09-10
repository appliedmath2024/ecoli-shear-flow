# E. coli Simulation in Shear Flow

Fortran source code for simulating a flagellated E. coli cell near a wall in a background flow, with a MATLAB script for rendering the simulated geometry and trajectory as a movie.

The supplied configuration uses four flagella, a wall-distance parameter of 0.8, and a shear-rate parameter of 10.0. Parameters are defined in the source files rather than a separate configuration file.

**Current execution requirement:** the supplied configuration starts from an existing simulation state (`frame = 80`). The required restart files are included in `Data/restart_data/`. See [Initialization and restart data](#initialization-and-restart-data) before running the simulation.

## Repository contents

```text
.
├── README.md
├── ecoli_shear_flow/
│   ├── Ecoli_OMP.f90          # Main program and simulation parameters
│   ├── commonval.f90          # Shared parameters and variables
│   ├── nrtype.f90             # Numerical type definitions
│   ├── interface_rod.f90      # Procedure interfaces
│   ├── *.f90                 # Geometry, force, and solver routines
│   └── makefile               # Build configuration
├── Data/
│   └── restart_data/ # Restart files for the supplied configuration
├── generate_ecoli_movie.m     # MATLAB visualization script
└── A_Movies/
    └── Movie_allFolders.mp4   # Existing example movie; also the default output
```

Generated simulation data and compiled executables are excluded from the supplied simulation directory. Many output files have a `.m` extension but contain numeric text, not MATLAB functions. `status.m` contains MATLAB variable assignments.

## Requirements

### Simulation

- A Linux/Unix-like build environment with `make`.
- Intel Fortran compiler `ifort`, as specified in the supplied makefile.
- OpenMP support and a linkable LAPACK library (`-llapack`).

The makefile uses Intel-specific options, including `-qopenmp`, `-xHost`, and `-mcmodel=medium`, and searches `/usr/lib/`. Adjust the compiler, flags, and library paths for your environment. Switching to a different compiler requires adapting the flags; compatibility has not been verified.

The main program explicitly requests 16 OpenMP threads with `OMP_set_NUM_THREADS(16)`. Edit that call to change the requested thread count. The `-xHost` option targets the build machine's instruction set, so rebuild on the target machine rather than distributing the executable as a portable binary.

### Visualization

- MATLAB with a graphical display available for figure capture.
- MPEG-4 support for the default MP4 output; a GIF output mode is also provided.

The script uses local functions, string arrays, `isfile`, `isfolder`, and reductions with the `'all'` option. The minimum supported MATLAB release has not been verified.

## Initialization and restart data

In `ecoli_shear_flow/Ecoli_OMP.f90`, the supplied settings are:

```fortran
data_dir = '../Data/restart_data/'
frame = 80
```

The restart path is relative to the simulation working directory. With the run commands below, it refers to `Data/restart_data/` in the repository root. That directory is included and contains the 11 files opened by the restart branch.

The following files are supplied for the current restart branch:

```text
coord_pt.m  triad1.m  triad2.m  triad3.m  body_pt.m
Fram.m      YCM.m     VCM.m     curv_tors.m  AD.m  AD2.m
```

The code opens all of these files. Geometry and director data must match the compiled dimensions and contain enough saved frames for the requested restart state. Change `data_dir` if the data are stored elsewhere.

The source also contains a fresh-initialization branch selected by `frame = 0`. To use it, edit the source and rebuild. This branch has not been validated as part of the publication preparation, and it does not reproduce the supplied restart state. No source settings have been changed automatically.

## Build and run

From the repository root, after configuring initialization:

```sh
cd ecoli_shear_flow
make
./ecoli-shortrod.out
```

`make` builds `ecoli-shortrod.out` and then removes generated `.mod` files. Rebuild after modifying simulation parameters.

Run from `ecoli_shear_flow` so that relative input and output paths resolve as documented. The program writes output files in its working directory and opens many of them with `status='replace'`, overwriting existing files with the same names. Use a separate working copy when retaining earlier results.

The default simulation duration is 3 seconds with a time step of `4e-8` seconds, corresponding to 75 million time steps. Runtime and memory use have not been benchmarked. For an initial local check, consider a shorter `Simulation_Time` and rebuild; such a run is not the full configured experiment.

## Main parameters

| Parameter | Supplied value | Source | Meaning |
| --- | --- | --- | --- |
| `Simulation_Time` | `3.0` | `Ecoli_OMP.f90` | Simulated duration in seconds |
| `dt` | `4e-8` | `Ecoli_OMP.f90` | Time step in seconds |
| `wall_switch` | `1` | `Ecoli_OMP.f90` | Enables the wall setting |
| `walldist` | `0.8` | `Ecoli_OMP.f90` | Wall-distance parameter |
| `sconst` | `10.0` | `Ecoli_OMP.f90` | Background-flow shear parameter |
| `frame` | `80` | `Ecoli_OMP.f90` | Restart selector; zero selects fresh initialization |
| `nrod` | `4` | `commonval.f90` | Number of flagella |
| `nPt` | `374` | `commonval.f90` | Number of discretization points per flagellum |

`commonval.f90` documents the units as microns, grams, and seconds. Changes to geometry or discretization must remain consistent with any restart data.

## Generate a movie

After generating simulation outputs, open `generate_ecoli_movie.m` in MATLAB and run it. Alternatively, from the repository root in MATLAB:

```matlab
generate_ecoli_movie
```

The script resolves paths relative to its own location and loads these files from `ecoli_shear_flow`:

- `status.m`: simulation metadata, including dimensions and output interval.
- `coord_pt.m`: flagellar centerline coordinates.
- `triad1.m`, `triad2.m`: material directors used to construct filament surfaces.
- `body_pt.m`, `body_face.m`: cell-body vertices and triangle connectivity.
- `AD.m`: data used to draw the projected trajectory.

The script reports missing inputs if simulation outputs are unavailable. It renders the flagella and cell body together with the wall and a trajectory projected onto the wall plane.

Editable settings in the script include:

| Setting | Default | Purpose |
| --- | --- | --- |
| `folder_list` | `"ecoli_shear_flow"` | Simulation folders relative to the script |
| `legend_name` | `{}` | Optional legend labels |
| `xz_switch` | `1` | Swaps x and z for the displayed view |
| `movie_type` | `2` | `1`: GIF; `2`: MP4 |
| `frame_stride` | `2` | Renders every second saved frame |
| `M.FrameRate` | `150` | MP4 playback frame rate |
| `M.Quality` | `75` | MP4 quality setting |

Output is saved as `A_Movies/Movie_allFolders.mp4` or `A_Movies/Movie_allFolders.gif`. The existing file of the same name is overwritten. The title reports simulation time based on `dt` and `nSkip`; playback time is different. GIF mode uses zero frame delay, whose playback behavior depends on the viewer.

## Validation status

The source files referenced by the makefile were checked to be present after removing unused source files and generated outputs. Compilation, numerical validation, and end-to-end movie generation have not been completed during publication preparation. An attempted MATLAB check could not run because MATLAB failed to start with a file-system error.

The restart files are included. Row-count checks found 81 saved frames in the centerline, director, body-coordinate, body-frame, and center-of-mass files, sufficient for the current `frame = 80` read loop. This is a structural check, not numerical or end-to-end validation. Before publishing a validated release, record the compiler and MATLAB versions, operating system, and a successful example run.

## Authors

- **Sookkyung Lim** — Department of Mathematical Sciences, University of Cincinnati, Cincinnati, OH 45221, USA
  Email: [limsg@ucmail.uc.edu](mailto:limsg@ucmail.uc.edu)
- **Yongsam Kim** — Department of Mathematics, Chung-Ang University, Dongjak-gu, Heukseok-dong, Seoul 06974, Republic of Korea
  Email: [kimy@cau.ac.kr](mailto:kimy@cau.ac.kr)
- **Wanho Lee** — National Institute for Mathematical Sciences, Daejeon 34047, Republic of Korea
  Email: [wlee@nims.re.kr](mailto:wlee@nims.re.kr)

## License and citation

**Citation and notification:** If you use or modify this code in your research, please cite the archived software release using its DOI. The associated manuscript is currently in preparation; publication details will be added when available. If possible, please let the authors know about your use of the code. We would be interested to hear about applications or extensions of the code.

A license has not yet been selected, and a `LICENSE` file is not included. Public availability alone does not specify permission to reuse or redistribute the code. License information should be finalized before release.

A release version and a DOI have not yet been assigned in this repository. The authors are listed above. A `CITATION.cff` file can be added with the author and release information; this section will be updated with the archived release DOI when available. No DOI or publication citation is claimed here.
