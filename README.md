# MRSI

Package for CRT MRSI reconstruction with low RAM requirement.

## Implemented steps

- Reading from Siemens DAT file
- FOV Shift
- Frequency Offset Correction
- Density Compensation
- Discrete Fourier Transform (from circles to cartesian)
- Slice FFT
- Coil combination with PATREFSCAN
- Regularization (L1 and L2)

## Install this package

1. [Download](https://julialang.org/downloads/) and install a recent version of Julia (1.9+ recommended)

2. Open the Julia REPL and type

    ```julia-repl
    julia> ] # enters the package manager
    (@v1.9) pkg> add https://github.com/korbinian90/MRSI.jl
    ```

    All dependencies will be installed automatically

## Quick Start

It is recommended to use [Julia with VS Code](https://code.visualstudio.com/docs/languages/julia#_getting-started).

```julia
using MRSI
dat_file = "path-to-dat-file/file.dat"
image = reconstruct(dat_file)

image = reconstruct(dat_file; old_headers=true) # for old sequences with non-standardized header entries

# only extract info from the scan
headers, info = read_scan_info(dat_file, :ONLINE)
```

Saving to NIfTI

```julia
import Pkg; Pkg.add("MriResearchTools") # installs package, many dependencies, takes a minute or two
using MriResearchTools

savenii(abs.(combined), "path-to-file/mag.nii")
savenii(angle.(combined), "path-to-file/phase.nii")
```

Saving as MAT file

```julia
import Pkg; Pkg.add("MATLAB") # first time install of package
using MATLAB # automatically finds and uses installed MATLAB
write_matfile(filename; mat_var_name=combined, n_frequency=info[:n_frequency], ...)  # writes all variables given in the keyword argument list to a MAT file
```

## Options
To see the options of the `reconstruct` functions, use the help menu of the julia REPL

```julia-repl
julia> using MRSI
julia> ?reconstruct
```

## Work with the source code

This downloads the package to the `.julia/dev` directory, usually in `home` directory. The source code can be easily modified and used.

```julia
import Pkg; Pkg.dev(url="https://github.com/korbinian90/MRSI.jl")
using Revise # makes sure that changed code is immediately available in running session
using MRSI
# <Perform some modification in MRSI/src>
dat_file = "path-to-dat-file/file.dat"
image = reconstruct(dat_file) # uses changed package
```

## Questions

Please don't hesitate to write to me, or create an issue on this repository.
