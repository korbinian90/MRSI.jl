# MRSI

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://korbinian90.github.io/MRSI.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://korbinian90.github.io/MRSI.jl/dev/)
[![Build Status](https://github.com/korbinian90/MRSI.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/korbinian90/MRSI.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/korbinian90/MRSI.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/korbinian90/MRSI.jl)

# Install this package
1. [Download](https://julialang.org/downloads/) and install a recent version of Julia (1.6+)
2. Open the Julia REPL and type
    ```
    julia> ] # enters the package manager
    (@v1.8) pkg> add https://github.com/korbinian90/MRSI.jl
    ```  
    All dependencies will be installed automatically

# Quick Start
It is recommended to use [Julia with VS Code](https://code.visualstudio.com/docs/languages/julia#_getting-started).

```julia
using MRSI
dat_file = "path-to-dat-file/file.dat"
image = reconstruct(dat_file)

# only extract info from the scan
headers, info = read_scan_info(dat_file, :ONLINE)
```
