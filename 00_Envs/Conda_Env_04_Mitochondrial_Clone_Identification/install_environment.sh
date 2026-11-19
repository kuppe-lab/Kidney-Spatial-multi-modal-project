#!/usr/bin/env

# Since the setup can go wrong in a million ways, we make the script fail and stop when any error occurs.
set -euo pipefail

environment_prefix="/beegfs/data/users/mgrasshoff/spatial_sigurd_py/spatial_sigurd_py_env2"
spatial_sigurd_commit="f072fb6d4bea2976b5d216e594429e2809da06be"

conda env create -y -f "environment.yml" -p "${environment_prefix}"

UV_LINK_MODE=copy "${environment_prefix}/bin/uv" pip install \
    --python "${environment_prefix}/bin/python" \
    --override "overrides.txt" \
    --constraint "constraints.txt" \
    "spatial-sigurd-py @ git+https://github.com/CostaLab/spatial_sigurd_py.git@${spatial_sigurd_commit}"

# The version of SpaGFT is stale. The authors forgot to update it. The meta data is correct though.
"${environment_prefix}/bin/python" -c "from importlib.metadata import version; print('spagft', version('spagft'))"

