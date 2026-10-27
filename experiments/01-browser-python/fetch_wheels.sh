#!/bin/sh
# The wheels the page installs into the browser's Python: tifffile from PyPI, gpxiltools built from its repository.
set -e
cd "$(dirname "$0")"
mkdir -p wheels
pip download --no-deps --only-binary=:all: --python-version 3.13 --platform any -d wheels tifffile
pip wheel --no-deps -w wheels "git+https://github.com/comparativechrono/Gpxil-tools"
