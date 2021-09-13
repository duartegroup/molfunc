
#!/bin/bash

for minor_version in {7..9}; do
    conda create -n test "python=3.$minor_version" twine --yes
    conda activate test
    python3 -m pip install --upgrade build
    python3 -m build
    conda deactivate
done
