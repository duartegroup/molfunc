
#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/ -v
done

# Bundle external shared libraries into the wheels
mkdir -p dist
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
    mv "$whl" dist/
done
