### Development and Distribution

To build the pip wheels, from the root directory assuming Docker+conda are installed 
on Mac OS:
```bash
source dev/build_mac_wheels.sh && docker run --rm -e PLAT=manylinux2014_x86_64 -v `pwd`:/io -w "/io" quay.io/pypa/manylinux2014_x86_64 bash dev/build_linux_wheels.sh
```

To upload to PyPi first test with
```bash
python3 -m twine upload --repository testpypi dist/*
```
then to the *real* PyPi
```bash
twine upload dist/*
```
