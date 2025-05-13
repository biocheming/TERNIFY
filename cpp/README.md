### COMPILE CPP VERSION TERNIFY

1. Set up the rdkit lib in the CMakeLists.txt

2. Run the commands

```
mkdir build
cd build
cmake ..
make -j 10
if [ -f "ternify" ];then
	cp ternify /usr/local/bin/
fi
cd ..
```

3. Run ternify
```
ternify tcs.inp
```
