# Quad Mesh Quantization Without a T-Mesh
Yoann Coudert--Osmont, David Desobry, Martin Heistermann, David Bommes, Nicolas Ray and Dmitry Sokolov


## Compilation
Prerequisites:
- CMake
- C++ compiler

Run in a terminal:
```
mkdir build/
cd build/
cmake ..
make
```

The program `Quantization` is generated in the `build/` subdirectory.

## Usage 
From the folder `./Quantization/build/` use the command line:
```
./Quantization [OPTION] input_seamless_map.obj
```

OPTION:
- -c: for coarse quantization
- -d `matrix file path`: output the matrix in the file given in parameter if `scale` is a numeric value. Otherwise if `scale` is equal to "a", an automatic scale value is computed.
- -h: show this help.
- -i: imprint the quantized map in the original mesh.
- -r: re-embed the quantized map in the original mesh.
- -o `output path`: change the file in which the output will be written. Default is "out.obj".
- -s `scale`: scale the input seamless map by `scale` if `scale` is a numeric value. Otherwise if `scale` is "a", an automatic scale value is computed.
- -sa `scale`: scale the autoscale value.

## License
Copyright (C) 2022, Coudert--Osmont Yoann \
SPDX-License-Identifier: AGPL-3.0-or-later 

If you make use of this code in scientific work we kindly ask you to cite our paper:
```
@inproceedings{coudert2024quad,
  title={Quad Mesh Quantization Without a T-Mesh},
  author={Coudert-Osmont, Yoann and Desobry, David and Heistermann, Martin and Bommes, David and Ray, Nicolas and Sokolov, Dmitry},
  booktitle={Computer Graphics Forum},
  volume={43},
  number={1},
  year={2024}
}
``` 