pyDBT
======

This repository is a python extension of the [DBT toolbox](https://github.com/LAVI-USP/DBT-Reconstruction) from LAVI-USP. 



## How to install?

 1. Download the toolbox or clone the directory:
    
    * ```git clone https://github.com/LAVI-USP/pyDBT.git```

 2. Go to parent directory:

    * ```cd pyDBT```

 3. Clone NVIDIA cuda-samples directory inside pyDBT:
    
    * ```git clone https://github.com/NVIDIA/cuda-samples```

 3. Install the package:

    * ```python3 setup.py install```

 3. If you have problems with `arch=sm_XX`, modify it in the `setup.py` accordingly to your NVIDIA-GPU architecture. This [link](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/) has some references.

 4. run the example:

    * ```cd pydbt && python3 example.py```
 
 5. The toolbox was tested on **Linux** (Ubuntu 20) x64, and **macOS** BigSur (Intel) machines, with python **3.7.5**.
     * You will need to install gcc for the compilation. You can install gcc 11 on macOS through [homebrew](https://formulae.brew.sh/formula/gcc).
     * Fell free to reach me if you want binaries for either Ubuntu or BigSur
 
 6. You can also run the [MATLAB version](https://github.com/LAVI-USP/DBT-Reconstruction) of the toolbox.

** Please report issues [here](https://github.com/LAVI-USP/pyDBT/issues). **

## Contribute? 

We are pleased with any contributions. Feel free to make any [pull requests](https://github.com/LAVI-USP/pyDBT/pulls) or send us an e-mail.


## Toolbox manual:

You can find the [manual](https://github.com/LAVI-USP/DBT-Reconstruction/wiki/Toolbox-Manual) from the MATLAB version, which is pretty much the same. I will create a specific one for the python version in the future.

## Contact:

If you have any questions or suggestion, please send us an e-mail:

- Rodrigo - rodrigo dot vimieiro at gmail dot com
- Marcelo - mvieira at sc dot usp dot br

## License:

The toolbox is licensed under the **GNU General Public License v3.0**. Please check the [licence file](https://github.com/LAVI-USP/pyDBT/blob/master/LICENSE).

## Reference:

If you use the toolbox, we will be very grateful if you refer to this [paper](https://doi.org/10.1007/978-981-13-2517-5_53):

> Vimieiro R.B., Borges L.R., Vieira M.A.C. (2019) Open-Source Reconstruction Toolbox for Digital Breast Tomosynthesis. In: Costa-Felix R., Machado J., Alvarenga A. (eds) XXVI Brazilian Congress on Biomedical Engineering. IFMBE Proceedings, vol 70/2. Springer, Singapore.

## Citations:

You can find [here](https://scholar.google.com.br/scholar?oi=bibs&hl=pt-BR&cites=3156269064066227282) the papers that have used the toolbox.

## Acknowledgments:

This work was supported by the São Paulo Research Foundation ([FAPESP](http://www.fapesp.br/) grant 2016/25750-0) and by the National Council for Scientific and Technological Development ([CNPq](http://www.cnpq.br/)). Nobody does anything alone, so we would like to thank the contribution of our lab members and the [Barretos Love Hospital](https://www.hcancerbarretos.com.br) for providing the images of DBT.

---

Laboratory of Computer Vision ([Lavi](http://iris.sel.eesc.usp.br/lavi/))  
Department of Electrical and Computer Engineering  
São Carlos School of Engineering, University of São Paulo  
São Carlos - Brazil
