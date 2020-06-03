 ```

$$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\   $$$$$$\  $$$$$$$\  
$$  __$$\\_$$  _|$$ | $$  |$$  __$$\ $$ ___$$\ $$  __$$\ 
$$ |  $$ | $$ |  $$ |$$  / $$ /  \__|\_/   $$ |$$ |  $$ |
$$$$$$$  | $$ |  $$$$$  /  \$$$$$$\    $$$$$ / $$ |  $$ |
$$  ____/  $$ |  $$  $$<    \____$$\   \___$$\ $$ |  $$ |
$$ |       $$ |  $$ |\$$\  $$\   $$ |$$\   $$ |$$ |  $$ |
$$ |     $$$$$$\ $$ | \$$\ \$$$$$$  |\$$$$$$  |$$$$$$$  |
\__|     \______|\__|  \__| \______/  \______/ \_______/ 

        Parallel Image-based Kinetic Solver (3D)
```



## What is PIKS3D?

PIKS3D is an open-source 3D parallel pore-scale rarefied gas-flow simulator. 
It solves the linearized gas-kinetic equation on a uniform Cartesian grid using the Discrete
Velocity Method (DVM). It can be run in parallel both OpenMP and MPI. The porous
structure can be arbitrary complex and and is input into the simulator as an binary
images of '0' (fluid) and '1' (solid). Below is an example showing the flow field in a sandstone simulated by the solver.

<p align="center"><a href="https://ibb.co/6yVgGFJ"><img src="https://i.ibb.co/GdBcGR3/gas-porous-media.png" alt="gas-porous-media" border="0" width="480"></a></p>


## How do I use PIKS3D?

See the [wiki](https://github.com/iPACT-Platform/PIKS3D/wiki) pages for installation and tutorials.
See the reference below for the theory and numerical method. The full-text PDFs are provided in the PIKS2D's directory [reference](https://github.com/iPACT-Platform/PIKS2D/reference).

## How do I cite PIKS3D?

* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Zhi-Hui Li, and Yonghao Zhang. “A Multi-Level Parallel Solver for Rarefied Gas Flows in Porous Media.” Computer Physics Communications 234 (January 1, 2019): 14–25. [DOI: 10.1016/j.cpc.2018.08.009](https://doi.org/10.1016/j.cpc.2018.08.009).
* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Jingsheng Ma, and Yonghao Zhang. “Pore-Scale Simulations of Rarefied Gas Flows in Ultra-Tight Porous Media.” Fuel 249 (August 1, 2019): 341–51. [DOI: 10.1016/j.fuel.2019.03.106](https://doi.org/10.1016/j.fuel.2019.03.106).


## License

The PIKS3D is licensed under the MIT license, see the file `LICENSE`.

## Who have funded PIKS3D
Development of PIKS3D received funding from Engineering and Physical Sciences Research Council, European Union’s Horizon 2020 Marie Skłodowska-Curie Individual Fellowship and the “Advanced Hybrid Method for Pore-Scale Simulation of Shale Gas Flows”, a global partnership project grant funded by [KFUPM](http://www.kfupm.edu.sa).
