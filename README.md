# FastQSL 2

To calculate the squashing factor $Q$, and other quatities related to the magnetic connectivity, at the bottom or at a cross section or in a box volume or on some seed points, given a 3D magnetic field on a Cartesian or spherical, uniformed or stretched grid.

This program can be downloaded with the command
```
git clone https://github.com/FastQSL2
```

Please address comments and suggestions to [Dr. Chen, Jun (陈俊)](mailto:chenjun@pmo.ac.cn)

If your markdown reader can not can not render the formulae in README\.md, please read README.html directly.

-----------------------------
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This program is licensed under a [CC BY-NC-SA 4.0 License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

-----------------------------
## Cite as

* Jun Chen*, Thomas Wiegelmann, Li Feng*, Bernhard Kliem, Chaowei Jiang, and Rui Liu. FastQSL 2: A Comprehensive Toolkit for Magnetic Connectivity Analysis.  2026, SCIENCE CHINA Physics, Mechanics & Astronomy, submitted

* [Peijin Zhang, Jun Chen*, Rui Liu and ChuanBing Wang. FastQSL: A Fast Computation Method for Quasi-separatrix Layers. 2022, The Astrophysical Journal, 937, 26](https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)

```bibtex
@ARTICLE{2022ApJ...937...26Z,
       author = {{Zhang}, PeiJin and {Chen}, Jun and {Liu}, Rui and {Wang}, ChuanBing},
        title = "{FastQSL: A Fast Computation Method for Quasi-separatrix Layers}",
      journal = {\apj},
     keywords = {Solar magnetic fields, GPU computing, 1503, 1969},
         year = 2022,
        month = sep,
       volume = {937},
       number = {1},
          eid = {26},
        pages = {26},
          doi = {10.3847/1538-4357/ac8d61},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ApJ...937...26Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
-----------------------------
## Software for Interface

### If use fastqsl\.pro
* install **IDL** https://www.nv5geospatialsoftware.com/Products/IDL 
  * setting an environmental variable of `$IDL_PATH` and placing fastqsl\.pro into a private path is suggested. If your IDL is installed at /usr/local/exelis/idl, just append such line to ~/.bashrc:
     ```bash
    export IDL_DIR=/usr/local/exelis/idl
    export IDL_PATH="$IDL_IDL/lib:+$HOME/your/private/pro/path"
    ``` 
* or install **GDL** https://gnudatalanguage.github.io/downloads.html
  * setting an environmental variable of `$GDL_PATH` is **necessary for write\_png**, and placing fastqsl\.pro into a private path is suggested. The default `$GDL_DIR` is depends on the distribution:
    * Ubuntu & Fedora:  /usr/share/gnudatalanguage
    * ArchLinux: /usr/lib/gdl
    * Gentoo: /usr/local/share/gdl
    * macOS: /opt/local/share/gnudatalanguage
    **Please append such lines to ~/.bashrc**, e.g. for Ubuntu
    ```bash
    export GDL_DIR=/usr/share/gnudatalanguage
    export GDL_PATH="$GDL_IDL/lib:+$HOME/your/private/pro/path"
    ```
  * If you also need [SSW](http://www.lmsal.com/solarsoft/) for some other analysis, please take a look at https://github.com/rbluosolar/sswgdl
### If use fastqsl\.py
* install **python** https://www.python.org/
  * **numpy** and **matplotlib** should be intalled
  * **scipy** is suggested to intall for reading *.sav from IDL in a demo
  * setting an environmental variable of `$PYTHONPATH` and placing fastqsl\.py into a private path is suggested, just append such line to ~/.bashrc
  ```bash
  export PYTHONPATH="$HOME/your/private/py/path:$PYTHONPATH"
  ```
-----------------------------
## Computation core with Fortran
### Compiler install 
* gfortran https://fortran-lang.org/learn/os_setup/install_gfortran/ or
* Intel® Fortran Compiler https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html
  * please append this line to ~/.bashrc
    ``` 
    source /opt/intel/oneapi/setvars.sh intel64
    ```
* Other compilers in https://fortran-lang.org/compilers/ should also work for FastQSL, while I have not test them, welcome for testing and sharing your experiences to me
* For checking whether your compiler is successfully installed, you can try https://github.com/el2718/sudoku
### Compilation
fields.f90 and trace_bline.f90 are included in fastqsl.f90, all *.mod produced from Compilation can be deleted

* For Linux and macOS (either by ifx/ifort or gfortran):
    ```bash
    ifx -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
    ```
    ```bash
    ifort -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
    ```
    ```bash
    gfortran -o fastqsl.x fastqsl.f90 -fopenmp -O3 -march=native
    ```
  * set -O3, -xHost, -ipo, -march=native for a better efficiency
* For Windows (either by ifort or gfortran):
  executing "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" in cmd first would be necessary
  ```bash
  ifort /o fastqsl.exe fastqsl.f90 /Qopenmp /O3 /QxHost /Qipo
  ``` 
  ```bash
  gfortran -o fastqsl.exe fastqsl.f90 -fopenmp -O3 -march=native
  ``` 
  * In Windows 11, also in some upgraded Windows 10, the pop-up window for fastqsl.exe can not be closed automatically, please uncomment this line in fastqsl.f90 to kill the pop-up window (delete !):
    ```
    ! call system('taskkill /im fastqsl.exe /f')
    ```
  * for ifx, the compilation should be the same as ifort, while I have not tested it
### Path of fastqsl.x
* please specify the path of fastqsl.x, 
  * in fastqsl\.pro, please correct the line of
    ```idl
    spawn, '/path/of/fastqsl.x'
    ```
  * in fastqsl\.py, please correct the line of  
    ```python
    subprocess.run(r'/path/of/fastqsl.x', shell=True)
    ```
* or move fastqsl.x to the `$PATH` (e.g. `/usr/local/bin/`) of the system and delete the text `/path/of/`
  * you can append this line to `~/.bashrc`
    ```
    export PATH=$HOME/bin:$PATH
    ```
    then `$HOME/bin` is a `$PATH` of the system
* For Windows, use fastqsl.exe instead of fastqsl.x
-----------------------------
## Parameters

The following introductions are wrote for fastqsl\.pro, the case is similar for fastqsl\.py, the tiny differences are elucidated. Please notice that the index order of an array in IDL and Fortran is Fortran-order (column-major), while the index order of an array in python is C-order (row-major), therefore the index order of an array for fastqsl\.py should be reversed.

The IDL language is case-insensitive, and the name of a keyword parameter can abbreviated to the shortest unambiguous string https://www.nv5geospatialsoftware.com/docs/Using_Keyword_Parameters.html . For example, **B_out** can be invoked by anyone of `,B_out=1`, `,/b_OUT`, `,B=1`, `,/b` to fastqsl\.pro. **Bx, By, Bz** are positional parameters (arguments) but not keyword parameters, therefore setting `,B=1` don't make unambiguous. Please do not apply these features to fastqsl\.py

### Magnetic field
  * **Bx, By, Bz**: 
    * For fastqsl\.pro, the typical dimesions of Bx, By and Bz are (nx,ny,nz). The dimesions of Bx can also be (3,nx,ny,nz), and then By and Bz should not be presented. For example, the last two lines of the following code produce identical images
      ```
      sz3d=size(By, /dim)
      nx=sz3d[0]
      ny=sz3d[1]
      nz=sz3d[2]
      Bvec=fltarr(3,nx,ny,nz)
      Bvec[0, *, *, *]=Bx
      Bvec[1, *, *, *]=By
      Bvec[2, *, *, *]=Bz
      fastqsl, Bx, By, Bz, /preview, fname='B3'
      fastqsl, Bvec,       /preview, fname='Bvec'
      ```
    * will be forcibly converted to 4-byte float arrays while writing 'bfield.bin'
    * some NaN values or magnetic nulls (where $\vec{B}=\vec{0}$) exist on grid is no matter for computation
### Coordinates

  * **xa, ya, za**: axis coordinates of magnetic field in 1D arrays
    * should not be presented If the **Bx, By, Bz** are set on a Cartesian uniformed grid
      ```
      stretchFlag= keyword_set(xa) and keyword_set(ya) and keyword_set(za)
      ```
    * the size of **xa, ya, za** must be consistant with the size of the 3D magnetic field
    * The values in **xa, ya, za** must be increased by order
  * **spherical**: whether the magnetic field is settled on a spherical grid. 
    * FastQSL uses longitude (in radian), latitude (in radian) and radius as the coordinates for spherical grid, denoted as $\varphi$, $\vartheta$, and $r$, respectively. The maximum range of $\varphi$ is $[0, 2 \pi]$, the maximum range of $\vartheta$ is $[-\pi/2, \pi/2]$. The relations between $\{\varphi, \vartheta, r\}$ and $\{x, y, z\}$ are
    $x=r \cos \vartheta \cos \varphi,$
    $y=r \cos \vartheta \sin \varphi,$
    $z=r \sin \vartheta.$
    * Default is 0 (Cartesian coordinates). If invoked, these keywords have such meanings:
      * **Bx, By, Bz**: longitudinal, latitudinal, radial component of the magnetic field, $B_\varphi, B_\vartheta, B_r$. 
        * **Be careful**: the index order of these arrays is `[i_longitude, i_latitude, i_radius]` for fastqsl\.pro, and is `[i_radius, i_latitude, i_longitude]` for fastqsl\.py.
      * **xa, ya, za**: axis coordinates of $\varphi, \vartheta, r$
      * **xreg, yreg, zreg**: output ranges of $\varphi, \vartheta, r$
        * will be rewrote as **lon_reg, lat_reg, r_reg** in returned **qsl**
      * **lon_delta, lat_delta, r_delta**: output grid spacing of $\varphi, \vartheta, r$
      * **arc_delta**: output grid spacing (in radian) of the arc on the great circle
          * works when csflag is invoked, the first curved axis is the arc on the great circle from point0 to point1, and the second axis is point0 -> point2, see appendix A3 of Chen (2026).
    * The classical spherical coordinates are $\{r, \theta, \varphi\}$, **which are not the spherical coordinates for FastQSL!**  The maximum range of $\theta$ is $[0, \pi]$. If you take a magnetic field on a grid with the classical spherical coordinates, two relations should be applied:
      * $\vartheta=\pi/2-\theta$
      * $B_\vartheta=-B_\theta$
    * For example, in a IDL script, if the arrays of $B_r(r, \theta, \varphi), B_\theta(r, \theta, \varphi), B_\varphi(r, \theta, \varphi), r, \theta, \varphi$ are `br, bp, bt, radius, theta, phi`, they should be converted to $B_\varphi(\varphi, \vartheta, r), B_\vartheta(\varphi, \vartheta, r), B_r(\varphi, \vartheta, r), \varphi, \vartheta, r$ in arrays of `B_lon, B_lat, B_r, lon_rad, lat_rad, radius`, the following code can give **q** at the bottom:
      ```
      b_lon  = reverse(transpose( bp, [2, 1, 0]), 2)
      b_lat  = reverse(transpose(-bt, [2, 1, 0]), 2)
      b_r    = reverse(transpose( br, [2, 1, 0]), 2)
      lon_rad= phi
      lat_rad= !pi/2. - reverse(theta)

      fastqsl, b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, /spherical, qsl=qsl
      ```
      In a python script, the corresponding code is
      ```python
      import numpy as np

      b_lon  = np.flip(np.transpose( bp, (2, 1, 0)), 1)
      b_lat  = np.flip(np.transpose(-bt, (2, 1, 0)), 1)
      b_r    = np.flip(np.transpose( br, (2, 1, 0)), 1)
      lon_rad= phi
      lat_rad= np.pi/2. - np.flip(theta)

      qsl=fastqsl(b_lon, b_lat, b_r, xa=lon_rad, ya=lat_rad, za=radius, spherical=True)
      ```
    * If you know Chinese, and have interest to see more properties of the coordinates, see the part of 例：经纬球坐标系 in 基矢量与张量\.pdf on https://github.com/el2718/thoughts/releases/tag/thoughts
  * **xperiod, yperiod, zperiod**: whether the $x, y, z$ -directions are periodical. For example, if **yperiod** is invoked, the period of y-direction is `ya[ny-1]-ya[0]` (if stretchFlag is invoked), the field line tracing will not stop at ymin, ymax. And the layer of `Bx[*,0,*]` will be copyed to `Bx[*,ny-1,*]` in fastqsl.x (also for `By, Bz`).
    * defaults are 0
    * can only work with Cartesian coordinates. For spherical coordinates, $\varphi$ -direction will be automatical periodical if `ya[ny-1]-ya[0]` is closed enough to $2 \pi$
    * Normally, FastQSL requires at least 3 layers for every direction to compute $\nabla \times \vec{B}$ and $\nabla \dfrac{\vec{B}}{B}$. If you need to analysis a field with the symmetry like 2D, e.g. $\dfrac{\partial \{B_x, B_y, B_z\}}{\partial y}=\{0,0,0\}$, just stack $\vec{B}(x,z)$ for 2 layers by
      ```
      size2d=size(Bx2d,/dim)
      nx=size2d(0)
      nz=size2d(1)
      Bx=fltarr(nx, 2, nz)
      By=fltarr(nx, 2, nz)
      Bz=fltarr(nx, 2, nz)
      for j=0, 1 do Bx[*,j,*]=Bx2d
      for j=0, 1 do By[*,j,*]=By2d
      for j=0, 1 do Bz[*,j,*]=Bz2d
      fastqsl, Bx, By, Bz, xreg=[0,nx-1], yreg=[0,0], zreg=[0,nz-1], qsl=qsl
      ```
      * Actually **yperiod** is invoked in the above code, the space period is 1
      * For spherical coordinates, ny, nz can not be 2, because $\vartheta, r$ -directions can not be periodical in any case
### Output domain
  * **xreg, yreg, zreg**: coordinates of the output region, in arrays of two elements
    * default is to include the whole 2D region at the bottom
    * If `(xreg[0] ne xreg[1]) and (yreg[0] ne yreg[1]) and (zreg[0] ne zreg[1]) and not csFlag`
      * the output domain is a box volume
      * invoke vFlag
    * If `(xreg[0] eq xreg[1]) or (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])`
      * the output domain is a cross section perpendicular to X or Y or Z axis
      * invoke cFlag
    * If `zreg=[0, 0]` (stretchFlag=0B) or `zreg=[za(0), za(0)]` (stretchFlag=1B)
      * the output domain is set at the bottom plane
      * invoke bflag
    * If csFlag is set, see below
      * invoke cFlag
    * in return, `qsl.xreg[1], qsl.yreg[1], qsl.zreg[1]` may slightly be changed, which is same as the final grid of `qsl.seed`, due to the trim at the final grid

  * **csFlag**: 
    * the output domain is the cross section is defined by three points:
      * the origin of output: 
        ```
        point0 = [xreg[0], yreg[0], zreg[0]] 
        ```
      * point0 -> point1 is the first axis, and 
        ```
        point1 = [xreg[1], yreg[1], zreg[0]]
        ```
      * point0 -> point2 is the second axis, and 
        ```
        point2 = [xreg[0], yreg[0], zreg[1]]
        ```  
    * `xreg[0], yreg[0], zreg[0]` does not have to be smaller than `xreg[1], yreg[1], zreg[1]` in this case
    * default is 0
  * **factor**: to bloat up the original resolution by a factor
    * For uniformed input grid, the grid spacing of output = 1/factor
    * default is 4
  * **delta**: grid spacing of output
    * For uniformed input grid, default is 1/factor
    * if keyword_set(delta), **factor** will be ignored
  * **lon_delta, lat_delta, r_delta, arc_delta**: see **spherical**
  * **seed**: launch points for tracing
    * if set `seed = 'original'` or `seed = 'original_bottom'` at the input; or **seed** is an array of coordinates with dimensions of (3) or (3,n1) or (3,n1,n2) or (3,n1,n2,n3), its units should be same as **xa, ya, za**
      * invoke sflag. And then FastQSL applied **the second way** to define the output domain by **seed**, and **xreg, yreg, zreg, csFlag, factor, delta, lon_delta, lat_delta, r_delta, arc_delta** will be ignored
      * if set `seed = 'original'`, then the output **seed** in **qsl** will be the original 3D grid of magnetic field. For example, for a field on an uniformed grid, if you only need $\nabla \times \vec{B}$ on the same input grid, just run
        ```
        IDL> fastqsl, Bx, By, Bz, seed='original', /CurlB, maxsteps=0, qsl=qsl
        ```
      * if set `seed = 'original_bottom'`, then the output **seed** in **qsl** will be the original 2D grid at the bottom of magnetic field
    * if set `seed = 1` at the input of fastqsl\.pro ( `, /seed` also makes seed eq 1), or set `seed = True` at the input of fastqsl\.py
      * not invoke sflag. The output grid is still described by **xreg, yreg, zreg, csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta**, and return as the array **seed** in **qsl**. For example, if **spherical** is invoked, and the output domain is set on the surface $\vartheta=\vartheta_0$ (i.e. `qsl.lat_reg[0] eq qsl.lat_reg[1]`), then `qsl.seed[0,i,j]` is same as `qsl.lon_reg[0] + i*qsl.lon_delta`, `qsl.seed[1,i,j]` is same as `qsl.lat_reg[0]`, and `qsl.seed[2,i,j]` is same as `qsl.r_reg[0] + j*qsl.r_delta`
### Tracing details
A magnetic field line is integrated by $\dfrac{\mathrm{d} \vec{r}(s)}{\mathrm{d} s}=\dfrac{\vec{B}}{B}$
  * **RK4Flag**:     to trace field line by RK4
    * default is 0 (use RKF45)
  * **step**:        step size of in tracing field lines ($\mathrm{d} s$) for RK4
    * default is 1.0
  * **tol**:         tolerance of a step in RKF45
    * default is 10.^(-4)
    * the unit of **step** and **tol** is the original grid spacing. This unit can vare from cell to cell within a stretched grid by a self-adaptive fashion with Equation (16-18) of [Zhang (2022)]((https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61))
  * **scottFlag**:  to integrate $\dfrac{\mathrm{d} \{\vec{U}(s),\vec{V}(s)\}}{\mathrm{d} s}=\{\vec{U}(s),\vec{V}(s)\} \cdot\nabla\dfrac{\vec{B}}{B}$ along with the field line tracing, and to give $Q$ and $Q_\perp$ by Equation (22) of [Scott_2017_ApJ_848_117](https://iopscience.iop.org/article/10.3847/1538-4357/aa8a64)
    * default is 0 (Method 3 of [Pariat_2012_A&A_541_A78](https://www.aanda.org/articles/aa/full_html/2012/05/aa18515-11/aa18515-11.html), some problematic sites are filled with [Scott (2017)]((https://iopscience.iop.org/article/10.3847/1538-4357/aa8a64)))
  * **inclineFlag**: to apply Equation (20) or (21) in [Zhang (2022)]((https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)):
  $\texttt{step}|_{S_0} = \textrm{max}([\texttt{step}_{\perp} \times |B_{n,0} / B|, \texttt{step}_\mathrm{min}])$
  $\texttt{tol}|_{S_0} = \texttt{tol}_{\perp} \times | B_{n,0} / B |^{1.5}$
for every field line launched from $S_0$, then **step** and **tol** at the input are actually $\texttt{step}_\perp$ and $\texttt{tol}_\perp$
    * invoke it can provide a better quality of **q** calculated with Method 3 of [Pariat (2012)](https://www.aanda.org/articles/aa/full_html/2012/05/aa18515-11/aa18515-11.html), while then FastQSL will takes slightly longer time for computation and makes slightly more points on **path** (then set a slightly larger **maxsteps** maybe necessary)
    * default is 0
    * can only works when scottFlag is not invoked
  * **maxsteps**:    maximum steps for stracing a field line at one direction
    * default is `4*(nx+ny+nz)` if traced by RKF45, is `long(4*(nx+ny+nz)/step)` if traced by RK4.
    * if we want field lines terminated at boundaries but not inside of the box, **maxsteps** should be large enough.
    * if **maxsteps** is too small, many 0 will appear in **rboundary**,  then **q** can not be given and results NaN, while **length, twist, q_perp** still have their values for one segment of field lines. For example, if `qsl.rboundary[i, j]` is 0 and `not stretchFlag and RK4flag and ~keyword_set(inclineFlag)`, then `qsl.length[i, j]` is approximately `2*maxsteps*step`
    * Sometimes we want disable the tracing and get **B, CurlB, seed** only, just run
      ```
      IDL> fastqsl, Bvec, /B, /CurlB, /seed, maxsteps=0, qsl=qsl
      ```
      * And then **q, rboundary, sign2d, tol, step, RK4Flag** will not exist in **qsl**. 
      * Even if **length_out, twist_out, rF_out, scottflag, path_out, loopB_out, loopCurlB_out** are set to 1 in the above command, they will be ignored, which means **length, twist, rFs, rFe, q_perp, path, loopB, loopCurlB** will not exist in **qsl**
  * **r_local**: the radius of the local sphere
    * default is 0.
    * If it is set to a positive value, **q_local** will be exported
### Efficiencicy and verbose
  * **nthreads**:    number of processors to engage the computation
    * default is OMP_GET_NUM_PROCS() - 2 (reset in fastqsl.x if the value is 0)
  * **silent**:      do not print anything if no mistake occured
    * default is 0
### Optional outputs
See **Products** for more details. All default values here are 0.
  * **B_out, CurlB_out, length_out, twist_out, path_out, loopB_out, loopCurlB_out**:  to export **B, CurlB, length, twist, path, loopB, loopCurlB**, respectively
    * **path_out** is not allowed to invoke with 3D output grid, due to the too huge memory occupation; if this happened, **path_out** will be ignored
    * **path_out** should be invoked first for invoking **loopB_out, loopCurlB_out**
  * **rF_out**:      to export **rFs, rFe**
  * **targetB_out**: to export **Bs, Be**
  * **targetCurlB_out**: to export **CurlBs, CurlBe**
### Output files
  * **odir**:        directory to save the results
    * default is `cdir+'fastqsl/'`, where cdir is the current directory
  * **fname**:        filename of the results
  * **save_file**: to produce `odir+fname+'.sav'` in fastqsl\.pro (`odir+fname+'.pkl'` in fastqsl\.py)
    * default is 0
  * **compress**:    to invoke the keyword compress of the IDL routine save, for saving storage
    * default is 0
    * not exist in fastqsl\.py
  * **preview**:  to produce PNG images for preview
    * default is 0. If you prefer the default value to be 1 in fastqsl\.pro, change this line
      ```
      preview         = keyword_set(preview) ; or n_elements(preview) eq 0
      ```
      to 
      ```
      preview         = keyword_set(preview) or n_elements(preview) eq 0
      ```   
  * **tmp_dir**:     the temporary directory for the data transmission between fastqsl.x and fastqsl\.pro
    * default is `cdir+'tmpFastQSL/'`
    * a virtual disk created with computer memory  provide an incredibly fast IO performance, set tmp_dir on such disk is suggested
      * For linux, /dev/shm/ is the directory from memory, one can set `tmp_dir='/dev/shm/tmpFastQSL/'`
      * For macOS, one way is detailed in https://lvv.me/posts/2025/09/25_ramdisk_on_macos/
      * For Windows, one choice is https://sourceforge.net/projects/imdisk-toolkit/
  * **keep_tmp**:    do not delete the temporary  binary files output from fastqsl.x
    * default is 0
-----------------------------
## Products
For fastqsl\.pro, the results is given by the structure **qsl**, can be returned by the keyword **qsl**, or can be saved as `odir+fname+'.sav'`. For example, the elements **q** can be accessed by `qsl.q`

For fastqsl\.py, the results is given by the [dictionary](https://docs.python.org/3.14/tutorial/datastructures.html#dictionaries) **qsl**, can be returned by the return of the function fastqsl, or can be saved as `odir+fname+'.pkl'`. For example, the element **q** can be accessed by `qsl['q']`

Possible elements in **qsl** are:
  * **csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta, RK4Flag, step, tol** can also appear, their meanings are the same as the input keywords
  * **seed**:    the coordinates of the output grid for the launch of tracing, its units are same as **xa, ya, za** if stretchFlag
  * **axis1**:  the coordinates $x, y$ ($\varphi, \vartheta$, if **spherical** is invoked) from point0 to point1
    * only appears when **csFlag** is invoked, then **axis1** is same as `qsl.seed[0:1, *, 0]`
  * **length**: $L= \int_\mathrm{path} \mathrm{d}l$, length of field lines launched from **seed**
  * **twist**: $T_w = \int_\mathrm{path} \dfrac{(\nabla \times \vec{B}) \cdot \vec{B}}{4\pi B^2} \mathrm{d}l$, can be used to measure how many turns two infinitesimally close field lines winding about each other. Eq. (16) of [Berger and Prior (2006) J. Phys. A: Math. Gen. 39 8321](https://iopscience.iop.org/article/10.1088/0305-4470/39/26/005); Also see [Liu_2016_ApJ_818_148](https://iopscience.iop.org/article/10.3847/0004-637X/818/2/148)
  * **q, q_perp**: squashing factor $Q$ and $Q_\perp$, see  [Titov_2002_JGRA_107_1164](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2001JA000278) and [Titov_2007_ApJ_660_863](https://iopscience.iop.org/article/10.1086/512671)
    * **q_perp** is only available when **scottFlag** invoked, [Pariat (2012)]((https://www.aanda.org/articles/aa/full_html/2012/05/aa18515-11/aa18515-11.html)) is not precise enough for $Q_\perp$
  * **q_local**: see Chen (2026), for locating where magnetic field lines bifurcate at, i.e.  (quasi-) separators 
  * **B, CurlB**:  $\vec{B}$, $\nabla \times \vec{B}$ on **seed**
    For example, sometimes we want to know the density, pressure, temperature distribution on a field line. The field lines is given by `*qsl.path[i]` from a previous run, and density, pressure, temperature are 3D arrays on the same grid of Bx, By, Bz. Then just run
    ```
    IDL> fastqsl, density, pressure, temperatrue, seed=*qsl.path[i], maxsteps=0, /B_out, qsl=qsl
    ```
    then `reform(qsl.B[0, *]), reform(qsl.B[1, *]), reform(qsl.B[2, *])` are actually the distributions of density, pressure, temperature on the field line
  * **sign2d**:   $\mathrm{sign}(B_z)|_{z=z_\mathrm{min}}$
    * only exist when the bottom plane is included
    * e.g. `slogq = alog10(qsl.q[*, *, 0] > 1.) * qsl.sign2d`
  * **rFs, rFe**:  coordinates of terminal foot points (r:remote, F:foot, s:start, e:end). A segment of a field line have two terminal points, at the start (or end) point, $\vec{B}$ (or $-\vec{B}$) points to the whole calculated path of the field line.  
    * If calculate at the bottom
      * if `qsl.sign2d[i, j]` is 1, then `qsl.rFs[*, i, j]` is identical to `qsl.seed[*, i, j]`, i.e.  the foot for launch; and `rFe[*, i, j]` is the target foot.
      * if `qsl.sign2d[i, j]` is -1, then `qsl.rFe[*, i, j]` is identical to `qsl.seed[*, i, j]`, and `rFs[*, i, j]` is the target foot.
      * If `qsl.sign2d[i, j]` is 0 and both `qsl.rFs[*, i, j]` and `qsl.rFe[*, i, j]` are not `qsl.seed[*, i, j]`, here must be on a bald patch ([Seehafer (1986)](https://link.springer.com/article/10.1007/BF00172044); [Titov (1993)](https://ui.adsabs.harvard.edu/abs/1993A%26A...276..564T/abstract) where $\vec{B}\cdot\nabla B_z|_\mathrm{PIL} > 0$
  * **rboundary**:  nature of the terminal points
    * rboundary is given by 10*rbs+rbe in fastqsl.x, therefore:
      ```
      rbs = qsl.rboundary / 10
      rbe = qsl.rboundary mod 10
      ```
      their values mark for where **rFs, rFe** are terminated:
      * 0 - inside the box. This value can appear in two cases:
        * the field line is longer than the length with the number of **maxsteps** can reach, then increasing **maxsteps** to a large enough number can terminate the field line at a boundary 
        * the field line is closed inside the box (e.g. a circular field line)
      * 1 - $z=z_\textrm{min}$ or $r=r_\textrm{min}$
      * 2 - $z=z_\textrm{max}$ or $r=r_\textrm{max}$
      * 3 - $y=y_\textrm{min}$ or $\vartheta=\vartheta_\textrm{min}$
      * 4 - $y=y_\textrm{max}$ or $\vartheta=\vartheta_\textrm{max}$
      * 5 - $x=x_\textrm{min}$ or $\varphi=\varphi_\textrm{min}$
      * 6 - $x=x_\textrm{max}$ or $\varphi=\varphi_\textrm{max}$
      * 7 - edge or corner
        * a special case is that the launch point is located outside the box
      * 8 - B is 0 or NaN
        * Strictly speaking, this is also included in the case of a field line terminated inside the box, but I want define this to a different category from 0

      So for a closed field line that both its two foots stand on the photosphere, its **rboundary** value is 11
    * If calculate at the bottom, 
      * rb_launch=1 for all launch points. i.e. if `qsl.sign2d[i,j]` is 1 (or -1), then `rbs[i, j]` (or `rbe[i, j]`) must be 1, and `rb_target[i, j] = rbe[i, j]` (or `rbs[i, j]`) for all target points
    * boundary_mark_colors.pdf is the color table for *_rbs.png, *_rbe.png *_rb_target.png
  * **Bs, Be**: $\vec{B}$ on **rFs, rFe**
      * [Priest and Demoulin (1995)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/95JA02740) use $N$ at the photosphere to locate QSLs; in [Titov (2007)](https://iopscience.iop.org/article/10.1086/512671),
      $Q = N^2 / \Delta$, and $\Delta =|B_\mathrm{n,launch}/B_\mathrm{n,target}|$ derived from $\nabla \cdot \vec{B}=0$
      If targetB_out is invoked, FastQSL can produce the image of $N$ and Bnr = $|B_\mathrm{n,launch}/B_\mathrm{n,target}|$ from **rboundary, Bs, Be**
  * **CurlBs, CurlBe**: $\nabla \times \vec{B}$  on **rFs, rFe**
  * **path**: path of field lines launched from **seed**
    * For example, if the output domain is 2D,
      * in fastqsl\.pro, `qsl.path` is a pointer array, `*qsl.path[i, j]` gives a field line with dimensions of (3, n), and
        * `(*qsl.path[i, j])[*, 0]` is identical to `qsl.rFs[*, i, j]`
        * `(*qsl.path[i, j])[*, n-1]` is identical to `qsl.rFe[*, i, j]`
        * `(*qsl.path[i, j])[*, qsl.index_seed[i, j]]` is identical to `qsl.seed[*, i, j]`
      * in fastqsl\.py, qsl['path'] is a list, qsl['path'][j][i] gives a field line with  dimensions of (n, 3), and
        * `qsl['path'][j][i][0, :]` is identical to `qsl['rFs'][j, i, :]`
        * `qsl['path'][j][i][-1, :]` is identical to `qsl['rFe'][j, i, :]`
        * `qsl['path'][j][i][qsl['index_seed'][j, i], :]` is identical to `qsl['seed'][j, i, :]`
    * Using RKF45 requires much less grid points a path than using RK4
    * A smaller tol (or step) requires more grid points on a path, while then the coordinate precision of the path is better
  * **loopB, loopCurlB**: $\vec{B}$, $\nabla \times \vec{B}$ on **path**
  * **index_seed**: the index in **path** for the launch points
-----------------------------
## Demos

### if use fastqsl\.pro
```
IDL> .r demo_charge4.pro
```
if you use Linux or macOS, and don't want to entry the interactive environment of IDL, you can create demo_charge4.sh with the content:
```bash
#! /bin/bash -f
idl <<EOF
.r demo_charge4.pro
EOF
```
and the user should have the execute permission of demo_charge4.sh:
```bash
chmod u+x demo_charge4.sh
```
then submit it in a terminal by
```bash
nohup ./demo_charge4.sh > verbose_demo.txt 2>&1 &
```
### if use fastqsl\.py
```bash
python3 demo_charge4.py
```
----------------------------
## Memory occupation 
In fastqsl.x, the most memory is occupied by: 
* a 3D magnetic field
* a 3D CurlB_field (if twistflag or CurlB_out or targetCurlB_out or loopCurlB_out)
* a dbdc_field (3 times as the occupation of the 3D magnetic field)
* data on a 2D slice
  * Even if the output domain is 3D, FastQSL processes the computaion layer by layer. Once a layer's computation is finished, the results are appended to associated *.bin files. The program then proceeds to the subsequent layer.
  * If **path_out** is invoked, **path** can occupied a quite large amount of memory

-----------------------------
## Derived Products
### Two Parameters Used for Modeling Solar Wind Speed

[Arge (2003)](https://pubs.aip.org/aip/acp/article-abstract/679/1/190/1010917/Improved-Method-for-Specifying-Solar-Wind-Speed?redirectedFrom=fulltext) found the Solar Wind Speed can be roughly modeled by
$v_\textrm{sw}=265+\dfrac{25}{f_s^{2/7}} \left(5-1.1\times \exp(1-(\theta_b/4)^2)\right)~\mathrm{km/s}, $
and two parameters in this formula are defined as:
* Magnetic field expansion factor
$
f_\mathrm{s}(\varphi, \vartheta, r)=
\left(\dfrac{R_0}{R_1}\right) ^2
\dfrac{B_r(\varphi_0, \vartheta_0, R_0)}
{B_r(\varphi_1, \vartheta_1, R_1)}, 
$
where $(\varphi_0, \vartheta_0, R_0)$ are the target coordinates traced from $(\varphi, \vartheta, r)$ to the inner boundary of $r=R_0$, and $(\varphi_1, \vartheta_1, R_1)$ are the target coordinates traced from $(\varphi, \vartheta, r)$ to the outer boundary of $r=R_1$.
* $\theta_b(\varphi, \vartheta, r)$, the minimum angular distance of an open-field footpoint from a coronal hole boundary.
* For a closed field line, its **rboundary** is 11, its $f_\mathrm{s}$ is set to 1000., its $\theta_b$ is set to 0., these default values can be adjusted in par2solarwind\.pro (par2solarwind\.py)


If you need this derived code, please visit https://github.com/el2718/par2solarwind
### Slip-Squashing Factors
Slip-squashing factors $Q_\mathrm{sf}$ and $Q_\mathrm{sb}$ ([Titov_2009_ApJ_693_1029](https://iopscience.iop.org/article/10.1088/0004-637X/693/1/1029)) are defined by two field line mappings and two boundary flow mappings between two instants; their large values define the surfaces that border of the reconnected or to-be-reconnected magnetic flux tubes for a given period of time during the magnetic evolution. 

For the case of static boundaries, we can compute the slip-squashing factors using the coordinate mapping provided by FastQSL. Following the initial coordinate mapping within the first magnetic field, the resulting mapped coordinates can be served as a seed grid for applying FastQSL to the second magnetic field. 

If you need this derived code, please visit https://github.com/el2718/slipq

-----------------------------
## History
* Jun 30, 2014 Rui Liu @ USTC, IDL edition
* Apr 21, 2015 Rui Liu and Jun Chen, deal with field lines pass through the boundary other than the bottom
* Apr 27, 2015 Rui Liu, calculate the squashing factor Q at a cross section
* Jun 15, 2015 Jun Chen, Fortran Edition; correct foot points with $\dfrac{\mathrm{d} \vec{r}}{\mathrm{d} r_i}=\dfrac{\vec{B}}{B_i}$
* Jul  8, 2015 Jun Chen, calculate the squashing factor Q in a box volume
* Oct 29, 2015 Jun Chen, deal with field lines touching the cut plane: use the plane quasi-perp to the field line  
* Nov  1, 2015 Jun Chen, fuse qcs and qfactor(z=0) in qfactor.f90
* Jun 22, 2016 Jun Chen, add the map of length
* Oct 30, 2017 Jun Chen, add the map of Bnr
* Aug 28, 2018 Jun Chen, supplement Q at maginal points
* May  1, 2021 PeiJin Zhang and Jun Chen, supplement RKF45 for tracing; previous classic RK4 is modified to 3/8-rule RK4
* May  1, 2021 Jun Chen, adapted to gfortran compiler
* Jun  1, 2021 Jun Chen, forcibly convert the input Bx, By, Bz to float arrays in IDL (Real(4) in Fortran)
* Jun 10, 2021 Jun Chen, add a keyword of maxsteps, suggested by Jiang, Chaowei
* Jun 11, 2021 Jun Chen, add the coordinates of mapping points to '*.sav' data, suggested by Jiang, Chaowei. 
* Jul  5, 2021 Jun Chen, switch the order of indexes of Bfield in trace_bline.f90 for a better efficiency
* Jul  9, 2021 Jun Chen, add a keyword scottFlag, add q_perp in output
* Dec 13, 2021 PeiJin Zhang and Jun Chen, adjust tol or step by incline
* Dec 25, 2021 Jun Chen, remove the reliance of Solar SoftWare
* Jan 30, 2022 Jun Chen, remove doppler_color_mix, due to the poor recognizability of green-white-yellow
* Feb 16, 2022 Jun Chen, reduce the memory occupation for 3D case in Fortran
* Apr 27, 2022 Jun Chen,
  * forcibly assign the minimum value of Q to 2, the theoretical minimum;
  * zreg[0] can be non-zero when calculate in a box volume;
  * the format of cut_str can be '(i0)' or '(f0.1)' or '(f0.2)' or '(f0.3)', according to the input
* Jun 10, 2022 Jun Chen, support stretched grid
* Oct 11, 2022 Jun Chen, support Windows
* Jan 17, 2023 Jun Chen, integrate doppler color in qfactor\.pro, doppler_color\.pro is not necessary anymore;
            for avioding an error in a remote server: % TVLCT: Unable to open X Windows display
* Jun 23, 2024 Jun Chen, change all txt files to binary files in tmp_dir, since a txt file could introduce errors to float values
* Jul  3, 2024 Jun Chen, support spherical coordinates
* Nov 11, 2024 Jun Chen, store output data in the structure QSL
* Nov 15, 2024 Jun Chen, add keywords of b_out, compress, redefine CurlB_out
* Nov 29, 2024 Jun Chen, add keywords of rF_out and length_out, rFs/rFe and length will not be saved by default
* Nov 30, 2024 Jun Chen, change the keyword twistflag to twist_out
* Dec  5, 2024 Jun Chen, add a keyword of seed
* Dec  6, 2024 Jun Chen, redefine rboundary
* Dec  8, 2024 Jun Chen, add keywords of path_out, loopB_out, loopCurlB_out
* Dec  9, 2024 Jun Chen, deal with B is NaN/0
* May  6, 2025 Jun Chen, deal the polar regions in spherical coordinates with two coordinates systems
* Jan, 6, 2026 Jun Chen,
  * add keywords of targetB_out, targetCurlB_out, qsl, silent, tmp_dir, keep_tmp
  * allow seed = 'original' or 'original_bottom'
  * remove the keyword of tmpB, RAMtmp
  * remove qsl.Bnr in output, this can be derived if targetB_out
  * rename the keyword nbridge to nthreads
  * rename qfactor to fastqsl
  * use os_sep=PATH_SEP() instead of '/'
  * kill the pop-up window for fastqsl.exe finally in Windows
* Jan, 12, 2026 Jun Chen, support python
* Jan, 13, 2026 Jun Chen, allow seed = 1 in fastqsl\.pro (True in fastqsl\.py) for exporting output grid
* Jan, 17, 2026 Jun Chen, remove the keyword no_preview, add a keyword of save_file, preview
* Jan, 20, 2026 Jun Chen, add a keyword of inclineFlag
* Jan, 24, 2026 Jun Chen, add keywords of xperiod, yperiod, zperiod
* Feb,  8, 2026 Jun Chen, add a keyword of r_local
