# FastQSL 2

To calculate the squashing factor $Q$, and other quatities related to the magnetic connectivity, at the bottom or at a cross section or in a box volume or on some seed points, given a 3D magnetic field on a Cartesian or spherical, uniformed or stretched grid.

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This module is licensed under a [CC BY-NC-SA 4.0 License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

-----------------------------
## Cite as

* [Jun Chen*, Thomas Wiegelmann, Li Feng*, Bernhard Kliem, Chaowei Jiang, and Rui Liu. FastQSL 2: A Comprehensive Toolkit for Magnetic Connectivity Analysis.  2026, SCIENCE CHINA Physics, Mechanics & Astronomy, submitted](https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)

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
* install **IDL** https://www.nv5geospatialsoftware.com/Products/IDL or
* install **GDL** https://gnudatalanguage.github.io
  * setting an environmental variable of GDL_PATH is necessary for write_png. 
    If you used one of the binary packages available for Linux, then it depends on the distribution:
    -Ubuntu & Fedora:  /usr/share/gnudatalanguage/lib
    -ArchLinux: /usr/lib/gdl
    -Gentoo: /usr/local/share/gdl
    -MacOS: /opt/local/share/gnudatalanguage/lib  
    **Please append such line to ~/.bashrc**, e.g. for Ubuntu
    export GDL_PATH=/usr/share/gnudatalanguage/lib
### If use fastqsl\.py
* install **python** https://www.python.org/
  * **numpy** and **matplotlib** should be intalled
-----------------------------
## Computation core with Fortran
### Compiler install 
* gfortran https://fortran-lang.org/learn/os_setup/install_gfortran/ or
* Intel® Fortran Compiler https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler-download.html
  * please append this line to ~/.bashrc   
  source /opt/intel/oneapi/setvars.sh intel64
* Other compilers in https://fortran-lang.org/compilers/ should also work for FastQSL, while I have not test them, welcome for testing and sharing your experiences to me
* For checking whether your compiler is successfully installed, you can try https://github.com/el2718/sudoku
### Compilation
* For Linux and macOS (either by ifx/ifort or gfortran):
  * set -O3, -xHost, -ipo, -march=native for a better efficiency
```bash
ifx -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
```
```bash
ifort -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
```
```bash
gfortran -o fastqsl.x fastqsl.f90 -fopenmp -O3 -march=native
```

* For Windows (either by ifort or gfortran):
  * for ifx, the compilation should be the same as ifort, while I have not tested it

executing "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" in cmd first would be necessary
```bash
ifort /o fastqsl.exe fastqsl.f90 /Qopenmp /O3 /QxHost /Qipo
``` 
```bash
gfortran -o fastqsl.exe fastqsl.f90 -fopenmp -O3 -march=native
``` 
### Path of fastqsl.x
* please specify the path of fastqsl.x, 
  * in fastqsl\.pro, please correct the line of 
    spawn, '/path/of/fastqsl.x'
  * in fastqsl\.py, please correct the line of  
    os.system(r'/path/of/fastqsl.x')
* or move fastqsl.x to the $PATH (e.g. /usr/local/bin/) of the system and delete the text /path/of/
* For Windows, use fastqsl.exe instead of fastqsl.x
-----------------------------
## Keywords

The following introductions are wrote for fastqsl\.pro, the case is similar for fastqsl\.py, the tiny differences are elucidated. Please notice that the array in IDL and Fortran uses Fortran-order (column-major), while the array in python uses C-order (row-major), therefore the index order of an array for fastqsl\.py should be reversed.

### Magnetic field
  * **Bx, By, Bz**: 
    * For fastqsl\.pro, the normal dimesions of Bx, By and Bz are (nx,ny,nz). If the dimesions of Bx are (3,nx,ny,nz), then By and Bz should not be presented, and then
    reform(Bx[0, \*, \*, \*]) is the real $B_x$
    reform(Bx[1, \*, \*, \*]) is the real $B_y$
    reform(Bx[2, \*, \*, \*]) is the real $B_z$
    * will be forcibly converted to 4-byte float arrays while writing 'bfield.bin'
    * some NaN values or nulls exist on grid is no matter for computation
### Coordinates for stretched grid
  * **xa, ya, za**: axis coordinates of magnetic field in 1D arrays
    * should not be presented For uniformed input grid
    stretchFlag= keyword_set(xa) and keyword_set(ya) and keyword_set(za)
    * the size of x{yz}a must be consistant with the size of the 3D magnetic field
    * The values in x{yz}a should be increased by order
  * **sphericalFlag**: whether the magnetic field is settled on a spherical grid. 
    * FastQSL uses longitude (in radian), latitude (in radian) and radius as the coordinates for spherical grid, denoted as $\varphi$, $\vartheta$, and $r$, respectively. The maximum range of $\varphi$ is $[0,\, 2 \pi]$, the maximum range of $\vartheta$ is $[-\pi/2,\, \pi/2]$. The relations between $\{\varphi,\,\vartheta ,\,r\}$ and $\{x,\,y,\,z\}$ are
    $
    x=r\,\cos \vartheta\,\cos \varphi,\\
    y=r\,\cos \vartheta\,\sin \varphi,\\
    z=r\,\sin \vartheta.
    $
    The classical spherical coordinates are $\{r,\,\theta,\,\varphi\}$, **which are not the spherical coordinates for FastQSL!**  The maximum range of $\theta$ is $[0, \pi]$. If you take a magnetic field on a grid with the classical spherical coordinates, two relations should be applied:
      * $\vartheta=\pi/2-\theta$
      * $B_\vartheta=-B_\theta$
    * Default is 0 (Cartesian coordinates). If invoked, these keywords have such meanings:
      * Bx, By, Bz: longitudinal, latitudinal, radial component of the magnetic field, $B_\varphi, B_\vartheta, B_r$. 
        * **Be careful**: the index order of these arrays is [i_longitude, i_latitude, i_radius] for  fastqsl\.pro, and is [i_radius, i_latitude, i_longitude] for fastqsl\.py
      * xa, ya, za: axis coordinates of $\varphi,\,\vartheta ,\,r$
      * xreg, yreg, zreg: output ranges of $\varphi,\,\vartheta ,\,r$
        * will be rewrote as **lon_reg, lat_reg, r_reg** in returned **qsl**
      * lon_delta, lat_delta, r_delta: output grid spacing of $\varphi,\,\vartheta ,\,r$
      * arc_delta: output grid spacing (in radian) of the arc on the great circle
          * works when csflag is invoked, the first curved axis is the arc on the great circle from point0 to point1, and the second axis is point0 -> point2
### Output domain
  * **xreg, yreg, zreg**: coordinates of the output region, in arrays of two elements
    * default is to include the whole 2D region at the bottom
    * If (xreg[0] ne xreg[1]) and (yreg[0] ne yreg[1]) and (zreg[0] ne zreg[1]) and not csFlag
      * the output domain is a box volume
      * invoke vFlag
    * If (xreg[0] eq xreg[1]) or (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])
      * the output domain is a cross section perpendicular to X or Y or Z axis
      * invoke cFlag
    * If zreg=[0, 0] (stretchFlag=0B) or zreg=[za(0), za(0)] (stretchFlag=1B)
      * the output domain is set at the bottom plane
      * invoke bflag
    * If csFlag is set, see below
      * invoke cFlag
    * in return, qsl.x{yz}reg[1] may be slightly changed, which is same as the final grid of qsl.seed, due to the trim at the final grid

  * **csFlag**: 
    * the output domain is the cross section is defined by three points:
      * the origin of output: 
        point0 = [xreg[0], yreg[0], zreg[0]] 
      * point0 -> point1 is the first axis, and 
        point1 = [xreg[1], yreg[1], zreg[0]]
      * point0 -> point2 is the second axis, and 
        point2 = [xreg[0], yreg[0], zreg[1]]  
    * default is 0
    * x{yz}reg[0] does not have to be smaller than x{yz}reg[1] in this case
  * **factor**: to bloat up the original resolution
    * For uniformed input grid, grid spacing of output = 1/factor
    * default is 4
  * **delta**: grid spacing of output
    * For uniformed input grid, default is 1/factor
    * if keyword_set(delta), the keyword of factor will be ignored
  * **lon_delta, lat_delta, r_delta, arc_delta**: see **sphericalFlag**
  * **seed**: launch points for tracing
    * if seed is 'original', or 'original_bottom'; or seed is an array of coordinates with dimensions of (3) or (3,n1) or (3,n1,n2) or (3,n1,n2,n3), its units should be same as x{yz}a
      * invoke sflag
      * if seed is 'original', then seed will be the original 3D grid
      * if seed is 'original_bottom', then seed will be the original 2D grid at the bottom
      * x{yz}reg, csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta will be ignored
      * if scottFlag, **q, q_perp** will be exported; otherwise even **q** can not be found
    * if seed eq 1 in fastqsl\.pro ( /seed also makes seed eq 1), or seed is True in fastqsl\.py
      * not invoke sflag. The output grid is still described by x{yz}reg, csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta, and return as the array **seed** in **qsl**
### Tracing details
  * **RK4Flag**:     to trace field line by RK4
    * default is 0 (use RKF45)
  * **step**:        step size in tracing field lines for RK4
    * default is 1.0
  * **tol**:         tolerance of a step in RKF45
    * default is 10.^(-4)
    * if not (scottFlag or keyword_set(seed)), step and tol will be adjusted by Equations (20) (21) in Zhang (2022)
  * **maxsteps**:    maximum steps for stracing a field line at one direction; suggested by Jiang, Chaowei
    * default is 10*(nx+ny+nz) if traced by RKF45, or long(10*(nx+ny+nz)/step) if traced by RK4.
    * if we want field lines terminated at boundaries but not inside of the box, maxsteps should be large enough.
    * if maxsteps is too small, many 0 will appear in **rboundary**,  then **q** can not be given and results NaN, while **length, twist, q_perp** still have values.
    * sometimes we want see **twist, q_perp** of a segments of field line at a fixed length.
             If not stretchFlag and RK4flag and (scottFlag or keyword_set(seed)), and if rboundary[i, j] is 0, 
             then length[i, j] is approximately 2\*maxsteps\*step

    * Sometimes we want get **B, CurlB, seed** only, and immediately. If maxsteps is set to 0, then:
      * length_out, twist_out, rF_out, scottflag, path_out, loopB_out, loopCurlB_out will forcibly be 0, which means **length, twist, rFs, rFe, q_perp, path, loopB, loopCurlB** will not be presented in **qsl**
      * **q, rboundary, sign2d, tol, step, RK4Flag** will not be presented in **qsl**
      * if even B_out, CurlB_out, launch_out are not invoked, then do nothing and return directly
### Output details
  * **odir**:        directory to save the results
    * default is cdir+'fastqsl/', where cdir is the current directory
  * **fstr**:        filename of the results
  * **preview**:     produce PNG images for preview
    * default is 0
  * **save_file**:   produce odir+fstr+'.sav' in fastqsl\.pro (odir+fstr+'.pkl' in fastqsl\.py)
    * default is 0
  * **compress**:    to invoke the keyword compress of save, for saving storage
    * default is 0
    * not exist in fastqsl\.py
  * **tmp_dir**:     the temporary directory for the data transmission between fastqsl.x and fastqsl\.pro
    * default is cdir+'tmpFastQSL/'
    * a virtual disk created with computer memory  provide an incredibly fast IO performance, set tmp_dir on such disk is suggested
      * For linux, /dev/shm/ is the directory from memory, one can set tmp_dir='/dev/shm/tmpFastQSL/'
      * For macOS, one way is detailed in https://lvv.me/posts/2025/09/25_ramdisk_on_macos/
      * For Windows, one choice is https://sourceforge.net/projects/imdisk-toolkit/
  * **keep_tmp**:    do not delete the temporary  binary files output from fastqsl.x
    * default is 0
### Rest details
  * **nthreads**:    number of processors to engage the computation
    * default is OMP_GET_NUM_PROCS() - 2 (set in fastqsl.x)
  * **silent**:      do not print anything if no mistake occured
    * default is 0
  * **scottFlag**:   to calculate $Q$ and $Q_\perp$ by the method of [Scott_2017_ApJ_848_117](https://iopscience.iop.org/article/10.3847/1538-4357/aa8a64)
    * default is 0 (method 3 of [Pariat_2012_A&A_541_A78](https://www.aanda.org/articles/aa/full_html/2012/05/aa18515-11/aa18515-11.html), some problematic sites are filled with Scott (2017))
### Optional outputs
See **Products** for more details. All default values here are 0. 
  * **B_out, CurlB_out, length_out, twist_out, path_out, loopB_out, loopCurlB_out**:  to export **B, CurlB, length, twist, path, loopB, loopCurlB**, respectively
    * path_out should be invoked first for invoking loopB_out, loopCurlB_out
  * **rF_out**:      to export **rFs, rFe**
  * **targetB_out**: to export **Bs, Be**
  * **targetCurlB_out**: to export **CurlBs, CurlBe**
-----------------------------
## Products
For fastqsl\.pro, the results is given by the structure **qsl**, can be returned by the keyword qsl, or can be saved as odir+fstr+'.sav'. For example, the elements **q** can be accessed by qsl.q.

For fastqsl\.py, the results is given by the dictionary **qsl**, can be returned by the return of the function fastqsl, or can be saved as odir+fstr+'.pkl'. For example, the element **q** can be accessed by qsl['q']

Possible elements in **qsl** are:
  * **csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta, RK4Flag, step, tol** can also appear, their meanings are the same as the input keywords
  * **seed**:    the coordinates of the output grid for the launch of tracing, its units are same as x{yz}a if stretchFlag
  * **arc**:      longitude (in radian) and latitude (in radian) of the arc on the great circle from point0 to point1. It appears when csflag and sphericalflag are invoked, so then arc is seed[0:1, \*, 0]
  * **q**:        squashing factor $Q$, see  [Titov_2002_JGRA_107_1164](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2001JA000278) and [Titov_2007_ApJ_660_863](https://iopscience.iop.org/article/10.1086/512671)
  * **q_perp**:   $Q_\perp$ in Titov (2007)
    * only available when scottFlag invoked, Pariat (2012) is not precise enough for $Q_\perp$
  * **twist**: $T_w = \int_L (\nabla \times \vec{B}) \cdot \vec{B}/(4\pi B^2)\, \mathrm{d}l$, can be used to measure how many turns two infinitesimally close field lines winding about each other. Eq. (16) of [Berger and Prior (2006) J. Phys. A: Math. Gen. 39 8321](https://iopscience.iop.org/article/10.1088/0305-4470/39/26/005); Also see [Liu_2016_ApJ_818_148](https://iopscience.iop.org/article/10.3847/0004-637X/818/2/148).

  * **sign2d**:   sign(Bz) at the bottom
    * only exist when the bottom plane is included
    * e.g. slogq = alog10(q > 1.)*sign2d
  * **length**:   length of field lines
  * **B, CurlB**:  $\vec{B}$, $\nabla \times \vec{B}$ on the output grid
    For example, sometimes we want to know the density, pressure, temperature distribution on a field line. The field lines is given by *qsl.path[i] from a previous run, and density, pressure, temperature are 3D arrays on the same grid of Bx, By, Bz. Then just run

```bash
IDL> fastqsl, density, pressure, temperatrue, seed=*qsl.path[i], maxsteps=0, /B_out, qsl=qsl
```

then reform(qsl.B[0, \*]), reform(qsl.B[1, \*]), reform(qsl.B[2, \*]) are actually the density, pressure, temperature distribution on the field line

  * **rFs, rFe**:  coordinates of terminal foot points (r:remote, F:foot, s:start, e:end), suggested by Jiang, Chaowei. A segment of a field line have two terminal points, at the start (or end) point, $\vec{B}$ (or $-\vec{B}$) points to the whole calculated path of the field line.  
    * If calculate at the bottom
      * if sign2d[i, j] is 1, then rFs[\*, i, j] is the foot for launch, i.e. seed[\*, i, j]; and rFe[\*, i, j] is the target foot.
      * if sign2d[i, j] is -1, then rFe[\*, i, j] is the foot for launch, and rFs[\*, i, j] is the target foot.
      * If sign2d[i, j] is 0 and both rFs[\*, i, j] and rFe[\*, i, j] are not seed[\*, i, j], here must be on a bald patch where 
      $\vec{B}\cdot\nabla B_z|_\mathrm{PIL} > 0$.
  * **rboundary**:  nature of the terminal points, see 'subroutine trim_size' in trace_bline.f90. 
  
    * rboundary is given by 10*rbs+rbe in fastqsl.x, therefore:
       rbs = rboundary / 10
       rbe = rboundary mod 10
       their values mark for where rFs, rFe are terminated:
      * 0 - inside the domain
      * 1 - zmin, r_min
      * 2 - zmax, r_max
      * 3 - ymin, lat_min
      * 4 - ymax, lat_max
      * 5 - xmin, lon_min
      * 6 - xmax, lon_max
      * 7 - an edge or a corner
      * 8 - B is 0 or NaN
    So for a closed field line that its both two foots stand on the photosphere, its rboundary is 11. 
    * If calculate at the bottom, rb_launch=1 for all launch points. i.e. where sign2d is 1 (or -1), then rbs (or rbe) must be 1. Then we can define rb_target for all target points. Where sign2d is 1 (or -1), rb_target=rbe (or rbs)
    * boundary_mark_colors.pdf is the color table for *_rbs.png, *_rbe.png *_rb_target.png.
  * **Bs, Be**: $\vec{B}$ on rFs, rFe
      * [Priest and Demoulin (1995)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/95JA02740) use $N$ at the photosphere to locate QSLs; in Titov (2007),
      $Q = N^2 / \Delta$, and $\Delta =|B_\mathrm{n,launch}/B_\mathrm{n,target}|$ derived from $\nabla \cdot \vec{B}=0$.
      If targetB_out is invoked, FastQSL can produce the image of $N$ and Bnr = $|B_\mathrm{n,launch}/B_\mathrm{n,target}|$ from rboundary, Bs and Be.
  * **CurlBs, CurlBe**: $\nabla \times \vec{B}$  on rFs, rFe
  * **path**: path of field lines launched from the output grid. 
    * For example, if the output domain is 2D,
      * In fastqsl\.pro, qsl.path is a pointer array, *qsl.path[i, j] gives a field line with size of 3xN, (\*qsl.path[i, j])[\*, k] is a grid on the field line
      * In fastqsl\.py, qsl['path'] is a list, qsl['path'][j][i] gives a field line with size of 3xN, qsl['path'][j][i][k, :] is a grid on the field line
    * using RKF45 requires much less grid points a path than using RK4; and a smaller tol (or step) requires more grid points on a path, while then the coordinate precision of the path is better
  * **loopB, loopCurlB**: $\vec{B}$, $\nabla \times \vec{B}$ on path
  * **index_launch**:    in index in path for launch points, e.g. (\*qsl.path[i, j])[\*, qsl.index_launch[i, j]] is qsl.seed[\*, i, j]
-----------------------------
## Derived Products
### Two Parameters in WSA Model
[Wang-Sheeley-Arge (WSA) model](https://pubs.aip.org/aip/acp/article-abstract/679/1/190/1010917/Improved-Method-for-Specifying-Solar-Wind-Speed?redirectedFrom=fulltext) is a famous model for solar wind speed upgraded from [Sheeley-Arge (WS) model](https://ui.adsabs.harvard.edu/abs/1990ApJ...355..726W/abstract), two parameters in WSA model  can be derived from the coordinate mapping of FastQSL:
* Magnetic field expansion factor
$
f_\mathrm{s}(\varphi,\,\vartheta,\,r)=
\left(\dfrac{R_\textrm{sun}}{R_1}\right) ^2
\dfrac{B_r(\varphi_\textrm{sun},\,\vartheta_\textrm{sun},\,R_\textrm{sun})}
{B_r(\varphi_1,\,\vartheta_1,\,R_1)}, 
$
where $(\varphi_\textrm{sun},\,\vartheta_\textrm{sun},\,R_\textrm{sun})$ are the target coordinates traced from $(\varphi,\,\vartheta,\,r)$ to the inner boundary of $r=R_\textrm{sun}$, and $(\varphi_1,\,\vartheta_1,\,R_1)$ are the target coordinates traced from $(\varphi,\,\vartheta,\,r)$ to the outer boundary of $r=R_1$.
* $\theta_b(\varphi,\,\vartheta,\,r)$, the minimum angular distance of an open-field footpoint from a coronal hole boundary.
* For a closed field line, its rboundary is 11, its $f_\mathrm{s}$ is set to 1000., its $\theta_b$ is set to 0., these default values can be adjusted in was_par\.pro (was_par\.py)

### Slip-Squashing Factors
Slip-squashing factors $Q_\mathrm{sf}$ and $Q_\mathrm{sb}$ ([Titov_2009_ApJ_693_1029](https://iopscience.iop.org/article/10.1088/0004-637X/693/1/1029)) are defined by two field line mappings and two boundary flow mappings between two instants; their large values define the surfaces that border of the reconnected or to-be-reconnected magnetic flux tubes for a given period of time during the magnetic evolution. 

For the case of static boundaries, we can compute the slip-squashing factors using the coordinate mapping provided by FastQSL. Following the initial coordinate mapping within the first magnetic field, the resulting mapped coordinates can be served as a seed grid for applying FastQSL in the second magnetic field. 

-----------------------------
## Demos
```bash
IDL> .r demo_charge4.pro
```
or
```bash
python3 demo_charge4.py
```
----------------------------
## Memory occupation 
In fastqsl.x, the most memory is occupied by: 
* a 3D magnetic field
* a 3D CurlB_field (if twistflag or CurlB_out or targetCurlB_out or loopCurlB_out)
* some 2D arrays
* dbdc_field (3 times as the occupation of the 3D magnetic field)
* some field lines (if path_out)
-----------------------------
## History
* Jun 30, 2014 Rui Liu @ USTC, IDL edition
* Apr 21, 2015 Rui Liu and Jun Chen, deal with field lines pass through the boundary other than the bottom
* Apr 27, 2015 Rui Liu, calculate the squashing factor Q at a cross section
* Jun 15, 2015 Jun Chen, Fortran Edition; correct foot points with $\dfrac{\mathrm{d}\, \vec{r}}{\mathrm{d} r_i}=\dfrac{\vec{B}}{B_i}$
* Jul  8, 2015 Jun Chen, calculate the squashing factor Q in a box volume
* Oct 29, 2015 Jun Chen, deal with field lines touching the cut plane: use the plane quasi-perp to the field line  
* Nov  1, 2015 Jun Chen, fuse qcs and qfactor(z=0) in qfactor.f90
* Jun 22, 2016 Jun Chen, add the map of length
* Oct 30, 2017 Jun Chen, add the map of Bnr
* Aug 28, 2018 Jun Chen, supplement Q at maginal points
* May  1, 2021 PeiJin Zhang and Jun Chen, supplement RKF45 for tracing; previous classic RK4 is modified to 3/8-rule RK4
* May  1, 2021 Jun Chen, adapted to gfortran compiler
* Jun  1, 2021 Jun Chen, forcibly convert the input Bx, By, Bz to float arrays in IDL (Real(4) in Fortran)
* Jun 10, 2021 Jun Chen, provide a keyword of maxsteps
* Jun 11, 2021 Jun Chen, add the coordinates of mapping points to '*.sav' data
* Jul  5, 2021 Jun Chen, switch the order of indexes of Bfield in trace_bline.f90 for a better efficiency
* Jul  9, 2021 Jun Chen, provide an option of the method of Scott_2017_ApJ_848_117, and q_perp as an output
* Dec 13, 2021 PeiJin Zhang and Jun Chen, adjust tol or step by incline
* Dec 25, 2021 Jun Chen, remove the reliance of Solar SoftWare
* Jan 30, 2022 Jun Chen, remove doppler_color_mix, due to the poor recognizability of green-white-yellow
* Feb 16, 2022 Jun Chen, reduce the memory occupation for 3D case in Fortran
* Apr 27, 2022 Jun Chen,
  (1) forcibly assign the minimum value of Q to 2, the theoretical minimum;
  (2) zreg[0] can be non-zero when calculate in a box volume;
  (3) the format of cut_str can be '(i0)' or '(f0.1)' or '(f0.2)' or '(f0.3)', according to the input
* Jun 10, 2022 Jun Chen, support stretched grid
* Oct 11, 2022 Jun Chen, support Windows
* Jan 17, 2023 Jun Chen, integrate doppler color in qfactor.pro, doppler_color.pro is not more necessary;
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
* Dec  9, 2024 Jun Chen, deal with B=NaN/0
* May  6, 2025 Jun Chen, deal the polar regions in spherical coordinates with two coordinates systems
* Jan, 6, 2026 Jun Chen,
  (1) add keywords of targetB_out, targetCurlB_out, qsl, silent, tmp_dir, keep_tmp
  (2) allow seed = 'original' or 'original_bottom'
  (3) remove the keyword of tmpB, RAMtmp
  (4) remove qsl.Bnr in output, this can be derived if targetB_out
  (5) rename the keyword nbridge to nthreads
  (6) rename qfactor to fastqsl
  (7) use os_sep=PATH_SEP() instead of '/'
  (8) kill the pop-up window for fastqsl.exe finally in Windows
* Jan, 12, 2026 Jun Chen, support python
* Jan, 13, 2026 Jun Chen, allow seed = 1 in fastqsl\.pro (True in fastqsl\.py) for exporting output grid
* Jan, 16, 2026 Jun Chen, remove the key_word  of no_preview, add key_words of preview, save_file
