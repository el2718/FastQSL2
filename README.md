# FastQSL 2

To calculate the squashing factor Q, and other quatities relate to the magnetic connectivity, at the bottom or at a cross section or in a box volume or on some seeds, given a 3D magnetic field with Cartesian or spherical, uniform or stretched grid.

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This module is licensed under a [CC BY-NC-SA 4.0 License][cc-by-nc-sa].

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

## Cite as

* [Jun Chen*, Thomas Wiegelmann, Li Feng*, Bernhard Kliem, Chaowei Jiang, and Rui Liu. FastQSL 2: A Comprehensive Toolkit for Magnetic Connectivity Analysis.  2026, SCIENCE CHINA Physics, Mechanics & Astronomy, submitted](https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)

* [Zhang, P., Chen, J.*, Liu, R. and Wang, C., FastQSL: A Fast Computation Method for Quasi-separatrix Layers. 2022, The Astrophysical Journal, 937, 26](https://iopscience.iop.org/article/10.3847/1538-4357/ac8d61)

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
## Computation core with Fortran
### Fortran compiler install 
* gfortran https://fortran-lang.org/learn/os_setup/install_gfortran/ or
* ifort https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran
  * note: please append this line to ~/.bashrc   
  source /opt/intel/oneapi/setvars.sh intel64
### Compilation
* For Linux and MacOS (either by ifx/ifort or gfortran):
```bash
ifx -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
```
```bash
ifort -o fastqsl.x fastqsl.f90 -fopenmp -O3 -xHost -ipo
```
```bash
gfortran -o fastqsl.x fastqsl.f90 -fopenmp -O3 -march=native
```
* For Windows (either by ifort or gfortran)::
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
  * in fastqsl.pro, please correct the line of 
    spawn, '/path/of/fastqsl.x'
  * in fastqsl.py, please correct the line of  
    os.system(r'/path/of/fastqsl.x')
* or move fastqsl.x to the $PATH (e.g. /usr/local/bin/) of the system and delete the text /path/of/
* For Windows, use fastqsl.exe instead of fastqsl.x
-----------------------------
## Interface

### If use fastqsl.pro
* install IDL https://www.nv5geospatialsoftware.com/Products/IDL or
* install GDL https://gnudatalanguage.github.io
  * setting an environmental variable of GDL_PATH is necessary for write_png. 
    If you used one of the binary packages available for Linux, then it depends on the distribution:
    -Ubuntu & Fedora:  /usr/share/gnudatalanguage/lib
    -ArchLinux: /usr/lib/gdl
    -Gentoo: /usr/local/share/gdl
    -MacOS: /opt/local/share/gnudatalanguage/lib  
    **Please append such line to ~/.bashrc** (e.g. for Ubuntu)
    export GDL_PATH=/usr/share/gnudatalanguage/lib
### If use fastqsl.py
* install python https://www.python.org/
  * **numpy** and **matplotlib** should be intalled
-----------------------------
## Keywords

The following introductions are wrote for fastqsl.pro, the case is similar for fastqsl.py.
please notice that the array in IDL and Fortran uses Fortran-order (column-major), while the array in python uses C-order (row-major), therefore the index order of an array with fastqsl.py should be reversed.

### Magnetic field and its coordinates
  * **Bx, By, Bz**:  3D magnetic field, will be forcibly converted to float arrays while writing 'bfield.bin', 
    some NaN values or nulls exist on grid is no matter for computation
    * For fastqsl.pro, the normal dimesions of Bx, By and Bz are (nx,ny,nz). If the dimesions of Bx are (3,nx,ny,nz), then By and Bz should not be presented, and

          reform(Bx[0,*,*,*]) is the real Bx
          reform(Bx[1,*,*,*]) is the real By
          reform(Bx[2,*,*,*]) is the real Bz

  * **xa, ya, za**:  coordinates of grid in 1D arrays for streched grid. 
    * The size of x{yz}a must be consistant with the size of the 3D magnetic field
    * the values in x{yz}a should be increased by order
  * **sphericalFlag**: the field is set in spherical coordinatesdefault is 0BIf invoked, these keywords have such meanings:
    * Bx, By, Bz: longitudinal, latitudinal, radial component of the magnetic field,
      the indexe order of these arrays is [i_longitude, i_latitude, i_radius]
    * xa, ya, za: grid coordinates of longitude (radian), latitude (radian), radius
    * xreg, yreg, zreg: output ranges of longitude (radian), latitude (radian), radius,
      will be rewrote as QSL.lon_reg/lat_reg/r_reg in QSL
    * lon_delta, lat_delta: output grid spacing of longitude (radian), latitude (radian)
    * r_delta: output grid spacing of radius
    * If csflag invoked, the first curved axis is the arc on the great circle from point0 to point1, and the second axis is point0 -> point2
    * arc_delta: output grid spacing (radian) of the arc on the great circle
### Output domain
  * **xreg, yreg, zreg**: coordinates of the output region, in arrays of two elementsdefault is to include the whole 2D region at the bottom.
    * If (xreg[0] ne xreg[1]) AND (yreg[0] ne yreg[1]) AND (zreg[0] ne zreg[1]) AND NOT csFlag
      calculate Q in a box volume
      invoke vFlag
    * If (xreg[0] eq xreg[1]) OR (yreg[0] eq yreg[1]) or (zreg[0] eq zreg[1])
      calculate Q in a cross section perpendicular to X or Y or Z axis
      invoke cFlag
    * If zreg=[0, 0] (stretchFlag=0B) or zreg=[za(0), za(0)] (stretchFlag=1B), calculate Q at the bottom
      invoke bflag
    * If csFlag is set, see below
      nvoke cFlag
    * in output, QSL.x{yz}reg may be slightly changed, which is same as the final grid of qsl.surface, due to the trim at the final grid

  * **csFlag**: to calculate Q in a cross section defined by three pointsdefault is 0B

          point0=[xreg[0],yreg[0],zreg[0]] the origin of output
          point1=[xreg[1],yreg[1],zreg[0]] point0 -> point1, first axis
          point2=[xreg[0],yreg[0],zreg[1]] point0 -> point2, second axis

    * x{yz}reg[0] does not have to be smaller than x{yz}reg[1] in this case
  * **factor**: to bloat up the original resolution, i.e. grid spacing of output = 1/factordefault is 4
  * **delta**: grid spacing of outputdefault is 1/factor (when stretchFlag is 0B);
    * if keyword_set(delta), the keyword of factor will be ignored
  * **lon_delta, lat_delta, r_delta, arc_delta**: see sphericalFlag
  * **seeds**: launch points for tracingseeds is an array of coordinates with dimensions of (3) or (3,n1) or (3,n1,n2) or (3,n1,n2,n3), its units should be same as x{yz}a

    * if keyword_set(seeds), then:
      (1) ignore keywords of x{yz}reg, *delta, csflag, surface_out
      (2) if scottflag, QSL.q/q_perp will be saved
      (3) if seeds eq 'original', then seeds will be the original 3D grid
      (4) if seeds eq 'original_bottom', then seeds will be the original 2D grid at the bottom
             
    * for example, QSL.surface in a previous run can be seeds for a next run
### Optional outputs
  * **surface_out**: to save the coordinates of surface for output as QSL.surfaceonly works for 2D outputdefault is 0B
  * **B_out**:       to save magnetic field on the output grid as QSL.Bdefault is 0B
  * **CurlB_out**:   to save \nabla \times B on the output grid as QSL.curlBdefault is 0B
  * **length_out**:  to save length of field lines on the output grid as QSL.lengthdefault is 0B		
  * **twist_out**:   to save twist number Tw on the output grid as QSL.twistdefault is 0B
  * **rF_out**:      to save rFs/rFe on the output grid as QSL.rFs/rFedefault is 0B
  * **targetB_out/targetCurlB_out**: to save B/CurlB on rFs/rFe as QSL.Bs/Be/CurlBs/CurlBedefault is 0B
  * **path_out**:    to save QSL.path launched from the output grid, see following OUTPUTSdefault is 0B
  * **loopB_out/loopCurlB_out**: to save B/CurlB on QSL.pathpath_out should be invoked first for invoking themdefault is 0B
### Tracing details
  * **RK4Flag**:     to trace field line by RK4default is 0B (RKF45)
  * **step**:        step size in tracing field lines for RK4default is 1.0
  * **tol**:         tolerance of a step in RKF45default is 10.^(-4), if ~(scottFlag or keyword_set(seeds)), step and tol will be adjusted by Equations (20) (21) in Zhang (2022)
  * **maxsteps**:    maxium steps for stracing a field line at one directionsuggested by Jiang, Chaowei
    * default is 10*(nx+ny+nz) if traced by RKF45, or long(10*(nx+ny+nz)/step) if traced by RK4.
    * if we want field lines terminate at boundaries but not inside of the box, maxsteps should be large enough.
    * if maxsteps is too small, many 0 will appear in rboundary, 
      then q can not be given and results NaN, while length, twist, q_perp still have values.
    * sometimes we want see twist, q_perp of a segments of field line at a fixed length.
             If ~stretchFlag and RK4flag and (scottFlag or keyword_set(seeds)), and if rboundary(i,j) is 0, 
             then length(i,j) is approximately 2*maxsteps*step

    * sometimes we want get qsl.b/curlb/surface only immediately. If maxsteps is set to 0, then:
      (1) length_out, twist_out, rF_out, scottflag, path_out, loopB_out, loopCurlB_out will forcibly be 0B, which means length, twist, rFs, rFe, q_perp, path, loopB, loopCurlB will not be presented in output
      (2) q, rboundary, sign2d, tol, step, rk4flag will not be saved
      (3) if even B_out, CurlB_out, surface_out are not invoked, then do nothing and return directly
### Output details
  * **odir**:        directory to save the results		
  * **fstr**:        filename of the results
  * **no_preview**:  don't produce PNG images for previewdefault is 0B
  * **no_save**:     don't produce odir+fstr+'.sav'default is 0B
  * **compress**:    to invoke the keyword compress of save, for saving storagedefault is 0B
    * not exist in fastqsl.py
  * **tmp_dir**:     the temporary directory for the data transmission between fastqsl.x and fastqsl.prodefault is cdir+'tmpFastQSL/'

    * For linux, /dev/shm/ is a directory that utilizes shared memory, 
      setting tmp_dir='/dev/shm/tmpFastQSL/' can make read and write operations extremely fast.
      For macOS, https://lvv.me/posts/2025/09/25_ramdisk_on_macos/
      For Windows, one choice is https://sourceforge.net/projects/imdisk-toolkit/
  * **keep_tmp**:    keep binary file output from fastqsl.xdefault is 0B
### Rest details
  * **nthreads**:    number of processors to engagedefault is OMP_GET_NUM_PROCS() - 2 (set in fastqsl.x)
  * **silent**:      do not print anything if no mistake occureddefault is 0B
  * **scottFlag**: to calculate Q and Q_perp by the method of Scott_2017_ApJ_848_117
    default is 0B (method 3 of Pariat_2012_A&A_541_A78, some problematic sites are filled with Scott (2017))
-----------------------------
## Results
the results can be returned by the keyword QSL, or saved as odir+fstr+'.sav' (odir+fstr+'.pkl' with fastqsl.py). 
in the structure of QSL:
  * seeds, csFlag, delta, lon_delta, lat_delta, r_delta, arc_delta, scottFlag, RK4Flag, step, tol can also appear in QSL, the meaning are the same the input keywords
  * **surface**:     the coordinates of the 2D output grid for the launch of tracing, its units are same as x{yz}a
  * **arc**:         longitude (radian) and latitude (radian) of the arc on the great circle from point0 to point1.

          It appears when csflag and sphericalflag are invoked. Then arc = surface(0:1, *, 0)
  * **q**:        squashing factor Q (Titov_2002_JGRA_107_1164Titov_2007_ApJ_660_863)
  * **q_perp**:   q_perp in Titov_2007_ApJ_660_863
                  only available when scottFlag invoked, Pariat (2012) is not precise enough for q_perp
  * **twist**:    see Liu_2016_ApJ_818_148Eq. (16) of Berger and Prior (2006) J. Phys. A**: Math. Gen. 39 8321
                 T_w = \int (\curlB \cdot B)/(4\pi B^2) ds
  * **sign2d**:   sign(Bz) at the bottome.g. slogq = alog10(q > 1.)*sign2donly exist if the bottom plane is included
  * **length**:   length of field lines
  * **B/CurlB**:  B/curlB on the output grid
                 
          For example, sometimes we want to know the density/pressure/temperature distribution on a field line, 
          density/pressure/temperature is a 3D array on the same grid of Bx, By, Bz.
          the field lines is *QSL.path[i] from a previous run, then just run

          IDL> fastqsl, density, pressure, temperatrue, seeds=*QSL.path[i], maxstep=0, /B_out, qsl=qsl

          then reform(qsl.B[0,*])/reform(qsl.B[1,*])/reform(qsl.B[2,*]) is actually
          the density/pressure/temperature distribution on the field line

  * **rFs/rFe**:  coordinates of terminal foot points (r:remoteF:foots/e:start/end)suggested by Jiang, Chaowei

          A segment of a field line have two terminal pointsat the start/end point, B/-B points to the whole path. 
          If calculate at the bottom, if sign2d(i,j) is 1/-1, then rFe(*,i,j)/rFs(*,i,j) is the target foot, and rFs(*,i,j)/rFe(*,i,j) is the foot for launch, i.e. surface(*,i,j). 
          If sign2d(i,j) is 0, and if both rFs(*,i,j) and rFe(*,i,j) are not surface(*,i,j), here must be bald patch
  * **rboundary**:   nature of the terminal points, see 'subroutine trim_size' in trace_bline.f90

          rboundary is given by 10*rbs+rbe in fastqsl.x, therefore:
          rbs= rboundary/10, rbe= rboundary mod 10, their values mark for where rFs/rFe are terminated:

          0 - inside the domain
          1 - zmin, r_min
          2 - zmax, r_max
          3 - ymin, lat_min
          4 - ymax, lat_max
          5 - xmin, lon_min
          6 - xmax, lon_max
          7 - an edge or a corner
          8 - B is 0/NaN

          Then a closed field line is where rboundary is 11. 

          If calculate at the bottom, rb_launch=1 for all launch pointsi.e. where sign2d is 1/-1, then rbs/rbe must be 1. 
          Then we can define rb_target for all target points. Where sign2d is 1/-1, rb_target=rbe/rbs

  * **Bs/Be/CurlBs/CurlBe**: B/CurlB on QSL.rFs/QSL.rFe
  * **path**: paths of field lines through the output grid
  
          QSL.path is a pointer array, *QSL.path[i] gives a field line with size of 3xN, (*QSL.path[i])[*,j] is a grid on the field line

  * **loopB/loopCurlB**: B/CurlB on QSL.path
  * **index_launch**:    in index in path for launch points, 

          e.g. (*QSL.path[i,j])[*,QSL.index_launch[i,j]] is seeds[*,i,j] or surface[*,i,j]
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
In fastqsl.x: 
* a 3D magnetic field + 
* a 3D curlB_field (if twistflag or CurlB_out or targetCurlB_out or loopCurlB_out) + 
* some 2D arrays + 
* dbdc_field (3 times as the occupation of the 3D magnetic field)
-----------------------------
## History
* Jun 30, 2014 Rui Liu @ USTC, IDL edition
* Apr 21, 2015 Rui Liu and Jun Chen, deal with field lines pass through the boundary other than the bottom
* Apr 27, 2015 Rui Liu, calculate the squashing factor Q at a cross section
* Jun 15, 2015 Jun Chen, Fortran Edition; correct foot points with  d vp/d r_i=Bp/B_i
* Jul  8, 2015 Jun Chen, calculate the squashing factor Q in a box volume
* Oct 29, 2015 Jun Chen, deal with field lines touching the cut plane: use the plane quasi-perp to the field line  
* Nov  1, 2015 Jun Chen, fuse qcs and qfactor(z=0) in qfactor.f90
* Jun 22, 2016 Jun Chen, add the map of length
* Oct 30, 2017 Jun Chen, add the map of Bnr
* Aug 28, 2018 Jun Chen, supplement Q at maginal points
* May  1, 2021 Peijing Zhang and Jun Chen, supplement RKF45 for tracing; previous classic RK4 is modified to 3/8-rule RK4
* May  1, 2021 Jun Chen, adapted to gfortran compiler
* Jun  1, 2021 Jun Chen, forcibly convert the input Bx, By, Bz to float arrays in IDL (Real(4) in Fortran)
* Jun 10, 2021 Jun Chen, provide a keyword of maxsteps
* Jun 11, 2021 Jun Chen, add the coordinates of mapping points to '*.sav' data
* Jul  5, 2021 Jun Chen, switch the order of indexes of Bfield in trace_bline.f90 for a better efficiency
* Jul  9, 2021 Jun Chen, provide an option of the method of Scott_2017_ApJ_848_117, and q_perp as an output
* Dec 13, 2021 Peijing Zhang and Jun Chen, adjust tol or step by incline
* Dec 25, 2021 Jun Chen, remove the reliance of Solar SoftWare
* Jan 30, 2022 Jun Chen, remove doppler_color_mix, due to the poor recognizability of green-white-yellow
* Feb 16, 2022 Jun Chen, reduce the memory occupation for 3D case in Fortran
* Apr 27, 2022 Jun Chen,
  (1) forcibly assign the minimum value of Q to 2, the theoretical minimum;
  (2) zreg[0] can be non-zero when calculate in a box volume;
  (3) the format of cut_str can be '(i0)' or '(f0.1)' or '(f0.2)' or '(f0.3)', according to the input
* May 10, 2022 Jun Chen, determine nq1, nq2, nq3, bflag, cflag, vflag in Fortran
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
* Dec  5, 2024 Jun Chen, add a keyword of seeds
* Dec  6, 2024 Jun Chen, redefine rboundary
* Dec  8, 2024 Jun Chen, add keywords of path_out, loopB_out, loopCurlB_out
* Dec  9, 2024 Jun Chen, deal with B=NaN/0
* May  6, 2025 Jun Chen, deal the polar regions in spherical coordinates with two coordinates systems
* Jan, 6, 2026 Jun Chen,
  (1) add keywords of targetB_out, targetCurlB_out, qsl, silent, tmp_dir, keep_tmp
  (2) allow seeds = 'original' or 'original_bottom'
  (3) remove the keyword of tmpB, RAMtmp
  (4) remove qsl.Bnr in output, this can be derived if targetB_out
  (5) rename the keyword nbridge to nthreads
  (6) rename qfactor to fastqsl
  (7) use os_sep=PATH_SEP() instead of '/'
  (8) kill the window for fastqsl.exe finally in Windows
* Jan, 12, 2026 Jun Chen, support python
