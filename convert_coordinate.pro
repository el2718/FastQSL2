function convert_coordinate, coordinate, v1, v2, v3, v4, $
v1out=v1out, v2out=v2out, v3out=v3out, v4out=v4out, $
mode=mode, tmp_dir=tmp_dir, nthreads=nthreads
;-----------------------------------------------------
if ~keyword_set(nthreads) then nthreads=0
if ~keyword_set(mode) then mode=0
if size(mode,/tname) eq 'STRING' then begin
    case mode of
		'xyz_to_lon_lat_r': mode=0
        'lon_lat_r_to_xyz': mode=1
		'xyz_to_lon2_lat2_r': mode=2
        'lon2_lat2_r_to_xyz': mode=3
        'lon_lat_r_to_lon2_lat2_r': mode= 4
        'lon2_lat2_r_to_lon_lat_r': mode= 5
    endcase
    if size(mode,/tname) eq 'STRING' then message, 'Something is wrong with mode'
endif
;-----------------------------------------------------
sz_coor=size(coordinate)
if sz_coor[1] ne 3 then message, 'Something is wrong with coordinate'
ndata = sz_coor[n_elements(sz_coor)-1]

r4flag = size(coordinate,/tname) ne 'DOUBLE'
if (r4flag and size(coordinate,/tname) ne 'FLOAT') then coordinate=float(coordinate)
;-----------------------------------------------------
present1= N_PARAMS() ge 2
present2= N_PARAMS() ge 3
present3= N_PARAMS() ge 4
present4= N_PARAMS() ge 5

if present1 then begin
    sz_v=size(v1)
    if sz_v[0] ne sz_coor[0] then message, 'Something is wrong with v1'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v1'
    if size(v1,/tname) ne size(coordinate,/tname) then begin
        if r4flag then v1=float(v1) else v1=double(v1)
    endif
endif

if present2 then begin
    sz_v=size(v2)
    if sz_v[0] ne sz_coor[0] then message, 'Something is wrong with v2'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v2'
    if size(v2,/tname) ne size(coordinate,/tname) then begin
        if r4flag then v2=float(v2) else v2=double(v2)
    endif
endif

if present3 then begin
    sz_v=size(v3)
    if sz_v[0] ne sz_coor[0] then message, 'Something is wrong with v3'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v3'
    if size(v3,/tname) ne size(coordinate,/tname) then begin
        if r4flag then v3=float(v3) else v3=double(v3)
    endif
endif

if present4 then begin
    sz_v=size(v4)
    if sz_v[0] ne sz_coor[0] then message, 'Something is wrong with v4'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v4'
    if size(v4,/tname) ne size(coordinate,/tname) then begin
        if r4flag then v4=float(v4) else v4=double(v4)
    endif
endif
;-----------------------------------------------------
os_sep=PATH_SEP()

cd, current = cdir
IF STRMID(cdir, STRLEN(cdir)-1) NE os_sep THEN cdir=cdir+os_sep

if keyword_set(tmp_dir) then begin
	IF STRMID(tmp_dir, STRLEN(tmp_dir)-1) NE os_sep THEN tmp_dir=tmp_dir+os_sep
endif else tmp_dir= cdir+'tmpFastQSL'+os_sep

old_tmp_dir=file_test(tmp_dir)
if ~old_tmp_dir then file_mkdir, tmp_dir
;-----------------------------------------------------
get_lun, unit
openw,  unit, tmp_dir+ 'head.bin'
writeu, unit, long([mode, nthreads, r4flag]), long64(ndata) 
close,  unit

openw,  unit, tmp_dir+'coordinate.bin'
writeu, unit, coordinate
close,  unit

if present1 then begin
openw,  unit, tmp_dir+'v1.bin'
writeu, unit, v1
close,  unit
endif

if present2 then begin
openw,  unit, tmp_dir+'v2.bin'
writeu, unit, v2
close,  unit
endif

if present3 then begin
openw,  unit, tmp_dir+'v3.bin'
writeu, unit, v3
close,  unit
endif

if present4 then begin
openw,  unit, tmp_dir+'v4.bin'
writeu, unit, v4
close,  unit
endif
;-----------------------------------------------------
cd, tmp_dir
; please specify the path
; spawn, '/path/of/convert_coordinate.x'
spawn, '~/Desktop/QSLS/update/convert_coordinate.x'
cd, cdir
;-----------------------------------------------------
openr,  unit, tmp_dir+'coordinate_out.bin'
coordinate_out=coordinate
readu, unit, coordinate_out
close,  unit

if present1 then begin
openw,  unit, tmp_dir+'v1out.bin'
v1out=v1
writeu, unit, v1out
close,  unit
endif

if present2 then begin
openw,  unit, tmp_dir+'v2out.bin'
v2out=v2
writeu, unit, v2out
close,  unit
endif

if present3 then begin
openw,  unit, tmp_dir+'v3out.bin'
v3out=v3
writeu, unit, v3out
close,  unit
endif

if present4 then begin
openw,  unit, tmp_dir+'v4out.bin'
v4out=v4
writeu, unit, v4out
close,  unit
endif
;-----------------------------------------------------
free_lun, unit, /force

if old_tmp_dir then begin
    file_delete, tmp_dir+['head.bin','coordinate.bin','coordinate_out.bin']
    if present1 then file_delete, tmp_dir+['v1.bin', 'v1out.bin']
    if present2 then file_delete, tmp_dir+['v2.bin', 'v2out.bin']
    if present3 then file_delete, tmp_dir+['v3.bin', 'v3out.bin']
    if present4 then file_delete, tmp_dir+['v4.bin', 'v4out.bin']
endif else file_delete, tmp_dir, /recursive

return, coordinate_out

end