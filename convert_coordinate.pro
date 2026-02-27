; coordinate: the array of coordinate to be tranformed, 
; its dimesions can be (3) or (3, n1) or (3, n1, n2) or (3, n1, n2, n3) 

; v1, v2, v3: the array of vectors to be tranformed, 
; their dimesions should be the same as coordinate

; mode:
; 0 or 'xyz_to_lon_lat_r'
; 1 or 'lon_lat_r_to_xyz'
; 2 or 'lon_lat_r_to_lon2_lat2_r' ; not finished
; 3 or 'lon2_lat2_r_to_lon_lat_r' ; not finished


function convert_coordinate, coordinate, v1, v2, v3, mode=mode, $
v1out=v1out, v2out=v2out, v3out=v3out
;-----------------------------------------------------
; print, mode, size(mode,/tname) eq 'STRING', mode eq 'xyz_to_lon_lat_r'
if ~keyword_set(mode) then mode=0
if size(mode,/tname) eq 'STRING' then begin
    case mode of
		'xyz_to_lon_lat_r': mode=0
        'lon_lat_r_to_xyz': mode=1
        ; 'lon_lat_r_to_lon2_lat2_r': mode= 2
        ; 'lon2_lat2_r_to_lon_lat_r': mode= 3
    endcase
    if size(mode,/tname) eq 'STRING' then message, 'Something is wrong with mode'
endif
if mode eq 0 then two_pi=!pi*2.
;-----------------------------------------------------

sz_coor=size(coordinate)
if sz_coor[1] ne 3 then message, 'Something is wrong with coordinate'

ngrid=sz_coor[n_elements(sz_coor)-1]/3
;-----------------------------------------------------

present1= N_PARAMS() ge 2
present2= N_PARAMS() ge 3
present3= N_PARAMS() ge 4

if present1 then begin
    sz_v=size(v1)
    if sz_v[0] ne sz_coor[0]     then message, 'Something is wrong with v1'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v1'
    v1out=v1
endif

if present2 then begin
    sz_v=size(v2)
    if sz_v[0] ne sz_coor[0]     then message, 'Something is wrong with v2'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v2'
    v2out=v2
endif

if present3 then begin
    sz_v=size(v3)
    if sz_v[0] ne sz_coor[0]     then message, 'Something is wrong with v3'
    dummy=where(sz_v ne sz_coor, count)
    if count ne 0 then message, 'Something is wrong with v3'
    v3out=v3
endif
;-----------------------------------------------------
coordinate_out=coordinate
coor_out=fltarr(3)

for i=0, ngrid-1 do begin
    coor_in=coordinate[i*3:i*3+2]
    if mode eq 0 then begin ; 'xyz_to_lon_lat_r'
        ; coor_out[2]=norm(coor_in) ; This approach is three times slower.
        coor_out[2]= SQRT(TOTAL(ABS(coor_in)^2))
        coor_out[1]= asin(coor_in[2]/coor_out[2])
        if coor_in[0] eq 0. and coor_in[1] eq 0. then begin 
            coor_out[0] = 0.
        endif else begin
            cos_lon= coor_in[0]/SQRT(TOTAL(ABS(coor_in[0:1])^2))
            if ~(cos_lon lt 1.) then begin
                coor_out[0] = 0.
            endif else if (cos_lon le -1.) then begin
                coor_out[0] = !pi
            endif else if (coor_in[1] ge 0.) then begin
                coor_out[0] =        acos(cos_lon)
            endif else begin
                coor_out[0] = two_pi-acos(cos_lon)
            endelse
        endelse
        if present1 then begin
            sin01=sin(coor_out[0:1])
            cos01=cos(coor_out[0:1])
        endif
    endif else if mode eq 1 then begin    ; 'lon_lat_r_to_xyz'
        sin01=sin(coor_in[0:1])
        cos01=cos(coor_in[0:1])
        coor_out=coor_in[2]*[cos01[1]*cos01[0], cos01[1]*sin01[0], sin01[1]]
    endif
    coordinate_out[i*3:i*3+2]=coor_out

    if present1 then begin
        ; stack e_lon, e_lat, e_r for 'xyz_to_lon_lat_r'
        matrix=[[         -sin01[0],           cos01[0],       0.],$
                [-sin01[1]*cos01[0], -sin01[1]*sin01[0], cos01[1]],$
                [ cos01[1]*cos01[0],  cos01[1]*sin01[0], sin01[1]]]
        if mode eq 1 then matrix=transpose(matrix)

        v1out[i*3:i*3+2]= reform(matrix ## v1[i*3:i*3+2])
    endif

    if present2 then v2out[i*3:i*3+2]= reform(matrix ## v2[i*3:i*3+2])
    if present3 then v3out[i*3:i*3+2]= reform(matrix ## v3[i*3:i*3+2])
endfor

return, coordinate_out
end