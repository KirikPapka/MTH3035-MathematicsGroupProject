
; Program to calculate level sets
; used by Weng and Taylor
; model top (equivalent to zztop)
;z_t=6000.0
z_t=20000.0
; no. of grid levels (equivalent to kkp)
;N_g=141
N_g=200
z_0=0.1
b_0=67.5
;a_grid=100.0
a_grid=500.0
; b_grid=4.0
b_grid=2.0
; Uniformly distributed levels set
big_z=findgen(N_g)/float(N_g-1)
;fact=alog((z_t+z_0)/z_0) + z_t/b_0
fact= alog((z_t+b_grid)/b_grid) + z_t/a_grid

big_z=big_z*fact




; a very fine grid for performing the calculations.
n_fine=30000

z_fine=findgen(n_fine)*z_t/float(n_fine-1)


dz=fltarr(N_g)
z=fltarr(N_g)

for l=0,N_g-1 do begin
  min_diff=1e6
  for k=0, n_fine-1 do begin
;    diff=abs( big_z(l) - alog((z_fine(k)+z_0)/z_0) - z_fine(k)/b_0 )
    diff=abs( big_z(l) - alog((z_fine(k)+b_grid)/b_grid) - z_fine(k)/a_grid) 
    if (diff lt min_diff) then begin
      min_diff=diff
      z(l)=z_fine(k)
    endif

  endfor ;k loop over n_fine  
  
endfor ;l loop over N_g

for l=0,N_g-2 do begin

dz(l)=z(l+1)-z(l)
endfor

end
