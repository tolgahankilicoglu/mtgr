FUNCTION mtgrf,te,logg,feh
;+
; NAME:
;       MTGRF (VER. 1.2); 
; PURPOSE:
;       This raw function DOES NOT INCLUDE ERROR CHECKING and was created for use 
;       only by MTGR and MTGRFUNC. Please use MTGRFUNC instead of this function!
;-
IF N_ELEMENTS(feh) EQ 0 THEN BEGIN
   a=0.769687
   b=0.409173
   c=0.611016
   d=-0.204335
   e=-0.175664
ENDIF ELSE BEGIN
   IF feh GE -0.4 AND feh LE 0.7 THEN BEGIN 
      a=0.037539*feh^2.+0.106185*feh+0.760313
      b=-0.074620*feh^3.-0.035236*feh^2.+0.139705*feh+0.397282
      c=0.122190*feh^3.-0.001445*feh^2.-0.387297*feh+0.645036
      d=0.017391*feh^3.-0.015925*feh^2.+ 0.021622*feh-0.205739
      e=-0.167246*feh-0.191036
   ENDIF ELSE BEGIN
      IF feh GE -0.5 AND feh LT -0.4 THEN BEGIN
         a=-0.213170*feh+0.641994
      ENDIF ELSE BEGIN
         a=0.080416*feh+0.788564
      ENDELSE
      b=0.051362*feh^2.+0.175915*feh+0.402455
      c=0.500361*feh^3.+1.059979*feh^2.+0.447992*feh+0.831917
      d=-0.040907*feh^2.+ 0.024522*feh-0.201473
      e=0.215200*feh-0.174482
   ENDELSE
ENDELSE
mass=10.^(a*(ALOG10(te)-4.)+b*10.^(c*(ALOG10(te)-4.))+d*(logg-4.)*10.^(e*(logg-4.)))
RETURN,mass
END
