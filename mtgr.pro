PRO mtgr,teff,logg,feh=feh,alfe=alfe,eteff=eteff,elogg=elogg,efeh=efeh,ealfe=ealfe,extrapolate=extrapolate,silent=silent
;+
; NAME:
;       MTGR (VER. 1.2)
; PURPOSE:
;       MTGR calculates the mass of intermediate-mass main-sequence stars
;       from their effective temperature and surface gravity.
;
; SCIENTIFIC BASE:
;       This calibration is based on the Mass - Effective Temperature - 
;       Surface Gravity Relation (MTGR) developed in Kilicoglu, T. 2021, 
;       Astronomy & Astrophysics, ###, #.
;
; CALLING SEQUENCE:
;       MTGR, teff, logg, [feh=, eteff=, elogg=, efeh=, /EXTRAPOLATE, 
;           /SILENT]
;
; INPUTS:
;       teff - Effective temperature of the star [K]. 
;       logg - Logarithm of the surface gravity of the star where the gravity
;                in CGS units.
; OPTIONAL INPUT KEYWORDS:
;       feh - Metal abundance with respect to the Sun: [Fe/H] or [M/H]. 
;       alfe - alpha element abundance with respect to the Sun: [alpha/Fe].
;       eteff - Uncertainty in Teff [K]. 
;       elogg - Uncertainty in Logg (dex). 
;       efeh - Uncertainty in metal abundance (dex).
;       ealfe - Uncertainty in alpha element abundance (dex).
;       /EXTRAPOLATE - allows MTGR to extrapolate for parameters outside
;               the allowed range. Warning: The uncertainties of the results
;               obtained using this switch CAN BE MUCH GREATER than that
;               calculated. 
;       /SILENT - If SILENT is set, warning and error messages will be
;               suppressed.
;
; OUTPUTS:
;       The estimated mass appears on the command prompt.
;
;
; EXAMPLES:
;       Usage without metallicity:
;
;       IDL> MTGR,9200.,3.90
;
;
;       Usage without metallicity and with uncertainties:
;
;       IDL> MTGR,9200.,3.90,eteff=200.,elogg=0.15
;
;
;       Usage with metallicity
;
;       IDL> MTGR,9200.,3.90,feh=0.12
;
;
;       Usage with metallicity and uncertainties
;
;       IDL> MTGR,9200.,3.90,feh=0.12,eteff=200.,elogg=0.15,efeh=0.08
;
;
;       Usage with metallicity, alpha abundance and uncertainties
;
;       IDL> MTGR,9200.,3.90,feh=0.12,alfe=0.2,eteff=200.,elogg=0.15,efeh=0.08,ealfe=0.1
;
;
; RESTRICTIONS:
;       The Teff and Logg validity range of the MTGR slightly changes with
;       metallicity. The approximate limits are as follows:
;       6400 K < Teff < 20000 K,
;       Logg > 3.44, and
;       -0.4 < [Fe/H] < 0.7
;
; FUNCTIONS CALLED:
;       MTGRF()
;
; REVISION HISTORY:
; VER 1.0:       Written                  T. Kilicoglu        June, 2021
; VER 1.1: alpha abundances added         T. Kilicoglu      August, 2021  
; VER 1.2: corr. for the allowed ranges   T. Kilicoglu      August, 2021  
;
; CONTACT:
;       tkilicoglu@ankara.edu.tr
;-
IF N_ELEMENTS(teff) EQ 0 || N_ELEMENTS(logg) EQ 0 THEN BEGIN
   IF N_ELEMENTS(silent) NE 1 THEN PRINT,'ERROR: Please supply Teff and Logg.'
ENDIF ELSE BEGIN
   IF (N_ELEMENTS(feh) GE 1 && (feh GE -1.0 && feh LE 0.7)) THEN BEGIN
      feh_list=[-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7]
      llim_teff_list=[7950,7870,7670,7610,7430,7090,7040.,6840.,6770.,6550.,6420.,6400.,6390.,6380.,6360.,6355.,6350.,6340.]
      hlim_teff_list=[26000.,25500.,25500.,25500.,24500.,24000.,22500.,21700,21200.,20400.,20200.,20000.,19000.,18500.,18300.,17800.,17300.,17000.]
      llim_logg_list=[3.76,3.73,3.70,3.67,3.64,3.61,3.56,3.53,3.50,3.47,3.44,3.44,3.41,3.39,3.37,3.34,3.34,3.33]
      hlim_logg_list=[4.622,4.603,4.581,4.558,4.535,4.509,4.475,4.447,4.428,4.399,4.397,4.387,4.365,4.346,4.341,4.325,4.305,4.289]
      IF feh NE 0.7 THEN i=where(feh_list EQ FLOOR(feh*10.)/10.) ELSE i=N_ELEMENTS(feh_list)-2
      llim_teff=llim_teff_list[i]+(feh-feh_list[i])*(llim_teff_list[i+1]-llim_teff_list[i])/(feh_list[i+1]-feh_list[i])
      hlim_teff=hlim_teff_list[i]+(feh-feh_list[i])*(hlim_teff_list[i+1]-hlim_teff_list[i])/(feh_list[i+1]-feh_list[i])
      llim_logg=llim_logg_list[i]+(feh-feh_list[i])*(llim_logg_list[i+1]-llim_logg_list[i])/(feh_list[i+1]-feh_list[i])
      hlim_logg=hlim_logg_list[i]+(feh-feh_list[i])*(hlim_logg_list[i+1]-hlim_logg_list[i])/(feh_list[i+1]-feh_list[i])
   ENDIF ELSE BEGIN
      llim_teff=6400.
      hlim_teff=20000.
      llim_logg=3.44
      hlim_logg=4.387
   ENDELSE
   IF (logg LT llim_logg-0.113 || logg GT hlim_logg+0.113) || (teff LT llim_teff*0.96875 || teff GT hlim_teff*1.25) || (N_ELEMENTS(feh) EQ 1 && (feh LT -1.0 || feh GT 0.7)) || (N_ELEMENTS(alfe) EQ 1 && (alfe LT -0.2 || alfe GT 0.6)) THEN BEGIN
      IF N_ELEMENTS(silent) NE 1 THEN BEGIN 
         PRINT,'ERROR: Teff, Logg, [M/H], or [alpha/Fe] is out of the allowed range.'
         PRINT,-1.0,0.7,FORMAT="('Allowed metallicity range:           : ',F5.2, ' <    [M/H]   < ',F4.2)"
         PRINT,-0.2,0.6,FORMAT="('Allowed alpha-element abundance range: ',F5.2, ' < [alpha/Fe] < ',F4.2)"
         IF (N_ELEMENTS(feh) EQ 0) || (N_ELEMENTS(feh) EQ 1 && (feh LT -1.0 || feh GT 0.7)) THEN BEGIN 
            PRINT,'Maximum allowed values with extrapolation for the solar initial composition:'
         ENDIF ELSE BEGIN
            PRINT,feh,FORMAT="('Maximum allowed values with extrapolation for [M/H] = ',F5.2,':')"
         ENDELSE
         PRINT,CEIL(llim_teff*0.96875),FLOOR(hlim_teff*1.25),FORMAT="(I5,' <    Teff    < ',I5)"
         PRINT,CEIL((llim_logg-0.113)*100.)/100.,FLOOR((hlim_logg+0.113)*100.)/100.,FORMAT="(F5.2, ' <    logg    < ',F4.2)"
      ENDIF
   ENDIF ELSE BEGIN
      IF ((logg LT llim_logg) || (teff LT llim_teff || teff GT hlim_teff) || (N_ELEMENTS(alfe) EQ 1 && (alfe LT 0.0 || alfe GT 0.4))) && N_ELEMENTS(extrapolate) EQ 0 THEN BEGIN
         IF N_ELEMENTS(silent) NE 1 THEN PRINT,'ERROR: Teff, Logg, or [alpha/Fe] is slightly out of the allowed range. Use /EXTRAPOLATE for extrapolation.'
      ENDIF ELSE BEGIN
         IF N_ELEMENTS(eteff) EQ 0 THEN eteff=0.
         IF N_ELEMENTS(elogg) EQ 0 THEN elogg=0.
         IF N_ELEMENTS(efeh) EQ 0 THEN efeh=0.
         IF N_ELEMENTS(ealfe) EQ 0 THEN ealfe=0.
         IF (eteff LT 0. || eteff GT 3000.) || (elogg LT 0. || elogg GT 1.) || (efeh LT 0. || efeh GT 1.) || (ealfe LT 0. || ealfe GT 1.) THEN BEGIN
            IF N_ELEMENTS(silent) NE 1 THEN PRINT,'ERROR: Please check the uncertainties.'
         ENDIF ELSE BEGIN
            IF N_ELEMENTS(feh) EQ 1 THEN mass=mtgrf(teff,logg,feh) ELSE mass=mtgrf(teff,logg)
            IF eteff GT 0. THEN BEGIN
               IF N_ELEMENTS(feh) EQ 1 THEN BEGIN
                  mte_low=mtgrf(teff-eteff,logg,feh)
                  mte_hgh=mtgrf(teff+eteff,logg,feh)
               ENDIF ELSE BEGIN
                  mte_low=mtgrf(teff-eteff,logg)
                  mte_hgh=mtgrf(teff+eteff,logg)
               ENDELSE
               dmte=(mte_low-mte_hgh)/2.
            ENDIF ELSE BEGIN
               dmte=0.
            ENDELSE
            IF elogg GT 0. THEN BEGIN
               IF N_ELEMENTS(feh) EQ 1 THEN BEGIN
                  mlogg_low=mtgrf(teff,logg-elogg,feh)
                  mlogg_hgh=mtgrf(teff,logg+elogg,feh)
               ENDIF ELSE BEGIN
                  mlogg_low=mtgrf(teff,logg-elogg)
                  mlogg_hgh=mtgrf(teff,logg+elogg)
               ENDELSE
               dmlogg=(mlogg_low-mlogg_hgh)/2.
            ENDIF ELSE BEGIN
               dmlogg=0.
            ENDELSE
            IF N_ELEMENTS(feh) EQ 1 && efeh GT 0. THEN BEGIN
               mfeh_low=mtgrf(teff,logg,feh-efeh)
               mfeh_hgh=mtgrf(teff,logg,feh+efeh)
               dmfeh=(mfeh_low-mfeh_hgh)/2.
            ENDIF ELSE BEGIN
               dmfeh=0.
            ENDELSE
            IF N_ELEMENTS(alfe) EQ 1 && ealfe GT 0. THEN BEGIN 
               dmalfe=(-0.17*ALOG10(mass)+0.255)*mass*ealfe
            ENDIF ELSE BEGIN
               dmalfe=0.
            ENDELSE
            IF N_ELEMENTS(alfe) EQ 1 THEN BEGIN
               mass=mass*(1.+alfe*(-0.17*ALOG10(mass)+0.255))
            ENDIF
            dm=SQRT(dmte^2.+dmlogg^2.+dmfeh^2.+dmalfe^2.)
            IF dmte+dmlogg+dmfeh+dmalfe NE 0. THEN PRINT,mass,dm,FORMAT='("M = ",F5.2," ± ",F5.2," M☉")' ELSE PRINT,mass,FORMAT='("M = ",F5.2," M☉")'
         ENDELSE
      ENDELSE
   ENDELSE
ENDELSE
END
