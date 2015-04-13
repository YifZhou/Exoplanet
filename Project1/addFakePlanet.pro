PRO addFakePlanet, PSFAmplitude = PSFAmplitude
  IF n_elements(PSFAmplitude) EQ 0 THEN PSFAmplitude = 2000
  PSF = readfits('median_unsat.fits')
  PSF = (PSF - min(PSF))/max(PSF) ;; normalize the PSF
  ims = readfits('center_im.fits')
  rotAngle = readfits('rotnth.fits')
  pixelScale = 0.009942 ;; pixel scale of narrow camera on KECK
  pixelSep = findgen(10)/9.0 * 3/pixelScale ;; the separation of the fake planets and the primary. In angular Spacep, sepration is set from 0 arcsec to 3 arcsec. 10 evenly distributed points are sampled.

  ;; calculate the x and y position of the peak of fake planet in
  ;; unrotated images and add the fake planet
  theta0 = [0, 60, 120, 180, 240, 300] ;; position angle of fake planet in rotated images
  FOR i = 0, (size(ims))[(size(ims))[0]] - 1 DO BEGIN
     FOR j = 0, n_elements(theta0) -1 DO BEGIN
        theta = (theta0[j] + rotAngle[i])/180. * !PI
        FOR k = 0, n_elements(pixelSep) - 1 DO BEGIN
           x = pixelSep[k] * cos(theta) + 612 ;; 612 is the center of the image
           y = pixelSep[k] * sin(theta) + 612
           ims[x-29.5:x+29.5, y-29.5:y+29.5, i] = ims[x-29.5:x+29.5, y-29.5:y+29.5, i] + PSF*PSFAmplitude
           ;; in IDL use the nearest integer when the index is
           ;; fractional, uncertaity can be introduced here
           ;; use 29.5 to match the size of PSF
        ENDFOR 
     ENDFOR
  ENDFOR
  outFn = 'fakePlanetCube_Amp=' + strn(PSFAmplitude) + '.fits'
  writefits, outFn, ims
END
