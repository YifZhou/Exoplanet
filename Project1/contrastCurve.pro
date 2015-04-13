PRO contrastCurve, image, PSF_Amp, outFN
  thetaList = [0, 60, 120, 180, 240, 300] ;; position angle
  pixelScale = 0.009942
  pixelSep = findgen(10)/9.0 * 3/pixelScale
  contrastCurveStack = fltarr(10, 6)
  meshgrid, findgen(1224), findgen(1224), xx, yy
  FOR i =0, 5 DO BEGIN ;; loop through each position angle
     x = cos(thetaList[i]/180.*!PI)*pixelSep + 1224./2.
     y = sin(thetaList[i]/180.*!PI)*pixelSep + 1224./2.
     FOR j =0, 9 DO BEGIN 
        dist = sqrt((xx - x[j])^2 + (yy - y[j])^2)
        inAnnuli = where((dist le 20) AND (dist GE 15))
        sky = sqrt(mean(image(inAnnuli)^2))
       ;; use aper to do apeture photometure
       ;; calculated sky value is noise
       contrastCurveStack[j, i] = PSF_Amp/sky ;;signal is known as the peak of PSF
    ENDFOR 
  ENDFOR
  Curve = mean(contrastCurveStack,dimension = 2)
  plot, pixelSep * pixelScale, Curve
  forprint, pixelSep * pixelScale, Curve, textout = outFN,/NoCOMMENT
END
