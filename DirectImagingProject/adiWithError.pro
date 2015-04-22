PRO adiwithError, infn
  cube = readfits(infn)
  rotnth= readfits('rotnth.fits')
  adiCube = fltarr((size(cube))[1:(size(cube))[0]])
  index = findgen((size(cube))[(size(cube))[0]])
  FOR i=0, (size(cube))[(size(cube))[0]] - 1 DO BEGIN
     cube[*, *, i] = fshift(cube[*, *, i], 5*randomn(seed), 5*randomn(seed))
  ENDFOR
  
  FOR i=0, (size(cube))[(size(cube))[0]] - 1 DO BEGIN
     adiCube[*, *, i] = rot(cube[*, *, i] - median(cube[*, *, where(index ne i)], dimension = 3), rotnth[i], 1., 1224./2, 1224./2, /piv, cub = -0.5)
     print, 'number ', i, ' image finished adi'
  ENDFOR
  median_im = median(adiCube, dimension = 3)
  writefits, 'ADIwithError_' + infn, median_im
END
