(TeX-add-style-hook
 "assignment1"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "paper=letter" "fontsize=11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("babel" "english") ("placeins" "section")))
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "fontenc"
    "fourier"
    "babel"
    "amsmath"
    "amsfonts"
    "amsthm"
    "bm"
    "graphicx"
    "placeins"
    "sectsty"
    "fancyhdr")
   (TeX-add-symbols
    '("horrule" 1))
   (LaTeX-add-labels
    "fig:median"
    "fig:median_fakeplanet"
    "fig:median_curve"
    "fig:simple_adi"
    "fig:adi_fakeplanet"
    "fig:adi_median_curve"
    "fig:adi_with_error"
    "fig:loci"
    "tab:ratio")))

