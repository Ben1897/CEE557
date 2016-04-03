(TeX-add-style-hook
 "setup"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt" "letterpaper" "titlepage")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("xcolor" "svgnames") ("geometry" "margin=1in")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "latexpython"
    "article"
    "art11"
    "setspace"
    "amsmath"
    "amsfonts"
    "amssymb"
    "hyperref"
    "listings"
    "cancel"
    "natbib"
    "verbatim"
    "lastpage"
    "xcolor"
    "xspace"
    "geometry"
    "graphicx"
    "epstopdf"
    "url"
    "subfig"
    "fancyhdr"
    ""
    "upgreek"
    "rotfloat"
    "pdfpages"
    "csvsimple"
    "tikz")
   (TeX-add-symbols
    '("figref" 1)
    '("eref" 1)
    '("unit" 1)
    '("val" 2)
    '("pdd" 2)
    '("pd" 2)
    '("background" 5)
    "reportTitle"
    "reportSubTitle"
    "reportDueDate"
    "reportClass"
    "reportClassShort"
    "reportClassTime"
    "reportClassInstructor"
    "reportAuthorName"
    "reportGroupName"
    "reportGroupMembers"
    "reportAuthorNameShort"
    "reportTitleShort"
    "reportGroupNameShort"
    "VR"
    "vr")))

