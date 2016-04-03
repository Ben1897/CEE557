(TeX-add-style-hook
 "CP2"
 (lambda ()
   (TeX-run-style-hooks
    "setup"
    "Part1"
    "Part2"
    "Appendix")
   (LaTeX-add-labels
    "eq:1"
    "eq:2a"
    "eq:2b"
    "eq:2c")))

