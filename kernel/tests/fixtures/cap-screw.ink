; A part whose record rides on its body as `spec`, a name the file also
; defines (the catalog's form).
(import "std")
(def spec (read-json "cap-screw.json"))
(body part (shape (cylinder-z 1.5 8)) (attr spec spec))
