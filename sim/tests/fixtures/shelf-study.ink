; The shelf (shelf.ink) under 5 kg of books: a study of the sim package.
(import "../../../kernel/tests/fixtures/shelf.ink")
(study books-5kg
  (shape (min (shape-of root/board) (shape-of root/left/bracket) (shape-of root/right/bracket)))
  (elastic (E 2100MPa) (nu .38) (density 1270kg/m3) (yield 45MPa))
  (fix wall (max (- y .5) (- z 44)))
  (load books (force 0 0 -49.05N) (region (- 62 z)))
  (gravity 0 0 -9.81m/s2) (res 6))
