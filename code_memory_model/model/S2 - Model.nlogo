; ANSI Shadow FONT http://patorjk.com/software/taag/#p=testall&f=Dancing%20Font&t=SETUP


globals
[
  ; INTERFACE PARAMETERS
    ; TRAIT CONFIGURATION
      ; ITV: continuous, one-sided width of the intraspecfic variation in behaviour types
      ; species-1-mean: continuous, mean alpha value of species 1 (= mean alpha value if n-species = 1)
      ; species-2-mean: continuous, mean alpha value of species 2

    ; COMMUNITY CONFIGURATION
      ; n-species: integer, number of species
      ; n-inds: integer, total number of individuals (NOT per species)


    ; LANDSCAPE
      ; set-this-seed: integer, seed for landscape generation - if 0 no seed is set
      ; patchiness: continouus, heterogeneity of resource distribution
      ; resource-cover: continouus, amount of resources

    ; MISC
      ; subfolder: string, where to save the output (relative to the working directory)
      ; max-ticks: integer, number of simulated time steps
      ; max-output: boolean, continouus output or only after simulation
      ; population-dynamics; boolean, enable population dynamics


  current-resources  ; counter for the number of resources that have been added to the landscape (used during initialization)

  ; parameter of the two-component memory algorithm adapted from van Moorter et al. (2009)
  rate-mem-ref       ; reference memory decay
  rate-mem-work      ; working memory decay

  ; parameters related to population dynamics and energetic costs (only needed if population dynamics are set to true)
  breeding-cost      ; increase in breeding costs per 100 time steps
  breeding-duration  ; number of time steps for the whole breeding process
  litter-size        ; number of offspring per litter
  maintenance-cost   ; maintenance costs imposed per time step


  ; output file
  file-name          ; path to save the shapefile

  ; lists with data of each animal. only used for coviability analysis (population-dynamics = true)
  list-speciesID     ; list to track species ID
  list-who           ; list to track animal ID
  list-ticks         ; list to track time step
  list-alpha         ; list to track alpha values
]



; PATCHES ARE SPATIALLY FIXED AND MAY PROVIDE RESOURCES
patches-own
[
  resource          ; boolean resource level (0 or 1)
  fertile           ; boolean for the capability of having the resource level 1
]





; ANIMALS ARE AGENTS THAT MOVE IN THE LANDSCAPE

breed[animals animal]
animals-own
[

  species-ID        ; identifier for different species

  ; elements of the memory
  mem-time          ; list with ticks since visit (only resources are memorized)
  mem-patch         ; list with memorized patches
  mem-resources     ; numeric as the number of gathered resources

  ; elements of the behavioural reaction norms
  alpha             ; behavioural trait / behaviour type that combines average behavioural expression and responsiveness
  PoD               ; persistence of direction (behavioural phenotype)
  mem-feedback      ; mem-feedback is based on the memory of patches where resources were found.
  species-mean

  ; elements of the get-memory-heading and get-PoD submodels
  mem-head          ; current heading
  memdx             ; x-direction of a memorized patch
  memdy             ; y-direction of a memorized patch
  head-x            ; list of x-directions to memorized patches weighted by the perceived utility
  head-y            ; list of y-directions to memorized patches weighted by the perceived utility
  head-x-sum        ; sum of all weighted x-directions to memorized patches
  head-y-sum        ; sum of all weighted y-directions to memorized patches

  ;population dynamics
  breeds            ; boolean, determines if the animal is currently breeding
  breeding-stage    ; integer, determines the time since the animal has started breeding
]






; ███████╗███████╗████████╗██╗   ██╗██████╗
; ██╔════╝██╔════╝╚══██╔══╝██║   ██║██╔══██╗
; ███████╗█████╗     ██║   ██║   ██║██████╔╝
; ╚════██║██╔══╝     ██║   ██║   ██║██╔═══╝
; ███████║███████╗   ██║   ╚██████╔╝██║
; ╚══════╝╚══════╝   ╚═╝    ╚═════╝ ╚═╝

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;; SETUP PROCEDURE ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup

  ca


  ;track change in trait distribution if population dynamics are enabled
  if population-dynamics [
    set list-speciesID ["Species"]
    set list-who ["who"]
    set list-ticks ["ticks"]
    set list-alpha ["alpha"]
  ]


  ; define a working directory
  set-current-directory working-dir


  ; define parameters of the two-component memory algorithm
  set rate-mem-ref 0.99
  set rate-mem-work 0.999

  ; define parameters related to population dynamics and energetic costs (only neccessary if population dynamics are true

  set breeding-cost 0.0075
  set breeding-duration 3600
  set litter-size 5
  set maintenance-cost 0.18


  ; create a unique file name and create an initial .csv-file
  let file-id random 999999999
  set file-name (word subfolder (remove ":" remove "" remove " " remove "." date-and-time) "_" "p-" patchiness "_rc-" resource-cover file-id)
  if max-output [let world-name (word subfolder (remove ":" remove "" remove " " remove "." date-and-time) "_" "p-" patchiness "_rc-" resource-cover file-id "_world.png")]

  file-open word file-name ".csv"
  file-print (word "XCOR" "," "YCOR" "," "ALPHA" "," "RESOURCESG" "," "SPECIESID" "," "PATCHINESS" "," "RESOURCECO" "," "TRAITMEAN" "," "POPSIZE" "," "DATETIME" "," "WHO" "," "ITV" "," "SPECMEAN" "," "SEED")
  file-close


  ; generate landscape and animal populations via submodels
  set current-resources 0 ; initial value for resources added to the landscape
  create-landscape

  ; export the landscape as a .png if neccessary
  if max-output [
    let world-name (word subfolder (remove ":" remove "" remove " " remove "." date-and-time) "_" "p-" patchiness "_rc-" resource-cover file-id "_world.png")
     export-view world-name]


  spawn-animals




  ; update plots
  set-current-plot "Trait Distribution"
  histogram [alpha] of turtles


  reset-ticks

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;; SETUP SUBMODELS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to create-landscape

  if set-this-seed != 0 [random-seed set-this-seed]

  ; set all resources to 0
  ask patches[
    set resource 0
    set pcolor white
  ]

  ; calculate the required number of resources and clumps.
  let max-resources (count patches) * resource-cover / 100
  let n-clumps floor ((100 - patchiness) / 100 * max-resources)

  if n-clumps <= 0 [set n-clumps 1]

  ; add one resource to each clump
  ask n-of n-clumps patches [fertilize]

  set current-resources n-clumps


  ; add resources the neighbors of existing resource patches
  while [current-resources < (max-resources)]
  [
    ask one-of patches with [resource = 1][
      if any? neighbors with [resource = 0] [
        set current-resources current-resources + 1
        ask (one-of neighbors with [resource = 0]) [fertilize]
      ]
    ]
  ]

  random-seed new-seed

  show (word sum [resource] of patches "of" max-resources "resources were created")
end


; subsubmodel to add resources to the landscape
to fertilize
  set resource 1
  set fertile 1
  set pcolor green
end



; submodel to spawn animals along with their species-ID and behaviour types (alpha)
to spawn-animals

  let species-added 0
  let inds-per-species floor n-inds / n-species

  ; set individuals per species with alpha sampled from a species-specific distribution
  while [species-added < n-species]
  [
    set species-added species-added + 1

    create-animals inds-per-species
    [
      ; sample individual alpha from a uniform distribution around the species-specific mean with ITV as the one-sided width of the distribution
      ifelse species-added = 1 [set species-mean species-1-mean][set species-mean species-2-mean]

      let valid-alpha false
      let propose-alpha -999
      while [not valid-alpha]
      [
        set propose-alpha species-mean + random-float ITV * 2 - ITV
        if propose-alpha <= 1 and propose-alpha >= 0 [set valid-alpha true]
      ]

      set alpha propose-alpha

      set species-ID species-added
    ]
  ]



  ; create empty memory and add animals add random locations
  ask animals
  [
    set mem-time (list)
    set mem-patch (list)
    set mem-resources 0
    if population-dynamics [set mem-resources 20]
    set xcor random-xcor set ycor random-ycor
    set breeds false
    set breeding-stage 0
  ]

end





; ██████╗  ██████╗
; ██╔════╝ ██╔═══██╗
; ██║  ███╗██║   ██║
; ██║   ██║██║   ██║
; ╚██████╔╝╚██████╔╝
;  ╚═════╝  ╚═════


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; GO PROCEDURE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to go

  ; measure duration of simulation run

  if (ticks < 1)
  [
    if population-dynamics [save-output]
    reset-timer
  ]


  ; animal submodels
  ask animals[
    ; as soon as an individual has memorized 2 patches, start using memory
    if length mem-time > 1
    [
      get-memory-heading
      get-PoD
    ]

    set-new-heading-and-move
    update-memory

    if population-dynamics [
      do-energetics
      if not breeds [attempt-breeding]
      create-offspring


    ]
  ]

  if (length remove-duplicates [species-ID] of turtles = 1 and population-dynamics) or count turtles = 0
    [
      save-output
      stop
    ]

  ; patches submodels
  grow-resources


  ; end simulation
  if ticks = max-ticks [
    save-output
    show timer
    stop
  ]

  if ticks > 1000 [
    ask turtles [pen-down]
  ]

  if population-dynamics and ticks mod 1000 = 0
  [
    set-current-plot "Trait Distribution"
    histogram [alpha] of turtles
  ]

  ; track movement via max-output after a spin-up-time of 1000 time steps
  if max-output and ticks > 1000
  [
    save-output
  ]


  if (population-dynamics) and (ticks mod 50 = 0)
  [
    ask animals [
        set list-speciesID lput species-ID list-speciesID
        set list-who lput who list-who
        set list-ticks lput ticks list-ticks
        set list-alpha lput alpha list-alpha
    ]
  ]

  tick

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;; GO SUBMODELS  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; MEMORY RELATED ;;;

to get-memory-heading


  ; intialize lists of utility-weighted memory directions
  set head-x (list)
  set head-y (list)

  ; calculate perceived utility - the formulation differs slightly from van Moorter et al. (2009), I use a decay function instead of an iterative approach
  (foreach but-last mem-patch but-last mem-time
    [ [?1 ?2] ->
      ; calculate perceived utility for memorized resource locations
      let dist distance ?1 ;distance to patch
      let X [pxcor] of ?1 ;x-coordinate of the patch
      let Y [pycor] of ?1 ;y-coordinate of the patch

      ; MODEL EQUATION I
      let util (1 - (rate-mem-work ^ (?2 + 1))) * (rate-mem-ref ^ (?2 + 1)) / (dist ^ 2) ;utility based on decay functions and distance

      ; on periodical landscapes the setting the heading via face xy does not account for the shortest path
      ; so the memory vector is always biased towards the center as it never accounts for the torus
      ; to circumvent this, we need to determine the shortest path on the toroid landscape manually
      let Xpos (list abs (X - max-pxcor - xcor) abs (X - xcor) abs (X + max-pxcor - xcor)) ; possible x-distances of the patch on periodical landscapes (left, mid, right)
      let Ypos (list abs (Y - max-pycor - ycor) abs (Y - ycor) abs (Y + max-pycor - ycor)) ; possible y-distances of the patch on periodical landscapes (bottom, mid, top)

      ; check which x and y  direction is the shortest on the periodical landscape
      let newx position min Xpos Xpos
      let newy position min Ypos Ypos

      if newx = 0 [set memdx X - max-pxcor - xcor]
      if newx = 1 [set memdx X - xcor]
      if newx = 2 [set memdx X + max-pxcor - xcor]

      if newy = 0 [set memdy Y - max-pycor - ycor]
      if newy = 1 [set memdy Y - ycor]
      if newy = 2 [set memdy Y + max-pycor - ycor]

      let magnitude sqrt (memdx ^ 2 + memdy ^ 2)

      ; MODEL EQUATION II for x and y component: weight the vector to the memorized location by its length
      set memdx memdx / magnitude
      set memdy memdy / magnitude

      ; MODEL EQUATION III for x and y component: calculate the xy-elements of the single memory vector
      set head-x lput (memdx * util) head-x
      set head-y lput (memdy * util) head-y
    ]
  )


  ; MODEL EQUATION IV for the x and y component: sum of all X and Y directions of all memories
  set head-y-sum sum head-y
  set head-x-sum sum head-x

  ; calculate  the heading of the combined attraction vector from memory
  set mem-head towardsxy (xcor + head-x-sum) (ycor + head-y-sum)

end

to get-PoD
  ; MODEL EQUATION V: asses the environmental gradient (memory-feedback)
  set mem-feedback  (sqrt (head-y-sum ^ 2 + head-x-sum ^ 2)) ^ .25

  ; MODEL EQUATION VI: determine PoD
  set PoD 1 - mem-feedback * alpha
  if PoD < 0 [set PoD 0]

end

to set-new-heading-and-move
  ; MODEL EQUATION VII: change the heading based on memory weighted by the persistence of direction
  let mean-head heading + ((1 - PoD) * (subtract-headings mem-head heading))
  set heading rvm mean-head 10
  fd 1
end

; add new patches to the memory, update the age of the memory, exploit resources
to update-memory
  ; check if the patch provides resources
  if [resource] of patch-here = 1
  [
    ; if yes, check if the patch is already part of the memory
    ifelse member? patch-here mem-patch
    [
      ; if yes update time
      set mem-time replace-item (position patch-here mem-patch) mem-time 0
    ]
    [
      ; if no, add patch to memory and set time since last visit to 0
      set mem-patch lput patch-here mem-patch
      set mem-time lput 0 mem-time
    ]

    ; exploit resources and increase the number of foraged resources by 1
    if ticks > 1000 or population-dynamics [set mem-resources mem-resources + 1]

    ask patch-here [
      set resource 0
      set pcolor orange
    ]
  ]

  if length mem-time > 0 [
    set mem-time map [ ?1 ->  ?1 + 1 ] mem-time
  ]
end

;;; LANDSCAPE RELATED ;;;

; regrow resources with a probability of 1 %
to grow-resources
  ask patches with [resource = 0 and fertile = 1][
    if random-float 1 > 0.99 [
      set resource 1
      set pcolor green
    ]
  ]
end

; submodel to save state variables and parameters
to save-output
  file-open word file-name ".csv"

  ask animals
  [
  file-print (word (precision xcor 2)  ","  (precision ycor 2) "," alpha "," mem-resources "," species-ID "," patchiness "," resource-cover "," species-1-mean "," n-inds "," ticks "," who "," ITV "," species-mean "," set-this-seed)
  ]

  file-close


  if population-dynamics
  [
    file-open word file-name "_alpha.csv"

    file-print  (word list-speciesID "," list-who ["who"] "," list-ticks ["ticks"] "," list-alpha ["alpha"])

    file-close


  ]

end


; Related to population dynamcis


; attempt breeding with a probability that depends on the number of available resources
to attempt-breeding
  if random-float 10000 < mem-resources
  [
    set breeds true
    set breeding-stage 0
  ]
end

; impose energetic costs due to maintenance and breeding
to do-energetics
      set mem-resources mem-resources - maintenance-cost

    if breeds
    [
      set breeding-stage breeding-stage + 1

      set mem-resources mem-resources - (breeding-cost * (breeding-stage / 100))
    ]

    if mem-resources < 0[die]
end


to create-offspring
  if breeding-stage > breeding-duration[

   set breeding-stage 0
   set breeds false
   let parental-alpha alpha

   hatch litter-size[
     set mem-time (list)
     set mem-patch (list)
     set mem-resources 10
     set breeds false
     set alpha species-mean + (random-float (ITV * 2) - ITV) ;random sampling from species-specific distribution
     if alpha > (species-mean + ITV) [set alpha species-mean + ITV]
     if alpha < (species-mean - ITV) [set alpha species-mean - ITV]
     if alpha > 1 [set alpha 1]
     if alpha < 0 [set alpha 0]
     set species-ID species-ID
    ]
   ]
end


; ██████╗ ███████╗██████╗  ██████╗ ██████╗ ████████╗███████╗██████╗ ███████╗
; ██╔══██╗██╔════╝██╔══██╗██╔═══██╗██╔══██╗╚══██╔══╝██╔════╝██╔══██╗██╔════╝
; ██████╔╝█████╗  ██████╔╝██║   ██║██████╔╝   ██║   █████╗  ██████╔╝███████╗
; ██╔══██╗██╔══╝  ██╔═══╝ ██║   ██║██╔══██╗   ██║   ██╔══╝  ██╔══██╗╚════██║
; ██║  ██║███████╗██║     ╚██████╔╝██║  ██║   ██║   ███████╗██║  ██║███████║
; ╚═╝  ╚═╝╚══════╝╚═╝      ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚══════╝



; van-Mises distribution
; this distribution can be used to implement correlated random walks (warum einfach, wenn es auch kompliziert geht):
; reimplementation of the rvm-function in the R-packages CircStats

to-report rvm [vm-mean k]
  set vm-mean vm-mean * pi / 180
  let a 1 + (1 + 4 * (k ^ 2)) ^ 0.5
  let b (a - (2 * a) ^ 0.5)/(2 * k)
  let r (1 + b ^ 2)/(2 * b)
  let calculated false

  while [not calculated]
  [
    let U1 random-float 1
    let z cos (180 * U1)
    let f (1 + r * z)/(r + z)
    let c k * (r - f)
    let U2 random-float 1

    if c * (2 - c) - U2 > 0 or (log (c / U2) e) + 1 - c > 0 [
      ifelse random-float 1 - 0.5 > 0 [
        let U3 -1
        let vm-rnd U3 * ((acos f) * pi / 180) + vm-mean
        set vm-rnd vm-rnd mod (2 * pi)
        set calculated true
        report (vm-rnd * 180 / pi)
      ]
      [
        let U3 1
        let vm-rnd U3 * ((acos f) * pi / 180) + vm-mean
        set vm-rnd vm-rnd mod (2 * pi)
        set calculated true
        report (vm-rnd * 180 / pi)
      ]
    ]

  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
501
80
884
464
-1
-1
1.5
1
10
1
1
1
0
1
1
1
0
249
0
249
1
1
1
ticks
30.0

SLIDER
17
378
269
411
resource-cover
resource-cover
1
100
60.0
1
1
%
HORIZONTAL

SLIDER
17
416
270
449
patchiness
patchiness
0
100
10.0
1
1
%
HORIZONTAL

BUTTON
21
28
198
93
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
249
27
455
90
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
18
532
271
565
n-inds
n-inds
1
100
750.0
1
1
NIL
HORIZONTAL

INPUTBOX
19
652
269
712
subfolder
NIL
1
0
String

SLIDER
18
487
271
520
n-species
n-species
1
2
2.0
1
1
NIL
HORIZONTAL

TEXTBOX
287
487
437
557
COMMUNITY CONFIGURATION\n\nthe number of individuals must be an integer multiple of the number of species
11
0.0
1

TEXTBOX
281
314
462
470
LANDSCAPE CONFIGURATION\n\nset-this-seed: landscape generation seed, if 0 no seed is set\n\nresource-cover: proportion of pixels with resources\n\npatchiness: 100 % minus the proportion of (resource-cover minus
11
0.0
0

TEXTBOX
294
600
444
698
OUTPUT CONFIGURATION\n\nsubfolder is the directory where the output is saved. The current directory of the model is the \"working directory\".
11
0.0
1

SLIDER
19
720
270
753
max-ticks
max-ticks
0
5000
2000.0
100
1
NIL
HORIZONTAL

PLOT
939
78
1139
228
Trait Distribution
NIL
NIL
-0.2
1.2
0.0
1.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" ""

SWITCH
19
604
269
637
max-output
max-output
0
1
-1000

SLIDER
16
337
269
370
set-this-seed
set-this-seed
0
100
0.0
1
1
NIL
HORIZONTAL

PLOT
932
268
1132
418
mean PoD
NIL
NIL
0.0
10.0
1.0
1.0
true
false
"" ""
PENS
"pen-1" 1.0 0 -7500403 true "" "plot mean [PoD] of animals"

SLIDER
21
165
273
198
ITV
ITV
0
.5
0.09
0.01
1
NIL
HORIZONTAL

INPUTBOX
18
761
261
821
working-dir
C:\\Users\\milles\\Documents\\
1
0
String

SWITCH
21
836
190
869
population-dynamics
population-dynamics
1
1
-1000

SLIDER
22
224
133
257
species-1-mean
species-1-mean
0
1
1.0
0.01
1
NIL
HORIZONTAL

SLIDER
141
225
264
258
species-2-mean
species-2-mean
0
1
0.25
.01
1
NIL
HORIZONTAL

TEXTBOX
282
140
432
280
COMMUNITY CONFIGURATION\n\nITV: degree of intraspecific variation in alpha\n\n\nspecies-1-mean / species-2-mean: Mean alpha level of species-specific distribution
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="fitness_trait_relationship_750inds" repetitions="40" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="species-1-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_fitness_trait_relationship_750inds/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
      <value value="70"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="750"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="fitness_trait_relationship_375inds" repetitions="40" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="species-1-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_fitness_trait_relationship_375inds/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
      <value value="70"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="375"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="trade_off_750inds" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="species-1-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_trade_off_750inds/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
      <value value="70"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="750"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="trade_off_375inds" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_trade_off_375inds/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
      <value value="70"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="species-1-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="375"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="coviability_analysis_heterogeneous" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200518/final_coviability_analysis_heterogeneous_noher/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="750"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.128"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.255" step="0.035" last="0.745"/>
    <steppedValueSet variable="species-2-mean" first="0.255" step="0.035" last="0.745"/>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="coviability_analysis_intermediate" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_coviability_analysis_intermediate/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="750"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.128"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="70"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.255" step="0.035" last="0.745"/>
    <steppedValueSet variable="species-2-mean" first="0.255" step="0.035" last="0.745"/>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="coviability_analysis_homogeneous" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200514/final_coviability_analysis_homogeneous/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-inds">
      <value value="750"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0.128"/>
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="2000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.255" step="0.035" last="0.745"/>
    <steppedValueSet variable="species-2-mean" first="0.255" step="0.035" last="0.745"/>
    <enumeratedValueSet variable="population-dynamics">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="coviability_analysis_pop_sim_homogeneous" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="n-inds">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="200000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200520/Coviability_analysis_homogeneous_noher/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0"/>
      <value value="0.128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.25" step="0.0625" last="0.75"/>
    <steppedValueSet variable="species-2-mean" first="0.25" step="0.0625" last="0.75"/>
  </experiment>
  <experiment name="coviability_analysis_pop_sim_heterogeneous" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="n-inds">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="200000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200520/Coviability_analysis_heterogeneous_noher/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0"/>
      <value value="0.128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.25" step="0.0625" last="0.75"/>
    <steppedValueSet variable="species-2-mean" first="0.25" step="0.0625" last="0.75"/>
  </experiment>
  <experiment name="coviability_analysis_pop_sim_intermediate" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="n-inds">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-ticks">
      <value value="200000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="working-dir">
      <value value="&quot;//oesa-erc/Oesa/Public/milles/proj_Memory_Model/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchiness">
      <value value="70"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cover">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="subfolder">
      <value value="&quot;simulations/200520/Coviability_analysis_intermediate_noher/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-output">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-species">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-this-seed">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ITV">
      <value value="0"/>
      <value value="0.128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population-dynamics">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="species-1-mean" first="0.25" step="0.0625" last="0.75"/>
    <steppedValueSet variable="species-2-mean" first="0.25" step="0.0625" last="0.75"/>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
