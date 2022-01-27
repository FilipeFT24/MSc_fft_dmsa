# 2D Code structure
## Top-level directory layout

## Grid generation
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_11_1.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_11_2.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_11_3.png" width="250" />
</p>
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_12_1.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_12_2.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_12_3.png" width="250" />
</p>

## Stencil generation
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_ST9_1.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_2.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_3.png" width="250" />
</p>
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_ST9_4.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_5.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_6.png" width="250" />
</p>
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_ST9_7.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_8.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_9.png" width="250" />
</p>
<p align="center">
  <img src="/[Figures]/Fig_1/Fig_ST9_10.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_11.png" width="250" />
  <img src="/[Figures]/Fig_1/Fig_ST9_12.png" width="250" />
</p>

## 2D quadrature
<p align="center">
  <img src="/[Figures]/Fig_0/Fig_01_1.png" width="250" />
  <img src="/[Figures]/Fig_0/Fig_01_5.png" width="250" />
  <img src="/[Figures]/Fig_0/Fig_01_9.png" width="250" />
</p>
<p align="center">
  <img src="/[Figures]/Fig_0/Fig_02_1.png" width="250" />
  <img src="/[Figures]/Fig_0/Fig_02_5.png" width="250" />
  <img src="/[Figures]/Fig_0/Fig_02_9.png" width="250" />
</p>

```
'
├── [Tools - Data]            # >>> Save and load data    (for plotting purposes...).
│
├── [Tools - Numerical]       # >>> Numerical tools (mostly) used for code speed up/efficiency.
│
├── [Tools - Post-processing] # >>> Post-processing tools (for plotting purposes...).
│
└── src  
    │
    ├── Main                  # >>> Call 'Class' functions.
    │
    ├── A                     # >>> Call 'SubClass' functions.
    │   ├── A_1               #  >> Set inputs.
    │   ├── A_2               #  >> Generate grid and set boundaries.
    │   ├── A_3               #  >> Call 'SubSubClass' functions.
    │   │   ├── A_3_1         #   > Compute grid properties.
    │   │   └── A_3_2         #   > Call 'SubSubSubClass' functions.
    │   │       ├── A_3_2_1   #   > Compute stencil elements (w/o extension(s)).
    │   │       └── A_3_2_2   #   > Compute stencil elements (w/  extension(s)).
    │   └── A_Tools           #  >> 'Class' Tools.
    │
    ├── B                     # >>> Call 'SubClass' functions.
    │   ├── B_1               #  >> Call 'SubSubClass' functions.
    │   │   ├── B_1_1         #   > Compute analytic solution.
    │   │   └── B_1_2         #   > Set 1D/2D quadrature rules.
    │   │                     #   > Compute source term.
    │   ├── B_2               #  >> Call 'SubSubClass' functions.
    │   │   ├── B_2_1         #   > Compute face polynomial coefficients (recursively).
    │   │   │                 #   > Set weighting function.
    │   │   └── B_2_2         #   > Assemble matrices Df, Dwf, Pf, Tf, A and B.
    │   │                     #   > Compute source term.
    │   │                     #   > Compute error norms.
    │   └── B_3               #  >> 
    │
    └── C                     # >>> Call 'Fig_X' functions.
        ├── Fig_0             #  >> Fig_0.
        ├── Fig_1             #  >> Fig_1.
        ├── Fig_2             #  >> Fig_2.
        ├── Fig_Tools         #  >> Fig_Tools.
        └── ...               #  >> ...    
```
```
'
inp 
├── msh                       # >>> Field: Grid.
│   ├── lim                   #  >> Field: Limits.
│   │   ├── Xv_i              #   > LHS (x-direction).
│   │   ├── Xv_f              #   > RHS (x-direction).
│   │   ├── Yv_i              #   > LHS (y-direction).
│   │   └── Yv_f              #   > RHS (y-direction).
│   ├── h                     #  >> ...case of a     uniform grid: "Approximated" hydraulic diameter used to generate the grid. 
│   ├── Nv                    #  >> ...case of a non-uniform grid: Number of vertices used in the x and y-directions used to generate the grid: [Nv(:)] = [Nv(X),Nv(Y)].  
│   ├── pt                    #  >> Face polygon type.
│   │   ├── (1)               #   - Option (1): Triangles  ('v').
│   │   └── (2)               #   - Option (2): Squares    ('s').
│   ├── eg                    #  >> Mesh example type #1   (Check SubClass A_2).
│   ├── dm                    #  >> Mesh example type #2   (Check SubClass A_2).
│   └── su                    #  >> Clustering parameters  (...case of a non-uniform grid).
│       ├── Nf_X              #   > Clustering location    (x-direction).
│       ├── Nf_Y              #   > Clustering location    (y-direction).
│       ├── Ks_X              #   > Clustering factor      (x-direction).
│       └── Ks_Y              #   > Clustering factor      (y-direction).  
│
├── pr                        # >>> Field: Problem reconstruction.
│   ├── vx                    #  >> Convection coefficient (x-direction).
│   ├── vy                    #  >> Convection coefficient (y-direction).
│   ├── gx                    #  >> Diffusion  coefficient (x-direction).
│   └── gy                    #  >> Diffusion  coefficient (y-direction).                
│
└── fr                        # >>> Field: Flux reconstruction.
    ├── st                    #  >> Simulation type.
    │   ├── (1)               #   - Option (1): Implicit.
    │   ├── (2)               #   - Option (2): Explicit.
    │   └── (3)               #   - Option (3): Deferred-correction (DC).
    ├── wf                    #  >> Weighting function.
    │   ├── (1)               #   - Option (1): Unweighted.
    │   └── (2)               #   - Option (2): Weighted.
    ├── np                    #  >> Face  polynomial order.
    ├── ng                    #  >> 1D/2D quadrature order.
    ├── nt                    #  >> Face neighbouring type.
    │   ├── (1)               #   - Option (1): Vertex neighbour.          ('false')
    │   └── (2)               #   - Option (2): Face   neighbour.          ('true ')
    └── et                    #  >> Extension type.
        ├── (1)               #   - Option (1): Do nothing...              ('false')
        └── (2)               #   - Option (2): Perform "extra" extension. ('true ')               
```
```
'
msh
│ 
├── d                # >>> Field: Domain.
│   └── h_ref        #  >> Hydraulic diameter. 
│
├── c                # >>> Field: Cell.
│   ├── NC           #  >> Number of cells.
│   ├── xy_v         #  >> Cell 'i' vertices: [Xv(:),Yv(:)] = [xy_v{i}(:,1),xy_v{i}(:,2)].
│   ├── mean         #  >> Cell 'i' centroid: [Xv,Yv] = [mean(1,i),mean(2,i)].
│   ├── h            #  >> Cell 'i' hydraulic diameter.
│   ├── vol          #  >> Cell 'i' volume.
│   ├── c            #  >> Cell 'i' neighbouring cell indices. 
│   └── f            #  >> Field: Face.
│       ├── f        #   > Face indices of cell 'i'.
│       ├── xy_v     #   > Face 'j' vertices of cell 'i'.  
│       ├── mean     #   > Face 'j' centroid: [Xv,Yv] = [mean(1,j),mean(2,j)].     
│       ├── len      #   > Face 'j' length of cell 'i'.
│       ├── Nf       #   > Face 'j' outer (unit) normal of cell 'i'.
│       └── Sf       #   > Face 'j' outer (face) normal of cell 'i'.
│
├── f                # >>> Field: Face.
│   ├── NF           #  >> Number of faces.
│   ├── xy_v         #  >> Face 'i' vertices: [Xv(:),Yv(:)] = [xy_v{i}(:,1),xy_v{i}(:,2)].
│   ├── mean         #  >> Cell 'i' centroid: [Xv,Yv] = [mean(1,i),mean(2,i)].
│   └── c            #  >> Face 'i' direct neighbouring cells.
│                                ├── 1. length(c{i}) = 1. Boundary face.
│                                └── 2. length(c{i}) = 2. Bulk face.
│
├── bnd              # >>> Field: Boundary.
│   ├── c            #  >> Field: Cell.
│   │   ├── c{1,i}   #   > Face 'j' vertices of boundary cell 'i': [Xv(:),Yv(:)] = [c{1,i}{j}(:,1),c{1,i}{j}(:,2)]. 
│   │   └── c{2,i}   #   > Cell index of boundary cell 'i'.
│   │
│   └── f            #  >> Field: Face.  
│       ├── f{1,i}   #   > Face vertices of boundary face 'i'. 
│       ├── f{2,i}   #   > Face index of boundary face 'i'.
│       ├── f{3,i}   #   > Cell index of boundary face 'i'.
│       ├── f{4,i}   #   > Boundary (char). 
│       └── f{5,i}   #   > Outer normal of boundary face 'i'.                 
│
└── s                # >>> Field: Stencil.
    ├── c            #  >> Face 'i' stencil cell indices.                           
    ├── f            #  >> Face 'i' stencil face indices.
    │                            └──  c{j,i}: Layer 'j' of face 'i' stencil.
    ├── c_e          #  >> Face 'i' extended stencil cell indices.                           
    ├── f_e          #  >> Face 'i' extended stencil face indices.
    │                            └──  c{j,i}: (Extension) layer 'j' of face 'i' stencil.
    ├── xy_v_c       #  >> Face 'i' (cell)  stencil points: [Xv(:),Yv(:)] = [xy_v_l{j,i}(1,:),xy_v_l{j,i}(2,:)]. 
    ├── xy_v_f       #  >> Face 'i' (face)  stencil points: [Xv(:),Yv(:)] = [xy_v_t  {i}(1,:),xy_v_t  {i}(2,:)]. 
    ├── xy_v_t       #  >> Face 'i' (total) stencil points: [Xv(:),Yv(:)] = [xy_v_t  {i}(1,:),xy_v_t  {i}(2,:)]. 
    └── par          #  >> Face 'i' stencil parameters.
        ├── n_e      #   > Face 'i' number of extensions.
        │                        ├── ne(1,i): number of extensions (x-direction).
        │                        └── ne(2,i): number of extensions (y-direction)
        ├── n_x      #   > Adimensional          length (x-direction).
        ├── n_y      #   > Adimensional          length (y-direction).
        ├── l_x      #   > Limit                        (x-direction): [l_x_min(i),l_x_max(i)] = [l_x(i,1),l_x(i,2)].
        └── l_y      #   > Limit                        (y-direction): [l_y_min(i),l_y_max(i)] = [l_y(i,1),l_y(i,2)].             ​
```

