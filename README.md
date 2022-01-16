# 2D Code structure
## Top-level directory layout

```
'
├── [Tools - Data]
|   └── ...
|
├── [Tools - Numerical]
|   └── ...
|
├── [Tools - Post-processing]
|   ├── Fig_0
|   ├── Fig_1
|   ├── Fig_2
|   └── ...
│
└── src  
    |
    ├── Main                        # >>> Call 'Class' functions.
    |
    ├── Class_1                     # >>> Call 'SubClass' functions.
    |   ├── SubClass_1_1            #  >> Set inputs.
    |   ├── SubClass_1_2            #  >> Grid generation. 
    |   └── SubClass_1_3            #  >> Call 'SubSubClass' functions.
    |       ├── SubSubClass_1_3_1   #   > Compute grid properties.
    |       └── SubSubClass_1_3_2   #   > Compute stencil elements (w/o and w/ extension(s)).
    |
    ├── Class_2                     # >>> Call 'SubClass' functions.
    |   ├── SubClass_2_1            #  >> Call 'SubSubClass' functions.
    |   |   ├── SubSubClass_2_1_1   #   > Compute analytic solution.
    |   |   └── SubSubClass_2_1_2   #   > 1D/2D quadrature.
    |   ├── SubClass_2_2            #  >> Call 'SubSubClass' functions.
    |   |   ├── SubSubClass_2_2_1   #   > Assemble matrices Df, Dwf, Pf and Tf.
    |   |   └── SubSubClass_2_2_2   #   > Assemble matrices A, B and C.
    |   └── SubClass_2_3            #  >> Compute error norms.
    |
    ├── Class_3                     # >>> Call 'SubClass' functions.
    |   ├── SubClass_3_1            #  >> 
    |   └── ...                     #  >> 
    |
    └── Class_4                     # >>> Call functions on ../[Tools - Post-processing].
```
```
'
inp 
|
├── X                      
|
├── Y
|
└── Z
```
```
'
msh 
|
├── bnd              # >>> Field: Boundary.
|   ├── c            #  >> Field: Cell.
|   |   ├── c{1,i}   #   > Face 'j' vertices of boundary cell 'i': [Xv(:),Yv(:)] = [c{1,i}{j}(:,1),c{1,i}{j}(:,2)]. 
|   |   └── c{2,i}   #   > Cell index of boundary cell 'i'.
|   |
|   └── f            #  >> Field: Face.  
|       ├── f{1,i}   #   > Face vertices of boundary face 'i'. 
|       ├── f{2,i}   #   > Face index of boundary face 'i'.
|       ├── f{3,i}   #   > Cell index of boundary face 'i'.
|       ├── f{4,i}   #   > Boundary (char). 
|       └── f{5,i}   #   > Outer normal of boundary face 'i'.                  
|
├── c                # >>> Field: Cell.
|   ├── NC           #   > Number of cells.
|   ├── xy_v         #   > Cell 'i' vertices: [Xv(:),Yv(:)] = [xy_v{i}(:,1),xy_v{i}(:,2)].
|   ├── vol          #   > Cell 'i' volume.
|   ├── mean         #   > Cell 'i' centroid: [Xv,Yv] = [mean(1,i),mean(2,i)].
|   ├── h            #   > Cell 'i' reference length/hydraulic diameter. 
|   ├── h_ref        #   > Reference length/hydraulic diameter. 
|   ├── nb           #   > Cell 'i' neighbouring cell indices. 
|   └── f            #  >> Field: Face.
|       ├── faces    #   > Face indices of cell 'i'.
|       ├── xy_v     #   > Face 'j' vertices of cell 'i'.
|       ├── len      #   > Face 'j' length of cell 'i'.
|       ├── mean     #   > Face 'j' centroid of cell 'i'.
|       └── Nf       #   > Face 'j' outer normal of cell 'i'.
|
├── d                # >>> Field: Domain.
|   └── xy_v         #  >> Domain vertices: [Xv(:),Yv(:)] = [xy_v(:,1),xy_v(:,2)]. 
|
├── f                # >>> Field: Face.
|   ├── NF           #   > Number of faces.
|   ├── xy_v         #   > Face 'i' vertices: [Xv(:),Yv(:)] = [xy_v{i}(:,1),xy_v{i}(:,2)].
|   ├── mean         #   > Cell 'i' centroid: [Xv,Yv] = [mean(1,i),mean(2,i)].
|   ├── len          #   > Face 'i' length.
|   └── c            #   > Face 'i' direct neighbouring cells.
|                                ├── 1. length(c{i}) = 1. Boundary face.
|                                └── 2. length(c{i}) = 2: Bulk face.
|
└── s                # >>> Field: Stencil.
    ├── c            #   > Face 'i' stencil cell indices.                           
    ├── f            #   > Face 'i' stencil face indices.
    |                            └──  c{j,i}: Layer 'j' of face 'i' stencil.
    ├── c_e          #   > Face 'i' extended stencil cell indices.                           
    ├── f_e          #   > Face 'i' extended stencil face indices.
    |                            └──  c{j,i}: (Extension) layer 'j' of face 'i' stencil.
    ├── xy_v_l       #   > Face 'i' (level) stencil points (cell/face centroid): [Xv(:),Yv(:)] = [xy_v_l{j,i}(1,:),xy_v_l{j,i}(2,:)]. 
    ├── xy_v_t       #   > Face 'i' (total) stencil points (cell/face centroid): [Xv(:),Yv(:)] = [xy_v_t  {i}(1,:),xy_v_t  {i}(2,:)]. 
    └── par          #  >> Face 'i' stencil parameters.
        ├── n_e      #   > Face 'i' number of extensions.
        |                        ├── ne(1,i): number of extensions (x-direction).
        |                        └── ne(2,i): number of extensions (y-direction)
        ├── n_x      #   > Adimensional length (x-direction).
        ├── n_y      #   > Adimensional length (y-direction).
        ├── h_x      #   > Reference    length (x-direction).
        ├── h_y      #   > Reference    length (y-direction). 
        ├── l_x      #   > Limit               (x-direction): [l_x_min(i),l_x_max(i)] = [l_x(i,1),l_x(i,2)].
        └── l_y      #   > Limit               (y-direction): [l_y_min(i),l_y_max(i)] = [l_y(i,1),l_y(i,2)].             ​
```

