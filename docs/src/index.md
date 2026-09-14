```@meta
CurrentModule = DihedralGroups
```

# DihedralGroups

`DihedralGroups.jl` implements the finite dihedral group

```math
D_d = \langle r, s \mid r^d = s^2 = 1,\; srs = r^{-1} \rangle
```

— the symmetry group of a regular `d`-gon — as a Julia type, together with the
group operation, the action on the vertices `0:d-1`, and iteration over the
whole group.

```@docs
DihedralGroups
```

## Quick start

```julia
julia> using DihedralGroups

julia> g = r(4, 1);        # a quarter turn

julia> h = s(4);           # a reflection

julia> g * h               # g, then h
D4: r^3·s

julia> g^4                 # back to the identity
D4: Id

julia> 2^g                 # g acts on the vertices 0:d-1
3

julia> collect(dihedralgroup(3))
6-element Vector{Dihedral{3}}:
 D3: Id
 D3: r^1
 D3: r^2
 D3: s
 D3: r^1·s
 D3: r^2·s
```

## Index

```@index
```
