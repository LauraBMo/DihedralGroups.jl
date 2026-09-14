```@meta
CurrentModule = DihedralGroups
```

# API

## The group element

```@docs
Dihedral
```

## Generators

```@docs
r
s
```

## Constructors

```@docs
dhrotation
dhreflect
```

## Group operations

These extend the corresponding `Base` functions, so they are listed under their
`Base` names. The signature must be written to match the method exactly: `{d}`
appears below wherever the method itself is declared with `where {d}`, and is
absent for `one(g::Dihedral)`, which is not.

```@docs
Base.:*(::Dihedral{d}, ::Dihedral{d}) where {d}
Base.inv(::Dihedral{d}) where {d}
Base.:^(::Dihedral{d}, ::Int) where {d}
Base.:^(::Int, ::Dihedral{d}) where {d}
Base.one(::Type{Dihedral{d}}) where {d}
Base.one(::Dihedral)
Base.show(::IO, ::Dihedral{d}) where {d}
```

## Queries

```@docs
order
rotation_angle
isreflect
```

## Iteration

```@docs
dihedralgroup
DihedralGroup
```
