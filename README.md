# DihedralGroups [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://LauraBMo.github.io/DihedralGroups.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://LauraBMo.github.io/DihedralGroups.jl/dev/) [![Build Status](https://github.com/LauraBMo/DihedralGroups.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/LauraBMo/DihedralGroups.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/LauraBMo/DihedralGroups.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/LauraBMo/DihedralGroups.jl)

Finite dihedral groups as a Julia type.

`D_d` is the symmetry group of a regular `d`-gon: the rotations `rᵏ` together
with the reflections `rᵏ·s`, for `0 ≤ k < d`. The group has order `2d` and the
presentation

```
⟨ r, s | rᵈ = s² = id, s·r·s = r⁻¹ ⟩
```

This package gives you the group operation, the action on the vertices `0:d-1`,
and iteration over every element of the group. `d ≥ 1`; there are no other
dependencies.

## Installation

```julia
using Pkg
Pkg.add("DihedralGroups")
```

## Usage

```julia
using DihedralGroups

g = r(4, 1)          # a quarter turn      -> D4: r^1
h = s(4)             # a reflection        -> D4: s

g * h                # D4: r^3·s  — g then h
g^4                  # D4: Id  — g has order 4
inv(g)               # D4: r^3

2^g                  # 3 — g sends vertex 2 to vertex 3
0^h                  # 0 — h fixes vertex 0

order(g)             # 4
rotation_angle(g)    # 1
isreflect(h)         # true

collect(dihedralgroup(3))   # all 6 elements of D3
```

Elements of the same group multiply, invert and compare; the elements of `D_d`
form a group under `*`, with `one(Dihedral{d})` as the identity.

## Documentation

See [the manual](https://LauraBMo.github.io/DihedralGroups.jl/stable/) for the
full API.
