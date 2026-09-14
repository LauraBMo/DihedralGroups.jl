"""
    DihedralGroups

Dihedral groups `D_d` — the symmetries of a regular `d`-gon — as a Julia type,
together with the group operation, the action on the vertices `0:d-1`, and
iteration over the whole group.

An element is either a rotation `rᵏ` or a reflection `rᵏ·s`, where `0 ≤ k < d`.
The group has order `2d` and the presentation

```
⟨ r, s | rᵈ = s² = id, s·r·s = r⁻¹ ⟩
```

`d` must be at least 1.

# Quick start
```jldoctest
julia> r(4, 1) * s(4)     # a quarter turn, then a reflection
D4: r^3·s

julia> 2^r(4, 1)          # rotate vertex 2 by one step
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

# See also
- [`Dihedral`](@ref) — the element type.
- [`r`](@ref) and [`s`](@ref) — the two generators.
- [`dihedralgroup`](@ref) — the iterator over a whole group.
"""
module DihedralGroups

export Dihedral, dhreflect, dhrotation, dihedralgroup, order,
       rotation_angle, isreflect, r, s

# The domain check for `d`, shared by every entry point that builds a group or
# an element, so the error message cannot drift between them.
@inline function _check_d(d)
    d ≥ 1 || throw(ArgumentError("the dihedral group D_d is only defined for d ≥ 1, got d = $d"))
    return nothing
end

"""
    Dihedral{d} <: Any

Dihedral group element representing a symmetry of a regular `d`-sided polygon.

# Fields
- `phi::Int`: Rotation component (0 ≤ phi < d)
- `reflect::Bool`: Reflection flag (true for reflection elements)

# Mathematical Representation
- Pure rotations: rᵏ where k = phi
- With reflections: rᵏ s where k = phi and s sends 0 to d-1

# Group Structure
The dihedral group D_d has order 2d with presentation:
⟨ r, s | rᵈ = s² = id, srs = r⁻¹ ⟩

# Examples
```jldoctest
julia> Dihedral{4}(1, false)  # 90° rotation
D4: r^1

julia> Dihedral{4}(0, true)   # Pure reflection
D4: s
```
"""
struct Dihedral{d}
    phi::Int
    reflect::Bool

    """
    Dihedral{d}(phi, reflect)

    Construct a dihedral group element with automatic modulo reduction of rotation index.

    Throws an `ArgumentError` if `d < 1`.
    """
    function Dihedral{d}(phi, reflect) where d
        _check_d(d)
        new{d}(mod(phi, d), reflect)
    end
end

"""
    Dihedral(d) -> Function

Group element constructor for dihedral group D_d.

Throws an `ArgumentError` if `d < 1`.

# Returns
- Constructor function `(i, reflect) -> Dihedral{d}(i, reflect)`

# Examples
```jldoctest
julia> pentagon = Dihedral(5);

julia> pentagon(2, false)  # 144° rotation
D5: r^2
```
"""
function Dihedral(d)
    _check_d(d)
    return (i, reflect) -> Dihedral{d}(i, reflect)
end

"""
    dhrotation(d::Int, i::Int) -> Dihedral{d}

Construct pure rotation element.

# Arguments
- `d`: Rotational symmetry order (`d ≥ 1`)
- `i`: Rotation steps (counterclockwise)

# Mathematical Form
rⁱ (rotation angle = (i/d) × 360°)

# Examples
```jldoctest
julia> dhrotation(6, 3)  # 180° rotation
D6: r^3
```
"""
dhrotation(d::Int, i::Int) = Dihedral{d}(i, false)

"""
    dhreflect(d::Int, i::Int) -> Dihedral{d}

Construct reflection element with rotational component.

# Arguments
- `d`: Polygon side count
- `i`: Pre-reflection rotation steps

# Mathematical Form
rⁱ s (rotation then reflection)

# Examples
```jldoctest
julia> dhreflect(4, 1)  # 90° rotation then reflection
D4: r^1·s
```
"""
dhreflect(d::Int, i::Int) = Dihedral{d}(i, true)

"""
    rotation_angle(g::Dihedral) -> Int

Rotation component of group element.

# Returns
- Rotation index k where 0 ≤ k < d

# Examples
```jldoctest
julia> rotation_angle(dhrotation(8, 3))
3
```
"""
rotation_angle(g::Dihedral) = g.phi

"""
    isreflect(g::Dihedral) -> Bool

Determine if element contains reflection.

# Returns
- `true` for reflection elements, `false` for pure rotations

# Examples
```jldoctest
julia> isreflect(s(5))
true
```
"""
isreflect(g::Dihedral) = g.reflect

"""
    r(d, phi = 1) -> Dihedral{d}

The rotation generator of `D_d`, and the constructor for any rotation element.

`r(d)` is the generator `r` itself: a single counterclockwise step of `360°/d`.
`r(d, phi)` is `rᵖʰⁱ`, a rotation by `phi` steps, reduced modulo `d`.

# Examples
```jldoctest
julia> r(4)
D4: r^1

julia> r(4, 3)      # three steps, i.e. 270°
D4: r^3

julia> r(4, 5)      # reduced modulo d
D4: r^1

julia> r(4, 1)^2
D4: r^2
```

# See also
- [`s`](@ref) — the reflection generator.
- [`dhrotation`](@ref) — the equivalent constructor.
"""
r(d, phi = 1) = Dihedral{d}(phi, false)

"""
    s(d, phi = 0) -> Dihedral{d}

The reflection generator of `D_d`, and the constructor for any reflection
element.

`s(d)` is the reflection `s` itself: it fixes vertex `0` and sends vertex `i` to
`d - i`. `s(d, phi)` is `rᵖʰⁱ·s`, that same reflection precomposed with a rotation
by `phi` steps.

Every reflection is its own inverse, so `s(d, phi)^2 == one(Dihedral{d})`.

# Examples
```jldoctest
julia> s(4)
D4: s

julia> s(4, 1)
D4: r^1·s

julia> s(4, 1)^2
D4: Id

julia> 0^s(5)       # s fixes vertex 0
0
```

# See also
- [`r`](@ref) — the rotation generator.
- [`dhreflect`](@ref) — the equivalent constructor.
"""
s(d, phi = 0) = Dihedral{d}(phi, true)

"""
    g::Dihedral{d} * h::Dihedral{d} -> Dihedral{d}

Group multiplication operation.

`g * h` applies `g` first and then `h`. Written out on the vertices, that makes
the action a left action:

```math
(i^g)^h = i^{g * h}
```

# Composition Rules
rᵃ * rᵇ = rᵃ⁺ᵇ
rᵃ * (rᵇ s) = rᵇ⁻ᵃ s
(rᵃ s) * rᵇ = rᵃ⁺ᵇ s
(rᵃ s) * (rᵇ s) = rᵇ⁻ᵃ

# Examples
```jldoctest
julia> r(4) * s(4)  # A quarter turn, then a reflection
D4: r^3·s

julia> s(4) * r(4)  # A reflection, then a quarter turn
D4: r^1·s
```
"""
function Base.:*(g::Dihedral{d}, h::Dihedral{d}) where {d}
    a = rotation_angle(g)
    b = rotation_angle(h)
    # `g * h` applies `g` first and then `h`, so the action below is a left
    # action: (i^g)^h == i^(g * h).
    return Dihedral{d}(isreflect(h) ? b - a : a + b, isreflect(g) ⊻ isreflect(h))
end

"""
    one(::Type{Dihedral{d}}) -> Dihedral{d}

Group identity element.

# Returns
- Identity element r⁰ (no rotation, no reflection)

# Examples
```jldoctest
julia> one(Dihedral{4})
D4: Id
```
"""
Base.one(::Type{Dihedral{d}}) where {d} = r(d, 0)

"""
    one(g::Dihedral{d}) -> Dihedral{d}

Identity element in same group.

# Examples
```jldoctest
julia> one(dhreflect(3,1))
D3: Id
```
"""
Base.one(g::Dihedral) = one(typeof(g))

"""
    inv(g::Dihedral{d}) -> Dihedral{d}

Group inverse element.

# Inverse Rules
(rᵏ)⁻¹ = rᵈ⁻ᵏ
(rᵏ s)⁻¹ = rᵏ s

# Examples
```jldoctest
julia> inv(r(5, 2))  # Inverse of 144° rotation
D5: r^3

julia> inv(s(4))     # Reflection is self-inverse
D4: s
```
"""
function Base.inv(g::Dihedral{d}) where {d}
   isreflect(g) ? g : r(d, d - rotation_angle(g))
end

"""
    i::Int ^ g::Dihedral{d} -> Int

Group action on polygon vertices.

# Action Rules
Pure rotation: i ↦ (i + k) mod d
Reflection: i ↦ (d - i + k) mod d

This is a left action: `(i^g)^h == i^(g * h)`.

# Arguments
- `i`: Vertex index (0 to d-1)
- `g`: Group element

# Examples
```jldoctest
julia> 1^r(4, 1)  # Rotate vertex 1 by 90°
2

julia> 3^s(5)     # Reflect vertex 3
2
```
"""
function Base.:^(i::Int, g::Dihedral{d}) where {d}
    j = isreflect(g) ? (d - i) : i
    return mod(j + rotation_angle(g), d)
end

"""
    g::Dihedral(n::Int) -> Int

Apply group action to vertex (alternative syntax).

# Examples
```jldoctest
julia> s = dhreflect(5, 2);

julia> s(3)  # Apply to vertex 3
4
```
"""
(g::Dihedral)(n::Int) = n^g

"""
    order(g::Dihedral{d}) -> Int

Element order (smallest k > 0 where gᵏ = identity).

# Order Calculation
Identity: 1
Pure rotation: d / gcd(k, d)
Reflection: 2

# Examples
```jldoctest
julia> order(r(8, 2))  # 90° rotation in octagon
4

julia> order(dhreflect(5, 0))
2
```
"""
function order(g::Dihedral{d}) where {d}
    o = 2 # g has a reflection
    if !isreflect(g)
        k = rotation_angle(g)
        o = iszero(k) ? 1 : div(d, gcd(k, d))
    end
    return o
end

"""
    g::Dihedral{d} ^ n::Int -> Dihedral{d}

Group element exponentiation.

# Algorithm
1. Compute element order m = order(g)
2. Reduce exponent modulo order: n_red = n mod m
3. Compute gⁿ_ʳᵉᵈ by repeated multiplication

# Examples
```jldoctest
julia> r(4, 1)^4  # Full rotation cycle
D4: Id

julia> s(3)^2      # Reflection squared
D3: Id
```
"""
function Base.:^(g::Dihedral{d}, n::Int) where {d}
    iszero(n) && return one(g)
    o = order(g)
    actual_n = mod(n, o)
    foldl((acc, _) -> acc * g, 1:actual_n; init=one(g))
end

"""
    show(io::IO, g::Dihedral{d})

Display element in algebraic notation.

# Representation
Identity: "Id"
Pure rotation: "r^k"
Pure reflection: "s"
Reflection+rotation: "r^k·s"

# Examples
```jldoctest
julia> Dihedral{3}(0, true)
D3: s

julia> Dihedral{4}(2, false)
D4: r^2
```
"""
function Base.show(io::IO, g::Dihedral{d}) where {d}
    print(io, "D$d: ")
    k = rotation_angle(g)
    if iszero(k)
        if isreflect(g)
            print(io, "s")
        else
            print(io, "Id")
        end
    else
        print(io, isreflect(g) ? "r^$k·s" : "r^$k")
    end
end

# Enable broadcasting
Base.broadcastable(g::Dihedral) = Ref(g)

"""
    DihedralGroup{d}

The iterator type over all `2d` elements of the dihedral group `D_d`, as returned
by [`dihedralgroup`](@ref).

`DihedralGroup` is not exported; construct one with `dihedralgroup(d)`.

# Element order
1. Rotations: r⁰, r¹, ..., rᵈ⁻¹
2. Reflections: r⁰s, r¹s, ..., rᵈ⁻¹s

# Examples
```jldoctest
julia> collect(dihedralgroup(3))
6-element Vector{Dihedral{3}}:
 D3: Id
 D3: r^1
 D3: r^2
 D3: s
 D3: r^1·s
 D3: r^2·s

julia> length(dihedralgroup(3))
6
```

# See also
- [`dihedralgroup`](@ref) — the constructor for this type.
"""
struct DihedralGroup{d} end

"""
    dihedralgroup(d) -> DihedralGroup{d}

Construct the iterator over all `2d` elements of `D_d` (`d ≥ 1`).

The rotations `r⁰ … rᵈ⁻¹` come first, then the reflections `r⁰s … rᵈ⁻¹s`.

Throws an `ArgumentError` if `d < 1`.

# Examples
```jldoctest
julia> collect(dihedralgroup(3))
6-element Vector{Dihedral{3}}:
 D3: Id
 D3: r^1
 D3: r^2
 D3: s
 D3: r^1·s
 D3: r^2·s

julia> length(dihedralgroup(5))
10
```

# See also
- [`DihedralGroup`](@ref) — the iterator type this returns.
"""
function dihedralgroup(d)
    _check_d(d)
    return DihedralGroup{d}()
end

# Iterator implementation
Base.eltype(_::DihedralGroup{d}) where d = Dihedral{d}
Base.length(_::DihedralGroup{d}) where d = 2*d

Base.iterate(D::DihedralGroup) = (one(eltype(D)), 1)
function Base.iterate(_::DihedralGroup{d}, status) where d
    (status > 2*d-1) && return nothing
    # Rotations r⁰ … rᵈ⁻¹ first, then reflections r⁰s … rᵈ⁻¹s.
    out = status > d-1 ? s(d) * r(d, status - d) : r(d, status)
    return out, status + 1
end

end
