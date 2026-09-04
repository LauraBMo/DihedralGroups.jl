# DihedralGroups.jl - Complete Implementation
module DihedralGroups

export Dihedral, dhreflect, dhrotation, dihedralgroup, order,
       rotation_angle, isreflect, r, s

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
    """
    Dihedral{d}(phi, reflect) where d = new(mod(phi, d), reflect)
end

"""
    Dihedral(d) -> Function

Group element constructor for dihedral group D_d.

# Returns
- Constructor function `(i, reflect) -> Dihedral{d}(i, reflect)`

# Examples
```jldoctest
julia> pentagon = Dihedral(5);

julia> pentagon(2, false)  # 144° rotation
D5: r^2
```
"""
Dihedral(d) = (i, reflect) -> Dihedral{d}(i, reflect)

"""
    dhrotation(d::Int, i::Int) -> Dihedral{d}

Construct pure rotation element.

# Arguments
- `d`: Polygon side count (≥3)
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

r(d, phi = 1) = Dihedral{d}(phi, false)
s(d, phi = 0) = Dihedral{d}(phi, true)

"""
    g::Dihedral{d} * h::Dihedral{d} -> Dihedral{d}

Group multiplication operation.

# Composition Rules
rᵃ * rᵇ = rᵃ⁺ᵇ
rᵃ * (rᵇ s) = rᵃ⁺ᵇ s
(rᵃ s) * rᵇ = rᵃ⁻ᵇ s
(rᵃ s) * (rᵇ s) = rᵃ⁻ᵇ

# Examples
```jldoctest
julia> r(4) * s(4)  # Rotation then reflection
D4: r^1·s

julia> s(4) * r(4)  # Reflection then rotation
D4: r^3·s
```
"""
function Base.:*(g::Dihedral{d}, h::Dihedral{d}) where {d}
    new_phi = rotation_angle(g) + rotation_angle(h)
    new_reflect = h.reflect
    if isreflect(g)
       new_phi = rotation_angle(g) - rotation_angle(h)
       new_reflect = !h.reflect
    end
    return Dihedral{d}(new_phi, new_reflect)
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

Iterator over all elements of dihedral group D_d.

# Element Order
1. Rotations: r⁰, r¹, ..., rᵈ⁻¹
2. Reflections: r⁰s, r¹s, ..., rᵈ⁻¹s

# Examples
```jldoctest
julia> collect(DihedralGroup{3})
6-element Vector{Dihedral{3}}:
 D3: Id
 D3: r^1
 D3: r^2
 D3: s
 D3: r^1·s
 D3: r^2·s
```
"""
struct DihedralGroup{d} end

"""
    dihedralgroup(d) -> DihedralGroup{d}

Construct group iterator for d ≥ 3.

# Examples
```jldoctest
julia> for g in dihedralgroup(3)
           println(g)
       end
D3: Id
D3: r^1
D3: r^2
D3: s
D3: r^1·s
D3: r^2·s
```
"""
dihedralgroup(d) = DihedralGroup{d}()

# Iterator implementation
Base.eltype(_::DihedralGroup{d}) where d = Dihedral{d}
Base.length(_::DihedralGroup{d}) where d = 2*d

Base.iterate(D::DihedralGroup) = (one(eltype(D)), 1)
function Base.iterate(_::DihedralGroup{d}, status) where d
    (status > 2*d-1) && return nothing
    out = r(d, status)
    if status > d-1
        out *= s(d)
    end
    return out, status + 1
end

end
