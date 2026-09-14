# using Revise ## Use when some test do not pass.
using DihedralGroups
using Test
using Aqua

@testset "DihedralGroups" begin

    @testset "Aqua" begin
        Aqua.test_all(DihedralGroups)
    end

    @testset "D₄ (Square)" begin
        d = 4
        r = Dihedral{d}(1, false)  # r^1
        s = Dihedral{d}(0, true)   # s
        id = one(Dihedral{d}) # r^0

        @testset "r & s" begin
            # `g * h` applies g then h, so s * r^phi is the reflection r^phi·s.
            @test all(s * (r^phi) == Dihedral{d}(phi, true) for phi in 0:(d-1))
        end

        @testset "Basic Properties" begin
            @test order(r) == 4
            @test order(s) == 2
            @test order(r*r) == 2
            @test inv(r) == Dihedral{d}(3, false)
            @test inv(s) == s
        end

        @testset "Group Operations" begin
            @test r * inv(r) == id
            @test s * s == id
            @test (r*s)^2 == id
            @test r*s != s*r
        end

        @testset "Acting on integers" begin
            I = 0:(d-1)
            @test I.^(s) == [0, 3, 2, 1]
            @test I.^(r) == [1, 2, 3, 0]
            @test I.^(r^2) == [2, 3, 0, 1]
            @test I.^(r^3) == [3, 0, 1, 2]
        end

    end

    @testset "D₅ (Pentagon)" begin
        d = 5
        r = Dihedral{d}(1, false)  # r^1
        s = Dihedral{d}(0, true)   # s
        id = one(Dihedral{d})

        @testset "r & s" begin
            # `g * h` applies g then h, so s * r^phi is the reflection r^phi·s.
            @test all(s * (r^phi) == Dihedral{d}(phi, true) for phi in 0:(d-1))
        end

        @testset "Element Orders" begin
            @test order(r) == 5
            @test order(r^2) == 5
            @test order(r*s) == 2
        end

        @testset "Inverse Operations" begin
            @test inv(r*s) == r*s
            @test inv(r^2) == Dihedral{d}(3, false)
        end

        @testset "Acting on integers" begin
            I = 0:(d-1)
            @test I.^(s) == [0, 4, 3, 2, 1]
            @test I.^(r) == [1, 2, 3, 4, 0]
            @test I.^(r^2) == [2, 3, 4, 0, 1]
            @test I.^(r^3) == [3, 4, 0, 1, 2]
            @test I.^(r^4) == [4, 0, 1, 2, 3]
        end
    end

    @testset "Order for D3..9" begin
        d = 9
        @test all(3:d) do i # Test order
            id = one(Dihedral{i}) # r^0
            all(all(s^j != id for j in 1:(order(s)-1)) for s in dihedralgroup(i))
        end
        @test all(3:d) do i # Test order
            id = one(Dihedral{i}) # r^0
            all(s^(order(s)) == id for s in dihedralgroup(i))
        end
    end

    @testset "Powers by integers" begin
        d = 24
        r = Dihedral{d}(1, false)  # r^1
        s = Dihedral{d}(0, true)   # s
        id = one(Dihedral{d}) # r^0
        @test all(3:d) do k
            rotations = r.^(1:(d-1))
            rotation_angle.((rotations).^(k)) == [mod(k*i, d) for i in 1:(d-1)]
        end
    end

    @testset "Group axioms over a whole group" begin
        for d in 1:12
            G = collect(dihedralgroup(d))
            id = one(Dihedral{d})
            @test length(G) == length(dihedralgroup(d)) == 2d
            @test eltype(dihedralgroup(d)) === Dihedral{d}
            @test all(g isa Dihedral{d} for g in G)
            @test length(unique(G)) == 2d                          # all distinct
            @test id in G
            @test all(g * id == g == id * g for g in G)            # identity
            @test all(g * inv(g) == id == inv(g) * g for g in G)   # inverses
            @test all(inv(inv(g)) == g for g in G)                 # involution
            @test all(g * h in G for g in G, h in G)               # closure
            @test all((g * h) * k == g * (h * k)
                      for g in G, h in G, k in G)                  # associativity
        end
    end

    @testset "Generators and element order" begin
        for d in 1:12
            R = r(d)
            S = s(d)
            id = one(Dihedral{d})
            @test order(id) == 1
            @test order(R) == d
            @test order(S) == 2
            @test all(order(g) == 2 for g in dihedralgroup(d) if isreflect(g))
            @test all(g^order(g) == id for g in dihedralgroup(d))
        end
    end

    @testset "Integer action is a group action" begin
        # `g * h` applies g first, so the action is a LEFT action.
        for d in 1:12
            for g in dihedralgroup(d), h in dihedralgroup(d), i in 0:(d-1)
                @test (i^g)^h == i^(g * h)
            end
            # Every element permutes the vertices.
            @test all(sort([i^g for i in 0:(d-1)]) == collect(0:(d-1))
                      for g in dihedralgroup(d))
        end
    end

    @testset "Modulo reduction and constructors" begin
        @test Dihedral{4}(5, false) == r(4, 1)
        @test Dihedral{4}(-1, false) == r(4, 3)
        @test Dihedral{4}(9, true) == s(4, 1)
        @test dhrotation(6, 3) == r(6, 3)
        @test dhreflect(6, 3) == s(6, 3)
        @test dhrotation(6, 9) == r(6, 3)      # reduced on construction

        pentagon = Dihedral(5)                 # curried constructor
        @test pentagon(2, false) == dhrotation(5, 2)
        @test pentagon(1, true) == dhreflect(5, 1)
    end

    @testset "Powers by zero and negatives" begin
        for d in 1:12
            R = r(d)
            id = one(Dihedral{d})
            @test R^0 == id
            @test R^(-1) == inv(R)
            @test R^(-3) == inv(R)^3
            @test all(R^k == R^(k + d) for k in 0:(d-1))
            @test all(g^0 == id for g in dihedralgroup(d))
            @test all(g^(-1) == inv(g) for g in dihedralgroup(d))
        end
    end

    @testset "one, rotation_angle and isreflect" begin
        @test one(Dihedral{4}) == r(4, 0)
        @test one(r(4, 3)) == one(Dihedral{4})
        @test !isreflect(one(Dihedral{4}))
        @test rotation_angle(one(Dihedral{4})) == 0
        @test rotation_angle(Dihedral{4}(7, true)) == 3
        @test isreflect(s(4, 2))
    end

    @testset "show" begin
        @test sprint(show, Dihedral{4}(0, false)) == "D4: Id"
        @test sprint(show, Dihedral{4}(2, false)) == "D4: r^2"
        @test sprint(show, Dihedral{4}(0, true)) == "D4: s"
        @test sprint(show, Dihedral{4}(3, true)) == "D4: r^3·s"
    end

    @testset "Broadcasting" begin
        @test Base.broadcastable(r(5)) isa Ref
        @test (0:4) .^ r(5) == [1, 2, 3, 4, 0]
    end

    @testset "Small d" begin
        @test r(1) == one(Dihedral{1})
        @test length(dihedralgroup(1)) == 2
        @test length(dihedralgroup(2)) == 4
        @test collect(dihedralgroup(1)) ==
              [Dihedral{1}(0, false), Dihedral{1}(0, true)]
        @test collect(dihedralgroup(2)) ==
              [Dihedral{2}(0, false), Dihedral{2}(1, false),
               Dihedral{2}(0, true), Dihedral{2}(1, true)]
    end

    @testset "d = 0 is rejected" begin
        # Every element construction funnels through the Dihedral{d} inner
        # constructor, so that is the single place the domain is checked.
        @test_throws ArgumentError Dihedral{0}(0, false)
        @test_throws ArgumentError Dihedral{0}(1, true)
        @test_throws ArgumentError r(0)
        @test_throws ArgumentError s(0)
        @test_throws ArgumentError dhrotation(0, 1)
        @test_throws ArgumentError dhreflect(0, 1)
        @test_throws ArgumentError Dihedral(0)
        @test_throws ArgumentError one(Dihedral{0})
        @test_throws ArgumentError dihedralgroup(0)
        @test_throws ArgumentError length(dihedralgroup(0))
    end
end
