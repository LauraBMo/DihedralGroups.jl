# using Revise ## Use when some test do not pass.
using DihedralGroups
using Test

using DihedralGroups: Dihedral

@testset "Dihedral Group Tests" begin
    @testset "D₄ (Square)" begin
        d = 4
        r = Dihedral{d}(1, false)  # r^1
        s = Dihedral{d}(0, true)   # s
        id = one(Dihedral{d}) # r^0

        @testset "r & s" begin
            @test all((r^phi)*s == Dihedral{d}(phi, true) for phi in 0:(d-1))
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
            @test all((r^phi)*s == Dihedral{d}(phi, true) for phi in 0:(d-1))
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
end
