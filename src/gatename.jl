# GateName definition
# -------------------
struct GateName{Name} end

"""
GateName is a parameterized type which allows making strings into Julia types for the
purpose of representing gate names.  The main use of GateName is overloading the `gate`
method which generates operators for Qubit and vQubit sites.

To make a GateName type, you can use the string macro notation: `GateName"MyTag"`.
To make an GateName value or object, you can use the notation: `GateName("mygate")`.
"""
GateName(s::AbstractString) = GateName{Symbol(s)}()
GateName(s::Symbol) = GateName{s}()
name(::GateName{N}) where {N} = N

macro GateName_str(s)
    return GateName{Symbol(s)}
end
