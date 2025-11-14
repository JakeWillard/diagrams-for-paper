
struct FieldLine <: Mark
    mu :: Float64
end


struct FieldSheet <: Mark
    V :: Float64
end


struct WorldLine <: Mark
    V :: Float64
end


struct FiniteDifferenceLines <: Mark
    dx :: Float64 
    dy :: Float64
    center :: TextMark
    top :: TextMark
    bottom :: TextMark
    side :: TextMark
end


struct ZoomInRight <: Mark
    mu :: Float64
    scale :: Float64
end