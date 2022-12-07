using FFTW

struct Spectrum{S, F}
    spec :: S
    freq :: F
end

function power_spectrum_1d_isotropic(var, x, y; condition = nothing)
    
    if !(condition isa Nothing)
        return Spectrum(nothing, nothing)
    end

    Nx = length(x)
    Ny = length(y)
    Nfx = Int64(Nx)
    Nfy = Int64(Ny)
    
    freq_x    = zeros(Float64, Nfx, Nfy)
    freq_y    = zeros(Float64, Nfy, Nfx)
   
    dx = x[2] - x[1]
    for j in 1:Nfy
        freq_x[:, j] = fftshift(fftfreq(Nfx, 1.0 / dx)) .* 2.0 .* π
    end
    dy = y[2] - y[1]
    for i in 1:Nfx
        freq_y[i, :] = fftshift(fftfreq(Nfy, 1.0 / dy)) .* 2.0 .* π
    end
    
    four_ful  = fftshift(fft(var)) / Nfx / Nfy
    spectra2D = abs.(four_ful .* four_ful)
    
    # Convert to 1D (shell integral over bins)
    freq = sqrt.(freq_x.^2 + freq_y.^2)

    wavef = range(0.0, maximum(freq) + 0.1, length = Int(min(Nfx/2, Nfy/2)+1))
    wave = zeros(length(wavef)-1)
    for i in 1:length(wave)
        wave[i] = 0.5*(wavef[i] + wavef[i+1])
    end

    spectra = zeros(Float64, length(wave))
    for i in 1:Nfx, j in 1:Nfy
        for k in 2:length(wave)
            if (freq[i, j] < wave[k]) && (freq[i, j] >= wave[k-1])
                spectra[k] += spectra2D[i, j]
            end
        end
    end

    return Spectrum(spectra, wave)
end

function power_spectrum_1d_x(var, x, y; condition = nothing)
    Nx = length(x)
    Ny = length(y)
    Nfx = Int64(Nx)
    
    spectra = zeros(Float64, Int(Nfx/2))
    
    dx = x[2] - x[1]
    freqs = fftfreq(Nfx, 1.0 / dx) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfx/2)] .* 2.0 .* π
    
    for j in 1:Ny
        if  (condition isa Nothing) || condition(π, y[j])
            fourier      = fft(var[:, j]) / Nfx
            spectra[1]  += fourier[1] .* conj(fourier[1])

            for m in 2:Int(Nfx/2)
                spectra[m] += 2.0 * fourier[m] * conj(fourier[m])  / Ny # factor 2 for neg freq contribution
            end
        end
    end
    return Spectrum(spectra, freqs)
end

    
function power_spectrum_1d_y(var, x, y; condition = nothing)
    Nx = length(x)
    Ny = length(y)
    Nfy = Int64(Ny)
    
    spectra = zeros(Float64, Int(Nfy/2))
    
    dy = y[2] - y[1]
    freqs = fftfreq(Nfy, 1.0 / dy) # 0,+ve freq,-ve freqs (lowest to highest)
    freqs = freqs[1:Int(Nfy/2)] .* 2.0 .* π
    
    for i in 1:Nx
        if  (condition isa Nothing) || condition(x[i], π/2)

            fourier      = fft(var[i, :]) / Nfy
            spectra[1]  += fourier[1] .* conj(fourier[1])

            for m in 2:Int(Nfy/2)
                spectra[m] += 2.0 * fourier[m] * conj(fourier[m]) / Nx # factor 2 for neg freq contribution
            end
        end
    end
    return Spectrum(spectra, freqs)
end

function immersed_condition(x, y)
    @inline toplft(x, y) = (((x > π/2) & (x < 3π/2)) & (((y > π/3) & (y < 2π/3)) | ((y > 4π/3) & (y < 5π/3))))
    @inline botlft(x, y) = (((x > π/2) & (x < 3π/2)) & (((y < -π/3) & (y > -2π/3)) | ((y < -4π/3) & (y > -5π/3))))
    @inline toprgt(x, y) = (((x < -π/2) & (x > -3π/2)) & (((y > π/3) & (y < 2π/3)) | ((y > 4π/3) & (y < 5π/3))))
    @inline botrgt(x, y) = (((x < -π/2) & (x > -3π/2)) & (((y < -π/3) & (y > -2π/3)) | ((y < -4π/3) & (y > -5π/3))))

    return !(toplft(x, y) | toprgt(x, y) | botlft(x, y) | botrgt(x, y))
end
