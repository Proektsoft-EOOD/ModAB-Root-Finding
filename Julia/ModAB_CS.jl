function modab_CS(f, x1::Real, x2::Real, y::Real=0.0; xtol::Float64=1e-14, ytol::Float64=1e-16, maxIter::Int=200)
    if x1 > x2
        x1, x2 = x2, x1
    end
    epsy = ytol * max(1, abs(y))
    y1 = f(x1) - y
    abs(y1) <= epsy && return x1
    y2 = f(x2) - y
    abs(y2) <= epsy && return x2
    epsx = xtol * (x2 - x1)
    side = 0
    bisection = true
    C = 16 # safety factor for threshold corresponding to 4 iterations = 2^4
    threshold = x2 - x1  # Threshold to fall back to bisection if AB fails to shrink the interval enough
    # calculate k on each bisection step with account for local function properties and symmetry
    for i in 1:maxIter
        local x3 = bisection ? (x1 + x2) / 2 : (x1 * y2 - y1 * x2) / (y2 - y1)
        if x2 - x1 <= epsx # x-convergence check
            return x3
        end
        local y3
        if bisection
            y3 = f(x3) - y  # Function value at midpoint
            ym = (y1 + y2) / 2 # Ordinate of chord at midpoint
            r = 1 - abs(ym / (y2 - y1)) # Symmetry factor
            k = r * r # Deviation factor
            # Check if the function is close enough to linear
            if abs(ym - y3) < k * (abs(y3) + abs(ym))
                bisection = false
                threshold = (x2 - x1) * C
            end
        else
            if x3 <= x1
                x3, y3 = x1, y1
            elseif x3 >= x2
                x3 , y3 = x2, y2
            else
                y3 = f(x3) - y
            end    
            threshold /= 2
        end

        if abs(y3) <= epsy # y-convergence check
            return x3
        end

        if sign(y1) == sign(y3)
            if side == 1
                m = 1 - y3 / y1
                y2 *= m <= 0 ? 0.5 : m
            elseif !bisection
                side = 1
            end
            x1, y1 = x3, y3
        else
            if side == -1
                m = 1 - y3 / y2
                y1 *=  m <= 0 ? 0.5 : m
            elseif !bisection
                side = -1
            end
            x2, y2 = x3, y3
        end
        if x2 - x1 > threshold # AB failed to shrink the interval enough
            bisection = true
            side = 0
        end
    end
    return NaN
end