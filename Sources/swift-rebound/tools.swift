import Foundation

let TINY = 1e-308

func tools_orbit_to_particle_err(G: Double, primary: Particle, m: Double, a: Double, e: Double, inc: Double, Omega: Double, omega: Double, f: Double) -> (particle: Particle, err: Int) {
    if (e == 1) {
        // Can't initialize a radial orbit with orbital elements.
        return (ParticleNAN(), 1)
    }
    if (e < 0) {
        // Eccentricity must be greater than or equal to zero.
        return (ParticleNAN(), 2)
    }
    if (e > 1) {
        if (a > 0) {
            // Bound orbit (a > 0) must have e < 1. 
            return (ParticleNAN(), 3)
        }
    }
    else {
        if (a < 0){
            // Unbound orbit (a < 0) must have e > 1.
			return (ParticleNAN(), 4)
        }
    }
    if (e*cos(f) < -1) {
        // Unbound orbit can't have f set beyond the range allowed by the asymptotes set by the parabola.
        return (ParticleNAN(), 5)
    }
    if (primary.m < TINY) {
        // Primary has no mass
        return (ParticleNAN(), 6)
    }
    return (Particle(), 0)

}