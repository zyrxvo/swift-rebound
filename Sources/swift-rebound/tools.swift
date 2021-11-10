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
    let p = Particle(m: m)
    let r = a*(1 - e*e) / (1 + e*cos(f))
    let v0 = sqrt(G*(m+primary.m)/a/(1.0 - e*e)) // in this form it works for elliptical and hyperbolic orbits

	let cO = cos(Omega)
	let sO = sin(Omega)
	let co = cos(omega)
	let so = sin(omega)
	let cf = cos(f)
	let sf = sin(f)
	let ci = cos(inc)
	let si = sin(inc)
	
	// Murray & Dermott Eq 2.122
	p.x = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
	p.y = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
	p.z = primary.z + r*(so*cf+co*sf)*si

	// Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
	p.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
	p.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
	p.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)
	
	p.ax = 0; 	p.ay = 0; 	p.az = 0;

    return (p, 0)

}

func tools_particle_to_orbit_err(G: Double, p: Particle, primary: Particle) -> (orbit: Orbit, err: Int) {

    return (Orbit(), 0)
}