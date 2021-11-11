import Foundation

let TINY = 1e-308
let M_PI = 3.141592653589793238462643383279502884197
let MIN_REL_ERROR = 1.0e-12
// Close to smallest relative floating point number, used for orbit calculation

let MIN_INC = 1.0e-8
// Below this inclination, the broken angles pomega and theta equal the corresponding
// unbroken angles to within machine precision, so a practical boundary for planar orbits

let MIN_ECC = 1.0e-8
// Below this eccentricity, corrections at order e^2 are below machine precision, so we use
// stable expressions accurate to O(e) for the mean longitude below for near-circular orbits.
// returns acos(num/denom), using disambiguator to tell which quadrant to return.  
// will return 0 or pi appropriately if num is larger than denom by machine precision
// and will return 0 if denom is exactly 0.0

func tools_mod_twopi(f: Double) -> Double {
    let twopi = 2.0*M_PI;
    return fmod(twopi + fmod(f, twopi), twopi);
}

func tools_orbit_to_particle(G: Double, primary: Particle, m: Double, orbit: Orbit) throws -> Particle {
    return try tools_orbital_elements_to_particle(G: G, primary: primary, m: m, a: orbit.a, e: orbit.e, inc: orbit.inc, Omega: orbit.Omega, omega: orbit.omega, f: orbit.f)
}

func tools_orbital_elements_to_particle(G: Double, primary: Particle, m: Double, a: Double, e: Double, inc: Double, Omega: Double, omega: Double, f: Double) throws -> Particle {
    guard e != 1 else { throw OrbitError.radialOrbit }
    guard e >= 0 else { throw OrbitError.negativeEccentricity }
    if (e > 1) {
        guard a < 0 else { throw OrbitError.boundOrbitError }
    }
    else {
        guard a > 0 else { throw OrbitError.unboundOrbitError }
    }
    guard (e*cos(f) >= -1) else { throw OrbitError.invalidTrueAnomaly }
    guard (primary.m >= TINY) else { throw OrbitError.masslessPrimary }

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

    return p

}

func acos2(num: Double, denom: Double, disambiguator: Double) -> Double {
	var val: Double
	let cosine = num/denom
	if(cosine > -1.0 && cosine < 1.0){
		val = acos(cosine)
		if(disambiguator < 0.0){
			val = -val
		}
	}
	else{
		val = (cosine <= -1.0) ? M_PI : 0.0
	}
	return val;
}

func tools_particle_to_orbit(G: Double, p: Particle, primary: Particle) throws -> Orbit {
    let o = Orbit()
    
    guard (primary.m >= TINY) else { throw OrbitError.masslessPrimary }
    
	var mu,dx,dy,dz,dvx,dvy,dvz,vsquared,vcircsquared,vdiffsquared : Double
	var hx,hy,hz,vr,rvr,muinv,ex,ey,ez,nx,ny,n,ea : Double
	mu = G*(p.m+primary.m)
	dx = p.x - primary.x
	dy = p.y - primary.y
	dz = p.z - primary.z
	dvx = p.vx - primary.vx
	dvy = p.vy - primary.vy
	dvz = p.vz - primary.vz
	o.d = sqrt ( dx*dx + dy*dy + dz*dz )
	
	vsquared = dvx*dvx + dvy*dvy + dvz*dvz
	o.v = sqrt(vsquared)
	vcircsquared = mu/o.d	
	o.a = -mu/( vsquared - 2.0*vcircsquared ) // semi major axis
    
    o.rhill = o.a*cbrt(p.m/(3.0*primary.m))
	
	hx = (dy*dvz - dz*dvy) //angular momentum vector
	hy = (dz*dvx - dx*dvz)
	hz = (dx*dvy - dy*dvx)
	o.h = sqrt ( hx*hx + hy*hy + hz*hz ) // abs value of angular momentum

	vdiffsquared = vsquared - vcircsquared;
    guard (o.d > TINY) else { throw OrbitError.coincidesWithPrimary }

	vr = (dx*dvx + dy*dvy + dz*dvz)/o.d	
	rvr = o.d*vr
	muinv = 1.0/mu

	ex = muinv*( vdiffsquared*dx - rvr*dvx )
	ey = muinv*( vdiffsquared*dy - rvr*dvy )
	ez = muinv*( vdiffsquared*dz - rvr*dvz )
 	o.e = sqrt( ex*ex + ey*ey + ez*ez ) // eccentricity
	o.n = o.a/fabs(o.a)*sqrt(fabs(mu/(o.a*o.a*o.a))) // mean motion (negative if hyperbolic)
	o.P = 2*M_PI/o.n // period (negative if hyperbolic)

	o.inc = acos2(num: hz, denom: o.h, disambiguator: 1.0) // cosi = dot product of h and z unit vectors.  Always in [0,pi], so pass dummy disambiguator
							    // will = 0 if h is 0.0

	nx = -hy;// vector pointing along the ascending node = zhat cross h
	ny =  hx
	n = sqrt( nx*nx + ny*ny )

	// Omega, pomega and theta are measured from x axis, so we can always use y component to disambiguate if in the range [0,pi] or [pi,2pi]
	o.Omega = acos2(num: nx, denom: n, disambiguator: ny) // cos Omega is dot product of x and n unit vectors. Will = 0 if i=0.0

    if(o.e < 1.0){
	    ea = acos2(num: 1.0-o.d/o.a, denom: o.e, disambiguator: vr) // from definition of eccentric anomaly.  If vr < 0, must be going from apo to peri, so ea = [pi, 2pi] so ea = -acos(cosea)
	    o.M = ea - o.e*sin(ea)			// mean anomaly (Kepler's equation)
    }
    else{
        ea = acosh((1.0-o.d/o.a)/o.e)
        if(vr < 0.0){                    // Approaching pericenter, so eccentric anomaly < 0.0
            ea = -ea
        }
        o.M = o.e*sinh(ea) - ea          // Hyperbolic Kepler's equation
    }

	// in the near-planar case, the true longitude is always well defined for the position, and pomega for the pericenter if e!= 0
	// we therefore calculate those and calculate the remaining angles from them
	if(o.inc < MIN_INC || o.inc > M_PI - MIN_INC){	// nearly planar.  Use longitudes rather than angles referenced to node for numerical stability.
		o.theta = acos2(num: dx, denom: o.d, disambiguator: dy)		// cos theta is dot product of x and r vectors (true longitude). 
        o.pomega = acos2(num: ex, denom: o.e, disambiguator: ey)		// cos pomega is dot product of x and e unit vectors.  Will = 0 if e=0.0

		if(o.inc < M_PI/2.0){
			o.omega = o.pomega - o.Omega
			o.f = o.theta - o.pomega
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega + o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.0*o.e*sin(o.f) // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
		}
		else{
			o.omega = o.Omega - o.pomega
			o.f = o.pomega - o.theta
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega - o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.0*o.e*sin(o.f) // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
		}
	    
	}
	// in the non-planar case, we can't calculate the broken angles from vectors like above.  omega+f is always well defined, and omega if e!=0
	else{
		let wpf = acos2(num: nx*dx + ny*dy, denom: n*o.d, disambiguator: dz)	// omega plus f.  Both angles measured in orbital plane, and always well defined for i!=0.0
		o.omega = acos2(num: nx*ex + ny*ey, denom: n*o.e, disambiguator: ez)
		if(o.inc < M_PI/2.0){
			o.pomega = o.Omega + o.omega
			o.f = wpf - o.omega
			o.theta = o.Omega + wpf
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega + o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta - 2.0*o.e*sin(o.f) // M-f from Murray & Dermott Eq 2.93. This way l->theta smoothly as e->0
            }
		}
		else{
			o.pomega = o.Omega - o.omega
			o.f = wpf - o.omega
			o.theta = o.Omega - wpf
            if(o.e > MIN_ECC){              // pomega well defined
			    o.l = o.pomega - o.M
            }
            else{                           // when e << 1 and pomega ill defined, use l = theta+(M-f). M-f is O(e) so well behaved
                o.l = o.theta + 2.0*o.e*sin(o.f) // M-f from Murray & Dermott Eq 2.93 (retrograde changes sign). This way l->theta smoothly as e->0
            }
		}
	}
    
    if let sim = p.sim {
        o.T = sim.t - o.M/fabs(o.n) // time of pericenter passage (M = n(t-T).  Works for hyperbolic with fabs and n defined as above).
    }
    else {
        o.T = Double.nan // if particle isn't in simulation yet, can't get time.  You can still manually apply the equation below using o.M and o.n
    }

    // move some of the angles into [0,2pi) range
    o.f = tools_mod_twopi(f: o.f)
    o.l = tools_mod_twopi(f: o.l)
    o.M = tools_mod_twopi(f: o.M)
    o.theta = tools_mod_twopi(f: o.theta)
    o.omega = tools_mod_twopi(f: o.omega)

    return o
}
