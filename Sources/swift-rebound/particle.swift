import Foundation

class Particle: Equatable  {
    var m, x, y, z, vx, vy, vz, ax, ay, az, a, P, e, inc, Omega, omega, pomega, f, M, E, l, theta, T: Double
    var primary : Particle?
    var sim : Simulation?
    
    init(m: Double? = nil, x: Double? = nil, y: Double? = nil, z: Double? = nil, vx: Double? = nil, vy: Double? = nil, vz: Double? = nil, ax: Double? = nil, ay: Double? = nil, az: Double? = nil, a: Double? = nil, P: Double? = nil, e: Double? = nil, inc: Double? = nil, Omega: Double? = nil, omega: Double? = nil, pomega: Double? = nil, f: Double? = nil, M: Double? = nil, E: Double? = nil, l: Double? = nil, theta: Double? = nil, T: Double? = nil, primary: Particle? = nil, G: Double? = 1.0) throws {
        
        self.m = m ?? 0.0 // If no mass is given, set the mass to 0.
        
        let cart = [x, y, z, vx, vy, vz]
        let orbi = [a,P,e,inc,Omega,omega,pomega,f,M,E,l,theta,T]
        
        // Check if user has mixed cartesian coordinates with orbital elements.
        var isCart = false
        for el in cart {
            isCart = isCart || el != nil
        }
        var isOrbi = false
        for el in orbi {
            isOrbi = isOrbi || el != nil
        }
        // Only throws an error if isCart is true AND isOrbit is true.
        guard !isCart || !isOrbi else { throw OrbitError.mixedCartesianAndKepler }
        
        if isOrbi {
            print("Orbital elements were passed given.")
            if primary == nil {
                throw OrbitError.missingPrimary
            }
            else {
                self.primary = primary
            }
            
            guard a != nil || P != nil else { throw OrbitError.missing_a_or_P }
            guard a == nil || P == nil else { throw OrbitError.mixed_a_and_P }
            
            if a == nil {
                self.P = P ?? 0.0
                let P2 = pow(self.P,2)
                let mu = (G ?? 1.0) * (self.primary!.m + self.m)
                let twopi2 = pow(2.0*PI, 2)
                let base = P2 * mu/twopi2
                let exponent = 1.0/3.0
                self.a = pow(base, exponent)
            }
            else {
                self.a = a ?? 0.0
                let mu = (G ?? 1.0) * (self.primary!.m + self.m)
                self.P = (2.0*PI) * sqrt( pow(self.a, 3)/mu )
            }
            // Normal orbital parameters.
            self.e = e ?? 0.0
            self.inc = inc ?? 0.0
            self.Omega = Omega ?? 0.0 // We require that Omega be passed if you want to specify longitude of node.
            
            // Need omega for conversion function. It can be specified either directly or through pomega indirectly.
            guard omega == nil || pomega == nil else { throw OrbitError.mixed_omega_and_pomega }
            if omega == nil && pomega == nil {
                // If neither omega or pomega were passed defualt omega to 0.
                self.omega = 0.0
            }
            // Use omega if it was passed. If pomega was passed, use it to update omega.
            self.omega = omega ?? 0.0
            self.pomega = pomega ?? 0.0
            if pomega != nil {
                if cos(self.inc) > 0 {
                    // inc is in range [-pi/2, pi/2] (prograde), so pomega = Omega + omega
                    self.omega = self.pomega - self.Omega
                }
                else {
                    // For retrograde orbits, pomega = Omega - omega.
                    self.omega = self.Omega - self.pomega
                }
            }
                        
            // Can specify longitude through any of these, but only accept one. Need f for conversion function.
            let longitudes = [f,M,E,l,theta,T]
            var nils = 0
            for angle in longitudes {
                nils = (angle == nil) ? nils + 1 : nils
            }
            // Check that only one phase angle was given.
            guard nils >= 5 else { throw OrbitError.mixedAngles }
            
            // Use f if it was passed. If a different angle was passed update f.
            self.f = f ?? 0.0
            self.M = M ?? 0.0
            self.E = E ?? 0.0
            self.l = l ?? 0.0
            self.theta = theta ?? 0.0
            self.T = T ?? 0.0
            
            // If it was passed, use theta to compute f because it's the next easiest.
            if theta != nil {
                if cos(self.inc) > 0 {
                    // For prograde orbits, theta = Omega + omega + f.
                    self.f = self.theta - self.Omega - self.omega
                }
                else {
                    // For retrograde, theta = Omega - omega - f.
                    self.f = self.Omega - self.omega - self.theta
                }
            }
            // Next, if passed, use l to compute f.
            if l != nil {
                if cos(self.inc) > 0 {
                    // For prograde orbits, l = Omega + omega + M.
                    self.M = self.l - self.Omega - self.omega
                }
                else {
                    // For retrograde, l = Omega - omega - M.
                    self.M = self.Omega - self.omega - self.l
                }
                self.f = tools_M_to_f(e: self.e, M: self.M)
            }
            // Next, if passed, use T to compute f.
            if T != nil {
                // Works for both elliptical and hyperbolic orbits
                // TODO: has accuracy problems for M=n*(t-T) << 1
                let n = pow(self.sim!.G * (primary!.m + self.m)/abs( pow(self.a, 3.0) ), 0.5)
                self.M = n*(self.sim!.t - self.T)
                self.f = tools_M_to_f(e: self.e, M: self.M)
            }
            // Next, if passed, use M to compute f.
            if M != nil {
                self.f = tools_M_to_f(e: self.e, M: self.M)
            }
            // Finally, if passed, use E to compute f.
            if E != nil {
                self.f = tools_E_to_f(e: self.e, E: self.E)
            }
            
            // With all of the necessary orbital elements computed, compute the cartesian coordinates.
            let p = try tools_orbital_elements_to_particle(G: (G ?? 1.0), primary: self.primary!, m: self.m, a: self.a, e: self.e, inc: self.inc, Omega: self.Omega, omega: self.omega, f: self.f)
            self.x = p.x
            self.y = p.y
            self.z = p.z
            self.vx = p.vx
            self.vy = p.vy
            self.vz = p.vz
            
        }
        else {
            // If cartesian coordinates were passed assign the appropriate values.
            self.x = x ?? 0.0
            self.y = y ?? 0.0
            self.z = z ?? 0.0
            self.vx = vx ?? 0.0
            self.vy = vy ?? 0.0
            self.vz = vz ?? 0.0
            
            // Assign the remaining property values to 0. The correct values will be computed when needed.
            self.a = a ?? 0.0
            self.P = P ?? 0.0
            self.e = e ?? 0.0
            self.inc = inc ?? 0.0
            self.Omega = Omega ?? 0.0
            self.omega = omega ?? 0.0
            self.pomega = pomega ?? 0.0
            self.f = f ?? 0.0
            self.M = M ?? 0.0
            self.E = E ?? 0.0
            self.l = l ?? 0.0
            self.theta = theta ?? 0.0
            self.T = T ?? 0.0
        }
        self.ax = ax ?? 0.0
        self.ay = ay ?? 0.0
        self.az = az ?? 0.0
    }
    
    static func ==(lhs: Particle, rhs: Particle) -> Bool {
        return lhs.m == rhs.m && lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.vx == rhs.vx && lhs.vy == rhs.vy && lhs.vz == rhs.vx && lhs.ax == rhs.ax && lhs.ay == rhs.ay && lhs.az == rhs.az && lhs.a == rhs.a && lhs.P == rhs.P && lhs.e == rhs.e && lhs.inc == rhs.inc && lhs.Omega == rhs.Omega && lhs.omega == rhs.omega && lhs.pomega == rhs.pomega && lhs.f == rhs.f && lhs.M == rhs.M && lhs.E == rhs.E && lhs.l == rhs.l
    }

}
extension Particle: CustomStringConvertible {
    var description: String {
        // Make sure Cartesian coordinates are defined.
        return "< m=\(m), x=\(x), y=\(y), z=\(z), vx=\(vx), vy=\(vy), vz=\(vz) >"
    }
    
}

enum OrbitError: Error {
    case radialOrbit
    case negativeEccentricity
    case boundOrbitError
    case unboundOrbitError
    case invalidTrueAnomaly
    case masslessPrimary
    case missingPrimary
    case coincidesWithPrimary
    case mixedCartesianAndKepler
    case missing_a_or_P
    case mixed_a_and_P
    case mixed_omega_and_pomega
    case mixedAngles
}

extension OrbitError: CustomStringConvertible {
    var description: String {
        switch self {
        case .radialOrbit:
            return "Can't initialize a radial orbit with orbital elements."
        case .negativeEccentricity:
            return "Eccentricity must be greater than or equal to zero."
        case .boundOrbitError:
            return "Bound orbit (a > 0) must have e < 1."
        case .unboundOrbitError:
            return "Unbound orbit (a < 0) must have e > 1."
        case .invalidTrueAnomaly:
            return "Unbound orbit can't have true anomaly, f, set beyond the range allowed by the asymptotes set by the parabola."
        case .masslessPrimary:
            return "Primary has no mass."
        case .missingPrimary:
            return "Need to specify a primary when initializing particle with orbital elements."
        case .coincidesWithPrimary:
            return "Particle is on top of primary."
        case .mixedCartesianAndKepler:
            return "You cannot pass cartesian coordinates and orbital elements (and/or primary) at the same time."
        case .missing_a_or_P:
            return "You need to pass either a semimajor axis, a, or orbital period, P, to initialize the particle using orbital elements."
        case .mixed_a_and_P:
            return "You can pass either the semimajor axis, a, or orbital period, P, but not both."
        case .mixed_omega_and_pomega:
            return "Can't pass both omega and pomega."
        case .mixedAngles:
            return "Can only pass one longitude/anomaly in the set [f, M, E, l, theta, T]."
        }
    }
    
}


/// A class containing orbital parameters for a particle.
///
/// When using the various swift-rebound functions using Orbits, all angles are in radians.
/// In swift-rebound the reference direction is the positive x direction, the reference plane is the xy plane.
class Orbit {
    var d, v, h, P, n, a, e, inc, Omega, omega, pomega, f, M, l, theta, T, rhill: Double
    
    /// Orbit initializer
    /// - Parameters:
    ///   - d: radial distance from reference
    ///   - v: velocity relative to central object's velocity
    ///   - h: specific angular momentum
    ///   - P: orbital period (negative if hyperbolic)
    ///   - n: mean motion    (negative if hyperbolic)
    ///   - a: semimajor axis
    ///   - e: eccentricity
    ///   - inc: inclination
    ///   - Omega: longitude of ascending node
    ///   - omega: argument of pericenter
    ///   - pomega: longitude of pericenter
    ///   - f: true anomaly
    ///   - M: mean anomaly
    ///   - l: mean longitude = Omega + omega + M
    ///   - theta: true longitude = Omega + omega + f
    ///   - T: time of pericenter passage
    ///   - rhill: Hill radius ( =a*pow(m/(3M),1./3.) )
    init(d: Double = 0.0, v: Double = 0.0, h: Double = 0.0, P: Double = 0.0, n: Double = 0.0, a: Double = 0.0, e: Double = 0.0, inc: Double = 0.0, Omega: Double = 0.0, omega: Double = 0.0, pomega: Double = 0.0, f: Double = 0.0, M: Double = 0.0, l: Double = 0.0, theta: Double = 0.0, T: Double = 0.0, rhill: Double = 0.0){
        self.d = d; self.v = v; self.h = h; self.P = P; self.n = n; self.a = a; self.e = e; self.inc = inc; self.Omega = Omega; self.omega = omega; self.pomega = pomega; self.f = f; self.M = M; self.l = l; self.theta = theta; self.T = T; self.rhill = rhill;
    }

}

extension Orbit: CustomStringConvertible {
    var description: String {
        // Make sure Orbital elements are defined.
        return "< a=\(a), e=\(e), inc=\(inc), Omega=\(Omega), omega=\(omega), f=\(f)>"
    }
    
}
