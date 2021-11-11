
class Particle {
    var m, x, y, z, vx, vy, vz, ax, ay, az, a, P, e, inc, Omega, omega, pomega, f, M, E, l: Double
    var sim : Simulation?
    
    init(m: Double = 0.0, x: Double = 0.0, y: Double = 0.0, z: Double = 0.0, vx: Double = 0.0, vy: Double = 0.0, vz: Double = 0.0, ax: Double = 0.0, ay: Double = 0.0, az: Double = 0.0, a: Double = 0.0, P: Double = 0.0, e: Double = 0.0, inc: Double = 0.0, Omega: Double = 0.0, omega: Double = 0.0, pomega: Double = 0.0, f: Double = 0.0, M: Double = 0.0, E: Double = 0.0, l: Double = 0.0){
        self.m = m; self.x = x; self.y = y; self.z = z; self.vx = vx; self.vy = vy; self.vz = vz; self.ax = ax; self.ay = ay; self.az = az; self.a = a; self.P = P; self.e = e; self.inc = inc; self.Omega = Omega; self.omega = omega; self.pomega = pomega; self.f = f; self.M = M; self.E = E; self.l = l
    }

}
extension Particle: CustomStringConvertible {
    var description: String {
        // Make sure Cartesian coordinates are defined.
        return "< m=\(m), x=\(x), y=\(y), z=\(z), vx=\(vx), vy=\(vy), vz=\(vz) >"
    }
    
}

enum OrbitError : Error {
    case radialOrbit
    case negativeEccentricity
    case boundOrbitError
    case unboundOrbitError
    case invalidTrueAnomaly
    case masslessPrimary
    case coincidesWithPrimary
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
        case .coincidesWithPrimary:
            return "Particle is on top of primary."
        }
    }
    
}

class Orbit {
    /// """
    /// A class containing orbital parameters for a particle.

    /// When using the various swift-rebound functions using Orbits, all angles are in radians. 
    /// In swift-rebound the reference direction is the positive x direction, the reference plane is the xy plane.

    /// Attributes
    /// ----------
    /// d       : float           
    ///     radial distance from reference 
    /// v       : float         
    ///     velocity relative to central object's velocity
    /// h       : float           
    ///     specific angular momentum
    /// P       : float           
    ///     orbital period (negative if hyperbolic)
    /// n       : float           
    ///     mean motion    (negative if hyperbolic)
    /// a       : float           
    ///     semimajor axis
    /// e       : float           
    ///     eccentricity
    /// inc     : float           
    ///     inclination
    /// Omega   : float           
    ///     longitude of ascending node
    /// omega   : float           
    ///     argument of pericenter
    /// pomega  : float           
    ///     longitude of pericenter
    /// f       : float           
    ///     true anomaly
    /// M       : float           
    ///     mean anomaly
    /// E       : float           
    ///     eccentric anomaly (requires solving Kepler's equation - only calculated when needed)
    /// l       : float           
    ///     mean longitude = Omega + omega + M
    /// theta   : float           
    ///     true longitude = Omega + omega + f
    /// T       : float
    ///     time of pericenter passage
    /// rhill   : float
    ///     Hill radius ( =a*pow(m/(3M),1./3.) )
    /// """  
    var d, v, h, P, n, a, e, inc, Omega, omega, pomega, f, M, l, theta, T, rhill: Double
    
    init(d: Double = 0.0, v: Double = 0.0, h: Double = 0.0, P: Double = 0.0, n: Double = 0.0, a: Double = 0.0, e: Double = 0.0, inc: Double = 0.0, Omega: Double = 0.0, omega: Double = 0.0, pomega: Double = 0.0, f: Double = 0.0, M: Double = 0.0, l: Double = 0.0, theta: Double = 0.0, T: Double = 0.0, rhill: Double = 0.0){
        self.d = d; self.v = v; self.h = h; self.P = P; self.n = n; self.a = a; self.e = e; self.inc = inc; self.Omega = Omega; self.omega = omega; self.pomega = pomega; self.f = f; self.M = M; self.l = l; self.theta = theta; self.T = T; self.rhill = rhill;
    }

}

class OrbitNAN: Orbit {  
    init(){
        super.init()
        super.d = Double.nan; super.v = Double.nan; super.h = Double.nan; super.P = Double.nan; super.n = Double.nan; super.a = Double.nan; super.e = Double.nan; super.inc = Double.nan; super.Omega = Double.nan; super.omega = Double.nan; super.pomega = Double.nan; super.f = Double.nan; super.M = Double.nan; super.l = Double.nan; super.theta = Double.nan; super.T = Double.nan; super.rhill = Double.nan;
    }

}

extension Orbit: CustomStringConvertible {
    var description: String {
        // Make sure Cartesian coordinates are defined.
        return "< a=\(a), e=\(e), inc=\(inc), Omega=\(Omega), omega=\(omega), f=\(f)>"
    }
    
}
