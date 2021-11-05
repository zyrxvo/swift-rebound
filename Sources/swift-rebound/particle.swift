
class Particle {
    var m, x, y, z, vx, vy, vz, a, P, e, inc, Omega, omega, pomega, f, M, E, l: Double
    
    init(m: Double = 0.0, x: Double = 0.0, y: Double = 0.0, z: Double = 0.0, vx: Double = 0.0, vy: Double = 0.0, vz: Double = 0.0, a: Double = 0.0, P: Double = 0.0, e: Double = 0.0, inc: Double = 0.0, Omega: Double = 0.0, omega: Double = 0.0, pomega: Double = 0.0, f: Double = 0.0, M: Double = 0.0, E: Double = 0.0, l: Double = 0.0){
        self.m = m; self.x = x; self.y = y; self.z = z; self.vx = vx; self.vy = vy; self.vz = vz; self.a = a; self.P = P; self.e = e; self.inc = inc; self.Omega = Omega; self.omega = omega; self.pomega = pomega; self.f = f; self.M = M; self.E = E; self.l = l
    }

}

class ParticleNAN: Particle {
    init(){
        super.init()
        super.m = Double.nan; super.x = Double.nan; super.y = Double.nan; super.z = Double.nan; super.vx = Double.nan; super.vy = Double.nan; super.vz = Double.nan; super.a = Double.nan; super.P = Double.nan; super.e = Double.nan; super.inc = Double.nan; super.Omega = Double.nan; super.omega = Double.nan; super.pomega = Double.nan; super.f = Double.nan; super.M = Double.nan; super.E = Double.nan; super.l = Double.nan
    }
}

extension Particle: CustomStringConvertible {
    var description: String {
        // Make sure Cartesian coordinates are defined.
        return "particle object < m=\(m), x=\(x), y=\(y), z=\(z), vx=\(vx), vy=\(vy), vz=\(vz) >"
    }
    
}