class Simulation {
    var t : Double = 0.0
    var G : Double = 1.0
    var dt : Double = 0.001
    var dt_last_done : Int = 0
    var steps_done : Int = 0
    var N : Int {
        get {
            return self.particles.count
        }
    }
    var particles : [Particle] = []
    
    init(){
    }
    
    public func add(particle: Particle? = nil, m: Double? = nil, x: Double? = nil, y: Double? = nil, z: Double? = nil, vx: Double? = nil, vy: Double? = nil, vz: Double? = nil, ax: Double? = nil, ay: Double? = nil, az: Double? = nil, a: Double? = nil, P: Double? = nil, e: Double? = nil, inc: Double? = nil, Omega: Double? = nil, omega: Double? = nil, pomega: Double? = nil, f: Double? = nil, M: Double? = nil, E: Double? = nil, l: Double? = nil, theta: Double? = nil, T: Double? = nil, primary: Particle? = nil, G: Double? = 1.0) {
        if let p = particle {
            self.particles.append(p)
        }
        else {
            do {
                var prim : Particle?
                if self.N > 0 {
                    // If there are already particles in the Simulation determine the primary particle. Either the primary particle was passed in, or the default primary (the first particle added) is used.
                    if let p = primary { prim = p }
                    else { prim = self.particles[0] }
                }
                let p = try Particle(m: m, x: x, y: y, z: z, vx: vx, vy: vy, vz: vz, ax: ax, ay: ay, az: az, a: a, P: P, e: e, inc: inc, Omega: Omega, omega: omega, pomega: pomega, f: f, M: M, E: E, l: l, theta: theta, T: T, primary: prim, G: G)
                self.particles.append(p)
            }
            catch {
                print(error)
            }
        }
    }
}

extension Simulation: CustomStringConvertible {
    var description: String {
        // Make sure Cartesian coordinates are defined.
        var s = ""
        s += "---------------------------------\n"
        //        s += "REBOUND version:     \t%s\n" %__version__
        //        s += "REBOUND built on:    \t%s\n" %__build__
        s += "Number of particles: \t\(self.N)\n"
        //        s += "Selected integrator: \t" + self.integrator + "\n"
        s += "Simulation time:     \t\(String(format: "%.16f", self.t))\n"
        s += "Current timestep:    \t\(String(format: "%f", self.dt))\n"
        if self.N > 0 {
            s += "---------------------------------\n"
            for p in self.particles {
                s += "< m=\(p.m), x=\(p.x), y=\(p.y), z=\(p.z), vx=\(p.vx), vy=\(p.vy), vz=\(p.vz) >\n"
            }
        }
        s += "---------------------------------"
        return s
    }
    
}
