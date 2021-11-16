import Foundation

struct Vector {
    var x : Double = 0.0
    var y : Double = 0.0
    var z : Double = 0.0
    
    init(vec: [Double]){
        self.x = vec[0]
        self.y = vec[1]
        self.z = vec[2]
    }
    init(x: Double, y: Double, z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
    
    static func +(lhs: Vector, rhs: Vector) -> Vector {
        return Vector(x: lhs.x + rhs.x, y: lhs.y + rhs.y, z: lhs.z + rhs.z)
    }
    
    static func +=(lhs : inout Vector, rhs: Vector) {
        lhs.x = lhs.x + rhs.x
        lhs.y = lhs.y + rhs.y
        lhs.z = lhs.z + rhs.z
    }
    
    static func -(lhs: Vector, rhs: Vector) -> Vector {
        return Vector(x: lhs.x - rhs.x, y: lhs.y - rhs.y, z: lhs.z - rhs.z)
    }
    
    static func -=(lhs : inout Vector, rhs: Vector) {
        lhs.x = lhs.x - rhs.x
        lhs.y = lhs.y - rhs.y
        lhs.z = lhs.z - rhs.z
    }
    
    static func *(lhs: Vector, rhs: Vector) -> Double {
        return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z
    }
    
    static func *(scalar: Double, vec: Vector) -> Vector {
        return Vector(x: scalar*vec.x, y: scalar*vec.y, z: scalar*vec.z)
    }
    
    static func *(vec: Vector, scalar: Double) -> Vector {
        return Vector(x: scalar*vec.x, y: scalar*vec.y, z: scalar*vec.z)
    }
    
    static func *=(lhs : inout Vector, scalar: Double) {
        lhs.x = lhs.x * scalar
        lhs.y = lhs.y * scalar
        lhs.z = lhs.z * scalar
    }
}

class Integrator {
    var name : String {
        get { return "none" }
    }
    init(){}

    func calculateAccelerationPair(G: Double, p1: Particle, p2: Particle) -> Vector {
        let r = Vector(vec: [p1.x - p2.x, p1.y - p2.y, p1.z - p2.z])
        let rmag = sqrt( r*r )
        // Don't include p1.m because this is computing only the acceleration.
        let factor = -(G * p2.m)/(rmag * rmag * rmag)
        return factor*r
    }

    func calculateAcceleration(G: Double, particles: [Particle]) {
        for i in 0..<particles.count {
            var total_accel = Vector(vec: [0.0, 0.0, 0.0])
            for j in 0..<particles.count {
                if i != j {
                    let accel_pair = calculateAccelerationPair(G: G, p1: particles[i], p2: particles[j])
                    total_accel += accel_pair
                }
            }
            particles[i].ax = total_accel.x
            particles[i].ay = total_accel.y
            particles[i].az = total_accel.z
        }
    }
    func step(s: Simulation) {}
}
extension Integrator: CustomStringConvertible {
    public var description: String {
        return self.name
    }
    
}

class LeapFrog: Integrator {
    override var name: String {
        get { return "Leap Frog" }
    }
    
    override func step(s: Simulation) {
        drift(s: s)
        kick(s: s)
        drift(s: s)
        s.t += s.dt
        s.steps_done += s.steps_done
    }
    
    func drift(s: Simulation) {
        for p in s.particles {
            p.x += p.vx * 0.5 * s.dt
            p.y += p.vy * 0.5 * s.dt
            p.z += p.vz * 0.5 * s.dt
        }
    }

    func kick(s: Simulation) {
        self.calculateAcceleration(G: s.G, particles: s.particles)
        for p in s.particles {
            p.vx += p.ax * s.dt
            p.vy += p.ay * s.dt
            p.vz += p.az * s.dt
        }
    }
}
