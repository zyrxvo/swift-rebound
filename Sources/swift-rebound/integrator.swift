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
    
    static func *(lhs: Vector, rhs: Vector) -> Double {
        return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z
    }
    
    static func *(scalar: Double, vec: Vector) -> Vector {
        return Vector(x: scalar*vec.x, y: scalar*vec.y, z: scalar*vec.z)
    }
    
    static func *(vec: Vector, scalar: Double) -> Vector {
        return Vector(x: scalar*vec.x, y: scalar*vec.y, z: scalar*vec.z)
    }
}

class Integrator {
    init(){}

    func calculateForcePair(G: Double, p1: Particle, p2: Particle) -> Vector {
        let r = Vector(vec: [p1.x - p2.x, p1.y - p2.y, p1.z - p2.z])
        let rmag = sqrt( r*r )
        let factor = -(G * p1.m * p2.m)/(rmag * rmag * rmag)
        return factor*r
    }

    func calculateAcceleration(G: Double, particles: [Particle]) {
        // Preallocate an NxNx3 array.
//        var forces = Array(repeating: Array(repeating: Array(repeating: 0.0, count: 3), count: particles.count), count: particles.count)
        for i in 0..<particles.count {
            var total_force = Vector(vec: [0.0, 0.0, 0.0])
            for j in 0..<particles.count {
                if i != j {
                    let force_pair = calculateForcePair(G: G, p1: particles[i], p2: particles[j])
                    total_force = total_force + force_pair
                }
            }
            particles[i].ax = total_force.x
            particles[i].ay = total_force.y
            particles[i].az = total_force.z
        }
    }
    func step() {}
}

class LeapFrog: Integrator {
    override func step() {
        drift()
        kick()
        drift()
    }
    func drift() {

    }

    func kick() {

    }
}
