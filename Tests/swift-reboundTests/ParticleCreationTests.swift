@testable import swift_rebound
import XCTest

class ParticleCreationTests: XCTestCase {
    var primary : Particle!
    var orbit : Orbit!
    let G = 1.0

    override func setUp() {
        super.setUp()
        do {
            primary = try Particle(m: 1)
        }
        catch {
            print(error)
        }
       
        orbit = Orbit()
    }

    override func tearDown() {
        primary = nil
        orbit = nil
        super.tearDown()
    }

}
