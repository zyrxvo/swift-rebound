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

    func test_particle_initialization_with_cartesian_and_orbital_elements() throws {
        let expectedError = OrbitError.mixedCartesianAndKepler
        var err: OrbitError?
        var cart : [Double?]
        var orbi : [Double?]
        // Test every combination of cartesian coordiantes and orbital elements using only one from each.
        for i in 0...5 {
            for j in 0...12 {
                cart = [nil, nil, nil, nil, nil, nil]
                orbi = [nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil]
                cart[i] = 0.0
                orbi[j] = 0.0
                XCTAssertThrowsError(try Particle(x: cart[0], y: cart[1], z: cart[2], vx: cart[3], vy: cart[4], vz: cart[5],
                                                  a: orbi[0], P: orbi[1], e: orbi[2], inc: orbi[3], Omega: orbi[4], omega: orbi[5],
                                                  pomega: orbi[6], f: orbi[7], M: orbi[8], E: orbi[9], l: orbi[10], theta: orbi[11], T: orbi[12])) {
                    thrownError in err = thrownError as? OrbitError
                }
                XCTAssertEqual(expectedError, err)
                XCTAssertEqual(expectedError.description, err?.description)
            }
        }
    }
}
