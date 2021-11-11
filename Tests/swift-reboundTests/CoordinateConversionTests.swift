@testable import swift_rebound
import XCTest

class CoordinateConversionTests: XCTestCase {
    var particle : Particle!
    var orbit : Orbit!

    override func setUp() {
        super.setUp()
        particle = Particle()
        orbit = Orbit()
    }

    override func tearDown() {
        particle = nil
        orbit = nil
        super.tearDown()
    }

}