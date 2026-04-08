#include "r.h"

double r_b(const double z) {
    switch (bodyParams.bodyType)
    {
        case BodyType::Cylindrical:
            return tan(PI / 12.0);

        case BodyType::Cone:
            return tan(PI / 12.0) * z;

        case BodyType::Parabolic:
            return tan(PI / 12.0) * sqrt(2*z - 1);

        case BodyType::DoubleCone: {
            constexpr double SLOPE = 0.01;      // coefficient of cone-to-cylinder approximation
            return tan(PI / 12.0) + SLOPE * (z - 1.0);
        }

        default:
            return tan(PI / 12.0) * z;
    }
}

double r_b_z(const double z) {
    double dz = 1e-12;
    return (r_b(z + dz) - r_b(z - dz)) / dz * 0.5;
}
