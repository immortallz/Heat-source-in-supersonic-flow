#include "r.h"

double r_b(double z, BodyType bodyType) {
    switch (bodyType)
    {
        case BodyType::Cylindrical:
            return tan(Pi / 12.0);

        case BodyType::Cone:
            return tan(Pi / 12.0) * z;

        case BodyType::Parabolic:
            return tan(Pi / 12.0) * sqrt(2*z - 1);

        case BodyType::DoubleCone: {
            double k = 0.01;
            return tan(Pi / 12.0) + k * (z - 1.0);
        }

        default:
            return tan(Pi / 12.0) * z;
    }
}

double r_b_z(double z, BodyType bodyType) {
    double dz = 1e-12;
    return (r_b(z + dz, bodyType) - r_b(z - dz, bodyType)) / dz * 0.5;
}
