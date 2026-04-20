#include "r.h"
#include <iostream>
#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

FlowParams flowParams;
BodyParams bodyParams;
NumericalParams numericalParams;
HeatSource heatSource;

std::string get_body_type_path() {
    switch(bodyParams.bodyType) {
        case BodyType::Cylindrical: return "cylinder";
        case BodyType::Cone: return "cone";
        case BodyType::Parabolic: return "parabolic";
        case BodyType::DoubleCone: return "double_cone";
    }
    return "unknown";
}

void force_and_momentum_by_heat_location(int N_z, int N_r, double L = 0.02, double Q = 1.0/20.0)
{
    numericalParams.log_force_and_momentum_only = true;
    std::string basePath =
        "heat_source_variation/"
        + get_body_type_path()
        + "_L_" + std::to_string(bodyParams.bodyLength).substr(0, 3)
        + "/Q_" + std::to_string(Q).substr(0, 4) + "/";
    std::ofstream
        Fy_out(basePath + "Fy_out_heat.txt"),
        Mz_out(basePath + "Mz_out_heat.txt"),
        Q_out(basePath + "Q_out_heat.txt"),
        z_out(basePath + "z_out_heat.txt"),
        r_out(basePath + "r_out_heat.txt");
    std::vector<double> result;
    
    double
        h_z = (bodyParams.bodyLength - bodyParams.transitionPoint) / N_z, h_r,
        z = bodyParams.transitionPoint, r;

    for (int i = 0; i < N_z + 1; i++) {
        r = r_b(z);
        h_r = (tan(PI / 12.0)*2.0 - r_b(z)) / N_r; // hardcoded upper limit = tan(PI / 12.0)*2.0
        for (int j = 0; j < N_r + 1; j++) {
            HeatSource src {r, 0, z, L, Q};
            result = solver(src);
            z_out << z << " ";
            r_out << r << " ";
            Fy_out << result[0] << " ";
            Mz_out << result[1] << " ";
            Q_out << result[2] << " ";
            r += h_r;
            if (!std::isnan(result[0])) {
                std::cout << "solver completed successfully for ";
            } else {
                std::cout << "solver returned NaN for ";
            }
            std::cout << get_body_type_path()
                + ", Q = " + std::to_string(Q).substr(0, 4)
                + ", (i, j) = (" + std::to_string(i) + ", " + std::to_string(j) + ")"
                << std::endl;
        }
        z_out << "\n";
        r_out << "\n";
        Fy_out << "\n";
        Mz_out << "\n";
        Q_out << "\n";
        z += h_z;
    }
}

std::string to_string(const BodyType type) {
    switch(type) {
        case BodyType::Cylindrical: return "Cylindrical";
        case BodyType::Cone: return "Cone";
        case BodyType::Parabolic: return "Parabolic";
        case BodyType::DoubleCone: return "DoubleCone";
    }
    return "Unknown";
}

void save_params()
{
    json flowJson, bodyJson, numericalJson, heatSourceJson;
    flowJson["Mach_inf"] = flowParams.Mach_inf;
    flowJson["p_inf"] = flowParams.p_inf;
    flowJson["rho_inf"] = flowParams.rho_inf;

    bodyJson["transitionPoint"] = bodyParams.transitionPoint;
    bodyJson["bodyLength"] = bodyParams.bodyLength;
    bodyJson["bodyType"] = to_string(bodyParams.bodyType);

    numericalJson["N"] = numericalParams.N;
    numericalJson["M"] = numericalParams.M;
    numericalJson["num_step_percent"] = numericalParams.num_step_percent;
    numericalJson["files_count"] = numericalParams.files_count;
    numericalJson["flux_scheme"] = numericalParams.flux_scheme;

    heatSourceJson["x"] = heatSource.x;
    heatSourceJson["y"] = heatSource.y;
    heatSourceJson["z"] = heatSource.z;
    heatSourceJson["L"] = heatSource.L;
    heatSourceJson["Q"] = heatSource.Q;

    std::ofstream
        flowFile("parameters/flowParams.json"),
        bodyFile("parameters/bodyParams.json"),
        numericalFile("parameters/numericalParams.json"),
        heatSourceFile("parameters/heatSourceParams.json");

    flowFile << flowJson.dump(4); // красиво с отступами
    bodyFile << bodyJson.dump(4);
    numericalFile << numericalJson.dump(4);
    heatSourceFile << heatSourceJson.dump(4);
}

int main()
{
    flowParams.Mach_inf = sqrt(2);
    flowParams.p_inf = 101330;
    flowParams.rho_inf = 1.2255;
    flowParams.a_inf = sqrt(GAMMA * flowParams.p_inf / flowParams.rho_inf);
    flowParams.V_inf = flowParams.Mach_inf * flowParams.a_inf;
    flowParams.is_adiabatic = false;

    bodyParams.transitionPoint = 1.0; // not recommended to change, works bad on other values yet
    bodyParams.bodyLength = bodyParams.transitionPoint + 5.0;
    bodyParams.bodyType = BodyType::Cylindrical;

    constexpr int numerical_mesh_multiplier = 16;

    numericalParams.N = 100 * numerical_mesh_multiplier;
    numericalParams.M = 300 * numerical_mesh_multiplier;
    numericalParams.CFL = 0.9;
    numericalParams.num_step_percent = 1000;
    numericalParams.files_count = 10000;
    numericalParams.flux_scheme = FluxScheme::MacCormack;
    numericalParams.log_force_and_momentum_only = false;

    double r_b_value = r_b(bodyParams.transitionPoint);
    std::cout << "r_b_value = " << r_b_value << std::endl;
    heatSource.x = 1.2 * r_b_value; // distance from body symmetry axis to heat source i.e. radial distance
    heatSource.L = 0.01 * r_b_value;
    heatSource.y = 0.0; // do not change this because we assume source is located in theta=0 plane
    heatSource.z = 4.3; // distance from beginning of the body to heat source
    heatSource.Q = 1.0 / 20.0;

    // tan(pi / 12) ~ 0.26794919243

    const double start = omp_get_wtime();

    // const std::vector<double> res = solver(heatSource);
    //
    // std::cout << "Lifting force: " << res[0] << std::endl;
    // std::cout << "Rotation momentum: " << res[1] << std::endl;
    // std::cout << "Q = " << res[2] << std::endl;

    numericalParams.N = 200;
    numericalParams.M = numericalParams.N * 3;
    numericalParams.num_step_percent = 1;
    numericalParams.files_count = 1;

    bodyParams.bodyType = BodyType::Cylindrical;

    // force_and_momentum_by_heat_location(30, 30, 0.02, 0.25);
    force_and_momentum_by_heat_location(1, 1, 0.02, 1.0 / 20.0);

    // bodyParams.bodyType = BodyType::Parabolic;
    //
    // force_and_momentum_by_heat_location(30, 30, 0.02, 0.05);

    save_params();

    const double finish = omp_get_wtime();
    const double seconds = finish - start;
    std::cout << "Elapsed time: " << seconds << "s" << std::endl;

    return 0;
}
