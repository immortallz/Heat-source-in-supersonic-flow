#include "r.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

FlowParams flowParams;
BodyParams bodyParams;
NumericalParams numericalParams;
HeatSource heatSource;


void test()
{
    // Пример использования базового класса и производных
    // G_array G = {1, 2, 3, 4, 5};

    // G_array F;
    // F = G;
    // F.print();

    // std::cout << G[0] << " " << G[4] << std::endl;
    // std::cout << G.get_rho(3) << std::endl;
    // std::cout << G.get_p(4) << std::endl;
    // std::cout << G.get_u() << std::endl;

    E_array E = {1,2,3,4,5};
    F_array F = {5,4,3,2,1};
    G_array G = {9,9,9,9,0};
    R_array R = {1,1,2,2,2};
    double rho = 1.2, p = 100000, u = 100, v = 200, w = 300, r = 3;
    rho = 1.615063, p = 149129.877770, u = 93.149968, v = 1.000000, w = 621.925338;

    E.data[0] = rho*u;
    E.data[1] = rho*u*u + p;
    E.data[2] = rho*u*v;
    E.data[3] = rho*u*w;
    E.data[4] = u * (Gamma * p / (Gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
    E = r * E;

    F.data[0] = rho*v;
    F.data[1] = rho*u*v;
    F.data[2] = rho*v*v + p;
    F.data[3] = rho*v*w;
    F.data[4] = v * (Gamma * p / (Gamma - 1) + rho*(u*u + v*v + w*w)*0.5);

    G.data[0] = rho*w;
    G.data[1] = rho*u*w;
    G.data[2] = rho*v*w;
    G.data[3] = rho*w*w + p;
    G.data[4] = w * (Gamma * p / (Gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
    G = r * G;

    R.data[0] = 0;
    R.data[1] = rho*v*v + p;
    R.data[2] = -rho*u*v;
    R.data[3] = 0;
    R.data[4] = 2.28;

    G.print();
    
    // G = 0.1*E + 1*G + 0.03*F;
    // G.print();
    // G = 1.0*G - 0.1*E - F*0.03;

    // (2*E).print();
    // (a*E + b*R).print();
    
    G.print();
    std::cout << G.get_rho(r) << std::endl;
    std::cout << G.get_p(r) << std::endl;
    std::cout << G.get_u() << std::endl;
    std::cout << G.get_v() << std::endl;
    std::cout << G.get_w() << std::endl;

    std::cout << "\n" << std::endl;

    E.print();
    get_E(G, r).print();
    F.print();
    get_F(G, r).print();
    R.print();
    get_R(G, r, 2.28).print();

    std::cout << "F_mirrored check:" << std::endl;
    F.print();
    F_array F_mirrored = (-1) * F;
    F_mirrored[2] = (-1) * F_mirrored[2];
    F_mirrored.print();
}

void force_and_momentum_by_heat_location(int N_z, int N_r, double L = 0.02, double Q = 0.2)
{
    std::ofstream
        Fy_out("Fy_out_heat.txt"),
        Mz_out("Mz_out_heat.txt"),
        Q_out("Q_out_heat.txt"),
        z_out("z_out_heat.txt"),
        r_out("r_out_heat.txt");
    std::vector<double> result;
    
    double
        h_z = (bodyParams.bodyLength - bodyParams.transitionPoint) / N_z, h_r,
        z = bodyParams.transitionPoint, r;

    for (int i = 0; i < N_z + 1; i++) {
        r = r_b(z);
        h_r = (tan(Pi / 12.0)*z - r_b(z)) / N_r;
        for (int j = 0; j < N_r + 1; j++) {
            HeatSource src {r, 0, z, L, Q};
            result = solver(src);
            z_out << z << " ";
            r_out << r << " ";
            Fy_out << result[0] << " ";
            Mz_out << result[1] << " ";
            Q_out << result[2] << " ";
            r += h_r;
        }
        z_out << "\n";
        r_out << "\n";
        Fy_out << "\n";
        Mz_out << "\n";
        Q_out << "\n";
        z += h_z;
    }
}

std::string to_string(BodyType type) {
    switch(type) {
        case BodyType::Cylindrical: return "Cylindrical";
        case BodyType::Cone: return "Cone";
        case BodyType::Parabolic: return "Parabolic";
        case BodyType::DoubleCone: return "DoubleCone";
    }
    return "Unknown";
}

void save_params(
    const FlowParams& flowParams,
    const BodyParams& bodyParams,
    const NumericalParams& numericalParams)
{
    json flowJson, bodyJson, numericalJson;
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

    std::ofstream
        flowFile("flowParams.json"),
        bodyFile("bodyParams.json"),
        numericalFile("numericalParams.json");

    flowFile << flowJson.dump(4); // красиво с отступами
    bodyFile << bodyJson.dump(4);
    numericalFile << numericalJson.dump(4);
}

int main()
{
    flowParams.Mach_inf = 3;
    flowParams.p_inf = 101330;
    flowParams.rho_inf = 1.2255;

    flowParams.a_inf = sqrt(Gamma * flowParams.p_inf / flowParams.rho_inf);
    flowParams.V_inf = flowParams.Mach_inf * flowParams.a_inf;

    flowParams.is_adiabatic = false;

    bodyParams.transitionPoint = 1.0; // not recommended to change, works bad on other values yet
    bodyParams.bodyLength = 2.0;
    bodyParams.bodyType = BodyType::DoubleCone;

    numericalParams.N = 100;
    numericalParams.M = 200;
    numericalParams.num_step_percent = 100;
    numericalParams.files_count = 100;

    heatSource.x = 0.5; // distance from body symmetry axis to heat source i.e. radial distance
    heatSource.y = 0.0; // do not change this because we assume source is located in theta=0 plane
    heatSource.z = 1.1; // distance from beginning of the body to heat source
    heatSource.L = 0.02;
    heatSource.Q = 1.0 / 3.0;

    // tan(pi / 12) ~ 0.26794919243

    // clock_t start = clock();
    // std::vector<double> res = solver(heatSource);
    // std::cout << "Fy = " << res[0] << std::endl;
    // std::cout << "Mz = " << res[1] << std::endl;
    // std::cout << "Q (total) = " << res[2] << std::endl;
    // clock_t finish = clock();
    // double seconds = double(finish - start) / CLOCKS_PER_SEC;
    // std::cout << "Elapsed time: " << seconds << "s" << std::endl;

    // save_params(flowParams, bodyParams, numericalParams);

    numericalParams.N = 200;
    numericalParams.M = 600;
    numericalParams.num_step_percent = 1;
    numericalParams.files_count = 1;

    force_and_momentum_by_heat_location(20, 15, 0.02, 0.2);
    force_and_momentum_by_heat_location(20, 15, 0.02, 0.05);

    return 0;
}