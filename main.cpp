#include "r.h"

Params params;

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

void force_and_momentum_by_heat_location()
{
    std::vector <double>
        z_q_array = {
            1.     , 1.01017, 1.02045, 1.03006, 1.04059, 1.05043, 1.06038,
            1.07045, 1.08064, 1.09009, 1.10051, 1.11073, 1.12067, 1.13084,
            1.14065, 1.15077, 1.16008, 1.17043, 1.1809 , 1.19052, 1.20025,
            1.21009, 1.22003, 1.23008, 1.24024, 1.25052, 1.26091, 1.27036,
            1.28097, 1.29062, 1.30038, 1.31023, 1.32018, 1.33024, 1.34041,
            1.35068, 1.36105, 1.37037, 1.38096, 1.39046, 1.40005, 1.41096,
            1.42075, 1.43063, 1.44061, 1.45068, 1.46085, 1.47112, 1.48019,
            1.49065, 1.50121, 1.51054, 1.52129, 1.53079, 1.54037, 1.55004,
            1.56118, 1.57103, 1.58096, 1.59097, 1.60108, 1.61126, 1.62007,
            1.63041, 1.64084, 1.65135, 1.66042, 1.67105, 1.68019, 1.6909 ,
            1.70011, 1.71088, 1.72017, 1.73107, 1.7405 , 1.75158, 1.76115,
            1.7708 , 1.78052, 1.79032, 1.80019, 1.81014, 1.82017, 1.83028,
            1.84046, 1.85073, 1.86108, 1.8715 , 1.88026, 1.89084, 1.9015 ,
            1.91046, 1.92128, 1.93036, 1.94135, 1.95056, 1.96171, 1.97107,
            1.98049, 1.99188, 2.00144
        },
        x_q_array = {
            0.3698846 , 0.37365026, 0.37744542, 0.3809791 , 0.38483902,
            0.38843316, 0.39205698, 0.39571268, 0.39940015, 0.40280976,
            0.40656018, 0.4102268 , 0.41378238, 0.41741284, 0.4209029 ,
            0.42449549, 0.42779042, 0.43144557, 0.43513374, 0.43851581,
            0.44192718, 0.44536828, 0.44883836, 0.45233899, 0.45587006,
            0.45943357, 0.46302833, 0.46629096, 0.46994773, 0.47326703,
            0.47661519, 0.47999053, 0.48339403, 0.48682665, 0.4902898 ,
            0.49378136, 0.49730277, 0.50045971, 0.50404031, 0.50724909,
            0.51048278, 0.51415292, 0.51744154, 0.52075687, 0.52409881,
            0.52746732, 0.53086331, 0.53428672, 0.53730562, 0.54078234,
            0.54428778, 0.54737907, 0.55093846, 0.55407795, 0.55724028,
            0.56042683, 0.564096  , 0.56733275, 0.57059363, 0.57387812,
            0.57718847, 0.58052242, 0.58340108, 0.58678202, 0.59018922,
            0.59362276, 0.59658696, 0.60006916, 0.60307425, 0.60660508,
            0.60965306, 0.61323283, 0.61632309, 0.61995332, 0.62308544,
            0.62675627, 0.62991561, 0.63308711, 0.63626889, 0.63946426,
            0.64267283, 0.64589842, 0.64914297, 0.65240795, 0.65569401,
            0.65900274, 0.66233378, 0.66568529, 0.66849528, 0.67188428,
            0.67529244, 0.6781478 , 0.68159074, 0.68447412, 0.68795298,
            0.69086657, 0.69438342, 0.69733067, 0.70029326, 0.70386937,
            0.70686717
        };
    std::ofstream Fy_out("Fy_out_heat.txt"), Mz_out("Mz_out_heat.txt");
    std::vector<double> result;
    int N = int(x_q_array.size());
    for(int i = 0; i < 70; i++){
        std::cout << i << " iteration:" << std::endl;
        result = solver(x_q_array[i], z_q_array[i]);
        Fy_out << result[0] << " ";
        Mz_out << result[1] << " ";
        std::cout << "==================================\n==================================" << std::endl;
    }
}

int main()
{
    params.Mach_inf = 3;
    params.p_inf = 101330;
    params.rho_inf = 1.2255;

    params.a_inf = sqrt(Gamma * params.p_inf / params.rho_inf);
    params.V_inf = params.Mach_inf * params.a_inf;

    params.is_adiabatic = false;

    params.bodyType = BodyType::DoubleCone;

    // double Mach_inf = 3, p_inf = 101330, rho_inf = 1.2255;
    // double a_inf = sqrt(Gamma * p_inf / rho_inf);
    // double V_inf = Mach_inf*a_inf;
    // bool is_adiabatic = false;

    clock_t start = clock();
    // double z0 = 1 / tan(Pi / 12.0);
    solver(0.5, 1 + 0.1);
    clock_t finish = clock();
    double seconds = double(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time: " << seconds << "s" << std::endl;
    // test();

    return 0;
}