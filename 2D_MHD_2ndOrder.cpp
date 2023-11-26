#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "/opt/homebrew/Cellar/libomp/17.0.5/include/omp.h"

using namespace std;

double adiabatic_constant = 2.0;
double cfl_constant = 0.45;
double domain_size = 1.0;
int cell_number = 200;
double final_time = 0.2;
double bias = 0.0;
double c_h = 0.1; // hyperbolic wave speed
double c_p = sqrt(0.36); // parabolic damping term

vector<double> flux_x(vector<double> consVar)
{
    double density = consVar[0];
    double vel_x = consVar[1] / density;
    double vel_y = consVar[2] / density;
    double vel_z = consVar[3] / density;
    double magnetic_x = consVar[4];
    double magnetic_y = consVar[5];
    double magnetic_z = consVar[6];
    double energy = consVar[7];
    double phi = consVar[8];
    double pressure = (adiabatic_constant - 1.0) * (energy - (0.5 * density * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z))) - (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z))));
    double pressure_star = pressure + (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z)));
    
    vector<double> ret(9);
    
    ret[0] = density * vel_x;
    ret[1] = density * (vel_x * vel_x) + pressure_star - (magnetic_x * magnetic_x);
    ret[2] = density * (vel_x * vel_y) - (magnetic_x * magnetic_y);
    ret[3] = density * (vel_x * vel_z) - (magnetic_x * magnetic_z);
    ret[4] = (magnetic_x * vel_x) - (magnetic_x * vel_x) + phi;
    ret[5] = (magnetic_y * vel_x) - (magnetic_x * vel_y);
    ret[6] = (magnetic_z * vel_x) - (magnetic_x * vel_z);
    ret[7] = ((energy + pressure_star) * vel_x) - (magnetic_x * ((magnetic_x * vel_x) + (magnetic_y * vel_y) + (magnetic_z * vel_z)));
    ret[8] = (c_h * c_h) * magnetic_x;
    
    return ret;
}

vector<double> flux_y(vector<double> consVar)
{
    double density = consVar[0];
    double vel_x = consVar[1] / density;
    double vel_y = consVar[2] / density;
    double vel_z = consVar[3] / density;
    double magnetic_x = consVar[4];
    double magnetic_y = consVar[5];
    double magnetic_z = consVar[6];
    double energy = consVar[7];
    double phi = consVar[8];
    double pressure = (adiabatic_constant - 1.0) * (energy - (0.5 * density * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z))) - (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z))));
    double pressure_star = pressure + (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z)));
    
    vector<double> ret(9);
    
    ret[0] = density * vel_y;
    ret[1] = density * (vel_y * vel_x) - (magnetic_y * magnetic_x);
    ret[2] = density * (vel_y * vel_y) + pressure_star - (magnetic_y * magnetic_y);
    ret[3] = density * (vel_y * vel_z) - (magnetic_y * magnetic_z);
    ret[4] = (magnetic_x * vel_y) - (magnetic_y * vel_x);
    ret[5] = (magnetic_y * vel_y) - (magnetic_y * vel_y) + phi;
    ret[6] = (magnetic_z * vel_y) - (magnetic_y * vel_z);
    ret[7] = ((energy + pressure_star) * vel_y) - (magnetic_y * ((magnetic_x * vel_x) + (magnetic_y * vel_y) + (magnetic_z * vel_z)));
    ret[8] = (c_h * c_h) * magnetic_y;
    
    return ret;
}

vector<double> cons(double density, double vel_x, double vel_y, double vel_z, double magnetic_x, double magnetic_y, double magnetic_z, double energy, double phi)
{
    vector<double> ret(9);
    
    ret[0] = density;
    ret[1] = density * vel_x;
    ret[2] = density * vel_y;
    ret[3] = density * vel_z;
    ret[4] = magnetic_x;
    ret[5] = magnetic_y;
    ret[6] = magnetic_z;
    ret[7] = energy;
    ret[8] = phi;
    
    return ret;
}

vector<double> source(vector<double> consVar)
{
    vector<double> ret(9);
    
    for (int i = 0; i < 8; i++)
    {
        ret[i] = 0.0;
    }
    
    ret[8] = -((c_h * c_h) / (c_p * c_p)) * consVar[8];
    
    return ret;
}

vector<double> addVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] + vector2[i];
    }
    
    return ret;
}

vector<double> subtractVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] - vector2[i];
    }
    
    return ret;
}

vector<double> multiplyVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] * vector2[i];
    }
    
    return ret;
}

vector<double> divideVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] / vector2[i];
    }
    
    return ret;
}

vector<double> multiplyVector(double scalar, vector<double> vector1)
{
    vector<double> ret(vector1.size());

    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = scalar * vector1[i];
    }
    
    return ret;
}

vector<double> consVar_left_halfEvolve_x(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_left, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_x(consVar_left), flux_x(consVar_right))));
}

vector<double> consVar_right_halfEvolve_x(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_right, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_x(consVar_left), flux_x(consVar_right))));
}

vector<double> consVar_left_halfEvolve_y(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_left, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_y(consVar_left), flux_y(consVar_right))));
}

vector<double> consVar_right_halfEvolve_y(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_right, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_y(consVar_left), flux_y(consVar_right))));
}

vector<double> slope(vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right)
{
    return addVectors(multiplyVector((0.5 * (1.0 + bias)), subtractVectors(consVar, consVar_left)), multiplyVector((0.5 * (1.0 - bias)), subtractVectors(consVar_right, consVar)));
}

vector<double> slope_limiter(vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right)
{
    /*
    vector<double> r = divideVectors(
                  subtractVectors(slope(consVar_left, consVar, consVar_right), slope(consVar_left_left, consVar_left, consVar)),
                  subtractVectors(slope(consVar, consVar_right, consVar_right_right), slope(consVar_left, consVar, consVar_right)));
     */
    vector<double> r = divideVectors(subtractVectors(consVar, consVar_left), subtractVectors(consVar_right, consVar));
    
    vector<double> xiR(r.size());
    
    for (int i = 0; i < xiR.size(); i++)
    {
        xiR[i] = 2.0 / (1.0 - bias + (1.0 + bias) * r[i]);
    }
    
    vector<double> xi(r.size());
    
    for (int i = 0; i < xi.size(); i++)
    {
        if (r[i] <= 0.0)
        {
            xi[i] = 0.0;
        }
        else if (r[i] >= 0.0 && r[i] <= 0.5)
        {
            xi[i] = 2.0 * r[i];
        }
        else if (r[i] >= 0.5 && r[i] <= 1.0)
        {
            xi[i] = 1.0;
        }
        else if (r[i] >= 1.0)
        {
            xi[i] = min(min(r[i], xiR[i]), 2.0);
        }
    }
    
    return xi;
}

vector<double> extrapolate_boundary_left(vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right)
{
    return subtractVectors(consVar,
                           multiplyVector(0.5,
                                          multiplyVectors(
                                                          slope_limiter(consVar_left, consVar, consVar_right),
                                                          slope(consVar_left, consVar, consVar_right))));
}

vector<double> extrapolate_boundary_right(vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right)
{
    return addVectors(consVar,
                           multiplyVector(0.5,
                                          multiplyVectors(
                                                          slope_limiter(consVar_left, consVar, consVar_right),
                                                          slope(consVar_left, consVar, consVar_right))));
}

vector<double> force_LF(vector<double> flux_left, vector<double> flux_right, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(multiplyVector(0.5, addVectors(flux_left, flux_right)), multiplyVector((0.5 * (domain_size / cell_number) / dt), subtractVectors(consVar_left, consVar_right)));
}

vector<double> force_RI(vector<double> flux_left, vector<double> flux_right, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(multiplyVector(0.5, addVectors(consVar_left, consVar_right)), multiplyVector((0.5 * dt / (domain_size / cell_number)), subtractVectors(flux_left, flux_right)));
}

vector<double> force_x(vector<double> flux_left, vector<double> flux_right, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return multiplyVector(0.5, addVectors(force_LF(flux_left, flux_right, consVar_left, consVar_right, dt), flux_x(force_RI(flux_left, flux_right, consVar_left, consVar_right, dt))));
}

vector<double> force_y(vector<double> flux_left, vector<double> flux_right, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return multiplyVector(0.5, addVectors(force_LF(flux_left, flux_right, consVar_left, consVar_right, dt), flux_y(force_RI(flux_left, flux_right, consVar_left, consVar_right, dt))));
}


vector<double> updateSLIC_x(vector<double> consVar_left_left, vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right, vector<double> consVar_right_right, double dt)
{
    vector<double> consVar_left_boundaryL = extrapolate_boundary_left(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryR = extrapolate_boundary_right(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryL_evolved = consVar_left_halfEvolve_x(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    vector<double> consVar_left_boundaryR_evolved = consVar_right_halfEvolve_x(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    
    vector<double> consVar_boundaryL = extrapolate_boundary_left(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryR = extrapolate_boundary_right(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryL_evolved = consVar_left_halfEvolve_x(consVar_boundaryL, consVar_boundaryR, dt);
    vector<double> consVar_boundaryR_evolved = consVar_right_halfEvolve_x(consVar_boundaryL, consVar_boundaryR, dt);
    
    vector<double> consVar_right_boundaryL = extrapolate_boundary_left(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryR = extrapolate_boundary_right(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryL_evolved = consVar_left_halfEvolve_x(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    vector<double> consVar_right_boundaryR_evolved = consVar_right_halfEvolve_x(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)),
                                              subtractVectors(
                                                              force_x(flux_x(consVar_left_boundaryR_evolved), flux_x(consVar_boundaryL_evolved), consVar_left_boundaryR_evolved, consVar_boundaryL_evolved, dt),
                                                              force_x(flux_x(consVar_boundaryR_evolved), flux_x(consVar_right_boundaryL_evolved), consVar_boundaryR_evolved, consVar_right_boundaryL_evolved, dt))));
}

vector<double> updateSLIC_y(vector<double> consVar_left_left, vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right, vector<double> consVar_right_right, double dt)
{
    vector<double> consVar_left_boundaryL = extrapolate_boundary_left(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryR = extrapolate_boundary_right(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryL_evolved = consVar_left_halfEvolve_y(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    vector<double> consVar_left_boundaryR_evolved = consVar_right_halfEvolve_y(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    
    vector<double> consVar_boundaryL = extrapolate_boundary_left(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryR = extrapolate_boundary_right(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryL_evolved = consVar_left_halfEvolve_y(consVar_boundaryL, consVar_boundaryR, dt);
    vector<double> consVar_boundaryR_evolved = consVar_right_halfEvolve_y(consVar_boundaryL, consVar_boundaryR, dt);
    
    vector<double> consVar_right_boundaryL = extrapolate_boundary_left(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryR = extrapolate_boundary_right(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryL_evolved = consVar_left_halfEvolve_y(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    vector<double> consVar_right_boundaryR_evolved = consVar_right_halfEvolve_y(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)),
                                              subtractVectors(
                                                              force_y(flux_y(consVar_left_boundaryR_evolved), flux_y(consVar_boundaryL_evolved), consVar_left_boundaryR_evolved, consVar_boundaryL_evolved, dt),
                                                              force_y(flux_y(consVar_boundaryR_evolved), flux_y(consVar_right_boundaryL_evolved), consVar_boundaryR_evolved, consVar_right_boundaryL_evolved, dt))));
}

double time_step(vector<vector<double> > pressure, vector<vector<double> > density, vector<vector<double> > vel_x, vector<vector<double> > magnetic_x, vector<vector<double> > magnetic_y, vector<vector<double> > magnetic_z)
{
    double u_max = 0.0;

    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {
            double expr = ((adiabatic_constant * pressure[i][j]) + ((magnetic_x[i][j] * magnetic_x[i][j]) + (magnetic_y[i][j] * magnetic_y[i][j]) + (magnetic_z[i][j] * magnetic_z[i][j]))) / density[i][j];
            double c_f = sqrt(
                              0.5 * (expr + sqrt((expr * expr) - (4 * (adiabatic_constant * pressure[i][j] / density[i][j]) * (magnetic_x[i][j] * magnetic_x[i][j] / density[i][j])))));
            double speed = abs(vel_x[i][j]) + abs(c_f);
            
            if (speed > u_max)
            {
                u_max = speed;
            }
        }
        
    }
    
    return cfl_constant * (domain_size / cell_number) / u_max;
}

vector<vector<double> > pad(vector<vector<double> > quantity)
{
    vector<vector<double> > paddedQuantity(cell_number + 4, vector<double>(cell_number + 4));
    
    for (int i = 0; i < cell_number; i++)
    {
        paddedQuantity[i + 2][0] = quantity[i][0];
        paddedQuantity[i + 2][1] = quantity[i][0];
        
        for (int k = 0; k < cell_number; k++)
        {
            paddedQuantity[i + 2][k + 2] = quantity[i][k];
        }
        
        paddedQuantity[i + 2][cell_number + 2] = quantity[i][cell_number - 1];
        paddedQuantity[i + 2][cell_number + 3] = quantity[i][cell_number - 1];
        
    }
    
    for (int j = 0; j < cell_number; j++)
    {
        paddedQuantity[0][j + 2] = quantity[0][j];
        paddedQuantity[1][j + 2] = quantity[0][j];
        
        for (int k = 0; k < cell_number; k++)
        {
            paddedQuantity[k + 2][j + 2] = quantity[k][j];
        }
        
        paddedQuantity[cell_number + 2][j + 2] = quantity[cell_number - 1][j];
        paddedQuantity[cell_number + 3][j + 2] = quantity[cell_number - 1][j];
    }
    
    return paddedQuantity;
}

vector<vector<double> > divergence(vector<vector<double> > magnetic_x, vector<vector<double> > magnetic_y, vector<vector<double> > magnetic_z)
{
    vector<vector<double> > divB(cell_number, vector<double>(cell_number));
    
    vector<vector<double> > magnetic_x_padded = pad(magnetic_x);
    vector<vector<double> > magnetic_y_padded = pad(magnetic_y);
    vector<vector<double> > magnetic_z_padded = pad(magnetic_z);
    
    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {
            divB[i][j] = (magnetic_x_padded[i + 3][j + 2] - magnetic_x_padded[i + 2][j + 2]) / (domain_size / cell_number) + (magnetic_y_padded[i + 2][j + 3] - magnetic_y_padded[i + 2][j + 3]) / (domain_size / cell_number);
        }
    }
    
    return divB;
}

vector<double> update_x(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)), subtractVectors(force_x(flux_x(consVar_left), flux_x(consVar), consVar_left, consVar, dt), force_x(flux_x(consVar), flux_x(consVar_right), consVar, consVar_right, dt))));
}

vector<double> update_y(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)), subtractVectors(force_y(flux_y(consVar_left), flux_y(consVar), consVar_left, consVar, dt), force_y(flux_y(consVar), flux_y(consVar_right), consVar, consVar_right, dt))));
}

vector<double> rk4(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, vector<double> consVar_top, vector<double> consVar_bot, double dt)
{
    vector<double> k1 = multiplyVector(dt, source(consVar));
    vector<double> consVar1 = addVectors(consVar, multiplyVector(0.5, k1));
    consVar1 = update_x(consVar1, consVar_left, consVar_right, dt * 0.25);
    consVar1 = update_y(consVar1, consVar_top, consVar_bot, dt * 0.5);
    consVar1 = update_x(consVar1, consVar_left, consVar_right, dt * 0.25);
    
    vector<double> k2 = multiplyVector(dt, source(consVar1));
    vector<double> consVar2 = addVectors(consVar, multiplyVector(0.5, k2));
    consVar2 = update_x(consVar2, consVar_left, consVar_right, dt * 0.25);
    consVar2 = update_y(consVar2, consVar_top, consVar_bot, dt * 0.5);
    consVar2 = update_x(consVar2, consVar_left, consVar_right, dt * 0.25);
    
    vector<double> k3 = multiplyVector(dt, source(consVar2));
    vector<double> consVar3 = addVectors(consVar, k3);
    consVar3 = update_x(consVar3, consVar_left, consVar_right, dt * 0.5);
    consVar3 = update_y(consVar3, consVar_top, consVar_bot, dt);
    consVar3 = update_x(consVar3, consVar_left, consVar_right, dt * 0.5);
    
    vector<double> k4 = multiplyVector(dt, source(consVar3));
    
    return addVectors(consVar, multiplyVector((1.0 / 6.0), addVectors(k1, addVectors(multiplyVector(2.0, k2), addVectors(multiplyVector(2.0, k3), k4)))));
}

int main()
{
    vector<vector<double> > density(cell_number, vector<double>(cell_number));
    vector<vector<double> > pressure(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_x(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_y(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_z(cell_number, vector<double>(cell_number));
    vector<vector<double> > magnetic_x(cell_number, vector<double>(cell_number));
    vector<vector<double> > magnetic_y(cell_number, vector<double>(cell_number));
    vector<vector<double> > magnetic_z(cell_number, vector<double>(cell_number));
    vector<vector<double> > phi(cell_number, vector<double>(cell_number));
    
    double t = 0.0;
    int n = 1;
    
    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {
            if ((sqrt(((i - (cell_number * 0.5)) * (i - (cell_number * 0.5))) + ((j - (cell_number * 0.5)) * (j - (cell_number * 0.5)))) * (domain_size / cell_number)) < 0.25)
            {
                density[i][j] = 1.0;
                vel_x[i][j] = 0.0;
                vel_y[i][j] = 0.0;
                vel_z[i][j] = 0.0;
                magnetic_x[i][j] = 0.75;
                magnetic_y[i][j] = 1.0;
                magnetic_z[i][j] = 0.0;
                pressure[i][j] = 1.0;
                phi[i][j] = 0.0;
            }
            else
            {
                density[i][j] = 0.125;
                vel_x[i][j] = 0.0;
                vel_y[i][j] = 0.0;
                vel_z[i][j] = 0.0;
                magnetic_x[i][j] = 0.75;
                magnetic_y[i][j] = -1.0;
                magnetic_z[i][j] = 0.0;
                pressure[i][j] = 0.1;
                phi[i][j] = 0.0;
            }
        }
    }
    
    while (t < final_time)
    {
        double dt = time_step(pressure, density, vel_x, magnetic_x, magnetic_y, magnetic_z);
        cout << "n = " << n << "; t = " << t << "; dt = " << dt << "\n";
        if (n <= 5)
        {
            dt *= 0.2;
        }
        t += dt;
        n += 1;
        
        vector<vector<double> > divB = divergence(magnetic_x, magnetic_y, magnetic_z);
        double divB_max = 0.0;
        
        for (int i = 0; i < cell_number; i++)
        {
            for (int j = 0; j < cell_number; j++)
            {
                if (abs(divB[i][j]) > divB_max)
                {
                    divB_max = abs(divB[i][j]);
                }
            }
        }
        
        cout << "divB_max = " << divB_max << "\n";
        
        // x - direction - 1/2
        
        vector<vector<double> > density_padded = pad(density);
        vector<vector<double> > pressure_padded = pad(pressure);
        vector<vector<double> > vel_x_padded = pad(vel_x);
        vector<vector<double> > vel_y_padded = pad(vel_y);
        vector<vector<double> > vel_z_padded = pad(vel_z);
        vector<vector<double> > magnetic_x_padded = pad(magnetic_x);
        vector<vector<double> > magnetic_y_padded = pad(magnetic_y);
        vector<vector<double> > magnetic_z_padded = pad(magnetic_z);
        vector<vector<double> > phi_padded = pad(phi);
        
        vector<vector<double> > density_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > pressure_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_x_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_y_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_z_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > magnetic_x_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > magnetic_y_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > magnetic_z_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > phi_new(cell_number + 4, vector<double>(cell_number + 4));
        
        vector<vector<double> > energy_tot(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > energy_new(cell_number + 4, vector<double>(cell_number + 4));
        
        for (int i = 0; i < cell_number + 4; i++)
        {
            for (int j = 0; j < cell_number + 4; j++)
            {
                double energy_int = pressure_padded[i][j] / (density_padded[i][j] * (adiabatic_constant - 1.0));
                
                energy_tot[i][j] = (0.5 * density_padded[i][j] * ((vel_x_padded[i][j] * vel_x_padded[i][j]) + (vel_y_padded[i][j] * vel_y_padded[i][j]) + (vel_z_padded[i][j] * vel_z_padded[i][j]))) + (0.5 * ((magnetic_x_padded[i][j] * magnetic_x_padded[i][j]) + (magnetic_y_padded[i][j] * magnetic_y_padded[i][j]) + (magnetic_z_padded[i][j] * magnetic_z_padded[i][j]))) + (density_padded[i][j] * energy_int);
            }
            
        }
        
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], magnetic_x_padded[i][j], magnetic_y_padded[i][j], magnetic_z_padded[i][j], energy_tot[i][j], phi_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i - 1][j], vel_x_padded[i - 1][j], vel_y_padded[i - 1][j], vel_z_padded[i - 1][j], magnetic_x_padded[i - 1][j], magnetic_y_padded[i - 1][j], magnetic_z_padded[i - 1][j], energy_tot[i - 1][j], phi_padded[i - 1][j]);
                vector<double> consVar_left_left = cons(density_padded[i - 2][j], vel_x_padded[i - 2][j], vel_y_padded[i - 2][j], vel_z_padded[i - 2][j], magnetic_x_padded[i - 2][j], magnetic_y_padded[i - 2][j], magnetic_z_padded[i - 2][j], energy_tot[i - 2][j], phi_padded[i - 2][j]);
                vector<double> consVar_right = cons(density_padded[i + 1][j], vel_x_padded[i + 1][j], vel_y_padded[i + 1][j], vel_z_padded[i + 1][j], magnetic_x_padded[i + 1][j], magnetic_y_padded[i + 1][j], magnetic_z_padded[i + 1][j], energy_tot[i + 1][j], phi_padded[i + 1][j]);
                vector<double> consVar_right_right = cons(density_padded[i + 2][j], vel_x_padded[i + 2][j], vel_y_padded[i + 2][j], vel_z_padded[i + 2][j], magnetic_x_padded[i + 2][j], magnetic_y_padded[i + 2][j], magnetic_z_padded[i + 2][j], energy_tot[i + 2][j], phi_padded[i + 2][j]);
                
                vector<double> newConsVar = updateSLIC_x(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt * 0.5);
                
                density_new[i][j] = newConsVar[0];
                vel_x_new[i][j] = newConsVar[1] / density_new[i][j];
                vel_y_new[i][j] = newConsVar[2] / density_new[i][j];
                vel_z_new[i][j] = newConsVar[3] / density_new[i][j];
                magnetic_x_new[i][j] = newConsVar[4];
                magnetic_y_new[i][j] = newConsVar[5];
                magnetic_z_new[i][j] = newConsVar[6];
                energy_new[i][j] = newConsVar[7];
                pressure_new[i][j] = (adiabatic_constant - 1.0) * (energy_new[i][j] - (0.5 * density_new[i][j] * ((vel_x_new[i][j] * vel_x_new[i][j]) + (vel_y_new[i][j] * vel_y_new[i][j]) + (vel_z_new[i][j] * vel_z_new[i][j]))) - (0.5 * ((magnetic_x_new[i][j] * magnetic_x_new[i][j]) + (magnetic_y_new[i][j] * magnetic_y_new[i][j]) + (magnetic_z_new[i][j] * magnetic_z_new[i][j]))));
                phi_new[i][j] = newConsVar[8];
            }
        }
        
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                density[i - 2][j - 2] = density_new[i][j];
                pressure[i - 2][j - 2] = pressure_new[i][j];
                vel_x[i - 2][j - 2] = vel_x_new[i][j];
                vel_y[i - 2][j - 2] = vel_y_new[i][j];
                vel_z[i - 2][j - 2] = vel_z_new[i][j];
                magnetic_x[i - 2][j - 2] = magnetic_x_new[i][j];
                magnetic_y[i - 2][j - 2] = magnetic_y_new[i][j];
                magnetic_z[i - 2][j - 2] = magnetic_z_new[i][j];
                phi[i - 2][j - 2] = phi_new[i][j];
            }
        }
        
        // y direction - 1
        
        density_padded = pad(density);
        pressure_padded = pad(pressure);
        vel_x_padded = pad(vel_x);
        vel_y_padded = pad(vel_y);
        vel_z_padded = pad(vel_z);
        magnetic_x_padded = pad(magnetic_x);
        magnetic_y_padded = pad(magnetic_y);
        magnetic_z_padded = pad(magnetic_z);
        phi_padded = pad(phi);
        
        for (int i = 0; i < cell_number + 4; i++)
        {
            for (int j = 0; j < cell_number + 4; j++)
            {
                double energy_int = pressure_padded[i][j] / (density_padded[i][j] * (adiabatic_constant - 1.0));
                
                energy_tot[i][j] = (0.5 * density_padded[i][j] * ((vel_x_padded[i][j] * vel_x_padded[i][j]) + (vel_y_padded[i][j] * vel_y_padded[i][j]) + (vel_z_padded[i][j] * vel_z_padded[i][j]))) + (0.5 * ((magnetic_x_padded[i][j] * magnetic_x_padded[i][j]) + (magnetic_y_padded[i][j] * magnetic_y_padded[i][j]) + (magnetic_z_padded[i][j] * magnetic_z_padded[i][j]))) + (density_padded[i][j] * energy_int);
            }
            
        }
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], magnetic_x_padded[i][j], magnetic_y_padded[i][j], magnetic_z_padded[i][j], energy_tot[i][j], phi_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i][j - 1], vel_x_padded[i][j - 1], vel_y_padded[i][j - 1], vel_z_padded[i][j - 1], magnetic_x_padded[i][j - 1], magnetic_y_padded[i][j - 1], magnetic_z_padded[i][j - 1], energy_tot[i][j - 1], phi_padded[i][j - 1]);
                vector<double> consVar_left_left = cons(density_padded[i][j - 2], vel_x_padded[i][j - 2], vel_y_padded[i][j - 2], vel_z_padded[i][j - 2], magnetic_x_padded[i][j - 2], magnetic_y_padded[i][j - 2], magnetic_z_padded[i][j - 2], energy_tot[i][j - 2], phi_padded[i][j - 2]);
                vector<double> consVar_right = cons(density_padded[i][j + 1], vel_x_padded[i][j + 1], vel_y_padded[i][j + 1], vel_z_padded[i][j + 1], magnetic_x_padded[i][j + 1], magnetic_y_padded[i][j + 1], magnetic_z_padded[i][j + 1], energy_tot[i][j + 1], phi_padded[i][j + 1]);
                vector<double> consVar_right_right = cons(density_padded[i][j + 2], vel_x_padded[i][j + 2], vel_y_padded[i][j + 2], vel_z_padded[i][j + 2], magnetic_x_padded[i][j + 2], magnetic_y_padded[i][j + 2], magnetic_z_padded[i][j + 2], energy_tot[i][j + 2], phi_padded[i][j + 2]);
                
                vector<double> newConsVar = updateSLIC_y(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt);
                
                density_new[i][j] = newConsVar[0];
                vel_x_new[i][j] = newConsVar[1] / density_new[i][j];
                vel_y_new[i][j] = newConsVar[2] / density_new[i][j];
                vel_z_new[i][j] = newConsVar[3] / density_new[i][j];
                magnetic_x_new[i][j] = newConsVar[4];
                magnetic_y_new[i][j] = newConsVar[5];
                magnetic_z_new[i][j] = newConsVar[6];
                energy_new[i][j] = newConsVar[7];
                pressure_new[i][j] = (adiabatic_constant - 1.0) * (energy_new[i][j] - (0.5 * density_new[i][j] * ((vel_x_new[i][j] * vel_x_new[i][j]) + (vel_y_new[i][j] * vel_y_new[i][j]) + (vel_z_new[i][j] * vel_z_new[i][j]))) - (0.5 * ((magnetic_x_new[i][j] * magnetic_x_new[i][j]) + (magnetic_y_new[i][j] * magnetic_y_new[i][j]) + (magnetic_z_new[i][j] * magnetic_z_new[i][j]))));
                phi_new[i][j] = newConsVar[8];
            }
        }
        
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                density[i - 2][j - 2] = density_new[i][j];
                pressure[i - 2][j - 2] = pressure_new[i][j];
                vel_x[i - 2][j - 2] = vel_x_new[i][j];
                vel_y[i - 2][j - 2] = vel_y_new[i][j];
                vel_z[i - 2][j - 2] = vel_z_new[i][j];
                magnetic_x[i - 2][j - 2] = magnetic_x_new[i][j];
                magnetic_y[i - 2][j - 2] = magnetic_y_new[i][j];
                magnetic_z[i - 2][j - 2] = magnetic_z_new[i][j];
                phi[i - 2][j - 2] = phi_new[i][j];
            }
        }
        
        
        // x - direction - 1/2
        
        density_padded = pad(density);
        pressure_padded = pad(pressure);
        vel_x_padded = pad(vel_x);
        vel_y_padded = pad(vel_y);
        vel_z_padded = pad(vel_z);
        magnetic_x_padded = pad(magnetic_x);
        magnetic_y_padded = pad(magnetic_y);
        magnetic_z_padded = pad(magnetic_z);
        phi_padded = pad(phi);
        
        for (int i = 0; i < cell_number + 4; i++)
        {
            for (int j = 0; j < cell_number + 4; j++)
            {
                double energy_int = pressure_padded[i][j] / (density_padded[i][j] * (adiabatic_constant - 1.0));
                
                energy_tot[i][j] = (0.5 * density_padded[i][j] * ((vel_x_padded[i][j] * vel_x_padded[i][j]) + (vel_y_padded[i][j] * vel_y_padded[i][j]) + (vel_z_padded[i][j] * vel_z_padded[i][j]))) + (0.5 * ((magnetic_x_padded[i][j] * magnetic_x_padded[i][j]) + (magnetic_y_padded[i][j] * magnetic_y_padded[i][j]) + (magnetic_z_padded[i][j] * magnetic_z_padded[i][j]))) + (density_padded[i][j] * energy_int);
            }
            
        }
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], magnetic_x_padded[i][j], magnetic_y_padded[i][j], magnetic_z_padded[i][j], energy_tot[i][j], phi_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i - 1][j], vel_x_padded[i - 1][j], vel_y_padded[i - 1][j], vel_z_padded[i - 1][j], magnetic_x_padded[i - 1][j], magnetic_y_padded[i - 1][j], magnetic_z_padded[i - 1][j], energy_tot[i - 1][j], phi_padded[i - 1][j]);
                vector<double> consVar_left_left = cons(density_padded[i - 2][j], vel_x_padded[i - 2][j], vel_y_padded[i - 2][j], vel_z_padded[i - 2][j], magnetic_x_padded[i - 2][j], magnetic_y_padded[i - 2][j], magnetic_z_padded[i - 2][j], energy_tot[i - 2][j], phi_padded[i - 2][j]);
                vector<double> consVar_right = cons(density_padded[i + 1][j], vel_x_padded[i + 1][j], vel_y_padded[i + 1][j], vel_z_padded[i + 1][j], magnetic_x_padded[i + 1][j], magnetic_y_padded[i + 1][j], magnetic_z_padded[i + 1][j], energy_tot[i + 1][j], phi_padded[i + 1][j]);
                vector<double> consVar_right_right = cons(density_padded[i + 2][j], vel_x_padded[i + 2][j], vel_y_padded[i + 2][j], vel_z_padded[i + 2][j], magnetic_x_padded[i + 2][j], magnetic_y_padded[i + 2][j], magnetic_z_padded[i + 2][j], energy_tot[i + 2][j], phi_padded[i + 2][j]);
                
                vector<double> newConsVar = updateSLIC_x(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt * 0.5);
                
                density_new[i][j] = newConsVar[0];
                vel_x_new[i][j] = newConsVar[1] / density_new[i][j];
                vel_y_new[i][j] = newConsVar[2] / density_new[i][j];
                vel_z_new[i][j] = newConsVar[3] / density_new[i][j];
                magnetic_x_new[i][j] = newConsVar[4];
                magnetic_y_new[i][j] = newConsVar[5];
                magnetic_z_new[i][j] = newConsVar[6];
                energy_new[i][j] = newConsVar[7];
                pressure_new[i][j] = (adiabatic_constant - 1.0) * (energy_new[i][j] - (0.5 * density_new[i][j] * ((vel_x_new[i][j] * vel_x_new[i][j]) + (vel_y_new[i][j] * vel_y_new[i][j]) + (vel_z_new[i][j] * vel_z_new[i][j]))) - (0.5 * ((magnetic_x_new[i][j] * magnetic_x_new[i][j]) + (magnetic_y_new[i][j] * magnetic_y_new[i][j]) + (magnetic_z_new[i][j] * magnetic_z_new[i][j]))));
                phi_new[i][j] = newConsVar[8];
            }
        }
        
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                density[i - 2][j - 2] = density_new[i][j];
                pressure[i - 2][j - 2] = pressure_new[i][j];
                vel_x[i - 2][j - 2] = vel_x_new[i][j];
                vel_y[i - 2][j - 2] = vel_y_new[i][j];
                vel_z[i - 2][j - 2] = vel_z_new[i][j];
                magnetic_x[i - 2][j - 2] = magnetic_x_new[i][j];
                magnetic_y[i - 2][j - 2] = magnetic_y_new[i][j];
                magnetic_z[i - 2][j - 2] = magnetic_z_new[i][j];
                phi[i - 2][j - 2] = phi_new[i][j];
            }
        }
        
        // Source update
        
        density_padded = pad(density);
        pressure_padded = pad(pressure);
        vel_x_padded = pad(vel_x);
        vel_y_padded = pad(vel_y);
        vel_z_padded = pad(vel_z);
        magnetic_x_padded = pad(magnetic_x);
        magnetic_y_padded = pad(magnetic_y);
        magnetic_z_padded = pad(magnetic_z);
        phi_padded = pad(phi);
        
#pragma omp parallel for
        for (int i = 2; i < cell_number + 2; i++)
        {
            for (int j = 2; j < cell_number + 2; j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], magnetic_x_padded[i][j], magnetic_y_padded[i][j], magnetic_z_padded[i][j], energy_tot[i][j], phi_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i - 1][j], vel_x_padded[i - 1][j], vel_y_padded[i - 1][j], vel_z_padded[i - 1][j], magnetic_x_padded[i - 1][j], magnetic_y_padded[i - 1][j], magnetic_z_padded[i - 1][j], energy_tot[i - 1][j], phi_padded[i - 1][j]);
                vector<double> consVar_right = cons(density_padded[i + 1][j], vel_x_padded[i + 1][j], vel_y_padded[i + 1][j], vel_z_padded[i + 1][j], magnetic_x_padded[i + 1][j], magnetic_y_padded[i + 1][j], magnetic_z_padded[i + 1][j], energy_tot[i + 1][j], phi_padded[i + 1][j]);
                vector<double> consVar_top = cons(density_padded[i][j - 1], vel_x_padded[i][j - 1], vel_y_padded[i][j - 1], vel_z_padded[i][j - 1], magnetic_x_padded[i][j - 1], magnetic_y_padded[i][j - 1], magnetic_z_padded[i][j - 1], energy_tot[i][j - 1], phi_padded[i][j - 1]);
                vector<double> consVar_bot = cons(density_padded[i][j + 1], vel_x_padded[i][j + 1], vel_y_padded[i][j + 1], vel_z_padded[i][j + 1], magnetic_x_padded[i][j + 1], magnetic_y_padded[i][j + 1], magnetic_z_padded[i][j + 1], energy_tot[i][j + 1], phi_padded[i][j + 1]);
                
                consVar = rk4(consVar, consVar_left, consVar_right, consVar_top, consVar_bot, dt);
                
                phi_new[i][j] = consVar[8];
            }
        }
        
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                density[i - 2][j - 2] = density_new[i][j];
                pressure[i - 2][j - 2] = pressure_new[i][j];
                vel_x[i - 2][j - 2] = vel_x_new[i][j];
                vel_y[i - 2][j - 2] = vel_y_new[i][j];
                vel_z[i - 2][j - 2] = vel_z_new[i][j];
                magnetic_x[i - 2][j - 2] = magnetic_x_new[i][j];
                magnetic_y[i - 2][j - 2] = magnetic_y_new[i][j];
                magnetic_z[i - 2][j - 2] = magnetic_z_new[i][j];
                phi[i - 2][j - 2] = phi_new[i][j];
            }
        }
    }
    
    ofstream density_file, pressure_file, velocity_x_file, velocity_y_file, velocity_z_file, magnetic_x_file, magnetic_y_file, magnetic_z_file, energy_int_file, divergence_file;
    density_file.open("density.dat");
    pressure_file.open("pressure.dat");
    velocity_x_file.open("velocity_x.dat");
    velocity_y_file.open("velocity_y.dat");
    velocity_z_file.open("velocity_z.dat");
    magnetic_x_file.open("magnetic_x.dat");
    magnetic_y_file.open("magnetic_y.dat");
    magnetic_z_file.open("magnetic_z.dat");
    energy_int_file.open("energy.dat");
    divergence_file.open("divergence.dat");
    
    vector<vector<double> > divB = divergence(magnetic_x, magnetic_y, magnetic_z);
    
    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {
            density_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << density[i][j] << "\n";
            pressure_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << pressure[i][j] << "\n";
            velocity_x_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << vel_x[i][j] << "\n";
            velocity_y_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << vel_y[i][j] << "\n";
            velocity_z_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << vel_z[i][j] << "\n";
            magnetic_x_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << magnetic_x[i][j] << "\n";
            magnetic_y_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << magnetic_y[i][j] << "\n";
            magnetic_z_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << magnetic_z[i][j] << "\n";
            energy_int_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << pressure[i][j] / (density[i][j] * (adiabatic_constant - 1.0)) << "\n";
            divergence_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << phi[i][j] << "\n";
        }
        density_file << "\n";
        pressure_file << "\n";
        velocity_x_file << "\n";
        velocity_y_file << "\n";
        velocity_z_file << "\n";
        magnetic_x_file << "\n";
        magnetic_y_file << "\n";
        magnetic_z_file << "\n";
        energy_int_file << "\n";
        divergence_file << "\n";
    }
    
    density_file.close();
    pressure_file.close();
    velocity_x_file.close();
    velocity_y_file.close();
    velocity_z_file.close();
    magnetic_x_file.close();
    magnetic_y_file.close();
    magnetic_z_file.close();
    energy_int_file.close();
    divergence_file.close();
    
    return 0;
}
