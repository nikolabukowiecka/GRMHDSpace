#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
#include "/opt/homebrew/Cellar/libomp/17.0.5/include/omp.h"

double adiabatic_constant = 5.0 / 3.0;
double cfl_constant = 0.9;
double domain_size = 1.0;
int cell_number = 100;
double final_time = 0.2;
double bias = 0.0;

vector<double> prim(double density_rel, double momentum_rel_x, double momentum_rel_y, double momentum_rel_z, double energy_rel)
{
    double momentum = (momentum_rel_x * momentum_rel_x) + (momentum_rel_y * momentum_rel_y) + (momentum_rel_z * momentum_rel_z);
    
    double c0 = (energy_rel + density_rel) / sqrt((energy_rel + density_rel) * (energy_rel + density_rel) - momentum);
    double c = density_rel / sqrt((energy_rel + density_rel) * (energy_rel + density_rel) - momentum);
    
    if ((energy_rel + density_rel) * (energy_rel + density_rel) - momentum < pow(10.0, -8.0))
    {
        c0 = (energy_rel + density_rel) / sqrt(pow(10.0, -8.0));
        c = density_rel / sqrt(pow(10.0, -8.0));
    }
    
    double a0 = - 1.0 / (adiabatic_constant * adiabatic_constant);
    double a1 = - (2.0 * c) * (adiabatic_constant - 1.0) / (adiabatic_constant * adiabatic_constant);
    double a2 = ((adiabatic_constant - 2.0) / adiabatic_constant) * ((c0 * c0) - 1.0) + 1.0 - (c * c) * ((adiabatic_constant - 1.0) / adiabatic_constant) * ((adiabatic_constant - 1.0) / adiabatic_constant);
    double a4 = (c0 * c0) - 1.0;
    double eta = 2.0 * c * (adiabatic_constant - 1.0) / adiabatic_constant;
    
    double xi = 1.0;
    int n = 0;
    while(n < 1000)
    {
        double fn = (a4 * (xi * xi * xi) * (xi - eta)) + (a2 * (xi * xi)) + (a1 * xi) + a0;
        double fn_derivative = a1 + (xi * ((2.0 * a2) + (a4 * xi * ((4.0 * xi) - (3.0 * eta)))));
        double xi_plus = xi - (fn / fn_derivative);
        
        if (abs(xi - xi_plus) < pow(10.0, -8.0))
        {
            n = 1000;
        }
        else
        {
            xi = xi_plus;
            n++;
        }
    }
        
    double lorentz = 0.5 * c0 * xi * (1.0 + sqrt(1.0 + (4.0 * ((adiabatic_constant - 1.0) / adiabatic_constant) * ((1.0 - (c * xi)) / ((c0 * c0) * (xi * xi))))));
    double enthalpy = 1.0 / (c * xi);
    
    vector<double> ret(7);
    ret[0] = density_rel / lorentz;
    ret[1] = momentum_rel_x / (ret[0] * enthalpy * (lorentz * lorentz));
    ret[2] = momentum_rel_y / (ret[0] * enthalpy * (lorentz * lorentz));
    ret[3] = momentum_rel_z / (ret[0] * enthalpy * (lorentz * lorentz));
    //ret[4] = enthalpy * ret[0] * ((adiabatic_constant - 1.0) / adiabatic_constant) - ret[0];
    ret[4] = ret[0] * enthalpy * (lorentz * lorentz) - energy_rel - density_rel;
    ret[5] = lorentz;
    ret[6] = enthalpy;
    
    if (ret[0] < pow(10.0, -8.0))
    {
        ret[0] = pow(10.0, -8.0);
    }
    if (ret[4] < pow(10.0, -8.0))
    {
        ret[4] = pow(10.0, -8.0);
    }
    
    return ret;
}

vector<double> flux_x(vector<double> consVar)
{
    double density_rel = consVar[0];
    double momentum_x = consVar[1];
    double momentum_y = consVar[2];
    double momentum_z = consVar[3];
    double energy_rel = consVar[4];
    
    vector<double> primVar = prim(density_rel, momentum_x, momentum_y, momentum_z, energy_rel);

    double density = primVar[0];
    double vel_x = primVar[1];
    double vel_y = primVar[2];
    double vel_z = primVar[3];
    double pressure = primVar[4];
    double lorentz = primVar[5];
    double enthalpy = primVar[6];
    
    /*
    double lorentz = 1.0 / sqrt(1.0 - ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)));
    double enthalpy = 1.0 + ((pressure / density) * (adiabatic_constant / (adiabatic_constant - 1.0)));
     */
    
    vector<double> ret(5);
    
    ret[0] = density * lorentz * vel_x;
    ret[1] = (density * enthalpy * (lorentz * lorentz) * vel_x * vel_x) + pressure;
    ret[2] = (density * enthalpy * (lorentz * lorentz) * vel_y * vel_x);
    ret[3] = (density * enthalpy * (lorentz * lorentz) * vel_z * vel_x);
    ret[4] = ((density * enthalpy * (lorentz * lorentz)) - pressure - (density * lorentz)) * vel_x + (pressure * vel_x);
    
    return ret;
}

vector<double> flux_y(vector<double> consVar)
{
    double density_rel = consVar[0];
    double momentum_x = consVar[1];
    double momentum_y = consVar[2];
    double momentum_z = consVar[3];
    double energy_rel = consVar[4];
    
    vector<double> primVar = prim(density_rel, momentum_x, momentum_y, momentum_z, energy_rel);

    double density = primVar[0];
    double vel_x = primVar[1];
    double vel_y = primVar[2];
    double vel_z = primVar[3];
    double pressure = primVar[4];
    double lorentz = primVar[5];
    double enthalpy = primVar[6];
    
    /*
    double lorentz = 1.0 / sqrt(1.0 - ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)));
    double enthalpy = 1.0 + ((pressure / density) * (adiabatic_constant / (adiabatic_constant - 1.0)));
     */
    
    vector<double> ret(5);
    
    ret[0] = density * lorentz * vel_y;
    ret[1] = (density * enthalpy * (lorentz * lorentz) * vel_x * vel_y);
    ret[2] = (density * enthalpy * (lorentz * lorentz) * vel_y * vel_y) + pressure;
    ret[3] = (density * enthalpy * (lorentz * lorentz) * vel_z * vel_y);
    ret[4] = ((density * enthalpy * (lorentz * lorentz)) - pressure - (density * lorentz)) * vel_y + (pressure * vel_y);
    
    return ret;
}


vector<double> cons(double density, double vel_x, double vel_y, double vel_z, double pressure)
{
    vector<double> ret(5);
    double lorentz = 1.0 / sqrt(1.0 - ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)));
    
    if (((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)) > 1.0 - pow(10.0, -8.0))
    {
        lorentz = 1.0 / sqrt(1.0 - pow(10.0, -8.0));
    }
    if (density < pow(10.0, -8.0))
    {
        density = pow(10.0, -8.0);
    }
    if (pressure < pow(10.0, -8.0))
    {
        pressure = pow(10.0, -8.0);
    }
    
    double enthalpy = 1.0 + ((pressure / density) * (adiabatic_constant / (adiabatic_constant - 1.0)));
                             
    ret[0] = density * lorentz;
    ret[1] = density * enthalpy * (lorentz * lorentz) * vel_x;
    ret[2] = density * enthalpy * (lorentz * lorentz) * vel_y;
    ret[3] = density * enthalpy * (lorentz * lorentz) * vel_z;
    ret[4] = (density * enthalpy * (lorentz * lorentz)) - pressure - (density * lorentz);
    
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
        if (abs(vector1[i]) < pow(10.0, -5.0))
        {
            vector1[i] = pow(10.0, -5.0);
        }
        if (abs(vector2[i]) < pow(10.0, -5.0))
        {
            vector2[i] = pow(10.0, -5.0);
        }
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

vector<double> consVar_left_halfEvolve_y(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_left, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_y(consVar_left), flux_y(consVar_right))));
}

vector<double> consVar_right_halfEvolve_x(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_right, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux_x(consVar_left), flux_x(consVar_right))));
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
    
    /*
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
     */
    
    vector<double> xi(r.size());
    for (int i = 0; i < xi.size(); i++)
    {
        if (r[i] < 0.0)
        {
            xi[i] = 0.0;
        }
        else
        {
            xi[i] = min((2.0 * r[i]) / (1.0 + r[i]), xiR[i]);
        }
    }
    
    /*
    vector<double> xi(r.size());
    for (int i = 0; i < xi.size(); i++)
    {
        if (r[i] <= 0.0)
        {
            xi[i] = 0.0;
        }
        else if (r[i] > 0.0 && r[i] <= 1.0)
        {
            xi[i] = r[i];
        }
        else
        {
            xi[i] = min(1.0, xiR[i]);
        }
    }
     */
    
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

vector<double> update_x(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, vector<double> fluxVect, vector<double> fluxVect_left, vector<double> fluxVect_right, double dt)
{
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)), subtractVectors(force_x(fluxVect_left, fluxVect, consVar_left, consVar, dt), force_x(fluxVect, fluxVect_right, consVar, consVar_right, dt))));
}

vector<double> update_y(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, vector<double> fluxVect, vector<double> fluxVect_left, vector<double> fluxVect_right, double dt)
{
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)), subtractVectors(force_y(fluxVect_left, fluxVect, consVar_left, consVar, dt), force_y(fluxVect, fluxVect_right, consVar, consVar_right, dt))));
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



double time_step(vector<vector<double> > pressure, vector<vector<double> > density, vector<vector<double> > vel_x, vector<vector<double> > vel_y)
{
    double u_max = 0.0;

    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {
            double num = (adiabatic_constant * pressure[i][j]) / density[i][j];
            double den = 1.0 + ((pressure[i][j] / density[i][j]) * (adiabatic_constant / (adiabatic_constant - 1.0)));
            double a = sqrt(num / den);
            double speed = max(abs(vel_x[i][j]), abs(vel_y[i][j])) + abs(a);
            
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

int main()
{
    /*
    double density = 1.0;
    double vel_x = 0.05;
    double vel_y = 0.0;
    double vel_z = 0.0;
    double pressure = 0;
    
    vector<double> consVar = cons(density, vel_x, vel_y, vel_z, pressure);
    vector<double> primsVar = prim(consVar[0], consVar[1], consVar[2], consVar[3], consVar[4]);
    
    cout << primsVar[0] << " " << primsVar[1] << " " << primsVar[2] << " " << primsVar[3] << " " << primsVar[4];
     */
    
    vector<vector<double> > density(cell_number, vector<double>(cell_number));
    vector<vector<double> > pressure(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_x(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_y(cell_number, vector<double>(cell_number));
    vector<vector<double> > vel_z(cell_number, vector<double>(cell_number));
    
    double t = 0.0;
    int n = 1;
    
    for (int i = 0; i < cell_number; i++)
    {
        for (int j = 0; j < cell_number; j++)
        {

            double r = sqrt((((i * (domain_size / cell_number)) - 0.5) * ((i * (domain_size / cell_number)) - 0.5)) + (((j * (domain_size / cell_number)) - 0.5) * ((j * (domain_size / cell_number)) - 0.5)));
            
            if (r <= 0.1)
            {
                density[i][j] = 10.0;
                pressure[i][j] = 13.3;
            }
            else if (r > 0.1)
            {
                density[i][j] = 1.0;
                pressure[i][j] = 0.0;
            }
            
            vel_x[i][j] = 0.0;
            vel_y[i][j] = 0.0;
            vel_z[i][j] = 0.0;
            
        }
    }
    
    while (t < final_time)
    {
        double dt = time_step(pressure, density, vel_x, vel_y);
        cout << "n = " << n << "; t = " << t << "; dt = " << dt << "\n";
        if (n <= 5)
        {
            dt *= 0.2;
        }
        t += dt;
        n += 1;
        
        // x direction (1/2 step)
        vector<vector<double> > density_padded = pad(density);
        vector<vector<double> > pressure_padded = pad(pressure);
        vector<vector<double> > vel_x_padded = pad(vel_x);
        vector<vector<double> > vel_y_padded = pad(vel_y);
        vector<vector<double> > vel_z_padded = pad(vel_z);
        
        vector<vector<double> > density_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > pressure_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_x_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_y_new(cell_number + 4, vector<double>(cell_number + 4));
        vector<vector<double> > vel_z_new(cell_number + 4, vector<double>(cell_number + 4));
        
#pragma omp parallel for
        
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], pressure_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i - 1][j], vel_x_padded[i - 1][j], vel_y_padded[i - 1][j], vel_z_padded[i - 1][j], pressure_padded[i - 1][j]);
                vector<double> consVar_left_left = cons(density_padded[i - 2][j], vel_x_padded[i - 2][j], vel_y_padded[i - 2][j], vel_z_padded[i - 2][j], pressure_padded[i - 2][j]);
                vector<double> consVar_right = cons(density_padded[i + 1][j], vel_x_padded[i + 1][j], vel_y_padded[i + 1][j], vel_z_padded[i + 1][j], pressure_padded[i + 1][j]);
                vector<double> consVar_right_right = cons(density_padded[i + 2][j], vel_x_padded[i + 2][j], vel_y_padded[i + 2][j], vel_z_padded[i + 2][j], pressure_padded[i + 2][j]);
                
                vector<double> newConsVar = updateSLIC_x(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt * 0.5);
                
                
                vector<double> primsNew = prim(newConsVar[0], newConsVar[1], newConsVar[2], newConsVar[3], newConsVar[4]);
                
                density_new[i][j] = primsNew[0];
                vel_x_new[i][j] = primsNew[1];
                vel_y_new[i][j] = primsNew[2];
                vel_z_new[i][j] = primsNew[3];
                pressure_new[i][j] = primsNew[4];
                
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
            }
        }
        
        // y direction (1 step)
        
        density_padded = pad(density);
        pressure_padded = pad(pressure);
        vel_x_padded = pad(vel_x);
        vel_y_padded = pad(vel_y);
        vel_z_padded = pad(vel_z);
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], pressure_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i][j - 1], vel_x_padded[i][j - 1], vel_y_padded[i][j - 1], vel_z_padded[i][j - 1], pressure_padded[i][j - 1]);
                vector<double> consVar_left_left = cons(density_padded[i][j - 2], vel_x_padded[i][j - 2], vel_y_padded[i][j - 2], vel_z_padded[i][j - 2], pressure_padded[i][j - 2]);
                vector<double> consVar_right = cons(density_padded[i][j + 1], vel_x_padded[i][j + 1], vel_y_padded[i][j + 1], vel_z_padded[i][j + 1], pressure_padded[i][j + 1]);
                vector<double> consVar_right_right = cons(density_padded[i][j + 2], vel_x_padded[i][j + 2], vel_y_padded[i][j + 2], vel_z_padded[i][j + 2], pressure_padded[i][j + 2]);
                
                vector<double> newConsVar = updateSLIC_y(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt);
                
                vector<double> primsNew = prim(newConsVar[0], newConsVar[1], newConsVar[2], newConsVar[3], newConsVar[4]);
                
                density_new[i][j] = primsNew[0];
                vel_x_new[i][j] = primsNew[1];
                vel_y_new[i][j] = primsNew[2];
                vel_z_new[i][j] = primsNew[3];
                pressure_new[i][j] = primsNew[4];
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
            }
        }
        
        // x direction (1/2 step)
        
        density_padded = pad(density);
        pressure_padded = pad(pressure);
        vel_x_padded = pad(vel_x);
        vel_y_padded = pad(vel_y);
        vel_z_padded = pad(vel_z);
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            for (int j = 2; j < (cell_number + 2); j++)
            {
                vector<double> consVar = cons(density_padded[i][j], vel_x_padded[i][j], vel_y_padded[i][j], vel_z_padded[i][j], pressure_padded[i][j]);
                vector<double> consVar_left = cons(density_padded[i - 1][j], vel_x_padded[i - 1][j], vel_y_padded[i - 1][j], vel_z_padded[i - 1][j], pressure_padded[i - 1][j]);
                vector<double> consVar_left_left = cons(density_padded[i - 2][j], vel_x_padded[i - 2][j], vel_y_padded[i - 2][j], vel_z_padded[i - 2][j], pressure_padded[i - 2][j]);
                vector<double> consVar_right = cons(density_padded[i + 1][j], vel_x_padded[i + 1][j], vel_y_padded[i + 1][j], vel_z_padded[i + 1][j], pressure_padded[i + 1][j]);
                vector<double> consVar_right_right = cons(density_padded[i + 2][j], vel_x_padded[i + 2][j], vel_y_padded[i + 2][j], vel_z_padded[i + 2][j], pressure_padded[i + 2][j]);
                
                vector<double> newConsVar = updateSLIC_x(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt * 0.5);
                
                vector<double> primsNew = prim(newConsVar[0], newConsVar[1], newConsVar[2], newConsVar[3], newConsVar[4]);
                
                density_new[i][j] = primsNew[0];
                vel_x_new[i][j] = primsNew[1];
                vel_y_new[i][j] = primsNew[2];
                vel_z_new[i][j] = primsNew[3];
                pressure_new[i][j] = primsNew[4];
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
            }
        }
    }
        
    
    ofstream density_file, pressure_file, velocity_file;
    density_file.open("density.dat");
    pressure_file.open("pressure.dat");
    velocity_file.open("velocity.dat");
    
        for (int i = 0; i < cell_number; i++)
        {
            for (int j = 0; j < cell_number; j++)
            {
                density_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << density[i][j] << "\n";
                pressure_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << pressure[i][j] << "\n";
                velocity_file << i * (domain_size / cell_number) << " " << j * (domain_size / cell_number) << " " << vel_x[i][j] << "\n";
            
            }
            density_file << "\n";
            pressure_file << "\n";
            velocity_file << "\n";
        }
        
        density_file.close();
        pressure_file.close();
        velocity_file.close();
    
    return 0;
}
