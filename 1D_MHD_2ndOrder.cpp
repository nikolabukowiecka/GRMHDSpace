#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "/opt/homebrew/Cellar/libomp/17.0.5/include/omp.h"

using namespace std;

double adiabatic_constant = 5.0 / 3.0;
double cfl_constant = 0.9;
double domain_size = 1.0;
int cell_number = 800;
double final_time = 0.2;
double bias = 0.0;

vector<double> flux(vector<double> consVar)
{
    double density = consVar[0];
    double vel_x = consVar[1] / density;
    double vel_y = consVar[2] / density;
    double vel_z = consVar[3] / density;
    double magnetic_x = consVar[4];
    double magnetic_y = consVar[5];
    double magnetic_z = consVar[6];
    double energy = consVar[7];
    double pressure = (adiabatic_constant - 1.0) * (energy - (0.5 * density * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z))) - (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z))));
    double pressure_star = pressure + (0.5 * ((magnetic_x * magnetic_x) + (magnetic_y * magnetic_y) + (magnetic_z * magnetic_z)));
    
    vector<double> ret(8);
    
    ret[0] = density * vel_x;
    ret[1] = density * (vel_x * vel_x) + pressure_star - (magnetic_x * magnetic_x);
    ret[2] = density * (vel_x * vel_y) - (magnetic_x * magnetic_y);
    ret[3] = density * (vel_x * vel_z) - (magnetic_x * magnetic_z);
    ret[4] = (magnetic_x * vel_x) - (magnetic_x * vel_x);
    ret[5] = (magnetic_y * vel_x) - (magnetic_x * vel_y);
    ret[6] = (magnetic_z * vel_x) - (magnetic_x * vel_z);
    ret[7] = ((energy + pressure_star) * vel_x) - (magnetic_x * ((magnetic_x * vel_x) + (magnetic_y * vel_y) + (magnetic_z * vel_z)));
    
    return ret;
}

vector<double> cons(double density, double vel_x, double vel_y, double vel_z, double magnetic_x, double magnetic_y, double magnetic_z, double energy)
{
    vector<double> ret(8);
    
    ret[0] = density;
    ret[1] = density * vel_x;
    ret[2] = density * vel_y;
    ret[3] = density * vel_z;
    ret[4] = magnetic_x;
    ret[5] = magnetic_y;
    ret[6] = magnetic_z;
    ret[7] = energy;
    
    return ret;
}

vector<double> addVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
#pragma omp parallel for
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] + vector2[i];
    }
    
    return ret;
}

vector<double> subtractVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
#pragma omp parallel for
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] - vector2[i];
    }
    
    return ret;
}

vector<double> multiplyVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
#pragma omp parallel for
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] * vector2[i];
    }
    
    return ret;
}

vector<double> divideVectors(vector<double> vector1, vector<double> vector2)
{
    vector<double> ret(vector1.size());
    
#pragma omp parallel for
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = vector1[i] / vector2[i];
    }
    
    return ret;
}

vector<double> multiplyVector(double scalar, vector<double> vector1)
{
    vector<double> ret(vector1.size());
    
#pragma omp parallel for
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = scalar * vector1[i];
    }
    
    return ret;
}

vector<double> consVar_left_halfEvolve(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_left, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux(consVar_left), flux(consVar_right))));
}

vector<double> consVar_right_halfEvolve(vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return addVectors(consVar_right, multiplyVector(0.5 * (dt / (domain_size / cell_number)), subtractVectors(flux(consVar_left), flux(consVar_right))));
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
    
#pragma omp parallel for
    for (int i = 0; i < xiR.size(); i++)
    {
        xiR[i] = 2.0 / (1.0 - bias + (1.0 + bias) * r[i]);
    }
    
    vector<double> xi(r.size());
    
#pragma omp parallel for
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

vector<double> force(vector<double> flux_left, vector<double> flux_right, vector<double> consVar_left, vector<double> consVar_right, double dt)
{
    return multiplyVector(0.5, addVectors(force_LF(flux_left, flux_right, consVar_left, consVar_right, dt), flux(force_RI(flux_left, flux_right, consVar_left, consVar_right, dt))));
}

vector<double> update(vector<double> consVar, vector<double> consVar_left, vector<double> consVar_right, vector<double> fluxVect, vector<double> fluxVect_left, vector<double> fluxVect_right, double dt)
{
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)), subtractVectors(force(fluxVect_left, fluxVect, consVar_left, consVar, dt), force(fluxVect, fluxVect_right, consVar, consVar_right, dt))));
}

vector<double> updateSLIC(vector<double> consVar_left_left, vector<double> consVar_left, vector<double> consVar, vector<double> consVar_right, vector<double> consVar_right_right, double dt)
{
    vector<double> consVar_left_boundaryL = extrapolate_boundary_left(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryR = extrapolate_boundary_right(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_boundaryL_evolved = consVar_left_halfEvolve(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    vector<double> consVar_left_boundaryR_evolved = consVar_right_halfEvolve(consVar_left_boundaryL, consVar_left_boundaryR, dt);
    
    vector<double> consVar_boundaryL = extrapolate_boundary_left(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryR = extrapolate_boundary_right(consVar_left, consVar, consVar_right);
    vector<double> consVar_boundaryL_evolved = consVar_left_halfEvolve(consVar_boundaryL, consVar_boundaryR, dt);
    vector<double> consVar_boundaryR_evolved = consVar_right_halfEvolve(consVar_boundaryL, consVar_boundaryR, dt);
    
    vector<double> consVar_right_boundaryL = extrapolate_boundary_left(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryR = extrapolate_boundary_right(consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_boundaryL_evolved = consVar_left_halfEvolve(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    vector<double> consVar_right_boundaryR_evolved = consVar_right_halfEvolve(consVar_right_boundaryL, consVar_right_boundaryR, dt);
    
    return addVectors(consVar, multiplyVector((dt / (domain_size / cell_number)),
                                              subtractVectors(
                                                              force(flux(consVar_left_boundaryR_evolved), flux(consVar_boundaryL_evolved), consVar_left_boundaryR_evolved, consVar_boundaryL_evolved, dt),
                                                              force(flux(consVar_boundaryR_evolved), flux(consVar_right_boundaryL_evolved), consVar_boundaryR_evolved, consVar_right_boundaryL_evolved, dt))));
    /*
    vector<double> consVar_left_boundary = extrapolate_boundary_left(consVar_left_left, consVar_left, consVar);
    vector<double> consVar_left_next_boundary = extrapolate_boundary_left(consVar_left, consVar, consVar_right);
    vector<double> consVar_right_boundary = extrapolate_boundary_right(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right);
    vector<double> consVar_right_next_boundary = extrapolate_boundary_right(consVar_left, consVar, consVar_right, consVar_right_right, consVar_right_right_right);
    
    vector<double> consVar_left_half = consVar_left_halfEvolve(consVar_left_boundary, consVar_right_boundary, dt);
    vector<double> consVar_left_next_half = consVar_left_halfEvolve(consVar_left_next_boundary, consVar_right_next_boundary, dt);
    vector<double> consVar_right_half = consVar_right_halfEvolve(consVar_left_boundary, consVar_right_boundary, dt);
    */
}

double time_step(vector<double> pressure, vector<double> density, vector<double> vel_x, vector<double> magnetic_x, vector<double> magnetic_y, vector<double> magnetic_z)
{
    double u_max = 0.0;

#pragma omp parallel for
    for (int i = 0; i < cell_number; i++)
    {
        double expr = ((adiabatic_constant * pressure[i]) + ((magnetic_x[i] * magnetic_x[i]) + (magnetic_y[i] * magnetic_y[i]) + (magnetic_z[i] * magnetic_z[i]))) / density[i];
        double c_f = sqrt(
                          0.5 * (expr + sqrt((expr * expr) - (4 * (adiabatic_constant * pressure[i] / density[i]) * (magnetic_x[i] * magnetic_x[i] / density[i])))));
        double speed = abs(vel_x[i]) + abs(c_f);
        
        if (speed > u_max)
        {
            u_max = speed;
        }
    }
    
    return cfl_constant * (domain_size / cell_number) / u_max;
}

vector<double> pad(vector<double> quantity)
{
    vector<double> paddedQuantity(cell_number + 4);
    paddedQuantity[0] = quantity[0];
    paddedQuantity[1] = quantity[0];
    
#pragma omp parallel for
    for (int i = 0; i < cell_number; i++)
    {
        paddedQuantity[i + 2] = quantity[i];
    }
    
    paddedQuantity[cell_number + 2] = quantity[cell_number - 1];
    paddedQuantity[cell_number + 3] = quantity[cell_number - 1];
    
    return paddedQuantity;
}

int main()
{
    vector<double> density(cell_number);
    vector<double> pressure(cell_number);
    vector<double> vel_x(cell_number);
    vector<double> vel_y(cell_number);
    vector<double> vel_z(cell_number);
    vector<double> magnetic_x(cell_number);
    vector<double> magnetic_y(cell_number);
    vector<double> magnetic_z(cell_number);
    
    double t = 0.0;
    int n = 1;
    
    for (int i = 0; i <= (cell_number * 0.5); i++)
    {
        density[i] = 1.08;
        vel_x[i] = 1.2;
        vel_y[i] = 0.01;
        vel_z[i] = 0.5;
        magnetic_x[i] = 2.0 / sqrt(4.0 * M_PI);
        magnetic_y[i] = 3.6 / sqrt(4.0 * M_PI);
        magnetic_z[i] = 2.0 / sqrt(4.0 * M_PI);
        pressure[i] = 0.95;
    }
#pragma omp parallel for
    for (int i = ((cell_number * 0.5) + 1); i < cell_number; i++)
    {
        density[i] = 1.0;
        vel_x[i] = 0.0;
        vel_y[i] = 0.0;
        vel_z[i] = 0.0;
        magnetic_x[i] = 2.0 / sqrt(4.0 * M_PI);
        magnetic_y[i] = 4.0 / sqrt(4.0 * M_PI);
        magnetic_z[i] = 2.0 / sqrt(4.0 * M_PI);
        pressure[i] = 1.0;
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
        
        vector<double> density_padded = pad(density);
        vector<double> pressure_padded = pad(pressure);
        vector<double> vel_x_padded = pad(vel_x);
        vector<double> vel_y_padded = pad(vel_y);
        vector<double> vel_z_padded = pad(vel_z);
        vector<double> magnetic_x_padded = pad(magnetic_x);
        vector<double> magnetic_y_padded = pad(magnetic_y);
        vector<double> magnetic_z_padded = pad(magnetic_z);
        
        vector<double> density_new(cell_number + 4);
        vector<double> pressure_new(cell_number + 4);
        vector<double> vel_x_new(cell_number + 4);
        vector<double> vel_y_new(cell_number + 4);
        vector<double> vel_z_new(cell_number + 4);
        vector<double> magnetic_x_new(cell_number + 4);
        vector<double> magnetic_y_new(cell_number + 4);
        vector<double> magnetic_z_new(cell_number + 4);
        
        vector<double> energy_tot(cell_number + 4);
        vector<double> energy_new(cell_number + 4);
        
#pragma omp parallel for
        for (int i = 0; i < cell_number + 4; i++)
        {
            double energy_int = pressure_padded[i] / (density_padded[i] * (adiabatic_constant - 1.0));
            
            energy_tot[i] = (0.5 * density_padded[i] * ((vel_x_padded[i] * vel_x_padded[i]) + (vel_y_padded[i] * vel_y_padded[i]) + (vel_z_padded[i] * vel_z_padded[i]))) + (0.5 * ((magnetic_x_padded[i] * magnetic_x_padded[i]) + (magnetic_y_padded[i] * magnetic_y_padded[i]) + (magnetic_z_padded[i] * magnetic_z_padded[i]))) + (density_padded[i] * energy_int);
        }
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            vector<double> consVar = cons(density_padded[i], vel_x_padded[i], vel_y_padded[i], vel_z_padded[i], magnetic_x_padded[i], magnetic_y_padded[i], magnetic_z_padded[i], energy_tot[i]);
            vector<double> consVar_left = cons(density_padded[i - 1], vel_x_padded[i - 1], vel_y_padded[i - 1], vel_z_padded[i - 1], magnetic_x_padded[i - 1], magnetic_y_padded[i - 1], magnetic_z_padded[i - 1], energy_tot[i - 1]);
            vector<double> consVar_left_left = cons(density_padded[i - 2], vel_x_padded[i - 2], vel_y_padded[i - 2], vel_z_padded[i - 2], magnetic_x_padded[i - 2], magnetic_y_padded[i - 2], magnetic_z_padded[i - 2], energy_tot[i - 2]);
            vector<double> consVar_right = cons(density_padded[i + 1], vel_x_padded[i + 1], vel_y_padded[i + 1], vel_z_padded[i + 1], magnetic_x_padded[i + 1], magnetic_y_padded[i + 1], magnetic_z_padded[i + 1], energy_tot[i + 1]);
            vector<double> consVar_right_right = cons(density_padded[i + 2], vel_x_padded[i + 2], vel_y_padded[i + 2], vel_z_padded[i + 2], magnetic_x_padded[i + 2], magnetic_y_padded[i + 2], magnetic_z_padded[i + 2], energy_tot[i + 2]);
            
            vector<double> newConsVar = updateSLIC(consVar_left_left, consVar_left, consVar, consVar_right, consVar_right_right, dt);
            
            density_new[i] = newConsVar[0];
            vel_x_new[i] = newConsVar[1] / density_new[i];
            vel_y_new[i] = newConsVar[2] / density_new[i];
            vel_z_new[i] = newConsVar[3] / density_new[i];
            magnetic_x_new[i] = newConsVar[4];
            magnetic_y_new[i] = newConsVar[5];
            magnetic_z_new[i] = newConsVar[6];
            energy_new[i] = newConsVar[7];
            pressure_new[i] = (adiabatic_constant - 1.0) * (energy_new[i] - (0.5 * density_new[i] * ((vel_x_new[i] * vel_x_new[i]) + (vel_y_new[i] * vel_y_new[i]) + (vel_z_new[i] * vel_z_new[i]))) - (0.5 * ((magnetic_x_new[i] * magnetic_x_new[i]) + (magnetic_y_new[i] * magnetic_y_new[i]) + (magnetic_z_new[i] * magnetic_z_new[i]))));
        }
        
#pragma omp parallel for
        for (int i = 2; i < (cell_number + 2); i++)
        {
            density[i - 2] = density_new[i];
            pressure[i - 2] = pressure_new[i];
            vel_x[i - 2] = vel_x_new[i];
            vel_y[i - 2] = vel_y_new[i];
            vel_z[i - 2] = vel_z_new[i];
            magnetic_x[i - 2] = magnetic_x_new[i];
            magnetic_y[i - 2] = magnetic_y_new[i];
            magnetic_z[i - 2] = magnetic_z_new[i];
        }
    }
    
    ofstream density_file, pressure_file, velocity_x_file, velocity_y_file, magnetic_y_file, energy_int_file;
    density_file.open("density.dat");
    pressure_file.open("pressure.dat");
    velocity_x_file.open("velocity_x.dat");
    velocity_y_file.open("velocity_y.dat");
    magnetic_y_file.open("magnetic_y.dat");
    energy_int_file.open("energy.dat");
    
    for (int i = 0; i < cell_number; i++)
    {
        density_file << i * (domain_size / cell_number) << " " << density[i] << "\n";
        pressure_file << i * (domain_size / cell_number) << " " << pressure[i] << "\n";
        velocity_x_file << i * (domain_size / cell_number) << " " << vel_x[i] << "\n";
        velocity_y_file << i * (domain_size / cell_number) << " " << vel_y[i] << "\n";
        magnetic_y_file << i * (domain_size / cell_number) << " " << magnetic_y[i] << "\n";
        energy_int_file << i * (domain_size / cell_number) << " " << pressure[i] / (density[i] * (adiabatic_constant - 1.0)) << "\n";
    }
    
    density_file.close();
    pressure_file.close();
    velocity_x_file.close();
    velocity_y_file.close();
    magnetic_y_file.close();
    energy_int_file.close();
    
    return 0;
}
