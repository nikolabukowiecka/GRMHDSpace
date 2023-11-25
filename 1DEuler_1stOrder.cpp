#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

double adiabatic_constant = 1.4;
double cfl_constant = 0.9;
double domain_size = 1.0;
int cell_number = 800;
double final_time = 0.15;

vector<double> flux(vector<double> consVar)
{
    double density = consVar[0];
    double vel_x = consVar[1] / density;
    double vel_y = consVar[2] / density;
    double vel_z = consVar[3] / density;
    double energy = consVar[4];
    double pressure = (adiabatic_constant - 1.0) * (energy - 0.5 * density * ((vel_x * vel_x) + (vel_y * vel_y) + (vel_z * vel_z)));
    
    vector<double> ret(5);
    
    ret[0] = density * vel_x;
    ret[1] = density * (vel_x * vel_x) + pressure;
    ret[2] = density * (vel_x * vel_y);
    ret[3] = density * (vel_x * vel_z);
    ret[4] = vel_x * (energy + pressure);
    
    return ret;
}

vector<double> cons(double density, double vel_x, double vel_y, double vel_z, double energy)
{
    vector<double> ret(5);
    
    ret[0] = density;
    ret[1] = density * vel_x;
    ret[2] = density * vel_y;
    ret[3] = density * vel_z;
    ret[4] = energy;
    
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

vector<double> multiplyVector(double scalar, vector<double> vector1)
{
    vector<double> ret(vector1.size());
    
    for (int i = 0; i < vector1.size(); i++)
    {
        ret[i] = scalar * vector1[i];
    }
    
    return ret;
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

double time_step(vector<double> pressure, vector<double> density, vector<double> vel_x)
{
    double u_max = 0.0;

    for (int i = 0; i < cell_number; i++)
    {
        double a = sqrt(adiabatic_constant * pressure[i] / density[i]);
        double speed = abs(vel_x[i]) + abs(a);
        
        if (speed > u_max)
        {
            u_max = speed;
        }
    }
    
    return cfl_constant * (domain_size / cell_number) / u_max;
}

vector<double> pad(vector<double> quantity)
{
    vector<double> paddedQuantity(cell_number + 2);
    paddedQuantity[0] = quantity[0];
    
    for (int i = 0; i < cell_number; i++)
    {
        paddedQuantity[i + 1] = quantity[i];
    }
    
    paddedQuantity[cell_number + 1] = quantity[cell_number - 1];
    
    return paddedQuantity;
}

int main()
{
    vector<double> density(cell_number);
    vector<double> pressure(cell_number);
    vector<double> vel_x(cell_number);
    vector<double> vel_y(cell_number);
    vector<double> vel_z(cell_number);
    
    double t = 0.0;
    int n = 1;
    
    for (int i = 0; i <= (cell_number / 2); i++)
    {
        density[i] = 1.0;
        vel_x[i] = -2.0;
        vel_y[i] = 0.0;
        vel_z[i] = 0.0;
        pressure[i] = 0.4;
    }
    for (int i = ((cell_number / 2) + 1); i < cell_number; i++)
    {
        density[i] = 1.0;
        vel_x[i] = 2.0;
        vel_y[i] = 0.0;
        vel_z[i] = 0.0;
        pressure[i] = 0.4;
    }
    
    while (t < final_time)
    {
        double dt = time_step(pressure, density, vel_x);
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
        
        vector<double> density_new(cell_number + 2);
        vector<double> pressure_new(cell_number + 2);
        vector<double> vel_x_new(cell_number + 2);
        vector<double> vel_y_new(cell_number + 2);
        vector<double> vel_z_new(cell_number + 2);
        
        vector<double> energy_tot(cell_number + 2);
        vector<double> energy_new(cell_number + 2);
        for (int i = 0; i < cell_number + 2; i++)
        {
            double energy_int = pressure_padded[i] / (density_padded[i] * (adiabatic_constant - 1.0));
            
            energy_tot[i] = (0.5 * density_padded[i] * ((vel_x_padded[i] * vel_x_padded[i]) + (vel_y_padded[i] * vel_y_padded[i]) + (vel_z_padded[i] * vel_z_padded[i]))) + (density_padded[i] * energy_int);
        }
        
        for (int i = 1; i < (cell_number + 1); i++)
        {
            vector<double> consVar = cons(density_padded[i], vel_x_padded[i], vel_y_padded[i], vel_z_padded[i], energy_tot[i]);
            vector<double> consVar_left = cons(density_padded[i - 1], vel_x_padded[i - 1], vel_y_padded[i - 1], vel_z_padded[i - 1], energy_tot[i - 1]);
            vector<double> consVar_right = cons(density_padded[i + 1], vel_x_padded[i + 1], vel_y_padded[i + 1], vel_z_padded[i + 1], energy_tot[i + 1]);
            
            vector<double> newConsVar = update(consVar, consVar_left, consVar_right, flux(consVar), flux(consVar_left), flux(consVar_right), dt);
            
            density_new[i] = newConsVar[0];
            vel_x_new[i] = newConsVar[1] / density_new[i];
            vel_y_new[i] = newConsVar[2] / density_new[i];
            vel_z_new[i] = newConsVar[3] / density_new[i];
            energy_new[i] = newConsVar[4];
            pressure_new[i] = (adiabatic_constant - 1.0) * (energy_new[i] - (0.5 * density_new[i] * ((vel_x_new[i] * vel_x_new[i]) + (vel_y_new[i] * vel_y_new[i]) + (vel_z_new[i] * vel_z_new[i]))));
        }
        
        for (int i = 1; i < (cell_number + 1); i++)
        {
            density[i - 1] = density_new[i];
            pressure[i - 1] = pressure_new[i];
            vel_x[i - 1] = vel_x_new[i];
            vel_y[i - 1] = vel_y_new[i];
            vel_z[i - 1] = vel_z_new[i];
        }
    }
    
    ofstream density_file, pressure_file, velocity_file, energy_int_file;
    density_file.open("density.dat");
    pressure_file.open("pressure.dat");
    velocity_file.open("velocity.dat");
    energy_int_file.open("energy.dat");
    
    for (int i = 0; i < cell_number; i++)
    {
        density_file << i * (domain_size / cell_number) << " " << density[i] << "\n";
        pressure_file << i * (domain_size / cell_number) << " " << pressure[i] << "\n";
        velocity_file << i * (domain_size / cell_number) << " " << vel_x[i] << "\n";
        energy_int_file << i * (domain_size / cell_number) << " " << pressure[i] / (density[i] * (adiabatic_constant - 1.0)) << "\n";
    }
    
    density_file.close();
    pressure_file.close();
    velocity_file.close();
    energy_int_file.close();
    
    return 0;
}
