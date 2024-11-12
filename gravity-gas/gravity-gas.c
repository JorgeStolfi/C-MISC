/* Last edited on 2018-05-29 16:40:40 by jstolfi */
/* Simulation of gas atmosphere in a star. */

typedef struct params_t
  { int32_t np;    // Number of particles.
    double dtMin;  // Min time step.
    double dtMax;  // Max time step.
    double rSyst;  // Radius of region of interest.
    // Star data:
    double rStar;  // Radius of nominal star surface.
    double QStar;  // Total charge of star.
    // Particle data:
    double QPart;  // Charge unit of particles.
    double rPart;  // Radius of particles.
    double EPart;  // Mean energy of particle at star surface.
  } params_t;

int main (int argc, char **argv)
  { 
    params_t *par = get_parameters();
    double t = 0;
    for (int32_t kt = 0; kt < niter; kt++)
      { /* Simulate another time step: */
        compute_accelerations(np, chr, pos, acc);
        update_velocities_and_positions(np, pos, vel, acc);
        t = t + dt;
      }
    return 0;

  }
  
params_t *get_parameters(void)
  { params_t *par = malloc(sizeof(params_t));
    par->np = 100;
    par->dtMax = 1.0e-3;
    par->dtMin = par->dtMax/1000;
    par->rStar = 1.0;
    par->QStar = 1.0e+3;
    par->QPart = 1.0e-3;
    par->rPart = 0.01;
    par->vPart = 
    par->rSyst = 10* par->rStar
    return par;
  }
  
void compute_accelerations(int32_t np, double chr[], double pos[], double acc[])
  { 
    for (int32_t kp = 0; kp < np; kp++)
      { /* Compute acceleration of particle {kp}:  */
        acc[kp] = compute_particle_accel(np, chr[kp], pos[kp]);
      }
  }

void detect_next_collision
  ( int32_t np, 
    double pos[], 
    double vel[], 
    double acc[], 
    double *dt, 
    int32_t *ipP, 
    int32_t *jpP
  ) 
  {
    int32_t iNext = -1;
    int32_t jNext = -1;
    double dtNext = +INF;
    for (int32_t ip = 0; ip < np; ip++)
      { for (int32_t jp = i+1; jp < np; jp++)
          { dtij = detect_two_particle_collision
              ( pos[ip], vel[ip], acc[ip],
                pos[jp], vel[jp], acc[jp]
        
  }
  
void update_velocities_and_positons()
  { 
    detect_next_collision(np, pos, vel, acc, dt, i, j);

    for (int32_t kp = 0; kp < np; kp++)
      { 
        /* Simulate motion with time step: */
        vel[kp] = vel[kp] + dt * acc[kp];
        pos[kp] = pos[kp] + dt * vel[kp] + 0.5*dt*dt*acc[kp];
        int32_t jp = col[kp];
        if (jp >= 0)
          { /* Simulate collision with particle {jp}: */
            vel[kp] = simul_collision(kp, jp);
          }
     }
  }
  
void update_velocities()
  { for (int32_t kp = 0; kp < np; kp++)
      { int32_t jp = col[kp];
        if (jp >= 0)
          { /* Simulate collision with particle {jp}: */
            vel[kp] = simul_collision(kp, jp);
          }
        /* Simulate motion with time step: */
        vel[kp] = vel[kp] + dt * acc[kp];
      }
  }

        
