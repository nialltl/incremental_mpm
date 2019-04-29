using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;
using UnityEngine;
using UnityEngine.Jobs;
using System.Collections.Generic;

using UnityEngine.Profiling;
using Unity.Collections.LowLevel.Unsafe;
using System.Runtime.InteropServices;

public class MLS_MPM_Fluid_Multithreaded : MonoBehaviour {
    struct Particle {
        public float2 x; // position
        public float2 v; // velocity
        public float2x2 C; // affine momentum matrix
        public float mass;
        public float padding; // unused
    }

    struct Cell {
        public float2 v; // velocity
        public float mass;
        public float padding; // unused
    }
    
    const int grid_res = 64;
    const int num_cells = grid_res * grid_res;
    const int division = 128; // job system batch size, chosen experimentally
    
    // simulation
    const float dt = 0.2f;
    const float sim_iterations = (int)(1.0f / dt);
    
    const float gravity = -0.3f;
    // fluid parameters
    const float rest_density = 4.0f;
    const float dynamic_viscosity = 0.1f;
    // equation of state
    const float eos_stiffness = 10.0f;
    const float eos_power = 4;
    
    NativeArray<Particle> ps;
    NativeArray<Cell> grid;

    float2[] weights = new float2[3];

    int num_particles;
    
    SimRenderer sim_renderer;
    
    // interaction
    const float mouse_radius = 10;
    bool mouse_down = false;
    float2 mouse_pos;
    
    void Start () {
        // populate our array of particles from the samples given, set their initial state
        List<float2> temp_positions = new List<float2>();
        const float spacing = 0.5f;
        const int box_x = 32, box_y = 32;
        const float sx = grid_res / 2.0f, sy = grid_res / 2.0f;
        for (float i = sx - box_x / 2; i < sx + box_x / 2; i += spacing) {
            for (float j = sy - box_y / 2; j < sy + box_y / 2; j += spacing) {
                var pos = math.float2(i, j);

                temp_positions.Add(pos);
            }
        }

        // round number of particles down to the nearest power of 2. taken from 
        int po2_amnt = 1; while (po2_amnt <= temp_positions.Count) po2_amnt <<= 1;
        num_particles = po2_amnt >> 1;

        ps = new NativeArray<Particle>(num_particles, Allocator.Persistent);

        for (int i = 0; i < num_particles; ++i) {
            Particle p = new Particle();
            p.x = temp_positions[i];
            p.v = 0;
            p.C = 0;
            p.mass = 1.0f;
            ps[i] = p;
        }

        grid = new NativeArray<Cell>(num_cells, Allocator.Persistent);

        for (int i = 0; i < num_cells; ++i) {
            var cell = new Cell();
            cell.v = 0;
            grid[i] = cell;
        }

        sim_renderer = GameObject.FindObjectOfType<SimRenderer>();
        sim_renderer.Initialise(num_particles, Marshal.SizeOf(new Particle()));
    }

    private void Update() {
        HandleMouseInteraction();

        for (int i = 0; i < sim_iterations; ++i) {
            Simulate();
        }

        sim_renderer.RenderFrame(ps);
    }
    
    void HandleMouseInteraction() {
        mouse_down = false;
        if (Input.GetMouseButton(0)) {
            mouse_down = true;
            var mp = Camera.main.ScreenToViewportPoint(Input.mousePosition);
            mouse_pos = math.float2(mp.x * grid_res, mp.y * grid_res);
        }
    }
    
    void Simulate() {
        Profiler.BeginSample("ClearGrid");
        new Job_ClearGrid() {
            grid = grid
        }.Schedule(num_cells, division).Complete();
        Profiler.EndSample();

        // P2G, first round
        Profiler.BeginSample("P2G 1");
        new Job_P2G_1() {
            ps = ps,
            grid = grid,
            num_particles = num_particles
        }.Schedule().Complete();
        Profiler.EndSample();
        
        // P2G, second round      
        Profiler.BeginSample("P2G 2");  
        new Job_P2G_2() {
            ps = ps,
            grid = grid,
            num_particles = num_particles
        }.Schedule().Complete();
        Profiler.EndSample();
        
        Profiler.BeginSample("Update grid");
        new Job_UpdateGrid() {
            grid = grid
        }.Schedule(num_cells, division).Complete();
        Profiler.EndSample();
        
        Profiler.BeginSample("G2P");
        new Job_G2P() {
            ps = ps,
            grid = grid,
            mouse_down = mouse_down,
            mouse_pos = mouse_pos
        }.Schedule(num_particles, division).Complete();
        Profiler.EndSample();
    }

    #region Jobs
    [BurstCompile]
    struct Job_ClearGrid : IJobParallelFor {
        public NativeArray<Cell> grid;

        public void Execute(int i) {
            var cell = grid[i];

            // reset grid scratch-pad entirely
            cell.mass = 0;
            cell.v = 0;

            grid[i] = cell;
        }
    }
    
    [BurstCompile]
    unsafe struct Job_P2G_1 : IJob {
        public NativeArray<Particle> ps;
        [ReadOnly] public NativeArray<Cell> grid;
        [ReadOnly] public int num_particles;
        
        public void Execute() {
            var weights = stackalloc float2[3];

            for (int i = 0; i < num_particles; ++i) {
                var p = ps[i];

                uint2 cell_idx = (uint2)p.x;
                float2 cell_diff = (p.x - cell_idx) - 0.5f;
                weights[0] = 0.5f * math.pow(0.5f - cell_diff, 2);
                weights[1] = 0.75f - math.pow(cell_diff, 2);
                weights[2] = 0.5f * math.pow(0.5f + cell_diff, 2);

                float2x2 C = p.C;
                
                for (uint gx = 0; gx < 3; ++gx) {
                    for (uint gy = 0; gy < 3; ++gy) {
                        float weight = weights[gx].x * weights[gy].y;

                        uint2 cell_x = math.uint2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                        float2 cell_dist = (cell_x - p.x) + 0.5f;
                        float2 Q = math.mul(C, cell_dist);

                        // MPM course, equation 172
                        float mass_contrib = weight * p.mass;

                        int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                        Cell cell = grid[cell_index];

                        // mass and momentum update
                        cell.mass += mass_contrib;
                        cell.v += mass_contrib * (p.v + Q);
                        
                        grid[cell_index] = cell;
                    }
                }
            }
        }
    }

    [BurstCompile]
    unsafe struct Job_P2G_2 : IJob {
        public NativeArray<Cell> grid;
        [ReadOnly] public NativeArray<Particle> ps;
        [ReadOnly] public int num_particles;

        // we now have 2 P2G phases as we need to ensure we have scattered particle masses to the grid,
        // in order to get our density estimate at each frame
        
        public void Execute() {
            var weights = stackalloc float2[3];

            for (int i = 0; i < num_particles; ++i) {
                var p = ps[i];

                uint2 cell_idx = (uint2)p.x;
                float2 cell_diff = (p.x - cell_idx) - 0.5f;
                weights[0] = 0.5f * math.pow(0.5f - cell_diff, 2);
                weights[1] = 0.75f - math.pow(cell_diff, 2);
                weights[2] = 0.5f * math.pow(0.5f + cell_diff, 2);

                // estimating particle volume by summing up neighbourhood's weighted mass contribution
                // MPM course, equation 152 
                float density = 0.0f;
                uint gx, gy;
                for (gx = 0; gx < 3; ++gx) {
                    for (gy = 0; gy < 3; ++gy) {
                        float weight = weights[gx].x * weights[gy].y;
                        int cell_index = (int)(cell_idx.x + gx - 1) * grid_res + (int)(cell_idx.y + gy - 1);
                        density += grid[cell_index].mass * weight;
                    }
                }

                float volume = p.mass / density;

                // end goal, constitutive equation for isotropic fluid: 
                // stress = -pressure * I + viscosity * (velocity_gradient + velocity_gradient_transposed)

                // Tait equation of state. i clamped it as a bit of a hack.
                // clamping helps prevent particles absorbing into each other with negative pressures
                float pressure = math.max(-0.1f, eos_stiffness * (math.pow(density / rest_density, eos_power) - 1));
                
                float2x2 stress = math.float2x2(
                    -pressure, 0, 
                    0, -pressure
                );

                // velocity gradient - CPIC eq. 17, where deriv of quadratic polynomial is linear
                float2x2 dudv = p.C;
                float2x2 strain = dudv;
                
                float trace = strain.c1.x + strain.c0.y;
                strain.c0.y = strain.c1.x = trace;

                float2x2 viscosity_term = dynamic_viscosity * strain;
                stress += viscosity_term;
                
                var eq_16_term_0 = -volume * 4 * stress * dt;

                for (gx = 0; gx < 3; ++gx) {
                    for (gy = 0; gy < 3; ++gy) {
                        float weight = weights[gx].x * weights[gy].y;

                        uint2 cell_x = math.uint2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                        float2 cell_dist = (cell_x - p.x) + 0.5f;

                        int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                        Cell cell = grid[cell_index];
                        
                        // fused force + momentum contribution from MLS-MPM
                        float2 momentum = math.mul(eq_16_term_0 * weight, cell_dist);
                        cell.v += momentum;

                        grid[cell_index] = cell;
                    }
                }
            }
        }
    }
    
    [BurstCompile]
    struct Job_UpdateGrid : IJobParallelFor {
        public NativeArray<Cell> grid;

        public void Execute(int i) {
            var cell = grid[i];

            if (cell.mass > 0) {
                // convert momentum to velocity, apply gravity
                cell.v /= cell.mass;
                cell.v += dt * math.float2(0, gravity);

                // boundary conditions
                int x = i / grid_res;
                int y = i % grid_res;
                if (x < 2 || x > grid_res - 3) { cell.v.x = 0; }
                if (y < 2 || y > grid_res - 3) { cell.v.y = 0; }

                grid[i] = cell;
            }
        }
    }


    [BurstCompile]
    unsafe struct Job_G2P : IJobParallelFor {
        public NativeArray<Particle> ps;
        [ReadOnly] public NativeArray<Cell> grid;
        [ReadOnly] public bool mouse_down;
        [ReadOnly] public float2 mouse_pos;
        
        public void Execute(int i) {
            Particle p = ps[i];

            // reset particle velocity. we calculate it from scratch each step using the grid
            p.v = 0;

            uint2 cell_idx = (uint2)p.x;
            float2 cell_diff = (p.x - cell_idx) - 0.5f;

            var weights = stackalloc float2[] {
                0.5f * math.pow(0.5f - cell_diff, 2),
                0.75f - math.pow(cell_diff, 2), 
                0.5f * math.pow(0.5f + cell_diff, 2)
            };
            
            // constructing affine per-particle momentum matrix from APIC / MLS-MPM.
            // see APIC paper (https://web.archive.org/web/20190427165435/https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
            // below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
            // where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
            float2x2 B = 0;
            for (uint gx = 0; gx < 3; ++gx) {
                for (uint gy = 0; gy < 3; ++gy) {
                    float weight = weights[gx].x * weights[gy].y;

                    uint2 cell_x = math.uint2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                    int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                    
                    float2 dist = (cell_x - p.x) + 0.5f;
                    float2 weighted_velocity = grid[cell_index].v * weight;

                    var term = math.float2x2(weighted_velocity * dist.x, weighted_velocity * dist.y);

                    B += term;

                    p.v += weighted_velocity;
                }
            }

            p.C = B * 4;
            
            // advect particles
            p.x += p.v * dt;

            // safety clamp to ensure particles don't exit simulation domain
            p.x = math.clamp(p.x, 1, grid_res - 2);

            if (mouse_down) {
                var dist = p.x - mouse_pos;
                if (math.dot(dist, dist) < mouse_radius * mouse_radius) {
                    var force = math.normalizesafe(dist, 0) * 1.0f;
                    p.v += force;
                }
            }

            // boundaries
            float2 x_n = p.x + p.v;
            const float wall_min = 3;
            float wall_max = (float)grid_res - 4;
            if (x_n.x < wall_min) p.v.x += wall_min - x_n.x;
            if (x_n.x > wall_max) p.v.x += wall_max - x_n.x;
            if (x_n.y < wall_min) p.v.y += wall_min - x_n.y;
            if (x_n.y > wall_max) p.v.y += wall_max - x_n.y;

            // no need for the deformation gradient update here,
            // as we never use it in our constitutive equation

            ps[i] = p;
        }
    }
    #endregion

    void OnDestroy() {
        ps.Dispose();
        grid.Dispose();
    }
}

