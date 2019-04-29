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

public class MLS_MPM_NeoHookean_Multithreaded : MonoBehaviour {
    struct Particle {
        public float2 x; // position
        public float2 v; // velocity
        public float2x2 C; // affine momentum matrix
        public float mass;
        public float volume_0; // initial volume
    }

    struct Cell {
        public float2 v; // velocity
        public float mass;
        public float padding; // unused
    }
    
    const int grid_res = 64;
    const int num_cells = grid_res * grid_res;

    // batch size for the job system. just determined experimentally
    const int division = 16;
    
    // simulation parameters
    const float dt = 0.1f; // timestep
    const float iterations = (int)(1.0f / dt);
    const float gravity = -0.3f;

    // Lamé parameters for stress-strain relationship
    const float elastic_lambda = 10.0f;
    const float elastic_mu = 20.0f;
    
    NativeArray<Particle> ps; // particles
    NativeArray<Cell> grid;
    
    // deformation gradient. stored as a separate array to use same rendering code for all demos, but feel free to store this field in the particle struct instead
    NativeArray<float2x2> Fs;

    float2[] weights = new float2[3];

    int num_particles;
    List<float2> temp_positions;
    
    SimRenderer sim_renderer;
    
    // interaction
    const float mouse_radius = 10;
    bool mouse_down = false;
    float2 mouse_pos;

    void spawn_box(int x, int y, int box_x = 8, int box_y = 8) {
        const float spacing = 0.5f;
        for (float i = - box_x / 2; i < box_x / 2; i += spacing) {
            for (float j = - box_y / 2; j < box_y / 2; j += spacing) {
                var pos = math.float2(x + i, y + j);

                temp_positions.Add(pos);
            }
        }
    }
    
    void Start () {
        // populate our array of particles
        temp_positions = new List<float2>();
        spawn_box(grid_res / 2, grid_res / 2, 32, 32);
        num_particles = temp_positions.Count;

        ps = new NativeArray<Particle>(num_particles, Allocator.Persistent);
        Fs = new NativeArray<float2x2>(num_particles, Allocator.Persistent);

        // initialise particles
        for (int i = 0; i < num_particles; ++i) {
            Particle p = new Particle();
            p.x = temp_positions[i];
            p.v = 0;
            p.C = 0;
            p.mass = 1.0f;
            ps[i] = p;

            // deformation gradient initialised to the identity
            Fs[i] = math.float2x2(
                1, 0, 
                0, 1
            );
        }

        grid = new NativeArray<Cell>(num_cells, Allocator.Persistent);

        for (int i = 0; i < num_cells; ++i) {
            var cell = new Cell();
            cell.v = 0;
            grid[i] = cell;
        }

        // ---- begin precomputation of particle volumes
        // MPM course, equation 152 

        // launch a P2G job to scatter particle mass to the grid
        new Job_P2G() {
            ps = ps,
            Fs = Fs,
            grid = grid,
            num_particles = num_particles
        }.Schedule().Complete();
        
        for (int i = 0; i < num_particles; ++i) {
            var p = ps[i];

            // quadratic interpolation weights
            float2 cell_idx = math.floor(p.x);
            float2 cell_diff = (p.x - cell_idx) - 0.5f;
            weights[0] = 0.5f * math.pow(0.5f - cell_diff, 2);
            weights[1] = 0.75f - math.pow(cell_diff, 2);
            weights[2] = 0.5f * math.pow(0.5f + cell_diff, 2);

            float density = 0.0f;
            // iterate over neighbouring 3x3 cells
            for (int gx = 0; gx < 3; ++gx) {
                for (int gy = 0; gy < 3; ++gy) {
                    float weight = weights[gx].x * weights[gy].y;

                    // map 2D to 1D index in grid
                    int cell_index = ((int)cell_idx.x + (gx - 1)) * grid_res + ((int)cell_idx.y + gy - 1);
                    density += grid[cell_index].mass * weight;
                }
            }

            // per-particle volume estimate has now been computed
            float volume = p.mass / density;
            p.volume_0 = volume;

            ps[i] = p;
        }

        // ---- end precomputation of particle volumes

        // boilerplate rendering code handled elsewhere
        sim_renderer = GameObject.FindObjectOfType<SimRenderer>();
        sim_renderer.Initialise(num_particles, Marshal.SizeOf(new Particle()));
    }

    private void Update() {
        HandleMouseInteraction();

        for (int i = 0; i < iterations; ++i) {
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
        Profiler.BeginSample("P2G");
        new Job_P2G() {
            ps = ps,
            Fs = Fs,
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
            Fs = Fs,
            mouse_down = mouse_down,
            mouse_pos = mouse_pos,
            grid = grid
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
    unsafe struct Job_P2G : IJob {
        public NativeArray<Cell> grid;
        [ReadOnly] public NativeArray<Particle> ps;
        [ReadOnly] public NativeArray<float2x2> Fs;
        [ReadOnly] public int num_particles;
        
        public void Execute() {
            var weights = stackalloc float2[3];

            for (int i = 0; i < num_particles; ++i) {
                var p = ps[i];
                
                float2x2 stress = 0;

                // deformation gradient
                var F = Fs[i];

                var J = math.determinant(F);

                // MPM course, page 46
                var volume = p.volume_0 * J;

                // useful matrices for Neo-Hookean model
                var F_T = math.transpose(F);
                var F_inv_T = math.inverse(F_T);
                var F_minus_F_inv_T = F - F_inv_T;

                // MPM course equation 48
                var P_term_0 = elastic_mu * (F_minus_F_inv_T);
                var P_term_1 = elastic_lambda * math.log(J) * F_inv_T;
                var P = P_term_0 + P_term_1;

                // cauchy_stress = (1 / det(F)) * P * F_T
                // equation 38, MPM course
                stress = (1.0f / J) * math.mul(P, F_T);

                // (M_p)^-1 = 4, see APIC paper and MPM course page 42
                // this term is used in MLS-MPM paper eq. 16. with quadratic weights, Mp = (1/4) * (delta_x)^2.
                // in this simulation, delta_x = 1, because i scale the rendering of the domain rather than the domain itself.
                // we multiply by dt as part of the process of fusing the momentum and force update for MLS-MPM
                var eq_16_term_0 = -volume * 4 * stress * dt;

                // quadratic interpolation weights
                uint2 cell_idx = (uint2)p.x;
                float2 cell_diff = (p.x - cell_idx) - 0.5f;
                weights[0] = 0.5f * math.pow(0.5f - cell_diff, 2);
                weights[1] = 0.75f - math.pow(cell_diff, 2);
                weights[2] = 0.5f * math.pow(0.5f + cell_diff, 2);

                // for all surrounding 9 cells
                for (uint gx = 0; gx < 3; ++gx) {
                    for (uint gy = 0; gy < 3; ++gy) {
                        float weight = weights[gx].x * weights[gy].y;
                        
                        uint2 cell_x = math.uint2(cell_idx.x + gx - 1, cell_idx.y + gy - 1);
                        float2 cell_dist = (cell_x - p.x) + 0.5f;
                        float2 Q = math.mul(p.C, cell_dist);

                        // scatter mass and momentum to the grid
                        int cell_index = (int)cell_x.x * grid_res + (int)cell_x.y;
                        Cell cell = grid[cell_index];

                        // MPM course, equation 172
                        float weighted_mass = weight * p.mass;
                        cell.mass += weighted_mass;

                        // APIC P2G momentum contribution
                        cell.v += weighted_mass * (p.v + Q);

                        // fused force/momentum update from MLS-MPM
                        // see MLS-MPM paper, equation listed after eqn. 28
                        float2 momentum = math.mul(eq_16_term_0 * weight, cell_dist);
                        cell.v += momentum;

                        // total update on cell.v is now:
                        // weight * (dt * M^-1 * p.volume * p.stress + p.mass * p.C)
                        // this is the fused momentum + force from MLS-MPM. however, instead of our stress being derived from the energy density,
                        // i use the weak form with cauchy stress. converted:
                        // p.volume_0 * (dΨ/dF)(Fp)*(Fp_transposed)
                        // is equal to p.volume * σ

                        // note: currently "cell.v" refers to MOMENTUM, not velocity!
                        // this gets converted in the UpdateGrid step below.

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

                // 'slip' boundary conditions
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
        public NativeArray<float2x2> Fs;
        [ReadOnly] public NativeArray<Cell> grid;

        [ReadOnly] public bool mouse_down;
        [ReadOnly] public float2 mouse_pos;
        
        public void Execute(int i) {
            Particle p = ps[i];

            // reset particle velocity. we calculate it from scratch each step using the grid
            p.v = 0;

            // quadratic interpolation weights
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

                    // APIC paper equation 10, constructing inner term for B
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
            
            // mouse interaction
            if (mouse_down) {
                var dist = p.x - mouse_pos;
                if (math.dot(dist, dist) < mouse_radius * mouse_radius) {
                    float norm_factor = (math.length(dist) / mouse_radius);
                    norm_factor = math.pow(math.sqrt(norm_factor), 8);
                    var force = math.normalize(dist) * norm_factor * 0.5f;
                    p.v += force;
                }
            }

            // deformation gradient update - MPM course, equation 181
            // Fp' = (I + dt * p.C) * Fp
            var Fp_new = math.float2x2(
                1, 0, 
                0, 1
            );
            Fp_new += dt * p.C;
            Fs[i] = math.mul(Fp_new, Fs[i]);

            ps[i] = p;
        }
    }

    #endregion

    private void OnDestroy() {
        ps.Dispose();
        grid.Dispose();
        Fs.Dispose();
    }
}

