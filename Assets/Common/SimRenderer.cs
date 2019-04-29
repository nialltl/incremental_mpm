using UnityEngine;
using System.Collections;
using Unity.Collections;

public class SimRenderer : MonoBehaviour {
    [SerializeField] Mesh instance_mesh;
    [SerializeField] Material instance_material;

    // Compute buffers used for indirect mesh drawing
    ComputeBuffer point_buffer;
    ComputeBuffer args_buffer;
    uint[] args = new uint[5] { 0, 0, 0, 0, 0 };
    
    Bounds bounds;

    public void Initialise(int buffer_size, int buffer_element_size) {
        // Create a compute buffer that holds (# of particles), with a byte offset of (size of particle) for the GPU to process
        point_buffer = new ComputeBuffer(buffer_size, buffer_element_size, ComputeBufferType.Default);
        
        // ensure the material has a reference to this consistent buffer so the shader can access it
        instance_material.SetBuffer("particle_buffer", point_buffer);
        
        // indirect arguments for mesh instances
        args_buffer = new ComputeBuffer(1, args.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        uint numIndices = (uint)instance_mesh.GetIndexCount(0);
        args[0] = numIndices;
        args[1] = (uint)point_buffer.count;
        args_buffer.SetData(args);

        // define rendering bounds for DrawMeshInstancedIndirect - basically just arbitrarily big enough not to get clipped
        bounds = new Bounds(Vector3.zero, new Vector3(100, 100, 100));
    }

    public void RenderFrame<T>(NativeArray<T> ps) where T : struct {
        // fill compute buffer with all data of our particles NativeArray, gets passed to our instanced material
        point_buffer.SetData(ps);

        Graphics.DrawMeshInstancedIndirect(instance_mesh, 0, instance_material, bounds, args_buffer);
    }
    
    void OnDisable() {
        if (args_buffer != null) args_buffer.Release();
        if (point_buffer != null) point_buffer.Release();
    }
}