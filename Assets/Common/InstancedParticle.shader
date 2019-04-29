  Shader "Instanced/InstancedParticle" {
    Properties {
        _Colour ("Colour", COLOR) = (1, 1, 1, 1)
        _Size ("Size", float) = 0.035
    }

    SubShader {
        Pass {
            Tags { "LightMode"="ForwardBase" "Queue" = "Transparent" "RenderType" = "Transparent" "IgnoreProjector" = "True" }
            ZWrite Off
            Blend SrcAlpha OneMinusSrcAlpha

            CGPROGRAM

            #pragma glsl
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_fwdbase nolightmap nodirlightmap nodynlightmap novertexlight
            #pragma target 4.5
            
            #include "UnityCG.cginc"

            // matches the structure of our data on the CPU side
            struct Particle {
                float2 x;
                float2 v;
                float4 C;
                float mass;
                float padding;
            };

            struct v2f {
                float4 pos : SV_POSITION;
            };

            float _Size;
            fixed4 _Colour;

            StructuredBuffer<Particle> particle_buffer;

            v2f vert (appdata_full v, uint instanceID : SV_InstanceID) {
                // take in data from the compute buffer, filled with data each frame in SimRenderer
                // offsetting and scaling it from the (0...grid_res, 0...grid_res) resolution of our sim into a nicer range for rendering
                float4 data = float4((particle_buffer[instanceID].x.xy - float2(32, 32)) * 0.1, 0, 1.0);
                
                // Scaling vertices by our base size param (configurable in the material) and the mass of the particle
                float3 localPosition = v.vertex.xyz * (_Size * data.w);
                float3 worldPosition = data.xyz + localPosition;
				
                // project into camera space
                v2f o;
                o.pos = mul(UNITY_MATRIX_VP, float4(worldPosition, 1.0f));
                
                return o;
            }

            fixed4 frag (v2f i) : SV_Target {
                return _Colour;
            }

            ENDCG
        }
    }
}