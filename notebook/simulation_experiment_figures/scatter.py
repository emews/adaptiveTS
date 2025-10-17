# Simplified scene for debugging - v3

# Camera positioned directly in front of the object
cam_loc = (0.5, 0.5, -3)

scene=f"""
#version 3.7;
global_settings {{ assumed_gamma 1.0 }}

#include "colors.inc"

camera {{
  location <{cam_loc[0]}, {cam_loc[1]}, {cam_loc[2]}>
  look_at  <0.5, 0.5, 0.5>
  sky      <0, 1, 0>
  right    <1.333, 0, 0>
  angle    45
}}
// Light source placed behind the camera
light_source {{ <0.5, 0.5, -5> color White }}
background {{ color SkyBlue }}

// A single large sphere at the center of the look_at point
sphere {{ <0.5, 0.5, 0.5>, 0.5 pigment {{ color Red }} }}

"""

with open("notebook/simulation_experiment_figures/scatter.pov","w") as f:
    f.write(scene)

print("notebook/simulation_experiment_figures/scatter.pov written (simplified test scene v3)")
