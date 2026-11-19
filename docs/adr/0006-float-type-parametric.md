# Float type is a parameter

Date: 2026-11-19. Status: accepted.

Every array and kernel is parametric in the float type `T`, default `Float64`. The development GPU (RTX 5070) runs Float64 at about 1/64 of its Float32 rate; cluster GPUs (A100, H100, MI250) run Float64 at full or half rate. Blocking tests run in Float64. Each module also has a Float32 smoke test, so local GPU work stays fast.

Write literal constants as `T(0.5)` or `one(T)/2` so that no kernel promotes to Float64 by accident.
