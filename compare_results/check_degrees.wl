e = Get["/scratch/mella/pipeline_to_copy/dist/out/phys_c1_w6_ord18.m"];
Print["ord18: max deg x=", Max[Exponent[#, x] & /@ e], "  max deg y=", Max[Exponent[#, y] & /@ e]];
e2 = Get["/scratch/mella/pipeline_to_copy/dist/out/phys_c1_w6_ord19.m"];
Print["ord19: max deg x=", Max[Exponent[#, x] & /@ e2], "  max deg y=", Max[Exponent[#, y] & /@ e2]];
