# Scalar-vs-matrix PSD benchmark

This report compares direct PSD matrix variables against scalar-built PSD matrices.

Modeling styles:
- `matrix_dense`: `@variable(model, X[1:n,1:n], PSD)`
- `matrix_zeroeq`: same, plus many `X[i,j] == 0`
- `scalar_dense`: scalar variables assembled into a symmetric matrix, then `M in PSDCone()`
- `scalar_zeroeq`: same scalar-built matrix, plus many scalar `== 0` constraints
- `scalar_structural_zero`: only create scalars for allowed nonzeros, fill forbidden entries with literal zeros, then `M in PSDCone()`

## COSMO, n=40, zero_fraction=0.8, objective=linear

| formulation | median solve time (s) | statuses | raw times |
| :--- | ---: | :--- | :--- |
| matrix_dense | 0.0195 | OPTIMAL | `[0.019515037536621094, 0.0224759578704834, 0.019036054611206055]` |
| matrix_zeroeq | 0.0492 | OPTIMAL | `[0.06154608726501465, 0.049223899841308594, 0.041337013244628906]` |
| scalar_dense | 0.0223 | OPTIMAL | `[0.018245935440063477, 0.03403306007385254, 0.022305965423583984]` |
| scalar_zeroeq | 0.0539 | OPTIMAL | `[0.06934905052185059, 0.05387687683105469, 0.03965497016906738]` |
| scalar_structural_zero | 0.0602 | OPTIMAL | `[0.1351299285888672, 0.06023406982421875, 0.048326969146728516]` |

## COSMO, n=80, zero_fraction=0.8, objective=feasibility

| formulation | median solve time (s) | statuses | raw times |
| :--- | ---: | :--- | :--- |
| matrix_dense | 0.1104 | OPTIMAL | `[0.11394596099853516, 0.11039400100708008, 0.109375]` |
| matrix_zeroeq | 0.2315 | OPTIMAL | `[0.23083019256591797, 0.2315359115600586, 0.23221397399902344]` |
| scalar_dense | 0.1105 | OPTIMAL | `[0.11053705215454102, 0.1120610237121582, 0.10994911193847656]` |
| scalar_zeroeq | 0.2302 | OPTIMAL | `[0.23031091690063477, 0.22975921630859375, 0.23016691207885742]` |
| scalar_structural_zero | 0.1315 | OPTIMAL | `[0.13146114349365234, 0.142503023147583, 0.10260820388793945]` |

## MosekTools, n=40, zero_fraction=0.8, objective=linear

| formulation | median solve time (s) | statuses | raw times |
| :--- | ---: | :--- | :--- |
| matrix_dense | 0.0022 | OPTIMAL | `[0.0023391246795654297, 0.0021200180053710938, 0.002206087112426758]` |
| matrix_zeroeq | 0.0202 | OPTIMAL | `[0.022310972213745117, 0.01669907569885254, 0.020190000534057617]` |
| scalar_dense | 0.0278 | OPTIMAL | `[0.02890491485595703, 0.027830839157104492, 0.026405811309814453]` |
| scalar_zeroeq | 0.0346 | OPTIMAL | `[0.03604888916015625, 0.03462505340576172, 0.032195091247558594]` |
| scalar_structural_zero | 0.0340 | OPTIMAL | `[0.03585505485534668, 0.03396892547607422, 0.03266310691833496]` |

## MosekTools, n=80, zero_fraction=0.8, objective=feasibility

| formulation | median solve time (s) | statuses | raw times |
| :--- | ---: | :--- | :--- |
| matrix_dense | 0.0041 | OPTIMAL | `[0.0040531158447265625, 0.004045009613037109, 0.004532814025878906]` |
| matrix_zeroeq | 0.1059 | OPTIMAL | `[0.11158013343811035, 0.10590386390686035, 0.10573005676269531]` |
| scalar_dense | 0.3285 | OPTIMAL | `[0.3285188674926758, 0.34251904487609863, 0.3275721073150635]` |
| scalar_zeroeq | 0.3506 | OPTIMAL | `[0.3502500057220459, 0.3531169891357422, 0.35059309005737305]` |
| scalar_structural_zero | 0.3530 | OPTIMAL | `[0.3504939079284668, 0.3530099391937256, 0.3535299301147461]` |

## Takeaway

Declaring scalar variables first and then setting many of them to zero does **not** buy speed by itself. In these tests it behaves about the same as the matrix-variable version for COSMO, and it is often worse for Mosek. The main reason is unchanged: the PSD cone dimension is still the same, and explicit zero equalities just add more linear constraints.

Even the `scalar_structural_zero` version does not magically fix this. It removes some scalar variables, but the PSD cone is still an `n × n` cone, so the win is limited unless you reformulate the problem to expose smaller cones or chordal structure.
