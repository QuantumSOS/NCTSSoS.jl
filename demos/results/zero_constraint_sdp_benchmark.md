# Zero-constraint SDP benchmark

Two formulations were compared:

- `dense`: PSD matrix variable with no explicit zero equalities
- `zero_eq`: same PSD matrix variable, plus many constraints `X[i,j] == 0`

That means the PSD cone dimension stays the same. We are only changing the number of linear equalities.

## Feasibility problem

### COSMO

| n | zero fraction | zero equalities | dense median (s) | zero-eq median (s) | zero-eq / dense | dense status | zero-eq status |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| 40 | 0.5 | 390 | 0.0086 | 0.0131 | 1.53x | OPTIMAL | OPTIMAL |
| 40 | 0.8 | 624 | 0.0086 | 0.0167 | 1.94x | OPTIMAL | OPTIMAL |
| 80 | 0.5 | 1580 | 0.1151 | 0.1844 | 1.60x | OPTIMAL | OPTIMAL |
| 80 | 0.8 | 2528 | 0.1155 | 0.2412 | 2.09x | OPTIMAL | OPTIMAL |
| 120 | 0.5 | 3570 | 0.5608 | 0.9149 | 1.63x | OPTIMAL | OPTIMAL |
| 120 | 0.8 | 5712 | 0.5638 | 1.1993 | 2.13x | OPTIMAL | OPTIMAL |

Geometric-mean slowdown from adding zero equalities: 1.80x

Raw solve times by seed:

- n=40, zero_fraction=0.5: dense=[0.008684873580932617, 0.008589982986450195, 0.008530855178833008], zero_eq=[0.013161897659301758, 0.01314687728881836, 0.013102054595947266]
- n=40, zero_fraction=0.8: dense=[0.008612871170043945, 0.008586883544921875, 0.008656978607177734], zero_eq=[0.01627206802368164, 0.016811132431030273, 0.016736984252929688]
- n=80, zero_fraction=0.5: dense=[0.11546897888183594, 0.11482596397399902, 0.11512303352355957], zero_eq=[0.18443012237548828, 0.1851949691772461, 0.18387508392333984]
- n=80, zero_fraction=0.8: dense=[0.11562395095825195, 0.11453700065612793, 0.11554503440856934], zero_eq=[0.24252891540527344, 0.24118304252624512, 0.2405409812927246]
- n=120, zero_fraction=0.5: dense=[0.560776948928833, 0.5633690357208252, 0.5595378875732422], zero_eq=[0.914160966873169, 0.9149439334869385, 0.9160370826721191]
- n=120, zero_fraction=0.8: dense=[0.5634341239929199, 0.566906213760376, 0.5638289451599121], zero_eq=[1.201275110244751, 1.1992759704589844, 1.1980469226837158]

### MosekTools

| n | zero fraction | zero equalities | dense median (s) | zero-eq median (s) | zero-eq / dense | dense status | zero-eq status |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| 40 | 0.5 | 390 | 0.0010 | 0.0042 | 4.26x | OPTIMAL | OPTIMAL |
| 40 | 0.8 | 624 | 0.0010 | 0.0076 | 7.95x | OPTIMAL | OPTIMAL |
| 80 | 0.5 | 1580 | 0.0059 | 0.0495 | 8.38x | OPTIMAL | OPTIMAL |
| 80 | 0.8 | 2528 | 0.0056 | 0.1228 | 21.78x | OPTIMAL | OPTIMAL |
| 120 | 0.5 | 3570 | 0.0098 | 0.2553 | 26.15x | OPTIMAL | OPTIMAL |
| 120 | 0.8 | 5712 | 0.0103 | 0.7367 | 71.78x | OPTIMAL | OPTIMAL |

Geometric-mean slowdown from adding zero equalities: 15.04x

Raw solve times by seed:

- n=40, zero_fraction=0.5: dense=[0.001110076904296875, 0.0009779930114746094, 0.0009410381317138672], zero_eq=[0.004146099090576172, 0.004164934158325195, 0.004585981369018555]
- n=40, zero_fraction=0.8: dense=[0.00092315673828125, 0.0009570121765136719, 0.0009639263153076172], zero_eq=[0.007603883743286133, 0.007775068283081055, 0.007421016693115234]
- n=80, zero_fraction=0.5: dense=[0.004045009613037109, 0.005908012390136719, 0.006565093994140625], zero_eq=[0.04952096939086914, 0.049440860748291016, 0.052696943283081055]
- n=80, zero_fraction=0.8: dense=[0.005639076232910156, 0.006467103958129883, 0.005081892013549805], zero_eq=[0.1264050006866455, 0.12172794342041016, 0.12280416488647461]
- n=120, zero_fraction=0.5: dense=[0.009353876113891602, 0.011749982833862305, 0.009763002395629883], zero_eq=[0.25530505180358887, 0.2515280246734619, 0.25789785385131836]
- n=120, zero_fraction=0.8: dense=[0.01026296615600586, 0.011940956115722656, 0.008986949920654297], zero_eq=[0.7366650104522705, 0.7429800033569336, 0.7299678325653076]

## Random linear objective

### COSMO

| n | zero fraction | zero equalities | dense median (s) | zero-eq median (s) | zero-eq / dense | dense status | zero-eq status |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| 40 | 0.5 | 390 | 0.0206 | 0.0804 | 3.91x | OPTIMAL | OPTIMAL |
| 40 | 0.8 | 624 | 0.0219 | 0.0523 | 2.38x | OPTIMAL | OPTIMAL |
| 60 | 0.5 | 885 | 0.0839 | 0.1512 | 1.80x | OPTIMAL | OPTIMAL |
| 60 | 0.8 | 1416 | 0.0812 | 0.1250 | 1.54x | OPTIMAL | OPTIMAL |

Geometric-mean slowdown from adding zero equalities: 2.25x

Raw solve times by seed:

- n=40, zero_fraction=0.5: dense=[0.02055501937866211, 0.023324966430664062, 0.01968693733215332], zero_eq=[0.0803520679473877, 0.03381490707397461, 0.1409299373626709]
- n=40, zero_fraction=0.8: dense=[0.021946191787719727, 0.023561954498291016, 0.019665002822875977], zero_eq=[0.06514596939086914, 0.05229687690734863, 0.04482078552246094]
- n=60, zero_fraction=0.5: dense=[0.09877705574035645, 0.06991791725158691, 0.0839240550994873], zero_eq=[0.15116500854492188, 0.2479081153869629, 0.1506180763244629]
- n=60, zero_fraction=0.8: dense=[0.09879088401794434, 0.0680840015411377, 0.08115005493164062], zero_eq=[0.136152982711792, 0.11983585357666016, 0.12498593330383301]

### MosekTools

| n | zero fraction | zero equalities | dense median (s) | zero-eq median (s) | zero-eq / dense | dense status | zero-eq status |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |
| 40 | 0.5 | 390 | 0.0022 | 0.0129 | 5.90x | OPTIMAL | OPTIMAL |
| 40 | 0.8 | 624 | 0.0021 | 0.0196 | 9.14x | OPTIMAL | OPTIMAL |
| 60 | 0.5 | 885 | 0.0045 | 0.0456 | 10.13x | OPTIMAL | OPTIMAL |
| 60 | 0.8 | 1416 | 0.0047 | 0.0653 | 13.92x | OPTIMAL | OPTIMAL |

Geometric-mean slowdown from adding zero equalities: 9.34x

Raw solve times by seed:

- n=40, zero_fraction=0.5: dense=[0.002398967742919922, 0.0021889209747314453, 0.002048015594482422], zero_eq=[0.014641046524047852, 0.011039972305297852, 0.012911796569824219]
- n=40, zero_fraction=0.8: dense=[0.0021209716796875, 0.002146005630493164, 0.0021588802337646484], zero_eq=[0.02220892906188965, 0.016403913497924805, 0.019622087478637695]
- n=60, zero_fraction=0.5: dense=[0.004502058029174805, 0.004573822021484375, 0.00442194938659668], zero_eq=[0.044883012771606445, 0.04873800277709961, 0.045590877532958984]
- n=60, zero_fraction=0.8: dense=[0.004694223403930664, 0.004799842834472656, 0.004436016082763672], zero_eq=[0.07048702239990234, 0.06031298637390137, 0.0653371810913086]

## Takeaway

In these JuMP benchmarks, explicit zero equalities did not speed up solves. They usually slowed them down, sometimes by a lot. That matches the solver mechanics: the PSD block is still the same size, but the model has more linear constraints to process.

If you want real speedups from sparsity, you generally need a formulation that actually reduces cone size or exploits chordal/sparse structure, not just `X[i,j] == 0` piled on top of a dense PSD variable.
