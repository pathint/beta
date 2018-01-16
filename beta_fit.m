(* Functions *) 
(* Frequencies to Proportions *)
ProcData[file_] := 
  If[EvenQ[Length[#]], 
     Table[{i, N[(#[[i + 1]] + #[[-(i + 1)]])/Total[#]]}, {i, 0, 
       Length[#]/2 - 1}], 
     Append[Table[{i, N[(#[[i + 1]] + #[[-(i + 1)]])/Total[#]]}, {i, 
        0, (Length[#] - 1)/2 - 1}], {(Length[#] - 1)/2, 
       N[#[[(Length[#] + 1)/2]]/Total[#]]}]] &@
   Flatten[Import[file, "Table"]];
CombineRule[values__] := 
  Module[{nonzeros, minNonZeros = 2}, 
   nonzeros = DeleteCases[List[values], 0];
   If[Length[nonzeros] >= minNonZeros, Median[nonzeros], 0]];
CombineExpressionProfiles[samples_] := MapThread[CombineRule, samples];

(* Do the Integral *)
BetaFit[m_, alpha_] := 
  If[EvenQ[m], 
     Append[#, {m/2, 
       Integrate[
        PDF[BinomialDistribution[m, p], m/2] PDF[
          BetaDistribution[alpha, alpha], p], {p, 0, 1}]}], #] &@
   Table[{i, 
     Integrate[(PDF[BinomialDistribution[m, p], i] + 
         PDF[BinomialDistribution[m, p], m - i]) PDF[
        BetaDistribution[alpha, alpha], p], {p, 0, 1}]}, {i, 
     0, (m - 1)/2}];

(* Example Usage *)

result = ProcData["GSE92332_result_" <> ToString[#] <> ".txt"] & /@
    Range[1, 15];

model = BetaFit[36, 0.242]
ListPlot[{result[[1]], model}, Joined -> False, Filling -> Axis, 
 PlotRange -> All]

