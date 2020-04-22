cd /polycomb
num=2020_244
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.5 0.5 1.0 0.2 2020_244python3 upload.py polycomb $num results/$num
