cd /polycomb
num=2020_76
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.3 1.0 1.0 0.2 2020_76
python3 upload.py polycomb $num results/$num
