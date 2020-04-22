cd /polycomb
num=2020_15
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.2 0.4 1.0 4.0 0.15 2020_15
python3 upload.py polycomb $num results/$num
