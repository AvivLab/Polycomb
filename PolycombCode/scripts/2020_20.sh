cd /polycomb
num=2020_20
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.3 1.0 1.0 0.15 2020_20
python3 upload.py polycomb $num results/$num
