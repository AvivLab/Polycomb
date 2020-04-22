cd /polycomb
num=2020_53
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.5 1.0 1.0 0.1 2020_53
python3 upload.py polycomb $num results/$num
