cd /polycomb
num=2020_113
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 2.0 6.0 0.1 2020_113
python3 upload.py polycomb $num results/$num
