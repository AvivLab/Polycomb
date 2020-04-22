cd /polycomb
num=2020_5
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.1 0.4 1.0 6.0 0.15 2020_5
python3 upload.py polycomb $num results/$num
