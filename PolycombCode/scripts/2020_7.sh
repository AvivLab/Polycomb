cd /polycomb
num=2020_7
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.4 1.0 6.0 0.15 2020_7
python3 upload.py polycomb $num results/$num
