cd /polycomb
num=2020_115
echo "Running job for $num"
julia src/phenotypicPliancyAnalyses.jl 0.3 0.4 2.0 6.0 0.1 2020_115
python3 upload.py polycomb $num results/$num
